"""
04_build_inputs.py
==================
Purpose : Turn the downloaded data into the exact inputs NetMedPy needs:
            disease -> set of genes           (genetic evidence only, score >= cutoff)
            drug    -> set of target genes    (approved drugs)
            known drug-disease pairs          (approved indications = the benchmark)
          for each STRING network (conf700, conf900), in two versions:
            'full'   : everything present in that network
            'shared' : restricted to proteins present in ALL networks (fair comparison)
          and gene-symbol-labelled network files for NetMedPy.
Input   : data/processed/*  (script 03), data/networks/*.pkl, results/physical_shared_proteins.csv
Output  : data/networks/physical_conf{700,900}_symbols.pkl   (nodes are gene symbols)
          data/processed/inputs_<network>_<mode>.pkl
          results/input_coverage.csv
          results/known_pairs_named.csv   (every benchmark pair with drug name + targets, to eyeball)
          results/dropped_drugs.csv       (approved drugs excluded, and why)
"""

import os
import pickle
from pathlib import Path

import networkx as nx
import pandas as pd
import pyarrow.dataset as pads

# ------------------------------------------------------------------
# 1. SETTINGS
# ------------------------------------------------------------------
BASE_DIR = Path(os.environ.get("NDR_BASE", (
    r"C:\Users\Vyomini\Documents\MU Engineering_Btech [CB]"
    r"\Projects\Course_Projects\Sem7_Projects"
    r"\Sys Bio & Network Model Project"
    r"\Network_Drug_Repurposing")))
OT_DIR, NET_DIR = BASE_DIR / "data/raw/opentargets", BASE_DIR / "data/networks"
PROC_DIR, RES_DIR = BASE_DIR / "data/processed", BASE_DIR / "results"

SCORE_CUTOFF = 0.5          # genetic association score cut-off for disease genes
MIN_GENES = 20              # a disease needs at least this many genes in the network
NETWORKS = ["conf400", "conf700", "conf900"]
STUDY_DISEASES = {          # exact Open Targets name -> group
    "type 2 diabetes mellitus": "A_cardiometabolic",
    "heart failure": "A_cardiometabolic",
    "hypercholesterolemia": "A_cardiometabolic",
    "rheumatoid arthritis": "B_immune",
    "asthma": "B_immune",
    "psoriasis": "B_immune",
}
GENETIC_TYPES = {"genetic_association"}

# ------------------------------------------------------------------
# 2. LOAD PROCESSED TABLES + LEAKAGE GUARD
# ------------------------------------------------------------------
assoc = pd.read_parquet(PROC_DIR / "genetic_associations.parquet")
assert set(assoc["aggregationValue"].unique()) <= GENETIC_TYPES, \
    "STOP: non-genetic evidence found in disease-gene table (leakage)"
print("Leakage guard passed. Evidence types in disease genes:", sorted(assoc["aggregationValue"].unique()))

drug_targets = pd.read_csv(PROC_DIR / "drug_targets.csv").dropna(subset=["symbol"])
approved = pd.read_csv(PROC_DIR / "approved_indications.csv")

# Salt/formulation IDs (e.g. SITAGLIPTIN FUMARATE) are children of a parent compound
# (SITAGLIPTIN). Merge them to the parent so (1) the same molecule is not counted twice
# and (2) a drug whose target is recorded on the parent/child is not lost.
drug_info = pads.dataset(str(OT_DIR / "drug_molecule"), format="parquet").to_table(
    columns=["id", "name", "drugType", "parentId"]).to_pandas().set_index("id")
canon = {d: (p if isinstance(p, str) and p else d)
         for d, p in zip(drug_info.index, drug_info["parentId"])}
n_children = sum(1 for d, c in canon.items() if d != c)
approved["drugId"] = approved["drugId"].map(lambda d: canon.get(d, d))
approved = approved.drop_duplicates()
drug_targets["drugId"] = drug_targets["drugId"].map(lambda d: canon.get(d, d))
drug_targets = drug_targets.drop_duplicates()
print(f"Drug IDs harmonised to parent compounds ({n_children:,} child IDs in Open Targets)")

disease = pads.dataset(str(OT_DIR / "disease"), format="parquet").to_table(
    columns=["id", "name", "descendants"]).to_pandas()
name_to_ids, desc = {}, {}
for did, nm, ds in zip(disease["id"], disease["name"], disease["descendants"]):
    name_to_ids.setdefault(str(nm).lower(), []).append(did)
    desc[did] = list(ds) if ds is not None else []

# ------------------------------------------------------------------
# 3. DISEASE GENES AND KNOWN PAIRS (independent of any network)
# ------------------------------------------------------------------
disease_genes_all, known_rows = {}, []
for name in STUDY_DISEASES:
    ids = name_to_ids.get(name.lower(), [])
    if not ids:
        print(f"  ! disease '{name}' not found - skipped"); continue
    dset = set(ids)
    for i in ids:
        dset.update(desc.get(i, []))                      # include subtypes
    a = assoc[assoc["diseaseId"].isin(dset)].dropna(subset=["symbol"])
    best = a.groupby("symbol")["associationScore"].max()
    disease_genes_all[name] = set(best[best >= SCORE_CUTOFF].index)
    for drug in set(approved[approved["diseaseId"].isin(dset)]["drugId"]):
        known_rows.append((drug, name))

known_all = pd.DataFrame(known_rows, columns=["drugId", "disease"]).drop_duplicates()
drug_targets_all = drug_targets.groupby("drugId")["symbol"].apply(set).to_dict()
dname = drug_info["name"].to_dict()
dtype = drug_info["drugType"].to_dict()
dropped_rows, named_frames = [], []
print(f"\nStudy diseases found: {len(disease_genes_all)}   known approved pairs: {len(known_all)}")

# ------------------------------------------------------------------
# 4. NETWORKS: relabel nodes by gene symbol, then build inputs per network/mode
# ------------------------------------------------------------------
shared_symbols = set(pd.read_csv(RES_DIR / "physical_shared_proteins.csv")["preferred_name"].dropna().str.upper())
coverage = []

for net in NETWORKS:
    with open(NET_DIR / f"physical_{net}_lcc.pkl", "rb") as f:
        g = pickle.load(f)
    mapping = {n: str(d["preferred_name"]).upper() for n, d in g.nodes(data=True)}
    if len(set(mapping.values())) != len(mapping):
        raise SystemExit(f"STOP: duplicate gene symbols in {net}; relabelling would merge proteins.")
    gs = nx.relabel_nodes(g, mapping)                    # nodes are now gene symbols
    with open(NET_DIR / f"physical_{net}_symbols.pkl", "wb") as f:
        pickle.dump(gs, f, protocol=pickle.HIGHEST_PROTOCOL)
    nodes = set(gs.nodes())

    for mode, allowed in (("full", nodes), ("shared", nodes & shared_symbols)):
        dgenes = {d: s & allowed for d, s in disease_genes_all.items()}
        dgenes = {d: s for d, s in dgenes.items() if len(s) >= MIN_GENES}       # keep usable diseases
        # drugs: approved for a kept disease, with >=1 target in the allowed protein set
        kept_pairs = known_all[known_all["disease"].isin(dgenes)]
        dtargets = {}
        for drug in set(kept_pairs["drugId"]):
            t = drug_targets_all.get(drug, set()) & allowed
            if t:
                dtargets[drug] = t
        kept_pairs = kept_pairs[kept_pairs["drugId"].isin(dtargets)].reset_index(drop=True)

        # --- transparency: which approved drugs were lost, and why ---
        wanted = set(known_all[known_all["disease"].isin(dgenes)]["drugId"])
        for drug in sorted(wanted - set(dtargets)):
            why = ("no target recorded in Open Targets" if drug not in drug_targets_all
                   else "its targets are not in this network")
            dropped_rows.append({"network": net, "mode": mode, "drugId": drug,
                                 "drug_name": dname.get(drug), "drugType": dtype.get(drug), "reason": why})
        kp = kept_pairs.copy()
        kp["drug_name"] = kp["drugId"].map(dname)
        kp["drugType"] = kp["drugId"].map(dtype)
        kp["targets_in_network"] = kp["drugId"].map(lambda d: ";".join(sorted(dtargets[d])))
        kp.insert(0, "mode", mode); kp.insert(0, "network", net)
        named_frames.append(kp)

        with open(PROC_DIR / f"inputs_{net}_{mode}.pkl", "wb") as f:
            pickle.dump({"disease_genes": dgenes, "drug_targets": dtargets, "known_pairs": kept_pairs,
                         "drug_names": {d: dname.get(d) for d in dtargets},
                         "score_cutoff": SCORE_CUTOFF, "network": net, "mode": mode}, f)

        for d in disease_genes_all:
            coverage.append({"network": net, "mode": mode, "disease": d,
                             "genes_total": len(disease_genes_all[d]),
                             "genes_used": len(dgenes.get(d, [])),
                             "disease_kept": d in dgenes,
                             "known_drugs_used": int((kept_pairs["disease"] == d).sum())})
        print(f"[{net}/{mode}] diseases kept: {len(dgenes)}  drugs: {len(dtargets)}  "
              f"known pairs: {len(kept_pairs)}  all possible pairs: {len(dgenes) * len(dtargets)}")

cov = pd.DataFrame(coverage)
cov.to_csv(RES_DIR / "input_coverage.csv", index=False)
pd.set_option("display.width", 200)
print("\n" + "=" * 70 + "\nINPUT COVERAGE (score cut-off", SCORE_CUTOFF, ")\n" + "=" * 70)
print(cov.to_string(index=False))
pd.concat(named_frames).to_csv(RES_DIR / "known_pairs_named.csv", index=False)
dropped = pd.DataFrame(dropped_rows)
dropped.to_csv(RES_DIR / "dropped_drugs.csv", index=False)
d7 = dropped[(dropped["network"] == "conf700") & (dropped["mode"] == "full")]
print("\nApproved drugs excluded in conf700/full:", len(d7))
print(d7["reason"].value_counts().to_string())
print("Examples:", ", ".join(map(str, d7["drug_name"].head(12))))
ex = pd.concat(named_frames)
ex = ex[(ex["network"] == "conf700") & (ex["mode"] == "full")].head(10)
print("\nSample of benchmark pairs (check these make biological sense):")
print(ex[["drug_name", "disease", "targets_in_network"]].to_string(index=False))
print("\nScript 04 finished.")