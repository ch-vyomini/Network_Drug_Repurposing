"""
03_coverage_table.py
====================
Purpose : For ~12 candidate diseases, count what data we actually have, so the
          study diseases are chosen from evidence (not by hope).
For each disease it reports:
  - disease genes from GENETIC evidence only (never drug evidence), at 3 score cut-offs,
    and how many of them exist in our STRING networks
  - approved drugs (maxClinicalStage == APPROVAL) and how many have a target in the network
  - target diversity (are all drugs hitting the same one protein?)
Input   : data/raw/opentargets/*  (script 02),  data/networks/*.pkl and
          results/physical_shared_proteins.csv  (script 01)
Output  : results/disease_coverage.csv
          data/processed/genetic_associations.parquet   (filtered, small)
          data/processed/drug_targets.csv, approved_indications.csv
Leakage control: only the evidence types in GENETIC_TYPES are used. The script
          prints every evidence type present and refuses to run if a name is missing.
"""

import os
import pickle
from pathlib import Path

import pandas as pd
import pyarrow.compute as pc
import pyarrow.dataset as pads

# ------------------------------------------------------------------
# 1. SETTINGS
# ------------------------------------------------------------------
BASE_DIR = Path(os.environ.get("NDR_BASE", (
    r"C:\Users\Vyomini\Documents\MU Engineering_Btech [CB]"
    r"\Projects\Course_Projects\Sem7_Projects"
    r"\Sys Bio & Network Model Project"
    r"\Network_Drug_Repurposing")))

OT_DIR      = BASE_DIR / "data" / "raw" / "opentargets"
NETWORK_DIR = BASE_DIR / "data" / "networks"
RESULTS_DIR = BASE_DIR / "results"
PROC_DIR    = BASE_DIR / "data" / "processed"
PROC_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Evidence types allowed for disease genes. 'known_drug' must NEVER be here.
GENETIC_TYPES = ["genetic_association"]
FORBIDDEN_TYPES = {"known_drug"}
assert not (set(GENETIC_TYPES) & FORBIDDEN_TYPES), "Drug evidence must not be used for disease genes"

SCORE_CUTOFFS = [0.1, 0.3, 0.5]          # association-score cut-offs to compare
APPROVED_LABEL = "APPROVAL"              # value of maxClinicalStage for approved drugs

# Candidate diseases (exact Open Targets names, case-insensitive). Edit freely.
CANDIDATES = {
    "A_cardiometabolic": ["type 2 diabetes mellitus", "hypertension", "coronary artery disease",
                          "hypercholesterolemia", "heart failure", "obesity"],
    "B_immune": ["rheumatoid arthritis", "psoriasis", "Crohn disease", "ulcerative colitis",
                 "asthma", "systemic lupus erythematosus", "multiple sclerosis"],
}

# ------------------------------------------------------------------
# 2. NETWORK PROTEIN SYMBOLS (700 network, and the shared set = 900 LCC)
# ------------------------------------------------------------------
def symbols_in_graph(path):
    with open(path, "rb") as f:
        g = pickle.load(f)
    return {str(d.get("preferred_name")).upper() for _, d in g.nodes(data=True) if d.get("preferred_name")}

net700 = symbols_in_graph(NETWORK_DIR / "physical_conf700_lcc.pkl")
shared = set(pd.read_csv(RESULTS_DIR / "physical_shared_proteins.csv")["preferred_name"].dropna().str.upper())
print(f"Network symbols: conf700 LCC = {len(net700):,}   shared (conf900 LCC) = {len(shared):,}")

# ------------------------------------------------------------------
# 3. LOAD OPEN TARGETS TABLES (only the columns we need)
# ------------------------------------------------------------------
def load(name, columns, filt=None):
    d = pads.dataset(str(OT_DIR / name), format="parquet")
    return d.to_table(columns=columns, filter=filt).to_pandas()

# 3a. evidence types present (leakage check)
types = pc.unique(pads.dataset(str(OT_DIR / "association_by_datatype_direct"), format="parquet")
                  .to_table(columns=["aggregationValue"]).column("aggregationValue")).to_pylist()
print("\nEvidence types in the association table:", sorted(map(str, types)))
missing = [t for t in GENETIC_TYPES if t not in types]
if missing:
    raise SystemExit(f"STOP: expected genetic evidence type(s) {missing} not found. "
                     f"Send me the list above and I will fix GENETIC_TYPES.")
print("Using ONLY:", GENETIC_TYPES, "| excluded drug evidence:", sorted(FORBIDDEN_TYPES & set(types)))

# 3b. genetic gene-disease associations (filtered at read time = low memory)
filt = (pc.field("aggregationType") == "datatypeId") & pc.field("aggregationValue").isin(GENETIC_TYPES)
assoc = load("association_by_datatype_direct",
             ["diseaseId", "targetId", "aggregationValue", "associationScore", "evidenceCount"], filt)
print(f"Genetic association rows kept: {len(assoc):,}")

# 3c. target Ensembl -> symbol, disease names/descendants, drugs
targets = load("target", ["id", "approvedSymbol"]).rename(columns={"id": "targetId", "approvedSymbol": "symbol"})
ens2sym = dict(zip(targets["targetId"], targets["symbol"].str.upper()))
assoc["symbol"] = assoc["targetId"].map(ens2sym)

disease = load("disease", ["id", "name", "descendants"])
name_to_ids = {}
for did, nm in zip(disease["id"], disease["name"]):
    name_to_ids.setdefault(str(nm).lower(), []).append(did)
descendants = {did: list(desc) if desc is not None else [] for did, desc in zip(disease["id"], disease["descendants"])}

ind = load("clinical_indication", ["diseaseId", "drugId", "maxClinicalStage"])
ind = ind[ind["maxClinicalStage"] == APPROVED_LABEL][["diseaseId", "drugId"]].drop_duplicates()

moa = load("drug_mechanism_of_action", ["chemblIds", "targets", "actionType", "targetType"])
rows = []
for ids, tg, act, ttype in zip(moa["chemblIds"], moa["targets"], moa["actionType"], moa["targetType"]):
    for drug in (ids if ids is not None else []):
        for t in (tg if tg is not None else []):
            rows.append((drug, t, ens2sym.get(t), act, ttype))
drug_targets = pd.DataFrame(rows, columns=["drugId", "targetId", "symbol", "actionType", "targetType"]).drop_duplicates()
print(f"Approved indication pairs: {len(ind):,}   drug-target pairs: {len(drug_targets):,}")

# save processed tables for script 04
assoc.to_parquet(PROC_DIR / "genetic_associations.parquet", index=False)
drug_targets.to_csv(PROC_DIR / "drug_targets.csv", index=False)
ind.to_csv(PROC_DIR / "approved_indications.csv", index=False)

# ------------------------------------------------------------------
# 4. COVERAGE TABLE
# ------------------------------------------------------------------
dt_net = drug_targets.dropna(subset=["symbol"])
rows = []
for group, names in CANDIDATES.items():
    for nm in names:
        ids = name_to_ids.get(nm.lower(), [])
        if not ids:
            near = [n for n in name_to_ids if nm.lower().split()[0] in n][:5]
            print(f"  ! '{nm}' not found exactly. Similar names: {near}")
            continue
        dset = set(ids)
        for i in ids:
            dset.update(descendants.get(i, []))

        # genes: best score per gene across the disease and its subtypes
        a = assoc[assoc["diseaseId"].isin(dset)].dropna(subset=["symbol"])
        best = a.groupby("symbol")["associationScore"].max()
        row = {"group": group, "disease": nm, "n_ontology_ids": len(dset), "genes_any": len(best)}
        for c in SCORE_CUTOFFS:
            gs = set(best[best >= c].index)
            row[f"genes>={c}_in_net700"] = len(gs & net700)
            row[f"genes>={c}_in_shared"] = len(gs & shared)

        # drugs: approved for the disease or a subtype, with a target in the network
        drugs = set(ind[ind["diseaseId"].isin(dset)]["drugId"])
        sub = dt_net[dt_net["drugId"].isin(drugs)]
        in700 = sub[sub["symbol"].isin(net700)]
        inshared = sub[sub["symbol"].isin(shared)]
        row["approved_drugs"] = len(drugs)
        row["drugs_target_in_net700"] = in700["drugId"].nunique()
        row["drugs_target_in_shared"] = inshared["drugId"].nunique()
        row["distinct_targets_net700"] = in700["symbol"].nunique()
        if len(in700):
            per_target = in700.groupby("symbol")["drugId"].nunique()
            row["top_target_share"] = round(per_target.max() / in700["drugId"].nunique(), 2)
            row["top_target"] = per_target.idxmax()
        rows.append(row)

cov = pd.DataFrame(rows)
cov.to_csv(RESULTS_DIR / "disease_coverage.csv", index=False)
pd.set_option("display.width", 250, "display.max_columns", 40)
print("\n" + "=" * 70 + "\nDISEASE COVERAGE TABLE\n" + "=" * 70)
print(cov.to_string(index=False))
print("\nSaved results/disease_coverage.csv and data/processed/*. Script 03 finished.")