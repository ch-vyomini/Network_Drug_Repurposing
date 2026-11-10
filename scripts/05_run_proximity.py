"""
05_run_proximity.py
===================
Purpose : Apply NetMedPy network proximity (Guney et al. 2016) to every drug x disease
          pair, for ALL the runs we need, one after another (hands-off).
          We do NOT implement the algorithm - we call netmedpy.screening().

What goes in (made by script 04)
  graph        : STRING physical network, nodes = gene symbols
  drug_targets : {drug id -> set of target gene symbols}          (the set T)
  disease_genes: {disease -> set of genetic disease genes}        (the set S)
What netmedpy does
  1. all_pair_distances : shortest-path length between every pair of proteins (cached on disk)
  2. screening          : for each (drug, disease): average, over each drug target, of the
                          distance to the NEAREST disease gene (= 'closest' measure d_c),
                          compared with N_ITER random target/disease sets of the same size and
                          degree  ->  z-score = (d_c - mean_random) / sd_random
What comes out
  results/proximity_<network>_<mode>_<null>_iter<N>.csv : one row per pair:
      z_score, raw_distance, p_value, known (0/1), drug_name, ...

The RUNS list says what to compute. A run whose result file already exists is SKIPPED,
so you can stop and restart safely. Windows note: keep everything under
'if __name__ == "__main__"' (NetMedPy uses worker processes).
"""

import gc
import os
import pickle
import time
from pathlib import Path

import pandas as pd
import netmedpy

# ------------------------------------------------------------------
# 1. SETTINGS
# ------------------------------------------------------------------
BASE_DIR = Path(os.environ.get("NDR_BASE", (
    r"C:\Users\Vyomini\Documents\MU Engineering_Btech [CB]"
    r"\Projects\Course_Projects\Sem7_Projects"
    r"\Sys Bio & Network Model Project"
    r"\Network_Drug_Repurposing")))

# (network, mode). 'full' = everything in that network (primary analysis);
# 'shared' = only proteins present in all three networks (fair threshold comparison).
# conf900/shared is identical to conf900/full, so it is listed once.
RUNS = [("conf900", "shared"),
        ("conf700", "full"),
        ("conf700", "shared"),
        ("conf400", "shared")]

NULL_MODEL = "log_binning"   # degree-binned random sets, as in the paper (bins of >= 100 nodes)
BIN_SIZE   = 100
N_ITER     = 1000            # random reference sets per pair (the paper used 1,000)
N_PROCS    = max(1, min(4, (os.cpu_count() or 2) - 1))   # lower to 2 if memory runs out
QUICK_TEST = os.environ.get("NDR_QUICK", "0") == "1"      # 1 disease, few drugs, 100 iterations
PROXIMAL_Z = -0.15           # the paper's cut-off for calling a pair 'proximal'

NET_DIR, PROC_DIR, RES_DIR = BASE_DIR / "data/networks", BASE_DIR / "data/processed", BASE_DIR / "results"


def auc_lower_is_positive(z, y):
    """AUC for 'known' pairs having LOWER z than unknown pairs (rank-sum formula)."""
    ranks = pd.Series(-z.values).rank()              # high rank = very low z
    n1, n0 = int(y.sum()), int((1 - y).sum())
    return (ranks[y.values == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def run_one(network, mode):
    n_iter = 100 if QUICK_TEST else N_ITER
    tag = f"{network}_{mode}_{NULL_MODEL}_iter{n_iter}" + ("_QUICK" if QUICK_TEST else "")
    out = RES_DIR / f"proximity_{tag}.csv"
    if out.exists():
        print(f"\n[{network}/{mode}] result already exists -> skipped ({out.name})")
        return
    t0 = time.time()
    print("\n" + "=" * 70 + f"\nRUN: {network} / {mode}   (processes: {N_PROCS})\n" + "=" * 70)

    with open(NET_DIR / f"physical_{network}_symbols.pkl", "rb") as f:
        graph = pickle.load(f)
    with open(PROC_DIR / f"inputs_{network}_{mode}.pkl", "rb") as f:
        inp = pickle.load(f)
    drug_targets, disease_genes = inp["drug_targets"], inp["disease_genes"]
    known = set(map(tuple, inp["known_pairs"][["drugId", "disease"]].values))

    if QUICK_TEST:
        first = sorted(disease_genes)[0]
        disease_genes = {first: disease_genes[first]}
        keep = sorted({d for d, ds in known if ds == first})[:6] + sorted(drug_targets)[:6]
        drug_targets = {d: drug_targets[d] for d in dict.fromkeys(keep)}
        print(">>> QUICK TEST: 1 disease,", len(drug_targets), "drugs, 100 iterations")

    print(f"Network: {graph.number_of_nodes():,} proteins, {graph.number_of_edges():,} edges")
    print(f"Pairs to score: {len(drug_targets)} drugs x {len(disease_genes)} diseases "
          f"= {len(drug_targets) * len(disease_genes):,}")
    for name, sets in (("drug targets", drug_targets.values()), ("disease genes", disease_genes.values())):
        missing = set().union(*sets) - set(graph.nodes())
        assert not missing, f"{name} not in network: {list(missing)[:5]}"

    # distances between all protein pairs (computed once per network, then cached)
    dist_file = PROC_DIR / f"dist_{network}.npz"
    if dist_file.exists():
        D = netmedpy.load_distances(str(dist_file))
        print("Loaded cached distance matrix:", dist_file.name)
    else:
        print("Computing all-pairs shortest paths (one-off; minutes)...")
        D = netmedpy.all_pair_distances(graph, distance="shortest_path", n_processors=N_PROCS)
        netmedpy.save_distances(D, str(dist_file))
        print(f"  done in {time.time() - t0:.0f}s, matrix {D.get_matrix().shape}")

    t1 = time.time()
    res = netmedpy.screening(
        drug_targets, disease_genes, graph, D, score="proximity",
        properties=["z_score", "raw_amspl", "p_value_single_tail"],
        null_model=NULL_MODEL, bin_size=BIN_SIZE, n_iter=n_iter, n_procs=N_PROCS)
    print(f"\nScreening took {time.time() - t1:.0f}s for {len(drug_targets) * len(disease_genes):,} pairs")

    long = res["z_score"].stack().rename("z_score").reset_index()
    long.columns = ["drug", "disease", "z_score"]
    long["raw_distance"] = res["raw_amspl"].stack().values
    long["p_value"] = res["p_value_single_tail"].stack().values
    long["known"] = [int((d, ds) in known) for d, ds in zip(long["drug"], long["disease"])]
    long["drug_name"] = long["drug"].map(inp.get("drug_names", {}))
    long["n_targets"] = long["drug"].map(lambda d: len(inp["drug_targets"][d]))
    long["n_disease_genes"] = long["disease"].map(lambda s: len(inp["disease_genes"][s]))
    long.to_csv(out, index=False)

    y = long["known"]
    print(f"Saved {out.name}  ({len(long):,} pairs, {int(y.sum())} known)")
    if y.sum() and (1 - y).sum():
        k, u = long[y == 1]["z_score"], long[y == 0]["z_score"]
        print(f"Mean z   known: {k.mean():.2f}   unknown: {u.mean():.2f}   (lower = closer)")
        print(f"Proximal (z <= {PROXIMAL_Z}) known: {(k <= PROXIMAL_Z).mean():.0%}   unknown: {(u <= PROXIMAL_Z).mean():.0%}")
        print(f"Quick AUC (known ranked above unknown): {auc_lower_is_positive(long['z_score'], y):.3f}")
    print(f"Run time {time.time() - t0:.0f}s")
    del D, graph
    gc.collect()


def main():
    t_all = time.time()
    for network, mode in RUNS:
        run_one(network, mode)
    print(f"\nALL RUNS DONE in {(time.time() - t_all) / 60:.1f} min.  Script 05 finished.")


if __name__ == "__main__":
    main()