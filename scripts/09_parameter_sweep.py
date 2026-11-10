"""
09_parameter_sweep.py
=====================
Purpose : Change NetMedPy's OWN parameters on the SAME dataset and see how the results change
          (this tests understanding of the algorithm, not just of the data).
Dataset : the primary analysis (conf700 / full): 246 drugs x 6 diseases = 1,476 pairs (from script 04).
Baseline: results/proximity_conf700_full_log_binning_iter1000.csv (script 05) =
          closest distance (T->S), null_model='log_binning', bin_size=100, n_iter=1000.
Variants (ONE parameter changed at a time):
  null_model = degree_match   exact-degree matching instead of degree bins
  null_model = uniform        random sets with NO degree matching  (what the paper's correction is for)
  n_iter     = 100            only 100 random sets instead of 1,000  (more Monte-Carlo noise)
  symmetric  = True           distance = average of T->S and S->T instead of only T->S
  bin_size   = 30             finer degree bins
For each variant: Spearman correlation of z with the baseline, % 'proximal' calls and their agreement (kappa),
  pooled AUC (95% drug-level bootstrap CI) and the paired AUC difference vs the baseline.
Output  : results/parameter_sweep_summary.csv, results/sweep_<variant>.csv, figures/parameter_sweep.png
Re-run  : safe (finished variants are skipped).   Windows: keep under 'if __name__ == "__main__"'.
"""
import os
import pickle
import re
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import netmedpy
from scipy.stats import rankdata, spearmanr

BASE_DIR = Path(os.environ.get("NDR_BASE", Path(__file__).resolve().parents[1]))
NET_DIR, PROC_DIR, RES_DIR = BASE_DIR / "data/networks", BASE_DIR / "data/processed", BASE_DIR / "results"
FIG_DIR = RES_DIR / "figures"; FIG_DIR.mkdir(parents=True, exist_ok=True)
NETWORK, MODE = "conf700", "full"
N_PROCS = max(1, min(4, (os.cpu_count() or 2) - 1))
PROXIMAL_Z = -0.15
N_BOOT = int(os.environ.get("NDR_BOOT", 1000))
FAST_ITER = os.environ.get("NDR_ITER")                     # testing only
BASE = dict(null_model="log_binning", bin_size=100, n_iter=1000, symmetric=False)
VARIANTS = {
    "null_degree_match": dict(null_model="degree_match"),
    "null_uniform": dict(null_model="uniform"),
    "n_iter_100": dict(n_iter=100),
    "symmetric_true": dict(symmetric=True),
    "bin_size_30": dict(bin_size=30),
}


def auc(y, s):
    n1 = int(y.sum()); n0 = len(y) - n1
    r = rankdata(s)
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def kappa(a, b):
    po = np.mean(a == b); pe = a.mean() * b.mean() + (1 - a.mean()) * (1 - b.mean())
    return (po - pe) / (1 - pe)


def main():
    rng = np.random.default_rng(42)
    base = pd.read_csv(RES_DIR / f"proximity_{NETWORK}_{MODE}_log_binning_iter1000.csv").set_index(["drug", "disease"])
    with open(NET_DIR / f"physical_{NETWORK}_symbols.pkl", "rb") as f:
        graph = pickle.load(f)
    with open(PROC_DIR / f"inputs_{NETWORK}_{MODE}.pkl", "rb") as f:
        inp = pickle.load(f)
    D = netmedpy.load_distances(str(PROC_DIR / f"dist_{NETWORK}.npz"))
    y = base["known"].values.astype(int)
    drugs = base.index.get_level_values("drug").values
    disease = base.index.get_level_values("disease").values
    groups = list(pd.Series(np.arange(len(base))).groupby(drugs).indices.values())
    z0 = base["z_score"].values
    print(f"Baseline {NETWORK}/{MODE}: {len(base):,} pairs, {y.sum()} known; baseline AUC = {auc(y, -z0):.3f}")

    zs = {"baseline": z0}
    for name, change in VARIANTS.items():
        out = RES_DIR / f"sweep_{name}.csv"
        if out.exists():
            zs[name] = pd.read_csv(out).set_index(["drug", "disease"]).reindex(base.index)["z_score"].values
            print(f"[{name}] already done -> loaded"); continue
        params = {**BASE, **change}
        if FAST_ITER and "n_iter" not in change:
            params["n_iter"] = int(FAST_ITER)
        t0 = time.time()
        print(f"\n[{name}] running with {params} ...")
        res = netmedpy.screening(inp["drug_targets"], inp["disease_genes"], graph, D, score="proximity",
                                 properties=["z_score", "raw_amspl"], n_procs=N_PROCS, **params)
        long = res["z_score"].stack().rename("z_score").reset_index(); long.columns = ["drug", "disease", "z_score"]
        long["raw_distance"] = res["raw_amspl"].stack().values
        long.to_csv(out, index=False)
        zs[name] = long.set_index(["drug", "disease"]).reindex(base.index)["z_score"].values
        print(f"[{name}] done in {time.time() - t0:.0f}s")

    # ---------- compare every variant with the baseline
    # Some z-scores can be undefined (NaN): e.g. with exact degree matching, a target whose degree is unique has no
    # other node to be swapped with, so the random sets have zero spread. Such pairs are dropped for that comparison
    # (the baseline is evaluated on exactly the same pairs) and counted in 'n_undefined'.
    rows = []
    for k, z in zs.items():
        m = np.isfinite(z) & np.isfinite(z0)
        yy, zz, bb = y[m], z[m], z0[m]
        dd_ = disease[m]
        # drug-level bootstrap on the valid pairs
        gr = list(pd.Series(np.arange(m.sum())).groupby(drugs[m]).indices.values())
        bt, bb_ = np.empty(N_BOOT), np.empty(N_BOOT)
        for i in range(N_BOOT):
            idx = np.concatenate([gr[j] for j in rng.integers(0, len(gr), len(gr))])
            bt[i], bb_[i] = auc(yy[idx], -zz[idx]), auc(yy[idx], -bb[idx])
        a = auc(yy, -zz); a0 = auc(yy, -bb)
        lo, hi = np.percentile(bt, [2.5, 97.5]); dl, dh = np.percentile(bt - bb_, [2.5, 97.5])
        within = [spearmanr(zz[dd_ == d], bb[dd_ == d]).statistic for d in np.unique(dd_)]
        cz, cb = zz <= PROXIMAL_Z, bb <= PROXIMAL_Z
        rows.append({"variant": k, "n_pairs_used": int(m.sum()), "n_undefined": int((~m).sum()),
                     "auc": a, "ci_low": lo, "ci_high": hi, "auc_diff_vs_baseline": a - a0,
                     "diff_ci_low": dl, "diff_ci_high": dh,
                     "spearman_vs_baseline": spearmanr(zz, bb).statistic,
                     "spearman_within_disease_median": np.median(within),
                     "pct_proximal": cz.mean(), "call_agreement": np.mean(cz == cb), "kappa": kappa(cz, cb),
                     "mean_z_known": zz[yy == 1].mean(), "mean_z_unknown": zz[yy == 0].mean()})
    summ = pd.DataFrame(rows)
    summ.to_csv(RES_DIR / "parameter_sweep_summary.csv", index=False)
    pd.set_option("display.width", 220, "display.max_columns", 30)
    print("\n" + "=" * 78 + "\nPARAMETER SWEEP (one parameter changed at a time; baseline = log_binning, bin 100, 1000 iterations)\n" + "=" * 78)
    print(summ[["variant", "n_undefined", "auc", "ci_low", "ci_high", "auc_diff_vs_baseline", "diff_ci_low", "diff_ci_high"]].round(3).to_string(index=False))
    print()
    print(summ[["variant", "spearman_vs_baseline", "spearman_within_disease_median", "pct_proximal", "call_agreement", "kappa",
                "mean_z_known", "mean_z_unknown"]].round(3).to_string(index=False))

    # ---------- figure: z of each variant against the baseline
    names = list(VARIANTS)
    fig, ax = plt.subplots(1, len(names), figsize=(4.2 * len(names), 4.2))
    for a_, n in zip(np.atleast_1d(ax), names):
        a_.scatter(z0[y == 0], zs[n][y == 0], s=4, alpha=0.25, color="#8da0cb")
        a_.scatter(z0[y == 1], zs[n][y == 1], s=6, alpha=0.5, color="#d4a017")
        lim = [np.nanmin([np.nanmin(z0), np.nanmin(zs[n])]), np.nanmax([np.nanmax(z0), np.nanmax(zs[n])])]
        a_.plot(lim, lim, "k--", lw=0.7)
        a_.set_xlabel("z, baseline"); a_.set_ylabel(f"z, {n}")
        ok = np.isfinite(zs[n]) & np.isfinite(z0); a_.set_title(f"rho = {spearmanr(z0[ok], zs[n][ok]).statistic:.2f}", fontsize=9)
    plt.tight_layout(); plt.savefig(FIG_DIR / "parameter_sweep.png", dpi=150); plt.close()
    print("\nFigure saved:", FIG_DIR / "parameter_sweep.png\nScript 09 finished.")


if __name__ == "__main__":
    main()