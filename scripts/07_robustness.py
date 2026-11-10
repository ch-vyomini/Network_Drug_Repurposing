"""
07_robustness.py
================
Question : If the STRING confidence threshold changes (400 / 700 / 900), do the proximity results stay the same?
Data     : the three 'shared' runs from script 05 - SAME 228 drugs x 6 diseases (1,368 pairs) and the same
           proteins in every network, so any difference comes from the network's structure, not from
           proteins appearing or disappearing.
Measures :
  A. Spearman rank correlation of z-scores between thresholds (all pairs, and within each disease)
  B. Overlap of the top-K most proximal drugs per disease (K = 10, 20) vs what chance would give
  C. Agreement of the 'proximal' call (z <= -0.15): % agreement and Cohen's kappa
  D. Does AUC really differ between thresholds?  (paired, drug-level bootstrap on the AUC difference)
  E. Coverage lost when going to the shared protein set (from script 04)
Input    : results/proximity_conf{400,700,900}_shared_*.csv, results/input_coverage.csv
Output   : results/robustness_*.csv, results/figures/robustness.png
"""

import itertools
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

BASE_DIR = Path(os.environ.get("NDR_BASE", (
    r"C:\Users\Vyomini\Documents\MU Engineering_Btech [CB]"
    r"\Projects\Course_Projects\Sem7_Projects"
    r"\Sys Bio & Network Model Project"
    r"\Network_Drug_Repurposing")))
RES_DIR = BASE_DIR / "results"
FIG_DIR = RES_DIR / "figures"; FIG_DIR.mkdir(parents=True, exist_ok=True)
NETS = ["conf400", "conf700", "conf900"]
SUFFIX = "log_binning_iter1000"
PROXIMAL_Z = -0.15
TOP_KS = [10, 20]
N_BOOT = int(os.environ.get("NDR_BOOT", 1000))
rng = np.random.default_rng(42)


def auc(y, s):
    n1 = int(y.sum()); n0 = len(y) - n1
    r = rankdata(s)
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


# ------------------------------------------------------------------
# 1. LOAD AND ALIGN THE THREE RUNS (same drug-disease pairs, same order)
# ------------------------------------------------------------------
runs = {}
for n in NETS:
    d = pd.read_csv(RES_DIR / f"proximity_{n}_shared_{SUFFIX}.csv").set_index(["drug", "disease"])
    runs[n] = d
base = runs[NETS[0]].index
for n in NETS:
    assert runs[n].index.equals(base) or set(runs[n].index) == set(base), f"pairs differ in {n}"
Z = pd.DataFrame({n: runs[n]["z_score"].reindex(base) for n in NETS})
known = runs[NETS[0]]["known"].reindex(base).values.astype(int)
disease = base.get_level_values("disease")
drug = base.get_level_values("drug")
print(f"Aligned pairs: {len(Z):,}  ({known.sum()} known)  drugs: {drug.nunique()}  diseases: {disease.nunique()}")

pairs = list(itertools.combinations(NETS, 2))

# ------------------------------------------------------------------
# A. RANK CORRELATION
# ------------------------------------------------------------------
rows = []
for a, b in pairs:
    rho_all = spearmanr(Z[a], Z[b]).statistic
    per = [spearmanr(Z[a][disease == d], Z[b][disease == d]).statistic for d in sorted(set(disease))]
    rows.append({"pair": f"{a} vs {b}", "spearman_all_pairs": rho_all,
                 "spearman_within_disease_median": np.median(per),
                 "spearman_within_disease_min": np.min(per), "spearman_within_disease_max": np.max(per)})
    print(f"A. {a} vs {b}:  rho (all pairs) = {rho_all:.3f};  within-disease median {np.median(per):.3f} "
          f"(range {np.min(per):.3f} to {np.max(per):.3f})")
pd.DataFrame(rows).to_csv(RES_DIR / "robustness_rank_correlation.csv", index=False)

# ------------------------------------------------------------------
# B. TOP-K OVERLAP PER DISEASE
# ------------------------------------------------------------------
rows = []
n_drugs = drug.nunique()
for K in TOP_KS:
    for a, b in pairs:
        ov = []
        for d in sorted(set(disease)):
            m = disease == d
            ta = set(Z[a][m].sort_values().head(K).index.get_level_values("drug"))
            tb = set(Z[b][m].sort_values().head(K).index.get_level_values("drug"))
            ov.append(len(ta & tb) / K)
            rows.append({"K": K, "pair": f"{a} vs {b}", "disease": d, "overlap_fraction": ov[-1]})
        print(f"B. top-{K} drugs, {a} vs {b}: mean overlap {np.mean(ov):.0%}   (chance would be about {K / n_drugs:.0%})")
pd.DataFrame(rows).to_csv(RES_DIR / "robustness_topk.csv", index=False)

# ------------------------------------------------------------------
# C. AGREEMENT OF THE 'PROXIMAL' CALL
# ------------------------------------------------------------------
rows = []
for a, b in pairs:
    ca, cb = (Z[a] <= PROXIMAL_Z).values, (Z[b] <= PROXIMAL_Z).values
    po = np.mean(ca == cb)
    pe = ca.mean() * cb.mean() + (1 - ca.mean()) * (1 - cb.mean())
    kappa = (po - pe) / (1 - pe)
    rows.append({"pair": f"{a} vs {b}", "agreement": po, "cohens_kappa": kappa,
                 "proximal_a": ca.mean(), "proximal_b": cb.mean()})
    print(f"C. 'proximal' call {a} vs {b}: agreement {po:.0%}, kappa {kappa:.2f}  "
          f"(proximal in {ca.mean():.0%} vs {cb.mean():.0%} of pairs)")
pd.DataFrame(rows).to_csv(RES_DIR / "robustness_calls.csv", index=False)

# ------------------------------------------------------------------
# D. PAIRED BOOTSTRAP OF THE AUC DIFFERENCE BETWEEN THRESHOLDS
# ------------------------------------------------------------------
groups = list(pd.Series(np.arange(len(Z))).groupby(drug).indices.values())
S = {n: -Z[n].values for n in NETS}
boot = {n: np.empty(N_BOOT) for n in NETS}
for i in range(N_BOOT):
    idx = np.concatenate([groups[j] for j in rng.integers(0, len(groups), len(groups))])
    for n in NETS:
        boot[n][i] = auc(known[idx], S[n][idx])
rows = []
for a, b in pairs:
    obs = auc(known, S[a]) - auc(known, S[b])
    lo, hi = np.percentile(boot[a] - boot[b], [2.5, 97.5])
    verdict = "differ" if (lo > 0 or hi < 0) else "no clear difference"
    rows.append({"comparison": f"AUC {a} - AUC {b}", "diff": obs, "ci_low": lo, "ci_high": hi, "verdict": verdict})
    print(f"D. AUC({a}) - AUC({b}) = {obs:+.3f}  95% CI [{lo:+.3f}, {hi:+.3f}]  -> {verdict}")
pd.DataFrame(rows).to_csv(RES_DIR / "robustness_auc_differences.csv", index=False)

# ------------------------------------------------------------------
# E. COVERAGE LOST WHEN RESTRICTING TO THE SHARED PROTEINS
# ------------------------------------------------------------------
cov_file = RES_DIR / "input_coverage.csv"
if cov_file.exists():
    cov = pd.read_csv(cov_file)
    t = cov[cov["mode"] == "full"].pivot(index="disease", columns="network", values="genes_used")
    t.insert(0, "genes_total", cov.drop_duplicates("disease").set_index("disease")["genes_total"])
    t["kept_in_shared_%"] = (100 * t["conf900"] / t["conf400"]).round(0)
    t.to_csv(RES_DIR / "coverage_loss.csv")
    print("\nE. Disease genes available per network (full) - the shared set equals conf900:")
    print(t.to_string())

# ------------------------------------------------------------------
# FIGURE
# ------------------------------------------------------------------
fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
for k, (a, b) in enumerate([("conf400", "conf900"), ("conf700", "conf900")]):
    ax[k].scatter(Z[b][known == 0], Z[a][known == 0], s=6, alpha=0.3, color="#8da0cb", label="unknown pairs")
    ax[k].scatter(Z[b][known == 1], Z[a][known == 1], s=9, alpha=0.6, color="#d4a017", label="known pairs")
    ax[k].set_xlabel(f"z-score, {b}"); ax[k].set_ylabel(f"z-score, {a}")
    ax[k].set_title(f"rho = {spearmanr(Z[a], Z[b]).statistic:.2f}"); ax[k].legend(fontsize=8)
tk = pd.read_csv(RES_DIR / "robustness_topk.csv"); tk = tk[tk["K"] == 10]
for i, p in enumerate(sorted(tk["pair"].unique())):
    v = tk[tk["pair"] == p]["overlap_fraction"]
    ax[2].bar(i, v.mean(), yerr=v.std(), capsize=4, color="#66c2a5")
ax[2].axhline(10 / n_drugs, color="k", ls="--", lw=0.8, label="chance")
ax[2].set_xticks(range(3)); ax[2].set_xticklabels(sorted(tk["pair"].unique()), fontsize=8)
ax[2].set_ylabel("top-10 drug overlap (mean over diseases)"); ax[2].set_ylim(0, 1); ax[2].legend(fontsize=8)
plt.tight_layout(); plt.savefig(FIG_DIR / "robustness.png", dpi=150); plt.close()
print("\nFigure saved:", FIG_DIR / "robustness.png")
print("Script 07 finished.")