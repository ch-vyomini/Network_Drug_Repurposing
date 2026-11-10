"""
06_benchmark.py
===============
Question : Does network proximity recover APPROVED drug-disease pairs better than simpler
           scores, and better than chance?   ("recovery of known indications")

For every proximity run made by script 05 this compares FIVE scores. Higher score = predicted 'known':
  1. Proximity            : -z   (NetMedPy, degree-matched null)         <- the method
  2. Raw distance         : -d   (same distance, but NO random-set correction)
  3. Overlap fraction     : share of the drug's targets that are themselves disease genes   <- baseline
  4. Overlap hypergeometric: -log10 p of the target/disease-gene overlap (an enrichment test)
  5. Target degree        : mean number of interaction partners of the drug's targets (hubness)
Metric   : AUC = probability that a random known pair scores higher than a random unknown pair
           (0.5 = no better than a coin flip). 95% CIs come from a DRUG-LEVEL bootstrap (all pairs
           of one drug are resampled together, because they are not independent).
Controls : label shuffling (1) across all pairs and (2) within each disease, to see what AUC chance gives;
           AUC restricted to pairs with NO direct target/disease-gene overlap (does the network help
           beyond direct hits?); per-disease and per-group AUCs.
Caveat   : 'unknown' pairs are NOT proven negatives (they may be untested or off-label uses).
Input    : results/proximity_*.csv (script 05), data/processed/inputs_*.pkl, data/networks/*_symbols.pkl
Output   : results/benchmark_summary.csv, benchmark_differences.csv, benchmark_by_disease.csv,
           shuffle_control.csv, figures/benchmark_auc.png
"""

import os
import pickle
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import hypergeom, rankdata

BASE_DIR = Path(os.environ.get("NDR_BASE", (
    r"C:\Users\Vyomini\Documents\MU Engineering_Btech [CB]"
    r"\Projects\Course_Projects\Sem7_Projects"
    r"\Sys Bio & Network Model Project"
    r"\Network_Drug_Repurposing")))
NET_DIR, PROC_DIR, RES_DIR = BASE_DIR / "data/networks", BASE_DIR / "data/processed", BASE_DIR / "results"
FIG_DIR = RES_DIR / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

RUNS = [("conf700", "full"), ("conf900", "shared"), ("conf700", "shared"), ("conf400", "shared")]
PRIMARY = ("conf700", "full")                 # main analysis (all usable data, paper-sized network)
SUFFIX = "log_binning_iter1000"
N_BOOT = int(os.environ.get("NDR_BOOT", 1000))
N_SHUFFLE = int(os.environ.get("NDR_SHUFFLE", 1000))
SEED = 42
GROUPS = {"type 2 diabetes mellitus": "A_cardiometabolic", "heart failure": "A_cardiometabolic",
          "hypercholesterolemia": "A_cardiometabolic", "rheumatoid arthritis": "B_immune",
          "asthma": "B_immune", "psoriasis": "B_immune"}
SCORES = ["Proximity (-z)", "Raw distance (no null)", "Overlap fraction",
          "Overlap hypergeometric", "Target degree"]
rng = np.random.default_rng(SEED)


def auc(y, s):
    """AUC via rank-sum (ties get average ranks)."""
    n1 = int(y.sum()); n0 = len(y) - n1
    if n1 == 0 or n0 == 0:
        return np.nan
    r = rankdata(s)
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def build_scores(df, inp, graph):
    """Add the five score columns to the proximity table."""
    N = graph.number_of_nodes()
    deg = dict(graph.degree())
    frac, hyp, k_over, mdeg = [], [], [], []
    for drug, dis in zip(df["drug"], df["disease"]):
        T, S = inp["drug_targets"][drug], inp["disease_genes"][dis]
        k = len(T & S)
        k_over.append(k)
        frac.append(k / len(T))
        hyp.append(-hypergeom.logsf(k - 1, N, len(S), len(T)) / np.log(10) if k > 0 else 0.0)
        mdeg.append(np.mean([deg[t] for t in T]))
    out = pd.DataFrame({
        "Proximity (-z)": -df["z_score"].values,
        "Raw distance (no null)": -df["raw_distance"].values,
        "Overlap fraction": frac, "Overlap hypergeometric": hyp, "Target degree": mdeg})
    return out, np.array(k_over)


all_summary, all_diff, all_disease, all_shuffle = [], [], [], []
fig_data = {}

for network, mode in RUNS:
    f = RES_DIR / f"proximity_{network}_{mode}_{SUFFIX}.csv"
    if not f.exists():
        print(f"[{network}/{mode}] file missing -> skipped"); continue
    df = pd.read_csv(f)
    with open(PROC_DIR / f"inputs_{network}_{mode}.pkl", "rb") as fh:
        inp = pickle.load(fh)
    with open(NET_DIR / f"physical_{network}_symbols.pkl", "rb") as fh:
        graph = pickle.load(fh)
    sc, k_over = build_scores(df, inp, graph)
    y = df["known"].values.astype(int)
    run = f"{network}/{mode}"
    print("\n" + "=" * 72 + f"\nRUN {run}:  {len(df):,} pairs, {y.sum()} known ({y.mean():.0%})\n" + "=" * 72)

    # ---- direct overlap facts (the paper reported ~15% of known pairs overlap)
    print(f"Direct target/disease-gene overlap:  known pairs {np.mean(k_over[y == 1] > 0):.0%}"
          f"   unknown pairs {np.mean(k_over[y == 0] > 0):.0%}")

    # ---- AUC + drug-level bootstrap (same resamples for all scores => paired differences)
    groups = list(df.groupby("drug").indices.values())
    boot = np.empty((N_BOOT, len(SCORES)))
    for b in range(N_BOOT):
        idx = np.concatenate([groups[i] for i in rng.integers(0, len(groups), len(groups))])
        for j, name in enumerate(SCORES):
            boot[b, j] = auc(y[idx], sc[name].values[idx])
    obs = np.array([auc(y, sc[n].values) for n in SCORES])
    lo, hi = np.nanpercentile(boot, 2.5, axis=0), np.nanpercentile(boot, 97.5, axis=0)
    print(f"\n{'Score':<26}{'AUC':>7}   95% CI")
    for j, n in enumerate(SCORES):
        print(f"{n:<26}{obs[j]:>7.3f}   [{lo[j]:.3f}, {hi[j]:.3f}]")
        all_summary.append({"run": run, "subset": "all pairs", "score": n, "auc": obs[j],
                            "ci_low": lo[j], "ci_high": hi[j], "n_pairs": len(df), "n_known": int(y.sum())})
    print("\nProximity minus each other score (paired, drug-level bootstrap):")
    for j in range(1, len(SCORES)):
        d = boot[:, 0] - boot[:, j]
        dl, dh = np.nanpercentile(d, [2.5, 97.5])
        verdict = "proximity better" if dl > 0 else ("proximity worse" if dh < 0 else "no clear difference")
        print(f"  vs {SCORES[j]:<24}{obs[0] - obs[j]:>+7.3f}   [{dl:+.3f}, {dh:+.3f}]   -> {verdict}")
        all_diff.append({"run": run, "comparison": f"Proximity - {SCORES[j]}", "diff": obs[0] - obs[j],
                         "ci_low": dl, "ci_high": dh, "verdict": verdict})

    # ---- subsets: no direct overlap, per group
    print("\nSubsets (AUC):")
    subsets = {"no direct overlap": k_over == 0,
               "group A cardiometabolic": df["disease"].map(GROUPS).eq("A_cardiometabolic").values,
               "group B immune": df["disease"].map(GROUPS).eq("B_immune").values}
    for sname, mask in subsets.items():
        a = [auc(y[mask], sc[n].values[mask]) for n in SCORES]
        print(f"  {sname:<25} n={mask.sum():>5} known={int(y[mask].sum()):>4}  "
              f"proximity {a[0]:.3f} | raw dist {a[1]:.3f} | overlap frac {a[2]:.3f} | degree {a[4]:.3f}")
        for n, v in zip(SCORES, a):
            all_summary.append({"run": run, "subset": sname, "score": n, "auc": v, "ci_low": np.nan,
                                "ci_high": np.nan, "n_pairs": int(mask.sum()), "n_known": int(y[mask].sum())})

    # ---- per disease
    for dis in sorted(df["disease"].unique()):
        m = (df["disease"] == dis).values
        row = {"run": run, "disease": dis, "n_pairs": int(m.sum()), "n_known": int(y[m].sum())}
        for n in SCORES:
            row[n] = auc(y[m], sc[n].values[m])
        all_disease.append(row)

    # ---- label-shuffle controls for the proximity score
    s_prox = sc["Proximity (-z)"].values
    dis_idx = [np.where(df["disease"].values == d)[0] for d in df["disease"].unique()]
    null_global, null_within = np.empty(N_SHUFFLE), np.empty(N_SHUFFLE)
    for b in range(N_SHUFFLE):
        null_global[b] = auc(rng.permutation(y), s_prox)
        y2 = y.copy()
        for ix in dis_idx:
            y2[ix] = rng.permutation(y[ix])
        null_within[b] = auc(y2, s_prox)
    for label, null in (("shuffled across all pairs", null_global), ("shuffled within each disease", null_within)):
        p = (1 + np.sum(null >= obs[0])) / (N_SHUFFLE + 1)
        print(f"\nChance control ({label}): AUC range {np.percentile(null, 2.5):.3f}-{np.percentile(null, 97.5):.3f}"
              f" (mean {null.mean():.3f});  observed {obs[0]:.3f};  p = {p:.4f}")
        all_shuffle.append({"run": run, "shuffle": label, "null_mean": null.mean(),
                            "null_2.5": np.percentile(null, 2.5), "null_97.5": np.percentile(null, 97.5),
                            "observed": obs[0], "p_value": p})
    fig_data[run] = dict(obs=obs, lo=lo, hi=hi, null=np.percentile(null_global, [2.5, 97.5]))

# ------------------------------------------------------------------
# SAVE TABLES
# ------------------------------------------------------------------
pd.DataFrame(all_summary).to_csv(RES_DIR / "benchmark_summary.csv", index=False)
pd.DataFrame(all_diff).to_csv(RES_DIR / "benchmark_differences.csv", index=False)
pd.DataFrame(all_disease).to_csv(RES_DIR / "benchmark_by_disease.csv", index=False)
pd.DataFrame(all_shuffle).to_csv(RES_DIR / "shuffle_control.csv", index=False)

prim = f"{PRIMARY[0]}/{PRIMARY[1]}"
if prim in fig_data:
    print("\n" + "=" * 72 + f"\nPER-DISEASE AUC, {prim}  (proximity vs overlap fraction)\n" + "=" * 72)
    pdz = pd.DataFrame([r for r in all_disease if r["run"] == prim])
    print(pdz[["disease", "n_pairs", "n_known", "Proximity (-z)", "Raw distance (no null)",
               "Overlap fraction", "Target degree"]].round(3).to_string(index=False))

# ------------------------------------------------------------------
# FIGURE
# ------------------------------------------------------------------
if fig_data:
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
    if prim in fig_data:
        d = fig_data[prim]
        x = np.arange(len(SCORES))
        ax[0].axhspan(d["null"][0], d["null"][1], color="grey", alpha=0.25, label="chance range (shuffled labels)")
        ax[0].bar(x, d["obs"], yerr=[d["obs"] - d["lo"], d["hi"] - d["obs"]], capsize=4,
                  color=["#d4a017"] + ["#8da0cb"] * 4)
        ax[0].axhline(0.5, color="k", lw=0.8, ls="--")
        ax[0].set_xticks(x); ax[0].set_xticklabels([s.replace(" (", "\n(") for s in SCORES], fontsize=8)
        ax[0].set_ylabel("AUC (recovery of approved indications)"); ax[0].set_ylim(0.3, 0.9)
        ax[0].set_title(f"Methods compared, {prim}"); ax[0].legend(loc="upper right", fontsize=8)
    runs = list(fig_data)
    vals = [fig_data[r]["obs"][0] for r in runs]
    err = [[fig_data[r]["obs"][0] - fig_data[r]["lo"][0] for r in runs], [fig_data[r]["hi"][0] - fig_data[r]["obs"][0] for r in runs]]
    ax[1].bar(range(len(runs)), vals, yerr=err, capsize=4, color="#d4a017")
    ax[1].axhline(0.5, color="k", lw=0.8, ls="--")
    ax[1].set_xticks(range(len(runs))); ax[1].set_xticklabels(runs, fontsize=8)
    ax[1].set_ylim(0.3, 0.9); ax[1].set_ylabel("Proximity AUC"); ax[1].set_title("Proximity across networks / modes")
    plt.tight_layout()
    plt.savefig(FIG_DIR / "benchmark_auc.png", dpi=150)
    plt.close()
    print("\nFigure saved:", FIG_DIR / "benchmark_auc.png")
print("\nScript 06 finished.")