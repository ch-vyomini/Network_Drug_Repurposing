"""
08_specificity.py
=================
Question : Is a drug closer to the disease it is APPROVED for than to the other diseases?
           (the paper's idea: gliclazide is close to type 2 diabetes but not to leukaemia)
Why this is the cleanest test: inside ONE drug the targets are fixed, so hubness and other drug
           properties cancel out - only the disease changes.
Analyses (for each run from script 05):
  1. DRUG-CENTRIC AUC : for each drug with >=1 known and >=1 unknown disease, AUC of 'known diseases have
                        lower z than unknown diseases', averaged over drugs. Chance = 0.5, obtained by shuffling the
                        known labels among the 6 diseases WITHIN each drug. 95% CI: bootstrap over drugs.
  2. MATRIX           : mean z of the drugs approved for disease X (rows) against each disease Y's genes (columns),
                        centred per column (removes 'this disease's gene set is simply central'). If proximity is
                        specific, the DIAGONAL should be the most negative in each row.
  3. PROMISCUITY      : in how many of the 6 diseases is each drug 'proximal' (z <= -0.15)?  Compared with what
                        independence would give; many all-or-none drugs = drug-level effects, not disease-specific.
Input    : results/proximity_<network>_<mode>_log_binning_iter1000.csv
Output   : results/specificity_summary.csv, specificity_drug_centric_<run>.csv, specificity_matrix_<run>.csv,
           specificity_promiscuity.csv, figures/specificity.png
"""

import os
from math import comb
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import rankdata

BASE_DIR = Path(os.environ.get("NDR_BASE", (
    r"C:\Users\Vyomini\Documents\MU Engineering_Btech [CB]"
    r"\Projects\Course_Projects\Sem7_Projects"
    r"\Sys Bio & Network Model Project"
    r"\Network_Drug_Repurposing")))
RES_DIR = BASE_DIR / "results"
FIG_DIR = RES_DIR / "figures"; FIG_DIR.mkdir(parents=True, exist_ok=True)
RUNS = [("conf700", "full"), ("conf900", "shared"), ("conf700", "shared"), ("conf400", "shared")]
PRIMARY = "conf700/full"
SUFFIX = "log_binning_iter1000"
PROXIMAL_Z = -0.15
N_PERM = int(os.environ.get("NDR_PERM", 1000))
N_BOOT = int(os.environ.get("NDR_BOOT", 1000))
rng = np.random.default_rng(42)


def drug_auc(scores, labels):
    """Per-row AUC. scores/labels: drugs x diseases; high score = predicted known. Rows need 1..n-1 known."""
    n = scores.shape[1]
    k = labels.sum(axis=1)
    R = rankdata(scores, axis=1)
    return ((R * labels).sum(axis=1) - k * (k + 1) / 2) / (k * (n - k))


summary, promisc_rows, mats = [], [], {}
for network, mode in RUNS:
    f = RES_DIR / f"proximity_{network}_{mode}_{SUFFIX}.csv"
    if not f.exists():
        print(f"[{network}/{mode}] missing -> skipped"); continue
    run = f"{network}/{mode}"
    df = pd.read_csv(f)
    Z = df.pivot(index="drug", columns="disease", values="z_score")
    K = df.pivot(index="drug", columns="disease", values="known").reindex(index=Z.index, columns=Z.columns).astype(int)
    names = df.drop_duplicates("drug").set_index("drug")["drug_name"].reindex(Z.index)
    print("\n" + "=" * 72 + f"\nRUN {run}:  {Z.shape[0]} drugs x {Z.shape[1]} diseases\n" + "=" * 72)

    # ---------- 1. drug-centric AUC
    keep = (K.sum(axis=1) >= 1) & (K.sum(axis=1) <= Z.shape[1] - 1)
    S, L = -Z[keep].values, K[keep].values
    a = drug_auc(S, L)
    obs = a.mean()
    boot = np.array([a[rng.integers(0, len(a), len(a))].mean() for _ in range(N_BOOT)])
    k_row = L.sum(axis=1); n_dis = L.shape[1]
    null = np.empty(N_PERM)
    for i in range(N_PERM):                         # shuffle known-labels among diseases WITHIN each drug
        Lp = (rankdata(rng.random(L.shape), axis=1) <= k_row[:, None]).astype(int)
        null[i] = drug_auc(S, Lp).mean()
    p = (1 + np.sum(null >= obs)) / (N_PERM + 1)
    lo, hi = np.percentile(boot, [2.5, 97.5])
    print(f"1. Drug-centric AUC (n={keep.sum()} drugs): {obs:.3f}  95% CI [{lo:.3f}, {hi:.3f}]"
          f"   chance range {np.percentile(null, 2.5):.3f}-{np.percentile(null, 97.5):.3f}   p = {p:.4f}")
    print(f"   share of drugs whose known disease is the closest or tied: {(a > 0.999).mean():.0%};"
          f"  share with AUC < 0.5 (known disease FARTHER than others): {(a < 0.5).mean():.0%}")
    summary.append({"run": run, "n_drugs": int(keep.sum()), "drug_centric_auc": obs, "ci_low": lo, "ci_high": hi,
                    "chance_2.5": np.percentile(null, 2.5), "chance_97.5": np.percentile(null, 97.5), "p_value": p})
    pd.DataFrame({"drug": Z[keep].index, "drug_name": names[keep].values, "n_known_diseases": k_row,
                  "auc": a}).sort_values("auc").to_csv(RES_DIR / f"specificity_drug_centric_{run.replace('/', '_')}.csv", index=False)

    # ---------- 2. matrix: drugs approved for X (rows) vs disease Y (columns)
    colmean = Z.mean(axis=0)
    M = pd.DataFrame({y: [Z.loc[K[x] == 1, y].mean() for x in Z.columns] for y in Z.columns}, index=Z.columns)
    Mc = M - colmean
    diag_rank = [int(rankdata(Mc.loc[x].values)[list(Mc.columns).index(x)]) for x in Mc.index]
    print("\n2. Mean z of drugs approved for the ROW disease, against each COLUMN disease (centred per column;"
          "\n   more negative = closer). Diagonal should be the lowest in its row:")
    print(Mc.round(2).to_string())
    print("   rank of the diagonal within its row (1 = best, 6 = worst):",
          dict(zip(Mc.index, diag_rank)), f"-> mean {np.mean(diag_rank):.2f} (chance 3.5)")
    Mc.to_csv(RES_DIR / f"specificity_matrix_{run.replace('/', '_')}.csv")
    mats[run] = Mc

    # ---------- 3. promiscuity
    n_close = (Z <= PROXIMAL_Z).sum(axis=1)
    pr = (Z <= PROXIMAL_Z).values.mean(); nd = Z.shape[1]
    obs_c = n_close.value_counts().reindex(range(nd + 1), fill_value=0)
    exp_c = [len(Z) * comb(nd, j) * pr ** j * (1 - pr) ** (nd - j) for j in range(nd + 1)]
    print(f"\n3. Number of diseases each drug is 'proximal' to (z <= {PROXIMAL_Z}); overall proximal rate {pr:.0%}")
    print("   diseases:", list(range(nd + 1)))
    print("   observed:", obs_c.tolist())
    print("   expected if independent:", [round(e) for e in exp_c])
    print(f"   drugs proximal to ALL {nd}: {(n_close == nd).sum()}  |  to NONE: {(n_close == 0).sum()}")
    for j in range(nd + 1):
        promisc_rows.append({"run": run, "n_diseases_proximal": j, "observed": int(obs_c[j]), "expected_if_independent": exp_c[j]})

pd.DataFrame(summary).to_csv(RES_DIR / "specificity_summary.csv", index=False)
pd.DataFrame(promisc_rows).to_csv(RES_DIR / "specificity_promiscuity.csv", index=False)

# ---------- figure: matrix (primary), matrix (400), promiscuity (primary)
if PRIMARY in mats:
    show = [PRIMARY] + (["conf400/shared"] if "conf400/shared" in mats else [])
    fig, ax = plt.subplots(1, len(show) + 1, figsize=(5.4 * (len(show) + 1), 4.8))
    for i, r in enumerate(show):
        m = mats[r]; v = np.abs(m.values).max()
        im = ax[i].imshow(m.values, cmap="RdBu", vmin=-v, vmax=v)
        ax[i].set_xticks(range(len(m))); ax[i].set_xticklabels([c[:14] for c in m.columns], rotation=60, ha="right", fontsize=7)
        ax[i].set_yticks(range(len(m))); ax[i].set_yticklabels([c[:18] for c in m.index], fontsize=7)
        ax[i].set_title(f"Centred mean z, {r}\n(rows: drugs approved for; cols: disease genes)", fontsize=8)
        ax[i].set_xlabel("disease whose genes are tested", fontsize=8)
        plt.colorbar(im, ax=ax[i], fraction=0.046)
    pp = pd.DataFrame(promisc_rows); pp = pp[pp["run"] == PRIMARY]
    w = 0.4; xs = np.arange(len(pp))
    ax[-1].bar(xs - w / 2, pp["observed"], w, label="observed", color="#d4a017")
    ax[-1].bar(xs + w / 2, pp["expected_if_independent"], w, label="expected if independent", color="#8da0cb")
    ax[-1].set_xlabel("number of diseases the drug is proximal to"); ax[-1].set_ylabel("number of drugs")
    ax[-1].set_title(f"Drug promiscuity, {PRIMARY}", fontsize=9); ax[-1].legend(fontsize=8)
    plt.tight_layout(); plt.savefig(FIG_DIR / "specificity.png", dpi=150); plt.close()
    print("\nFigure saved:", FIG_DIR / "specificity.png")
print("\nSummary across runs:")
print(pd.DataFrame(summary).round(3).to_string(index=False))
print("\nScript 08 finished.")