# Network-based drug repurposing: does network proximity recover approved indications?

CB 5313 Systems Biology & Network Modeling, Topic 9 (Guney et al., 2016; tool: NetMedPy)

## 1. Question

Guney et al. (2016) propose that a drug is likely to treat a disease when its targets lie close to the disease's genes in the human protein-protein interaction network. "Closeness" ("network proximity") is the average distance from each drug target to its nearest disease gene, expressed as a z-score against degree-matched random protein sets. This project applies **NetMedPy**'s implementation of that algorithm to **independent, current data** (not the package's tutorial data) and asks:

1. Does proximity recover **approved drug-disease pairs** better than chance and better than a simple **direct target-overlap** baseline?
2. How **robust** is the result to the STRING confidence threshold (400 / 700 / 900)?
3. Is a drug **specific** to the disease it is approved for, rather than close to every disease?
4. How sensitive are the results to **NetMedPy's own parameters** (null model, number of random sets, bin size, distance definition)?

## 2. Data (downloaded by scripts; nothing large is committed)

| What | Source | Version | Script |
| ----------------------------- | ------------------------------------------------------------------------------------------------------------------------------------ | ------------- | ---------------------------- |
| Protein interaction network | STRING, *physical* links, human (9606) | v12.5 | `00_download_string.py` |
| Disease genes | Open Targets `association_by_datatype_direct`, evidence type `genetic_association` **only** (no drug evidence, to avoid circularity) | release 26.09 | `02_download_opentargets.py` |
| Drug targets | Open Targets `drug_mechanism_of_action` (+ `target`, `drug_molecule`) | 26.09 | `02` |
| Known indications (benchmark) | Open Targets `clinical_indication`, `maxClinicalStage == APPROVAL` | 26.09 | `02` |

If the STRING download in script 00 fails, download the same files by hand from [STRING downloads](https://string-db.org/cgi/download) (Homo sapiens, v12.5): the *physical links* file and the *protein info* file, and place them in `data/raw/string/`.

## 3. Repository contents

```text
README.md
requirements.txt
.gitignore
scripts/
  00_download_string.py
  01_build_networks.py
  02_download_opentargets.py
  03_coverage_table.py
  04_build_inputs.py
  05_run_proximity.py
  06_benchmark.py
  07_robustness.py
  08_specificity.py
  09_parameter_sweep.py
results/
  selected .csv tables
  figures/*.png
data/                         # local only; excluded from the repository
```

The `data/` directory contains downloaded and processed data and is not committed to GitHub. Scripts 00-05 create the data needed for the analyses.

## 4. How to run

Run these commands from the repository's root directory:

```bash
pip install -r requirements.txt
python scripts/00_download_string.py        # STRING files (or place them in data/raw/string/ by hand)
python scripts/01_build_networks.py         # STRING networks at 400/700/900, largest connected component, shared proteins
python scripts/02_download_opentargets.py   # Open Targets datasets (about 0.8 GB; can take about an hour)
python scripts/03_coverage_table.py         # coverage of candidate diseases -> how the 6 study diseases were chosen
python scripts/04_build_inputs.py           # disease genes, drug targets, known pairs for NetMedPy
python scripts/05_run_proximity.py          # NetMedPy proximity, 4 runs (about 20 min on a laptop)
python scripts/06_benchmark.py              # AUC vs baselines, bootstrap CIs, label-shuffle controls
python scripts/07_robustness.py             # threshold robustness
python scripts/08_specificity.py            # disease specificity
python scripts/09_parameter_sweep.py        # NetMedPy parameters changed one at a time (about 45 min)
```

Outputs go to `results/`. Scripts 05 and 09 print many "... finished" lines: that is NetMedPy's own progress output.

## 5. Methods in brief

- **Disease genes:** genetic-association score >= 0.5, including subtype terms; a disease is kept if it has >= 20 genes in the network. Six diseases: type 2 diabetes, heart failure, hypercholesterolemia (cardiometabolic); rheumatoid arthritis, asthma, psoriasis (immune).
- **Drugs:** approved drugs with >= 1 target in the network; salt/formulation IDs merged into their parent compound.
- **Networks:** STRING physical links with combined score >= 400, 700, 900; largest connected component of each. `full` = everything present in that network; `shared` = only proteins present in all three (fair threshold comparison).
- **Proximity (NetMedPy):** `all_pair_distances` (shortest paths) and `screening` with the "closest" distance (target -> nearest disease gene), `null_model="log_binning"`, `bin_size=100`, `n_iter=1000`; "proximal" means z <= -0.15 (the paper's cut-off).
- **Benchmark:** AUC of known (approved) vs unknown pairs. Five scores compared: proximity, raw distance (no null), overlap fraction, overlap hypergeometric p, mean target degree. Drug-level bootstrap confidence intervals; label-shuffle controls.

## 6. Results (numbers from `results/`)

**Pooled AUC for recovering approved pairs (95% CI).** Shuffled-label chance is about 0.46-0.54; every proximity AUC is above it (p <= 0.001).

| Run | Proximity | Raw distance | Overlap fraction | Target degree | Proximity minus overlap |
| -------------------- | -------------------- | ------------ | -------------------- | ------------- | ----------------------- |
| 700 / full (primary) | 0.572 [0.537, 0.607] | 0.584 | 0.567 [0.543, 0.592] | 0.527 | +0.005 [-0.023, +0.034] |
| 900 / shared | 0.599 [0.565, 0.632] | 0.609 | 0.573 | 0.528 | +0.025 [+0.001, +0.051] |
| 700 / shared | 0.583 [0.544, 0.619] | 0.598 | 0.573 | 0.529 | +0.010 [-0.021, +0.038] |
| 400 / shared | 0.660 [0.626, 0.695] | 0.660 | 0.573 | 0.518 | +0.087 [+0.059, +0.114] |

- Direct overlap between a drug's targets and the disease's genes: 21-23% of known pairs vs 8-9% of unknown pairs.
- AUC on pairs with **no** direct overlap: 0.507 (700/full), 0.540 (900), 0.514 (700/shared), 0.622 (400).
- Per disease (700/full): asthma 0.50, heart failure 0.32, hypercholesterolemia 0.85 (only 17 known drugs), psoriasis 0.65, rheumatoid arthritis 0.76, type 2 diabetes 0.65.
- Per group (700/full): cardiometabolic 0.55, immune 0.61.

**Robustness across thresholds (shared set, 1,368 pairs).**

| | 400 vs 700 | 400 vs 900 | 700 vs 900 |
| ---------------------------------------- | ----------------------- | ----------------------- | ----------------------- |
| Spearman rho of z (all pairs) | 0.62 | 0.58 | 0.80 |
| Top-10 drug overlap (chance about 4%) | 55% | 40% | 62% |
| Top-20 drug overlap (chance about 9%) | 74% | 72% | 82% |
| Agreement of the "proximal" call (kappa) | 67% (0.34) | 64% (0.27) | 79% (0.58) |
| AUC difference | +0.077 [+0.045, +0.112] | +0.062 [+0.032, +0.092] | -0.016 [-0.039, +0.008] |

Restricting to the shared proteins keeps 64-79% of each disease's genes.

**Specificity (within each drug).**

| Run | Drug-centric AUC [95% CI] | p vs shuffled | Mean rank of diagonal (chance 3.5) |
| ------------ | ------------------------- | ------------- | ---------------------------------- |
| 700 / full | 0.517 [0.474, 0.563] | 0.21 | 2.33 |
| 900 / shared | 0.642 [0.607, 0.679] | 0.001 | 1.67 |
| 700 / shared | 0.534 [0.484, 0.577] | 0.06 | 2.00 |
| 400 / shared | 0.632 [0.592, 0.672] | 0.001 | 1.33 |

Drugs proximal to **all 6** diseases / to **none** (observed vs expected if independent): 700/full 40 vs 6 and 26 vs 3; 900/shared 18 vs 5 and 21 vs 2; 700/shared 42 vs 8 and 21 vs 1; 400/shared 15 vs 6 and 3 vs 2.

**NetMedPy parameter sweep (700 / full, one parameter changed at a time; baseline = log_binning, bin_size 100, 1,000 iterations, closest distance T -> S).**

| Variant | AUC [95% CI] | AUC change vs baseline [95% CI] | Spearman rho of z vs baseline | "Proximal" calls: agreement (kappa) | % pairs proximal | Run time |
| ----------------------------------------- | -------------------- | ------------------------------- | ----------------------------- | ----------------------------------- | ---------------- | -------- |
| baseline | 0.572 [0.537, 0.607] | - | 1.000 | - | 53.5% | 208 s |
| null_model = degree_match | 0.574 [0.542, 0.607] | +0.003 [-0.001, +0.006] | 0.993 | 97% (0.94) | 56.1% | 610 s |
| null_model = uniform (no degree matching) | 0.568 [0.536, 0.601] | -0.004 [-0.014, +0.007] | 0.921 | 87% (0.73) | 61.2% | 817 s |
| n_iter = 100 | 0.573 [0.541, 0.604] | +0.001 [-0.002, +0.004] | 0.995 | 98% (0.97) | 53.6% | 98 s |
| symmetric = True (T->S and S->T averaged) | 0.597 [0.565, 0.628] | +0.025 [+0.013, +0.036] | 0.911 | 87% (0.73) | 54.9% | 866 s |
| bin_size = 30 | 0.573 [0.540, 0.607] | +0.002 [-0.000, +0.004] | 0.998 | 99% (0.98) | 53.6% | 402 s |

## 7. Interpretation

Overall, network proximity recovered approved indications better than chance in these analyses, but it did **not consistently outperform direct target overlap**. In the primary 700/full run, the difference was small and its confidence interval included zero. The stronger result at 400/shared suggests that including more lower-confidence interactions may help connect drug targets to disease genes, although it does not prove that the 400 network is universally more accurate.

Performance also differed between diseases. Hypercholesterolemia and rheumatoid arthritis had higher AUCs, while asthma was near chance and heart failure was below chance in the 700/full analysis. These differences may reflect the number and quality of known indications, disease-gene coverage, and the fact that genetic risk genes are not necessarily the same genes targeted by treatments. The hypercholesterolemia result should be interpreted cautiously because it included only 17 known drugs.

Threshold comparisons showed that rankings were moderately to strongly correlated, especially between 700 and 900, but the binary "proximal" calls agreed less consistently. A continuous ranking can remain similar even when pairs close to the z-score cut-off switch between proximal and non-proximal. The large number of drugs close to all or none of the six diseases also suggests that proximity alone does not guarantee disease specificity.

The parameter sweep showed that reducing the number of random sets to 100 or changing the bin size to 30 had very little effect on rankings or AUC in this dataset. The uniform null model changed rankings more, while the symmetric distance definition produced the largest AUC increase. Symmetric distance measures a different relationship by averaging target-to-disease and disease-to-target distances, so its higher AUC should be treated as an alternative analysis rather than proof that it is biologically better.

Taken together, the results support network proximity as a useful complementary measure, particularly for pairs without direct target overlap, but not as a stand-alone test of drug efficacy. Unknown pairs are not confirmed negatives, and the analysis is limited by network coverage, disease-gene selection, and the six diseases studied. These computational results identify patterns in the data; they do not establish clinical effectiveness.

## 8. Limitations

- "Unknown" drug-disease pairs are **not proven negatives** (untested or off-label uses exist); the benchmark measures *recovery of approved indications*.
- 19% of pairs are positives (the paper had about 2%) because diseases with many approved drugs were chosen; only 6 diseases, 3 per group.
- Many drugs share identical targets (e.g. insulin analogues, gliflozins); the drug-level bootstrap treats them as independent, so confidence intervals are probably somewhat too narrow.
- 76 approved drugs (about 23%) could not be used: no target recorded, or targets absent from the network.
- Disease-gene definition (genetic evidence, score >= 0.5, with subtypes) differs from the paper (OMIM + GWAS); results may depend on it.
- STRING's lower-confidence physical edges may include literature-derived evidence, which can favour well-studied proteins; not tested.
- The "proximal" cut-off (z <= -0.15) is lenient: about half of all pairs pass it.
- Proximity gives no direction of effect (helpful vs harmful) and is not clinical evidence.

## 9. Reproducibility and AI-use disclosure

- Versions: STRING v12.5, Open Targets 26.09, NetMedPy 0.2.0.1, Python 3.11.18. Random seed 42 in scripts 06-09.
- **AI use (course policy requires disclosure):** AI assistants (Claude and ChatGPT) were used to discuss project design, help draft and debug analysis scripts, and edit parts of the README, including the interpretation in Section 7. I ran the scripts on my data and reviewed the outputs. The final analysis and conclusions should be checked against the saved results and the course's AI-use policy.
