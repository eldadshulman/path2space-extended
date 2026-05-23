# tissue_preservation — Figure 2G

Reproduces the fresh-frozen vs FFPE generalization analysis of Shulman et al.,
*Cell* 2026.

Path2Space is trained only on fresh-frozen (FF) breast-cancer slides — yet most
clinical archives are FFPE (formalin-fixed paraffin-embedded). This component
tests whether predictions hold up on FFPE tissue, and asks which slide-level
properties account for whatever preservation effect remains.

## The two analyses

| Analysis | Question | Notebook section |
|---|---|---|
| **Figure 2G** | Per-section gene-wise PCC stratified by FF vs FFPE — do the distributions overlap? Does mean PCC differ after adjusting for cohort? | Violin plot + cohort-adjusted MixedLM |
| **Nuclear morphology** | Of eight CellPose-derived nuclear features (area, equivalent diameter, circularity, solidity, eccentricity, aspect ratio, nearest-neighbor distance, local nuclear density), which best explain prediction accuracy after preservation is held fixed? | Per-feature MixedLM + BH-FDR |

## Contents

- `lib/ff_ffpe_stats.py` — `ff_vs_ffpe_test` (cohort-adjusted MixedLM), `morphology_associations` (per-feature MixedLM with BH-FDR), `MORPHOLOGY_FEATURES`.
- `notebooks/01_ff_vs_ffpe.ipynb` — Figure 2G violin, cohort-adjusted test, eight-feature morphology MixedLM, eccentricity-vs-PCC scatter.
- `data/per_sample_gene_pcc.parquet` — per (sample, gene) PCC for the two external cohorts with both FF and FFPE sections (HEST + HTAN), restricted to genes present in every sample of the cohort.
- `data/per_sample_summary.parquet` — per-sample mean PCC, eight CellPose morphology features (sample-level means), cohort, preservation method.

## Method

Per-sample mean PCC and the eight nuclear-morphology features (sample-level
means over CellPose-segmented nuclei) are the bundled inputs. Stage A — the
per-spot CellPose segmentation and per-spot morphology computation — needs the
full H&E images and is not bundled. Per-sample aggregates are.

The cohort-adjusted FF-vs-FFPE test fits
`mean_pcc ~ ff + (1 | cohort)` on the per-sample table restricted to the
external (non-training) cohorts. The morphology test fits
`mean_pcc ~ feature + preservation + (1 | cohort)` separately for each of the
eight features, with Benjamini-Hochberg FDR across the eight tests.

## Bundled data

| File | Rows | Notes |
|---|---|---|
| `per_sample_gene_pcc.parquet` | 163,865 | per (sample, gene) PCC for HEST (5 samples × 11,614 common genes) + HTAN (9 samples × 11,755 common genes) |
| `per_sample_summary.parquet` | 40 | per-sample mean PCC + eight CellPose morphology features + cohort + preservation method |

Upstream inputs that are not bundled (too large or not redistributable):

- The per-sample H&E tiles and the per-spot CellPose segmentation outputs that produce the morphology features.
- The full Path2Space per-spot predicted and measured expression matrices used to compute the per-(sample, gene) PCCs.

## Reproduction notes

- **Scope of the morphology MixedLM.** This component fits the morphology test at the per-sample level (n = 40) on bundled aggregated features, which makes the eccentricity coefficient (β ≈ −1.4) and its sign reproducible from the bundled data, but the per-feature p-values are noisier than at a per-spot scale. The eccentricity association is the strongest negative effect across the eight features in this per-sample fit; the paper's reported `q = 6.5 × 10⁻⁷` corresponds to an analysis with substantially more degrees of freedom than the per-sample table supports.
- **External cohort definition.** "External validation" here means the three non-training cohorts: HEST, HTAN, pierre_martinez. Only HEST and HTAN contribute both FF and FFPE sections; pierre_martinez is FF only and is included in the MixedLM fit but contributes no FF/FFPE contrast.
- **statsmodels version required.** `requirements.txt` pins `statsmodels>=0.14.6`. Earlier versions (including 0.14.4) return `NaN` p-values for the eccentricity MixedLM fit — the powell optimizer settles on a non-positive-definite Hessian; `lbfgs` converges to the same MLE with usable p-values, and 0.14.6+ uses it consistently. If you see the eccentricity row with `NaN` p-value, check your statsmodels version first.

## Running

```bash
pip install -r requirements.txt
jupyter nbconvert --to notebook --execute --inplace notebooks/01_ff_vs_ffpe.ipynb
```

## Citation

Shulman et al., *Cell* 2026. doi:10.1016/j.cell.2026.04.023
