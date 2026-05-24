# benchmarking — Path2Space vs. competing ST prediction methods

Reproduces the benchmarking analysis of Shulman et al., *Cell* 2026 — the
Figure 3 panels and Table S3. Compares Path2Space against 16 published
spatial-transcriptomics expression-prediction methods on the Bassiouni Visium
cohort and three external cohorts, plus the legacy 785-HVG benchmarking
dataset.

## Scope: evaluation only, not training

The training and prediction for the 16 competing methods was performed by the
authors of the [**TRIPLEX**](https://www.nature.com/articles/s41467-024-54472-y)
benchmark and is **hosted upstream at Nature Communications**, not here. This
component does not redistribute the competing methods' model weights, training
scripts, or prediction outputs at slide level — only the per-(slide, gene)
Pearson correlations between each method's predictions and the measured
expression on the QC-passing slides, aggregated to per-(cohort, gene)
medians. Path2Space training and prediction live in the companion repo,
[`path2space-companion`](https://github.com/eldadshulman/path2space-companion).

## What's reproduced

| Notebook | Panel | Methods | Genes | Cohorts |
|---|---|---|---|---|
| `01_benchmark_visium_bassiouni.ipynb` | 988-HVG Bassiouni benchmark | 17 (16 competitors + Path2Space) × 2 norms | 988 HVGs | Bassiouni (cross-val), HEST + HTAN + Martinez (external) |
| `01_benchmark_visium_bassiouni.ipynb` | Genome-wide top-4 | 4 top competitors + Path2Space (smoothed/unsmoothed) | ~13,000 | same four cohorts |
| `02_benchmark_legacy_dataset.ipynb`   | Legacy 785-HVG benchmark   | Two **disjoint-by-design** panels: 7 raw-count methods + Path2Space (`no_library_size`), 4 CPM methods + Path2Space log(CPM) (`library_size`). Path2Space is the only method evaluated under both normalizations. | 785 HVGs | the legacy benchmark's cross-validation splits |

Statistical comparison: per-(cohort, validation type) Mann-Whitney U test of
the matched per-gene PCC vectors of Path2Space vs. each competitor.

## Contents

- `lib/benchmark.py` — `model_order_by_median`, `mannwhitney_vs_reference`
  (per-(group) Mann-Whitney U of every model against a reference, on matched
  per-gene PCC vectors), `median_pcc_table`.
- `notebooks/01_benchmark_visium_bassiouni.ipynb` — 988-HVG panel (both
  normalizations) plus the genome-wide top-4 panel.
- `notebooks/02_benchmark_legacy_dataset.ipynb` — legacy 785-HVG benchmark
  (two normalizations with disjoint method sets).
- `data/bench_988_bassiouni.parquet` — bundled input for notebook 01's 988-HVG
  panel (~0.9 MB, aggregated to one PCC per (norm, model, validation_type,
  cohort, gene) by averaging per-patient PCCs of the QC-passing slides).
- `data/bench_genome_wide_top4.parquet` — bundled input for the genome-wide
  panel (~2.3 MB).
- `data/bench_785_legacy.parquet` — bundled input for the legacy benchmark
  (~65 KB).
- `data/hvg_988_bassiouni.csv` — the 988 HVG list (derived from the bundled
  data).
- `data/hvg_785_legacy.csv` — the 785 HVG list (from the legacy benchmark).
- `data/median_pcc_*.pkl`, `data/mannwhitney_*.pkl` — summary outputs the
  notebooks write.

## Method

The per-method per-(slide, gene) PCC tables (~7.5 MB each, 33 method
variants across two normalizations) live upstream and are not bundled. The
aggregation collapses them to one PCC per (norm, model, validation_type,
cohort, gene) by:

1. **QC filtering** — keep only the 40 QC-passing slides
   (`qc/data/metadata/st_samples_retained.csv`).
2. **Mean per patient** — within each (norm, model, validation_type, cohort,
   patient, gene), take the mean PCC across that patient's slides.
3. **Mean per cohort** — within each (norm, model, validation_type, cohort,
   gene), take the mean across patients.

The bundled parquet is the output of step 3.

## HVG selection (upstream)

For the 988-HVG benchmark, HVGs were selected per the protocol in the
canonical R notebook `joo_benchmarking/scr/HVG_SVG_R.ipynb` (R kernel,
scran + Seurat):

1. Per slide, top 1,000 HVGs via `scran::modelGeneVar` + `getTopHVGs`.
2. Union across QC-passing slides of the Bassiouni cohort.
3. Drop genes expressed in < 70 % of spots of those slides.

This yields the canonical 988-gene list used by the 17 competing methods'
prediction tables (the intersection of genes present across all methods).
The 785 genes of the legacy benchmark were defined upstream by the original
Nature Communications study and are reused as-is.

The Path2Space prediction tables contain **one additional gene**,
`AL627309.1`, that is not present in any of the 17 competitors' outputs.
Because cross-method comparison requires the same gene set, this gene is
filtered out of `bench_988_bassiouni.parquet`; the parquet and the bundled
`hvg_988_bassiouni.csv` both contain exactly the 988 canonical HVGs.

## HVG-selection (R) — required packages

The `HVG_SVG_R.ipynb` notebook itself is not run by this component, but is
referenced for the selection method. If you want to re-run it, it needs an
R 4.3+ environment with `scran`, `Seurat`, and `reticulate` for the Python ↔ R
bridge.

## Reproduction notes

- **Legacy 785-HVG normalization split is by design, not by reporting gap.**
  The original Nature Communications benchmark applied each method using the
  normalization its own authors specified — methods designed to consume
  library-size-normalized counts were evaluated on the `library_size` panel;
  methods designed to consume raw counts on the `no_library_size` panel. The
  two method sets are therefore **disjoint by design**, not overlapping with
  one subset reporting less. There is no method present in both panels except
  Path2Space, which we trained under both normalizations to demonstrate
  robustness.
- **SGN is excluded** from all comparisons in notebook 01, matching the
  published Figure 3 exclusion. SGN's `no_cpm` median is anomalously high
  (0.42 vs Path2Space 0.39 in external validation) while its `cpm` median
  collapses to ≈ 0 — a pattern not seen for any other method and consistent
  with a normalization or prediction artifact in SGN's `no_cpm` outputs.
  The raw SGN rows remain in `data/bench_988_bassiouni.parquet` for any
  reader who wants to inspect them; the notebook filters SGN out at load
  time.
- **Path2Space variants.** The Path2Space smoothed prediction
  (`Path2Space_log_converted` upstream → `Path2Space` after renaming) and
  the unsmoothed variant (`Path2Space_unsmoothed_log_converted` →
  `Path2Space (unsmoothed)`) are both bundled. The notebooks use the
  smoothed variant as the Mann-Whitney reference; the unsmoothed variant is
  shown alongside the competitors for context but excluded from the
  Mann-Whitney comparisons.
- **Cohort-name aliases.** The upstream pkls use `TNBC` and
  `pierre_martinez` for two cohorts; the bundled parquet renames these to
  `Bassiouni et al.` and `Martinez et al.` to match the figure labels.

## Running

```bash
pip install -r requirements.txt
jupyter nbconvert --to notebook --execute --inplace notebooks/01_benchmark_visium_bassiouni.ipynb
jupyter nbconvert --to notebook --execute --inplace notebooks/02_benchmark_legacy_dataset.ipynb
```

## Citation

Shulman et al., *Cell* 2026. doi:10.1016/j.cell.2026.04.023
