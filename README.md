# isotope_fw_metric_analysis

[日本語版はこちら](README_ja.md)

Global analysis of how latitude, temperature and temperature seasonality shape the structure of freshwater food webs, using stable-isotope similarity networks built from the IsoFresh database (606 food webs with at least 10 taxa).

## Summary of findings

- Toward the poles, connectance, nestedness and link strength decline (standardised effects -0.2 to -0.3), consistently in lakes and rivers and after controlling for sampling effort.
- About half of the decline in connectance is a direct climate effect; the rest is mediated by taxon richness and prey niche area (piecewise SEM). Nestedness and link strength act mainly through richness.
- Connectance rises by ~0.004 per deg C, matching the sign and magnitude predicted by Petchey, Brose & Rall (2010) for freshwater fish activation energies (Rall et al. 2012).
- Modularity (Q) increases toward the poles; standardised modularity (z-Q) shows no net trend because a richness-mediated decline and a direct climate increase cancel.

See `reports/report_concise.md` (Japanese) for the full write-up.

## Repository layout

| Path | Contents |
|---|---|
| `R/` | Analysis scripts, numbered in run order |
| `results/` | Model outputs (CSV / text). `server_rerun/` reproduces the sampling-effort analyses on the server (identical results) |
| `results/theory_results/` | Theory-driven tests (Petchey 2010, niche packing, fish generalism) |
| `figures/` | Manuscript and supplementary figures (PNG + PDF) |
| `reports/` | Report drafts (Markdown and Word) |
| `data/` | Place IsoFresh CSVs here (not tracked) |

## Scripts

| Script | Purpose |
|---|---|
| `01_network_analysis_full_pipeline.R` | Builds isotope-similarity networks, computes metrics and null-model z-scores, extracts WorldClim/HydroLAKES covariates, runs the mixed-model battery, bootstrap mediation and piecewise SEMs. Writes `IsoFresh_Metrics_v2.RData` and `figures_v2/` |
| `02_make_manuscript_figures.R` | Main and supplementary figures from the pipeline outputs |
| `03_sampling_effort_sensitivity.R` | Climate effects with log10 summed `N_Sample` (sampling effort) as covariate; subset of well-sampled webs |
| `04_sem_sampling_effort.R` | piecewise SEM with sampling effort in every equation (latitude, temperature, seasonality) |
| `05_mediation_bootstrap_sampling_effort.R` | Study-level cluster bootstrap (999) of indirect effects, sampling effort included |
| `06_sem_modularity_zQ.R` | SEM decomposition for Q and z-Q |
| `07_sem_path_diagrams.R` | SEM path diagrams (Fig8, FigS10, FigS11) |
| `08_theory_tests.R` | Petchey et al. 2010 test, temperature x body-size, link-species scaling x latitude, niche packing/expansion, fish trophic height |
| `09_theory_sem.R` | SEM including fish trophic height as a mediator |
| `10_theory_figures.R` | Figures for the theory tests (Fig5-Fig7) |

## How to run

```bash
mkdir work && cp data/IsoFresh_*.csv work/ && cd work
Rscript ../R/01_network_analysis_full_pipeline.R          # heavy step (network metrics + null models)
for s in 02 03 04 05 06 07 08 09 10; do ISOFRESH_RAW=. Rscript ../R/${s}_*.R; done
```

Useful environment variables: `ISOFRESH_RAW` (folder with the IsoFresh CSVs), `ISOFRESH_DATA` (folder with `IsoFresh_Metrics_v2.RData`), `ISOFRESH_N_SIMUL` (null-model replicates, default 199), `ISOFRESH_N_BOOT` (bootstrap replicates, default 999), `ISOFRESH_CORES`.

## Environment

Main analyses: R 4.5.1 on Ubuntu 24.04 (lme4 2.0.6, lmerTest 3.2.1, piecewiseSEM 2.3.1, rnetcarto 0.2.6, igraph 2.2.1, vegan 2.7.2, terra 1.9.46, rfishbase 5.0.1). Theory tests were run with R 4.4.2 on macOS.

## Methods in brief

- Nodes: taxa with mean d13C and d15N per food web. Edge weight: 1 - (Euclidean distance / maximum distance within the web). Links: pairs above one global similarity threshold (0.66, mean of all pooled similarities).
- Metrics: connectance L / [S(S-1)/2]; weighted NODF and modularity (rnetcarto) standardised against a symmetric weight-permutation null model; mean link strength; mean closeness.
- Climate: |latitude|; WorldClim 2.1 BIO1 (mean annual temperature) and BIO4/100 (temperature seasonality, SD of monthly means, deg C).
- Statistics: linear mixed models with `(1 | Region/Source_DOI)`, BH-FDR; piecewise SEM; study-level cluster bootstrap.
