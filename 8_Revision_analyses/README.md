# 8 · Revision analyses (PLOS Computational Biology resubmission)

> Analyses added in revision, answering the reviewer points on threshold robustness, class-frequency confounding, and annotation circularity. Their headline outputs are the **per-association stability column in S8 Table** and the S16 sensitivity figure; full result tables ship with the Zenodo data deposit (**DOI [10.5281/zenodo.21866340](https://doi.org/10.5281/zenodo.21866340)**).

**Pipeline position:** cross-cutting — consumes outputs of sections 2–3. [⌂ overview](../README.md)

All scripts read their inputs from a working directory passed via the `SENS_DIR` environment variable (see each header). Inputs are the pipeline artifacts named in the [overview manifest](../README.md#data-artifacts--zenodo-manifest) plus, for exact replication of the published baseline, the archived working files listed in each script.

## `threshold_sensitivity/` — MRMI-cutoff sensitivity (answers R1-24 / R2-14)

Replicates the network-based annotation (per-TF MRMI neighborhoods at decay D ≥ 0.005) and re-runs it at D ≥ {0.005, 0.02, 0.05, 0.10, 0.20} (MRMI ≤ 266–81).

| # | Script | What it does |
|---|---|---|
| 0 | `0_install_topgo.R` | one-time: install topGO into a scratch library |
| 1 | `1_validate_neighborhoods.R` | prove the harness reproduces every published neighborhood size (2,915 TFs, exact) |
| 2 | `2_validate_pvals.R` | prove the PWY p-value replication (GeneOverlap ≡ phyper, < 1e-12) |
| 3 | `3_sweep_PWY_sensitivity.R` | CornCyc PWY enrichment at all five thresholds |
| 4 | `4_sweep_GO_sensitivity.R` | topGO BP enrichment (paper's verbatim calls) at all five thresholds |
| 5 | `5_export_metrics.R` | yield / retention / Jaccard metrics + stability scores per arm |
| 6 | `6_compare_thresholds.R` | stratified retention, dropout decomposition, robust core |
| 7 | `7_process_GO_results.R` | GO metrics + S8 parent-level stability via the recorded term→parent map |

Key results: annotation yield is stable across the 40-fold threshold range (71→92 TFs at FDR ≤ 0.05 while median neighborhood shrinks 227→54); turnover concentrates in marginal calls; ~50% of catalog associations are significant in ≥3 of 5 thresholds in both arms (PWY 50.2%; GO 52.0% of the 5,469 associations assessable under the current GO release — terms obsolete in today's ontology are flagged NA, not scored).

## `gan_snp_normalization/` — class-frequency normalization (answers R2-07)

`1_snp_freq_normalization.R` rebuilds the GAN (155,058 edges; 23,943 sources; 18,958 targets — matching the paper exactly) and normalizes gene-class contributions by class SNP counts. Result: per-SNP interaction rates are flat across classes (0.106–0.148); enzyme/kinase dominance in raw counts reflects class SNP frequency, not preferential connectivity.

## `paralog_null/` — shared-association null model (answers R2-17)

`1_paralog_null_control.R` replicates the 932 published paralog-pair Jaccard values exactly, then tests each similarity bin against neighborhood-size-matched random TF pairs (20 per pair, fixed seed) plus a 20,000-pair global null. Result: only bin IX exceeds the null (3.55×, empirical P ≈ 0.0005); bins II–VIII are indistinguishable from random and bin I falls below it.

## `fig4A_tiers/` — Fig 4A grouping diagnostics and terciles (answers R2-13 / R2-14)

| # | Script | What it does |
|---|---|---|
| 1 | `1_K_diagnostics_rebuild.R` | rebuild the TF × dataset condition-enrichment matrix from the GSEA results and run elbow / silhouette / gap diagnostics |
| 2 | `2_K_diagnostics_original.R` | the same diagnostics on the archived original matrix + reproducibility check of the published k-means memberships (71.5% agreement, unseeded) |
| 3 | `3_tercile_assignments.R` | the deterministic tercile tiers that replace the k-means clusters (cutoffs 45.5 / 60.9%; 277/276/276 TFs) |

## Environment

Run with R 4.6.1: data.table 1.18.2.1, cluster 2.1.8.2, topGO 2.64.0, GO.db 3.23.1 (Bioconductor), ggplot2. Note that GO-arm significance calls depend on the installed GO.db release; the original analysis used the 2023 release (see S3 Text of the manuscript).
