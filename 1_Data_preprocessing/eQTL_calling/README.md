# eQTL calling — genotype panel → eight-model scan → filtered eQTL table

> Upstream of [`../SNPs_eQTL/`](../SNPs_eQTL), which classifies and cleans what this produces.
> **Manuscript:** Methods, *"Identification of eQTL"* and *"eQTL identification and classification"*; S1 and S4 Tables.

**Pipeline position:** [⌂ overview](../../README.md) · [`1_Data_preprocessing`](../README.md) · **next ➡** [`../SNPs_eQTL/`](../SNPs_eQTL)

These scripts were written by co-author Jonas Rodriguez in 2021 and ran on the UW-Madison CHTC
HTCondor pool. They are published here as delivered, with three classes of edit only:
cluster usernames and one host address replaced with `<chtc-user>` / `<beast-host>` placeholders,
and clearly-marked `EDITORIAL NOTE` headers on four files. **No analysis logic was changed.**

## Run order

| # | Script | What it does |
|---|---|---|
| 1 | `01.0_genotype_data_download_and_prep.sh` | Download Panzea HapMap3.2.1 + the WiDiv-942 panel (Mazaheri et al. 2019); per source keep biallelic SNPs at MAF ≥ 0.05, set heterozygous calls to missing; `bcftools concat --allow-overlaps --rm-dups all` (WiDiv listed first, so its record wins at shared sites); re-filter the merged panel (`F_MISSING<0.1`, biallelic, MAF ≥ 0.05) → 304-line VCF |
| 2 | `01.1_fix_REF_ALT_var.R` | Fix REF/ALT against the AGPv4 FASTA before the TASSEL hapmap→VCF conversion (called from step 1) |
| 3 | `02.1_make_gct_files_mazaheri.R` | Seedling expression (Mazaheri et al. 2019) → GCT |
| 4 | `02.2_make_gct_files_kremling.R` | Seven Kremling et al. 2018 tissues (GRoot, Gshoot, Kern, L3Base, L3Tip, LMAD, LMAN) → GCT |
| 5 | `02.3_make_individual_beds.R` | Per-tissue PLINK2 genotype subsets (MAF/missingness re-applied per tissue, which is why the SNP count is a range) |
| 6 | `02.0_tensorQTL_data_prep.sh` | GTEx-pipeline expression prep: keep genes with ≥ 6 reads **and** ≥ 0.1 TPM in ≥ 20% of samples, TMM-normalise; then 60 PEER factors per tissue (`broadinstitute/gtex_eqtl:V8`) |
| 7 | `02.4_make_annotation_table.R` | Gene annotation table + syntenic-gene labels (Schnable pan-grass set) |
| 8 | `03.1_format_files_for_rMVP.R` | `MVP.Data(fileKin=TRUE, filePC=TRUE)` → memory-mapped genotypes, kinship matrix, SNP PCs |
| 9 | `05.0_eQTL_permutations.sh` | tensorQTL **cis** permutations (±100 kb, no covariates) → per-gene nominal P thresholds |
| 10 | `05.0_eQTL_permutations_trans.py` | tensorQTL **trans** permutations (`nperms=10000`) → per-tissue threshold = 5th percentile of the minimum-P null ⚠ see header note |
| 11 | `05.1_eQTL_covariates.R` | Assemble the covariate matrix: **5 SNP PCs + top 25 PEER factors**; also writes the zeroed "null" covariate file used by step 10 |
| 12 | `06.0_eQTL_data_split.R` | Split each tissue into 150 marker blocks × 40 gene chunks for the HTCondor array |
| 13 | `06.1_eQTL_analysis_genome_wide_condor.R` | **The eight-model scan** (see below) ⚠ delivered copy is a debug snapshot — see header note |
| 14 | `06.3_identify_jobs_to_rerun.R` | Diff the returned tarballs against the expected 150 × 40 × 8 grid → redo list |
| 15 | `06.2_new_eQTL_analysis_complie_results.R` | **Strict** compile → `00_filtered_cis_trans_all_eQTL_results.txt` (feeds the eGRN) — run before step 16 |
| 16 | `06.2_new_eQTL_analysis_complie_results_liberal.R` | **Liberal** compile, trans only → the GAN input |
| 17 | `06.6_eQTL_results_compute_pve_results_liberal.R` | Join gene coordinates, the synteny label and MAF; add a per-SNP PVE column |
| — | `05.4_..._condor_complie_results.sh` | Reference only: the rsync/scp round trip that pulled results off the submit node |

## The eight candidate linear models

Fitted in rMVP. Models 2–8 all carry the marker-based kinship matrix as a random polygenic effect,
with variance components by Haseman–Elston regression (`vc.method = "HE"`).

| # | Result column | Fixed-effect covariates | Type |
|---|---|---|---|
| 1 | `naive_model_pval` | none (intercept + SNP dosage), `method = "GLM"` | fixed only |
| 2 | `k_model_pval` | none | mixed |
| 3 | `k_q_model_pval` | 5 SNP PCs | mixed |
| 4 | `k_q_p_model1_pval` | 5 SNP PCs + PEER 1–5 | mixed |
| 5 | `k_q_p_model2_pval` | 5 SNP PCs + PEER 1–10 | mixed |
| 6 | `k_q_p_model3_pval` | 5 SNP PCs + PEER 1–15 | mixed |
| 7 | `k_q_p_model4_pval` | 5 SNP PCs + PEER 1–20 | mixed |
| 8 | `k_q_p_model5_pval` | 5 SNP PCs + PEER 1–25 | mixed |

`n_models_support` = how many of the eight cleared the applicable permutation threshold.
**cis and trans are thresholded and filtered differently** — this matters when reading the outputs:

| | Threshold | Retained |
|---|---|---|
| **cis** | per **gene** (tensorQTL cis permutations, ±100 kb) | all 8 models for the eGRN set (`../SNPs_eQTL/Clean_00eQTLs.sh`) |
| **trans** | per **tissue** (5th percentile of the 10,000-permutation minimum-P null) | ≥ 2 models for the GAN |

The archived pilot [`archive/eQTL_calling/05.3_...R`](../../archive/eQTL_calling) states the same eight
models most legibly, but names them differently (`k_q_p_model1..6`, where its `model1` has no PEER term).
Use the column names above.

## Outputs → where they go

| Artifact | Produced by | Consumed by |
|---|---|---|
| `00_filtered_cis_trans_all_eQTL_results.txt` | step 15 | `../SNPs_eQTL/Clean_00eQTLs.sh` → eGRN |
| `..._with_pve_and_synteny_liberal.txt` | step 17 | the trans classification chain → GAN |

## Known gaps

- The **HTCondor submit description, job wrapper and `00_header.txt` were not preserved**, and the
  delivered `06.1` is a debug snapshot (step 13 header note). The models, thresholds and filters are
  fully documented; the exact job-submission harness is not recoverable.
- Two conversion steps between `00_filtered_cis_trans_all_eQTL_results.txt` and the intermediates read
  by `../SNPs_eQTL/Set_CleanFiles_eQTL.R` were one-off bedtools runs whose invocations were not kept —
  see the provenance note in [`../README.md`](../README.md).
- Software versions were not pinned. Known: R 4.0.2 on the execute nodes, PEER via
  `broadinstitute/gtex_eqtl:V8`, Ensembl Plants B73 RefGen_v4 release 50 (annotation builds v4.46,
  v4.48 and v4.50 all appear across the wider pipeline). rMVP was installed from GitHub HEAD in
  April 2021 and tensorQTL from a then-current release; neither commit was recorded.
