# Multi-Omic Network Integration in Maize

Code and processed-data description for **"Prioritizing Maize Metabolic Gene Regulators through Multi-Omic Network Integration"** (Gomez-Cano et al.).

This repository documents the primary data, processed data, and code used to build four transcription-factor (TF)–target network layers from maize expression, protein–DNA interaction, and eQTL data; integrate them; and predict, prioritize, and evaluate TF functions.

- **Paper:** _[journal / DOI — pending]_
- **Processed outputs** (network edge-lists, gene embeddings, the TF→function catalog, prioritization scores): **Zenodo _[DOI — pending]_**
- **License:** MIT (see `LICENSE`)

## Repository structure (pipeline order)

| Section | Contents |
|---------|----------|
| `1_Data_preprocessing/` | Raw & processed data prep by modality — `PDI/` (ChIP-/DAP-seq download, trimming, mapping, peak calling, peak→target assignment), `SNPs_eQTL/` (SNP→bed, eQTL identification/annotation), `Expression_coexpression/` (expression matrices, weighted-PCC) |
| `2_Network_construction/` | The four TF–target layers — `RFN_expression/` (RF-inferred regulatory network), `GAN_trans_eQTL/`, `GRN_PDI/`, `eGRN_cis_eQTL/` |
| `3_Functional_annotation/` | TF-annotation strategies — `common_target/`, `common_function/`, `network_based_embedding/` (PecanPy embeddings + MI/MR distance + clustering), `annotation_PWY_GO/` |
| `4_Evaluation_knockouts/` | Benchmarking predictions against TF-knockout DEGs |
| `5_Evaluation_random_networks/` | Benchmarking against randomized networks |
| `6_Prioritization_and_conditions/` | rZ/URS prioritization of regulators; GSEA condition-mapping |
| `7_Paralog_redundancy/` | Embedding-based prediction of paralog redundancy/divergence |
| `figures/` | Scripts that render the manuscript figures |
| `environment/` | Software/package versions to reproduce the analyses |
| `archive/` | Superseded script versions, kept for provenance |

> **Terminology note:** the expression-based layer built with random-forest regression (following Zhou et al. 2020) is referred to here as the **RF-inferred regulatory network (RFN)** — previously labeled "co-expression network (CEN)". It predicts regulatory (TF→target) relationships from expression and does not, by itself, imply correlation-based co-expression or causal regulation.

## Requirements

- **R** — topGO, GeneOverlap, GOSemSim, Rrvgo, DECIPHER, Parmigene, FGSEA, wCorr, igraph, Rsubread, DESeq2 (see `environment/`).
- **Python** — PecanPy, numpy, pandas, scikit-learn.
- **Command-line** — bcftools, Bowtie2, SAMtools, Trimmomatic, GEM, FastQC.

Exact versions are pinned in `environment/`.

## How to run

Sections are ordered; each contains a `README` with the specific commands. In brief:
`1_Data_preprocessing` → `2_Network_construction` → `3_Functional_annotation` → (`4`, `5`) evaluation → `6` prioritization → `7` paralogs. Figure scripts live in `figures/`.

## Citation

If you use this code or the associated data, please cite the manuscript _[citation — pending]_ and the Zenodo archive _[DOI — pending]_.
