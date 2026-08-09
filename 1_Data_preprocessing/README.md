# 1 · Data preprocessing

> Modality-specific preparation of every raw/public input. Outputs feed [`2_Network_construction/`](../2_Network_construction). **Manuscript:** Methods.

**Pipeline position:** [⌂ overview](../README.md) · **next ➡** [`2_Network_construction`](../2_Network_construction)

## Subsections (independent; run any order)

### `PDI/` — ChIP-/DAP-seq, full raw chain
Numbered steps 1→8: download (SRA) → trim (Trimmomatic + `Adapter.fastq`) → map to B73 RefGen_v4 (Bowtie2) → MAPQ30-filter + dedup (Picard) → GEM peak calling → peak→target-gene assignment by TSS distance → TF→target lists. `peak_coverage_qc/` computes per-peak CPMs (featureCounts) for reproducibility QC; `protoplast/` trims/counts the protoplast ChIP-seq. Dataset metadata: `Table_S2.txt`.
**Step-by-step commands & provenance:** [`PDI/README.md`](PDI/README.md).

### Open chromatin (ACR) — external data, no raw processing
The ACR layer uses published B73 scATAC ACRs (GEO **GSE155178**), overlapped with the peak sets for QC — nothing to run here.

### `SNPs_eQTL/` — SNP and eQTL preparation
| # | Script | What it does |
|---|---|---|
| 1 | `1_From_SNP_2_bed.sh` | SNP tables → BED |
| 2 | `2_annotate_trasn_cis_eQTL.sh` | annotate eQTLs as cis / trans |
| 2.1 | `2.1_FarmToBed.sh` | FarmCPU output → BED |
| 3 | `Set_CleanFiles_eQTL.R` + `Clean_00eQTLs.sh` | final cleaning of the eQTL sets |

### `Expression_coexpression/` — expression inputs
| Where | Notebook / script | What it does |
|---|---|---|
| `expression/` | `Get.Expressed.Genes.list.ipynb` | per-network expressed-gene lists |
| `wPCC_DB/` | `FromMatrix_2_Gene.ipynb` | split the weighted-PCC (wPCC) matrix into the per-gene coexpression database used downstream (uses the expressed-gene list from `expression/`) |

## Outputs → where they go

| Artifact | Produced by | Consumed by |
|---|---|---|
| `Net.Dis2TSS.txt` (peak→TSS assignments) + TF→target lists | `PDI/` steps 6–8 | `2_…/GRN_PDI/`, `2_…/eGRN_cis_eQTL/` |
| clean cis- / trans-eQTL BED sets | `SNPs_eQTL/` | `2_…/eGRN_cis_eQTL/`, `2_…/GAN_trans_eQTL/` |
| expressed-gene lists | `expression/` | `2_…/RFN_expression/`, filters throughout |
| wPCC per-gene database | `wPCC_DB/` | layer weighting in `3_…/network_based_embedding/`, `6_…` (GSEA), `7_…` (paralogs) |
