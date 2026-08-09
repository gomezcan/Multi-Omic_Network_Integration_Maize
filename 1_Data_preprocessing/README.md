# 1 · Data preprocessing
Raw and processed input preparation, by modality. Outputs feed `2_Network_construction/`.
- **`PDI/`** — ChIP-/DAP-seq, full raw chain: download (`1_…`), trim (`2.x_…` + `Adapter.fastq`), map to B73v4 (`3.x_…` Bowtie2), MAPQ30-filter + dedup (`4.x_…`), GEM peak calling (`5_…`), peak→target-gene assignment by TSS distance (`6–8_…`); `peak_coverage_qc/` = per-peak CPM QC (featureCounts); `protoplast/` = protoplast ChIP-seq trimming + read counts. Metadata: `Table_S2.txt`.
- **Open chromatin (ACR)** — no raw processing here: the ACR layer uses published B73 scATAC ACRs (GEO **GSE155178**), overlapped with the peak sets for QC.
- **`SNPs_eQTL/`** — SNP→bed (`1_From_SNP_2_bed.sh`, `2.1_FarmToBed.sh`), eQTL annotation/cleaning (`2_annotate_trasn_cis_eQTL.sh`, `Set_CleanFiles_eQTL.R`, `Clean_00eQTLs.sh`).
- **`Expression_coexpression/`** — expressed-gene lists (`expression/`) and weighted-PCC database (`wPCC_DB/`).
