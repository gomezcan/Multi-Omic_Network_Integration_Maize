# 1 · Data preprocessing
Raw and processed input preparation, by modality. Outputs feed `2_Network_construction/`.
- **`PDI/`** — ChIP-/DAP-seq: download (`1_Get_fastq_file.sh`), trim (`2.1/2.2_trimmomatic_*_job.sh`); `protoplast/` = protoplast ChIP-seq trimming + read counts. Metadata: `Table_S2.txt`.
- **`SNPs_eQTL/`** — SNP→bed (`1_From_SNP_2_bed.sh`, `2.1_FarmToBed.sh`), eQTL annotation/cleaning (`2_annotate_trasn_cis_eQTL.sh`, `Set_CleanFiles_eQTL.R`, `Clean_00eQTLs.sh`).
- **`Expression_coexpression/`** — expressed-gene lists (`expression/`) and weighted-PCC database (`wPCC_DB/`).
