# 3 · TF functional annotation
Three integration strategies plus the network embeddings.
- **`common_target/`** — common-interaction annotation (`Fig_CommonTarg.v2.R`).
- **`common_function/`** — per-layer GO enrichment (`GOs_CoExp.R`, `GOs_PDI.R`, `GOs_teQTL.R`).
- **`annotation_PWY_GO/`** — PWY/GO annotation (`Fig_Function_PWYs_GOs.v3.R` + `Source_*`).
- **`network_based_embedding/`** — core method: integrate layers (`Fig_integrationV3.R`, `Fig_Integration_CLR.R`), PecanPy embeddings (`1_Set_Pecanpy_weighted.R`, `2_pecanpy_weighted_job.sh`, `3_OptimumK_pepanpy.ipynb`), MI/MR distance + ClusterONE (`DistanceCalculation/`), t-SNE (`tsne/`).
- **`communities/`** — network community detection.
