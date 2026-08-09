# Figures — manuscript → script map

Most panels are rendered by the `Fig_*.R` script inside the section that produces the underlying result. ([⌂ overview](../README.md))

| Manuscript item | Script / location |
|---|---|
| **Fig 1** — network construction & integration | `2_Network_construction/*/Fig_*.R`, `3_Functional_annotation/network_based_embedding/Fig_integrationV3.R`, `…/tsne/Fig_tsne.R` |
| **Fig 2** — benchmarking (knockouts + random networks) | `4_Evaluation_knockouts/Fig_MethodsComparison.R`, `5_Evaluation_random_networks/4_*` |
| **Fig 3** — prioritization by process | `6_Prioritization_and_conditions/prioritization_rZ_URS/Fig_TFdescription.R` |
| **Fig 4** — condition mapping (GSEA) | `6_Prioritization_and_conditions/condition_mapping_GSEA/Fig_GSEAv2.R` |
| **Fig 5** — paralog redundancy | `7_Paralog_redundancy/Fig_Paralogs.v2.R` |
| Network drawings | `3_Functional_annotation/common_target/Plot_Net_igraph.ipynb` |
| **Table S2** — ChIP/DAP dataset metadata | `1_Data_preprocessing/PDI/Table_S2.txt` (input) |
| **Table S4** — peak sets / QC | `1_Data_preprocessing/PDI/` steps 5–7 + `peak_coverage_qc/` |
| **Table S8** — TF→function catalog ⭐ | `3_Functional_annotation/` strategies + `4_Evaluation_knockouts/Fig_MethodsComparison.R` (`Summary.Total.Annotation.txt`) |
