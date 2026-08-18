# 2 · Network construction — the four TF→target layers

> One script per layer turns the preprocessed inputs into an edge list. The four edge lists are the central data artifacts of the paper (Zenodo items 1–4). **Manuscript:** Fig 1, Methods.

**Pipeline position:** **⬅ prev** [`1_Data_preprocessing`](../1_Data_preprocessing) · [⌂ overview](../README.md) · **next ➡** [`3_Functional_annotation`](../3_Functional_annotation)

## Layers

| Layer | Folder | Script(s) | Input (from step 1) | Output edge list ⭐ |
|---|---|---|---|---|
| **RFN** — RF-inferred regulatory network (expression; formerly "CEN") | `RFN_expression/` | `Fig_Coexpression.v3.R` | expression matrices + expressed-gene lists | `CoExp_NetworkFinal.10_11_2021.txt` |
| **GRN** — protein–DNA interactions | `GRN_PDI/` | `Fig_PDI.v4.R` (+ `Fig_PDI_newDAP.R` for the new DAP-seq batch) | `Net.Dis2TSS` peak→target assignments, scATAC-overlap z-scores, Table S4 peak QC | `Only_PDI_NetworkFinal.10_14_2022.txt` |
| **eGRN** — cis-eQTL-supported PDIs | `eGRN_cis_eQTL/` | `Fig_ciseQTL.R` | PDI targets × clean cis-eQTL set | `CisE_PDI_NetworkFinal.10_14_2022.txt` |
| **GAN** — trans-eQTL gene association | `GAN_trans_eQTL/` | `Fig_transeQTL.v3.R` | `Clean_trans.eQTL.v2.txt` + `Clean_trans.eQTLp.v2.txt` from step 1's `SNPs_eQTL/Set_CleanFiles_eQTL.R` | `teQTL_NetworkFinal.10_11_2021.txt` |

RF regression follows Zhou et al. 2020 (scikit-learn `RandomForestRegressor`).

## Outputs → where they go

All four edge lists ⭐ are read by `3_Functional_annotation/` (integration + annotation), `5_Evaluation_random_networks/` (rewiring), and `6_Prioritization_and_conditions/` (per-TF summaries).

> **Path note:** scripts reference the original working layout (`Data/Annotations/…`, `../Fig_PDI/…`). Annotation inputs (TF dictionary, syntenic genes, CornCyc pathways, Y1H, …) lived in `Data/Annotations/`; the old `Fig_*` folders map to the sections of this repo (see the [overview](../README.md#step-by-step)).
