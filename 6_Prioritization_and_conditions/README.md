# 6 · Prioritization & condition mapping

> Rank candidate regulators per metabolic process (reciprocal-Z score, upstream-regulator score) and map TF–GO associations to biological conditions via GSEA. **Manuscript:** Figs 3 & 4.

**Pipeline position:** **⬅ prev** [`4_Evaluation_knockouts`](../4_Evaluation_knockouts) / [`5_…`](../5_Evaluation_random_networks) · [⌂ overview](../README.md) · **next ➡** [`7_Paralog_redundancy`](../7_Paralog_redundancy)

## `prioritization_rZ_URS/` — regulator ranking (Fig 3)

| Script | What it does | Inputs | Output |
|---|---|---|---|
| `Fig_TFdescription.R` | per-TF description across the four layers; process-level prioritization (rZ / URS) | the four edge lists ⭐ (step 2), `Summary.Total.Annotation.txt` (step 4), expression DB | prioritization scores ⭐; Fig 3 panels |

## `condition_mapping_GSEA/` — condition mapping (Fig 4)

| Script | What it does | Inputs | Output |
|---|---|---|---|
| `Fig_GSEAv2.R` | FGSEA of TF–GO associations across expression conditions | `TFGO_4337_net.rds`, `MaizeSyntenicGenes_GOparent.rds` (step 3), wPCC DB (step 1), `ExpressionSamples_annotation.txt` | `GSEA_results/GSEA_GOs.<TF>.txt`; Fig 4 panels |

> **Path note:** `../Fig_MethodsComparison/` = section 4; `../Fig_PecanpyPart2/` = section 3 `network_based_embedding/`.
