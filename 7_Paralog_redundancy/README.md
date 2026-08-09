# 7 · Paralog redundancy / divergence

> Predict functional redundancy vs divergence of TF paralog pairs by comparing network-embedding similarity with amino-acid-sequence distance and coexpression, against a random-pair null. **Manuscript:** Fig 5.

**Pipeline position:** **⬅ prev** [`6_Prioritization_and_conditions`](../6_Prioritization_and_conditions) · [⌂ overview](../README.md)

## Inputs
Paralog pairs (`Paralogs_Schable_Lab.csv`), embedding MR/MI distances (`InputClusterONE_…_syntenic.txt`, step 3), full final network (`Full_Final_network.11022022.txt`), wPCC DB (step 1).

## Steps

| # | Where | Script(s) | What it does |
|---|---|---|---|
| 1 | `AAsequenceDistance/` | `1_Keep_longest_pep.R` → `Extract_seq{,.len}.py` → `2_Pep_distance_Calculation.{R,sh}` | longest peptide per gene → pairwise AA distance (DECIPHER Hamming) |
| 2 | `wPCC_Distance/` | `1_Get_wPCC_from_files.sh` | wPCC coexpression distance per paralog pair |
| 3 | `Random_SCC_DB/` | `Spearman_RandomPairs.{R,sh}` | Spearman-correlation null from random gene pairs |
| 4 | (root) | `Fig_Paralogs.v2.R` | combine embedding similarity × AA distance × wPCC vs null → redundancy/divergence calls; Fig 5 panels |

## Outputs
Paralog redundancy/divergence predictions → Fig 5.
