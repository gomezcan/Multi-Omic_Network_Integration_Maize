# 5 · Evaluation against random networks

> Degree-preserving rewiring of each network layer, re-annotation, and comparison of observed vs random GO recovery — separately for each of the three annotation strategies. **Manuscript:** Fig 2.

**Pipeline position:** **⬅ prev** [`3_Functional_annotation`](../3_Functional_annotation) · [⌂ overview](../README.md) · **next ➡** [`6_Prioritization_and_conditions`](../6_Prioritization_and_conditions)

## Inputs
The four edge lists ⭐ (step 2), the integrated network / embedding pipeline (step 3, re-run on rewired networks), GO annotations (GAMER).

## Steps (run in order; three parallel strategy tracks)

| # | Scripts | What happens |
|---|---|---|
| 0–1 | `0_RewireNet.R`, `1_RewireNet.FullMR.R` | rewire the networks (igraph, degree-preserving); full-MR variant for the embedding track |
| random nets + GOs | `Get_RandomNet_and_GOs.R`, `Get_RandomCoExpNet_and_GOs.R`, `Get_RandomteQTLNet_and_GOs.R` | generate random networks per layer and their GO annotations |
| track "network-based" | `1_1_MapParentGOs_Netbase.{R,sh}`, `4_SetInput_wNet_Pencanpy.random.R` | map GO parents; rebuild embedding input from rewired networks |
| track "common-function" | `2_GO_targets.R`, `2_1_GSS_from_randomNets.R`, `2_2_MapParentGOs_CommonFunct.{R,sh}` | re-annotate targets on random networks; GO semantic similarity (GSS) |
| track "common-target" | `3_GO_Common_targets.R`, `3_1_MapParentGOs_Common_targets.{R,sh}` | same for the common-target strategy |
| GSS null | `Get_GO_SS_random.{R,parallel.R}`, `Get_GO_SS_random.teQTL.parallel.R` | random GO semantic-similarity distributions (GOSemSim) |
| 4 | `4_RandomNets_analysis.R` + `4_1_Obs_vs_Random_GSS_results/4_1_GSS_Obs_vs_random_{nbase,cFunc,cTarg}.R` | observed vs random GSS per strategy → Fig 2 |

## Outputs
Observed-vs-random GO semantic-similarity tables per strategy → Fig 2.
