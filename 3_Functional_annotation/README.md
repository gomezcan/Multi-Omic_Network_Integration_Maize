# 3 · Integration, embedding & TF functional annotation

> The core of the method: the four layers are merged into one weighted network, nodes are embedded with PecanPy, embedding distances define clusters, and TFs are annotated by **three strategies** (network-based, common-target, common-function) → the TF→function catalog (**Table S8**). **Manuscript:** Figs 1 & 3, Methods.

**Pipeline position:** **⬅ prev** [`2_Network_construction`](../2_Network_construction) · [⌂ overview](../README.md) · **next ➡** [`4_Evaluation_knockouts`](../4_Evaluation_knockouts), [`5_Evaluation_random_networks`](../5_Evaluation_random_networks), [`6_Prioritization_and_conditions`](../6_Prioritization_and_conditions)

## Inputs
The four edge lists ⭐ from step 2 (`CoExp_…`, `Only_PDI_…`, `CisE_PDI_…`, `teQTL_NetworkFinal…`), the wPCC database (step 1), and the annotation bundle (TF dictionary, syntenic genes, CornCyc, GO/GAMER).

## A · `network_based_embedding/` — integrate → embed → cluster (run in order)

| # | Script / notebook | What it does | Output |
|---|---|---|---|
| 0 | `Fig_integrationV3.R` (support: `wPCC_m_GAN_Random.parallel.R`; variant: `Fig_Integration_CLR.R` = CLR benchmark) | integrate/weight the four layers (wPCC-based edge weights; GAN random background) | integration stats (Fig 1) |
| 1 | `1_Set_Pecanpy_weighted.R` | merge the four edge lists into one weighted union network | `uniqFullNets_weighted.txt` ⭐ |
| 2 | `2_pecanpy_weighted_job.sh` | PecanPy node2vec on the weighted network | `Pecanpy_uFNetsW.Dim50_WL80_nW10.txt` ⭐ (50 dims, walk length 80, 10 walks) |
| 3 | `3_OptimumK_pepanpy.ipynb` 📓 | **decision record**: fuzzy c-means + knee-point (kneed) selection of the cluster number *K* on the embedding | chosen *K* |
| — | `CosMatrix_pecanpy.ipynb` 📓 | cosine-similarity matrix of the embeddings | cosine matrix |
| — | `Selecting_N_clusters.ipynb` 📓 | companion K-selection analysis (earlier variant of #3) | — |
| 4 | `DistanceCalculation/1_MI_and_MR_Distance.R` | mutual-information + mutual-rank distance between embedded genes (Parmigene) | `MR_edgesDB_…/MR_MI.pecanpy.<gene>.txt` |
| 5 | `DistanceCalculation/2_create_ClusterOneInput.sh` (+ `Filter_pVal.R`) | assemble/filter the MR-MI edge database | `InputClusterONE_…_syntenic.txt` |
| 6 | `DistanceCalculation/3_ClusterOne_job.sh` | retained for provenance only: the published analysis does not use ClusterONE; each TF's "cluster" is its MRMI neighborhood (D ≥ 0.005) taken directly from the step-5 edge list | — |
| 7 | `Fig_pecanpyV6.R`, `Fig_PecanpyPart2.R` | per-cluster GO/PWY enrichment; TF–GO network assembly | `BP_results_targets/`, `TFGO_4337_net.rds`, `MaizeSyntenicGenes_GOparent.rds` |
| viz | `tsne/` (`tsne_Adj_matrix_job.R`, `Fig_tsne.R`, `OptimumK_tsne.ipynb` 📓, `Get_tsneDataExample.R`) | t-SNE projection of the embedding + K-selection check in t-SNE space | Fig 1 panels |

## B · Annotation strategies (independent, compared in step 4)

| Strategy | Folder | Script(s) | Idea |
|---|---|---|---|
| **Network-based** | (A above) | cluster enrichments from steps 6–7 | TF inherits functions enriched in its embedding cluster |
| **Common-target** | `common_target/` | `Fig_CommonTarg.v2.R` (+ `Plot_Net_igraph.ipynb` 📓 = network drawing from `Full_Final_network.11022022.txt`) | TFs sharing targets share functions |
| **Common-function** | `common_function/` | `GOs_CoExp.R`, `GOs_PDI.R`, `GOs_teQTL.R` | per-layer GO enrichment of each TF's targets |
| **PWY/GO catalog** | `annotation_PWY_GO/` | `Fig_Function_PWYs_GOs.v3.R` (+ `Source_…v2.R`) | consolidate CornCyc-pathway & GO annotations → per-strategy enrichment tables |

## C · `communities/` — alternative community detection (supporting)
`Get_net_communities.py`, `Communities_Based_cmean.R`, `Community_Prediction.R`, `Identify_communities_by_tsne.R`, `testing_community_identification.ipynb` 📓 (igraph exploration).

## Outputs → where they go

| Artifact | Consumed by |
|---|---|
| `uniqFullNets_weighted.txt` ⭐, PecanPy embeddings ⭐ | steps 4–7 (`7_Paralog_redundancy` reads the embedding distances directly) |
| per-strategy PWY/GO enrichment tables (`*_PWY_enrichment.txt`, `*_GO_enrichment.txt`) | `4_Evaluation_knockouts/Fig_MethodsComparison.R` → assembled into `Summary.Total.Annotation.txt` → **Table S8** ⭐ |
| `TFGO_4337_net.rds`, `MaizeSyntenicGenes_GOparent.rds` | `6_…/condition_mapping_GSEA/` |

📓 = Jupyter notebook — parameter-selection / visualization decision records, kept so reviewers can audit the choices (e.g., number of clusters).

> **Path note:** scripts reference the original layout — `../Fig_PDI/`, `../Fig_Coexpression/`, `../Fig_transeQTL/` = section 2 folders; `../Fig_pecanpy/` = this folder; `Data/Annotations/` = the annotation bundle.
