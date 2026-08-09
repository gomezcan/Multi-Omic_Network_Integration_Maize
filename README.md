# Multi-Omic Network Integration in Maize

Code and processed-data description for **"Prioritizing Maize Metabolic Gene Regulators through Multi-Omic Network Integration"** (Gomez-Cano et al.).

Four transcription-factor (TF)→target network layers are built from maize expression, protein–DNA-interaction (PDI), and eQTL data, then integrated into a single weighted network whose node embeddings are used to annotate, evaluate, prioritize, and compare metabolic gene regulators.

- **Paper:** _[journal / DOI — pending]_
- **Processed outputs** (edge lists, embeddings, TF→function catalog, prioritization scores): **Zenodo _[DOI — pending]_** — see the [data manifest](#data-artifacts--zenodo-manifest) below
- **License:** MIT (see `LICENSE`)

---

## Pipeline at a glance

Folders `1_ … 7_` are the pipeline, in run order. Every artifact marked ⭐ is (or will be) in the Zenodo deposit.

```mermaid
flowchart TD
    subgraph IN["Inputs (raw / public)"]
        direction LR
        iEXP["Expression panel<br/>(coexpression atlas)"]
        iPDI["Public ChIP- & DAP-seq<br/>(metadata: Table S2)"]
        iQTL["cis- / trans-eQTL sets"]
        iACR["scATAC ACRs<br/>(GEO GSE155178)"]
    end

    subgraph S1["1 · Data preprocessing"]
        direction LR
        p1["PDI: trim → map → filter<br/>→ GEM peaks → peak-to-TSS targets"]
        p2["SNPs & eQTL:<br/>SNP-to-bed → annotate → clean"]
        p3["Expression: expressed-gene<br/>lists · wPCC database"]
    end

    subgraph S2["2 · Network construction — four TF→target layers"]
        direction LR
        n1["RFN ⭐<br/>(RF-inferred, expression)"]
        n2["GRN ⭐<br/>(PDI)"]
        n3["eGRN ⭐<br/>(cis-eQTL-supported PDI)"]
        n4["GAN ⭐<br/>(trans-eQTL)"]
    end

    subgraph S3["3 · Integration, embedding & TF annotation"]
        i1["Integrate layers → weighted union network ⭐"]
        i2["PecanPy node embeddings ⭐<br/>(Dim50 · WL80 · nW100)"]
        i3["MI/MR distance → ClusterONE clusters"]
        i4["TF→function annotation, 3 strategies:<br/>network-based · common-target · common-function<br/>→ TF→function catalog (Table S8) ⭐"]
    end

    s4["4 · Evaluation vs<br/>TF-knockout DEGs"]
    s5["5 · Evaluation vs<br/>random networks"]
    s6["6 · Prioritization (rZ / URS)<br/>+ GSEA condition mapping ⭐"]
    s7["7 · Paralog redundancy<br/>/ divergence"]

    iPDI --> p1
    iACR --> p1
    iQTL --> p2
    iEXP --> p3

    p3 --> n1
    p1 --> n2
    p1 --> n3
    p2 --> n3
    p2 --> n4

    n1 & n2 & n3 & n4 --> i1
    i1 --> i2 --> i3 --> i4

    i4 --> s4
    i4 --> s5
    i4 --> s6
    i2 --> s7
    p3 --> s7

    classDef zen stroke-width:2.5px,stroke-dasharray:0;
    class n1,n2,n3,n4,i1,i2,i4 zen;
```

## Step by step

| Step | Folder | What happens | Key outputs (→ used by) | Manuscript |
|---|---|---|---|---|
| **1** | [`1_Data_preprocessing/`](1_Data_preprocessing) | Modality-specific input prep: full raw ChIP-/DAP-seq chain (download → trim → Bowtie2 → MAPQ30+dedup → GEM peaks → peak→TSS targets), eQTL cleaning/annotation, expressed-gene lists + weighted-PCC (wPCC) database | `Net.Dis2TSS.txt`, clean cis/trans eQTL sets, expressed-gene lists, wPCC DB (→ 2, 3) | Methods |
| **2** | [`2_Network_construction/`](2_Network_construction) | One script per layer builds the four TF→target edge lists | `CoExp_NetworkFinal…txt` ⭐, `Only_PDI_NetworkFinal…txt` ⭐, `CisE_PDI_NetworkFinal…txt` ⭐, `teQTL_NetworkFinal…txt` ⭐ (→ 3, 5, 6) | Fig 1 |
| **3** | [`3_Functional_annotation/`](3_Functional_annotation) | Integrate the four layers, embed nodes with PecanPy, cluster (MI/MR + ClusterONE), annotate TFs by three strategies | `uniqFullNets_weighted.txt` ⭐, PecanPy embeddings ⭐, cluster memberships, PWY/GO enrichments → **TF→function catalog (Table S8)** ⭐ (→ 4, 5, 6) | Figs 1, 3 |
| **4** | [`4_Evaluation_knockouts/`](4_Evaluation_knockouts) | Benchmark the three annotation strategies against TF-knockout DEGs | method-comparison stats; `Summary.Total.Annotation.txt` (→ 6) | Fig 2 |
| **5** | [`5_Evaluation_random_networks/`](5_Evaluation_random_networks) | Rewire networks, re-annotate, compare observed vs random GO recovery | observed-vs-random GO semantic-similarity results | Fig 2 |
| **6** | [`6_Prioritization_and_conditions/`](6_Prioritization_and_conditions) | Rank regulators per process (reciprocal-Z, upstream-regulator score); map TF–GO associations to conditions via GSEA | prioritization scores ⭐, `GSEA_results/` (→ figures) | Figs 3, 4 |
| **7** | [`7_Paralog_redundancy/`](7_Paralog_redundancy) | Predict paralog redundancy/divergence from embedding similarity vs sequence distance | paralog divergence predictions | Fig 5 |
| — | [`figures/`](figures) | Figure-by-figure map to the scripts that render each panel | — | all |
| — | [`environment/`](environment) | Software & package versions | — | Methods |
| — | `archive/` | Superseded script versions, kept for provenance | — | — |

## From the manuscript to the code

| You are reading… | Go to |
|---|---|
| Fig 1 (networks & integration) | `2_Network_construction/` + `3_Functional_annotation/network_based_embedding/` |
| Fig 2 (benchmarking) | `4_Evaluation_knockouts/` + `5_Evaluation_random_networks/` |
| Fig 3 (prioritization by process) | `6_Prioritization_and_conditions/prioritization_rZ_URS/` |
| Fig 4 (condition mapping) | `6_Prioritization_and_conditions/condition_mapping_GSEA/` |
| Fig 5 (paralogs) | `7_Paralog_redundancy/` |
| Table S2 (ChIP/DAP datasets) | `1_Data_preprocessing/PDI/Table_S2.txt` (input metadata) |
| Table S4 (peak sets / reproducibility) | produced along `1_Data_preprocessing/PDI/` steps 5–7 + `peak_coverage_qc/` |
| Table S8 (TF→function catalog) | produced in `3_Functional_annotation/` (annotation strategies) |
| Methods, preprocessing | `1_Data_preprocessing/` (per-modality READMEs) |
| Methods, network inference | `2_Network_construction/README.md` |
| Methods, embedding & clustering | `3_Functional_annotation/network_based_embedding/` |

A finer per-panel map lives in [`figures/README.md`](figures/README.md).

## Data artifacts — Zenodo manifest

The processed artifacts below connect the pipeline to the Zenodo deposit _[DOI — pending]_. "Born in" = the section whose script writes the file.

| # | Artifact (deposit name) | Content | Born in | Consumed by |
|---|---|---|---|---|
| 1 | `CoExp_NetworkFinal.10_11_2021.txt` | RFN edge list (RF-inferred, expression) | 2 · `RFN_expression/` | 3, 5, 6 |
| 2 | `Only_PDI_NetworkFinal.10_14_2022.txt` | GRN edge list (PDI) | 2 · `GRN_PDI/` | 3, 5, 6 |
| 3 | `CisE_PDI_NetworkFinal.10_14_2022.txt` | eGRN edge list (cis-eQTL-supported PDI) | 2 · `eGRN_cis_eQTL/` | 3, 5, 6 |
| 4 | `teQTL_NetworkFinal.10_11_2021.txt` | GAN edge list (trans-eQTL) | 2 · `GAN_trans_eQTL/` | 3, 5, 6 |
| 5 | `uniqFullNets_weighted.txt` | Integrated weighted union network | 3 · `network_based_embedding/` (`1_Set_Pecanpy_weighted.R`) | PecanPy (step 2) |
| 6 | `Pecanpy_uFNetsW.Dim50_WL80_nW100.txt` | Node embeddings (50 dims, walk length 80, 100 walks) | 3 · `network_based_embedding/` (`2_pecanpy_weighted_job.sh`) | 3 (MI/MR), 7 |
| 7 | TF→function catalog (**Table S8**) | TF ↔ predicted pathway/GO functions, all strategies | 3 · annotation outputs (+ `Summary.Total.Annotation.txt`, assembled in 4) | 6 |
| 8 | Prioritization scores (rZ / URS) | Ranked regulators per metabolic process | 6 · `prioritization_rZ_URS/` | Fig 3 |
| 9 | Repo snapshot | This repository, archived | — | — |

> The exact deposit list and file names are frozen at deposit time (Goal 3); this table is the working manifest.

## Reproducing

1. Work through the numbered folders in order; **every section README documents its inputs, scripts in run order, and outputs.**
2. **Path note:** scripts are kept verbatim from the analysis environment for provenance; internal paths still reference the original by-figure working directories (e.g., `../Fig_PDI/`, `Data/Annotations/`). Each section README states which repo folder the old path corresponds to; adjust paths (or symlink) rather than editing scripts.
3. Software versions: [`environment/software_versions.md`](environment/software_versions.md).

### Requirements

- **R** — topGO, GeneOverlap, GOSemSim, Rrvgo, DECIPHER, Parmigene, FGSEA, wCorr, igraph, Rsubread, DESeq2
- **Python** — PecanPy, numpy, pandas, scikit-learn, fcmeans, kneed
- **Command-line** — bcftools, Bowtie2, SAMtools, Picard, Trimmomatic, GEM, bedtools, FastQC, ClusterONE

> **Terminology note:** the expression-based layer built with random-forest regression (following Zhou et al. 2020) is the **RF-inferred regulatory network (RFN)** — previously labeled "co-expression network (CEN)". It predicts regulatory (TF→target) relationships from expression and does not, by itself, imply correlation-based co-expression or causal regulation.

## Citation

If you use this code or the associated data, please cite the manuscript _[citation — pending]_ and the Zenodo archive _[DOI — pending]_.
