# Software & package versions

Versions used for the analyses in Gomez-Cano et al. Pin these (or use the container/renv once added) to reproduce results.

## Reference resources
- Maize reference genome: **B73 AGPv4**
- Metabolic pathways: **CornCyc** (MaizeGDB)
- GO annotations: **GAMER** (maize-GAMER)
- Syntenic gene set: maize–*Sorghum bicolor* (Schnable et al. 2019)

## Command-line tools
| Tool | Version | Used in |
|------|---------|---------|
| FastQC | 0.11.5 | PDI read QC |
| Trimmomatic | _[pin version]_ | PDI read trimming |
| Bowtie2 | 2.3.5.1 (module loaded in `3.x_MappingBowtie2_*` scripts; Methods draft said 2.3.4.1 — reconcile) | PDI read mapping |
| SAMtools | 1.9 | PDI alignment filtering |
| Picard | _[pin version]_ | PDI duplicate removal (MarkDuplicates) |
| GEM | 3.4 | PDI peak calling |
| bedtools | _[pin version]_ | PDI peak→TSS target assignment |
| bcftools | 1.7 | SNP dataset concatenation |
| PecanPy | _[pin version]_ | network embeddings |

## R packages
| Package | Version | Used in |
|---------|---------|---------|
| Rsubread | 1.32.2 | peak CPM counting |
| GeneOverlap | 1.30.0 | PWY enrichment |
| topGO | 2.46.0 | GO enrichment |
| GOSemSim | 2.20.0 | GO semantic similarity |
| Rrvgo | 1.6 | GO parent mapping |
| DECIPHER | 2.22 | paralog AA Hamming distance |
| Parmigene | _[pin version]_ | mutual information (embeddings) |
| FGSEA | 1.20 (1.18.0 for condition-GSEA) | gene-set enrichment |
| wCorr | 1.9.1 | weighted PCC |
| igraph | 1.2.4.1 | random-network rewiring |
| DESeq2 | _[pin version]_ | knockout DEGs |

## Python
| Package | Version | Used in |
|---------|---------|---------|
| scikit-learn (RandomForestRegressor) | _[pin version]_ | RFN construction (per Zhou et al. 2020) |
| PecanPy | _[pin version]_ | node embeddings |
| numpy / pandas | _[pin versions]_ | data handling |

> Items marked _[pin version]_ were not stated in the original Methods — recover from the analysis environment (`sessionInfo()` / `pip freeze` / `conda list`) before submission. A frozen `renv.lock` and/or `environment.yml` should be added here.
