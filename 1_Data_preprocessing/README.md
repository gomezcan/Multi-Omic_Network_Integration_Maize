# 1 · Data preprocessing

> Modality-specific preparation of every raw/public input. Outputs feed [`2_Network_construction/`](../2_Network_construction). **Manuscript:** Methods.

**Pipeline position:** [⌂ overview](../README.md) · **next ➡** [`2_Network_construction`](../2_Network_construction)

## Subsections (independent; run any order)

### `PDI/` — ChIP-/DAP-seq, full raw chain
Numbered steps 1→8: download (SRA) → trim (Trimmomatic + `Adapter.fastq`) → map to B73 RefGen_v4 (Bowtie2) → MAPQ30-filter + dedup (Picard) → GEM peak calling → peak→target-gene assignment by TSS distance → TF→target lists. `peak_coverage_qc/` computes per-peak CPMs (featureCounts) for reproducibility QC; `protoplast/` trims/counts the protoplast ChIP-seq. Dataset metadata: `Table_S2.txt`.
**Step-by-step commands & provenance:** [`PDI/README.md`](PDI/README.md).

### Open chromatin (ACR) — external data, no raw processing
The ACR layer uses published B73 scATAC ACRs (GEO **GSE155178**), overlapped with the peak sets for QC — nothing to run here.

### `eQTL_calling/` — eQTL identification (upstream of `SNPs_eQTL/`)
The genotype panel, expression/PEER preparation, permutation thresholds and the genome-wide
**eight-model** rMVP scan across the eight tissues, ending in the filtered eQTL results table that
`SNPs_eQTL/` consumes. cis and trans associations are thresholded and filtered differently
(per-gene vs. per-tissue thresholds; 8-of-8 models for the eGRN set vs. ≥ 2 for the GAN).
**Run order, the model table and known gaps:** [`eQTL_calling/README.md`](eQTL_calling/README.md).

### `SNPs_eQTL/` — SNP and eQTL preparation
| # | Script | What it does |
|---|---|---|
| 1 | `1_From_SNP_2_bed.sh` | GEM `narrowPeak` peak calls → per-TF summit BEDs (builds the PDI reference `All.Summit_*.bed` used by steps 2 and 3b) |
| 2 | `2_annotate_trasn_cis_eQTL.sh` | annotate an eQTL BED by gene-body overlap (`intersectBed`) and by distance to the nearest PDI peak summit (`closestBed`) |
| 3a | `Set_CleanFiles_eQTL.R` | classify eQTLs into the five S2A Fig classes → `Clean_*.v2.txt` sets (→ GAN, `2_…/GAN_trans_eQTL/Fig_transeQTL.v3.R`) |
| 3b | `Clean_00eQTLs.sh` | cis branch: eQTLs supported by all 8 models and ≤50 kb from the target TSS, matched to PDI summits ≤20 bp away → `cis_eQTL.pdi.network.txt` (→ eGRN, `2_…/eGRN_cis_eQTL/Fig_ciseQTL.R`) |

> **Provenance note.** Both branches start from the eQTL results table produced upstream on the
> original compute cluster by the genome-wide eight-model rMVP scan (Methods). The intermediate
> inputs read by `Set_CleanFiles_eQTL.R` (`Final_cis_trans_all_eQTL_10012021.bed`,
> `trans.eQTL*_target.txt`, `cis.eQTLt_target.txt`, `ciseQTL_noFiter.txt`) came from one-off
> `2_annotate_trasn_cis_eQTL.sh`-style bedtools runs whose exact invocations were not preserved;
> source genes were defined as gene body or ≤2 kb upstream of the TSS (S2A Fig). Genome resources
> referenced but not shipped: `All.Summit_10.2020.bed`, `Zea_mays.B73_RefGen_v4.46.bed`,
> `Zm.v4.synteny.genes.txt`.

### `Expression_coexpression/` — expression inputs
| Where | Notebook / script | What it does |
|---|---|---|
| `expression/` | `Get.Expressed.Genes.list.ipynb` | per-network expressed-gene lists |
| `wPCC_DB/` | `FromMatrix_2_Gene.ipynb` | split the weighted-PCC (wPCC) matrix into the per-gene coexpression database used downstream (uses the expressed-gene list from `expression/`) |

## Outputs → where they go

| Artifact | Produced by | Consumed by |
|---|---|---|
| `Net.Dis2TSS.txt` (peak→TSS assignments) + TF→target lists | `PDI/` steps 6–8 | `2_…/GRN_PDI/`, `2_…/eGRN_cis_eQTL/` |
| `00_filtered_cis_trans_all_eQTL_results.txt` (eight-model eQTL scan) | `eQTL_calling/` | `SNPs_eQTL/` |
| clean cis- / trans-eQTL BED sets | `SNPs_eQTL/` | `2_…/eGRN_cis_eQTL/`, `2_…/GAN_trans_eQTL/` |
| expressed-gene lists | `expression/` | `2_…/RFN_expression/`, filters throughout |
| wPCC per-gene database | `wPCC_DB/` | layer weighting in `3_…/network_based_embedding/`, `6_…` (GSEA), `7_…` (paralogs) |
