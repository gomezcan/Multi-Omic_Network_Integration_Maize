# Processing of raw PDI data (ChIP-/DAP-seq)

Public ChIP-seq and DAP-seq datasets (metadata: `Table_S2.txt`) are processed from raw FASTQ to TF→target-gene assignments in the numbered order below. The final outputs (`Net.Dis2TSS.txt` / TF→target lists) are the inputs to `2_Network_construction/GRN_PDI/`. Protoplast ChIP-seq has its own trimming/counting in `protoplast/`.

### 1. Download FASTQ files from SRA
```
./1_Get_fastq_file.sh Table_S2.txt
```

### 2. Trim adapters and low-quality reads (Trimmomatic)
`Adapter.fastq` (TruSeq3 adapters) must sit in the working directory.
```
./2.1_trimmomatic_SE_job.sh      # single-end
./2.2_trimmomatic_PE_job.sh      # paired-end
```

### 3. Map to B73 RefGen_v4 (Bowtie2 → sorted BAM)
Bowtie2 2.3.5.1 against the `Index_B73v4.dna_bowtie2` index; PE mapping adds `--no-mixed --no-discordant`; SAM→sorted BAM via SAMtools.
```
./3.1_MappingBowtie2_SE.sh
./3.2_MappingBowtie2_PE.sh
```

### 4. Filter and deduplicate alignments
Keep MAPQ ≥ 30 (removes multi-mappers), then Picard `MarkDuplicates REMOVE_DUPLICATES=true` → `DeDup.Q30.*.bam`.
```
./4.1_CleanMapped_Reads_SE.sh
./4.2_CleanMapped_Reads_PE.sh
```

### 5. Peak calling (GEM/GPS)
GEM (`gem.jar`) with `--relax`, k-mer range 6–15, default read distribution, B73 genome/chromosome files. The job script contains two modes: Mode 1 (gDNA control) and Mode 2 (no control; the mode used for the final peak sets). Run per sample, including the DAP methylation states (`_deM` / `_Met`).
```
./5_PeakCalling_GEM.job.sh
```

### 6–8. Peaks → target genes
```
./6_Peaks_to_Targets_disTSS.sh   # peak midpoint → closest TSS (bedtools closest vs
                                 #   TSS_B73_RefGen_v4.46.bed, signed distance)
                                 #   → Dis2TSS.<TF>.txt → concatenated Net.Dis2TSS.txt
./7_CountPeak_Targets.sh         # per-sample summary: peaks and target genes (|dist| ≤ 2 kb)
./8_GetNet_FromTargets.sh        # assemble the TF→target edge list from targets files
```

### QC: read coverage in peaks (`peak_coverage_qc/`)
```
./1_Get_saf.sh                   # GEM narrowPeak → SAF annotation
./2_TPMs_Peaks_by_Sample.R       # Rsubread::featureCounts per sample vs All.peaks.saf
                                 #   → per-peak CPM/TPM tables (peak reproducibility QC)
```

### Provenance
Steps 3–8, `Adapter.fastq`, and `peak_coverage_qc/` were recovered verbatim from the original processing archive (eg-server `MaizeENCODE/` project); only filenames were normalized to this numbered pipeline:

| Repo file | Original (`MaizeENCODE/…`) |
|---|---|
| `3.1_MappingBowtie2_SE.sh` | `Scripts/MappingBowtie2_SingleEnd.sh` |
| `3.2_MappingBowtie2_PE.sh` | `Scripts/MappingBowtie2_PairEnd.sh` |
| `4.1_CleanMapped_Reads_SE.sh` | `Scripts/CleanMapped_Reads_single.sh` |
| `4.2_CleanMapped_Reads_PE.sh` | `Scripts/CleanMapped_Reads_paired.sh` |
| `5_PeakCalling_GEM.job.sh` | `Scripts/PeakCallin_deMET.job.sh` |
| `6_Peaks_to_Targets_disTSS.sh` | `Targets/1_Peaks.2.Targets_disTSS.sh` |
| `7_CountPeak_Targets.sh` | `Targets/2_CountPeak_Targets.sh` |
| `8_GetNet_FromTargets.sh` | `Targets/GetNet_FromTargets.files.sh` |
| `peak_coverage_qc/1_Get_saf.sh` | `QC/PeaksCoverage/1_Get_saf.sh` |
| `peak_coverage_qc/2_TPMs_Peaks_by_Sample.R` | `QC/PeaksCoverage/2_TPMs_Peaks_by_Sample.R` |
| `Adapter.fastq` | `Scripts/Adapter.fastq` |

Genome resources referenced by the scripts (not shipped here): Bowtie2 index `Index_B73v4.dna_bowtie2`, `TSS_B73_RefGen_v4.46.bed`, `GenomeSize_B73.sizes`, GEM `Read_Distribution_default.txt`, and per-chromosome FASTA (`ChrsB73/`).
