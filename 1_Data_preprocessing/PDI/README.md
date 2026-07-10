# Processing of raw PDI data

### 1. Download fasts file using metadata file: Table_S2
```
./Get_sra_file.sh Table_S2.txt
```

### 2. Clean adapters and low quality reads 
```
# Single-read reads
cd RawData/SE/
./2.1_trimmomatic_SE_job.sh

# Paired-end reads
cd RawData/PE/
./2.2_trimmomatic_PE_job.sh
```

