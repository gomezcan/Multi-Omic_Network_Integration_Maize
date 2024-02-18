# Processing of raw PDI data

### 1. Download fasts file using metadata file: Table_S2
```
./Get_sra_file.sh Table_S2.txt
```

### 2. Clean adapters and low quality reads 
```
# Paired-end reads
cd RawData/PE/
./trimmomatic_PE_job.sh

# Single-read reads
cd RawData/SE/
./trimmomatic_SE_job.sh
```

