## Make annotation table to merge with eQTL results later
library(data.table)
library(tidyverse)
library(rtracklayer)

## Read in gff3 file downloaded from maizeGDB then filter rows with a gene id
setwd("~/eqtl_pdi_coexpr/data/06_annotation_data/")
gff = import('Zea_mays.B73_RefGen_v4.50.gff3.gz')
gff_filtered = as_tibble(gff) %>% 
  drop_na(gene_id) %>% 
  select(seqnames, start, end, type, gene_id, description) %>% 
  filter(!grepl("B73",seqnames)) %>% 
  filter(!grepl("Mt",seqnames)) %>% 
  filter(!grepl("Pt",seqnames)) %>% 
  mutate(seqnames = as.numeric(seqnames)) %>% 
  arrange(seqnames, start) %>% 
  setNames(.,c("gene_chrom","gene_TSS","gene_TTS","gene_type","gene_id","description")) %>% 
  mutate(description = ifelse(is.na(description), "not_annotated",description))

## Identify syntenic genes based on Schnable 2019: 10.6084/m9.figshare.7926674.v1
syntenic = fread("sorghum3_intell_plusteff.csv") %>% 
  select(maize1_v4,maize2_v4) %>% 
  mutate(check = maize1_v4==maize2_v4) %>% 
  filter(!check == TRUE) %>% 
  select(-check) %>% 
  pivot_longer(cols = c(maize1_v4,maize2_v4)) %>% 
  filter(!value=="No Gene") %>% 
  setNames(.,c("syntenic_subgenome","gene_id")) %>% 
  mutate(syntenic_subgenome = ifelse(syntenic_subgenome=="maize1_v4","sub_genome1","sub_genome2")) %>% 
  mutate(gene_synteny = "syntenic") %>% 
  distinct(.)

## Combine syntenic information with gff file
gff_filtered_with_syntenic = gff_filtered %>% 
  left_join(., syntenic) %>% 
  mutate(gene_synteny = ifelse(is.na(gene_synteny), "non_syntenic",gene_synteny))

fwrite(gff_filtered_with_syntenic,"Zea_mays.B73_RefGen_v4.48.gff3.with.desc.synteny.txt",sep = "\t", row.names = FALSE, col.names = TRUE)

