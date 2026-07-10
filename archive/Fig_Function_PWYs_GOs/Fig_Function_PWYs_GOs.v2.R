library(extrafont)
library(scales)
library(tidyverse)
library(data.table)
library(ggVennDiagram)
library(purrr)
library(gplots)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(viridis)
library(patchwork)
library(reshape2)
#
library(rrvgo)
library(topGO)
library(GOSemSim)
library(enrichplot)
library(GeneOverlap)
library(ComplexHeatmap)
library(circlize)
library(parallel)
library(org.Zmays.eg.db)


##########################################################
######                  Functions                   ######
##########################################################

source("Source_Fig_Function_PWYs_GOs.v2.R")

##################################################
##########        Annotations       ##############
##################################################

# saf <- as_tibble(read.table("Data/eQTL_data/Zea_mays.B73_RefGen_v4.46.saf", stringsAsFactors = F))
# saf1 <- subset(saf, V5=="+")[,c(1,2,3)]
# saf2 <- subset(saf, V5=="-")[,c(1,2,4)]
# colnames(saf1) <- c("GeneID", "chrAnn", "TSS")
# colnames(saf2) <- c("GeneID", "chrAnn", "TSS")
# #
# saf <- rbind(saf1, saf2)

# PDI
PDI <- unique(fread("../Fig_PDI/Only_PDI_NetworkFinal.10_14_2022.txt")[,c(2,3)])
colnames(PDI)[1] <- "TF"

# PDI eQTL
PDIeQTL <- unique(fread("../Fig_PDI/CisE_PDI_NetworkFinal.10_14_2022.txt")[,c(2,3)])
colnames(PDIeQTL)[1] <- "TF"

# PDI <- rbind(PDI, PDIeQTL)
# CoExp
CoExp <- unique(fread("../Fig_Coexpression/CoExp_NetworkFinal.10_11_2021.txt"))
CoExp <- unique(CoExp[,2:3])
colnames(CoExp)[1] <- "TF"

# teQTL
teQTL <- unique(fread("../Fig_transeQTL/teQTL_NetworkFinal.10_11_2021.txt"))
colnames(teQTL)[1] <- "TF"

# Genes in synteny
GenesMaize <- unique(as.character(read.table("Data/Annotations/Zm.v4.synteny.genes.txt", h=T)$gene_id))
length(GenesMaize)

# TF names
TFdic <- as_tibble(read.table("Data/Annotations/TF_Id_Name.txt", h=T, stringsAsFactors = F))

# TF and CoReg
TF_CoR <- as_tibble(read.table("Data/Annotations/TF_CoR_Mazie.txt", h=T, stringsAsFactors = F)) 

### Y1H network
# Y1H <- as_tibble(read.table("Data/Annotations/Y1H.data.txt", h=T, stringsAsFactors = F))[,4:6]
# ReplaceName(Y1H$TF.v4)

# CornCYC
CornCYC <- as_tibble(read.table("Data/Annotations/corn_pathways.0210325.reduced.txt", h=T, stringsAsFactors = F))
CornCYC <- CornCYC[!(CornCYC$Gene.id %in% "unknown"),]
CornCYC <- subset(CornCYC, Gene.id %in% Syntenic)
CornCYC[,"Class"] <- "Enzyme"
colnames(CornCYC)[3] <- "GeneID"
CornCYC$Pathway.id <- gsub("-", "_", CornCYC$Pathway.id)
#
CornCYC$Pathway.name <- gsub("</i>", "", gsub("<i>", "", CornCYC$Pathway.name))
CornCYC$Pathway.name <- gsub("UDP-&alpha;-D-xylose", "UDP.alpha.D_xylose", CornCYC$Pathway.name)
CornCYC$Pathway.name <- gsub("<sup>", ".", CornCYC$Pathway.name)
CornCYC$Pathway.name <- gsub("</sup>", "", CornCYC$Pathway.name)

CornCYC.list <- split(CornCYC$GeneID, CornCYC$Pathway.id)
CornCYC_size <- as_tibble(as.data.frame(table(unique(CornCYC[,c(1,3)])$Pathway.id), stringsAsFactors = F))
colnames(CornCYC_size) <- c("PWY", "nPWY")

All_TFs <- unique(c(unique(PDI$Source), CoExp$Source, TF_CoR$GeneID))

## syntenic genes 
Syntenic <- as_tibble(read.table("Data/Annotations/Zm.v4.synteny.genes.txt", h=T, stringsAsFactors = F, sep = '\t'))
Syntenic <- subset(Syntenic, gene_synteny=='syntenic')$gene_id


# GOs term annotations
background <- readMappings("synteny.ID_TopGO_V4_GRAMER.txt")
background_IDs <- as.character(unique(names(background)))
background_list <- unique(as.character(unlist(background)))

##################################################

####################################################################
########       Enrichment of PWY in targets by TF       ############
####################################################################

# ## PYW Enrichment 
PDI <- unique(rbind(PDI, PDIeQTL))
GRN_PWY <- Enrichmet_classes(unique(rbind(PDI, PDIeQTL)))
GRN_PWY$PWY <- as.character(GRN_PWY$PWY)
write.table(GRN_PWY, "PWY_results/GRN_PWY_enrichment.txt", sep = '\t', row.names = F)

#
eGRN_PWY <- Enrichmet_classes(PDIeQTL)
eGRN_PWY$PWY <- as.character(eGRN_PWY$PWY)
write.table(eGRN_PWY, "PWY_results/eGRN_PWY_enrichment.txt", sep = '\t', row.names = F)

#
CEN_PWY <- Enrichmet_classes(CoExp)
write.table(CEN_PWY, "PWY_results/CEN_PWY_enrichment.txt", sep = '\t', row.names = F)
CEN_PWY$PWY <- as.character(CEN_PWY$PWY)

#
GAN_PWY <- Enrichmet_classes(subset(teQTL, TF %in% All_TFs))
GAN_PWY$PWY <- as.character(GAN_PWY$PWY)
write.table(GAN_PWY, "PWY_results/GAN_PWY_enrichment.txt", sep = '\t', row.names = F)

GRN_PWY <-  fread("PWY_results/GRN_PWY_enrichment.txt")
eGRN_PWY <-  fread("PWY_results/eGRN_PWY_enrichment.txt")
CEN_PWY <-  fread("PWY_results/CEN_PWY_enrichment.txt")
GAN_PWY <-  fread("PWY_results/GRN_PWY_enrichment.txt")

# Add FDR by TF
GRN_PWY %>% 
  #filter(n.targ > 0) %>%
  group_by(TF) %>%
  mutate("FDR" = p.adjust(Pval, method = 'fdr')) %>%
  filter(Pval <= 0.05) %>% 
  arrange(TF) -> GRN_PWY

eGRN_PWY %>% 
  #filter(n.targ > 0) %>%
  group_by(TF) %>%
  mutate("FDR" = p.adjust(Pval, method = 'fdr')) %>%
  filter(Pval <= 0.05) %>% 
  arrange(TF) -> eGRN_PWY

CEN_PWY %>% 
  #filter(n.targ > 0) %>%
  group_by(TF) %>%
  mutate("FDR" = p.adjust(Pval, method = 'fdr')) %>%
  filter(Pval <= 0.05) %>% 
  arrange(TF)  -> CEN_PWY

GAN_PWY %>% 
  #filter(n.targ > 0) %>%
  group_by(TF) %>%
  mutate("FDR" = p.adjust(Pval, method = 'fdr')) %>%
  filter(Pval <= 0.05) %>% 
  arrange(TF)  -> GAN_PWY

# TF names
TFdic <- as_tibble(read.table("Data/Annotations/TF_Id_Name.txt", h=T, stringsAsFactors = F))

GRN_PWY[,"TFname"] <- ReplaceName(GRN_PWY$TF)
eGRN_PWY[,"TFname"] <- ReplaceName(eGRN_PWY$TF)
CEN_PWY[,"TFname"] <- ReplaceName(CEN_PWY$TF)
GAN_PWY[,"TFname"] <- ReplaceName(GAN_PWY$TF)


####################################################################


####################################################################
########       Enrichment of GOs in targets by TF       ############
####################################################################

## test GOS in CEN network
CENdone <- list.files(path = 'BP_results/', pattern = "^GOsBP_CEN_*")
CENdone <- gsub(".txt","", gsub("GOsBP_CEN_", "", CENdone))
length(CENdone)

GenesList <- unique(CoExp$TF)
GenesList <- GenesList[!(GenesList %in% CENdone)]
Lgenes <- length(GenesList)

w=5  # Size of range to test
print(".. Ready to start ..")
GOsDB_CEN <- list()

for (i in seq(0, Lgenes, w)){
  max=Lgenes
  Start=i+1
  end=i+w
  
  if (end<max){
    listtotest <- GenesList[Start:end]
    cat(" working on:", Start:end, "\n")
    GOsDB_CEN <- c(GOsDB_CEN, mclapply(listtotest, function(x) SuperGO(x, "CEN"), mc.cores=w))
  }
  else{
    r <- paste(Start, max, sep = "-")
    listtotest <- GenesList[Start:max]
    cat(" working on:", Start:max, "\n")
    w <- max - (Start-1)
    GOsDB_CEN <- c(GOsDB_CEN, mclapply(listtotest, function(x) SuperGO(x, "CEN"), mc.cores=w))
  }
}

length(GOsDB_CEN)

  
### test GOS in GAN network
GANdone <- list.files(path = 'BP_results/', pattern = "^GOsBP_GAN_*")
GANdone <- gsub(".txt","", gsub("GOsBP_GAN_", "", GANdone))

GenesList <- unique(subset(teQTL, TF %in% All_TFs)$TF)
GenesList <- GenesList[!(GenesList %in% GANdone)]
Lgenes <- length(GenesList)

w=50  # Size of range to test
print(".. Ready to start ..")
GOsDB_GAN <- list()

for (i in seq(0, Lgenes, w)){
  max=Lgenes
  Start=i+1
  end=i+w
  
  if (end<max){
    listtotest <- GenesList[Start:end]
    cat(" working on:", Start:end, "\n")
    GOsDB_GAN <- c(GOsDB_GAN, mclapply(listtotest, function(x) SuperGO(x, "GAN"), mc.cores=w))
  }
  else{
    r <- paste(Start, max, sep = "-")
    listtotest <- GenesList[Start:max]
    cat(" working on:", Start:max, "\n")
    w <- max - (Start-1)
    GOsDB_GAN <- c(GOsDB_GAN, mclapply(listtotest, function(x) SuperGO(x, "GAN"), mc.cores=w))
  }
}

#mask <- unlist(lapply(GOsDB_GAN, function(x) is.data.frame(x)))
#files2read <- paste0("BP_results/", list.files(path = "BP_results/", pattern = "^GOsBP_GAN_*"))
#GOsDB_GAN <- lapply(files2read, fread)

### Test GOS in GRN network
GRNdone <- list.files(path = 'BP_results/', pattern = "^GOsBP_GRN_*")
GRNdone <- gsub(".txt","", gsub("GOsBP_GRN_", "", GRNdone))

GenesList <- unique(PDI$TF)
Lgenes <- length(GenesList)
GenesList <- GenesList[!(GenesList %in% GRNdone)]
Lgenes <- length(GenesList)


w=20  # Size of range to test
print(".. Ready to start ..")
GOsDB_GRN <- list()
for (i in seq(0, Lgenes, w)){
  max=Lgenes
  Start=i+1
  end=i+w
  
  if (end<max){
    listtotest <- GenesList[Start:end]
    cat(paste0(" .. ", Start,":",end," .."))
    GOsDB_GRN <- c(GOsDB_GRN, mclapply(listtotest, function(x) SuperGO(x, "GRN"), mc.cores=w))
  }
  else{
    r <- paste(Start, max, sep = "-")
    listtotest <- GenesList[Start:max]
    w <- max - (Start-1)
    cat(paste0(" .. ", Start,":",max," .."))
    GOsDB_GRN <- c(GOsDB_GRN, mclapply(listtotest, function(x) SuperGO(x, "GRN"), mc.cores=w))
  }
}

#mask <- unlist(lapply(GOsDB_GRN, function(x) is.data.frame(x)))
#GOsDB_GRN <- GOsDB_GRN[mask]

### Test GOS in eGRN network
eGRNdone <- list.files(path = 'BP_results/', pattern = "^GOsBP_eGRN_*")
eGRNdone <- gsub(".txt","", gsub("GOsBP_eGRN_", "", eGRNdone))
GenesList <- unique(PDIeQTL$TF)
GenesList <- GenesList[!(GenesList %in% eGRNdone)]
Lgenes <- length(GenesList)

w=30  # Size of range to test
print(".. Ready to start ..")
GOsDB_eGRN <- list()
for (i in seq(0, Lgenes, w)){
  max=Lgenes
  Start=i+1
  end=i+w
  
  if (end<max){
    listtotest <- GenesList[Start:end]
    cat(paste0(" .. ", Start,":",end," .."))
    
    GOsDB_eGRN <- c(GOsDB_eGRN, mclapply(listtotest, function(x) SuperGO(x, "eGRN"), mc.cores=w))
  }
  else{
    r <- paste(Start, max, sep = "-")
    listtotest <- GenesList[Start:max]
    w <- max - (Start-1)
    cat(paste0(" .. ", Start,":",max," .."))
    GOsDB_eGRN <- c(GOsDB_eGRN, mclapply(listtotest, function(x) SuperGO(x, "eGRN"), mc.cores=w))
  }
}

#mask <- unlist(lapply(GOsDB_eGRN, function(x) is.data.frame(x)))
#GOsDB_eGRN <- GOsDB_eGRN[mask]

#mask <- unlist(lapply(GOsDB_CEN, function(x) is.data.frame(x)))
#GOsDB_CEN <- GOsDB_CEN[mask]

####
GOsDB_GRN <- list.files(path = 'BP_results/', pattern = "^GOsBP_GRN_*")
GOsDB_eGRN <- list.files(path = 'BP_results/', pattern = "^GOsBP_eGRN_*")
GOsDB_CEN <- list.files(path = 'BP_results/', pattern = "^GOsBP_CEN_*")
GOsDB_GAN <- list.files(path = 'BP_results/', pattern = "^GOsBP_GAN_*")

length(GOsDB_GRN)
length(GOsDB_eGRN)
length(GOsDB_CEN)
length(GOsDB_GAN)

# read files
GOsDB_GRN <- lapply(GOsDB_GRN, function(x) fread(paste0("BP_results/",x)))
GOsDB_eGRN <- lapply(GOsDB_eGRN, function(x) fread(paste0("BP_results/",x)))
GOsDB_CEN <- lapply(GOsDB_CEN, function(x) fread(paste0("BP_results/",x)))
GOsDB_GAN <- lapply(GOsDB_GAN, function(x) fread(paste0("BP_results/",x)))

## Combine DF results 
GOsDB_GRN <- rbindlist(GOsDB_GRN, idcol = F)
GOsDB_eGRN <- rbindlist(GOsDB_eGRN, idcol = F)
GOsDB_CEN <- rbindlist(GOsDB_CEN, idcol = F)
GOsDB_GAN <- rbindlist(GOsDB_GAN, idcol = F)

##  add FDR and filter
GOsDB_GRN <- GOsDB_GRN %>% 
  #filter(Significant >= 1) %>% # Remove GOs reported without targets
  group_by(Mutant) %>%
  mutate("FDR" = p.adjust(classic, method = 'fdr')) %>%
  filter(FDR <= 0.1) 

GOsDB_eGRN <- GOsDB_eGRN %>% 
  #filter(Significant >= 1) %>% # Remove GOs reported without targets
  group_by(Mutant) %>%
  mutate("FDR" = p.adjust(classic, method = 'fdr')) %>%
  filter(FDR <= 0.1) 

GOsDB_CEN <- GOsDB_CEN %>% 
  #filter(Significant >= 1) %>% # Remove GOs reported without targets
  group_by(Mutant) %>%
  mutate("FDR" = p.adjust(classic, method = 'fdr')) %>%
  filter(FDR <= 0.1) 

GOsDB_GAN <- GOsDB_GAN %>% 
  #filter(Significant >= 1) %>% # Remove GOs reported without targets
  group_by(Mutant) %>%
  mutate("FDR" = p.adjust(classic, method = 'fdr')) %>%
  filter(FDR <= 0.1) 

####################################################################

# https://genomevolution.org/r/1nkew
####################################################################
######       TFs with multi-network annotation           ###########
####################################################################

# Annotate genes by network
TFs_with_PWYs <- list(GRN =  unique(GRN_PWY$TF), 
                      eGRN = unique(eGRN_PWY$TF),
                      CEN = unique(CEN_PWY$TF),
                      GAN = unique(GAN_PWY$TF))

TFs_with_GOs <- list(GRN = unique(GOsDB_GRN$Mutant), 
                      eGRN = unique(GOsDB_eGRN$Mutant),
                      CEN =  unique(GOsDB_CEN$Mutant),
                      GAN =  unique(GOsDB_GAN$Mutant))


### Number of TFs annotated
TFsAnnotated <- rbind(data.frame(Network=names(unlist(lapply(TFs_with_PWYs, length))),
                           TFs=as.numeric(unlist(lapply(TFs_with_PWYs, length))),
                           Annotation="PWYs"),
                      data.frame(Network=names(unlist(lapply(TFs_with_GOs, length))),
                                 TFs=as.numeric(unlist(lapply(TFs_with_GOs, length))),
                                 Annotation="GOs")) %>% as_tibble()


# List of TFs with predictions in multiple networks
## PWYs
TFs_with_PWYs_DB <- venn(TFs_with_PWYs)
TFs_with_PWYs_DB <- as.list(attr(TFs_with_PWYs_DB, "intersections"))
TFs_with_PWYs_DB <- as.data.table(plyr::ldply(TFs_with_PWYs_DB, data.table))
TFs_with_PWYs_DB <- subset(TFs_with_PWYs_DB, !(.id %in% c("CEN", "eGRN", "GRN", "GAN")))
TFs_with_PWYs_DB <- unique(TFs_with_PWYs_DB)
8+8+830+8+51+49+1+9+1+1

## GOs
TFs_with_GOs_DB <- venn(TFs_with_GOs)
TFs_with_GOs_DB <- as.list(attr(TFs_with_GOs_DB, "intersections"))
TFs_with_GOs_DB <- as.data.table(plyr::ldply(TFs_with_GOs_DB, data.table))
TFs_with_GOs_DB <- subset(TFs_with_GOs_DB, !(.id %in% c("CEN", "eGRN", "GRN", "GAN")))
TFs_with_GOs_DB

####################################################################


####################################################################
########            PWYs similarities                   ############
####################################################################
source("Source_Fig_Function_PWYs_GOs.v2.R")

# 1. Calculate all PWYs vs PWYs  
# Defined format
CornCYC_net <- CornCYC[,c(1,3)]
colnames(CornCYC_net) <- c("TF", "Target")

# Calculate common genes and its significance
PWY_similarity <- Enrichmet_classes(CornCYC_net)
PWY_similarity <- PWY_similarity[,-c(5)]
colnames(PWY_similarity)[c(1, 2, 4)] <- c('PWY1', 'PWY2', 'CommonGenes')
PWY_similarity

# Keep Significant overlapping 
# PWY_similarity <- subset(PWY_similarity, CommonGenes >= 1)
# PWY_similarity <- subset(PWY_similarity, Pval <= 0.01)

# 2 
## Defined TF pairs to compare by networks

TFs_in_CEN_GAN_PWY  <- intersect(CEN_PWY$TF, GAN_PWY$TF)
TFs_in_CEN_eGRN_PWY <- intersect(CEN_PWY$TF, eGRN_PWY$TF)
TFs_in_CEN_GRN_PWY <- intersect(CEN_PWY$TF, GRN_PWY$TF)
TFs_in_GAN_eGRN_PWY <- intersect(GAN_PWY$TF, eGRN_PWY$TF)
TFs_in_GAN_GRN_PWY <- intersect(GAN_PWY$TF, GRN_PWY$TF)
TFs_in_GRN_eGRN_PWY <- intersect(GRN_PWY$TF, eGRN_PWY$TF)


# 3
tf <- TFs_in_CEN_GAN_PWY[1]
## compare TFs by network pairs
# CEN and GAN
CEN_GAN_PWY <- lapply(TFs_in_CEN_GAN_PWY, Func.CEN_GAN_PWY)
names(CEN_GAN_PWY) <- TFs_in_CEN_GAN_PWY
mask <- unlist(lapply(CEN_GAN_PWY, function(x) is.data.frame(x)))
CEN_GAN_PWY <- CEN_GAN_PWY[mask]
CEN_GAN_PWY <- rbindlist(CEN_GAN_PWY, idcol = T)

# CEN and eGRN
CEN_eGRN_PWY <- lapply(TFs_in_CEN_eGRN_PWY, Func.CEN_eGRN_PWY)
names(CEN_eGRN_PWY) <- TFs_in_CEN_eGRN_PWY
mask <- unlist(lapply(CEN_eGRN_PWY, function(x) is.data.frame(x)))
CEN_eGRN_PWY <- CEN_eGRN_PWY[mask]
CEN_eGRN_PWY <- rbindlist(CEN_eGRN_PWY, idcol = T)

# CEN and GRN
CEN_GRN_PWY <- lapply(TFs_in_CEN_GRN_PWY, Func.CEN_GRN_PWY)
names(CEN_GRN_PWY) <- TFs_in_CEN_GRN_PWY
mask <- unlist(lapply(CEN_GRN_PWY, function(x) is.data.frame(x)))
CEN_GRN_PWY <- CEN_GRN_PWY[mask]
CEN_GRN_PWY <- rbindlist(CEN_GRN_PWY, idcol = T)

# GAN and eGRN
GAN_eGRN_PWY <- lapply(TFs_in_GAN_eGRN_PWY, Func.GAN_eGRN_PWY)
names(GAN_eGRN_PWY) <- TFs_in_GAN_eGRN_PWY
mask <- unlist(lapply(GAN_eGRN_PWY, function(x) is.data.frame(x)))
GAN_eGRN_PWY <- GAN_eGRN_PWY[mask]
GAN_eGRN_PWY <- rbindlist(GAN_eGRN_PWY, idcol = T)

# GAN and GRN
GAN_GRN_PWY <- lapply(TFs_in_GAN_GRN_PWY, Func.GAN_GRN_PWY)
names(GAN_GRN_PWY) <- TFs_in_GAN_GRN_PWY
mask <- unlist(lapply(GAN_GRN_PWY, function(x) is.data.frame(x)))
GAN_GRN_PWY <- GAN_GRN_PWY[mask]
GAN_GRN_PWY <- rbindlist(GAN_GRN_PWY, idcol = T)

# GRN and GRN
GRN_eGRN_PWY <- lapply(TFs_in_GRN_eGRN_PWY, Func.GRN_eGRN_PWY)
names(GRN_eGRN_PWY) <- TFs_in_GRN_eGRN_PWY
mask <- unlist(lapply(GRN_eGRN_PWY, function(x) is.data.frame(x)))
GRN_eGRN_PWY <- GRN_eGRN_PWY[mask]
GRN_eGRN_PWY <- rbindlist(GRN_eGRN_PWY, idcol = T)

colnames(CEN_GAN_PWY) <- c("TF", "PWY1", "PWY2", "Pval", "PWYsCommonGenes", "PWY1_size", "PWY2_size", "TargPWY1", "TargPWY2", "PWY1name", "PWY2name")
colnames(CEN_eGRN_PWY) <- c("TF", "PWY1", "PWY2", "Pval", "PWYsCommonGenes", "PWY1_size", "PWY2_size", "TargPWY1", "TargPWY2", "PWY1name", "PWY2name")
colnames(CEN_GRN_PWY) <- c("TF", "PWY1", "PWY2", "Pval", "PWYsCommonGenes", "PWY1_size", "PWY2_size", "TargPWY1", "TargPWY2", "PWY1name", "PWY2name")
colnames(GAN_eGRN_PWY) <- c("TF", "PWY1", "PWY2", "Pval", "PWYsCommonGenes", "PWY1_size", "PWY2_size", "TargPWY1", "TargPWY2", "PWY1name", "PWY2name")
colnames(GAN_GRN_PWY) <- c("TF", "PWY1", "PWY2", "Pval", "PWYsCommonGenes", "PWY1_size", "PWY2_size", "TargPWY1", "TargPWY2", "PWY1name", "PWY2name")
colnames(GRN_eGRN_PWY) <- c("TF", "PWY1", "PWY2", "Pval", "PWYsCommonGenes", "PWY1_size", "PWY2_size", "TargPWY1", "TargPWY2", "PWY1name", "PWY2name")

CEN_GAN_PWY[,"PWYs"] <- "CEN:GAN"
CEN_eGRN_PWY[,"PWYs"] <- "CEN:eGRN"
CEN_GRN_PWY[,"PWYs"] <- "CEN:GRN"
GAN_eGRN_PWY[,"PWYs"] <- "GAN:eGRN"
GAN_GRN_PWY[,"PWYs"] <- "GAN:GRN"
GRN_eGRN_PWY[,"PWYs"] <- "GRN:eGRN"

PWY_Enriched <- rbind(CEN_GAN_PWY, 
                      CEN_eGRN_PWY,
                      CEN_GRN_PWY,
                      GAN_eGRN_PWY, 
                      GAN_GRN_PWY,
                      GRN_eGRN_PWY)


write.table(PWY_Enriched, 
            "PWY_GO_results/CommonFunction_PWY_enrichment.txt", 
            row.names = F, quote = F, sep = '\t')


####################################################################


####################################################################
########            GOs similarities                        ########
####################################################################
# Pre-calculate semantic similarity
Zm.GOSemSim.BP <- GOSemSim::godata(org.Zmays.eg.db, keytype="GENENAME", ont = 'BP')
typeof(org.Zmays.eg.db)


# list of TFs enriched with at least a GO in each network

############
## Semantic similarity background using AnnotationHub
# library(AnnotationHub)
# ah = AnnotationHub()
# zm <- query(ah, c("Zea mays"))
# zm <- zm["AH96032"]
# zm_hub <- ah[[zm$ah_id]] # Extract info for specific maize id
# # Create Semantic similarity calculations: in this cases based on BP
# zmGO_BP <- godata(zm_hub, ont="BP")
############

# Semantic similarity for TFs with GOs in at least two networks
as.data.table(table(TFs_with_GOs_DB$.id))
source("Source_Fig_Function_PWYs_GOs.v2.R")

TFs_in_GRN_CEN_GO <- unique(intersect(GOsDB_GRN$Mutant, GOsDB_CEN$Mutant))
TFs_in_GRN_GAN_GO <- unique(intersect(GOsDB_GRN$Mutant, GOsDB_GAN$Mutant))
TFs_in_GRN_eGRN_GO<- unique(intersect(GOsDB_GRN$Mutant, GOsDB_eGRN$Mutant))
TFs_in_CEN_GAN_GO <- unique(intersect(GOsDB_CEN$Mutant, GOsDB_GAN$Mutant))
TFs_in_CEN_eGRN_GO<- unique(intersect(GOsDB_CEN$Mutant, GOsDB_eGRN$Mutant)) 
TFs_in_GAN_eGRN_GO<- unique(intersect(GOsDB_GAN$Mutant, GOsDB_eGRN$Mutant))



#######################

## Pair Comparisons: without reducing GOs to parents ##
###
# SS_GRN_CEN <- Get_SS_GRN_CEN(TFs_in_GRN_CEN_GO)
###

w=40  # Size of range to test
print(".. Ready to start ..")
SS_GRN_CEN <- list()
GenesList <- TFs_in_GRN_CEN_GO
Lgenes <- length(GenesList)

for (i in seq(0, Lgenes, w)){
  max=Lgenes
  Start=i+1
  end=i+w
  
  if (end<max){
    listtotest <- GenesList[Start:end]
    cat(paste0(" .. ", Start,":",end," .."))
    SS_GRN_CEN <- c(SS_GRN_CEN, mclapply(listtotest, function(x) Get_SS_GRN_CEN(x), mc.cores=w))
  }
  else{
    r <- paste(Start, max, sep = "-")
    listtotest <- GenesList[Start:max]
    w <- max - (Start-1)
    cat(paste0(" .. ", Start,":",max," .."))
    SS_GRN_CEN <- c(SS_GRN_CEN, mclapply(listtotest, function(x) Get_SS_GRN_CEN(x), mc.cores=w))
    
  }
}

SS_GRN_CEN <- rbindlist(SS_GRN_CEN, idcol = F)

###
# SS_GRN_GAN <- Get_SS_GRN_GAN(TFs_in_GRN_GAN_GO)
###
w=5  # Size of range to test
print(".. Ready to start ..")
SS_GRN_GAN <- list()
GenesList <- TFs_in_GRN_GAN_GO
Lgenes <- length(GenesList)

for (i in seq(0, Lgenes, w)){
  max=Lgenes
  Start=i+1
  end=i+w
  
  if (end<max){
    listtotest <- GenesList[Start:end]
    cat(paste0(" .. ", Start,":",end," .."))
    SS_GRN_GAN <- c(SS_GRN_GAN, mclapply(listtotest, function(x) Get_SS_GRN_GAN(x), mc.cores=w))
  }
  else{
    r <- paste(Start, max, sep = "-")
    listtotest <- GenesList[Start:max]
    w <- max - (Start-1)
    cat(paste0(" .. ", Start,":",max," .."))
    SS_GRN_GAN <- c(SS_GRN_GAN, mclapply(listtotest, function(x) Get_SS_GRN_GAN(x), mc.cores=w))
    
  }
}

SS_GRN_GAN <- rbindlist(SS_GRN_GAN, idcol = F)

###
# SS_GRN_eGRN <- Get_SS_GRN_eGRN(TFs_in_GRN_eGRN_GO)
###

w=40  # Size of range to test
print(".. Ready to start ..")
SS_GRN_eGRN <- list()
GenesList <- TFs_in_GRN_eGRN_GO
Lgenes <- length(GenesList)

for (i in seq(0, Lgenes, w)){
  max=Lgenes
  Start=i+1
  end=i+w
  
  if (end<max){
    listtotest <- GenesList[Start:end]
    cat(paste0(" .. ", Start,":",end," .."))
    SS_GRN_eGRN <- c(SS_GRN_eGRN, mclapply(listtotest, function(x) Get_SS_GRN_eGRN(x), mc.cores=w))
  }
  else{
    r <- paste(Start, max, sep = "-")
    listtotest <- GenesList[Start:max]
    w <- max - (Start-1)
    cat(paste0(" .. ", Start,":",max," .."))
    SS_GRN_eGRN <- c(SS_GRN_eGRN, mclapply(listtotest, function(x) Get_SS_GRN_eGRN(x), mc.cores=w))
    
  }
}

SS_GRN_eGRN <- rbindlist(SS_GRN_eGRN, idcol = F)

###
# SS_CEN_GAN <- Get_SS_CEN_GAN(TFs_in_CEN_GAN_GO)
###

w=40  # Size of range to test
print(".. Ready to start ..")
SS_CEN_GAN <- list()
GenesList <- TFs_in_CEN_GAN_GO
Lgenes <- length(GenesList)

for (i in seq(0, Lgenes, w)){
  max=Lgenes
  Start=i+1
  end=i+w
  
  if (end<max){
    listtotest <- GenesList[Start:end]
    cat(paste0(" .. ", Start,":",end," .."))
    SS_CEN_GAN <- c(SS_CEN_GAN, mclapply(listtotest, function(x) Get_SS_CEN_GAN(x), mc.cores=w))
  }
  else{
    r <- paste(Start, max, sep = "-")
    listtotest <- GenesList[Start:max]
    w <- max - (Start-1)
    cat(paste0(" .. ", Start,":",max," .."))
    SS_CEN_GAN <- c(SS_CEN_GAN, mclapply(listtotest, function(x) Get_SS_CEN_GAN(x), mc.cores=w))
    
  }
}

SS_CEN_GAN <- rbindlist(SS_CEN_GAN, idcol = F)
SS_CEN_GAN[(SS_CEN_GAN$SS >= 0.5),]

###
# SS_eGRN_CEN <- Get_SS_eGRN_CEN(TFs_in_CEN_eGRN_GO)
###

w=40  # Size of range to test
print(".. Ready to start ..")
SS_eGRN_CEN <- list()
GenesList <- TFs_in_CEN_eGRN_GO
Lgenes <- length(GenesList)

for (i in seq(0, Lgenes, w)){
  max=Lgenes
  Start=i+1
  end=i+w
  
  if (end<max){
    listtotest <- GenesList[Start:end]
    cat(paste0(" .. ", Start,":",end," .."))
    SS_eGRN_CEN <- c(SS_eGRN_CEN, mclapply(listtotest, function(x) Get_SS_eGRN_CEN(x), mc.cores=w))
  }
  else{
    r <- paste(Start, max, sep = "-")
    listtotest <- GenesList[Start:max]
    w <- max - (Start-1)
    cat(paste0(" .. ", Start,":",max," .."))
    SS_eGRN_CEN <- c(SS_eGRN_CEN, mclapply(listtotest, function(x) Get_SS_eGRN_CEN(x), mc.cores=w))
    
  }
}

SS_eGRN_CEN <- rbindlist(SS_eGRN_CEN, idcol = F)

###
SS_GAN_eGRN <- Get_SS_GAN_eGRN(TFs_in_GAN_eGRN_GO)
###
#######################

#######################
#### Global GSS analysis: Single value per pair of networks
#######################
SS_GRN_CEN
SS_GRN_GAN
SS_GRN_eGRN
SS_CEN_GAN
SS_eGRN_CEN
SS_GAN_eGRN

DF_SS <- rbind(SS_GRN_CEN, SS_GRN_GAN, SS_GRN_eGRN, 
               SS_CEN_GAN, SS_eGRN_CEN, SS_GAN_eGRN)

DF_SS[,"TF"] <- ReplaceName(DF_SS$Source)
DF_SS$SS >= 0.6

source("Source_Fig_Function_PWYs_GOs.v2.R")
#######################

#######################
## Read individual values: TFs with at Sig. GOs in at least tow nets
#######################

# Read GSS for individual GOs by TF
gossfiles <- list.files(path = 'GO_SS_data/', pattern = 'GSS_*')

GO_ssTF_DB <- lapply(gossfiles, Read_GO_ssTF)
GO_ssTF_DB <- rbindlist(GO_ssTF_DB, idcol = F)
GO_ssTF_DB[,"TFname"] <- ReplaceName(GO_ssTF_DB$TF)

GO_ssTF_DB <- left_join(GO_ssTF_DB, TFdic, by=c("TF"="TF.v4"))
GO_ssTF_DB$TF.Name[is.na(GO_ssTF_DB$TF.Name)] <- GO_ssTF_DB$TF[is.na(GO_ssTF_DB$TF.Name)]
colnames(GO_ssTF_DB)[12] <- "TFname"

# TFs order in plot
GO_ssTF_DB %>%
  dplyr::group_by(TFname) %>%
  dplyr::summarise(SS=mean(GSS), SigGOs=sum((GSS>=0.6)*1)) %>%
  arrange(-SS) -> GSS_order


# defined TFs named with number included
GSS_order[,"TFname2"] <- paste0(GSS_order$TFname, ' [', GSS_order$SigGOs,']')

# add name to main GSS DF
GO_ssTF_DB <- left_join(GO_ssTF_DB, GSS_order[,c(1,4)], by="TFname")


# Filter significant comparisons
GO_ssTF_DB <- subset(GO_ssTF_DB, GSS >= 0.6)
#######################

write.table(GO_ssTF_DB, "PWY_GO_results/CommonFunction_GO_enrichment.txt", row.names = F, quote = F, sep = '\t')

#######################
# Summary number of GOs tested by TFs in GSS_DB
#######################

#########
# Count GO terms tested
#########

## Part 1: 
GO_ssTF_DB %>%
  dplyr::select(GO1, GO2, TF, Nets) %>%
  dplyr::mutate(Net.1 = chop(Nets, '[:]',1), Net.2=chop(Nets, '[:]',2)) %>%
  dplyr::select(GO1, TF, Net.1) %>%
  dplyr::group_by(TF, Net.1) %>%
  dplyr::summarise(n=length(unique(GO1))) -> GOs_in_GSS_1
colnames(GOs_in_GSS_1)[2] <- 'Net'

## Part 2
GO_ssTF_DB %>%
  dplyr::select(GO1, GO2, TF, Nets) %>%
  dplyr::mutate(Net.1 = chop(Nets, '[:]',1), Net.2=chop(Nets, '[:]',2)) %>%
  dplyr::select(GO2, TF, Net.2) %>%
  dplyr::group_by(TF, Net.2) %>%
  dplyr::summarise(n=length(unique(GO2))) -> GOs_in_GSS_2
colnames(GOs_in_GSS_2)[2] <- 'Net'

## Combine part 1 and 2
rbind(GOs_in_GSS_1, GOs_in_GSS_2) %>%
  dplyr::arrange(TF) -> GOs_in_GSS

GOs_in_GSS %>% 
  dplyr::group_by(Net) %>%
  dplyr::summarise(nmean=mean(n))
#########  

#########
## Count GO term that pass the GSS filter
#########

## Part 1
GO_ssTF_DB %>%
  dplyr::filter(GSS >=  0.6) %>%
  dplyr::select(GO1, GO2, TF, Nets) %>%
  dplyr::mutate(Net.1 = chop(Nets, '[:]',1), Net.2=chop(Nets, '[:]',2)) %>%
  dplyr::select(GO1, TF, Net.1) %>%
  dplyr::group_by(TF, Net.1) %>%
  dplyr::summarise(n=length(unique(GO1))) -> Top_GOs_in_GSS_1

colnames(Top_GOs_in_GSS_1)[2] <- 'Net'

## Part 2
GO_ssTF_DB %>%
  dplyr::filter(GSS >=  0.6) %>%
  dplyr::select(GO1, GO2, TF, Nets) %>%
  dplyr::mutate(Net.1 = chop(Nets, '[:]',1), Net.2=chop(Nets, '[:]',2)) %>%
  dplyr::select(GO2, TF, Net.2) %>%
  dplyr::group_by(TF, Net.2) %>%
  dplyr::summarise(n=length(unique(GO2))) -> Top_GOs_in_GSS_2

colnames(Top_GOs_in_GSS_2)[2] <- 'Net'

## Combine part 1 and 2
rbind(Top_GOs_in_GSS_1, Top_GOs_in_GSS_2) %>%
  dplyr::arrange(TF) -> Top_GOs_in_GSS


#########

GOs_in_GSS      # tested
Top_GOs_in_GSS  # Past
PWY_Enriched

# Combine results from total and only significant GO counts by TFa and Net
GOs_in_GSS_Summary <- left_join(GOs_in_GSS, Top_GOs_in_GSS, by=c("TF","Net"))
colnames(GOs_in_GSS_Summary)[c(3,5)] <- c("Total.GOs", "Similar.GOs")

# Replace NAs by zeru
GOs_in_GSS_Summary$Similar.GOs[is.na(GOs_in_GSS_Summary$Similar.GOs)] <- 0

# Add "Similar.GOs" as percentage
GOs_in_GSS_Summary[,"Perc.Similar.GOs"] <- (GOs_in_GSS_Summary$Similar.GOs/GOs_in_GSS_Summary$Total.GOs)*100
#########

##
#
100*length(unique(c(GRN_PWY$TF, GOsDB_GRN$Mutant)))/length(unique(PDI$TF))
100*length(unique(c(eGRN_PWY$TF, GOsDB_eGRN$Mutant)))/length(unique(PDIeQTL$TF))
100*length(unique(c(CEN_PWY$TF, GOsDB_CEN$Mutant)))/length(unique(CoExp$TF))
100*length(unique(c(GAN_PWY$TF, GOsDB_GAN$Mutant)))/length(unique(subset(teQTL, TF %in% All_TFs)$TF))


ReduceGOs <- function(tf, net){
  
  if(net=="GRN"){
    # Defined 
    GOsDB = subset(GOsDB_GRN, Mutant == tf)
    scores <- setNames(-log10(GOsDB$FDR), GOsDB$GO.ID) 
    GO_vector = subset(GOsDB, Mutant == tf)$GO.ID
  }
  else if (net=='eGRN') {
    GOsDB = subset(GOsDB_eGRN, Mutant == tf)
    scores <- setNames(-log10(GOsDB$FDR), GOsDB$GO.ID) 
    GO_vector = subset(GOsDB, Mutant == tf)$GO.ID
  }
  else if (net=='CEN') {
    GOsDB = subset(GOsDB_CEN, Mutant == tf)
    scores <- setNames(-log10(GOsDB$FDR), GOsDB$GO.ID) 
    GO_vector = subset(GOsDB, Mutant == tf)$GO.ID
  }
  else if (net=='GAN') {
    GOsDB = subset(GOsDB_GAN, Mutant == tf)
    scores <- setNames(-log10(GOsDB$FDR), GOsDB$GO.ID) 
    GO_vector = subset(GOsDB, Mutant == tf)$GO.ID
  }
  
  if (length(GO_vector) >1) {
    # Semantic similarity
    simMatrix <- calculateSimMatrix(GO_vector,  orgdb=org.Zmays.eg.db,  ont="BP", 
                                    semdata=Zm.GOSemSim.BP,
                                    method="Wang")
    
    # Reduce term
    reducedTerms <- reduceSimMatrix(simMatrix, scores, keytype="GENENAME",
                                    threshold=0.7, orgdb=org.Zmays.eg.db)
    
    # List of parents
    parent <- unique(reducedTerms[,c("parent", "parentTerm")])
    
    ReduceDB <- subset(GOsDB, Mutant == tf & GO.ID %in% parent$parent)
    ReduceDB <- left_join(ReduceDB, parent, by=c("GO.ID"="parent"))
    ReduceDB$Term <- ReduceDB$parentTerm
    ReduceDB <- ReduceDB[,-c(9)]
    
    return(ReduceDB)
    
  }
  
  return(subset(GOsDB, Mutant == tf))
}
# map GOs to parent and reduce GO df
GOsDB_GRN_reduced <- lapply(unique(GOsDB_GRN$Mutant),   function(x) ReduceGOs(x,"GRN"))
GOsDB_eGRN_reduced <- lapply(unique(GOsDB_eGRN$Mutant), function(x) ReduceGOs(x,"eGRN"))
GOsDB_CEN_reduced <- lapply(unique(GOsDB_CEN$Mutant),   function(x) ReduceGOs(x,"CEN"))
GOsDB_GAN_reduced <- lapply(unique(GOsDB_GAN$Mutant),   function(x) ReduceGOs(x,"GAN"))


GOsDB_GRN_reduced <-  rbindlist(GOsDB_GRN_reduced, idcol = F)
GOsDB_eGRN_reduced <- rbindlist(GOsDB_eGRN_reduced, idcol = F)
GOsDB_CEN_reduced <-  rbindlist(GOsDB_CEN_reduced, idcol = F)
GOsDB_GAN_reduced <-  rbindlist(GOsDB_GAN_reduced, idcol = F,  fill = T)


##
########################
## Plots
########################
source("Source_Fig_Function_PWYs_GOs.v2.R")

### Part A
# Total TFs annotated by PWYs and GOs
TFsAnnotated %>%
  ggplot(aes(y=Network, x=TFs, fill=Annotation)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_pubclean() +
  geom_text(aes(label=scales::comma(TFs), x=(TFs - TFs*0.25)),
            position = position_dodge(0.9), size=3) +
  scale_x_log10(expand=c(0,0), limits = c(1, 4000), label=scales::comma) +
  annotation_logticks(sides = "b", color = 'black') +
  theme(axis.text=element_text(size=10), legend.position = 'bottom',
        legend.text=element_text(size=9),
        text = element_text(size=10, family="Times"))  -> Plo_TFsAnnotated


Plo_TFsAnnotated
### Part B
# TFs with PWYs enriched
Plot_TFs_venn_PWYs <- vennfuncInt(TFs_with_PWYs)

### Part C
# TFs with GOs enriched
Plot_TFs_venn_GOs <- vennfuncInt(TFs_with_GOs)

## Combine  b, c, and d plots
Plot_bcd <- {Plo_TFsAnnotated/Plot_TFs_venn_PWYs/Plot_TFs_venn_GOs}  +  
            plot_layout(ncol = 1)
Plot_bcd

pdf("Plots/Plot_S5bcd.pdf", width = 4, height = 6)
print(Plot_bcd)
dev.off()

length(TFs_in_CEN_GAN_PWY)
length(TFs_in_GRN_eGRN_PWY)
length(TFs_in_CEN_eGRN_PWY)
length(TFs_in_CEN_GRN_PWY)
length(TFs_in_GAN_eGRN_PWY)
length(TFs_in_GAN_GRN_PWY)

length(TFs_in_CEN_GAN_GO)
length(TFs_in_GRN_CEN_GO)
length(TFs_in_CEN_eGRN_GO)
length(TFs_in_GRN_eGRN_GO)
length(TFs_in_GRN_GAN_GO)
length(TFs_in_GAN_eGRN_GO)
########################


### Part D
colorGroups <- c(CEN = '#FFD700', GRN='#6A5ACD', GAN='#1E90FF', eGRN='#FF1493')

############
GOs_in_GSS[,"TFname"] <- ReplaceName(GOs_in_GSS$TF)
GOs_in_GSS$TFname <- factor(GOs_in_GSS$TFname, levels = GSS_order$TFname)

GOs_in_GSS %>% 
  ggplot(aes(y=TFname, x=n, fill=Net)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_x_log10(expand=c(0,0), limits = c(-1, 600), label=Scales) +
  annotation_logticks(sides = "b", color = 'black') +
  scale_fill_manual(values = alpha(colorGroups, 0.6)) +
  xlab(bquote("GO term enriched (FDR "<="0.1)")) + 
  ylab("TFs") + 
  theme(strip.text.x = element_text(size = 10), 
        strip.text.y = element_text(size = 8), 
        axis.text=element_text(size=10),
        text = element_text(size=10, family="Helvetica"),
        legend.position = 'bottom', legend.direction = 'horizontal') + 
  guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
  labs(fill='Network') -> Plot_TotalGOs
############

# GSS from combined all GOs by TF
##### 
ggplot(DF_SS, aes(y=Net, x=SS, color=Net2)) +
  geom_jitter(height = 0.1) +
  theme_pubclean() +
  xlab("TF GO Semantic similarity") + 
  ylab("") + 
  geom_vline(xintercept = 0.6, linetype="dashed") +
  theme(strip.text.x = element_text(size = 12), 
        axis.text=element_text(size=12),
        text = element_text(size=12, family="Helvetica"),
        legend.position = 'right', legend.direction="vertical") +
  labs(color='Net. Compared') -> Plot_GOSS
##### 

### Part E
# GSS from individual GOs by TF
GO_ssTF_DB$TFname2 <- factor(GO_ssTF_DB$TFname2, levels = GSS_order$TFname2)
write.table(subset(GO_ssTF_DB, GSS >= 0.6), "MainResults/TFs_and_GOpairs_02012023.txt", sep = "\t", quote = F, row.names = F)

ggplot(GO_ssTF_DB, aes(y=TFname2, x=GSS, fill=Nets)) +
  #geom_jitter(height = 0.1) +
  geom_boxplot(notch = T, outlier.size = 0.1) + 
  theme_pubclean() +
  xlab("GO Semantic similarity") + 
  ylab("TFs") + 
  geom_vline(xintercept = 0.6, linetype="dashed") +
  scale_x_continuous(expand = c(0,0)) +
  theme(strip.text.x = element_text(size = 10), 
        strip.text.y = element_text(size = 8), 
        axis.text=element_text(size=10),
        text = element_text(size=10, family="Helvetica"),
        legend.position = 'bottom', legend.direction = 'horizontal') + 
  guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
  labs(fill='Net. Compared') -> Plot_GOSS_ind

### Part F
GOs_in_GSS_Summary$TFname <- factor(GOs_in_GSS_Summary$TFname, levels = GSS_order$TFname)
write.table(GOs_in_GSS_Summary, "TFs_with_SigGSS_02012023.txt", sep = "\t", quote = F, row.names = F)

GOs_in_GSS_Summary %>% 
  ggplot(aes(y=TFname, x=Perc.Similar.GOs, fill=Net)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_manual(values = alpha(colorGroups, 0.6)) +
  xlab("GOs similar/GOs annotated") + 
  ylab("TFs") + 
  scale_x_continuous(expand = c(0,0)) +
  theme(strip.text.x = element_text(size = 10), 
        strip.text.y = element_text(size = 8), 
        axis.text=element_text(size=10),
        text = element_text(size=10, family="Helvetica"),
        legend.position = 'bottom', legend.direction = 'horizontal') + 
  guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
  labs(fill='Network') -> Plot_Total_similarGOs

Plot_Total_similarGOs


# Final Plot
Plot_CommonTFs <- {{Plo_TFsAnnotated+Plot_TFs_venn_PWYs+Plot_TFs_venn_GOs +
    plot_layout( widths = c(.6, 1, 1))}/{Plot_TotalGOs+Plot_GOSS_ind+Plot_Total_similarGOs + 
        plot_layout(guides = "collect") & theme(legend.position = 'bottom')}} +
    plot_layout(heights = c(.2, 1)) 
Plot_CommonTFs

tiff("Plots/Plot_CommonTFs_annotated.tiff", units="in", width=11, height=15, res=300)
print(Plot_CommonTFs)
dev.off()

########################
DF_SS
summary(GOs_in_GSS_Summary$Perc.Similar.GOs)

GOs_in_GSS_Summary %>%
  dplyr::group_by(Net) %>%
  dplyr::summarise(mena=mean(Perc.Similar.GOs))


####################################################################

PWY_similarity[,"PWY1name"] <- ReplaceNamePWY(PWY_similarity$TF)
PWY_similarity[,"PWY2n"] <- ReplaceNamePWY(PWY_similarity$PWY)



