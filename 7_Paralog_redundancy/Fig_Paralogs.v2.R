###################################################################################
#######                             Libraries                               #######
###################################################################################
library(patchwork)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(scales)
library(viridis)
library(RColorBrewer)
library(ggrepel)
library(data.table)
library(circlize)
library(factoextra)
library(reshape2)
library(GeneOverlap)
library(fgsea)
library(parallel)
library(scales)
library(hrbrthemes)
library(ggpointdensity)
library(ggVennDiagram)
library(gplots)

set.seed(42)
suppressMessages(library(igraph))

###################################################################################
#######                        Functions                                    #######
###################################################################################

ReplaceName <- function(ids){
  # for (i in 1:nrow(Top45)){
  #   ids <- gsub(Top45$V1[i], Top45$V2[i], ids)
  # }
  
  for (i in 1:nrow(TFdic)){
    ids <- gsub(TFdic$V2[i], TFdic$V1[i], ids)
  }
  return(ids)
}

ReplaceNamePWY <- function(ids){
  
  for (i in 1:nrow(CornCYC)){
    w <- paste0('\\<', CornCYC$Pathway.id[i], '\\>')
    ids <- gsub(w, CornCYC$Pathway.name[i], ids)
    #ids <- gsub("_", " ", ids)
  }
  return(ids)
}

###
RewireNet <- function(net) {
  ####  
  # This function take as df_net-like to convert it into an igraph 
  # object. rewire it, and return a random network.
  ####
  igraph_R <- graph_from_data_frame(net, directed = F) 
  igraph_R <- simplify(igraph_R)
  # Rewire network with similar degree
  igraph_R <- rewire(igraph_R, with = keeping_degseq(loops = FALSE, niter = vcount(igraph_R)*1000))
  # Get  igraph DF
  out <- as.data.table(as_data_frame(igraph_R, what = "edges"))
  
  colnames(out) <- c("GeneID1", "GeneID2")
  return(out)
}

chop=function(myStr,mySep,myField){
  ## chop a string by a separator and return specified field
  choppedString=sapply(strsplit(myStr,mySep),"[",myField)
  if(length(myField)>1){
    choppedString=apply(choppedString,2,function(x){paste0(x[!is.na(x)],collapse=mySep)})
  }
  return(choppedString)
}


jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

##############################################################################
##################         Read data input           #########################
##############################################################################
##
saf <- as_tibble(read.table("Data/eQTL_data/Zea_mays.B73_RefGen_v4.46.saf", stringsAsFactors = F))
saf1 <- subset(saf, V5=="+")[,c(1,2,3)]
saf2 <- subset(saf, V5=="-")[,c(1,2,4)]
colnames(saf1) <- c("GeneID", "chrAnn", "TSS")
colnames(saf2) <- c("GeneID", "chrAnn", "TSS")
#
saf <- rbind(saf1, saf2)


## Syntenic genes 
Syntenic <- as_tibble(read.table("../Data/Annotations/Zm.v4.synteny.genes.txt", h=T, stringsAsFactors = F, sep = '\t'))
Syntenic <- subset(Syntenic, gene_synteny=='syntenic')$gene_id
length(Syntenic)

# TFdb
TFdb <- as_tibble(read.table("Data/Annotations/TF_CoR_Mazie.txt", h=T, stringsAsFactors = F))

# All TFs
AllTFs <- fread("../Fig_RandomNets/All_TFs.txt", header = F)$V1

# TF names
TFdic <- as_tibble(read.table("Data/Annotations/TF_Id_Name.txt", h=F, stringsAsFactors = F))

# TF target net from Network-based
MRMI_full <- fread('../Fig_pecanpy/DistanceCalculation/InputClusterONE_Dim50_WL80_nW10_0.005_syntenic.txt')
MRMI_full  <- subset(MRMI_full, V2 %in% Syntenic)
MRMI_full  <- subset(MRMI_full, V1 %in% Syntenic)

# Full network before MR_mi calculation
FullNet <- fread("../Fig_CommonTarg/Full_Final_network.11022022.txt")



# Read list of TDG
Paralogs <- fread("Data/Annotations/Paralogs_Schable_Lab.csv", h=T, sep = ',')
Paralogs <- subset(Paralogs, maize1_v4 !="No Gene" & maize2_v4 !="No Gene")
Paralogs <- subset(Paralogs, maize1_v4 !="No Gene" & maize2_v4 !="No Gene")

# confirm syntenic genes
Paralogs <- subset(Paralogs, maize1_v4 %in% Syntenic)
Paralogs <- subset(Paralogs, maize2_v4 %in% Syntenic)

# Genes in MR_MI DB
MR_MI_DB_ids <- list.files(path = 'MR_edgesDB_syntenic/', pattern = "^MR_MI.pecanpy*")
MR_MI_DB_ids <- gsub(".txt", "", gsub("MR_MI.pecanpy.", "", MR_MI_DB_ids))

# Confirm genes with embbeding information (in network)
Paralogs <- subset(Paralogs, maize1_v4 %in% MR_MI_DB_ids)
Paralogs <- subset(Paralogs, maize2_v4 %in% MR_MI_DB_ids)

Paralogs <- unique(Paralogs[,2:3])
##############################################################################

##############################################################################
##################            Paralogs_DB            #########################
##############################################################################

# Create db of all possible pairs of paralogs
Paralogs_DB <- c(split(Paralogs$maize1_v4, Paralogs$maize2_v4),
                 split(Paralogs$maize2_v4, Paralogs$maize1_v4))

# mark genes with more than a paralogs
mask <- unlist(lapply(Paralogs_DB, length)) >= 2
table(unlist(lapply(Paralogs_DB, length)))

# Paralogs with a single copy
Paralogs_only_1 <- as.data.frame(Paralogs_DB[unlist(lapply(Paralogs_DB, length)) == 1])

t(Paralogs_only_1) %>%
  as.data.frame() %>%
  rownames_to_column("GeneID1") %>%
  as.data.table() -> Paralogs_only_1

colnames(Paralogs_only_1)[2] <- 'GeneID2'

# Genes with more than 2 paralogs
Paralogs_more_1 <- Paralogs_DB[mask]

# add gene from maize 1 to identify al possible comparisons
for( i in names(Paralogs_more_1)){
  Paralogs_more_1[[i]] <- c(Paralogs_more_1[[i]], i)
}

# Identify all possible combination of two genes for each group of paralogs
Paralogs_more_1 <- rbindlist(lapply(Paralogs_more_1, function(x) as.data.table(t(combn(x, 2)))), 
                             idcol = F)
colnames(Paralogs_more_1) <- c('GeneID1', "GeneID2")


# transform pairs into graph obj to reduce redundancy 
Paralogs_only_1 <- graph_from_data_frame(Paralogs_only_1, directed = F) 
Paralogs_only_1 <- simplify(Paralogs_only_1)

Paralogs_more_1 <- graph_from_data_frame(Paralogs_more_1, directed = F) 
Paralogs_more_1 <- simplify(Paralogs_more_1)

# igraph as DF
Paralogs_only_1 <- as.data.table(as_data_frame(Paralogs_only_1, what = "edges"))
Paralogs_more_1 <- as.data.table(as_data_frame(Paralogs_more_1, what = "edges"))

## 


Paralogs_DB <- rbind(Paralogs_only_1 %>% 
                       dplyr::rename('GeneID1'='to') %>%
                       dplyr::rename('GeneID2'='from'),
                     Paralogs_more_1 %>% 
                       dplyr::rename('GeneID1'='to') %>%
                       dplyr::rename('GeneID2'='from'))

Paralogs_DB <- unique(Paralogs_DB)

# keep only TFs
Paralogs_DB[,'isTF'] <- (Paralogs_DB$GeneID1 %in% AllTFs)*1 + (Paralogs_DB$GeneID2 %in% AllTFs)*1
Paralogs_DB <- subset(Paralogs_DB, isTF > 0)

# save names for pairs of genes
Paralogs_DB <- subset(Paralogs_DB, GeneID1 != GeneID2)
Paralogs_DB <- unique(Paralogs_DB[,1:2])


for (i in 1:nrow(Paralogs_DB)){
  file <- paste0("AAsequnceDistance/",Paralogs_DB$GeneID1[i], "_", Paralogs_DB$GeneID2[i],".txt")
  fwrite(data.table(gid=c(Paralogs_DB$GeneID1[i], Paralogs_DB$GeneID2[i])),
         file,
         row.names = F, col.names = F, quote = F) 
  print(i)
}

##############################################################################

##################################################################
########          MR of MI  from the Embedding            ######## 
##################################################################

Read_m1_MI_files <- function(id, DB){
  # gene id from maize copy 2
  id2 <- unique(DB[DB$GeneID1 == id,]$GeneID2)
  # MR_MI.pecanpy.Zm00001d036331.txt
  file <- paste0('MR_edgesDB_syntenic/MR_MI.pecanpy.',id,'.txt')
  #
  mr.mi <- fread(file)[,1:4]
  mr.mi <- subset(mr.mi, GeneID2 %in% id2)
  return(mr.mi)
}

#
ParalogsPair_MRdb <- lapply(unique(Paralogs_DB$GeneID1), function(x) Read_m1_MI_files(x, Paralogs_DB))
table(unlist(lapply(ParalogsPair_MRdb, nrow)))
#
#mask <- unlist(lapply(ParalogsPairMRDB, nrow)) == 3
#ParalogsPairMRDB[mask]

# Combine results
ParalogsPair_MRdb <- rbindlist(ParalogsPair_MRdb, idcol = F)

## Create random background
DFrandom_MR_DB <- list()

for(i in seq(1000)){
  # Calculate random net of paralogs
  DFrandom <- RewireNet(Paralogs_DB[,1:2])
  
  # calculate MR for random pairs
  DFrandom_MR <- lapply(unique(DFrandom$GeneID1), function(x) Read_m1_MI_files(x, DFrandom))
  DFrandom_MR <- rbindlist(DFrandom_MR, idcol = F)
  DFrandom_MR_DB[[i]] <- DFrandom_MR
}

# random distribution by set of paralogs pair
MeanDis_randomMR <- unlist(lapply(DFrandom_MR_DB, function(x) mean(x$MR)))

# random distribution by gene
MeanDis_randomMR_list <- rbindlist(DFrandom_MR_DB, idcol = T)

MeanDis_randomMR_list <- rbind(MeanDis_randomMR_list[,c(2,5)] %>% rename('GeneID'="GeneID1"),
                               MeanDis_randomMR_list[,c(3,5)] %>% rename('GeneID'="GeneID2"))

MeanDis_randomMR_list <- split(MeanDis_randomMR_list$MR, MeanDis_randomMR_list$GeneID)

# Calculate Z value to observed MR with random dist
ParalogsPair_MRdb[,"Z"] <- 0

for(i in 1:nrow(ParalogsPair_MRdb)){
  id <- ParalogsPair_MRdb$GeneID1[[i]]
  ParalogsPair_MRdb$Z[i] <- (ParalogsPair_MRdb$MR[i] - mean(MeanDis_randomMR_list[[id]]))/sd(MeanDis_randomMR_list[[id]])
}

hist(ParalogsPair_MRdb$Z, 100)

# Test Z value is lower and larges than expected 
ParalogsPair_MRdb[,"Pvalneg"] <- unlist(lapply(ParalogsPair_MRdb$Z, 
                                               function (x) pnorm(x, lower.tail=TRUE)))
ParalogsPair_MRdb[,"Pvalpos"] <- unlist(lapply(ParalogsPair_MRdb$Z, 
                                               function (x) pnorm(x, lower.tail=FALSE)))

ParalogsPair_MRdb[,'Pval'] <- 1

# If z negative keep left tail
ParalogsPair_MRdb$Pval[ParalogsPair_MRdb$Z < 0]  <- ParalogsPair_MRdb$Pvalneg[ParalogsPair_MRdb$Z < 0] 
# If z positive keep right tail
ParalogsPair_MRdb$Pval[ParalogsPair_MRdb$Z > 0]  <- ParalogsPair_MRdb$Pvalpos[ParalogsPair_MRdb$Z > 0] 

ParalogsPair_MRdb <- ParalogsPair_MRdb[,-c("Pvalpos", 'Pvalneg')]

ParalogsPair_MRdb <- left_join(ParalogsPair_MRdb, saf[,1:2], by=c('GeneID1'='GeneID'))
ParalogsPair_MRdb <- left_join(ParalogsPair_MRdb, saf[,1:2], by=c('GeneID2'='GeneID'))
colnames(ParalogsPair_MRdb) <- gsub('chrAnn.x', 'Chr_G1', colnames(ParalogsPair_MRdb))
colnames(ParalogsPair_MRdb) <- gsub('chrAnn.y', 'Chr_G2', colnames(ParalogsPair_MRdb))

# Paralogs
ParalogsPair_MRdb[,'Type'] <- (ParalogsPair_MRdb$Chr_G1 == ParalogsPair_MRdb$Chr_G2)*1
ParalogsPair_MRdb$Type <- gsub("0", "DiffChr", ParalogsPair_MRdb$Type)
ParalogsPair_MRdb$Type <- gsub("1", "SameChr", ParalogsPair_MRdb$Type)


##################################################################
########    Spearman correlation  of genes paralogs       ######## 
##################################################################

Get_spearman_id1_id2 <- function(id1, DB){
  # gene id from maize copy 2
  id2 <- DB[DB$GeneID1 == id1,]$GeneID2
  # MR_MI.pecanpy.Zm00001d036331.txt
  file1 <- paste0('MR_edgesDB_syntenic/MR_MI.pecanpy.',id1,'.txt')
  file2 <- paste0('MR_edgesDB_syntenic/MR_MI.pecanpy.',id2,'.txt')
  #
  mr.mi1 <- fread(file1)[,1:4]
  geneorder <- mr.mi1$GeneID2
  
  #
  mr.mi1 <- unlist(split(mr.mi1$MR, mr.mi1$GeneID2))
  
  # Read DFs from paralog paris: id2
  mr.mi2 <- lapply(file2, function(x) fread(x)[,1:4])
  names(mr.mi2) <- id2
  
  # keep id2 as vector with geneid and MR values
  mr.mi2 <- lapply(mr.mi2, function(x) unlist(split(x$MR, x$GeneID2)))
  names(mr.mi2) <- id2
  
  # Define vector 1: id1
  v1 <- mr.mi1[geneorder]
  
  # Define empty list to save results
  spearmanResults <- list()
  
  for(id in id2){
    # To keep gene order
    v2 <- mr.mi2[[id]][geneorder] 
    
    # spearman cor of MR profile
    val <- cor(v1, v2, method = 'spearman')
    spearmanResults <- c(spearmanResults, val)
  }
  names(spearmanResults) <- id2
  
  t(as.data.frame(spearmanResults)) %>%
    as.data.frame() %>%
    rownames_to_column("GeneID2") -> spearmanResults
  
  spearmanResults[,'GeneID1'] <- id1
  spearmanResults <- spearmanResults[,c("GeneID1", "GeneID2", "V1")]
  colnames(spearmanResults)[3] <- "SCC"
  
  return(spearmanResults)
}

Paralogs_DB
ParalogsPair_SCCdb <- lapply(unique(Paralogs_DB$GeneID1), 
                             function(x) Get_spearman_id1_id2(x, Paralogs_DB))
table(unlist(lapply(ParalogsPair_SCCdb, nrow)))

ParalogsPair_SCCdb <- rbindlist(ParalogsPair_SCCdb, idcol = F)

## Create random background


for(i in seq(1000)){
  # Calculate random net of paralogs
  DFrandom <- RewireNet(Paralogs_DB[,1:2])
  # write random pair to test by ""
  fname <- paste0("Random_SCC_DB/ParalogsRandom/RandomParalogsPairs_",i, ".txt")
  
  fwrite(DFrandom, fname, quote = F, row.names = F)
  cat(paste0(" .. donde: ", i, " ..\n"))
  # # calculate MR for random pairs
  # DFrandom_SCC <- lapply(unique(DFrandom$GeneID1), function(x) Get_spearman_id1_id2(x, DFrandom))
  # DFrandom_SCC <- rbindlist(DFrandom_SCC, idcol = F)
  # DFrandom_SCC_DB[[i]] <- DFrandom_SCC
  
}
## Read random background results

# Get Mean values
DFrandom_SCC_DB <- list.files(path = 'Random_SCC_DB/SCC_PatalogsRandom/', 
                              pattern = "^SSC_RandomParalogsPairs.*")

DFrandom_SCC_DB <- lapply(DFrandom_SCC_DB, 
                          function(x) fread(paste0("Random_SCC_DB/SCC_PatalogsRandom/", x)))

unlist(lapply(DFrandom_SCC_DB, function(x) nrow(x)))

# distribution of SCC Means
MeanDis_randomSCC <- unlist(lapply(DFrandom_SCC_DB, function(x) mean(x$SCC)))
# distribution of SCC Means  scale between 0 and 1
MeanDis_randomSCC <- unlist(lapply(DFrandom_SCC_DB, function(x) mean( (x$SCC+1)/2)))

# random distribution by gene
MeanDis_randomSCC_list <- rbindlist(DFrandom_SCC_DB, idcol = T)

MeanDis_randomSCC_list <- rbind(MeanDis_randomSCC_list[,c(2,4)] %>% rename('GeneID'="GeneID1"),
                                MeanDis_randomSCC_list[,c(3,4)] %>% rename('GeneID'="GeneID2"))

MeanDis_randomSCC_list <- split(MeanDis_randomSCC_list$SCC, MeanDis_randomSCC_list$GeneID)

# Calculate Z value to observed SCC with random dist
ParalogsPair_SCCdb[,"Z"] <- 0

for(i in 1:nrow(ParalogsPair_SCCdb)){
  id <- ParalogsPair_SCCdb$GeneID1[[i]]
  ParalogsPair_SCCdb$Z[i] <- (ParalogsPair_SCCdb$SCC[i] - mean(MeanDis_randomSCC_list[[id]]))/sd(MeanDis_randomSCC_list[[id]])
}

## Test Z value is lower and larges than expected 
ParalogsPair_SCCdb[,"Pvalneg"] <- unlist(lapply(ParalogsPair_SCCdb$Z, function (x) pnorm(x, lower.tail=TRUE)))
#
ParalogsPair_SCCdb[,"Pvalpos"] <- unlist(lapply(ParalogsPair_SCCdb$Z, function (x) pnorm(x, lower.tail=FALSE)))

ParalogsPair_SCCdb[,'Pval'] <- 1
# if z negative keep left tail
ParalogsPair_SCCdb$Pval[ParalogsPair_SCCdb$Z < 0]  <- ParalogsPair_SCCdb$Pvalneg[ParalogsPair_SCCdb$Z < 0] 
# if z positive keep right tail
ParalogsPair_SCCdb$Pval[ParalogsPair_SCCdb$Z > 0]  <- ParalogsPair_SCCdb$Pvalpos[ParalogsPair_SCCdb$Z > 0]

ParalogsPair_SCCdb <- ParalogsPair_SCCdb[,-c("Pvalpos", 'Pvalneg')]

ParalogsPair_SCCdb <- left_join(ParalogsPair_SCCdb, saf[,1:2], by=c('GeneID1'='GeneID'))
ParalogsPair_SCCdb <- left_join(ParalogsPair_SCCdb, saf[,1:2], by=c('GeneID2'='GeneID'))
colnames(ParalogsPair_SCCdb) <- gsub('chrAnn.x', 'Chr_G1', colnames(ParalogsPair_SCCdb))
colnames(ParalogsPair_SCCdb) <- gsub('chrAnn.y', 'Chr_G2', colnames(ParalogsPair_SCCdb))

#
ParalogsPair_SCCdb[,'Type'] <- (ParalogsPair_SCCdb$Chr_G1 == ParalogsPair_SCCdb$Chr_G2)*1
ParalogsPair_SCCdb$Type <- gsub("0", "DiffChr", ParalogsPair_SCCdb$Type)
ParalogsPair_SCCdb$Type <- gsub("1", "SameChr", ParalogsPair_SCCdb$Type)
##################################################################

##################################################################
####  Combine results from MR and SCC of genes paralogs   ######## 
##################################################################
library(ggpointdensity)

MRdb_SCCdb <- left_join(ParalogsPair_MRdb, 
                        ParalogsPair_SCCdb[,c('GeneID1', 'GeneID2', 'SCC', 'Z')], 
                        by=c("GeneID1","GeneID2"))

colnames(MRdb_SCCdb) <- gsub('Z.x', 'Z_MR', colnames(MRdb_SCCdb))
colnames(MRdb_SCCdb) <- gsub('Z.y', 'Z_SCC', colnames(MRdb_SCCdb))

MRdb_SCCdb <- as.data.table(MRdb_SCCdb)


##################################################################

##################################################################
####                Protein distances                      #######
##################################################################

Pep.dis.files <- list.files(path = 'AAsequnceDistance/', pattern = "^Dis.*")
Pep.dis.files <- Pep.dis.files[grepl(".txt", Pep.dis.files)]

PepDisDB <- lapply(Pep.dis.files, function(x) fread(paste0("AAsequnceDistance/", x)))
PepDisDB <- rbindlist(PepDisDB, idcol = F)

# split gene and pep ids
PepDisDB[,"Pep1"] <- chop(PepDisDB$GeneID1, '[_]', 2)
PepDisDB[,"Pep2"] <- chop(PepDisDB$GeneID2, '[_]', 2)

PepDisDB$GeneID1 <- chop(PepDisDB$GeneID1, '[_]', 1)
PepDisDB$GeneID2 <- chop(PepDisDB$GeneID2, '[_]', 1)

# create index to calculate average distance
PepDisDB[,'Index1'] <- paste0(PepDisDB$GeneID1, ".",PepDisDB$GeneID2)
PepDisDB[,'Index2'] <- paste0(PepDisDB$GeneID2, ".",PepDisDB$GeneID1)

PepDisDB %>%
  dplyr::select(GeneID1, GeneID2, Dis, Index1) %>%
  dplyr::group_by(Index1) %>%
  dplyr::mutate(GeneID1, GeneID2, sDis=sum(Dis), n=length(Dis)) %>% 
  dplyr::rename('Index'='Index1') %>%
  unique()  -> PepDisDB_1

PepDisDB %>%
  dplyr::select(GeneID1, GeneID2, Dis, Index2) %>%
  dplyr::group_by(Index2) %>%
  dplyr::mutate(GeneID1, GeneID2, sDis=sum(Dis), n=length(Dis)) %>% 
  dplyr::rename('Index'='Index2') %>%
  unique()  -> PepDisDB_2

# Get average and kep unique comparisons
PepDisDB_1[,'mDis'] <- PepDisDB_1$sDis/PepDisDB_1$n
PepDisDB_2[,'mDis'] <- PepDisDB_1$sDis/PepDisDB_1$n

# Combine all comparison
PepDisDB <- rbind(unique(PepDisDB_1[,-c(3)]), unique(PepDisDB_1[,-c(3)]))


# Remove redundancies
PepDisDB %>%
  dplyr::group_by(Index) %>%
  dplyr::mutate(GeneID1, GeneID2, Index, sDis=sum(sDis), n=sum(n)) %>% 
  unique()  -> PepDisDB

PepDisDB[,'Dis'] <- PepDisDB$sDis/PepDisDB$n


# PepDisDB[,'Index2'] <- paste0(PepDisDB$GeneID2, ".",PepDisDB$GeneID1)

# MRdb_SCCdb db
MRdb_SCCdb[,'Index1'] <- paste0(MRdb_SCCdb$GeneID1, ".", MRdb_SCCdb$GeneID2)
MRdb_SCCdb[,'Index2'] <- paste0(MRdb_SCCdb$GeneID2, ".", MRdb_SCCdb$GeneID1)

#
MRdb_SCCdb <- left_join(MRdb_SCCdb, PepDisDB[,c(3,6)], by=c('Index1'='Index'))
MRdb_SCCdb <- left_join(MRdb_SCCdb, PepDisDB[,c(3,6)], by=c('Index2'='Index'))

# Replace NA with 0 to combined both columns without duplicate
MRdb_SCCdb$Dis.x[is.na(MRdb_SCCdb$Dis.x)] <- 0
MRdb_SCCdb$Dis.y[is.na(MRdb_SCCdb$Dis.y)] <- 0

MRdb_SCCdb[,"DisHammin"] <- MRdb_SCCdb$Dis.x + MRdb_SCCdb$Dis.y

# Remove temporal columns
MRdb_SCCdb <- MRdb_SCCdb[,-c('Dis.x', 'Dis.y', "Index1", "Index2")]

# Re-defined Zscore based on observed values
MRdb_SCCdb[,"Z_SCC2"] <- scale(MRdb_SCCdb$SCC, center = T, scale = T)
MRdb_SCCdb[,"Z_MR2"] <-  scale(MRdb_SCCdb$MR, center = T, scale = T)


# Define zones based on z score
MRdb_SCCdb[,"ZoneSCC"] <- NA
MRdb_SCCdb[,"ZoneMR"] <- NA

MRdb_SCCdb[,"Index"] <- seq(1, nrow(MRdb_SCCdb))

# Defined zones
MRdb_SCCdb$ZoneSCC[subset(MRdb_SCCdb, Z_SCC > -0.5 & Z_SCC < 0.5)$Index] <- 0
MRdb_SCCdb$ZoneSCC[subset(MRdb_SCCdb, Z_SCC <= -0.5 )$Index] <- -1
MRdb_SCCdb$ZoneSCC[subset(MRdb_SCCdb, Z_SCC >= 0.5 )$Index] <- 1

# Defined zones
MRdb_SCCdb$ZoneMR[subset(MRdb_SCCdb, Z_MR > -0.5 & Z_MR < 0.5)$Index] <- 0
MRdb_SCCdb$ZoneMR[subset(MRdb_SCCdb, Z_MR <= -0.5 )$Index] <- -1
MRdb_SCCdb$ZoneMR[subset(MRdb_SCCdb, Z_MR >= 0.5 )$Index] <- 1

# 
MRdb_SCCdb[,"ZonesIdex"] <- paste0(MRdb_SCCdb$ZoneSCC, MRdb_SCCdb$ZoneMR)
MRdb_SCCdb$ZonesIdex <- gsub("-11","I", MRdb_SCCdb$ZonesIdex)
MRdb_SCCdb$ZonesIdex <- gsub("01","II", MRdb_SCCdb$ZonesIdex)
MRdb_SCCdb$ZonesIdex <- gsub("11","III", MRdb_SCCdb$ZonesIdex)
MRdb_SCCdb$ZonesIdex <- gsub("-10","IV", MRdb_SCCdb$ZonesIdex)
MRdb_SCCdb$ZonesIdex <- gsub("00","V", MRdb_SCCdb$ZonesIdex)
MRdb_SCCdb$ZonesIdex <- gsub("10","VI", MRdb_SCCdb$ZonesIdex)
MRdb_SCCdb$ZonesIdex <- gsub("-1-1","VII", MRdb_SCCdb$ZonesIdex)
MRdb_SCCdb$ZonesIdex <- gsub("0-1","VIII", MRdb_SCCdb$ZonesIdex)
MRdb_SCCdb$ZonesIdex <- gsub("1-1","IX", MRdb_SCCdb$ZonesIdex)

MRdb_SCCdb$ZonesIdex <- factor(MRdb_SCCdb$ZonesIdex, 
                               levels = c("I", "II","III",
                                          "IV", "V", "VI",
                                          "VII", "VIII", "IX"))

# Count zones freq
ZonesIndex <- as.data.table(table(MRdb_SCCdb[,c("Type", "ZonesIdex")]))

ZonesIndex$ZonesIdex <- factor(ZonesIndex$ZonesIdex, 
                               levels = rev(c("I", "II","III",
                                          "IV", "V", "VI",
                                          "VII", "VIII", "IX")))

# MR test for Zs scores
MRdb_SCCdb[, "MR_Zs"] <- sqrt(rank(MRdb_SCCdb$Z_MR2)*rank(-MRdb_SCCdb$Z_SCC2))



##################################################################
####              Top 10 examples distances               #######
##################################################################

MRdb_SCCdb

Get_TFnet <- function(tf){
  ## Predicted Targets and regulators of TF
  net_MR_in_1 <- unique(subset(MRMI_full,  V2 %in% tf[1])$V1)
  net_MR_in_2 <- unique(subset(MRMI_full,  V2 %in% tf[2])$V1)
  
  #
  net_MR_out_1 <- unique(subset(MRMI_full, V1 %in% tf[1])$V2)
  net_MR_out_2 <- unique(subset(MRMI_full, V1 %in% tf[2])$V2)
  
  ## Observed Targets and regulators of TF
  net_obs_in_1 <- unique(subset(FullNet,  Target %in% tf[1])$Source)
  net_obs_in_2 <- unique(subset(FullNet,  Target %in% tf[2])$Source)
  #
  net_obs_out_1 <- unique(subset(FullNet, Source %in% tf[1])$Target)
  net_obs_out_2 <- unique(subset(FullNet, Source %in% tf[2])$Target)
  
  # Interaction groups
  MR_degree <- list("G1"=unique(c(net_MR_in_1, net_MR_out_1)), 
                    "G2"=unique(c(net_MR_in_2, net_MR_out_2)))
  #names(MR_degree) <- tf
  
  Obs_degree <- list("G1"=unique(c(net_obs_in_1, net_obs_out_1)), 
                     "G2"=unique(c(net_obs_in_2, net_obs_out_2)))
  #names(Obs_degree) <- tf
  #
  MR_degree <- venn(MR_degree)
  Obs_degree <- venn(Obs_degree)
  #
  MR_degree <- as.list(attr(MR_degree, "intersections"))
  MR_degree <- plyr::ldply(MR_degree, data.table)
  MR_degree <- as.data.table(MR_degree)
  MR_degree <- as.data.table(table(MR_degree$.id))
  #
  Obs_degree <- as.list(attr(Obs_degree, "intersections"))
  Obs_degree <- plyr::ldply(Obs_degree, data.table)
  Obs_degree <- as.data.table(Obs_degree)
  Obs_degree <- as.data.table(table(Obs_degree$.id))
  
  # add types of networks
  MR_degree[,'Net'] <- "MR"
  Obs_degree[,'Net'] <- "Obs"
  
  # observed Targets and regulators of TF
  DF <- rbind(MR_degree, Obs_degree)
  # 
  DF[,'TFname'] <- ReplaceName(DF$V1)
  
  return(DF)
}

Example_1 <- Get_TFnet(c('Zm00001d043937','Zm00001d011537'))
TFdic <- as_tibble(read.table("Data/Annotations/TF_Id_Name.txt", h=F, stringsAsFactors = F))


# Examples to plot
Top10 <- subset(MRdb_SCCdb, Z_SCC2 > 1.503 & Z_MR2 < -1.503)
Tail10 <- subset(MRdb_SCCdb, Z_SCC2 < -1.4926 & Z_MR2 > 1.4926 & ZonesIdex=='I')

# Combine examples
DF_ExamplesInput <- rbind(Top10, Tail10)

ExamplesVec <- list()

for(i in 1:nrow(DF_Examples)){
  #
  ExamplesVec[[paste0(DF_Examples$GeneID1[i], '_',DF_Examples$GeneID2[i])]] <- c(DF_Examples$GeneID1[i],DF_Examples$GeneID2[i])
}


DF_Examples <- lapply(ExamplesVec, Get_TFnet)
DF_Examples <- rbindlist(DF_Examples, idcol = T)
DF_Examples[,'TFname'] <- ReplaceName(DF_Examples$.id)

# defined order
DF_ExamplesInput %>%
  arrange(MR_Zs) -> DF_ExamplesInput

NamesOrder <- paste(DF_ExamplesInput$GeneID1, "_", DF_ExamplesInput$GeneID2, sep = "")
NamesOrder <- ReplaceName(NamesOrder)

DF_Examples$TFname <- factor(DF_Examples$TFname, levels = NamesOrder)
DF_Examples$V1 <- factor(DF_Examples$V1, levels = rev(c('G1', 'G1:G2','G2')))

#
as.data.table(table(MRdb_SCCdb[,c('Type','ZonesIdex')]))

## Calculate Jaccard index
Get_jaccard <- function(tf){
  ## Predicted Targets and regulators of TF
  net_MR_in_1 <- unique(subset(MRMI_full,  V2 %in% tf[1])$V1)
  net_MR_in_2 <- unique(subset(MRMI_full,  V2 %in% tf[2])$V1)
  
  #
  net_MR_out_1 <- unique(subset(MRMI_full, V1 %in% tf[1])$V2)
  net_MR_out_2 <- unique(subset(MRMI_full, V1 %in% tf[2])$V2)
  
  
  ## Observed Targets and regulators of TF
  net_obs_in_1 <- unique(subset(FullNet,  Target %in% tf[1])$Source)
  net_obs_in_2 <- unique(subset(FullNet,  Target %in% tf[2])$Source)
  #
  net_obs_out_1 <- unique(subset(FullNet, Source %in% tf[1])$Target)
  net_obs_out_2 <- unique(subset(FullNet, Source %in% tf[2])$Target)
  
  
  # all interactions : MR
  net_1 <- unique(c(net_MR_in_1, net_MR_out_1))
  net_2 <- unique(c(net_MR_in_2, net_MR_out_2))
  
  # all interactions : Obs
  obsnet_1 <- unique(c(net_obs_in_1, net_obs_out_1))
  obsnet_2 <- unique(c(net_obs_in_2, net_obs_out_2))
  
  df <- tibble(Index=paste0(tf, collapse = '_'), 
               J_MR =jaccard(net_1, net_2), 
               J_Obs=jaccard(obsnet_1, obsnet_2))
  
  return(df)
}

# define TF pairs to test
AllPairs <- list()

for(i in 1:nrow(MRdb_SCCdb)){
  #
  AllPairs[[paste0(MRdb_SCCdb$GeneID1[i], '_',MRdb_SCCdb$GeneID2[i])]] <- c(MRdb_SCCdb$GeneID1[i],MRdb_SCCdb$GeneID2[i])
}


DFjaccard <- lapply(AllPairs, Get_jaccard)
# 
DFjaccard <- rbindlist(DFjaccard, idcol = F)

MRdb_SCCdb_Reduce <- MRdb_SCCdb[,c(1,2,5,6,11,12,18,19)]
MRdb_SCCdb_Reduce[,"Index"] <- paste0(MRdb_SCCdb_Reduce$GeneID1, "_", MRdb_SCCdb_Reduce$GeneID2)


# Add Jaccard values
DFjaccard <- left_join(DFjaccard[,-c(4,5)], MRdb_SCCdb_Reduce[,c(7,8, 9)], by="Index")

##################################################################


##################################################################
####                Co-expression distance                  ######
##################################################################

wPCC_files <- list.files(path = 'wPCC_Distance/', pattern = "^wPCC_Zm*")

DB_wPCC_pairs <-lapply(wPCC_files, function(name){
  x <- try(read.table(paste("wPCC_Distance/",name, sep=""), fill=T, head=F, stringsAsFactors=FALSE))
  if(inherits(x, "try-error"))
    return(NULL)
  else
    return(x)
})

#DB_wPCC_pairs <- lapply(wPCC_files, function(x) read.table(paste0("wPCC_Distance/", x), fill=T))

names(DB_wPCC_pairs) <- gsub('.txt','', gsub('wPCC_','', wPCC_files))

#
mask <- unlist(lapply(DB_wPCC_pairs, function(x) !is.null(x)))
DB_wPCC_pairs <- DB_wPCC_pairs[mask]

# Filter for empty interactions
mask <- unlist(lapply(DB_wPCC_pairs, function(x) ncol(x) > 1))
DB_wPCC_pairs <- DB_wPCC_pairs[mask]
DB_wPCC_pairs$Zm00001d001824_Zm00001d026628

# Filter nets without value
hist(unlist(lapply(DB_wPCC_pairs, function(x) nrow(x))))
DB_wPCC_pairs <- lapply(DB_wPCC_pairs, function(x) unique(x[!is.na(x$V4),]))
hist(unlist(lapply(DB_wPCC_pairs, function(x) nrow(x))))

nrow(DB_wPCC_pairs$Zm00001d001865_Zm00001d026594)

DB_wPCC_pairs <- rbindlist(DB_wPCC_pairs, idcol = T)
DB_wPCC_pairs$V1 <- chop(DB_wPCC_pairs$V1, '[.]', 1)

colnames(DB_wPCC_pairs) <- c('Index', 'wPCC.id', 'GeneID1', 'GeneID2', 'wPCC')

# average wPCC
DB_wPCC_pairs %>%
  dplyr::select(Index, wPCC) %>%
  dplyr::group_by(Index) %>%
  dplyr::summarise(wPCC=mean(wPCC)) -> DB_wPCCm

# Subset(ParalogsPair_MRdb, GeneID2=="Zm00001d036803")
MRdb_SCCdb_Reduce[,10]
MRdb_SCCdb_Reduce <- left_join(MRdb_SCCdb_Reduce[,-c(10)], DB_wPCCm, by="Index")
MRdb_SCCdb_Reduce$wPCC[is.na(MRdb_SCCdb_Reduce$wPCC)] <- 0

##################################################################


##################################################################
########                      Plots                       ######## 
##################################################################

###########
# Plot MR between gene pairs
###########

######
# V1: density plot
######
ggplot(ParalogsPair_MRdb, aes(x=MR/1000)) + 
  geom_density(alpha=.3, colour="#FF4500", fill="#FFCC99") +
  scale_x_continuous(expand = c(0,0), labels = comma) +
  theme_pubclean() + 
  ylab('Density') + 
  xlab(bquote("Gene pairs" ~ MR[MI] ~ "(1000x)")) +
  theme(
    #strip.text.y = element_text(size = 5, angle = 0), 
    #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.text=element_text(size=10), 
    text = element_text(size=10, family="Times", color='black')) +
  facet_grid(.~ Type) -> Plot_3a
######

######
# V2: MR distribution boxplot
######
ParalogsPair_MRdb

lower_threshold <-  subset(ParalogsPair_MRdb, Z <= -0.5)
upper_threshold <- subset(ParalogsPair_MRdb,  Z >= 0.5)

lower_threshold <-  subset(lower_threshold, Z==max(lower_threshold$Z))
upper_threshold <- subset(upper_threshold, Z==min(upper_threshold$Z))
                    
lower_threshold <- lower_threshold$MR/1000
upper_threshold <- upper_threshold$MR/1000

ggplot(ParalogsPair_MRdb, aes(y=MR/1000, x=Type)) + 
  geom_boxplot(alpha=.3, colour="#FF4500", fill="#FFCC99", notch = T) +
  geom_hline(yintercept = c(lower_threshold, upper_threshold), linetype="dashed") + 
  theme_pubclean() + 
  scale_y_reverse() +
  xlab('') + 
  ylab(bquote("Gene pair" ~ MR[MI] ~ "(1e3)")) +
  theme(
    #strip.text.y = element_text(size = 5, angle = 0), 
    #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.text=element_text(size=10), 
    text = element_text(size=10, family="Times", color='black')) + 
  stat_compare_means(method = 'wilcox.test', 
                     size=3, paired = F, 
                     label="p.signif", 
                     label.x = 1.5) -> Plot_3a

Plot_3a
######

#####
# Plot SCC between gene pair profiles 
#####
ParalogsPair_SCCdb[,"SCCscale"] <-  (ParalogsPair_SCCdb$SCC + 1)/2

lower_threshold_scc <-  subset(ParalogsPair_SCCdb, Z <= -0.5)
upper_threshold_scc <-  subset(ParalogsPair_SCCdb, Z >= 0.5)

lower_threshold_scc <-  subset(lower_threshold_scc, Z==max(lower_threshold_scc$Z))
upper_threshold_scc <- subset(upper_threshold_scc, Z==min(upper_threshold_scc$Z))

lower_threshold_scc <- lower_threshold_scc$SCC
upper_threshold_scc <- upper_threshold_scc$SCC

ggplot(ParalogsPair_SCCdb, aes(y=SCC, x=Type)) + 
  geom_boxplot(alpha=.3, colour="#FF007F", fill="#FF99FF", notch = T) +
  geom_hline(yintercept = c(lower_threshold_scc, upper_threshold_scc), linetype="dashed") + 
  #scale_y_continuous(limits = c(-0.2, 0.6)) + 
  theme_pubclean() + 
  xlab('') + 
  ylab(bquote("Gene pair" ~ SCC[MR])) +
  theme(
    #strip.text.y = element_text(size = 5, angle = 0), 
    #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.text=element_text(size=10), 
    text = element_text(size=10, family="Times", color='black')) + 
  stat_compare_means(method = 'wilcox.test',  size=3, paired = F, 
                     label="p.signif", 
                     label.x = 1.5) -> Plot_3b

Plot_3b
#####
  
#####
# Plot MR vs SCC between gene pair
#####

lower_threshold
MRdb_SCCdb
table(test$ZonesIdex)

MRdb_SCCdb %>%
  ggplot(aes(x=Z_SCC, y=Z_MR)) +
  #geom_point(alpha=0.5) + #geom_hex(bins = 80) +
  geom_pointdensity() +
  scale_color_viridis(labels=comma) +
  scale_y_reverse() +
  theme_pubclean() + 
  geom_hline(yintercept = c(-0.5, 0.5), linetype="dashed") + 
  geom_vline(xintercept = c(-0.5, 0.5), linetype="dashed") + 
  ylab(bquote("Gene pairs" ~ Z[MR])) + 
  xlab(bquote("Gene pairs" ~ Z[SCC])) +
  theme(strip.text.x = element_text(size = 10), 
        axis.text=element_text(size=10),
        legend.position = 'right',
        legend.direction='vertical',
        legend.key.size = unit(0.3, 'cm'),
        text = element_text(size=10)) -> Plot_MR_vs_SCC

Plot_MR_vs_SCC

Plot_3abc <- Plot_3a + Plot_3b + Plot_MR_vs_SCC + plot_layout(widths = c(0.2,0.2,1))
Plot_3abc

pdf("Plots/Plot_3bcd.pdf", width=5, height=2)
print(Plot_3abc)
dev.off()

cor(-MRdb_SCCdb$Z_MR, MRdb_SCCdb$Z_SCC, method = 'spearman')

#####

#####
# Heat map with exampkes from top and tail Zscore values
#####

# redefined top value to reduce color scale
DF_Examples[,'N2'] <- DF_Examples$N

DF_Examples$N2[DF_Examples$N2 >1000] <- 1000

ggplot(subset(DF_Examples, Net=='MR'), aes(y=V1, x=TFname, fill=N2, label=N)) +
  geom_tile() +
  geom_text(aes(color=N), size=2) +
  theme_bw() + 
  scale_fill_viridis(option = 'A') + 
  scale_color_viridis(option = 'D', direction = -1) + 
  scale_y_discrete(expand = c(0,0)) +
  ylab('') + xlab("") +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text = element_text(size = 10), 
        axis.text.x = element_blank(),
        plot.title = element_text(size = 10),
        legend.position  = 'bottom',
        legend.key.size = unit(0.3, 'cm'),
        text = element_text(size=10, family="Times")) +
  guides(color = 'none') +
  facet_grid(Net ~. ) -> Plot_4e1

ggplot(subset(DF_Examples, Net=='Obs'), aes(y=V1, x=TFname, fill=N2, label=N)) +
  geom_tile() +
  geom_text(aes(color=N), size=2) +
  theme_bw() + 
  scale_fill_viridis(option = 'A') + 
  scale_color_viridis(option = 'D', direction = -1) + 
  scale_y_discrete(expand = c(0,0)) +
  ylab('') + xlab("") +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text = element_text(size = 10), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size=4),
        plot.title = element_text(size = 10),
        legend.position  = 'bottom',
        legend.key.size = unit(0.3, 'cm'),
        text = element_text(size=10, family="Times")) +
  guides(color = 'none') +
  facet_grid(Net ~. ) -> Plot_4e2

Plot_4e <- Plot_4e1/Plot_4e2 + plot_layout(guides = "collect") & theme(legend.position  = 'bottom')

pdf('Plots/Plot_4e.pdf', width = 5.5, height = 3)
print(Plot_4e)
dev.off()
#####

#####
# Plot number of TF pairs by zones
#####

colnames(ZonesIndex)[1] <- 'Chr'
ZonesIndex$Chr <- gsub('DiffChr', 'Diff.', ZonesIndex$Chr)
ZonesIndex$Chr <- gsub('SameChr', 'Same', ZonesIndex$Chr)

# boxplot of zones freq
ZonesIndex %>%
  ggplot(aes(y=ZonesIdex, x=N, fill=Chr)) +
  geom_bar(stat="identity", alpha=0.5) + 
  theme_pubclean() + 
  ylab("Z-score regions") + 
  xlab("TF pairs") +
  theme(strip.text.x = element_text(size = 10), 
        axis.text=element_text(size=10),
        axis.text.x =element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = 'bottom',
        legend.direction='vertical',
        legend.key.size = unit(0.3, 'cm'),
        text = element_text(size=10)) -> Plot_e

Plot_e
pdf("Plots/Plot_3e.pdf", width=1.7, height=3)
print(Plot_e)
dev.off()
#####

#####
# Plot Jaccard by zones
#####

DFjaccard[,"J_MR2"] <- DFjaccard$J_MR

DFjaccard$J_MR2[DFjaccard$J_MR2 >0.07] <- 0.07

ggplot(DFjaccard, aes(y=J_MR2, x=ZonesIdex)) +
  geom_boxplot(notch = T,  fill='#87CEFA', color='#00008B') +
  theme_pubclean() +
  ylab("Jaccard Index") + 
  xlab("Z-score bins") +
  theme(strip.text.x = element_text(size = 10), 
        axis.text=element_text(size=10),
        #axis.text.x =element_text(angle = 0 , vjust = 0.5, hjust = 1),
        #legend.position = 'bottom',
        #legend.direction='vertical',
        #legend.key.size = unit(0.3, 'cm'),
        text = element_text(size=10)) +
  stat_compare_means(label = "p.signif", 
                     method = "t.test",
                     ref.group = "I") -> Plot_4g

Plot_4g


#####


#####
# Plot Z_SCC vs DisHammin
#####

MRdb_SCCdb_Reduce <- left_join(MRdb_SCCdb_Reduce[,-c(11)], DFjaccard[,c(1,2,3)], by='Index')

Labels1 <- subset(MRdb_SCCdb_Reduce, DisHammin <= 0.3 & ZonesIdex %in% c('I'))
Labels2 <- subset(MRdb_SCCdb_Reduce, DisHammin <= 0.3 & ZonesIdex %in% c('IX'))
#Labels3 <- subset(MRdb_SCCdb_Reduce, Z_SCC >=  4)

Labels <- rbind(Labels1, Labels2)
Labels[,'Index2'] <- ReplaceName(Labels$Index)

pos <- position_jitter(width = 0.1, seed = 2)

subset(MRdb_SCCdb_Reduce, ZonesIdex %in% c('I', 'IX')) %>%
  ggplot(aes(x=ZonesIdex, y=DisHammin, alpha=0.5)) +
  #geom_point(alpha=0.5) + #geom_hex(bins = 80) +
  #geom_pointdensity()+
  geom_jitter(position = pos, size=1) +
  geom_text_repel(data=Labels, aes(x=ZonesIdex, y=, label=Index2),
                  size=1.5, 
                  segment.size=0.3, color='#1E90FF',
                  point.padding = 0,
                  box.padding = 1,
                  max.overlaps = Inf,
                  position = pos) +
  scale_color_viridis(labels=comma, option = 'B') +
  theme_pubclean() + 
  #geom_vline(xintercept = c(-0.5, 0.5), linetype="dashed") +
  ylab('AA hammin distance') + 
  xlab('Z−score bins') +
  theme(strip.text.x = element_text(size = 10), 
        axis.text=element_text(size=10),
        legend.position = 'none',
        legend.direction='horizontal',
        legend.key.size = unit(0.4, 'cm'),
        text = element_text(size=10)) +
  stat_compare_means(label = "..p.signif..", method = "t.test") -> Plot_4h

Plot_4h
#ylab('wPCC') + 
#xlab(bquote("TF pairs " ~ Z[SCC])) +

Labels_zn_1 <- subset(MRdb_SCCdb_Reduce, DisHammin <= 0.25 & ZonesIdex == 'I')
Labels_zn_2 <- subset(MRdb_SCCdb_Reduce, wPCC < -0.1 & ZonesIdex == 'I')
Labels_zn_3 <- subset(MRdb_SCCdb_Reduce, wPCC > 0.9 & DisHammin >= 0.5 & ZonesIdex == 'I')

Labels_zn <- rbind(Labels_zn_1, Labels_zn_2, Labels_zn_3)
Labels_zn[,'Index2'] <- ReplaceName(Labels_zn$Index)

subset(MRdb_SCCdb_Reduce, ZonesIdex %in% c('I')) %>%
  ggplot(aes(x=wPCC, y=DisHammin)) +
  geom_pointdensity() +
  geom_text_repel(data=Labels_zn, aes(x=wPCC, y=DisHammin, label=Index2),
                  size=1.5, segment.size=0.3, color='#1E90FF',
                  point.padding = 0, 
                  box.padding = 0.5,
                  max.overlaps = Inf) +
  scale_color_viridis(labels=comma, option = 'B') +
  theme_pubclean() + 
  scale_x_continuous(limits = c(-0.4,1)) +
  xlab('wPCC') + 
  ylab("AA hammin distance") +
  theme(strip.text.x = element_text(size = 10), 
        axis.text=element_text(size=10),
        legend.position = 'bottom',
        legend.direction='horizontal',
        legend.key.size = unit(0.4, 'cm'),
        text = element_text(size=10)) -> Plot_4i
Plot_4i

Labels_zp_1 <- subset(MRdb_SCCdb_Reduce, DisHammin <= 0.30 & ZonesIdex == 'IX')
Labels_zp_2 <- subset(MRdb_SCCdb_Reduce, wPCC < -0.1 & ZonesIdex == 'IX')
Labels_zp_3 <- subset(MRdb_SCCdb_Reduce, wPCC > 0.9 & ZonesIdex == 'IX')

Labels_zp <- rbind(Labels_zp_1, Labels_zp_2, Labels_zp_3) # Labels_zp_2, Labels_zp_3
Labels_zp[,'Index2'] <- ReplaceName(Labels_zp$Index)

subset(MRdb_SCCdb_Reduce, ZonesIdex %in% c('IX')) %>%
  ggplot(aes(x=wPCC, y=DisHammin)) +
  geom_pointdensity() +
  geom_text_repel(data=Labels_zp, aes(x=wPCC, y=DisHammin, label=Index2),
                  size=1.5, segment.size=0.3, color='#1E90FF',
                  point.padding = 0, 
                  box.padding = 0.5,
                  max.overlaps = Inf) +
  scale_color_viridis(labels=comma, option = 'B') +
  theme_pubclean() + 
  scale_x_continuous(limits = c(-0.4,1)) +
  xlab('wPCC') + 
  ylab("AA hammin distance") +
  theme(strip.text.x = element_text(size = 10), 
        axis.text=element_text(size=10),
        legend.position = 'bottom',
        legend.direction='horizontal',
        legend.key.size = unit(0.4, 'cm'),
        text = element_text(size=10)) -> Plot_4j

Plot_4j

MRdb_SCCdb_Reduce %>%
  dplyr::filter(wPCC != 0) %>%
  dplyr::group_by(ZonesIdex) %>%
  dplyr::summarise(wPCC=mean((wPCC)))

subset(MRdb_SCCdb_Reduce, ZonesIdex %in% c('I','IX')) %>%
  ggplot(aes(x=ZonesIdex, y=wPCC))+
  geom_boxplot(notch = T) +
  stat_compare_means(label = "..p.signif..", method = "t.test")


# TF names
TFdic <- as_tibble(read.table("Data/Annotations/TF_Id_Name.txt", h=F, stringsAsFactors = F))


Plot_4ghi <- {{Plot_4g/Plot_4h}|{Plot_4i/Plot_4j}} + plot_layout(widths = c(0.5,1))

pdf("Plots/Plot_4ghij.pdf", width=7.5, height=5.5)
print(Plot_4ghi)
dev.off()

##################################################################
Z_SCC
wPCCClass1 <- subset(MRdb_SCCdb_Reduce, Z_SCC <= -0.5)
wPCCClass2 <- subset(MRdb_SCCdb_Reduce, Z_SCC > -0.5 & Z_SCC < 0.5)
wPCCClass3 <- subset(MRdb_SCCdb_Reduce, Z_SCC >= 0.5)

wPCCClass1[,'Class'] <- '-0.5'
wPCCClass2[,'Class'] <- '0'
wPCCClass3[,'Class'] <- '0.5'

wPCCClass <- rbind(wPCCClass1, wPCCClass2, wPCCClass3)

wPCCClass %>%
  ggplot(aes(x=Class, y=wPCC)) +
  geom_boxplot(notch = T) +
  stat_compare_means(label = "p.signif", 
                     method = "t.test",
                     ref.group = "-0.5")


