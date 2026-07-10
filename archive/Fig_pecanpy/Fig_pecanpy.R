library(tidyverse)
library(ggrepel)
library(ggpubr)
library(data.table)
#library(GeneOverlap)
library(viridis)
library(ComplexHeatmap)
#library(fgsea)
library(reshape2)
library(circlize)
library(ggVennDiagram)
library(scales)
library(purrr)
library(gplots)
library(ggplot2)
library(PCAtools)
library(factoextra)

#

##################################################
##########          Functions        #############
##################################################

ReplaceNamePWY <- function(ids){
  
  for (i in 1:nrow(CornCYC)){
    w <- paste0('\\<', CornCYC$Pathway.id[i], '\\>')
    ids <- gsub(w, CornCYC$Pathway.name[i], ids)
    #ids <- gsub("_", " ", ids)
  }
  return(ids)
}
#
ReplaceName <- function(ids){
  # for (i in 1:nrow(Top45)){
  #   ids <- gsub(Top45$V1[i], Top45$V2[i], ids)
  # }
  
  for (i in 1:nrow(TFdic)){
    ids <- gsub(TFdic$TF.v4[i], TFdic$TF.Name[i], ids)
  }
  return(ids)
}

##################################################
##########        Annotations       ##############
##################################################

#### top 45 ####
Top45 <- as_tibble(read.table("Data/Annotations/Top45.txt", h=F, stringsAsFactors = F))

# TF names
TFdic <- as_tibble(read.table("Data/Annotations/TF_Id_Name.txt", h=T, stringsAsFactors = F))


# Phenolic related genes
PheGenes <- as_tibble(read.table("Data/Annotations/LinaPheGenes2020.txt", h=T, sep = "\t", quote="", stringsAsFactors = F))

# TF and CoReg
TF_CoR <- as_tibble(read.table("Data/Annotations/TF_CoR_Mazie.txt", h=T, stringsAsFactors = F))

## Y1H network
Y1H <- as_tibble(read.table("Data/Annotations/Y1H.data.txt", h=T, stringsAsFactors = F))[,4:6]
ReplaceName(Y1H$TF.v4)

## syntenic genes 
Syntenic <- as_tibble(read.table("Data/Annotations/Zm.v4.synteny.genes.txt", h=T, stringsAsFactors = F, sep = '\t'))
Syntenic <- subset(Syntenic, gene_synteny=='syntenic')$gene_id

# CornCYC
CornCYC <- as_tibble(read.table("Data/Annotations/corn_pathways.0210325.reduced.txt", h=T, stringsAsFactors = F))
CornCYC <- CornCYC[!(CornCYC$Gene.id %in% "unknown"),]
CornCYC <- subset(CornCYC, Gene.id %in% Syntenic)
#CornCYC[,"Class"] <- "Enzyme"
colnames(CornCYC)[3] <- "GeneID"
CornCYC$Pathway.id <- gsub("-", "_", CornCYC$Pathway.id)

## syntenic genes 
Syntenic <- as_tibble(read.table("Data/Annotations/Zm.v4.synteny.genes.txt", h=T, stringsAsFactors = F, sep = '\t'))
Syntenic <- subset(Syntenic, gene_synteny=='syntenic')$gene_id

# PDI
PDI <- as_tibble(unique(fread("../Fig_PDI/PDI_NetworkFinal.10_11_2021.txt")[,c(2,3)]))
colnames(PDI)[1] <- "Source"

# CoExp
CoExp <- as_tibble(unique(fread("../Fig_Coexpression/CoExp_NetworkFinal.10_11_2021.txt")))
colnames(CoExp)[2] <- "Source"
CoExp <- unique(CoExp[,2:3])

# teQTL
teQTL <- as_tibble(unique(fread("../Fig_transeQTL/teQTL_NetworkFinal.10_11_2021.txt")))
colnames(teQTL)[1] <- "Source"

teQTLtf <- subset(teQTL, Source %in% unique(c(TF_CoR$GeneID, PDI$Source, CoExp$Source))) 

All_TFs <- unique(c(PDI$Source, CoExp$Source, TF_CoR$GeneID))

##################################################
##########        Pencanty          ##############
##################################################

peGAN <- as.data.frame(fread("pecanpy_GAN.txt" , skip=1, header=F))
peGRN <- as.data.frame(fread("pecanpy_GRN.txt" , skip=1, header=F))
peCEN <- as.data.frame(fread("pecanpy_CEN.txt" , skip=1, header=F))

row.names(peGAN) <- peGAN$V1
row.names(peGRN) <- peGRN$V1
row.names(peCEN) <- peCEN$V1

# peGAN <- peGAN[,-c(1)]
# peGRN <- peGRN[,-c(1)]
# peCEN <- peCEN[,-c(1)]

peGAN[1:5,1:5]
peGRN[1:5,1:5]
peCEN[1:5,1:5]


##################################################
######     dimension reduction: PCA       ########
##################################################

pcaGAN <- pca(peGAN[,-c(1)]) # removeVar = 0.01
pcaGRN <- pca(peGRN[,-c(1)])
pcaCEN <- pca(peCEN[,-c(1)])

screeplot(pcaGAN, axisLabSize = 18, titleLabSize = 22, components = 1:20)
screeplot(pcaGRN, axisLabSize = 18, titleLabSize = 22, components = 1:20)
screeplot(pcaCEN, axisLabSize = 18, titleLabSize = 22, components = 1:20)

##################################################

##################################################
######      Combine pecanpy       ########
##################################################

length(row.names(peGAN))
length(row.names(peGRN))
length(row.names(peCEN))

# scale by gene
peGAN[,-c(1)] <- apply(peGAN[,-c(1)], 2, scale)
peGRN[,-c(1)] <- apply(peGRN[,-c(1)], 2, scale)
peCEN[,-c(1)] <- apply(peCEN[,-c(1)], 2, scale)

colnames(peGAN)[-c(1)] <- paste0("GAN_", colnames(peGAN)[-c(1)])
colnames(peGRN)[-c(1)] <- paste0("GRN_", colnames(peGRN)[-c(1)])
colnames(peCEN)[-c(1)] <- paste0("CEN_", colnames(peCEN)[-c(1)])

peAll <- as_tibble(data.frame(gid=Reduce(intersect, list(row.names(peGAN), row.names(peGRN), row.names(peCEN)))))
peGRN_CEN <- as_tibble(data.frame(gid=Reduce(intersect, list(row.names(peGAN), row.names(peGRN), row.names(peCEN)))))


peAll <- left_join(peAll, peGRN, by=c("gid"="V1"))
peAll <- left_join(peAll, peCEN, by=c("gid"="V1"))
peAll <- left_join(peAll, peGAN, by=c("gid"="V1"))

peGRN_CEN <- left_join(peAll, peGRN, by=c("gid"="V1"))
peGRN_CEN <- left_join(peAll, peCEN, by=c("gid"="V1"))

write.table(peAll, "All_pecanoy_layer.txt", row.names = F, quote = F, sep = "\t")
write.table(peGRN_CEN, "GRN_CEN_pecanoy_layer.txt", row.names = F, quote = F, sep = "\t")

peAll[1:5,1:5]

#
#dim(peAll)
#heatmap.2(as.matrix(peGRN_CEN[1:100,-c(1)]),trace="none")


##
library(WGCNA)
#################################################################################
######################             WGCNA                   ######################
#################################################################################
options(stringsAsFactors = FALSE);
enableWGCNAThreads()

# Cosine distance
dcos <- fread("CosMatrix_All_layer.txt")
row.names(dcos) <- dcos$gid
dcos <- dcos[,-c(1)]


# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(t(peAll[,-c(1)]), powerVector = powers, verbose = 5)

# Plot the results:

### Step=by--step network construction 

# Turn adjacency into topological overlap 
TOM = TOMsimilarity(as.matrix(dcos)) # using cos as dis
dissTOM = 1-dcos


# Call the hierarchical clustering function 
geneTree = hclust(as.dist(dissTOM), method = "average")

# Plot the resulting clustering tree (dendrogram) 
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", 
     main = "Gene clustering on TOM-based dissimilarity", 
     labels = FALSE, hang = 0.04);

# We like large modules, so we set the minimum module size relatively high: 
minModuleSize = 30; 

# Module identification using dynamic tree cut: 
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, 
                            deepSplit = 2, pamRespectsDendro = FALSE, 
                            minClusterSize = minModuleSize)
table(dynamicMods)

# Convert labels to colors for plotting
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)


# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", 
                    dendroLabels = FALSE, hang = 0.03, 
                    addGuide = TRUE, guideHang = 0.05, 
                    main = "Gene dendrogram and module colors")

## Combined modules
# Calculate eigengenes
MEList = moduleEigengenes(t(peAll[,-c(1)]), colors = dynamicColors)
MEs = MEList$eigengenes 

# Calculate dissimilarity of module eigengenes 
MEDiss = 1-cor(MEs)

# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")

# Plot the result
# Size: 7x12
sizeGrWindow(7, 6) 
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")

# Threshold
MEDissThres = 0.25 
# Plot the cut line into the dendrogram 
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(t(peAll[,-c(1)]), dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors
# Eigengenes of the new merged modules
mergedMEs = merge$newMEs

sizeGrWindow(12, 9)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), 
                    c("Dynamic Tree Cut", "Merged dynamic"), 
                    dendroLabels = FALSE, hang = 0.03, 
                    addGuide = TRUE, guideHang = 0.05)

# Rename to moduleColors
moduleColors = mergedColors

#
# Select modules
modules = tibble(Module=unique(moduleColors), M=seq(1:length(unique(moduleColors))))

# Select module probes
probes = row.names(peAll[,-c(1)])

# Empty DF
WGCNA_Modules <- as_tibble(as.data.frame(matrix(0, nrow = 0, ncol = 2)))
colnames(WGCNA_Modules) <- c("GeneID", "Module")

for (m in modules$Module) {
  inModule = (moduleColors==m)
  tem <- tibble(GeneID=probes[inModule], Module=m)
  WGCNA_Modules <- rbind(WGCNA_Modules, tem)
}

head(as.data.frame(table(moduleColors)))
head(as.data.frame(table(WGCNA_Modules$Module)))


## Construct numerical labels corresponding to the colors
# colorOrder = c("grey", standardColors(104)) standardColors(1)
# moduleLabels = match(moduleColors, colorOrder)-1;
# MEs = mergedMEs;
#################################################################################

####################################
########   save Modules     ########
####################################



#WGCNA_Modules <- tibble(GeneID= names(merge$colors), Module=merge$colors)

WGCNA_Modules <- left_join(WGCNA_Modules, saf[,1:2], by="GeneID")
WGCNA_Modules <- left_join(WGCNA_Modules, modules, by="Module")
#WGCNA_Modules["TFn"] <- ReplaceName(WGCNA_Modules$GeneID)

tem <- as_tibble(as.data.frame(table(WGCNA_Modules[,3:4])))
tem <- subset(tem, Freq>0)

ggplot(subset(tem, M %in% c(28, 46, 64)), aes(x=chrAnn, y=M, size=Freq))+
  geom_point()

subset(WGCNA_Modules, GeneID %in% TFdic[TFdic$TF.Name %in% c("P1", "P2", "MYB31","WRKY53"),]$TF.v4)
subset(WGCNA_Modules, GeneID %in% TFdic[TFdic$TF.Name %in% c("MYB31","WRKY53"),]$TF.v4)

####################################################
########     Test with known regulators     ########
########     CornCYC Enrichment test        ########
####################################################

# 
library(GeneOverlap)

CmeanPy_Modules <- as_tibble(read.table("CmeanGroupsPy.txt"))
colnames(CmeanPy_Modules)[2] <- "M"

subset(WGCNA_Modules, GeneID %in% TFdic[TFdic$TF.Name %in% c("MYB31","WRKY53"),]$TF.v4)
subset(CmeanPy_Modules, GeneID %in% TFdic[TFdic$TF.Name %in% c("MYB31","WRKY53"),]$TF.v4)

# make CornCyc list and remove small PWY == 1 
CornCYC  <- subset(CornCYC, GeneID %in% Syntenic)
CornCYC.list <- split(CornCYC$GeneID, CornCYC$Pathway.id)

CornCYC_size <-  as.data.frame(t(as.data.frame(lapply(CornCYC.list, length))))
colnames(CornCYC_size) <- "Freq"
head(CornCYC_size)
#
CornCYCred  <- subset(CornCYC, !(Pathway.id %in% row.names(subset(CornCYC_size, Freq == 1))))
CornCYC.list <- split(CornCYCred$GeneID, CornCYCred$Pathway.id)


hist(table(WGCNA_Modules$M))


Enrichmet_classes <- function(network){
  ## Count TF targets in network
  # Count Total
  network <- unique(network[,c("GeneID", "M")])
  #network <- subset(network, Target %in% Syntenic)
  #
  Total_targets <- as_tibble(as.data.frame(table(network$M)))
  colnames(Total_targets) <- c('Module', 'nGenes') 
  Total_targets$Module <- as.character(Total_targets$Module)
  
  
  # list input: network
  network.list <- split(network$GeneID, network$M)
  
  ## Compare list of predicted targets vs annotated genes in query.vector 
  #
  print(". Pre-newGOM .")
  go.obj <- newGOM(network.list, CornCYC.list, genome.size=nrow(WGCNA_Modules)) # annotated genes in Genome v4
  
  Pval <- getMatrix(go.obj, name="pval")
  Common <- getMatrix(go.obj, name="intersection")
  
  print(". Post-newGOM .")
  
  ### Summary tables
  ## adjust p value
  Pval <- as.data.frame(Pval)
  Pval[,1:ncol(Pval)] <- apply(Pval[,1:ncol(Pval)], 2, p.adjust)
  
  Pval_table   <- as_tibble(reshape2::melt(as.matrix(Pval)))
  colnames(Pval_table) <- c('Module', 'PWY', 'padj')
  Pval_table[,1:2] <- apply(Pval_table[,1:2], 2, as.character)
  
  #print(Pval_table)
  #
  Common_table <- as_tibble(reshape2::melt(as.matrix(Common)))
  colnames(Common_table) <- c('Module', 'PWY', 'n.targ')
  Common_table[,1:2] <- apply(Common_table[,1:2], 2, as.character)
  #print(Common_table)
  
  # Add predicted target in class by Module
  Pval_table <- left_join(Pval_table, Common_table , by=c('Module', 'PWY'))
  
  # Add total predicted targets
  Pval_table <- left_join(Pval_table, Total_targets, by="Module")
  
  # Select significant TFs 
  Pval_table <- subset(Pval_table, padj <= 0.1)
  #Pval_table <- tibble(TF=c("test", "a", "b"), TF2="test2")
  
  #write.table(Pval_table, name, sep = "\t", row.names = F, quote = F)
  print(paste("... Done ...", sep = ""))
  return(Pval_table)
}


# NetworkEnrichment <- Enrichmet_classes(subset(WGCNA_Modules, M %in% c(46,94))) # test modules with MB31 and WRKY53
# NetworkEnrichment_cmean <- Enrichmet_classes(subset(CmeanPy_Modules, M %in% c(23,21))) # test modules with MB31 and WRKY53
# NetworkEnrichment <- left_join(NetworkEnrichment, unique(CornCYC[,1:2]), by=c("PWY"="Pathway.id"))


table(subset(CmeanPy_Modules, M==23)$GeneID %in% subset(PDI, Source=="Zm00001d006236")$Source) # Zm00001d006236: MYB31; M 23
####################################################

####################################################
########           adjacence                ########
####################################################


library(Rtsne)

Matrix_freq <- function(TF_target_DF) {
  
  # Count TF-target associations
  TF_target_DF <- as_tibble(as.data.frame(table(TF_target_DF), stringsAsFactors = F))
  
  # filter 
  TF_target_DF <- subset(TF_target_DF, Freq > 0) # remove TF-target with zero freq to reduce speed during dcast
  TF_target_DF <- reshape2::dcast(TF_target_DF, Source ~ Target)
  row.names(TF_target_DF) <- TF_target_DF$Source
  TF_target_DF <- TF_target_DF[,-c(1)]
  # 
  TF_target_DF[is.na(TF_target_DF)] <- 0
  
  # scale values
  TF_target_DF_M <- apply(TF_target_DF, 2, scale) # scale by targets: Columns on this format
  row.names(TF_target_DF_M) <- row.names(TF_target_DF)
  print(dim(TF_target_DF_M) == dim(TF_target_DF))
  
  return(list(Freq=TF_target_DF, Z=TF_target_DF_M))
  
}

PDI_adj <- Matrix_freq(PDI)
CoExp_adj <- Matrix_freq(CoExp)
teQTL_adj <- Matrix_freq(teQTL)

All_adj <- Matrix_freq(teQTL)


tsne_PDI = Rtsne(as.matrix(t(PDI_adj$Z)), 
                 check_duplicates=FALSE, pca=TRUE, 
                 perplexity=30, theta=0.2, dims=2,
                   num_threads=40)

tsne_CoExp = Rtsne(as.matrix(t(CoExp_adj$Z)), 
                 check_duplicates=FALSE, pca=TRUE, 
                 perplexity=30, theta=0.2, dims=2,
                 num_threads=40)

tsne_teQTL = Rtsne(as.matrix(t(teQTL_adj$Z)), 
                 check_duplicates=FALSE, pca=TRUE, 
                 perplexity=30, theta=0.2, dims=2,
                 num_threads=40)


tsne_PDI = as.data.frame(tsne_PDI$Y)
row.names(tsne_PDI) <- colnames(PDI_adj$Z)

tsne_CoExp = as.data.frame(tsne_CoExp$Y)
row.names(tsne_CoExp) <- colnames(CoExp_adj$Z)

tsne_teQTL = as.data.frame(tsne_teQTL$Y)
row.names(tsne_teQTL) <- colnames(teQTL_adj$Z)

colnames(tsne_PDI) <- paste0("GRN_", colnames(tsne_PDI))
colnames(tsne_CoExp) <- paste0("CEN_", colnames(tsne_CoExp))
colnames(tsne_teQTL) <- paste0("GAN_", colnames(tsne_teQTL))


tsne_All <- as_tibble(data.frame(gid=Reduce(intersect, list(row.names(tsne_PDI), row.names(tsne_CoExp), row.names(tsne_teQTL)))))

tsne_All <- left_join(tsne_All, tsne_PDI %>% mutate(gid = rownames(tsne_PDI)), by="gid")
tsne_All <- left_join(tsne_All, tsne_CoExp %>% mutate(gid = rownames(tsne_CoExp)), by="gid")
tsne_All <- left_join(tsne_All, tsne_teQTL %>% mutate(gid = rownames(tsne_teQTL)), by="gid")

tsne_All[,2:7] <- apply(tsne_All[,2:7], 2, scale)

write.table(tsne_All, "tsne_All.txt", sep = "\t", quote = F, row.names = F)

row.names(peAll) <- peAll$gid
