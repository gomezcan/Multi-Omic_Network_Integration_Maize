library(tidyverse)
library(ggrepel)
library(ggpubr)
library(data.table)
library(GeneOverlap)
library(viridis)
library(ComplexHeatmap)


ReplaceName <- function(ids){
  
  for (i in 1:nrow(Top45)){
    ids <- gsub(Top45$V1[i], Top45$V2[i], ids)
  }
  
  for (i in 1:nrow(TFdic)){
    ids <- gsub(TFdic$TF.v4[i], TFdic$TF.Name[i], ids)
  }
  return(ids)
}

##################################################
##########        Annotations       ##############
##################################################

#### top 45
Top45 <- as_tibble(read.table("Data/Annotations/Top45.txt", h=F, stringsAsFactors = F))

# TF names
TFdic <- as_tibble(read.table("Data/Annotations/TF_Id_Name.txt", h=T, stringsAsFactors = F))

# Phenolic related genes
PheGenes <- as_tibble(read.table("Data/Annotations/LinaPheGenes2020.txt", h=T, sep = "\t", quote="", stringsAsFactors = F))

# TF and CoReg
TF_CoR <- as_tibble(read.table("Data/Annotations/TF_CoR_Mazie.txt", h=T, stringsAsFactors = F)) 

## Y1H network
Y1H <- as_tibble(read.table("Data/Annotations/Y1H.data.txt", h=T, stringsAsFactors = F))[,4:6]

## syntenic genes 
Syntenic <- as_tibble(read.table("Data/Annotations/Zm.v4.synteny.genes.txt", h=T, stringsAsFactors = F, sep = '\t'))
Syntenic <- subset(Syntenic, gene_synteny=='syntenic')$gene_id

# CornCYC
CornCYC <- as_tibble(read.table("Data/Annotations/corn_pathways.0210325.reduced.txt", h=T, sep = '\t'))


####
ciseqtls <- as_tibble(read.table("../Fig_ciseQTL/cis_eQTL.pdi.network.txt"))
colnames(ciseqtls) <- c('chr', 'snp', 'snp.id', 'target', 'source', 'summit', 'TF', 'snp_summit')
#
ciseqtls[,"TF.Name"] <- ciseqtls$TF
ciseqtls$TF <- gsub("_deM", "", ciseqtls$TF) 
ciseqtls$TF <- gsub("_Met", "", ciseqtls$TF)
ciseqtls

# Remove non-syntenic genes 
ciseqtls_full <- ciseqtls
ciseqtls <- subset(ciseqtls, target %in% Syntenic)


# define query list
List_query <- rbind(tibble(Target="Phe", GeneID=unique(PheGenes$GeneID)),
                    tibble(Target="TFs", GeneID=unique(TF_CoR$GeneID)), 
                    tibble(Target="TF.phe", GeneID=unique(Y1H$TF.v4)[is.na(unique(Y1H$TF.v4)) == FALSE]))

List_query <- subset(List_query, GeneID %in% Syntenic) # reduce analysis to only Syntenic genes
List_query <- split(List_query$GeneID, List_query$Target)


#
# replace name by gene ID
for (i in 1:nrow(TFdic)){
  ciseqtls$TF <- gsub(TFdic$TF.Name[i], TFdic$TF.v4[i], ciseqtls$TF)
}

ciseqtls$TF <- gsub('bZIZm00001d0288503', 'Zm00001d0288503', ciseqtls$TF)
ciseqtls$source <- as.character(ciseqtls$source)




# list cis_eQTLs by tissue
ciseqtls_list <- list()

for (i in unique(ciseqtls$source)) {
  ciseqtls_list[[i]] <- unique(subset(ciseqtls, source==i)[,c(4,5,7)])
}

Enrichmet.targets.eqTLs <- function(network){
  ## Count TF targets in network
  # Count Total
  Total_targtes <- as_tibble(as.data.frame(table(network$TF))) # mark unique targets
  colnames(Total_targtes) <- c('TF', 'targets') 
  Total_targtes$TF <- as.character(Total_targtes$TF)
  print(Total_targtes)
  
  # list input: network
  #print(network)
  network.list <- split(as.character(network$target), network$TF)
  
  
  ## Compare list of predicted targets vs annoated genes in query.vector 
  #
  go.obj <- newGOM(network.list, List_query, genome.size=length(Syntenic)) # annotated genes in Genome v4
  
  Pval <- getMatrix(go.obj, name="pval")
  print(Pval)
  Common <- getMatrix(go.obj, name="intersection")
  print(Common)
    
  ### Summary tables
  ## adjust p value
  Pval <- as.data.frame(Pval)
  Pval$Phe <- p.adjust(Pval$Phe, method = 'fdr')
  Pval$TF.phe <- p.adjust(Pval$TF.phe, method = 'fdr')
  Pval$TFs <- p.adjust(Pval$TFs, method = 'fdr')
  
  Pval_table   <- as_tibble(reshape::melt(as.matrix(Pval))) 
  colnames(Pval_table) <- c('TF', 'class', 'padj')
  Pval_table[,1:2] <- apply(Pval_table[,1:2], 2, as.character)
  # 
  Common_table <- as_tibble(reshape::melt(as.matrix(Common)))
  colnames(Common_table) <- c('TF', 'class', 'n.targ')
  Common_table[,1:2] <- apply(Common_table[,1:2], 2, as.character)
  
  # Add predicted target in class by TF
  
  Pval_table <- left_join(Pval_table, Common_table , by=c('TF', 'class'))
  
  # Add total predicted targets
  Pval_table <- left_join(Pval_table, Total_targtes, by="TF")
  
  # Select significant TFs 
  Pval_table <- subset(Pval_table, padj <= 0.1)
  
  return(Pval_table)
}

ciseqtls_total <- unique(ciseqtls[,c("target", "TF")])

Enrichmet.targets.eqTLs(ciseqtls_total)

Enrichmet.targets.eqTLs(ciseqtls_list$Kern)

TF.Candidates.ceQTLs <- lapply(ciseqtls_list, Enrichmet.targets.eqTLs)

TF.Candidates.ceQTLs$seedling

###
top.phe <- c("Zm00001d047017", "Zm00001d006236", "Zm00001d021019")
ciseqtls.candidates <- unique(subset(ciseqtls, TF %in% top.phe & target %in% PheGenes$GeneID)[,c(4,5,7,9)])
ciseqtls.candidates$TF.Name <- ReplaceName(ciseqtls.candidates$TF.Name)

write.table(ciseqtls.candidates, "Top.Phe.targets.ceQTL.txt", sep = "\t", quote = F, row.names = F)

#
TF.Candidates.ceQTLs_table <- as_tibble(rbindlist(TF.Candidates.ceQTLs, idcol = TRUE))
TF.Candidates.ceQTLs_table$class <- as.character(TF.Candidates.ceQTLs_table$class)
colnames(TF.Candidates.ceQTLs_table)[1] <- 'Network'
unique(TF.Candidates.ceQTLs_table$Network)

length(unique(TF.Candidates.ceQTLs_table$TF))
subset(Top45, V1 %in% unique(TF.Candidates.ceQTLs_table$TF))[,2]


subset(TF.Candidates, TF %in% CoExp.phe.cand & class=='Phe' & source=='DAP_New')[,1]
subset(TF.Candidates, TF %in% CoExp.TFphe.cand & class=='TF.phe' & source=='DAP_New')[,1]

subset(TF.Candidates, TF %in% unique(TF.Candidates.ceQTLs_table$TF))


table(unique(TF.Candidates$TF) %in% unique(TF.Candidates.ceQTLs_table$TF))
table(unique(CoExp.phe.cand) %in% unique(TF.Candidates.ceQTLs_table$TF))
CoExp.phe.cand[unique(CoExp.phe.cand) %in% unique(TF.Candidates.ceQTLs_table$TF)] %in%
  
  
  #### Heatmap
  Matrix <- subset(TF.Candidates.ceQTLs_table, class=='TF.phe')[,c('Network', 'TF', 'padj')]
Matrix <- reshape2::dcast(Matrix, Network ~ TF, value.var="padj")
row.names(Matrix) <- as.character(Matrix[,1])
Matrix <- Matrix[,-c(1)]
Matrix[is.na(Matrix)] <- 1
Matrix <- t(-log10(as.matrix(Matrix)))

Hetmap_Matrix <- Heatmap(Matrix, name="-log10 (FDR)",
                         # column_km = 3,
                         # column_names_rot = 90,
                         # row_names_rot = 45,
                         cluster_rows = TRUE, cluster_columns = FALSE,
                         show_column_dend = FALSE, show_row_dend = FALSE, 
                         clustering_method_columns = "ward.D2",
                         col=viridis(50, direction = -1, option = "A"),
                         column_names_gp = gpar(fontsize = 10),
                         row_names_gp = gpar(fontsize = 4),
                         show_heatmap_legend = T,
                         # heatmap_ = unit(10),
                         heatmap_height = unit(8, 'cm'),
                         heatmap_width  = unit(6, 'cm'),
                         heatmap_legend_param = list(title_gp = gpar(fontsize = 10), 
                                                     labels_gp = gpar(fontsize = 10), 
                                                     direction = "horizontal"))
draw(Hetmap_Matrix, heatmap_legend_side = "bottom")
