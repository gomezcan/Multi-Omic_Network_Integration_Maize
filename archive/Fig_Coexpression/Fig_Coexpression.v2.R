library(tidyverse)
library(ggrepel)
library(ggpubr)
library(data.table)
library(GeneOverlap)
library(viridis)
#library(ComplexHeatmap)

#BiocManager::install("ComplexHeatmap")

##################################################
##########          Functions       ##############
##################################################

ReplaceName <- function(ids){
  
  for (i in 1:nrow(Top45)){
    ids <- gsub(Top45$V1[i], Top45$V2[i], ids)
  }
  
  for (i in 1:nrow(TFdic)){
    ids <- gsub(TFdic$TF.v4[i], TFdic$TF.Name[i], ids)
  }
  return(ids)
}

Enrichmet.targets <- function(network){
  ## Count TF targets in network
  # Count Total
  network <- subset(network, tgt.gid %in% Syntenic) 
  #
  Total_targtes <- as_tibble(as.data.frame(table(unique(network[,1:2])$reg.gid)))
  colnames(Total_targtes) <- c('TF', 'targets') 
  Total_targtes$TF <- as.character(Total_targtes$TF)
  
  # Count positives Positive correlations
  net.pos <- subset(network, pcc > 0)
  net.pos <- as_tibble(as.data.frame(table(unique(net.pos[,1:2])$reg.gid)))
  colnames(net.pos) <- c('TF', 'targ.pos') 
  net.pos$TF <- as.character(net.pos$TF)
  
  # Total 
  Total_targtes <- left_join(Total_targtes, net.pos, by='TF')
  Total_targtes['targ.neg'] <- Total_targtes$targets - Total_targtes$targ.pos 
  
  # list input: network
  network.list <- unique(network[,1:2])
  network.list <- split(network.list$tgt.gid, network.list$reg.gid)
  
  ## Compare list of predicted targets vs annoated genes in query.vector 
  #
  go.obj <- newGOM(network.list, List_query, genome.size=length(Syntenic)) # annotated genes in Genome v4
  
  Pval <- getMatrix(go.obj, name="pval")
  Common <- getMatrix(go.obj, name="intersection")
  
  ### Summary tables
  ## adjust p value
  Pval <- as.data.frame(Pval)
  Pval[,1:ncol(Pval)] <- apply(Pval[,1:ncol(Pval)], 2, p.adjust)
  
  Pval_table   <- as_tibble(reshape::melt(as.matrix(Pval))) 
  colnames(Pval_table) <- c('TF', 'class', 'padj')
  
  
  # 
  Common_table <- as_tibble(reshape::melt(as.matrix(Common)))
  colnames(Common_table) <- c('TF', 'class', 'n.targ')
  
  # Add predicted target in class by TF
  Pval_table <- left_join(Pval_table, Common_table , by=c('TF', 'class'))
  
  # Add total predicted targets
  Pval_table <- left_join(Pval_table, Total_targtes, by="TF")
  
  # Select significant TFs 
  Pval_table <- subset(Pval_table, padj <= 0.1)
  
  return(Pval_table)
}


outdegree.counter <- function(network){
  ## Count TF targets in network
  # Count Total
  network <- subset(network, tgt.gid %in% Syntenic) 
  #
  Total_targtes <- as_tibble(as.data.frame(table(unique(network[,1:2])$reg.gid)))
  colnames(Total_targtes) <- c('TF', 'targets') 
  Total_targtes$TF <- as.character(Total_targtes$TF)
  
  
  return(Total_targtes)
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
CornCYC$Gene.id <- gsub('ZM', 'Zm', CornCYC$Gene.id)
CornCYC$Gene.id <- gsub('D', 'd', CornCYC$Gene.id)
CornCYC <- subset(CornCYC, Gene.id %in% Syntenic)

# Coexpression data
# CoexDB <- as_tibble(read.table("ExpressionData/widiv304.grn.tsv.gz", h=T))
# CoexDB <- CoexDB[order(-CoexDB$score),][1:100000,]

CoexDB = readRDS("Data/Coexpression/rf.100k.rds")
#  add network names
names(CoexDB$tn) <- CoexDB$nid
names(CoexDB)
CoexDB <- CoexDB$tn
CoexDB$n16b

# Coexpression results from 304 lines used on eQTL analysis
CoexDB_304 <- as_tibble(read.table("Data/Coexpression/widiv304.grn.tsv.gz", h=T))
CoexDB_304 <- CoexDB_304[order(-CoexDB_304$score),]
CoexDB_304 <- CoexDB_304[1:100000,] # select top predictions

# add 
CoexDB[['n304']] <- CoexDB_304


# define query list
List_query <- rbind(tibble(Target="Phe", GeneID=unique(PheGenes$GeneID)),
                    tibble(Target="TFs", GeneID=unique(TF_CoR$GeneID)), 
                    tibble(Target="TF.phe", GeneID=unique(Y1H$TF.v4)[is.na(unique(Y1H$TF.v4)) == FALSE]))

List_query <- split(List_query$GeneID, List_query$Target)
CornCYC.list <- split(CornCYC$Gene.id, CornCYC$Pathway.id)
#
List_query <- c(List_query, CornCYC.list)




outdegree.counter


######################################################################################
######### Count number of interactions for each class genes set of interest  #########
######################################################################################

#TF.Candidates <- Enrichmet.targets(CoexDB$n16b)
TF.Candidates <- lapply(CoexDB, Enrichmet.targets)

TF.Candidates_table <- as_tibble(rbindlist(TF.Candidates, idcol = TRUE))
colnames(TF.Candidates_table)[1] <- 'Network'

write.table(TF.Candidates_table, "CoExpression_candidates.fullReport.06-10-21.txt", sep = '\t', row.names = F, quote = F)


top.phe <- c("Zm00001d047017", "Zm00001d006236", "Zm00001d021019")
CoexDB.candiates <- as_tibble(rbindlist(CoexDB, idcol = TRUE))

CoexDB.candiates.HSF24 <- subset(CoexDB.candiates, reg.gid=="Zm00001d032923")
View(as.data.frame(table(CoexDB.candiates.HSF24$.id)))


#subset(TF.Candidates, TF=="Zm00001d032923")
CoexDB.candiates <- subset(CoexDB.candiates, reg.gid %in% top.phe & tgt.gid %in% PheGenes$GeneID)
write.table(CoexDB.candiates, "Top.Phe.targets.CoExp.txt", sep = "\t", quote = F, row.names = F)


########  matrix input to heatmap
## Phe genes

MakeHeamap <- function(geneset, r.c.val.vector){
  
  Matrix <- subset(TF.Candidates_table, class==geneset)[,r.c.val.vector]
  print(head(Matrix))
  print(dim(Matrix))
   
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
                                    column_names_gp = gpar(fontsize = 2),
                                    row_names_gp = gpar(fontsize = 0.5),
                                    show_heatmap_legend = T,
                                    # heatmap_ = unit(10),
                                    heatmap_height = unit(8, 'cm'),
                                    heatmap_width  = unit(6, 'cm'),
                                    heatmap_legend_param = list(title_gp = gpar(fontsize = 10), 
                                                                labels_gp = gpar(fontsize = 10), 
                                                                direction = "horizontal"))
  return(Hetmap_Matrix)
}

TF.Candidates_table.Summary <- as_tibble(as.data.frame(table(TF.Candidates_table[,2:3])))


write.table(subset(TF.Candidates_table.Summary, Freq >0), "CoExpression_candidates.06-02-21.txt", sep = '\t', row.names = F, quote = F)

#length(unique(subset(TF.Candidates_table.Summary, class == 'TF.phe' & Freq>0)$TF))

Plot.summary.phe <- ggplot(subset(TF.Candidates_table.Summary, class == 'Phe' & Freq>0), 
                       aes(y=reorder(TF, Freq), x=class, fill=Freq))+
  geom_tile()  +
  #scale_x_continuous(breaks = seq(0,100,20), expand = c(0,0)) + 
  scale_fill_viridis(option = "D", direction = 1) + #
  theme_bw() +
  theme(legend.title = element_text(vjust = 0.8), 
        panel.border = element_blank(), 
        panel.grid.major = element_blank())

Plot.summary.TFphe <- ggplot(subset(TF.Candidates_table.Summary, class == 'TF.phe' & Freq>0), 
                           aes(y=reorder(TF, Freq), x=class, fill=Freq))+
  geom_tile()  +
  #scale_x_continuous(breaks = seq(0,100,20), expand = c(0,0)) + 
  scale_fill_viridis(option = "D", direction = 1) + #
  theme_bw() +
  theme(legend.title = element_text(vjust = 0.8), 
        panel.border = element_blank(), 
        panel.grid.major = element_blank())


Plot.summary.phe <- ggpar(Plot.summary.phe, font.ytickslab = 0.1, xlab = "", ylab='TFs')
Plot.summary.TFphe <- ggpar(Plot.summary.TFphe, font.ytickslab = 0.1, xlab = "", ylab='TFs')


# heatmap TF enriched on phe and TF-phe genes
Candidates.phe <- unique(as.character(subset(TF.Candidates_table.Summary, class == 'Phe' & Freq>0)$TF))
Candidates.TFphe <- unique(as.character(subset(TF.Candidates_table.Summary, class == 'TF.phe' & Freq>0)$TF))

table(unique(Y1H$TF.v4) %in% Candidates.phe)
table(unique(Y1H$TF.v4) %in% Candidates.TFphe)
#
Both_TFs <- subset(TF.Candidates_table.Summary, class == 'TF.phe' & Freq>0)
Both_TFs <- as.character(subset(Both_TFs, TF %in% Candidates.phe)$TF)
#
subset(Top45, V1 %in% Candidates.phe)
subset(Top45, V1 %in% Candidates.TFphe)
subset(Top45, V1 %in% Both_TFs)

Plot.summary.phe_TFphe <- ggplot(subset(TF.Candidates_table.Summary, class != 'TFs' & TF %in% Both_TFs), 
                             aes(y=reorder(TF, Freq), x=class, fill=Freq))+
  geom_tile()  +
  #scale_x_continuous(breaks = seq(0,100,20), expand = c(0,0)) + 
  scale_fill_viridis(option = "D", direction = 1) + #
  theme_bw() +
  theme(legend.title = element_text(vjust = 0.8), 
        panel.border = element_blank(), 
        panel.grid.major = element_blank())

Plot.summary.phe_TFphe <- ggpar(Plot.summary.phe_TFphe, font.ytickslab = 0.1, xlab = "", ylab='TFs')

#3x5
Hetmap_Matrix_Phe.pval <- MakeHeamap('Phe', c('Network','TF', 'padj'))
Hetmap_Matrix_TF.Phe.pval <- MakeHeamap('TF.phe', c('Network','TF', 'padj'))
Hetmap_Matrix_Phe.pval <- MakeHeamap('TFs', c('Network','TF', 'padj'))

draw(Hetmap_Matrix_Phe.pval, heatmap_legend_side = "bottom")
draw(Hetmap_Matrix_TF.Phe.pval, heatmap_legend_side = "bottom")
draw(Hetmap_Matrix_Phe.pval, heatmap_legend_side = "bottom")

## Phe genes

Matrix_TF.Phe.pval <- reshape2::dcast(subset(TF.Candidates_table, class=='TF.phe')[,c(1,2,4)],
                                   Network ~ TF, value.var="padj")
row.names(Matrix_TF.Phe.pval) <- Matrix_TF.Phe.pval$Network
Matrix_TF.Phe.pval <- Matrix_TF.Phe.pval[,-c(1)]
Matrix_TF.Phe.pval[is.na(Matrix_TF.Phe.pval)] <- 1

Input2 <- t(-log10(as.matrix(Matrix_TF.Phe.pval)))



# TFs: trans-eQTLs mapping to TF, which could have as targets either TF or Phe genes
Trans_eQTL_TF <- paste(subset(teQTLs, isTF==TRUE)$Annotation, subset(teQTLs, isTF==TRUE)$Target, sep = "_")
Trans_eQTL_Phe <- paste(subset(teQTLs, isTF==TRUE & Target %in% PheGenes$GeneID )$Annotation, 
                        subset(teQTLs, isTF==TRUE & Target %in% PheGenes$GeneID)$Target, sep = "_")

All.Trans_eQTL_Phe <- paste(subset(teQTLs, Target %in% PheGenes$GeneID)$Annotation, 
                            subset(teQTLs, Target %in% PheGenes$GeneID)$Target, sep = "_")

subset(as.data.frame(table(subset(teQTLs, Target %in% PheGenes$GeneID)$Annotation)), Freq>2)

# Total
Trans_eQTL_Gene <- as_tibble(as.data.frame(table(Trans_eQTL_Gene)))
# TFs: trans-eQTLs mapping to TF, which could have as targets either TF or Phe genes
Trans_eQTL_TF <- as_tibble(as.data.frame(table(Trans_eQTL_TF)))
# Phe Genes
Trans_eQTL_Phe <- as_tibble(as.data.frame(table(Trans_eQTL_Phe)))

teQTLs_In_Phe_Plot <- ggplot(as_tibble(as.data.frame(table(Trans_eQTL_Phe$Freq))),
                             aes(x=Var1, y=Freq))+
  geom_bar(stat="identity", fill="hotpink") +
  geom_text(aes(label=Freq), vjust=-0.3, size=3.5) +
  theme_pubclean() 

# Multiple eQTLs by Phe Genes
teQTLs_In_Phe_Plot <- ggpar(teQTLs_In_Phe_Plot,  ylab = "TFs/CoR",  
                            font.ytickslab = 12, font.xtickslab = 12, 
                            font.x = 14, font.y = 14, xlab = "eQTL")
teQTLs_In_Phe_Plot

# Multiple eQTLs by TF Genes: regulators of regulators
teQTLs_In_TF_Plot <- ggplot(as_tibble(as.data.frame(table(Trans_eQTL_TF$Freq))),
                            aes(x=Var1, y=Freq))+
  geom_bar(stat="identity", fill="#56B4E9")+
  geom_text(aes(label=Freq), vjust=-0.3, size=2) +
  theme_pubclean()

teQTLs_In_TF_Plot <- ggpar(teQTLs_In_TF_Plot,  ylab = "TFs/CoR",  
                           font.ytickslab = 12, font.xtickslab = 12, 
                           font.x = 14, font.y = 14, xlab = "eQTL")

ggarrange(teQTLs_In_Phe_Plot, teQTLs_In_TF_Plot, ncol = 1) 
######################################################################################

######################################################################################
################  Count trans eQTLs: Total Potential targets         #################
######################################################################################

Trans_eQTL_Target_Gene <- as_tibble(as.data.frame(table(unique(teQTLs[,c(2,5)])$Annotation)))

# trans-eQTLs mapping to TF, which could have as targets either TF or Phe genes
# Total
Trans_eQTL_TF_total.Target <- as_tibble(as.data.frame(table(unique(subset(teQTLs, isTF==TRUE)[,c(2,5)])$Annotation)))
# Phe                       
Trans_eQTL_Phe_Target <- as_tibble(as.data.frame(table(unique(subset(teQTLs, isTF==TRUE & Target %in% PheGenes$GeneID)[,c(2,5)])$Annotation))) 
# TFs
Trans_eQTL_TF_Target <- as_tibble(as.data.frame(table(unique(subset(teQTLs, isTF==TRUE & IsTF_Target == TRUE)[,c(2,5)])$Annotation)))


###
teQTLs_In_PheTarget_Plot <- ggplot(as_tibble(as.data.frame(table(Trans_eQTL_Phe_Target$Freq))),
                                   aes(x=Var1, y=Freq))+
  geom_bar(stat="identity", fill="hotpink") + geom_text(aes(label=Freq), vjust=-0.3, size=3.5) +
  theme_pubclean()

teQTLs_In_PheTarget_Plot <- ggpar(teQTLs_In_PheTarget_Plot,  ylab = "TFs/CoR",  
                                  font.ytickslab = 12, font.xtickslab = 12, 
                                  font.x = 14, font.y = 14, xlab = "Targets")

teQTLs_In_TFsTarget_Plot <- ggplot(as_tibble(as.data.frame(table(Trans_eQTL_TF_Target$Freq))),
                                   aes(x=Var1, y=Freq))+
  geom_bar(stat="identity", fill="#56B4E9") + geom_text(aes(label=Freq), vjust=-0.3, size=3.5) +
  theme_pubclean()
teQTLs_In_TFsTarget_Plot <- ggpar(teQTLs_In_TFsTarget_Plot,  ylab = "TFs/CoR",  
                                  font.ytickslab = 12, font.xtickslab = 12, 
                                  font.x = 14, font.y = 14, xlab = "Targets")

teQTLs_In_TFsTarget_Plot

ggarrange(teQTLs_In_PheTarget_Plot, teQTLs_In_TFsTarget_Plot, ncol =  1, align = 'h')

# Top values
subset(Trans_eQTL_Phe, Freq >= 10)         # eQTLs supporting an association
subset(Trans_eQTL_Phe_Target, Freq >= 4)  # TF-"target" associations
# Summary by unique TFs-target associations

#subset(Trans_eQTL_Phe, Trans_eQTL_Phe, Freq >= 15)         # eQTLs supporting an association
subset(Trans_eQTL_Phe_Target, Var1 %in% Top45$V1)  # TF-"target" associations


#####################################################################################

length(unique(subset(teQTLs, isTF==TRUE & IsTF_Target != TRUE)$Target))

# trans_eQTL_TF_Table <- as_tibble(as.data.frame(table(unique(subset(teQTLs, isTF==TRUE)[,c(2,5)])$Annotation)))
# trans_eQTL_TF.with_TF.target_Table <- as_tibble(as.data.frame(
#   table(unique(subset(teQTLs, isTF==TRUE & IsTF_Target==TRUE)[,c(2,5)])$Annotation)))
# 
# trans_eQTL_TF.with_Phen.target_Table <- as_tibble(as.data.frame(
#   table(unique(subset(teQTLs, isTF==TRUE & IsTF_Target==FALSE)[,c(2,5)])$Annotation)))
# #
# trans_eQTL_TF.with_TF.target_Table <- trans_eQTL_TF.with_TF.target_Table[order(-trans_eQTL_TF.with_TF.target_Table$Freq),]
# trans_eQTL_TF.with_Phen.target_Table <- trans_eQTL_TF.with_Phen.target_Table[order(-trans_eQTL_TF.with_Phen.target_Table$Freq),]
# 
# teQTLs_TFs_plot <- ggplot(as_tibble(as.data.frame(table(subset(trans_eQTL_TF.with_Phen.target_Table, Freq>0)$Freq))),
#                           aes(x=Var1, y=Freq))+
#   geom_bar(stat="identity", fill="#56B4E9") +
#   geom_text(aes(label=Freq), vjust=-0.3, size=3.5) +
#   theme_pubclean() 
#   
# teQTLs_TFs_plot <- ggpar(teQTLs_TFs_plot, ylab = "TFs", xlab = "Number of trans-eQTLs")


#############################################################
#######     Define Unique trans-QTLs association     ########
#############################################################

teQTLs["IsPheG"] <- teQTLs$Target %in% PheGenes$GeneID
teQTLs["IsAnnotationPheG"] <- teQTLs$Annotation %in% PheGenes$GeneID

write.table(teQTLs, "eQTL_Results/Total.transeQTLs.txt",  row.names = F, sep = '\t', quote = F)
write.table(subset(teQTLs, isTF==TRUE), "eQTL_Results/TF-Gene.transeQTLs.txt",  row.names = F, sep = '\t', quote = F)
write.table(subset(teQTLs, isTF==TRUE & IsTF_Target==TRUE), "eQTL_Results/TF-TF.transeQTLs.txt",  row.names = F, sep = '\t', quote = F)
write.table(subset(teQTLs, Target %in% PheGenes$GeneID), 
            "eQTL_Results/All-PhenG.transeQTLs.txt",  row.names = F, sep = '\t', quote = F)

subset(teQTLs, isTF==TRUE & Target %in% PheGenes$GeneID)

edged_teQTLs.Total <- unique(subset(teQTLs, isTF==TRUE & IsTF_Target==FALSE)[,c("Annotation","Target")])
edged_teQTLs_TF_TF <- unique(subset(teQTLs, isTF==TRUE & IsTF_Target==TRUE)[,c("Annotation","Target")])

edged_teQTLs <- unique(subset(teQTLs, isTF==TRUE & (Target %in% PheGenes$GeneID))[,c("Annotation","Target")])

colnames(edged_teQTLs)[1] <- "GeneID"
colnames(edged_teQTLs_TF_TF)[1] <- "GeneID"
colnames(edged_teQTLs.Total)[1] <- "GeneID"

edged_teQTLs <- left_join(edged_teQTLs, TFdic, by="GeneID")
edged_teQTLs_TF_TF <- left_join(edged_teQTLs_TF_TF, TFdic, by="GeneID")
#############################################################


##################################################
######            cis-eQTLs               ######
##################################################

# Total before PDI data
cis_eQTLs_pass.Total <-  Filter_eQTLs <- subset(cis_eQTLs_pass, snp_chr==gene_chr)
cis_eQTLs_pass.Total <- subset(cis_eQTLs_pass.Total, (abs((snp_pos-gene_start)) <= 100000) | (abs((snp_pos-gene_stop)) <= 100000))
cis_eQTLs_pass.Total <- cis_eQTLs_pass.Total[,c(1:3,10,12,13)]

write.table(cis_eQTLs_pass.Total, "eQTL_Results/Total.beforePDI.ciseQTLs.txt", row.names = F, quote = F, sep = '\t')

# cis annotated to Phe genes
cis_eQTLs_pass_Phe <-  subset(cis_eQTLs_pass, Target %in% as.character(PheGenes$GeneID))
cis_eQTLs_pass_Phe <-  subset(cis_eQTLs_pass_Phe, snp_chr==gene_chr)
cis_eQTLs_pass_Phe <- subset(cis_eQTLs_pass_Phe, (abs((snp_pos-gene_start)) <= 100000) | (abs((snp_pos-gene_stop)) <= 100000))


# eQTLs and its distance to peaks summits
ceQTLs <- subset(ceQTLs, V9 !='.')[,c(4,8,9)]      # remove cis eQTL no mapping to peaks
ceQTLs[,c(1:2)] <- apply(ceQTLs[,c(1:2)], 2, as.character)
colnames(ceQTLs) <- c("snp", "Annotation", "Dis")
ceQTLs$Annotation <- gsub("_deM", "", ceQTLs$Annotation)
ceQTLs$Annotation<- gsub("_Met", "", ceQTLs$Annotation)
#

### Filter eQTLs that are on the same Chr
# 1. snp-target same chr
colnames(eQTLs)
Filter_eQTLs <- subset(cis_eQTLs_pass, snp_chr==gene_chr)[,c(1,3,10,12,13)]
# 2. snp-target dis <= 100 kbs
Filter_eQTLs <- subset(Filter_eQTLs, 
                       (abs((snp_pos-gene_start)) <= 100000) | 
                         (abs((snp_pos-gene_stop)) <= 100000))
length(unique(Filter_eQTLs$snp))
colnames(Filter_eQTLs)[3] <- "Target"

## Filter eQTLs based on distances betwen SNP and peak summit
## 1. Get ceQTLs that pass filter and that are in close proximity to a peak summit ##
ceQTLs_Pass <- subset(ceQTLs, snp %in% Filter_eQTLs$snp)
## 2. get eQTL that in close proximity to a peak summit 
ceQTLs_Pass <- subset(ceQTLs_Pass, abs(Dis) <=20)
length(unique(ceQTLs_Pass$snp))

########## Dis plot  ##########
###############################
ceQTLs_Dis_Freq <- as_tibble(as.data.frame(table(ceQTLs_Pass$Dis)))
ceQTLs_Dis_Freq$Var1 <- (as.numeric(as.vector(ceQTLs_Dis_Freq$Var1)))
ceQTLs_Dis_Freq
# Plot distances distribution
ceQTLs_Dis_plot <- ggplot(ceQTLs_Dis_Freq, aes(x=Var1, y=Freq))+
  geom_bar(stat="identity", fill="#56B4E9") +
  scale_x_continuous(limits = c(-30,30)) +
  #geom_text(aes(label=Freq), vjust=-0.3, size=2) +
  theme_pubclean()
ceQTLs_Dis_plot
###############################

########## snp mapped to multiple summits  ##########
#####################################################

ProblematicSNP_Summit <- subset(as_tibble(as.data.frame(table(ceQTLs_Pass$snp))), Freq>1)
ProblematicSNP_Summit_plot <- ggplot(as.data.frame(table(ProblematicSNP_Summit$Freq)), 
                                     aes(x=Var1, y=Freq))+
  geom_bar(stat="identity", fill="#56B4E9") +
  #scale_x_continuous(limits = c(-100,100)) +
  geom_text(aes(label=Freq), vjust=-0.3, size=5) +
  theme_pubclean()
ProblematicSNP_Summit_plot
ProblematicSNP_Summit_plot <- ggpar(ProblematicSNP_Summit_plot,  ylab = "Counts",  
                                    font.ytickslab = 12, font.xtickslab = 12, 
                                    font.x = 14, font.y = 14, xlab = "SNP-Summit Redundant")
ProblematicSNP_Summit_plot

subset(ceQTLs_Pass, snp %in% subset(ProblematicSNP_Summit, Freq>=7)$Var1)

#####################################################

# 3. Remove Problematic cis-eQTL: snp mapped to multiple summits
#ceQTLs_Pass2 <- subset(ceQTLs_Pass, !(snp %in% ProblematicSNP_Summit$Var1))
ceQTLs_Pass2 <- ceQTLs_Pass
ceQTLs_Pass2 <- left_join(ceQTLs_Pass2, Filter_eQTLs[,c("snp", "Target", "snp_pos")], by="snp")
ceQTLs_Pass2["isTargetTF"] <- ceQTLs_Pass2$Target %in% TF_CoR
colnames(ceQTLs_Pass2)[2] <- "TF"

ceQTLs_Pass2.total <- ceQTLs_Pass2
length(unique(ceQTLs_Pass2.total$snp))

ceQTLs_Pass2_TFs <- subset(ceQTLs_Pass2, isTargetTF==TRUE)
length(unique(ceQTLs_Pass2_TFs$snp))
ceQTLs_Pass2 <- subset(ceQTLs_Pass2, Target %in% PheGenes$GeneID)
length(unique(ceQTLs_Pass2$snp))


## Case 1: Miltiple SNP supporting the same association
# Phe
Case1 <- as_tibble(as.data.frame(table(paste(ceQTLs_Pass2$Annotation, ceQTLs_Pass2$Target, sep='_'))))
# TF
Case1.TF <- as_tibble(as.data.frame(table(paste(ceQTLs_Pass2_TFs$Annotation, ceQTLs_Pass2_TFs$Target, sep='_'))))

## Case 2: Same SNP supporting the multiples association
# Phe
Case2 <- as_tibble(as.data.frame(table(ceQTLs_Pass2$snp)))
# TF
Case2.TF <- as_tibble(as.data.frame(table(ceQTLs_Pass2_TFs$snp)))
## Case 3: Same TFs associated woth multiples association
# Phe
Case3 <- as_tibble(as.data.frame(table(unique(ceQTLs_Pass2[,c("Annotation", "Target")])$Annotation)))
# TF
Case3.TF <- as_tibble(as.data.frame(table(unique(ceQTLs_Pass2_TFs[,c("Annotation", "Target")])$Annotation)))

Case1_plot <- ggplot(as.data.frame(table(Case1$Freq)),aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity", fill="hotpink") +
  geom_text(aes(label=Freq), vjust=+0.1, size=5) +
  scale_y_continuous(limits = c(0,300)) +
  theme_pubclean()

Case1_TF_plot <- ggplot(as.data.frame(table(Case1.TF$Freq)),aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity", fill="#56B4E9") +
  geom_text(aes(label=Freq), vjust=+0.1, size=5) +
  scale_y_continuous(limits = c(0,600)) +
  theme_pubclean()

Case2_plot <- ggplot(as.data.frame(table(Case2$Freq)),aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity", fill="hotpink") +
  geom_text(aes(label=Freq), vjust=+0.1, size=5) +
  scale_y_continuous(limits = c(0,300)) +
  theme_pubclean()

Case2_TF_plot <- ggplot(as.data.frame(table(Case2.TF$Freq)),aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity", fill="#56B4E9") +
  geom_text(aes(label=Freq), vjust=+0.1, size=5) +
  scale_y_continuous(limits = c(0,600)) +
  theme_pubclean()

Case3_plot <- ggplot(as.data.frame(table(Case3$Freq)),aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity", fill="hotpink") +
  geom_text(aes(label=Freq), vjust=+0.1, size=5) +
  theme_pubclean()

Case3_TF_plot <- ggplot(as.data.frame(table(Case3.TF$Freq)),aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity", fill="#56B4E9") +
  geom_text(aes(label=Freq), vjust=+0.1, size=5) +
  theme_pubclean()

Case1_TF_plot
Case2_TF_plot
Case3_TF_plot

Case1_plot <- ggpar(Case1_plot,  ylab = "TFs",  xlab = "cis-eQTL", font.ytickslab = 12, font.xtickslab = 12, font.x = 14, font.y = 14)
Case2_plot <- ggpar(Case2_plot,  ylab = "TFs",  xlab = "cis-eQTL", font.ytickslab = 12, font.xtickslab = 12,  font.x = 14, font.y = 14)
Case3_plot <- ggpar(Case3_plot,  ylab = "TFs",  xlab = "Targets", font.ytickslab = 12, font.xtickslab = 12,  font.x = 14, font.y = 14)

Case1_TF_plot <- ggpar(Case1_TF_plot,  ylab = "TFs",  xlab = "cis-eQTL", font.ytickslab = 12, font.xtickslab = 12, font.x = 14, font.y = 14)
Case2_TF_plot <- ggpar(Case2_TF_plot,  ylab = "TFs",  xlab = "cis-eQTL", font.ytickslab = 12, font.xtickslab = 12,  font.x = 14, font.y = 14)
Case3_TF_plot <- ggpar(Case3_TF_plot,  ylab = "TFs",  xlab = "Targets", font.ytickslab = 12, font.xtickslab = 12,  font.x = 14, font.y = 14)

ggarrange(ggarrange(Case1_plot, Case2_plot, Case3_plot, nrow = 3),
          ggarrange(Case1_TF_plot, Case2_TF_plot, Case3_TF_plot, nrow = 3),
          ncol = 2)

#
subset(Case1, Freq >= 2) # Case 1 examples
subset(ceQTLs_Pass2, Annotation=="EREB6" & Target=="Zm00001d020401")
subset(ceQTLs_Pass2, Annotation=="Zm00001d002762" & Target=="Zm00001d020401")
#
subset(Case2, Freq >= 2) # Case 2 examples
subset(ceQTLs_Pass2, snp=="4-192889967")
#
subset(Case3, Freq >= 6) # Case 2 examples
subset(ceQTLs_Pass2, Annotation %in% c("DOF23","EREB6"))

#############################################################
#######     Define Unique cis-QTLs association     ########
#############################################################

# save results
colnames(ceQTLs_Pass2.total)[2] <- "GeneID"
ceQTLs_Pass2.total <- left_join(ceQTLs_Pass2.total, TFdic, by="GeneID")
ceQTLs_Pass2.total <- ReplaceName.by.ID(ceQTLs_Pass2.total)


write.table(ceQTLs_Pass2.total, "eQTL_Results/Total.ciseQTLs.txt", row.names = F, quote = F, sep = '\t')
write.table(subset(ceQTLs_Pass2.total, isTargetTF==TRUE), "eQTL_Results/TF-TF.ciseQTLs.txt", row.names = F, quote = F, sep = '\t')

# ceQTLs_Pass2: phe targets
edged_ceQTLs <- unique(ceQTLs_Pass2[,c("Annotation", "Target")])
edged_ceQTLs.TF_TF <- unique(ceQTLs_Pass2_TFs[,c("Annotation", "Target")])

colnames(edged_ceQTLs)[1] <- 'GeneID'
colnames(edged_ceQTLs.TF_TF)[1] <- 'GeneID'

edged_ceQTLs <- left_join(edged_ceQTLs, TFdic, by="GeneID")
edged_ceQTLs.TF_TF <- left_join(edged_ceQTLs.TF_TF, TFdic, by="GeneID")

######## fix names and IDs ######
################################
##
ReplaceName.by.ID <- function(cis.edges){
  tem <- subset(cis.edges, is.na(TFname)==TRUE)
  tem2 <- subset(cis.edges, is.na(TFname)==FALSE)
  #
  tem$TFname <- tem$GeneID
  tem <- left_join(tem, TFdic, by="TFname")
  tem$GeneID.x <- tem$GeneID.y
  tem <- tem[,(colnames(tem) != "GeneID.y")]
  colnames(tem)[colnames(tem)=="GeneID.x"] <- "GeneID"
  
  cis.edges <- rbind(tem,tem2)
  return(cis.edges)
}

edged_ceQTLs <- ReplaceName.by.ID(edged_ceQTLs)
edged_ceQTLs.TF_TF <- ReplaceName.by.ID(edged_ceQTLs.TF_TF)
#################################


#############################################################
##############      cis and trans eQTLs         #############
#############################################################
## Phe
edged_ceQTLs["edged"] <- paste(edged_ceQTLs$GeneID, edged_ceQTLs$Target, sep="_")
edged_ceQTLs["Source"] <- "cis"
edged_teQTLs["edged"] <- paste(edged_teQTLs$GeneID, edged_teQTLs$Target, sep="_")
edged_teQTLs["Source"] <- "trans"

## TFs
edged_ceQTLs.TF_TF["edged"] <- paste(edged_ceQTLs.TF_TF$GeneID, edged_ceQTLs.TF_TF$Target, sep="_")
edged_ceQTLs.TF_TF["Source"] <- "cis"
edged_teQTLs_TF_TF["edged"] <- paste(edged_teQTLs_TF_TF$GeneID, edged_teQTLs_TF_TF$Target, sep="_")
edged_teQTLs_TF_TF["Source"] <- "trans"

edged_ceQTLs$Target <- as.character(edged_ceQTLs$Target)

split(edged_ceQTLs$Target, edged_ceQTLs$GeneID)

eQTLs_edges <- list(Phe_cis=edged_ceQTLs$edged, 
                    Phe_trans=edged_teQTLs$edged,
                    TF_cis=edged_ceQTLs.TF_TF$edged,
                    TF_trans=edged_teQTLs_TF_TF$edged)
eQTLs_edges

##################################################
######            Co-expression               ####
##################################################

MakeEgde.list <- function(df){
  edge <- paste(df$regulator, df$target, sep='_')  
}

# define egdes 
CoexDB["edge"] <- paste(CoexDB$reg.gid, CoexDB$tgt.gid, sep = "_")
CoexDB["edge"] <- paste(CoexDB$TF, CoexDB$Target, sep = "_")
CoexDB <- CoexDB[,-c(8)]

CoexDB_15k <- CoexDB[1:15000,]
CoexDB_PCC0.5 <- subset(CoexDB, abs(pcc) >=0.5)



# set list of edges
CoexDB_45.list <- lapply(CoexDB_45, MakeEgde.list)
CoexDB_45.list[["n304"]] <- CoexDB$egde


TF.Coex <- unique(bind_rows(CoexDB_45)$regulator)
TF.Coex.304 <- unique(c(as.character(CoexDB$reg.gid), as.character(CoexDB$reg.gid)))
TF.eQTLs <- unique(c(edged_ceQTLs$GeneID, 
                     edged_teQTLs$GeneID, 
                     edged_ceQTLs.TF_TF$GeneID, 
                     edged_teQTLs_TF_TF$GeneID))



CountCommon.edges <- function(list.tf){
  # Count common target and test if significant
  go.obj <- newGOM(list.tf, CoexDB_45.list, genome.size=39591*2785) # all vs all 
  Pval <- getMatrix(go.obj, name="pval")
  J <- getMatrix(go.obj, name="Jaccard")
  common <- getNestedList(go.obj, name="intersection")
  return(list(M.pval=Pval, J=J, CommonEg=common))
}


CommonEgdes_CoExpression <- CountCommon.edges(eQTLs_edges)
#CommonTagets$M.pval[CommonTagets$M.pval==0] <- min(CommonTagets$M.pval[CommonTagets$M.pval!=min(CommonTagets$M.pval)])

m.text <- as.data.frame(lapply(CommonEgdes_CoExpression$CommonEg, sapply, length))
m.text <- m.text[, colnames(CommonTagets$M.pval)]

# 
library(ComplexHeatmap)
library(viridis)

Heatmap(as.matrix(-log10(CommonEgdes_CoExpression$M.pval)), 
        heatmap_width = unit(20, "cm"), 
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.f", m.text[i, j]), x, y, gp = gpar(fontsize = 6, col = "gray60"))
        },
        name = "-log10(P.val)", 
        cluster_columns = TRUE, show_column_dend = TRUE,
        show_row_dend = TRUE, column_dend_reorder = TRUE,
        col=viridis(100, direction = -1, option = "A"),
        column_names_gp = gpar(fontsize = 10),
        row_names_gp = gpar(fontsize = 12),
        show_heatmap_legend = T)

# Final_Net <- rbind(edged_ceQTLs[,c("edged", "Source")], 
#                    edged_teQTLs[,c("edged", "Source")],
#                    CoexDB_Sig[,c("edged", "Source")])


#CoexDB_45.list <- lapply(CoexDB_45, MakeEgde.list)
#CoexDB_45.list[["Net304.15k"]] <- CoexDB_15k$egde
#CoexDB_45.list[["Net304.PCC0.5"]] <- CoexDB_PCC0.5$egde

# long df vesion
CoexDB_45_long <- as_tibble(rbindlist(CoexDB_45, idcol = "Source"))
colnames(CoexDB_45_long) <- c("Source", "TF", "Target", "score")
CoexDB_45_long["edge"] <- paste(CoexDB_45_long$TF, CoexDB_45_long$Target, sep="_")
#
CoexDB["Source"] <- "n304"
colnames(CoexDB)[1:2] <- c("TF", "Target")
#CoexDB_15k["Source"] <- "n304"

#
CoexDB_long_DB <- rbind(CoexDB_45_long[,c("TF", "Target", "edge", "Source")], 
                        CoexDB[, c("TF", "Target", "edge", "Source")])

unique(CoexDB_long_DB$Source)
write.table(CoexDB_long_DB, "eQTL_Results/CoExpDB.txt", sep = '\t', row.names = F)

# Frequence of egde supported by co-expression
TF_Cis_Table <- as_tibble(as.data.frame(table(subset(CoexDB_long_DB, edge %in% edged_ceQTLs.TF_TF$edged)$edge)))
TF_Trans_Table <- as_tibble(as.data.frame(table(subset(CoexDB_long_DB, edge %in% edged_teQTLs_TF_TF$edged)$edge)))
#
Phe_cis_Table  <- as_tibble(as.data.frame(table(subset(CoexDB_long_DB, edge %in% edged_ceQTLs$edged)$edge)))
Phe_Trans_Table <- as_tibble(as.data.frame(table(subset(CoexDB_long_DB, edge %in% edged_teQTLs$edged)$edge)))


TF_cis_CoEp_plot <- ggplot(as.data.frame(table(TF_Cis_Table$Freq)), aes(x=Var1, y=Freq) ) +
  geom_bar(stat="identity", fill="#56B4E9") +
  geom_text(aes(label=Freq), vjust=+0.1, size=5) +
  theme_pubclean()

TF_trans_CoEp_plot <- ggplot(as.data.frame(table(TF_Trans_Table$Freq)), aes(x=Var1, y=Freq) ) +
  geom_bar(stat="identity", fill="#56B4E9") +
  geom_text(aes(label=Freq), vjust=+0.1, size=5) +
  theme_pubclean()

Phe_cis_CoEp_plot <- ggplot(as.data.frame(table(Phe_cis_Table$Freq)), aes(x=Var1, y=Freq) ) +
  geom_bar(stat="identity", fill="hotpink") +
  geom_text(aes(label=Freq), vjust=+0.1, size=5) +
  theme_pubclean()

Phe_trans_CoEp_plot <- ggplot(as.data.frame(table(Phe_Trans_Table$Freq)), aes(x=Var1, y=Freq) ) +
  geom_bar(stat="identity", fill="hotpink") +
  geom_text(aes(label=Freq), vjust=+0.1, size=5) +
  theme_pubclean()

TF_cis_CoEp_plot <- ggpar(TF_cis_CoEp_plot,  submain = "Cis-eQTL (TFs)", xlab="Coexp. Support", font.tickslab = 12)
TF_trans_CoEp_plot <- ggpar(TF_trans_CoEp_plot, submain = "Trans-eQTL (TFs)", xlab="Coexp. Support", font.tickslab = 12)
#
Phe_cis_CoEp_plot <- ggpar(Phe_cis_CoEp_plot, submain = "Cis-eQTL (Phe)", xlab="Coexp. Support", font.tickslab = 12)
Phe_trans_CoEp_plot <- ggpar(Phe_trans_CoEp_plot, submain = "Trans-eQTL (Phe)", xlab="Coexp. Support", font.tickslab = 12)

ggarrange(TF_cis_CoEp_plot, TF_trans_CoEp_plot, nrow = 1, align = 'h')
ggarrange(Phe_cis_CoEp_plot, Phe_trans_CoEp_plot, nrow = 1, align = 'h')

subset(TF_Cis_Table, Freq>=5) 
subset(TF_Trans_Table, Freq>=6) 
#
subset(Phe_cis_Table, Freq>=1) 
subset(Phe_Trans_Table, Freq>1) 
View(Phe_Trans_Table)

UGTs <- as.character(read.table("UGTs.txt", h=F)$V1)

#### Final list of eQTLs
Phe_eQTLs <- rbind(edged_ceQTLs, edged_teQTLs)
TF_eQTLs <- rbind(edged_ceQTLs.TF_TF, edged_teQTLs_TF_TF)

eQTLs <- rbind(edged_ceQTLs, edged_teQTLs,
               edged_ceQTLs.TF_TF, edged_teQTLs_TF_TF)

table(Phe_eQTLs$Source)
table(TF_eQTLs$Source)

tem <- subset(as_tibble(as.data.frame(table(subset(CoexDB_long_DB, TF=="Zm00001d047671")$Target))), Freq>0)

PheGenes

CandiatesFromCoExp <- function(tf){
  tem <- as_tibble(as.data.frame(table((subset(CoexDB_long_DB, TF==tf)[,1:2])$Target))) 
  #
  tem <- subset(tem, Freq>1)
  tem <- tem[order(-tem$Freq),]
  tem["InPhe"] <- tem$Var1 %in% PheGenes$GeneID
  print(table(unique(tem[,c(1,3)])$InPhe))
  tem <- subset(tem, InPhe==TRUE)
  tem <- left_join(tem, PheGenes, by=c("Var1"="GeneID"))
  colnames(tem)[1] <- "GeneID"
  return(tem[,-c(3)])
}


CandiatesFromCoExp("Zm00001d047671")

write.table(CoexDB_long_DB, "ExpressionData/CoexpDB.LongFormat.txt", row.names = F, sep = "\t", quote = F) 
tem <- CandiatesFromCoExp("Zm00001d027435")
View(tem)

subset(Phe_eQTLs, GeneID=="Zm00001d047671")
dim(as.data.frame(table(subset(CoexDB_long_DB, TF=="Zm00001d047671")$Target)))

CoexDB_long_DB[grep("Zm00001d047671", CoexDB_long_DB$edge),]
Phe_Trans_Table[grep("Zm00001d047671", TF_Trans_Table$Var1),]

Phe_cis_Table
head(Phe_Trans_Table,30)

library(stringr)
str_split_fixed(before$type, "_and_", 2)


## trans

##################################################
######     Expression data description      ######
##################################################

row.names(ExpData) <- ExpData$gene
ExpData[1:5,1:5]
ExpData <- ExpData[,-c(1:5)]



ExpData <- ExpData[(row.names(ExpData) %in%  PheGenes$GeneID), ]
ExpData[1:5,1:5]

MedianData <- data.frame(colMeans(ExpData))

colnames(MedianData) <- "MeanData"

MeanEx_byline <- ggplot(MedianData, aes(x=MeanData)) +
  geom_histogram(bins = 100) + 
  geom_vline(xintercept = mean(MedianData$MeanData)) +
  geom_vline(xintercept = (mean(MedianData$MeanData) + sd(MedianData$MeanData)), linetype="dotted") +
  geom_vline(xintercept = (mean(MedianData$MeanData) - sd(MedianData$MeanData)), linetype="dotted") +
  theme_pubr()

ggpar(MeanEx_byline, xlab = "Average expression - Phe-Genes")

lowerthreshold <- mean(MedianData$MeanData) - sd(MedianData$MeanData)
upperthreshold <- mean(MedianData$MeanData) + sd(MedianData$MeanData)

MedianData_low <- subset(MedianData, MeanData <= lowerthreshold) 
MedianData_up <- subset(MedianData, MeanData >= upperthreshold) 

MedianData_low["Class"] <- "Low"
MedianData_up["Class"] <- "Up"
MedianData_Constrasting <- rbind(MedianData_up, MedianData_low)

table(MedianData_Constrasting$Class)

PlotConst.lines <- ggplot(MedianData_Constrasting, aes(x=Class, y=MeanData, fill=Class)) +
  geom_boxplot() +
  theme_pubr()

ggpar(PlotConst.lines, ylab = "Average expression - Phe-Genes", xlab = "Lines")

write.table( MedianData_Constrasting, "MedianData_Constrasting.txt", sep = "\t", quote = F)
##################################################


#############################################################
#############    Data from  Public domain   #################
#############################################################

PengData<- read.table("../Data_45_net/")

Public_Net_Phe <- unique(subset(Public_Net, V2 %in% PheGenes$GeneID))

TablePhe <- as_tibble(as.data.frame(table(Public_Net_Phe$V1)))


TablePhe <- TablePhe[order(-TablePhe$Freq),]
colnames(TablePhe) <- c("TF", "Phe_Targets")

write.table(TablePhe, "PhenolicNet.txt", row.names = F, quote = F, sep = "\t")
#############################################################



