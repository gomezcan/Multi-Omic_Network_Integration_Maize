library(rrvgo)
library(topGO)
library(GOSemSim)
library(enrichplot)
library(scales)
library(tidyverse)
library(ggrepel)
library(ggpubr)
library(data.table)
library(GeneOverlap)
library(viridis)
library(ComplexHeatmap)
#library(fgsea)
library(reshape2)
library(circlize)
library(data.table)
library(ggVennDiagram)
library(scales)
library(purrr)
library(gplots)
library(ggplot2)
library(parallel)
library(org.Zmays.eg.db)
library(patchwork)

##########################################################
######                  Functions                   ######
##########################################################

ReplaceNamePWY <- function(ids){
  # 
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

vennfunc <- function(list_x){
  venn <- Venn(list_x)
  data <- process_data(venn)
  colorGroups <- c(CEN = 'goldenrod1', GRN='steelblue1', GAN='darkorchid1')
  colfunc <- colorRampPalette(colorGroups)
  col <- colfunc(3)
  
  colorGroups <- c(CEN="gray100",GRN="gray99", GAN="gray98")
  colfunc2 <- colorRampPalette(colorGroups)
  col2 <- colfunc2(7)
  
  ggplot() +
    geom_sf(aes(fill=name), data = venn_region(data), show.legend = F) +
    geom_sf(aes(color=name), size = 1.5, data = venn_setedge(data), show.legend = F) +
    #
    geom_sf_text(aes(label = name), size=6, data = venn_setlabel(data)) +
    geom_sf_text(aes(label= scales::comma(count, accuracy = 1)), size=5, data = venn_region(data)) +
    #
    scale_fill_manual(values = col2) + # 
    scale_color_manual(values = alpha(col, .5)) +
    #
    theme_void() +
    theme(plot.margin = unit(c(0.5,1,1,0.1), "cm")) +
    xlim(-150,1000)
}

vennfuncInt <- function(list_x){
  venn <- Venn(list_x)
  data <- process_data(venn)
  
  colorGroups <- c(CEN = '#FF9933', GRN='#1E90FF', GAN='#FFD700', eGRN='#FF1493')
  colfunc <- colorRampPalette(colorGroups)
  col <- colfunc(4)
  
  colorGroups <- c(CEN="gray100",GRN="gray99", GAN="gray98", eGRN='gray97')
  colfunc2 <- colorRampPalette(colorGroups)
  col2 <- colfunc2(15)
  
  ggplot() +
    geom_sf(aes(fill=name), data = venn_region(data), show.legend = F) +
    geom_sf(aes(color=name), size = 1.5, data = venn_setedge(data), show.legend = F) +
    #
    geom_sf_text(aes(label = name), size=6, data = venn_setlabel(data)) +
    geom_sf_text(aes(label= scales::comma(count, accuracy = 1)), size=5, data = venn_region(data)) +
    #
    scale_fill_manual(values = col2) + # 
    scale_color_manual(values = alpha(col, .5)) +
    #
    theme_void() +
    theme(plot.margin = unit(c(0.5,0.1, 0.5, 0.1), "cm"),
          text = element_text(family="Helvetica")) # +
    #xlim(-150,1000)
}

## PYW Enrichment 
Enrichmet_classes <- function(network){
  ## Count TF targets in network
  # Count Total
  network <- subset(network, Target %in% Syntenic)
  #
  Total_targtes <- as_tibble(as.data.frame(table(unique(network[,c("TF", "Target")])$TF)))
  colnames(Total_targtes) <- c('TF', 'targets') 
  Total_targtes$TF <- as.character(Total_targtes$TF)
  
  # list input: network
  network.list <- unique(network[,c("TF", "Target")])
  network.list <- split(network.list$Target, network.list$TF)
  
  ## Compare list of predicted targets vs annotated genes in query.vector 
  #
  print(". Pre-newGOM .")
  go.obj <- newGOM(network.list, CornCYC.list, genome.size=length(Syntenic)) # annotated genes in Genome v4
  
  Pval <- getMatrix(go.obj, name="pval")
  Common <- getMatrix(go.obj, name="intersection")
  
  print(". Post-newGOM .")
  
  ### Summary tables
  ## adjust p value
  Pval <- as.data.frame(Pval)
  #Pval[,1:ncol(Pval)] <- apply(Pval[,1:ncol(Pval)], 2, p.adjust)
  
  Pval_table   <- as_tibble(reshape2::melt(as.matrix(Pval)))
  colnames(Pval_table) <- c('TF', 'PWY', 'Pval')
  
  #
  Common_table <- as_tibble(reshape2::melt(as.matrix(Common)))
  colnames(Common_table) <- c('TF', 'PWY', 'n.targ')
  
  # Add predicted target in class by TF
  Pval_table <- left_join(Pval_table, Common_table , by=c('TF', 'PWY'))
  
  # Add total predicted targets
  Pval_table <- left_join(Pval_table, Total_targtes, by="TF")
  
  #write.table(Pval_table, name, sep = "\t", row.names = F, quote = F)
  print(paste("... Done ...", sep = ""))
  return(Pval_table)
}

## chop a string by a separator and return specified field
chop=function(myStr,mySep,myField){
  
  choppedString=sapply(strsplit(myStr,mySep),"[",myField)
  if(length(myField)>1){
    choppedString=apply(choppedString,2,function(x){paste0(x[!is.na(x)],collapse=mySep)})
  }
  return(choppedString)
}

GetGO <- function(degs, mutant){
  
  # Use a list of DEGs and the name of the mutant (string)
  # to identify GOs enriched. Required to have a background predefined
  #
  # mutant: description of a category of source of the data
  # degs: set of genes to test enrichment
  print(length(degs))
  GeneList <- factor(as.integer(background_IDs %in% degs))
  names(GeneList) <- background_IDs
  GOdata_BP <- new("topGOdata", ontology = "BP", allGenes = GeneList, annot = annFUN.gene2GO, gene2GO = background)
  #GOdata_MF <- new("topGOdata", ontology = "MF", allGenes = GeneList, annot = annFUN.gene2GO, gene2GO = background)
  #GOdata_CC <- new("topGOdata", ontology = "CC", allGenes = GeneList, annot = annFUN.gene2GO, gene2GO = background)
  
  #### Define test ####
  test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test", nodeSize = 10)
  
  ### test enrichment 
  results_BP <- getSigGroups(GOdata_BP, test.stat)
  
  ### save pdf Graph
  #namepdf=paste("GOs_Plots/GO.BP_",mutant, "", sep = "")
  #printGraph(GOdata_BP, results_BP, firstSigNodes=20,  fn.prefix = namepdf, useInfo = "def", pdfSW = TRUE) #
  
  ######## Get Significant GOs ########  
  Res_DF_BP <- as_tibble(GenTable(GOdata_BP, classic = results_BP, topNodes = 1000, orderBy='Fis')) # save as dataframe
  Res_DF_BP["Mutant"] <- mutant # add Mutant column name
  Res_DF_BP$classic <- as.numeric(Res_DF_BP$classic)
  
  ##### get all GOs and their genes from the topGO result #####
  gs <- genesInTerm(GOdata_BP) # list genes by GO
  # 
  ANOTATION = lapply(gs,function(x) x[x %in% degs]) ## Get only my Differential expressed genes
  
  ### Get only the GO's located in the result of topGO in Res_DF
  DF_GO_Genes <- ANOTATION[Res_DF_BP$GO.ID] # list
  
  ## Transform it to a data frame.
  DF_GO_Genes = list_to_DF(DF_GO_Genes)
  DF_GO_Genes <- left_join(DF_GO_Genes, Res_DF_BP[,c(1,2,6)], by='GO.ID')   # left join to add GO info
  DF_GO_Genes["Mutant"] <- mutant # add Mutant column name
  
  write.table(DF_GO_Genes, 
              paste("GOsResults/Genes_GOBP_",mutant, ".txt", sep = ""), sep = '\t', quote = F,
              row.names = F)
  
  #return(list(GOs_DF=Res_DF_BP, Genes=DF_GO_Genes)) # Return list of GOs-Stats and GeneID-GOs
  return(Res_DF_BP) # Return list of GOs-Stats and GeneID-GOs
}

SuperGO_CommomTarg <- function(tf){
  ## used TF's targets to test GO terms enrichment
  # 1. Select tf's Targets
  # 2. Filter out Targets predicted by only a layers
  # 3. Add up targets predicted by at least two layers
  # 4. test enrichment
  
  # Get network by TF
  network <- unique(Reduced_InteractionDB[Reduced_InteractionDB$TF==tf,]$Target)
  network <- network[network %in% Syntenic]
  #
  Total_targets <- length(network)

  
  if(Total_targets > 10){
    # If targets largert than 
    # Include 3_later targets into paired comparison if required 
    out <- GetGO(network, tf)
    
    # Save gene results GOs  
    write.table(out, paste0("GOsResults/GOsBP_",tf, ".txt"), sep = '\t', row.names = F, quote = F)
    
    }
  else { print(paste0(" Salado: ", tf, " ..")) }
  
  return(out)
}

ReduceGOs <- function(tf){
  
  # Defined 
  GO_vector = GOsDB_commonT[GOsDB_commonT$TF == tf, ]$GO.ID
  scores <- setNames(-log10(GOsDB_commonT[GOsDB_commonT$TF == tf, ]$FDR), GO_vector)
  
  if (length(GO_vector) > 1) {
    # Semantic similarity
    simMatrix <- calculateSimMatrix(GO_vector,  orgdb=org.Zmays.eg.db,  ont="BP", 
                                    semdata=Zm.GOSemSim.BP,
                                    method="Wang")
    
    # Reduce term
    reducedTerms <- reduceSimMatrix(simMatrix, scores, keytype="GENENAME",
                                    threshold=0.7, orgdb=org.Zmays.eg.db)
    
    # List of parents
    parent <- unique(reducedTerms[,c("parent", "parentTerm")])
    
    ReduceDB <- subset(GOsDB_commonT, TF == tf & GO.ID %in% parent$parent)
    ReduceDB <- left_join(ReduceDB, parent, by=c("GO.ID"="parent"))
    ReduceDB$Term <- ReduceDB$parentTerm
    ReduceDB <- ReduceDB[,-c(9)]
    
    return(ReduceDB)
    
  }
  
  return(subset(GOsDB_commonT, TF == tf))
}

##################################################
##########        Annotations       ##############
##################################################

saf <- as_tibble(read.table("Data/eQTL_data/Zea_mays.B73_RefGen_v4.46.saf", stringsAsFactors = F))
saf1 <- subset(saf, V5=="+")[,c(1,2,3)]
saf2 <- subset(saf, V5=="-")[,c(1,2,4)]
colnames(saf1) <- c("GeneID", "chrAnn", "TSS")
colnames(saf2) <- c("GeneID", "chrAnn", "TSS")
#
saf <- rbind(saf1, saf2)

# PDI
PDI <- unique(fread("../Fig_PDI/Only_PDI_NetworkFinal.10_14_2022.txt")[,c(2,3)])
colnames(PDI)[1] <- "Source"
dim(PDI)

mean(as.data.frame(table(PDI$Source))$Freq)
mean(as.data.frame(table(PDIeQTL$Source))$Freq)
  

# PDI eQTL
PDIeQTL <- unique(fread("../Fig_PDI/CisE_PDI_NetworkFinal.10_14_2022.txt")[,c(2,3)])
colnames(PDIeQTL)[1] <- "Source"

PDI <- rbind(PDI, PDIeQTL)

# CoExp
CoExp <- unique(fread("../Fig_Coexpression/CoExp_NetworkFinal.10_11_2021.txt"))
colnames(CoExp)[2] <- "Source"
CoExp <- unique(CoExp[,2:3])

# teQTL
teQTL <- unique(fread("../Fig_transeQTL/teQTL_NetworkFinal.10_11_2021.txt"))
colnames(teQTL)[1] <- "Source"

# genes in synteny
GenesMaize <- unique(as.character(read.table("Data/Annotations/Zm.v4.synteny.genes.txt", h=T)$gene_id))
length(GenesMaize)

# TF names
TFdic <- as_tibble(read.table("Data/Annotations/TF_Id_Name.txt", h=T, stringsAsFactors = F))

## Y1H network
Y1H <- as_tibble(read.table("Data/Annotations/Y1H.data.txt", h=T, stringsAsFactors = F))[,4:6]
ReplaceName(Y1H$TF.v4)


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

##################################################

##################################################
######        Common interactions       ##########
##################################################

# Define list
InteractionDB <- list()

InteractionDB[["GRN"]] <- paste0(PDI$Source, "_", PDI$Target)
InteractionDB[["eGRN"]] <- paste0(PDIeQTL$Source, "_", PDIeQTL$Target)
InteractionDB[["CEN"]] <- paste0(CoExp$Source, "_", CoExp$Target)
InteractionDB[["GAN"]] <- paste0(teQTL$Source, "_", teQTL$Target)

InteractionDB <- lapply(InteractionDB, unique)
length(InteractionDB)

# Interaction groups
InteractionDB_venn <- venn(InteractionDB)
InteractionDB_venn <- as.list(attr(InteractionDB_venn, "intersections"))
InteractionDB_venn <- plyr::ldply(InteractionDB_venn, data.table)
InteractionDB_venn <- as.data.table(InteractionDB_venn)

InteractionDB_venn <-  data.table(str_split_fixed(InteractionDB_venn$V1, pattern = "_", n = 2), Edge=InteractionDB_venn$.id)
colnames(InteractionDB_venn)[1:2] <- c("Source", "Target")

fwrite(InteractionDB_venn, "Full_Final_network.11022022.txt", sep = '\t', row.names = F, col.names = T, quote = F)

All_TFs <- unique(c(unique(PDI$Source), CoExp$Source, TF_CoR$GeneID))
All_TFs <- unique(c(unique(PDI$Source), CoExp$Source, TF_CoR$GeneID))

##################################################


#################################################################
#######            Common TFs among layers           ############
#################################################################


# Total TFs by layer
Total_TFs_list = list(CEN=CoExp$Source,  GRN=PDI$Source, GAN=teQTL$Source, eGRN=PDIeQTL$Source)
Total_TFs_list <- lapply(Total_TFs_list, unique)

# TF in venn groups
ven_file <- venn(Total_TFs_list)
ven_file <- as.list(attr(ven_file, "intersections"))

#
Total_TFs_DB <- as.data.table(plyr::ldply(ven_file, data.table)) 
colnames(Total_TFs_DB) <- c("Network", "TF")

Total_TFs_DB <- Total_TFs_DB[Total_TFs_DB$TF %in% All_TFs, ]
Total_TFs_DB[,"TFname"] <- ReplaceName(Total_TFs_DB$TF)
Total_TFs_DB

fwrite(subset(Total_TFs_DB, Network=='CEN:GRN:GAN:eGRN')[,2:3], "TF_In_3_Layers.txt", quote = F, sep = "\t", row.names = F)

vennfuncInt(Total_TFs_list)

GRN_size <- as.data.table(table(PDI$Source)) %>%
  dplyr::rename(TF=V1)

eGRN_size <- as.data.table(table(PDIeQTL$Source)) %>%
  dplyr::rename(TF=V1)

CEN_size <- as.data.table(table(CoExp$Source)) %>%
  dplyr::rename(TF=V1)

GAN_size <- as.data.table(table(teQTL$Source)) %>%
  dplyr::rename(TF=V1) %>%
  dplyr::filter(TF %in% All_TFs)

TFs_DB_list <- list(GRN = GRN_size$TF,
                    eGRN= eGRN_size$TF,
                    CEN = CEN_size$TF,
                    GAN = GAN_size$TF)

TFs_DB_list_10 <- list(GRN = GRN_size[GRN_size$N>=10, ]$TF,
                       eGRN= eGRN_size[eGRN_size$N>=10,]$TF,
                       CEN = CEN_size[CEN_size$N>=10, ]$TF,
                       GAN = GAN_size[GAN_size$N>=10, ]$TF)



#################################################################


####################################################################
##########       TFs by set of common interactions      ############
####################################################################

nrow(unique(InteractionDB_venn[,1:2]))
nrow(unique(subset(InteractionDB_venn, Source %in% All_TFs)[,1:2]))

# Interactions Freq
Interaction_Freq <- data.table(table(subset(InteractionDB_venn, Source %in% All_TFs)$Edge)) %>%
  arrange(-N)
Interaction_Freq

# Count TFs on Venn groups
TFs_InVennInte  <- data.table(table(unique(subset(InteractionDB_venn, Source %in% All_TFs)[,c("Source", "Edge")])$Edge))
Targ_InVennInte <- data.table(table(unique(subset(InteractionDB_venn, Source %in% All_TFs)[,c("Target", "Edge")])$Edge))

# Summary of interactions, TFs, and targets by networks
Interaction_Freq <- left_join(Interaction_Freq, TFs_InVennInte, by='V1')
Interaction_Freq <- left_join(Interaction_Freq, Targ_InVennInte, by='V1')
colnames(Interaction_Freq) <- c("Networks", 'Interaction', "TFs", "Targets")

sum(Interaction_Freq$Interaction)
sum(subset(Interaction_Freq, !(Networks %in% c("CEN", "GRN", "GAN")))$Interaction)

####################################################################

####################################################################
########     Enrichment of PWY in common interactions   ############
####################################################################

# Defined interactions with at least two lines of evidence 
Reduced_InteractionDB <- subset(InteractionDB_venn, !(Edge %in% c("CEN", "GRN", "GAN")))
colnames(Reduced_InteractionDB)[1] <- "TF"
Reduced_InteractionDB

# Total common targets
Total_Target <- as.data.table(table(unique(Reduced_InteractionDB[,1:2])$TF))
fwrite(Total_Target, "Total_Common_Targets.txt", sep = "\t", row.names = F)
colnames(Total_Target) <- c("TF", "nCommonTarg")

# PWY Enrichment test
CommonT_Enrichment <- Enrichmet_classes(Reduced_InteractionDB)
CommonT_Enrichment$PWY <- as.character(CommonT_Enrichment$PWY)
CommonT_Enrichment

CommonT_Enrichment <- left_join(CommonT_Enrichment, CornCYC_size, by="PWY")
CommonT_Enrichment <- left_join(CommonT_Enrichment, Total_Target, by="TF")

write.table(CommonT_Enrichment, 
            "PWY_GO_results/CommonTarg_PWY_enrichment.txt", 
            row.names = F, quote = F, sep = '\t')


# Add FDR by TF
CommonT_Enrichment %>% 
  #filter(n.targ >= 1) %>%
  group_by(TF) %>%
  mutate("FDR" = p.adjust(Pval, method = 'fdr')) %>%
  filter(Pval <= 0.05) -> Final_CommonT_Enrichment

Final_CommonT_Enrichment

Final_CommonT_Enrichment[,"TFname"] <- ReplaceName(Final_CommonT_Enrichment$TF)
Final_CommonT_Enrichment[,"PWYname"] <- ReplaceNamePWY(Final_CommonT_Enrichment$PWY)

# re-load names
TFdic <- as_tibble(read.table("Data/Annotations/TF_Id_Name.txt", h=T, stringsAsFactors = F))

write.table(left_join(Final_CommonT_Enrichment, CornCYC_size, by="PWY"), 
            "TopPWYs_CommonT.03232023.txt", row.names = F, quote = F, sep = '\t')

write.table(Final_CommonT_Enrichment,"TopPWYs_CommonT.03232023.txt", row.names = F, quote = F, sep = '\t')

length(unique(Final_CommonT_Enrichment$TF))
length(unique(Final_CommonT_Enrichment$PWY))

####################
# Vertical heatmap #
####################
MakeHeatmapVertical <- function(df_PWY, h) {
  
  # Add PWY name
  df_PWY <- left_join(df_PWY, CornCYC_size, by="PWY")
  df_PWY[,"Name"]  <- paste0("[", df_PWY$nPWY,"] ", ReplaceNamePWY(df_PWY$PWY))
  df_PWY$TF <- ReplaceName(df_PWY$TF)
  
  # Add % of PWY targeted
  df_PWY[,"Per"] <- round( (df_PWY$n.targ/df_PWY$nPWY)*100, 2)
  
  
  # heatmap values
  df1 <- df_PWY[,c("TF", "Name", "Per")] %>%
    group_by(TF, Name) %>%
    summarise(Per=max(Per))
  
  Values <- reshape2::dcast(df1, TF ~ Name, value.var="Per")
  row.names(Values) <- Values$TF
  Values <- Values[,-c(1)]
  Values[is.na(Values)] <- 0
  Values <- as.matrix(Values)
  
  
  # heatmap text
  df2 <- df_PWY[,c("TF", "Name", "n.targ")] %>%
    group_by(TF, Name) %>%
    summarise(n.targ=max(n.targ))
  
  Text <- reshape2::dcast(df2, TF ~ Name, value.var="n.targ")
  row.names(Text) <- Text$TF
  Text <- Text[,-c(1)]
  Text[is.na(Text)] <- 0
  Text <- as.matrix(t(Text))
  
  cor_scale <- colorRamp2(seq(0,99,1), c("white", viridis(99, direction = -1, option = "B")))
  colnames(Values) <- gsub("_", " ", colnames(Values))
  
  hm_common <- Heatmap(t(Values), 
                       column_names_rot = 75,
                       
                       cell_fun = function(j, i, x, y, width, height, fill) {
                         if((Text)[i, j]>0){
                           grid.text(sprintf("%.f", (Text)[i, j]), x, y, 
                                     gp = gpar(fontsize = 6, col = "gray80"))
                         }
                         
                       },
                       name = "PWY Percentage", 
                       cluster_columns = TRUE, column_dend_reorder = TRUE,
                       show_row_dend = FALSE, show_column_dend = FALSE,
                       col=cor_scale,
                       width = unit(12, "cm"),
                       height = unit(h, "cm"),
                       column_names_gp = gpar(fontsize = 7),
                       row_names_gp = gpar(fontsize = 7),
                       show_heatmap_legend = T,
                       heatmap_legend_param = list(title_gp = gpar(fontsize = 10), 
                                                   labels_gp = gpar(fontsize = 10),
                                                   width = unit(15, "mm"),
                                                   direction = "horizontal")) 
  
  return(hm_common)
  
}

HeatmapCommonPWYs <- MakeHeatmapVertical(Final_CommonT_Enrichment, 15)
draw(HeatmapCommonPWYs, heatmap_legend_side = "bottom")  

tiff("Plots/Plot_HeatmapCommonPWYs.tiff", units="in", width=8, height=8, res=300)
print(draw(HeatmapCommonPWYs, heatmap_legend_side = "bottom")  )
dev.off()

# scale_y_log10(
#   breaks = scales::trans_breaks("log10", function(x) 10^x),
#   labels = scales::trans_format("log10", scales::math_format(10^.x))
# ) +
#   annotation_logticks(color = 'white', sides = 'left')
####################

####################################################################

####################################################################
########     Enrichment of GOs in common interactions   ############
####################################################################

GenesList <- unique(Reduced_InteractionDB$TF)
Lgenes <- length(GenesList)


w=40  # Size of range to test
print(".. Ready to start ..")
GOsDB_commonT <- list()

for (i in seq(0, Lgenes, w)){
  max=Lgenes
  Start=i+1
  end=i+w
  
  if (end<max){
    listtotest <- GenesList[Start:end]
    GOsDB_commonT <- c(GOsDB_commonT, mclapply(listtotest, SuperGO_CommomTarg, mc.cores=w))
  }
  else{
    r <- paste(Start, max, sep = "-")
    listtotest <- GenesList[Start:max]
    GOsDB_commonT <- c(GOsDB_commonT, mclapply(listtotest, SuperGO_CommomTarg, mc.cores=w))
  }
}


# Filter out empty DFs
GOsDB_commonT <- GOsDB_commonT[unlist(lapply(GOsDB_commonT, function(x) is.data.frame(x)))]
#GOsDB_commonT <- list.files(path = "GOsResults/", pattern = '^GOsBP_*')
#GOsDB_commonT <- lapply(GOsDB_commonT, function(x) fread(paste0("GOsResults/",x)))
# Combine results
GOsDB_commonT <- rbindlist(GOsDB_commonT, idcol = F)
colnames(GOsDB_commonT)[7] <- 'TF'
GOsDB_commonT <- left_join(GOsDB_commonT, Total_Target, by="TF")

# Results used on Fig_MethodsComparison
write.table(GOsDB_commonT, "PWY_GO_results/CommonTarg_GO_enrichment.txt",  row.names = F, quote = F, sep = '\t')


# adjust Pvals and filter 
GOsDB_commonT <- GOsDB_commonT %>% 
  group_by(TF) %>%
  mutate("FDR" = p.adjust(classic, method = 'fdr')) %>%
  filter(FDR <= 0.1) 

# Pre-calculate semantic similarity
Zm.GOSemSim.BP <- GOSemSim::godata(org.Zmays.eg.db, keytype="GENENAME", ont = 'BP')


# map GOs to parent and reduce GO df
GOsDB_commonT_reduced <- lapply(unique(GOsDB_commonT$TF), ReduceGOs)

# Combine results
GOsDB_commonT_reduced <- rbindlist(GOsDB_commonT_reduced, idcol = F, fill = T)
GOsDB_commonT_reduced

GOsDB_commonT_reduced[,"GOname"] <- paste0('[', GOsDB_commonT_reduced$Annotated, '] ',GOsDB_commonT_reduced$Term)
GOsDB_commonT_reduced[,"TFname"] <- ReplaceName(GOsDB_commonT_reduced$TF)

write.table(GOsDB_commonT_reduced, "TopGOs_CommonT.03232023.txt", row.names = F, quote = F, sep = '\t')

####################################################################

Final_CommonT_Enrichment
GOsDB_commonT_reduced

####################################################################
########                       Plots                    ############
####################################################################

######
# 0: Sup. Fig. 4
######
## Totals common TFs by layer with and without filtered by n.targets

Plot_CommonTFs <- vennfuncInt(TFs_DB_list)/vennfuncInt(TFs_DB_list_10)
pdf(file = "Plots/Plot_CommonTFs_FigS4.pdf", width = 5, height = 8) 
print(Plot_CommonTFs)
dev.off()

######
# 1
######
InteractionDB_venn

## Plot of total interactions by comparison 
Plot_vennInte <- vennfuncInt(InteractionDB)


tiff("Plots/Plot_VennInteractions.tiff", units="in", width=4, height=4, res=300)
print(Plot_vennInte)
dev.off()

FigS4abc <- vennfuncInt(TFs_DB_list)/vennfuncInt(TFs_DB_list_10)/vennfuncInt(InteractionDB)
pdf(file = "Plots/Plot_FigS4abc.pdf", width = 5, height = 12) 
print(FigS4abc)
dev.off()

FigS4abc

######

######
# 2
######
## TFs annotation Freq
#View(Final_CommonT_Enrichment)

# list
TFs_anno_Freq <-list(PWY=unique(Final_CommonT_Enrichment$TF),
                     GO=unique(GOsDB_commonT_reduced$TF))

# Interaction groups
TFs_anno_Freq <- venn(TFs_anno_Freq)
TFs_anno_Freq <- as.list(attr(TFs_anno_Freq, "intersections"))
TFs_anno_Freq <- as.data.table(plyr::ldply(TFs_anno_Freq, data.table))

# Frequency
TFs_anno_FreqTable <-  as.data.table(table(TFs_anno_Freq$.id))
colnames(TFs_anno_FreqTable) <- c("Class", "TFs")

##
# 
##
Plot2 <- ggplot(TFs_anno_FreqTable, aes(x=Class, y=TFs, fill=Class)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_pubclean() +
  geom_text(aes(label=scales::comma(TFs), y=TFs+2, x=Class),
            position = position_dodge(0.9)) +
  scale_y_continuous(expand=c(0,0), limits = c(0,240)) +
  scale_x_discrete(expand=c(0,0)) +
  scale_fill_viridis(discrete = T, alpha = 0.5) +
  xlab("Annotation") +
  theme(strip.text.x = element_text(size = 12), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text=element_text(size=12), 
        legend.position = 'none',
        text = element_text(size=12, family="Helvetica"))

Plot2

######

####################################################################


table(GOsDB_commonT$Mutant)


GOsDB_commonT_reduced

GOsHetamp <- MakeHeatmapVerticalGOs(GOsDB_commonT_reduced, 40)

MakeHeatmapVerticalGOs <- function(df_PWY, h) {
  # Add % of PWY targeted
  df_PWY[,"Per"] <- round( (df_PWY$Significant/df_PWY$Annotated)*100, 2)
  
  
  # heatmap values
  df1 <- df_PWY[,c("TFname", "GOname", "Per")] %>%
    group_by(TFname, GOname) %>%
    summarise(Per=max(Per))
  
  Values <- reshape2::dcast(df1, TFname ~ GOname, value.var="Per")
  row.names(Values) <- Values$TFname
  Values <- Values[,-c(1)]
  Values[is.na(Values)] <- 0
  Values <- as.matrix(Values)
  
  
  # heatmap text
  df2 <- df_PWY[,c("TFname", "GOname", "Significant")] %>%
    group_by(TFname, GOname) %>%
    summarise(n.targ=max(Significant))
  
  Text <- reshape2::dcast(df2, TFname ~ GOname, value.var="n.targ")
  row.names(Text) <- Text$TFname
  Text <- Text[,-c(1)]
  Text[is.na(Text)] <- 0
  Text <- as.matrix(t(Text))
  
  cor_scale <- colorRamp2(seq(0,99,1), c("white", viridis(99, direction = -1, option = "B")))
  colnames(Values) <- gsub("_", " ", colnames(Values))
  
  hm_common <- Heatmap(t(Values), 
                       column_names_rot = 75,
                       
                       cell_fun = function(j, i, x, y, width, height, fill) {
                         if((Text)[i, j]>0){
                           grid.text(sprintf("%.f", (Text)[i, j]), x, y, 
                                     gp = gpar(fontsize = 6, col = "gray80"))
                         }
                         
                       },
                       name = "PWY Percentage", 
                       cluster_columns = TRUE, column_dend_reorder = TRUE,
                       show_row_dend = FALSE, show_column_dend = FALSE,
                       col=cor_scale,
                       width = unit(20, "cm"),
                       height = unit(h, "cm"),
                       column_names_gp = gpar(fontsize = 7),
                       row_names_gp = gpar(fontsize = 7),
                       show_heatmap_legend = T,
                       heatmap_legend_param = list(title_gp = gpar(fontsize = 10), 
                                                   labels_gp = gpar(fontsize = 10),
                                                   width = unit(15, "mm"),
                                                   direction = "horizontal")) 
  
  return(hm_common)
  
}

#tiff("Plots/Plot_C5_GOBP.tiff", units="in", width=10, height=4, res=300)
#print(treemapPlot(reducedTerms))
#dev.off()

GOsDB_size <- unique(GOsDB_commonT_reduced[,1:3])

##################################################
# Zm00001d006236	MYB31
# Zm00001d013443	COL8
# Zm00001d047017	bHLH91

tf = "Zm00001d006236"

SummaryHeatMap_PWYGO <- function(tf){
  # Select pwys enriched in target TF
  target_peys <- subset(Final_CommonT_Enrichment, TF==tf)$PWY
  
  # Select GOs enriched in target TF
  target_gos <- unique(subset(GOsDB_commonT_reduced, Mutant==tf)$GO.ID)
  
  # Genes in PWY
  PWYs_genes <- CornCYC[CornCYC$Pathway.id %in% target_peys,]
  PWYs_genes <- unique(PWYs_genes[PWYs_genes$GeneID %in% Syntenic, 1:3])
  
  #
  GOS_genes <- fread(paste0("GOsResults/Genes_GOBP_", tf, ".txt"), sep = '\t')
  GOS_genes <- subset(GOS_genes, GO.ID %in% target_gos)
  
  # TF Targets
  TFtarg <- as.data.table(subset(Reduced_InteractionDB, TF==tf & Target %in% unique(PWYs_genes$GeneID)))
  
  # TF Targets: GOs
  TFtargGOs <- as.data.table(subset(Reduced_InteractionDB, TF==tf & Target %in% unique(GOS_genes$GeneID)))
  
  # map target class to pwy genes
  PWYs_genes <- left_join(PWYs_genes, TFtarg[,2:3], by=c('GeneID'='Target'))
  #PWYs_genes$Edge[is.na(PWYs_genes$Edge)] <- 'Not target'
  
  # map target class to GOs genes
  GOS_genes <- left_join(GOS_genes, TFtargGOs[,2:3], by=c('GeneID'='Target'))
  
  # Count Freq
  PWYs_genes <- as.data.frame(table(PWYs_genes[,c(1,4)]))
  
  # Count Freq in GOs
  GOS_genes <- as.data.table(as.data.frame(table(GOS_genes[,c(1,6)])))
  GOS_genes <- left_join(GOS_genes, GOsDB_size, by="GO.ID")
  GOS_genes <- subset(GOS_genes, Annotated <= 1000)
  GOS_genes[,"Name"]  <- paste0(GOS_genes$Term, " [", GOS_genes$Annotated,"]")
  
  # Add PWY name
  PWYs_genes <- left_join(PWYs_genes, CornCYC_size, by=c('Pathway.id'="PWY"))
  PWYs_genes[,"Name"]  <- paste0(ReplaceNamePWY(PWYs_genes$Pathway.id), " [", PWYs_genes$nPWY,"]")
  PWYs_genes$Name <- gsub("_", " ", PWYs_genes$Name)
  
  # Add % of PWY targeted
  PWYs_genes[,"Per"] <- round((PWYs_genes$Freq/PWYs_genes$nPWY)*100, 2)
  GOS_genes[,"Per"] <- round((GOS_genes$Freq/GOS_genes$Annotated)*100, 2)
  
  PWYs_genes$Per[PWYs_genes$Per > 40] <- 40
  GOS_genes$Per[GOS_genes$Per >   40] <- 40
    
  PWYs_genes %>%
    dplyr::group_by(Edge) %>%
    dplyr::filter(sum(Per) > 0) -> PWYs_genes
  
  GOS_genes %>%
    dplyr::group_by(Edge) %>%
    dplyr::filter(sum(Per) > 0) -> GOS_genes
  
  Plot <- ggplot(PWYs_genes, aes(x=Name, y=Edge, fill=Per, label=Freq)) + 
    geom_tile(show.legend = FALSE) + #
    geom_text(aes(color=Freq), show.legend = FALSE, size=3) + # 
    scale_fill_viridis(discrete=FALSE, limits = c(0, 40)) + # 
    scale_color_viridis(discrete=FALSE, direction = -1, name=NULL, option = 'B') +
    scale_x_discrete(expand = c(0,0)) + 
    scale_y_discrete(expand = c(0,0)) + 
    xlab("") + ylab("") + 
    theme_pubclean() +
    theme(strip.text.x = element_text(size = 12),  # axis.text.x=element_blank(),
          axis.text = element_text(size=10),
          axis.text.x = element_text(angle = 70, hjust = 1, vjust = 1),
          text = element_text(size=10),
          axis.ticks.x=element_blank(),
          legend.title.align = 0,
          legend.position = 'bottom', 
          legend.direction = "horizontal",
          legend.title = element_text(hjust = 0.5, vjust = .85)) + 
    labs(fill="PWY/GO Percentage", subtitle = paste0("PWY: ",ReplaceName(tf))) 
  
  
  Plot2 <- ggplot(GOS_genes, aes(x=Name, y=Edge, fill=Per, label=Freq)) + 
    geom_tile() +
    geom_text(aes(color=Freq), show.legend = FALSE,size=3) +
    scale_fill_viridis(discrete=FALSE, limits = c(0, 40)) + # 
    scale_color_viridis(discrete=FALSE, direction = -1, name=NULL, option = 'B') +
    scale_x_discrete(expand = c(0,0)) + 
    scale_y_discrete(expand = c(0,0)) + 
    xlab("") + ylab("") + 
    theme_pubclean() +
    theme(strip.text.x = element_text(size = 12),  # axis.text.x=element_blank(),
          axis.text = element_text(size=10),
          axis.text.x = element_text(angle = 70, hjust = 1, vjust =  1),
          text = element_text(size=10),
          axis.ticks.x=element_blank(),
          legend.text = element_text(size=10),
          legend.title.align = 0,
          legend.position = 'bottom', 
          legend.direction = "horizontal",
          legend.title = element_text(hjust = 0.5, vjust = .85)) + 
    labs(fill="PWY/GO Percentage", subtitle = paste0("GO: ", ReplaceName(tf))) 
  
    #plotf <- Plot / Plot2 +  plot_layout(heights  = c(1, 2)) 
    # guides = "collect" & theme(legend.position = 'bottom')
  
  return(list(pPWY=Plot, pGO=Plot2))
  
  
}

Zm00001d006236_MYB31$pPWY
# 10x9
Zm00001d006236_MYB31  <- SummaryHeatMap_PWYGO("Zm00001d006236")
Zm00001d013443_COL8   <- SummaryHeatMap_PWYGO("Zm00001d013443")  
Zm00001d047017_bHLH91 <- SummaryHeatMap_PWYGO("Zm00001d047017")

######
tiff("Plots/Plot_PWYS_GOs_Zm00001d006236_MYB31_pwy.tiff", units="in", width=6, height=4, res=300)
print(Zm00001d006236_MYB31$pPWY)
dev.off()

tiff("Plots/Plot_PWYS_GOs_Zm00001d006236_MYB31_gos.tiff", units="in", width=6, height=6, res=300)
print(Zm00001d006236_MYB31$pGO)
dev.off()

tiff("Plots/Plot_PWYS_GOs_Zm00001d013443_COL8_pwy.tiff", units="in", width=6, height=6, res=300)
print(Zm00001d013443_COL8$pPWY)
dev.off()

tiff("Plots/Plot_PWYS_GOs_Zm00001d013443_COL8_gos.tiff", units="in", width=12, height=6, res=300)
print(Zm00001d013443_COL8$pGO)
dev.off()

tiff("Plots/Plot_PWYS_GOs_Zm00001d047017_bHLH91_pwy.tiff", units="in", width=4, height=5.5, res=300)
print(Zm00001d047017_bHLH91$pPWY)
dev.off()

tiff("Plots/Plot_PWYS_GOs_Zm00001d047017_bHLH91_gos.tiff", units="in", width=11.5, height=6.5, res=300)
print(Zm00001d047017_bHLH91$pGO)
dev.off()
######

bHLH91_plot <- Zm00001d047017_bHLH91$pPWY + Zm00001d047017_bHLH91$pGO + 
  plot_layout(widths = c(.20, 1), nrow = 1, guides = "collect") & theme(legend.position = 'bottom')

MYB31_plot <- Zm00001d006236_MYB31$pPWY + Zm00001d006236_MYB31$pGO + 
  plot_layout(widths = c(.1, 1), nrow = 1, guides = "collect") & theme(legend.position = 'bottom')

COL8_plot <- Zm00001d013443_COL8$pPWY + Zm00001d013443_COL8$pGO + 
  plot_layout(widths = c(.18, 1), nrow = 1, guides = "collect") & theme(legend.position = 'bottom')

tiff("Plots/Plot_PWYS_GOs_Zm00001d047017_bHLH91.tiff", units="in", width=11.5, height=6.5, res=300)
print(bHLH91_plot)
dev.off()

tiff("Plots/Plot_PWYS_GOs_Zm00001d006236_MYB31.tiff", units="in", width=6, height=5.4, res=300)
print(MYB31_plot)
dev.off()

tiff("Plots/Plot_PWYS_GOs_Zm00001d013443_COL8.tiff", units="in", width=12, height=6.4, res=300)
print(COL8_plot)
dev.off()
##################################################

