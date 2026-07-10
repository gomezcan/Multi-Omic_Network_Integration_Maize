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
library(fgsea)
library(reshape2)
library(circlize)
library(data.table)
library(ggVennDiagram)
library(scales)
library(purrr)
library(gplots)
library(ggplot2)

##########################################################
######                  Functions                   ######
##########################################################

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


flattenMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    Pval  =(cormat)[ut],
    Common = pmat[ut]
  )
}

overlap_significance <- function(genes_all, gene_sets, iterations) {
  ####
  # To test the significance of multiple gene sets (>2)
  observed <- length(reduce(gene_sets, intersect))
  simulated <- map_dbl(seq_len(iterations), function(x) {
    sim <- map(lengths(gene_sets), ~sample(genes_all, .x))
    sim <- length(reduce(sim, intersect))
    return(sim)
  })
  pval <- (sum(simulated >= observed) + 1) / (iterations + 1)
  # list(pval=pval, simulated_values=simulated, observed=observed)
  return(c(observed, pval))
}

overlap_significance_tf3 <- function(tf){
  # Comparison bases on simulation of the null distribution
  
  # make list of targtes
  List_test <- list(CEN=CEN_list[[tf]], GRN=GRN_list[[tf]], GAN=GAN_list[[tf]]) 
  List_test <- List_test[lengths(List_test) != 0]
  #print(lapply(List_test, length))
  
  # get p-value
  results <- overlap_significance(GenesMaize, List_test, 10000) # commonTarget, Pval
  
  #out <- c(tf, results)
  return(results)
}

overlap_significance_tf2 <- function(tf){
  # Comparison bases on hypergeometric null distribution
  
  # make list of targets
  List_test <- list(CEN=CEN_list[[tf]], GRN=GRN_list[[tf]], GAN=GAN_list[[tf]]) 
  List_test <- List_test[lengths(List_test) != 0]
  #print(lapply(List_test, length))
  
  # get p-value bases on sampling
  go.obj <- newGOM(List_test, genome.size=45546) # all vs all 
  Pval <- getMatrix(go.obj, name="pval")
  common <- getMatrix(go.obj, name="intersection")
  results <- as.vector(c(common, Pval))
  
  return(results)
}

overlap_significance_tf2_cen_grn <- function(tf){
  # Comparison bases on hypergeometric null distribution
  
  # make list of targets
  List_test <- list(CEN=CEN_list[[tf]], GRN=GRN_list[[tf]]) 
  List_test <- List_test[lengths(List_test) != 0]
  #print(lapply(List_test, length))
  
  # get p-value bases on sampling
  go.obj <- newGOM(List_test, genome.size=45546) # all vs all 
  Pval <- getMatrix(go.obj, name="pval")
  common <- getMatrix(go.obj, name="intersection")
  results <- as.vector(c(common, Pval))
  
  return(results)
}

overlap_significance_tf2_cen_gan <- function(tf){
  # Comparison bases on hypergeometric null distribution
  
  # make list of targets
  List_test <- list(CEN=CEN_list[[tf]], GAN=GAN_list[[tf]]) 
  List_test <- List_test[lengths(List_test) != 0]
  #print(lapply(List_test, length))
  
  # get p-value bases on sampling
  go.obj <- newGOM(List_test, genome.size=45546) # all vs all 
  Pval <- getMatrix(go.obj, name="pval")
  common <- getMatrix(go.obj, name="intersection")
  results <- as.vector(c(common, Pval))
  
  return(results)
}

overlap_significance_tf2_gan_grn <- function(tf){
  # Comparison bases on hypergeometric null distribution
  
  # make list of targets
  List_test <- list(CEN=GAN_list[[tf]], GRN=GRN_list[[tf]]) 
  List_test <- List_test[lengths(List_test) != 0]
  #print(lapply(List_test, length))
  
  # get p-value bases on sampling
  go.obj <- newGOM(List_test, genome.size=45546) # all vs all 
  Pval <- getMatrix(go.obj, name="pval")
  common <- getMatrix(go.obj, name="intersection")
  results <- as.vector(c(common, Pval))
  
  return(results)
}


add_PCC.GAN <- function(tf){
  
  # GAN net
  net <- unique(subset(teQTL, Source==tf))
  colnames(net)[2] <- "GeneID"
  
  ## Cor wPCC file
  pcc_files <-  list.files(path = "wPCC_net_only_TFs/", pattern = paste0("*",tf,'.txt'))
  
  # PCC net names
  PCC_Nets <- sapply(strsplit(pcc_files, split='.', fixed=TRUE), `[`, 2) # nets
  
  # Read co-expression files
  c=1
  for (f in pcc_files){
    PCCf <- as_tibble(fread(paste("wPCC_net_only_TFs/",f, sep = "")))[,2:3]
    colnames(PCCf)[2] <- PCC_Nets[c]
    net <- left_join(net, PCCf, by="GeneID")
    c=c+1
  }
  return(net)
}

list_to_DF <- function(list){
  # This function get a list of GOs and GeneIds to produce a table
  Net <- as_tibble(as.data.frame(matrix(0, nrow = 0, ncol = 2)))
  colnames(Net) <- c("GO.ID", "GeneID")
  
  GOs <- names(list)
  
  for (n in GOs){
    tem <- tibble(GO.ID=n, GeneID=list[[n]])
    Net <- rbind(Net, tem)
  }
  Net <- unique(Net)
  return(Net)
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
  
  return(list(GOs_DF=Res_DF_BP, Genes=DF_GO_Genes)) # Return list of GOs-Stats and GeneID-GOs
}

SuperGO_CommomTarg <- function(tf){
  ## used TF's targets to test GO terms enrichment
  # 1. Select tf's Targets
  # 2. Filter out Targets predicted by only a layers
  # 3. Add up targets predicted by at least two layers
  # 4. test enrichment
  
  # Get network by TF
  network <- TargetsClasses[[tf]]
  
  network <- subset(network, GeneID %in% Syntenic)
  
  #
  Total_targets <- as_tibble(as.data.frame(table(unique(network)$Class), stringsAsFactors = F))
  colnames(Total_targets) <- c('Class', 'nClass') 
  
  #
  # Remove targs without multiple layers 
  network <- subset(network, !(Class %in% c("CEN", "GRN", "GAN")))
  
  if(nrow(network)>1){
    # If tfs does have common targets
    # Include 3_later targets into paired comparison if required 
    
    if (length(subset(network, Class =="CEN:GRN:GAN")$GeneID)> 0) {
      
      CEN_GRN_GAN <- subset(network, Class =="CEN:GRN:GAN")$GeneID
      
      # Reduce list to targets with only two layers 
      network <- subset(network, Class !="CEN:GRN:GAN")
      
      # List input: network
      network.list <- split(network$GeneID, network$Class)
      network.list <- lapply(network.list, function (x) c(x, CEN_GRN_GAN))
      
    }
    
    # list input: network
    network.list <- split(network$GeneID, network$Class)
    
    ## Compare list of predicted targets vs annotated genes in query.vector 
    #
    print(". Pre-GO.")
    for (c in names(network.list)){
      names(network.list)
      out <- GetGO(network.list[[c]], c)
      #
      
      colnames(out$GOs_DF)[7] <- "Class"
      out$GOs_DF[,"TF"] <- tf
      #
      colnames(out$Genes)[5] <- "Class"
      out$Genes[,"TF"] <- tf
      
      # save GOs
      out$GOs_DF <- subset(out$GOs_DF, classic <= 0.05)
      write.table(out$GOs_DF, paste0("CommonTargGOs/CommonTarget_GOs_",tf, ".txt"), sep = '\t', row.names = F, quote = F)
      
      # save gene annotation GOs  
      write.table(out$Genes, paste0("CommonTargGOs/CommonTarg_Annotation_",tf, ".txt"), sep = '\t', row.names = F, quote = F)
      
    }
    
    
  }
  
  else { print(paste0(" Salado: ", tf, " ..")) }
  
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
rm(saf1)
rm(saf2)

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
ReplaceName(Y1H$TF.v4)

## syntenic genes 
Syntenic <- as_tibble(read.table("Data/Annotations/Zm.v4.synteny.genes.txt", h=T, stringsAsFactors = F, sep = '\t'))
Syntenic <- subset(Syntenic, gene_synteny=='syntenic')$gene_id

# CornCYC
CornCYC <- as_tibble(read.table("Data/Annotations/corn_pathways.0210325.reduced.txt", h=T, stringsAsFactors = F))
CornCYC <- CornCYC[!(CornCYC$Gene.id %in% "unknown"),]
CornCYC <- subset(CornCYC, Gene.id %in% Syntenic)
CornCYC[,"Class"] <- "Enzyme"
colnames(CornCYC)[3] <- "GeneID"
CornCYC$Pathway.id <- gsub("-", "_", CornCYC$Pathway.id)
#
CornCYC$Pathway.name <- gsub("</i>", "", gsub("<i>", "", CornCYC$Pathway.name))
CornCYC.list <- split(CornCYC$GeneID, CornCYC$Pathway.id)
CornCYC_size <- as_tibble(as.data.frame(table(unique(CornCYC[,c(1,3)])$Pathway.id), stringsAsFactors = F))
colnames(CornCYC_size) <- c("PWY", "nPWY")

# GOs term annotations
background <- readMappings("../Fig_GOs/synteny.ID_TopGO_V4_GRAMER.txt")
background_IDs <- as.character(unique(names(background)))
background_list <- unique(as.character(unlist(background)))


###################################################################################################


###################################################################################################
############                          Networks input                            ###################
###################################################################################################

# PDI
PDI <- unique(fread("../Fig_PDI/Only_PDI_NetworkFinal.10_14_2022.txt")[,c(2,3)])
colnames(PDI)[1] <- "Source"

ePDI <- unique(fread("../Fig_PDI/Only_PDI_NetworkFinal.10_14_2022.txt")[,c(2,3)])
colnames(PDI)[1] <- "Source"


# CoExp
CoExp <- as_tibble(unique(fread("../Fig_Coexpression/CoExp_NetworkFinal.10_11_2021.txt")))
colnames(CoExp)[2] <- "Source"
CoExp <- unique(CoExp[,2:3])

# teQTL
teQTL <- as_tibble(unique(fread("../Fig_transeQTL/teQTL_NetworkFinal.10_11_2021.txt")))
colnames(teQTL)[1] <- "Source"

# teQTL associated with TFs
teQTLtf <- subset(teQTL, Source %in% unique(c(TF_CoR$GeneID, PDI$Source, CoExp$Source))) 

All_TFs <- unique(c(PDI$Source, CoExp$Source, TF_CoR$GeneID))

###################################################################################################


#############################################################################################
##############            Common targets between layers          ############################
#############################################################################################

GenesMaize <- unique(as.character(read.table("Data/Annotations/Zm.v4.synteny.genes.txt", h=T)$gene_id))
length(GenesMaize)

# Total TFs by layer
Total_TFs_list = list(CEN=CoExp$Source,  GRN=PDI$Source, GAN=teQTLtf$Source)
Total_TFs_list <- lapply(Total_TFs_list, unique)

# TF in venn groups
ven_file <- venn(Total_TFs_list)
ven_file <- as.list(attr(ven_file, "intersections"))

TotalsWithTargets <- unique(as.character(unlist(Total_TFs_list)))

Count_CommonTargs <- function(tf){
  # Comparison bases on simulation of the null distribution
  
  # Make list of targtes
  List_test <- list(CEN=CEN_list[[tf]], GRN=GRN_list[[tf]], GAN=GAN_list[[tf]]) 
  List_test <- List_test[lengths(List_test) != 0]

  tem_ven_file <- venn(List_test, show.plot = F)
  tem_ven_file <- as.list(attr(tem_ven_file, "intersections"))
  
  tem_ven_file <- list_to_DF(tem_ven_file)
  tem_ven_file["TF"] <- tf
  #out <- c(tf, results)
  return(tem_ven_file)
}

# Label targets based on layer used to make the prediction by TF
TargetsClasses <-  lapply(TotalsWithTargets, Count_CommonTargs)
names(TargetsClasses) <- TotalsWithTargets

# Count common targets by layer and TF
TargetsClasses_Freq <- lapply(TargetsClasses, function(x) table(x$Class))
TargetsClasses_Freq <- lapply(TargetsClasses_Freq, as.data.frame)
names(TargetsClasses_Freq) <- TotalsWithTargets


TargetsClasses_Freq <- as_tibble(rbindlist(TargetsClasses_Freq, idcol = T))
colnames(TargetsClasses_Freq) <- c("TF", "Class", "Targets")
TargetsClasses_Freq["TFname"] <- ReplaceName(TargetsClasses_Freq$TF)

TargetsClasses_Freq$Class <- factor(TargetsClasses_Freq$Class, 
                                    levels = rev(c("GRN", "CEN", "GAN",
                                               "CEN:GRN", "CEN:GAN", "GRN:GAN",
                                               "CEN:GRN:GAN")))


## Annotate Common targets amount all three layers
ListCommonTargets <- function(tf){
  
  tem <- TargetsClasses[[tf]]
  tem <- subset(tem, Class=="CEN:GRN:GAN")
  return(tem)
}

TopCommonTargets <- lapply(subset(TargetsClasses_Freq, Class=="CEN:GRN:GAN")$TF, ListCommonTargets)
TopCommonTargets <- as_tibble(rbindlist(TopCommonTargets))
TopCommonTargets[,"TFname"] <- ReplaceName(TopCommonTargets$TF)
TopCommonTargets <- left_join(TopCommonTargets, CornCYC[,1:3], by="GeneID")
View(TopCommonTargets)

###################################################################################################
#########                                 Summary and plots                               #########
###################################################################################################

########################################################
########       Plots targets groups             ########
########################################################

Plot_TagetsGroups_1 <- ggplot(subset(TargetsClasses_Freq, TF %in% c(ven_file$CEN, ven_file$GAN, ven_file$GRN)), 
                              aes(y=Class, x=Targets)) + 
  geom_jitter(size=1, alpha=0.5) + theme_pubclean() +
  scale_x_continuous(trans = log2_trans(), 
                     breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x)))


Plot_TagetsGroups_2 <- ggplot(subset(TargetsClasses_Freq, TF %in% c(ven_file$`CEN:GRN`)), 
                              aes(y=Class, x=Targets)) + 
  geom_jitter(size=1, alpha=0.5) + theme_pubclean() +
  scale_x_continuous(trans = log2_trans(), 
                     breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x)))


Plot_TagetsGroups_3 <- ggplot(subset(TargetsClasses_Freq, TF %in% c(ven_file$`GRN:GAN`)), 
                              aes(y=Class, x=Targets)) + 
  geom_jitter(size=1, alpha=0.5) + theme_pubclean() +
  scale_x_continuous(trans = log2_trans(), 
                     breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x)))

Plot_TagetsGroups_4 <- ggplot(subset(TargetsClasses_Freq, TF %in% c(ven_file$`CEN:GAN`)), 
                              aes(y=Class, x=Targets)) + 
  geom_jitter(size=1, alpha=0.5) + theme_pubclean() +
  scale_x_continuous(trans = log2_trans(), 
                     breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x)))

Plot_TagetsGroups_5 <- ggplot(subset(TargetsClasses_Freq, TF %in% c(ven_file$`CEN:GRN:GAN`)), 
                              aes(y=Class, x=Targets)) + 
  geom_jitter(size=1, alpha=0.5) + theme_pubclean() +
  scale_x_continuous(trans = log2_trans(), 
                     breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x)))


Plot_TagetsGroups_1 <- ggpar(Plot_TagetsGroups_1,  font.tickslab = 10,  font.x = 10, font.y = 10,  legend = "bottom", ylab = "")
Plot_TagetsGroups_2 <- ggpar(Plot_TagetsGroups_2,  font.tickslab = 10,  font.x = 10, font.y = 10,  legend = "bottom", ylab = "")
Plot_TagetsGroups_3 <- ggpar(Plot_TagetsGroups_3,  font.tickslab = 10,  font.x = 10, font.y = 10,  legend = "bottom", ylab = "")
Plot_TagetsGroups_4 <- ggpar(Plot_TagetsGroups_4,  font.tickslab = 10,  font.x = 10, font.y = 10,  legend = "bottom", ylab = "")
Plot_TagetsGroups_5 <- ggpar(Plot_TagetsGroups_5,  font.tickslab = 10,  font.x = 10, font.y = 10,  legend = "bottom", ylab = "")


ggarrange(Plot_TagetsGroups_1, 
          Plot_TagetsGroups_2,
          Plot_TagetsGroups_3,
          Plot_TagetsGroups_4,
          Plot_TagetsGroups_5, ncol = 5, align = 'h')

# size = 3x8.5
ggarrange((vennfunc(Total_TFs_list) + labs(subtitle = "   Total TFs")), 
          Plot_TagetsGroups, 
          widths = c(1, 1.7))

########################################################


########################################################
########       Enrichment of PWYs               ########
########################################################

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
    xlim(-80,950)
}

TFdic <- as_tibble(read.table("Data/Annotations/TF_Id_Name.txt", h=T, stringsAsFactors = F))

TFs_with_HotTargets <- subset(TargetsClasses_Freq,  TF %in% ven_file$`CEN:GRN:GAN`)
TFs_with_HotTargets <- subset(TFs_with_HotTargets, !(Class %in% c("CEN", "GRN", "GAN")))


PWYEnr_CommonTarg <- function(tf){
  ## used TF's targets to test enrichmed of PWYs
  # 1. Select tf's Targets
  # 2. Filter out Targets predicted by only a layers
  # 3. Add up targets predicted by at least two layers
  # 4. test enrichment
  
  # Get network by TF
  network <- TargetsClasses[[tf]]
  
  network <- subset(network, GeneID %in% Syntenic)
  
  #
  Total_targets <- as_tibble(as.data.frame(table(unique(network)$Class), stringsAsFactors = F))
  colnames(Total_targets) <- c('Class', 'nClass') 
  
  #
  # Remove targs without multiple layers 
  network <- subset(network, !(Class %in% c("CEN", "GRN", "GAN")))
  
  if(nrow(network)>1){
    # If tfs does have common targets
    # Include 3_later targets into paired comparison if required 
    if (length(subset(network, Class =="CEN:GRN:GAN")$GeneID)> 0) {
      CEN_GRN_GAN <- subset(network, Class =="CEN:GRN:GAN")$GeneID
      
      # reduce list to targets with only two layers 
      network <- subset(network, Class !="CEN:GRN:GAN")
      
      # list input: network
      network.list <- split(network$GeneID, network$Class)
      network.list <- lapply(network.list, function (x) c(x, CEN_GRN_GAN))
      
    }
    
    
    # list input: network
    network.list <- split(network$GeneID, network$Class)
    
    # N of classes to test
    n_classes2test <- length(names(network.list))
    
    
    if (n_classes2test > 1) {
      ## Compare list of predicted targets vs annotated genes in query.vector 
      #
      print(". Pre-newGOM .")
      go.obj <- newGOM(network.list, CornCYC.list, genome.size=length(Syntenic)) # Annotated genes in Genome v4
      
      Pval <- getMatrix(go.obj, name="pval")
      Common <- getMatrix(go.obj, name="intersection")
      
      print(". Post-newGOM .")
      ### Summary tables
      
      ## adjust p value
      Pval <- as.data.frame(Pval)
      
      #
      Pval_table   <- as_tibble(reshape2::melt(as.matrix(Pval)))
      colnames(Pval_table) <- c('Class', 'PWY', 'Pval')
      
      #
      Common_table <- as_tibble(reshape2::melt(as.matrix(Common)))
      colnames(Common_table) <- c('Class', 'PWY', 'n.targ')    
      
      
    }
    else {
      ## Compare list of predicted targets vs annotated genes in query.vector 
      #
      print(". Pre-newGOM .")
      go.obj <- newGOM(network.list, CornCYC.list, genome.size=length(Syntenic)) # Annotated genes in Genome v4
      
      Pval <- getMatrix(go.obj, name="pval")
      Common <- getMatrix(go.obj, name="intersection")
      
      
      Pval_table <- tibble(Class=sapply(strsplit(names(Pval) , split='.', fixed=TRUE), `[`, 2),
                           PWY=sapply(strsplit(names(Pval) , split='.', fixed=TRUE), `[`, 1),
                           Pval=Pval)
      
      Common_table <- tibble(Class=sapply(strsplit(names(Pval) , split='.', fixed=TRUE), `[`, 2),
                             PWY=sapply(strsplit(names(Common) , split='.', fixed=TRUE), `[`, 1),
                             n.targ=Common)
      
      
    }
    
    # Add predicted target in class by TF
    Pval_table <- left_join(Pval_table, Common_table , by=c('Class', 'PWY'))
    
    # Add total predicted targets
    Pval_table <- left_join(Pval_table, Total_targets, by="Class")
    
    # add PWY size
    Pval_table <- left_join(Pval_table, CornCYC_size, by="PWY")
    
    Pval_table[,"TF"] <- tf
    
    # Select significant TFs 
    Pval_table <- subset(Pval_table, Pval <= 0.05)
    #Pval_table <- tibble(TF=c("test", "a", "b"), TF2="test2")
    
    #write.table(Pval_table, name, sep = "\t", row.names = F, quote = F)
    print(paste("... Done ...", sep = ""))
    
    return(Pval_table)
    
  }
  else {
    print(paste0(" Salado: ", tf, " .."))
  }
  
}

# test PWY enrichment within common TFs: CEN:GRN:GAN
TFs_3_layers_PWYEnr <- lapply(ven_file$`CEN:GRN:GAN`, PWYEnr_CommonTarg)
TFs_3_layers_PWYEnr <- TFs_3_layers_PWYEnr[unlist(lapply(TFs_3_layers_PWYEnr, function(x) is_tibble(x) ==TRUE))]
  
# test PWY enrichment within common TFs: CEN:GAN
TFs_2_layers_PWYEnr_CENGAN <- lapply(ven_file$`CEN:GAN`, PWYEnr_CommonTarg)
TFs_2_layers_PWYEnr_CENGAN <- TFs_2_layers_PWYEnr_CENGAN[unlist(lapply(TFs_2_layers_PWYEnr_CENGAN, function(x) is_tibble(x) ==TRUE))]
TFs_2_layers_PWYEnr_CENGAN <- as_tibble(rbindlist(TFs_2_layers_PWYEnr_CENGAN, idcol = F))

# test PWY enrichment within common TFs: CEN:GRN
TFs_2_layers_PWYEnr_CENGRN <- lapply(ven_file$`CEN:GRN`, PWYEnr_CommonTarg)
TFs_2_layers_PWYEnr_CENGRN <- TFs_2_layers_PWYEnr_CENGRN[unlist(lapply(TFs_2_layers_PWYEnr_CENGRN, function(x) is_tibble(x) ==TRUE))]
TFs_2_layers_PWYEnr_CENGRN <- as_tibble(rbindlist(TFs_2_layers_PWYEnr_CENGRN, idcol = F))

# test PWY enrichment within common TFs: GRN:GAN
TFs_2_layers_PWYEnr_GRNGAN <- lapply(ven_file$`GRN:GAN`, PWYEnr_CommonTarg)
TFs_2_layers_PWYEnr_GRNGAN <- TFs_2_layers_PWYEnr_GRNGAN[unlist(lapply(TFs_2_layers_PWYEnr_GRNGAN, function(x) is_tibble(x) ==TRUE))]
TFs_2_layers_PWYEnr_GRNGAN <- as_tibble(rbindlist(TFs_2_layers_PWYEnr_GRNGAN, idcol = F))

TFs_3_layers_PWYEnr <- left_join(TFs_3_layers_PWYEnr, unique(CornCYC[,c(1,2)]),  by=c("PWY"="Pathway.id"))

TFs_3_layers_PWYEnr
write.table(TFs_3_layers_PWYEnr, "CommTarg_Network_3TF.CornCYC.04.21.2022.txt", row.names = F, sep='\t', quote = F)
write.table(TFs_2_layers_PWYEnr_CENGAN, "CommTarg_Network_CENGAN.CornCYC.04.21.2022.txt", row.names = F, sep='\t', quote = F)
write.table(TFs_2_layers_PWYEnr_CENGRN, "CommTarg_Network_CENGRN.CornCYC.04.21.2022.txt", row.names = F, sep='\t', quote = F)
write.table(TFs_2_layers_PWYEnr_GRNGAN, "CommTarg_Network_GRNGAN.CornCYC.04.21.2022.txt", row.names = F, sep='\t', quote = F)

TFs_3_layers_PWYEnr[TFs_3_layers_PWYEnr$TF =="Zm00001d050195",]


###############################################################

###############################################################
########              Enrichment of GOs                ########
###############################################################

# target list based on layers classes
TargetsClasses$Zm00001d010805



SuperGO_CommomTarg("Zm00001d031796")

##############################################################
## 1 test GO enrichment within common TFs: CEN:GRN:GAN layers
##############################################################
Lgenes <- length(ven_file$`CEN:GRN:GAN`)
GenesList <- ven_file$`CEN:GRN:GAN`

w=40 # Size of range to test
print(".. Ready to start ..")

for (i in seq(0, Lgenes, w)){
  max=Lgenes
  Start=i+1
  end=i+w
  if (end<max){
    listtotest <- GenesList[Start:end]
    print(length(listtotest))
    
    mclapply(listtotest, SuperGO_CommomTarg, mc.cores=w)
    
  }
  else{
    r <- paste(Start, max, sep = "-")
    listtotest <- GenesList[Start:max]
    print(length(listtotest))
    mclapply(listtotest, SuperGO_CommomTarg, mc.cores=w)
  }
}
##############################################################

##############################################################
## 2 test PWY enrichment within common TFs: CEN:GAN
##############################################################

Lgenes <- length(ven_file$`CEN:GAN`)
GenesList <- ven_file$`CEN:GAN`

w=40 # Size of range to test
print(".. Ready to start ..")

for (i in seq(0, Lgenes, w)){
  max=Lgenes
  Start=i+1
  end=i+w
  if (end<max){
    listtotest <- GenesList[Start:end]
    print(length(listtotest))
    
    mclapply(listtotest, SuperGO_CommomTarg, mc.cores=w)
    
  }
  else{
    r <- paste(Start, max, sep = "-")
    listtotest <- GenesList[Start:max]
    print(length(listtotest))
    mclapply(listtotest, SuperGO_CommomTarg, mc.cores=w)
  }
}

##############################################################

##############################################################
## 3 test PWY enrichment within common TFs: GRN:GAN
##############################################################

Lgenes <- length(ven_file$`GRN:GAN`)
GenesList <- ven_file$`GRN:GAN`

w=10 # Size of range to test
print(".. Ready to start ..")

for (i in seq(0, Lgenes, w)){
  max=Lgenes
  Start=i+1
  end=i+w
  if (end<max){
    listtotest <- GenesList[Start:end]
    print(length(listtotest))
    
    mclapply(listtotest, SuperGO_CommomTarg, mc.cores=w)
    
  }
  else{
    r <- paste(Start, max, sep = "-")
    listtotest <- GenesList[Start:max]
    print(length(listtotest))
    mclapply(listtotest, SuperGO_CommomTarg, mc.cores=w)
  }
}
##############################################################

##############################################################
## 4 test PWY enrichment within common TFs: CEN:GRN
##############################################################

Lgenes <- length(ven_file$`CEN:GRN`)
GenesList <- ven_file$`CEN:GRN`

w=9 # Size of range to test
print(".. Ready to start ..")

for (i in seq(0, Lgenes, w)){
  max=Lgenes
  Start=i+1
  end=i+w
  if (end<max){
    listtotest <- GenesList[Start:end]
    print(length(listtotest))
    
    mclapply(listtotest, SuperGO_CommomTarg, mc.cores=w)
    
  }
  else{
    r <- paste(Start, max, sep = "-")
    listtotest <- GenesList[Start:max]
    print(length(listtotest))
    mclapply(listtotest, SuperGO_CommomTarg, mc.cores=w)
  }
}

##############################################################



######################################################################################
################             test common Targets amount layers        ################
######################################################################################

# List of targets by TF
CEN_list <- split(CoExp$Target, CoExp$Source)
GRN_list <- split(PDI$Target, PDI$Source)
GAN_list <- split(teQTLtf$Target, teQTLtf$Source)

CEN_list <- lapply(CEN_list, unique)
GRN_list <- lapply(GRN_list, unique)
GAN_list <- lapply(GAN_list, unique)

# CEN_GRN
Venn_CEN_GRN <- lapply(ven_file$`CEN:GRN`, overlap_significance_tf2)
names(Venn_CEN_GRN) <- ven_file$`CEN:GRN`

Venn_CEN_GRN <- as.data.frame(t(as.data.frame(Venn_CEN_GRN, stringsAsFactors = F)))
Venn_CEN_GRN["Group"] <- "CEN~GRN" 
Venn_CEN_GRN["TF"] <- ven_file$`CEN:GRN`
dim(Venn_CEN_GRN)

# GRN_GAN
Venn_GRN_GAN <- lapply(ven_file$`GRN:GAN`, overlap_significance_tf2)
names(Venn_GRN_GAN) <- ven_file$`GRN:GAN`
Venn_GRN_GAN <- as.data.frame(t(as.data.frame(Venn_GRN_GAN, stringsAsFactors = F)))
Venn_GRN_GAN["Group"] <- "GRN~GAN" 
Venn_GRN_GAN["TF"] <- ven_file$`GRN:GAN`
dim(Venn_GRN_GAN)

# CEN_GAN
Venn_CEN_GAN <- lapply(ven_file$`CEN:GAN`, overlap_significance_tf2)
names(Venn_CEN_GAN) <- ven_file$`CEN:GAN`
Venn_CEN_GAN <- as.data.frame(t(as.data.frame(Venn_CEN_GAN, stringsAsFactors = F)))
Venn_CEN_GAN["Group"] <- "CEN~GAN"
Venn_CEN_GAN["TF"] <- ven_file$`CEN:GAN`
dim(Venn_CEN_GAN)

# CEN_GAN_GRN
Venn_CEN_GAN_GRN <- lapply(ven_file$`CEN:GRN:GAN`, overlap_significance_tf3)
names(Venn_CEN_GAN_GRN) <- ven_file$`CEN:GRN:GAN`
Venn_CEN_GAN_GRN <- as.data.frame(t(as.data.frame(Venn_CEN_GAN_GRN, stringsAsFactors = F)))
Venn_CEN_GAN_GRN["Group"] <- "CEN~GAN~GRN"
Venn_CEN_GAN_GRN["TF"] <- ven_file$`CEN:GRN:GAN`
dim(Venn_CEN_GAN_GRN)

# CEN_GAN_GRN_cengrn
Venn_CEN_GAN_GRNcen_grn <- lapply(ven_file$`CEN:GRN:GAN`, overlap_significance_tf2_cen_grn)
names(Venn_CEN_GAN_GRNcen_grn) <- ven_file$`CEN:GRN:GAN`
Venn_CEN_GAN_GRNcen_grn <- as.data.frame(t(as.data.frame(Venn_CEN_GAN_GRNcen_grn, stringsAsFactors = F)))
Venn_CEN_GAN_GRNcen_grn["Group"] <- "cen~grn"
Venn_CEN_GAN_GRNcen_grn["TF"] <- ven_file$`CEN:GRN:GAN`
dim(Venn_CEN_GAN_GRNcen_grn)

# CEN_GAN_GRN_cengan
Venn_CEN_GAN_GRNcen_gan <- lapply(ven_file$`CEN:GRN:GAN`, overlap_significance_tf2_cen_gan)
names(Venn_CEN_GAN_GRNcen_gan) <- ven_file$`CEN:GRN:GAN`
Venn_CEN_GAN_GRNcen_gan <- as.data.frame(t(as.data.frame(Venn_CEN_GAN_GRNcen_gan, stringsAsFactors = F)))
Venn_CEN_GAN_GRNcen_gan["Group"] <- "cen~gan"
Venn_CEN_GAN_GRNcen_gan["TF"] <- ven_file$`CEN:GRN:GAN`
dim(Venn_CEN_GAN_GRNcen_gan)

# CEN_GAN_GRN_gangrn
Venn_CEN_GAN_GRNgan_grn <- lapply(ven_file$`CEN:GRN:GAN`, overlap_significance_tf2_gan_grn)
names(Venn_CEN_GAN_GRNgan_grn) <- ven_file$`CEN:GRN:GAN`
Venn_CEN_GAN_GRNgan_grn <- as.data.frame(t(as.data.frame(Venn_CEN_GAN_GRNgan_grn, stringsAsFactors = F)))
Venn_CEN_GAN_GRNgan_grn["Group"] <- "grn~gan"
Venn_CEN_GAN_GRNgan_grn["TF"] <- ven_file$`CEN:GRN:GAN`
dim(Venn_CEN_GAN_GRNgan_grn)

#
dim(subset(Venn_CEN_GAN_GRN, V2 <=0.05))
dim(subset(Venn_CEN_GAN_GRN, V2 <=0.05))

######################################################################################

Venn_DF <- as_tibble(rbind(Venn_CEN_GRN, Venn_GRN_GAN, Venn_CEN_GAN, Venn_CEN_GAN_GRN))

# To fix scale maximum 
Venn_DF$V2[(Venn_DF$V2 <= min(Venn_CEN_GAN$V2))] <- min(Venn_CEN_GAN$V2) 

Venn_DF_CEN_GAN_GRNpairs <- as_tibble(rbind(Venn_CEN_GAN_GRNcen_grn, 
                                            Venn_CEN_GAN_GRNcen_gan, 
                                            Venn_CEN_GAN_GRNgan_grn))
# To fix scale maximum 
Venn_DF_CEN_GAN_GRNpairs$V2[(Venn_DF_CEN_GAN_GRNpairs$V2 <= 1e-20)] <- 1e-20

########################################################################
########             Significane of commons targets             ########
########################################################################

Venn_DF$Group <- factor(Venn_DF$Group, levels = rev(c('CEN_GRN', "GRN_GAN", "CEN_GAN", "CEN_GAN_GRN")))


Plot_2 <- ggplot(Venn_DF, aes(y=Group, x=-log10(V2), color = (V2<=0.05))) +
  geom_jitter(width = 0.2, size=0.5) +
  scale_color_manual(values=c("TRUE" = "#F8766D", "FALSE" = "#00BFC4"), 
                     labels=c(expression("P-value" <= 0.05), "P-value > 0.05")) +
  xlab(expression(-Log[10]~Pvalue)) + ylab("Venn diagram groups") +
  theme_pubclean()


Venn_DF_CEN_GAN_GRNpairs$Group <- factor(Venn_DF_CEN_GAN_GRNpairs$Group, levels = rev(c('grn~gan', "cen~gan", "cen~grn")))
Plot_3 <- ggplot(Venn_DF_CEN_GAN_GRNpairs, aes(y=Group, x=-log10(V2), color = (V2<=0.05))) +
  geom_jitter(width = 0.2, size=0.5) +
  scale_color_manual(values=c("TRUE" = "#F8766D", "FALSE" = "#00BFC4"), 
                     labels=c(expression("P-value" <= 0.05), "P-value > 0.05")) +
  xlab(expression(-Log[10]~Pvalue)) + ylab("") +
  theme_pubclean()

Plot_2 <- ggpar(Plot_2, font.tickslab = 10, font.x = 10, font.y = 10, legend = "bottom", legend.title = "")
Plot_3 <- ggpar(Plot_3, font.tickslab = 10, font.x = 10, font.y = 10, legend = "bottom", legend.title = "")

# venn figure: Total TFs 
Plot_1 <- vennfunc(Total_TFs_list) + theme(text = element_text(size=8))

ggarrange(Plot_1 , Plot_2, Plot_3, ncol = 3, 
          widths = c(1.5,1.5,1), 
          common.legend = T, legend = "bottom")

# table(subset(Venn_DF, V2 <=0.05)$Group)
# table(subset(Venn_DF_CEN_GAN_GRNpairs, V2 <=0.05)$Group)

################################################################

################################################################
########   Count common TFs and number of associations  ######## 
################################################################

#
PDI_counts <- as_tibble(as.data.frame(table(PDI$Source), stringsAsFactors = F))
colnames(PDI_counts) <- c("Source", "Targets")
#
CoExp_counts <- as_tibble(as.data.frame(table(CoExp$Source), stringsAsFactors = F))
colnames(CoExp_counts) <- c("Source", "Targets")
#
teQTL_counts <- as_tibble(as.data.frame(table(teQTL$Source), stringsAsFactors = F))
colnames(teQTL_counts) <- c("Source", "Targets")

teQTL_countsTFs <- subset(teQTL_counts, Source %in% c(PDI_counts$Source, CoExp_counts$Source, TF_CoR$GeneID))

write.table(teQTL_countsTFs[,1], "~/Projects/MaizeENCODE/Data_45_net/wPCC_net_only_TFs/TF_list.txt", sep = "\t", quote = F, row.names = F)
################################################################

#############################################################################################

######################################################################
###        test with TFs filter based on number of targets         ###
######################################################################

Total_TFs_list = list(CEN=CoExp_counts$Source, 
                      GRN=PDI_counts$Source,
                      GAN=teQTL_countsTFs$Source)

Total_TFs_list_10 = list(CEN=subset(CoExp_counts, Targets>=10)$Source, 
                         GRN=subset(PDI_counts, Targets>=10)$Source,  
                         GAN=subset(teQTL_countsTFs, Targets>=10)$Source)

# size= 5x5
ggarrange(vennfunc(Total_TFs_list) + labs(title="Total TFs/CoRegs"),
          vennfunc(Total_TFs_list_10) +  labs(title=expression("Total TFs/CoRegs (">=" 10 targets)")), 
          nrow=1)

######################################################################

###################################################################### 
########                Add PCC to GAN network                ######## 
######################################################################
# TF names
TFdic <- as_tibble(read.table("Data/Annotations/TF_Id_Name.txt", h=T, stringsAsFactors = F))
tail(TFdic)

# List of TFs with CoExp data
TF_in_wPCCDB <- as.character(read.table("wPCC_TFs_Done.txt", h=F)$V1)
TF_in_wPCCDB <- as_tibble(as.data.frame(table(TF_in_wPCCDB), stringsAsFactors = F)) 

# Number of CoExp Nets
TF_in_wPCCDB <- TF_in_wPCCDB[order(-TF_in_wPCCDB$Freq),] 
hist(TF_in_wPCCDB$Freq, 45)

# Get observed PCC mean by TF
GAN_wPCC <- lapply(TF_in_wPCCDB$TF_in_wPCCDB, add_PCC.GAN)
names(GAN_wPCC) <- TF_in_wPCCDB$TF_in_wPCCDB

# Get mean wPCC by TF
GAN_wPCCmean <- lapply(GAN_wPCC, Get_PCCmeanByNet)


# Get Z and Pval values comparing with random networks
Result_Ztest <- lapply(TF_in_wPCCDB$TF_in_wPCCDB, Ztes_for_list)
names(Result_Ztest) <- TF_in_wPCCDB$TF_in_wPCCDB

Result_Ztest <- as_tibble(rbindlist(Result_Ztest, idcol = T))

Result_Ztest_Freq <- as_tibble(as.data.frame(table(subset(Result_Ztest, Pval <=0.05)$.id), stringsAsFactors = F))
Result_Ztest_Freq <- Result_Ztest_Freq[order(-Result_Ztest_Freq$Freq),] 
Result_Ztest_Freq

tf="Zm00001d048603"


MakeBoxPlot_RandomPCC <- function(tf){
  ## Read random wPCC means
  rdb <- as_tibble(fread(paste0("GAN_R_CoExpData/MeanVals/wPCCm_GAN_Random_",tf, ".txt")))
  colnames(rdb)[2] <- "PCC"
  
  obs <- as.data.frame(GAN_wPCCmean[[tf]])
  obs["Net"]<- row.names(obs)
  colnames(obs)[1] <- "PCC"
  
  tem <- subset(Result_Ztest, .id==tf)[,c("Pval", "Nets")]
  
  obs <- left_join(obs, tem, by=c("Net"="Nets"))
  
  Plot <- ggplot(rdb, aes(y=Net, x=PCC)) +
    geom_boxplot(notch = T, outlier.size = 0.5) +
    geom_point(data=obs, aes(y=Net, x=PCC, color=Pval<=0.05), width = 0.01, 
               size=2, alpha=0.9) +
    theme_pubclean() +
    scale_color_discrete(limits = c("FALSE", "TRUE"), labels = c("P > 0.05 ", expression("P"<=" 0.5"))) +
    labs(subtitle=ReplaceName(tf))
  
  Plot <- ggpar(Plot, font.tickslab=10, font.x = 10, font.y = 10, legend = 'bottom', ylab = "Dataset", legend.title = "Obs. PCC")
  return(Plot)
  
}

# GRAS26 Size 6x3
MakeBoxPlot_RandomPCC("Zm00001d048603")
GAN_wPCC$Zm00001d048603[,c("GeneID","n18e_2", "n18e_1", "n18e", "n17a_2")]

Result_Ztest <- Result_Ztest[is.na(Result_Ztest$Pval)==FALSE,]
Result_Ztest$Pval[(Result_Ztest$Pval)<=1e-10] <- 1e-10

# length(unique(teQTLtf$Source))
#Result_Ztest_summary <- as_tibble(as.data.frame(table(Result_Ztest$Pval <= 0.05, Result_Ztest$Nets), stringsAsFactors = F))
Result_Ztest_summary <- as_tibble(as.data.frame(table(subset(Result_Ztest, Pval <= 0.05)$.id), stringsAsFactors = F))
colnames(Result_Ztest_summary) <- c("TF", "Datasets")

# Summary all TFs that pass by Network

Result_Ztest_summary

Plot_GAN_wpcc <- ggplot(as.data.frame(table(Result_Ztest_summary$Datasets)), 
                        aes(x=Freq, y=Var1)) +
  geom_bar(stat="identity", fill="gray80")+
  geom_text(aes(x=Freq+2, y=Var1, label=Freq)) +
  theme_pubclean() + 
  scale_x_continuous(expand = c(0, 0), limits = c(0,480)) +
  ylab("Dataset") + xlab("TF")

Plot_GAN_wpcc <- ggpar(Plot_GAN_wpcc, 
                       font.tickslab=10, 
                       font.x = 10, font.y = 10)
Plot_GAN_wpcc

# Plot_GAN_wpcc <- ggplot(Result_Ztest_summary, aes(y=Nets, x= -log10(Pval), color= Pval<=0.05)) +
#   geom_point(size=0.2, alpha=0.7) +
#   theme_pubclean() + 
#   scale_x_continuous(expand = c(0, 0))+
#   scale_color_discrete(limits = c("FALSE", "TRUE"), labels = c("> 0.05 ", expression(""<=" 0.5"))) +
#   xlab(expression(-Log[10]~Pvalue))





######################################################################

######################################################################
########                                                      ########
########   Summary/general number of PWY enrichment results   ######## 
########                                                      ########
######################################################################


###############################################################
# Count number of TFs by PWY 
###############################################################

CoExp_PWY_freq <- as_tibble(as.data.frame(table(unique(CoExp_PWY[,c(2,6)])$PWY)))
CoExp_PWY_freq["Class"] = 'CEN'
CoExp_PWY_freq["Rank"] <- rank(-CoExp_PWY_freq$Freq)
CoExp_PWY_freq <- CoExp_PWY_freq[order(CoExp_PWY_freq$Rank),]

PDI_PWY_freq <- as_tibble(as.data.frame(table(unique(PDI_PWY[,c(2,6)])$PWY)))
PDI_PWY_freq["Class"] = 'GRN'
PDI_PWY_freq["Rank"] <- rank(-PDI_PWY_freq$Freq)
PDI_PWY_freq <- PDI_PWY_freq[order(PDI_PWY_freq$Rank),]

teQTL_PWY_freq <- as_tibble(as.data.frame(table(unique(teQTL_PWY[,c(2,6)])$PWY)))
teQTL_PWY_freq["Class"] = 'GAN'
teQTL_PWY_freq["Rank"] <- rank(-teQTL_PWY_freq$Freq)
teQTL_PWY_freq <- teQTL_PWY_freq[order(teQTL_PWY_freq$Rank),]


PWY_Freq_by_net <- rbind(CoExp_PWY_freq, PDI_PWY_freq, teQTL_PWY_freq)
PWY_Freq_by_net$Class <- factor(PWY_Freq_by_net$Class, levels = c("CEN", "GRN", "GAN"))

PWY_Freq_by_net["PWYname"] <- ReplaceNamePWY(PWY_Freq_by_net$Var1)
PWY_Freq_by_net$PWYname <- gsub("_", " ", PWY_Freq_by_net$PWYname)

Total_tf_pwys_net <- unique(rbind(CoExp_PWY[,c(2,6,7)], 
                                  PDI_PWY[,c(2,6,7)], 
                                  subset(teQTL_PWY, TFid %in% All_TFs)[,c(2,6,7)]))
Total_tf_pwys_net
length(unique(Total_tf_pwys_net$PWY))
length(unique(Total_tf_pwys_net$TFid))

table(Total_tf_pwys_net$Net)

Plot_Pwys_net <- ggplot(PWY_Freq_by_net, aes(x = Freq, y=Class, fill=Class)) + 
  geom_boxplot(notch = T, alpha=0.5, outlier.shape = NA) +
  geom_jitter(data=subset(PWY_Freq_by_net, Rank <=6), aes(x = Freq, y=Class),
              alpha=0.5, height = 0.01) + 
  geom_text_repel(data=subset(PWY_Freq_by_net, Rank <=6), 
                  aes(x=Freq, y=Class, label=PWYname),
                  nudge_y= 0.2,
                  nudge_x= 0.2,
                  direction= "y",
                  force_pull=0,
                  hjust        = 0,
                  segment.size = 0.2,
                  size=2) +
  scale_fill_manual(values=c("goldenrod1", "darkorchid1", 'dodgerblue')) +
  theme_pubclean() +
  scale_x_sqrt(breaks=c(1,5,20,40,80, 120, 160), expand=c(0,0)) + 
  
  xlab('Counts') 
#xlab( expression(paste(Log[2], '(Counts + 1)', sep = " "))) +
#xlim(0, 13)

Plot_Pwys_net
###############################################################

###############################################################
# Count number of PWYs by TF    
###############################################################

CoExp_TF_freq <- as_tibble(as.data.frame(table(unique(CoExp_PWY[,c(2,6)])$TFid)))
CoExp_TF_freq["Class"] = 'CEN'
CoExp_TF_freq["Rank"] <- rank(-CoExp_TF_freq$Freq)
CoExp_TF_freq <- CoExp_TF_freq[order(CoExp_TF_freq$Rank),]

PDI_TF_freq <- as_tibble(as.data.frame(table(unique(PDI_PWY[,c(2,6)])$TFid)))
PDI_TF_freq["Class"] = 'GRN'
PDI_TF_freq["Rank"] <- rank(-PDI_TF_freq$Freq)
PDI_TF_freq <- PDI_TF_freq[order(PDI_TF_freq$Rank),]


teQTL_TF_freq <- as_tibble(as.data.frame(table(unique(subset(teQTL_PWY, TFid %in% All_TFs)[,c(2,6)])$TFid)))
teQTL_TF_freq["Class"] = 'GAN'
teQTL_TF_freq["Rank"] <- rank(-teQTL_TF_freq$Freq)
teQTL_TF_freq <- teQTL_TF_freq[order(teQTL_TF_freq$Rank),]

TF_Freq_by_net <- rbind(CoExp_TF_freq, PDI_TF_freq, teQTL_TF_freq)
TF_Freq_by_net$Class <- factor(TF_Freq_by_net$Class, levels = c("CEN", "GRN", "GAN"))

TF_Freq_by_net["TFname"] <- ReplaceName(TF_Freq_by_net$Var1)

Toptfs <- subset(TF_Freq_by_net, Rank <=6)

Plot_TFs_net <- ggplot(TF_Freq_by_net, aes(x = Freq, y=Class, fill=Class)) + 
  geom_boxplot(notch = T, alpha=0.5) +
  geom_jitter(data=Toptfs, aes(x=Freq, y=Class), alpha=0.5, height = 0.01, size=1) + 
  geom_text_repel(data=Toptfs, aes(x=Freq, y=Class, label=TFname),
                  nudge_y= 0.2,
                  nudge_x= 0.4,
                  direction= "y",
                  force_pull=0,
                  hjust        = 0,
                  segment.size = 0.2,
                  size=2) +
  scale_fill_manual(values=c("goldenrod1", "darkorchid1", 'dodgerblue')) +
  theme_pubclean() +
  scale_x_sqrt(breaks=c(1,5,10,20,30), expand=c(0,0)) + 
  xlab("Counts")
#xlab( expression(paste(Log[2], '(Counts + 1)', sep = " "))) +
#xlim(0, 6)

# Figure 2a-b
Plot_Pwys_net <- ggpar(Plot_Pwys_net, font.tickslab=14, font.x = 14, font.y = 14, legend = 'none', ylab = "Network")
Plot_TFs_net <- ggpar(Plot_TFs_net, font.tickslab=14, font.x = 14, font.y = 14, legend = 'none', ylab = "Network")

# https://realpython.com/k-means-clustering-python/

PWY_Freq_by_net %>% 
  group_by(Class) %>%
  summarise(M=mean(Freq))

TF_Freq_by_net %>% 
  group_by(Class) %>%
  summarise(M=mean(Freq))

#subset(TF_Freq_by_net,  Class=='CEN')
#ReplaceNamePWY(PDI_PWY[PDI_PWY$TFid=="Zm00001d052405",]$PWY)
#ReplaceNamePWY(CoExp_PWY[CoExp_PWY$TFid=="Zm00001d052405",]$PWY)
###############################################################

###############################################################
##  Count Common PWYs among TFs
###############################################################

# Venn diagramm of common PWTs 
ResPWY_list = list(CEN=unique(CoExp_PWY$PWY), 
                   GRN=unique(PDI_PWY$PWY),
                   GAN=unique(subset(teQTL_PWY, TFid %in% All_TFs)$PWY))


# Venn diagramm of common TFs
ResTF_list = list(CEN=unique(CoExp_PWY$TFid), 
                  GRN=unique(PDI_PWY$TFid),
                  GAN=unique(subset(teQTL_PWY, TFid %in% TF_CoR$GeneID)$TFid))


Plot_AB <- ggarrange(Plot_Pwys_net, Plot_TFs_net, nrow = 2, align = 'hv')
Plot_CD <- ggarrange(vennfunc(ResPWY_list),
                     vennfunc(ResTF_list),
                     nrow = 2, align = 'hv')

# size= 5x8
ggarrange(Plot_AB, Plot_CD, ncol = 2, widths = c(1.3, 1))

###############################################################

###############################################################
##          Common PWYs by network                
###############################################################

pwystop <- unique(CoExp_PWY$PWY)[unique(CoExp_PWY$PWY) %in% unique(PDI_PWY$PWY)][unique(CoExp_PWY$PWY)[unique(CoExp_PWY$PWY) %in% unique(PDI_PWY$PWY)] %in% unique(subset(teQTL_PWY, TFid %in% TF_CoR$GeneID)$PWY)]
tfs.cen.grn.top <- unique(CoExp_PWY$TFid)[unique(CoExp_PWY$TFid) %in% unique(PDI_PWY$TFid)]
tfs.cen.gan.top <- unique(CoExp_PWY$TFid)[unique(CoExp_PWY$TFid) %in% unique(subset(teQTL_PWY, TFid %in% TF_CoR$GeneID)$TFid)]

#
ReplaceNamePWY(pwystop)
ReplaceName(tfs.cen.grn.top)[order(ReplaceName(tfs.cen.grn.top))]
ReplaceName(tfs.cen.gan.top)[order(ReplaceName(tfs.cen.gan.top))]

## TFs in PWYs highly predicted
TFsinPDI_toppwy   <- unique(subset(PDI_PWY, PWY %in% subset(PWY_Freq, Freq==3)$PWY)$TFid)
TFsinCoExp_toppwy <- unique(subset(CoExp_PWY, PWY %in% subset(PWY_Freq, Freq==3)$PWY)$TFid)

#
CoExpDFTop <- subset(CoExp_PWY, PWY %in% pwystop)
PDIDFTop   <- subset(PDI_PWY, PWY %in% pwystop)
teQTLDFTop <- subset(teQTL_PWY, PWY %in% pwystop)
teQTLDFTop <- subset(teQTLDFTop, TFid %in% TF_CoR$GeneID)


MakeHeatmap <- function(df_PWY) {
  
  # add PWY name
  df_PWY[,"Name"]  <- paste0("[", left_join(df_PWY, CornCYC_size, by="PWY")$PWYSize, "] ", ReplaceNamePWY(df_PWY$PWY))
  df_PWY$TF <- ReplaceName(df_PWY$TF)
  
  # add % of PWY targeted
  df_PWY[,"Per"] <- round((df_PWY$n.targ/left_join(df_PWY, CornCYC_size, by="PWY")$PWYSize)*100,2)
  
  
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
  Text <- as.matrix(Text)
  
  cor_scale <- colorRamp2(seq(0,99,1), viridis(100, direction = -1, option = "A"))
  
  hm_common <- Heatmap(t(Values), 
                       #heatmap_width = unit(15, "cm"), 
                       #heatmap_height = unit(8, "cm"), 
                       column_names_rot = 90, 
                       cell_fun = function(j, i, x, y, width, height, fill) {
                         if(t(Text)[i, j]>0){
                           grid.text(sprintf("%.f", t(Text)[i, j]), x, y, gp = gpar(fontsize = 12, col = "gray70"))
                         }
                         
                       },
                       name = "PWY Percentage", 
                       cluster_columns = TRUE, column_dend_reorder = TRUE,
                       show_row_dend = FALSE, show_column_dend = FALSE,
                       col=cor_scale,
                       column_names_gp = gpar(fontsize = 8),
                       row_names_gp = gpar(fontsize = 8),
                       show_heatmap_legend = T,
                       heatmap_legend_param = list(title_gp = gpar(fontsize = 10), 
                                                   labels_gp = gpar(fontsize = 10),
                                                   width = unit(15, "mm"),
                                                   direction = "horizontal")) 
  return(hm_common)
  
}
MakeHeatmapVertical <- function(df_PWY, h) {
  
  # add PWY name
  df_PWY[,"Name"]  <- paste0(ReplaceNamePWY(df_PWY$PWY), " [", left_join(df_PWY, CornCYC_size, by="PWY")$PWYSize, "]")
  df_PWY$TF <- ReplaceName(df_PWY$TF)
  
  # add % of PWY targeted
  df_PWY[,"Per"] <- round((df_PWY$n.targ/left_join(df_PWY, CornCYC_size, by="PWY")$PWYSize)*100,2)
  
  
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
  print(df2)
  
  Text <- reshape2::dcast(df2, TF ~ Name, value.var="n.targ")
  row.names(Text) <- Text$TF
  Text <- Text[,-c(1)]
  Text[is.na(Text)] <- 0
  Text <- as.matrix(t(Text))
  
  cor_scale <- colorRamp2(seq(0,99,1), viridis(100, direction = -1, option = "A"))
  
  colnames(Values) <- gsub("_", " ", colnames(Values))
  
  hm_common <- Heatmap(Values, 
                       #heatmap_width = unit(15, "cm"), 
                       #heatmap_height = unit(8, "cm"), 
                       column_names_rot = 70, 
                       cell_fun = function(j, i, x, y, width, height, fill) {
                         if(t(Text)[i, j]>0){
                           grid.text(sprintf("%.f", t(Text)[i, j]), x, y, gp = gpar(fontsize = 8, col = "gray60"))
                         }
                         
                       },
                       name = "PWY Percentage", 
                       cluster_columns = TRUE, column_dend_reorder = TRUE,
                       show_row_dend = FALSE, show_column_dend = FALSE,
                       col=cor_scale,
                       width = unit(2, "cm"),
                       height = unit(h, "cm"),
                       column_names_gp = gpar(fontsize = 8),
                       row_names_gp = gpar(fontsize = 8),
                       show_heatmap_legend = T,
                       heatmap_legend_param = list(title_gp = gpar(fontsize = 10), 
                                                   labels_gp = gpar(fontsize = 10),
                                                   width = unit(15, "mm"),
                                                   direction = "horizontal")) 
  return(hm_common)
  
}

CoExp_heatmap <- MakeHeatmapVertical(CoExpDFTop, 20)
PDI_heatmap <- MakeHeatmapVertical(PDIDFTop, 3)
teQTL_heatmap <- MakeHeatmapVertical(teQTLDFTop, 1)

list_hm = CoExp_heatmap %v% PDI_heatmap %v% teQTL_heatmap
list_hm

# size 20x6
draw(list_hm, heatmap_legend_side = "bottom")  
###############################################################

##################################################
##  Common TFs based on common PWY's analysis   ##
##################################################

CommonTFs_in_Top_PWYs <- as.data.frame(table(c(unique(CoExpDFTop$TFid), unique(PDIDFTop$TFid), unique(teQTLDFTop$TFid))))
CommonTFs_in_Top_PWYs <- subset(CommonTFs_in_Top_PWYs, Freq>1)
ReplaceName(CommonTFs_in_Top_PWYs$Var1)


CommonTFs_in_Top5_pwy <- unique(subset(CoExp_PWY, PWY %in% subset(PWY_Freq, Freq==3)$PWY)$TFid[(subset(CoExp_PWY, PWY %in% subset(PWY_Freq, Freq==3)$PWY)$TFid %in% subset(PDI_PWY, PWY %in% subset(PWY_Freq, Freq==3)$PWY)$TFid)])
ReplaceName(CommonTFs_in_Top5_pwy)
ReplaceName("Zm00001d020492")

CoExpDFTop[CoExpDFTop$TFid %in% CommonTFs_in_Top5_pwy,]
PDIDFTop[PDIDFTop$TFid %in% CommonTFs_in_Top5_pwy,]


GetTF_PWY_Net <- function(tf, pwy) {
  # Annotate pwy genes
  pwy.targ <- unique(subset(CornCYC, Pathway.id == pwy)[,1:3])
  
  # Mark targets
  pwy.targ["CEN"] <- (pwy.targ$GeneID %in% subset(CoExp, Source==tf)$Target)*1
  pwy.targ["GRN"] <- (pwy.targ$GeneID %in% subset(PDI, Source==tf)$Target)*1
  pwy.targ["GAN"] <- (pwy.targ$GeneID %in% subset(teQTL, Source==tf)$Target)*1
  #
  return(pwy.targ)
}


write.table(GetTF_PWY_Net("Zm00001d020492", "PWY_5048"), "Zm00001d020492_PWY_5048.net.txt", 
            sep = "\t", row.names = F, quote = F)

write.table(GetTF_PWY_Net("Zm00001d006236", "PWY_6457"), "Zm00001d006236_PWY_6457.net.txt", 
            sep = "\t", row.names = F, quote = F)

write.table(GetTF_PWY_Net("Zm00001d006236", "PWY1F_467"), "Zm00001d006236_PWY1F_467.net.txt", 
            sep = "\t", row.names = F, quote = F)
##################################################

###############################################################
########  Details about Common TFs by network           
###############################################################

# 1. Identification of Common TFs between CoExp & PDI
CEN_GRN_ids <- unique(CoExp_PWY$TFid[CoExp_PWY$TFid %in% PDI_PWY$TFid])


# 2. Identification of Common TFs between CoExp & teQTL
CEN_GAN_ids <- unique(CoExp_PWY$TFid[CoExp_PWY$TFid %in% teQTL_PWY$TFid])

# 3. Count targets by significant PWY

Count_CEN_GRN_in_pwy <- function(tfid){
  
  # pdi targets
  pdi <- unique(subset(PDI, Source==tfid)$Target)
  
  # CoExp targets
  coexp <- unique(subset(CoExp, Source==tfid)$Target)
  
  # PWY significantly enriched
  pdi_PWY <- unique(subset(PDI_PWY, TFid==tfid)$PWY)
  coexp_PWY <- unique(subset(CoExp_PWY, TFid==tfid)$PWY)
  
  # Count targets by PWY
  counts_pdi <- unique(subset(CornCYC, Pathway.id %in% pdi_PWY & GeneID %in% pdi)[,c(1,3)])
  counts_coexp <- unique(subset(CornCYC, Pathway.id %in% coexp_PWY & GeneID %in% coexp)[,c(1,3)])
  #
  counts_pdi <- as_tibble(as.data.frame(table(counts_pdi$Pathway.id)))
  counts_coexp <- as_tibble(as.data.frame(table(counts_coexp$Pathway.id)))
  #
  counts_pdi["Net"] <- "GRN"
  counts_coexp["Net"] <- "CEN"
  
  # merge results
  counts <- rbind(counts_pdi, counts_coexp)
  colnames(counts) <- c("PWY", "n.targ", "Net")
  counts["TFid"] <- tfid
  counts["TFname"] <- ReplaceName(counts$TFid)
  counts["PWYname"] <- ReplaceNamePWY(counts$PWY)
  
  # add pwy size
  counts <- left_join(counts, CornCYC_size, by="PWY")
  
  # add pwy targ. percentage
  counts["n.targ_perc"] <- (counts$n.targ/counts$PWYSize)*100
  
  
  return(counts)
  
  
}

Count_CEN_GAN_in_pwy <- function(tfid){
  
  # pdi targets
  teqtl <- unique(subset(teQTL, Source==tfid)$Target)
  
  # CoExp targets
  coexp <- unique(subset(CoExp, Source==tfid)$Target)
  
  # PWY significantly enriched
  teqtl_PWY <- unique(subset(teQTL_PWY, TFid==tfid)$PWY)
  coexp_PWY <- unique(subset(CoExp_PWY, TFid==tfid)$PWY)
  
  # Count targets by PWY
  counts_teqtl <- unique(subset(CornCYC, Pathway.id %in% teqtl_PWY & GeneID %in% teqtl)[,c(1,3)])
  counts_coexp <- unique(subset(CornCYC, Pathway.id %in% coexp_PWY & GeneID %in% coexp)[,c(1,3)])
  #
  counts_teqtl <- as_tibble(as.data.frame(table(counts_teqtl$Pathway.id)))
  counts_coexp <- as_tibble(as.data.frame(table(counts_coexp$Pathway.id)))
  #
  counts_teqtl["Net"] <- "GAN"
  counts_coexp["Net"] <- "CEN"
  
  # merge results
  counts <- rbind(counts_teqtl, counts_coexp)
  colnames(counts) <- c("PWY", "n.targ", "Net")
  counts["TFid"] <- tfid
  counts["TFname"] <- ReplaceName(counts$TFid)
  counts["PWYname"] <- ReplaceNamePWY(counts$PWY)
  
  # add pwy size
  counts <- left_join(counts, CornCYC_size, by="PWY")
  
  # add pwy targ. percentage
  counts["n.targ_perc"] <- (counts$n.targ/counts$PWYSize)*100
  
  
  return(counts)
  
  
}

common39_TF <- lapply(CEN_GRN_ids, Count_CEN_GRN_in_pwy)
common7_TF <- lapply(CEN_GAN_ids, Count_CEN_GAN_in_pwy)

common39_TF <- as_tibble(rbindlist(common39_TF, idcol = F))
common7_TF <- as_tibble(rbindlist(common7_TF, idcol = F))

common39_TF <- common39_TF[order(common39_TF$TFname),]
common7_TF <- common7_TF[order(common7_TF$TFname),]


# bar plot
make_bar_Plot <- function(df, colors) {
  
  tfs <- sort(unique(df$TFname))
  
  out <- list()
  for (t in tfs) {
    df_tem <- subset(df, TFname ==t)
    
    Plot <- ggplot(df_tem, aes(y=PWYname, x=n.targ_perc, fill=Net))+
      geom_col(alpha=0.5, position = "dodge") +
      scale_x_continuous(expand=c(0,0), limits = c(0, 100)) + 
      theme_pubclean() + 
      theme(strip.text.x = element_text(size = 8), 
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      scale_fill_manual(values=colors) +
      ylab("") + xlab("PWY percentage")
    
    Plot <- ggpar(Plot, font.tickslab = 9, subtitle = t, legend.title = 'Network')
    out[[t]] <- Plot
  }
  
  return(out)
}

plot_l_common7_TF <- make_bar_Plot(common7_TF, c("goldenrod1", 'dodgerblue'))
plot_l_common39_TF <- make_bar_Plot(common39_TF, c("goldenrod1", 'darkorchid1'))

# size=12x15
Fig_S6 <- ggarrange(plotlist=plot_l_common7_TF, ncol = 2, nrow = 4, 
                    common.legend = T, align = 'v', legend = "bottom")

# size=35x20
Fig_S7 <- ggarrange(plotlist=plot_l_common39_TF, 
                    ncol = 3, 
                    nrow = 13, 
                    common.legend = T, align = 'v', legend = "bottom")


ggsave(filename = "Common_CEN_GRN_39TFs.pdf", plot = Fig_S7, 
       width = 25, height = 35, units = "in")


write.table(common39_TF, 'CEN_GRN_Commom.39TFs.txt', sep = '\t', quote = F, row.names = F)
write.table(common7_TF, 'CEN_GAN_Commom.7TFs.txt', sep = '\t', quote = F, row.names = F)


# 2. heat-map with common TFs in top 5 PWYs
CoExp_heatmapTopTFs <- MakeHeatmap(subset(CoExpDFTop, TFid %in% CommonTFs_in_Top5_pwy))
PDI_heatmapTopTFs <- MakeHeatmap(subset(PDIDFTop, TFid %in% CommonTFs_in_Top5_pwy))

draw(CoExp_heatmapTopTFs+PDI_heatmapTopTFs, heatmap_legend_side = "bottom")  


Top5_pwy_Candiates <- rbind(tibble(layer="PDI",TFid=subset(PDI_PWY, PWY %in% subset(PWY_Freq, Freq==3)$PWY)$TFid),
                            tibble(layer="teQTL",TFid=subset(teQTL_PWY, PWY %in% subset(PWY_Freq, Freq==3)$PWY)$TFid),
                            tibble(layer="CoExp",TFid=subset(CoExp_PWY, PWY %in% subset(PWY_Freq, Freq==3)$PWY)$TFid))
write.table(Top5_pwy_Candiates, "Top5_pwy_Candiates.txt", row.names = F, quote = F, sep = "\t")


################################
# Plots for Lab meeting
################################

CEN_GRN_name <- ReplaceName(CEN_GRN_ids)
TopFromVisualanalysis <- c("ARF18", "bHLH91", "bHLH43", "COL8", "GLK9",  "MYB31", "WRKY53")
CEN_GRN_idsTop <- CEN_GRN_ids[CEN_GRN_name %in% TopFromVisualanalysis]

CEN_GAN_name <- ReplaceName(CEN_GAN_ids)
TopFromVisualanalysis2 <- c("DOF19", "SBP29")
CEN_GAN_idsTop <- CEN_GAN_ids[CEN_GAN_name %in% TopFromVisualanalysis2]


Top_CEN_GRN_common_TF <- lapply(CEN_GRN_idsTop, Count_CEN_GRN_in_pwy)
Top_CEN_GAN_common_TF <- lapply(CEN_GAN_idsTop, Count_CEN_GAN_in_pwy)

Top_CEN_GRN_common_TF <- as_tibble(rbindlist(Top_CEN_GRN_common_TF, idcol = F))
Top_CEN_GAN_common_TF <- as_tibble(rbindlist(Top_CEN_GAN_common_TF, idcol = F))

Top_CEN_GRN_common_TF <- Top_CEN_GRN_common_TF[order(Top_CEN_GRN_common_TF$TFname),]
Top_CEN_GAN_common_TF <- Top_CEN_GAN_common_TF[order(Top_CEN_GAN_common_TF$TFname),]

make_bar_Plot <- function(df, colors) {
  
  tfs <- sort(unique(df$TFname))
  
  out <- list()
  for (t in tfs) {
    df_tem <- subset(df, TFname ==t)
    
    Plot <- ggplot(df_tem, aes(y=PWYname, x=n.targ_perc, fill=Net))+
      geom_col(alpha=0.5, position = "dodge") +
      scale_x_continuous(expand=c(0,0), limits = c(0, 100)) + 
      theme_pubclean() + 
      theme(strip.text.x = element_text(size = 10), 
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      scale_fill_manual(values=colors) +
      ylab("") + xlab("PWY percentage")
    
    Plot <- ggpar(Plot, font.tickslab = 10, subtitle = t, legend.title = 'Network')
    out[[t]] <- Plot
  }
  
  return(out)
}

plot_Top_CEN_GRN <- make_bar_Plot(Top_CEN_GRN_common_TF, c("goldenrod1", 'darkorchid1'))

plot_Top_CEN_GAN <- make_bar_Plot(Top_CEN_GAN_common_TF, c("goldenrod1", 'dodgerblue'))


plot_Top_CEN_GRN$DOF19 <- plot_Top_CEN_GAN$DOF19
plot_Top_CEN_GRN$SBP29 <- plot_Top_CEN_GAN$SBP29


ggarrange(plotlist=plot_Top_CEN_GRN, 
          ncol = 3, 
          nrow = 3, 
          common.legend = T, align = 'hv', legend = "bottom")

##########################################################

################################################
####       Heatmap TF vs pathway
################################################

