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
  # to test the significance of multiple gene sets (>2)
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

Get_PCCmeanByNet <- function(df){
  apply(df[,-c(1:2)], 2, mean, na.rm=TRUE)
}

Ztest <- function(value, randomlist){
  value <- as.numeric(value)
  randomlist <- as.numeric(randomlist)
  value[is.na(value)] <- 0
  randomlist <- randomlist[is.na(randomlist)==FALSE]
  
  #
  value <- abs(value)
  randomlist <- abs(randomlist)
  
  zv <- (value - mean(randomlist, na.rm=TRUE))/sd(randomlist)
  Pval <- pnorm(zv, 0,1, lower.tail = F)         ## is zv larger than Random values? 
  #
  Zresults <- c(round(zv,3), Pval, round(value,3))
  
  return(Zresults)
}

Ztes_for_list <- function(tf){
  ## Read random wPCC means
  rdb <- fread(paste0("GAN_R_CoExpData/MeanVals/wPCCm_GAN_Random_",tf, ".txt"))
  # rdb <- R_GAN_wPCC[[tf]]
  
  ## make a list by CoExp net (rows)
  #rdb.list <- MakelistByRows(rdb)
  rdb.list <- split(rdb$wPCCm, rdb$Net)
  rdb.list$n13a
  nets <- names(rdb.list)
  
  # obs val
  obs <- GAN_wPCCmean[[tf]][nets]
  
  #print(obs)
  out <- as.data.frame(t(mapply(Ztest, obs[nets],  rdb.list[nets])), stringAsFactor=F)
  colnames(out) <- c("Z", "Pval", "wPCCm")
  out[,"Nets"] <- nets
  out[,"FDR"] <- p.adjust(out$Pval, method = 'fdr')
  out <- as_tibble(out)
  
  #
  return(out)
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

CornCYC_size <- as_tibble(as.data.frame(table(unique(CornCYC[,c(1,3)])$Pathway.id), stringsAsFactors = F))
colnames(CornCYC_size) <- c("PWY", "PWYSize")
##################################################

###################################################################################################
###################                     Networks                         ##########################
###################################################################################################

# PDI
PDI <- as_tibble(unique(fread("../Fig_PDI/PDI_NetworkFinal.10_11_2021.txt")[,c(2,3)]))
colnames(PDI)[1] <- "Source"
#write.table(PDI, "../Fig_pecanpy/uniqGRN.10_11_2021.txt", row.names = F, sep = "\t", quote = F, col.names = F)

# CoExp
CoExp <- as_tibble(unique(fread("../Fig_Coexpression/CoExp_NetworkFinal.10_11_2021.txt")))
colnames(CoExp)[2] <- "Source"
CoExp <- unique(CoExp[,2:3])
#write.table(CoExp, "../Fig_pecanpy/uniqCEN.10_11_2021.txt", row.names = F, sep = "\t", quote = F,  col.names = F)

# teQTL
teQTL <- as_tibble(unique(fread("../Fig_transeQTL/teQTL_NetworkFinal.10_11_2021.txt")))
colnames(teQTL)[1] <- "Source"
#teQTL["val"] <- 1
#write.table(teQTL, "../Fig_pecanpy/uniqGAN.10_11_2021.txt", row.names = F, sep = "\t", quote = F,  col.names = F)
#teQTL <- teQTL[,1:2]
teQTLtf <- subset(teQTL, Source %in% unique(c(TF_CoR$GeneID, PDI$Source, CoExp$Source))) 


All_TFs <- unique(c(PDI$Source, CoExp$Source, TF_CoR$GeneID))

# Total TFs by layer
Total_TFs_list = list(CEN=CoExp$Source,  GRN=PDI$Source, GAN=teQTLtf$Source)
Total_TFs_list <- lapply(Total_TFs_list, unique)

# TF in venn groups
ven_file <- venn(Total_TFs_list)
ven_file <- as.list(attr(ven_file, "intersections"))
lapply(ven_file, length)



###################################################################################################
###################                   PWYs results                         ########################
###################################################################################################

#######################################################
## TF_PWY enrichment results from individual layers  ##
#######################################################
## PDI
PDI_PWY <- as_tibble(read.table("../Fig_PDI/PDI_NetworkFinal.CornCYC.04_18_2021.txt", h=T, stringsAsFactors = F))
PDI_PWY$TF <- ReplaceName(PDI_PWY$TF)
PDI_PWY[,"Net"] <- "PDI"

length(unique(PDI_PWY$TFid))
length(unique(PDI_PWY$PWY))

## CoExp
CoExp_PWY <- as_tibble(read.table("../Fig_Coexpression/CoExp_NetworkFinal.Full.CornCYC.04_18_2022.txt", h=T, stringsAsFactors = F))
CoExp_PWY[,"TFid"] <- CoExp_PWY$TF
CoExp_PWY$TF <- ReplaceName(CoExp_PWY$TF)
colnames(CoExp_PWY)[3] <- 'PWY'
CoExp_PWY[,"Net"] <- "CoExp"

## teQTL
teQTL_PWY <- as_tibble(read.table("../Fig_transeQTL/teQTL_NetworkFinal.CornCYC.04_18_2021.txt", h=T, stringsAsFactors = F))
teQTL_PWY[,"TFid"] <- teQTL_PWY$TF
teQTL_PWY <- subset(teQTL_PWY, TFid %in% All_TFs)
teQTL_PWY$TF <- ReplaceName(teQTL_PWY$TF)
teQTL_PWY[,"Net"] <- "teQTL"


#######################################################



#######################################################
##  TF_PWY enrichment results from common layers     ##
#######################################################
CommTarg_3TF_PWY <- as_tibble(fread("../Fig_CommonTarg/CommTarg_Network_3TF.CornCYC.04.21.2022.txt", h=T, stringsAsFactors = F))
CommTarg_3TF_PWY[,"TFid"] <- CommTarg_3TF_PWY$TF
CommTarg_3TF_PWY$TF <- ReplaceName(CommTarg_3TF_PWY$TF)
CommTarg_3TF_PWY[,"Net"] <- paste0("3TF_", CommTarg_3TF_PWY$Class)

#######################################################


#######################################################
##         TF_PWY enrichment Clustering              ##
#######################################################

# enrichment on clusters-tsne
tsne_PWY <- as_tibble(read.table("../Fig_tsne/Modulestsne.CornCYC.04.21.2022.txt", h=T))
colnames(tsne_PWY)[6] <- "nPWY"

# enrichment on clusters-pecanpy
peca_PWY <- as_tibble(read.table("../Fig_pecanpy/ModulesPecapy.CornCYC.04_19_2021.txt", h=T))
colnames(peca_PWY)[6] <- "nPWY"

peca_PWY

subset(tsne_PWY, n.targ>1)
subset(peca_PWY, n.targ>1)

tsneModules <- as_tibble(read.table("../Fig_tsne/Cluster_Round2_all_tsne.txt", h=T))[,c(1,3)]
PecaModules <- as_tibble(read.table("../Fig_pecanpy/Cluster_Round2_all_Pecapy.txt", h=T))[,c(2,4)]

# Subset to TF in groups of total TFs 
TFs_tM <- subset(tsneModules, GeneID %in% All_TFs)
TFs_pM <- subset(PecaModules, GeneID %in% All_TFs)


# Subset to TF in groups of 112 TFs 
TF112_tM <- subset(tsneModules, GeneID %in% ven_file$`CEN:GRN:GAN`)
TF112_pM <- subset(PecaModules, GeneID %in% ven_file$`CEN:GRN:GAN`)

# map  PWY enriched to TF by module 
TF112_tM <- left_join(TF112_tM, tsne_PWY, by="Module")
TF112_pM <- left_join(TF112_pM, peca_PWY, by="Module")

TFs_tM  <- left_join(TFs_tM, tsne_PWY, by="Module")
TFs_pM  <- left_join(TFs_pM, peca_PWY, by="Module")

TFs_tM <- subset(TFs_tM, n.targ>2)
TFs_pM <- subset(TFs_pM, n.targ>2)

length(unique(TFs_tM$GeneID))
length(unique(TFs_pM$GeneID))

TF112_tM <- TF112_tM[!(is.na(TF112_tM$PWY)),]
TF112_pM <- TF112_pM[!(is.na(TF112_pM$PWY)),]

TFs_tM <- TFs_tM[!(is.na(TFs_tM$PWY)),]
TFs_pM <- TFs_pM[!(is.na(TFs_pM$PWY)),]

as.data.frame(table(TFs_tM$GeneID))
TFs_tM

TF112_tM[,"Net"] <- "tsne"
TF112_pM[,"Net"] <- "Peca"

TF112_tM[,"TFid"] <- TF112_tM$GeneID
TF112_tM$GeneID <- ReplaceName(TF112_tM$GeneID)

TF112_pM[,"TFid"] <- TF112_pM$GeneID
TF112_pM$GeneID <- ReplaceName(TF112_pM$GeneID)

length(unique(TF112_tM$TF))
length(unique(TF112_pM$TF))

colnames(TF112_tM)[c(1)] <- "TF"
colnames(TF112_pM)[c(1)] <- "TF"

subset(TF112_tM, GeneID == "Zm00001d005737") 
subset(TF112_pM, GeneID == "Zm00001d005737")
subset(CoExp_PWY, TFid == "Zm00001d005737")

#######################################################

#######################################################

ColSelcected <- c("TF","PWY", "TFid", "Net")
Red_CoExp_PWY <- unique(CoExp_PWY[, ColSelcected])
Red_PDI_PWY <- unique(PDI_PWY[,ColSelcected])
Red_teQTL_PWY <- unique(teQTL_PWY[, ColSelcected])

Red_CT_3TF_PWY <- unique(CommTarg_3TF_PWY[, ColSelcected])
Red_TF112_tM <- unique(TF112_tM[, ColSelcected])
Red_TF112_pM <- unique(TF112_pM[, ColSelcected])

CEN_GRN_GAN_TF_PWY <- unique(rbind(Red_CoExp_PWY,
                                   Red_PDI_PWY,
                                   Red_teQTL_PWY))

CEN_GRN_GAN_TF_PWY[,"Index"] <- paste0(CEN_GRN_GAN_TF_PWY$TFid, "_", CEN_GRN_GAN_TF_PWY$PWY)


CEN_GRN_GAN_TF_PWY_Freq <-  as_tibble(as.data.frame(table(CEN_GRN_GAN_TF_PWY$Index)))

Top_CEN_GRN_GAN_TF_PWY <- subset(CEN_GRN_GAN_TF_PWY, Index %in% subset(CEN_GRN_GAN_TF_PWY_Freq, Freq==3)$Var1)

Top_CEN_GRN_GAN_TF_PWY <- left_join(Top_CEN_GRN_GAN_TF_PWY, unique(CornCYC[,c(1,2)]), by=c("PWY"="Pathway.id"))
View(Top_CEN_GRN_GAN_TF_PWY)


## Total list of interactions
Total_TF_PWY <- unique(rbind(Red_CT_3TF_PWY, 
                             Red_CoExp_PWY,
                             Red_PDI_PWY,
                             Red_teQTL_PWY,
                             Red_TF112_tM,
                             Red_TF112_pM))


Total_TF_PWY <- left_join(Total_TF_PWY, unique(CornCYC[,c(1,2)]), by=c("PWY"="Pathway.id"))

Total_TF_PWY[, "Index"] <- paste0(Total_TF_PWY$TFid, "_", Total_TF_PWY$PWY)

Total_TF_PWY.list <- split(Total_TF_PWY$Index, Total_TF_PWY$Net)

CENGRNGAN_list <- split(CEN_GRN_GAN_TF_PWY$Index, CEN_GRN_GAN_TF_PWY$Net)

# TF in venn groups
ven_CENGRNGAN <- venn(CENGRNGAN_list)
ven_CENGRNGAN <- as.list(attr(ven_CENGRNGAN, "intersections"))
lapply(ven_CENGRNGAN, length)
names(ven_file)

# get p-value bases on sampling
PossibleComparisons <- round((length(unique(Total_TF_PWY$TF)) * length(unique(Total_TF_PWY$PWY)))/2)
All_PWYs_go.obj <- newGOM(Total_TF_PWY.list, Total_TF_PWY.list, genome.size=PossibleComparisons) # all vs all 
drawHeatmap(All_PWYs_go.obj)

Pval <- getMatrix(All_PWYs_go.obj, name="pval")
common <- getMatrix(All_PWYs_go.obj, name="intersection")

diag(Pval) <- 1
Pval_log <- -log10(Pval)
max(Pval_log != Inf, 2)

Pval_log[Pval_log == Inf ] <- max(Pval_log[Pval_log != Inf])

Pval_log[Pval_log >= 5 ] <- 5

hm1 <- Heatmap(Pval_log[colnames(Pval_log),],
        #heatmap_width = unit(15, "cm"), 
        #heatmap_height = unit(8, "cm"), 
        #column_names_rot = 90, 
        name = "-log10 Pval", 
        cluster_columns = TRUE, column_dend_reorder = TRUE,
        show_row_dend = FALSE, show_column_dend = FALSE,
        col=viridis(10, direction = 1),
        column_names_gp = gpar(fontsize =12),
        row_names_gp = gpar(fontsize = 12))

hm2 <- Heatmap(common,
               #heatmap_width = unit(15, "cm"), 
               #heatmap_height = unit(8, "cm"), 
               #column_names_rot = 90, 
               name = "-log10 Pval", 
               cluster_columns = TRUE, column_dend_reorder = TRUE,
               show_row_dend = FALSE, show_column_dend = FALSE,
               col=viridis(100, direction = 1),
               column_names_gp = gpar(fontsize =12),
               row_names_gp = gpar(fontsize = 12))

hm2

#######################################################



###################################################################################################

###################################################################################################
#########                                 Summary and plots                               #########
###################################################################################################

#############################################################################################
#######################            Common prediction          ###############################
#############################################################################################

GenesMaize <- as.character(read.table("Data/Annotations/Zm.v4.synteny.genes.txt", h=T)$gene_id)




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
########   Summary/general number of PWY enrichment results   ######## 
###################################################################### 

###############################################################
########   Count number of TFs by PWY 
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

# Venn diagram of common PWTs 
ResPWY_list = list(CEN=unique(CoExp_PWY$PWY), 
                   GRN=unique(PDI_PWY$PWY),
                   GAN=unique(subset(teQTL_PWY, TFid %in% All_TFs)$PWY))


# Venn diagram of common TFs
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

wc
ReplaceNamePWY(subset(CoExp_PWY, TFid=="Zm00001d050195" & n.targ > 1 & Pval <= 0.01)$PWY)  # bZIP38
ReplaceNamePWY(subset(PDI_PWY,   TFid=="Zm00001d050195" & n.targ > 1 & Pval <= 0.01)$PWY)  # 
ReplaceNamePWY(subset(teQTL_PWY,   TFid=="Zm00001d050195" & n.targ >= 1 & Pval <= 0.05)$PWY)  # 


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
                       width = unit(4, "cm"),
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


#######################################################
##      Compare common TF_PWY predictions            ##
#######################################################

##########
# GRN, CEN, and GAN comparisons
GRN_index <- paste0(PDI_PWY$TFid, "_",PDI_PWY$PWY)
CEN_index <- paste0(CoExp_PWY$TFid, "_",CoExp_PWY$PWY)
GAN_index <- paste0(teQTL_PWY$TFid, "_",teQTL_PWY$PWY)

TF_PWY_list <- list(GRN=GRN_index, CEN=CEN_index, GAN=GAN_index)

TFs_list <- list(GRN=unique(PDI_PWY$TFid), 
                 CEN=unique(CoExp_PWY$TFid),
                 GAN=unique(teQTL_PWY$TFid))

PWYs_list <- list(GRN=unique(PDI_PWY$PWY), 
                  CEN=unique(CoExp_PWY$PWY),
                  GAN=unique(teQTL_PWY$PWY))

# get p-value bases on sampling
Plot_Common_Classic_predictions <- ggarrange(vennfunc(TF_PWY_list), 
                                             vennfunc(PWYs_list),
                                             vennfunc(TFs_list), ncol = 3, 
                                             widths = c(1.5, 1.5, 1.5))


Plot_Common_Classic_predictions

tiff("Plots/Plot_Common_Classic_predictions.tiff", units="in", width=6, height=2.5, res=300)
print(Plot_Common_Classic_predictions)
dev.off()


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
    geom_sf(aes(color=name), size = 1, data = venn_setedge(data), show.legend = F) +
    #
    geom_sf_text(aes(label = name), size=3, data = venn_setlabel(data)) +
    geom_sf_text(aes(label= scales::comma(count, accuracy = 1)), size=3, data = venn_region(data)) +
    #
    scale_fill_manual(values = col2) + # 
    scale_color_manual(values = alpha(col, .5)) +
    #
    theme_void() +
    theme(plot.margin = unit(c(0.5,1,1,0.1), "cm")) #+
  #xlim(-150,1000)
}

TF_PWY_venn_obj <- venn(TF_PWY_list, show.plot = F, intersections = T)

TF_PWY_venn_obj
TF_PWY_venn_obj <- as.list(attr(TF_PWY_venn_obj, "intersections"))
TF_PWY_venn_obj$`GRN:CEN:GAN`

# Refine

Common_TF_PWY_GRN <- MakeHeatmapVertical(PDI_PWY[GRN_index %in% TF_PWY_venn_obj$`GRN:CEN:GAN`,], 4)
Common_TF_PWY_CEN <- MakeHeatmapVertical(CoExp_PWY[CEN_index %in% TF_PWY_venn_obj$`GRN:CEN:GAN`,], 4)
Common_TF_PWY_GAN <- MakeHeatmapVertical(teQTL_PWY[GAN_index %in% TF_PWY_venn_obj$`GRN:CEN:GAN`,], 4)

mask <- CEN_index %in% TF_PWY_venn_obj$`GRN:CEN:GAN`

CoExp_PWY[mask,]

dim(CoExp_PWY)

list_hm2 = Common_TF_PWY_GRN %v% Common_TF_PWY_CEN %v% Common_TF_PWY_GAN
list_hm2

# size 8x8
tiff("Plots/Plot_Heatmap_Common_TF_PWYs.tiff", units="in", width=4, height=8, res=300)
print(draw(list_hm2, heatmap_legend_side = "bottom"))
dev.off()

######################################################################################
################            Count common Targets amount layers 
######################################################################################

# list of targets by TF
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


Venn_DF <- as_tibble(rbind(Venn_CEN_GRN, Venn_GRN_GAN, Venn_CEN_GAN, Venn_CEN_GAN_GRN))
# to fix scale maximum 
Venn_DF$V2[(Venn_DF$V2 <= min(Venn_CEN_GAN$V2))] <- min(Venn_CEN_GAN$V2) 

Venn_DF_CEN_GAN_GRNpairs <- as_tibble(rbind(Venn_CEN_GAN_GRNcen_grn, 
                                            Venn_CEN_GAN_GRNcen_gan, 
                                            Venn_CEN_GAN_GRNgan_grn))
# to fix scale maximum 
Venn_DF_CEN_GAN_GRNpairs$V2[(Venn_DF_CEN_GAN_GRNpairs$V2 <= 1e-20)] <- 1e-20


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

######################################################################################

###############################################################
##        Common TFs based on common PWY's analysis          ##
###############################################################

CommonTFs_in_Top_PWYs <- as.data.frame(table(c(unique(CoExpDFTop$TFid), unique(PDIDFTop$TFid), unique(teQTLDFTop$TFid))))
CommonTFs_in_Top_PWYs <- subset(CommonTFs_in_Top_PWYs, Freq>1)
ReplaceName(CommonTFs_in_Top_PWYs$Var1)


CommonTFs_in_Top5_pwy <- unique(subset(CoExp_PWY, PWY %in% subset(PWY_Freq, Freq==3)$PWY)$TFid[(subset(CoExp_PWY, PWY %in% subset(PWY_Freq, Freq==3)$PWY)$TFid %in% subset(PDI_PWY, PWY %in% subset(PWY_Freq, Freq==3)$PWY)$TFid)])
ReplaceName(CommonTFs_in_Top5_pwy)
ReplaceName("Zm00001d020492")

CoExpDFTop[CoExpDFTop$TFid %in% CommonTFs_in_Top5_pwy,]
PDIDFTop[PDIDFTop$TFid %in% CommonTFs_in_Top5_pwy,]


GetTF_PWY_Net <- function(tf, pwy) {
  # Annotate PWY genes
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

###############################################################

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

