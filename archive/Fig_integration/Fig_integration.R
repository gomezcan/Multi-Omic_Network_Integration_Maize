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

#
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
                           grid.text(sprintf("%.f", t(Text)[i, j]), x, y, gp = gpar(fontsize = 12, col = "gray80"))
                         }
                         
                       },
                       name = "Percentage PWY", 
                       cluster_columns = TRUE, column_dend_reorder = TRUE,
                       show_row_dend = FALSE, show_column_dend = FALSE,
                       col=cor_scale,
                       column_names_gp = gpar(fontsize = 12),
                       row_names_gp = gpar(fontsize = 12),
                       show_heatmap_legend = T,
                       heatmap_legend_param = list(title_gp = gpar(fontsize = 10), 
                                                   labels_gp = gpar(fontsize = 10),
                                                   width = unit(15, "mm"),
                                                   direction = "horizontal")) 
  return(hm_common)
  
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

CornCYC_size <- as_tibble(as.data.frame(table(unique(CornCYC[,c(1,3)])$Pathway.id), stringsAsFactors = F))
colnames(CornCYC_size) <- c("PWY", "PWYSize")

###### TF_PWY enrichment results  ######
## PDI
PDI_PWY <- as_tibble(read.table("PDI_NetworkFinal.CornCYC.10_11_2021.txt", h=T, stringsAsFactors = F))
PDI_PWY$TF <- ReplaceName(PDI_PWY$TF)
PDI_PWY[,"Net"] <- "PDI"

## CoExp
CoExp_PWY <- as_tibble(read.table("CoExp_NetworkFinal.Full.CornCYC.10_11_2021.txt", h=T, stringsAsFactors = F))
CoExp_PWY[,"TFid"] <- CoExp_PWY$TF
CoExp_PWY$TF <- ReplaceName(CoExp_PWY$TF)
colnames(CoExp_PWY)[3] <- 'PWY'
CoExp_PWY[,"Net"] <- "CoExp"

## teQTL
teQTL_PWY <- as_tibble(read.table("teQTL_NetworkFinal.CornCYC.10_11_2021.txt", h=T, stringsAsFactors = F))
teQTL_PWY[,"TFid"] <- teQTL_PWY$TF
teQTL_PWY$TF <- ReplaceName(teQTL_PWY$TF)
teQTL_PWY[,"Net"] <- "teQTL"

#
ColSelcected <- c("TF","PWY","padj", "n.targ", "targets", "TFid", "Net")

#
CoExp_PWY <- unique(CoExp_PWY[, colnames(CoExp_PWY) %in% ColSelcected])
PDI_PWY <- unique(PDI_PWY[, colnames(PDI_PWY) %in% ColSelcected])
teQTL_PWY <- unique(teQTL_PWY[, colnames(teQTL_PWY) %in% ColSelcected])

## Total list of interactions
Total_TF_PWY <- unique(rbind(CoExp_PWY[, c("TFid", "PWY")],
                      PDI_PWY[, c("TFid", "PWY")],
                      teQTL_PWY[, c("TFid", "PWY")]))

unique(CoExp_PWY[,c("TFid", "PWY")])
length(unique(CoExp_PWY$TFid))
length(unique(CoExp_PWY$PWY))

unique(PDI_PWY[,c("TFid", "PWY")])
length(unique(PDI_PWY$TFid))
length(unique(PDI_PWY$PWY))

unique(teQTL_PWY[,c("TFid", "PWY")])
length(unique(teQTL_PWY$TFid))
length(unique(teQTL_PWY$PWY))


length(unique(Total_TF_PWY$TFid))
length(unique(Total_TF_PWY$PWY))



######################################################################################
######### Summary and plots
######################################################################################
# Load library
library(VennDiagram)

# Chart
venn.diagram(
  x = list(unique(CoExp_PWY$TFid), unique(PDI_PWY$TFid), unique(teQTL_PWY$TFid)),
  category.names = c("CoExp" , "PDI" , "teQTL"),
  filename = 'TF_venn_diagramm.png',
  output=TRUE
)

venn.diagram(
  x = list(unique(CoExp_PWY$PWY), unique(PDI_PWY$PWY), unique(teQTL_PWY$PWY)),
  category.names = c("CoExp" , "PDI" , "teQTL"),
  filename = 'PWY_venn_diagramm.png',
  output=TRUE
)

## count PWY enrichment frequency among all three networks  ##
PWY_Freq <- as_tibble(as.data.frame(table(c(unique(CoExp_PWY$PWY), unique(PDI_PWY$PWY), unique(teQTL_PWY$PWY))), stringsAsFactors = F))
colnames(PWY_Freq)[1] <- c("PWY")
PWY_Freq[,"PWYname"] <- ReplaceNamePWY(PWY_Freq$PWY)



# TFs in PWYs highly predicted
TFsinPDI_toppwy   <- unique(subset(PDI_PWY, PWY %in% subset(PWY_Freq, Freq==3)$PWY)$TFid)
TFsinCoExp_toppwy <- unique(subset(CoExp_PWY, PWY %in% subset(PWY_Freq, Freq==3)$PWY)$TFid)
#
CoExpDFTop <- subset(CoExp_PWY, PWY %in% subset(PWY_Freq, Freq==3)$PWY)
PDIDFTop   <- subset(PDI_PWY, PWY %in% subset(PWY_Freq, Freq==3)$PWY)
teQTLDFTop <- subset(teQTL_PWY, PWY %in% subset(PWY_Freq, Freq==3)$PWY)

CoExp_heatmap <- MakeHeatmap(CoExpDFTop)
PDI_heatmap <- MakeHeatmap(PDIDFTop)
teQTL_heatmap <- MakeHeatmap(teQTLDFTop)

#
CommonTFs_in_Top5_pwy <- unique(subset(CoExp_PWY, PWY %in% subset(PWY_Freq, Freq==3)$PWY)$TFid[(subset(CoExp_PWY, PWY %in% subset(PWY_Freq, Freq==3)$PWY)$TFid %in% subset(PDI_PWY, PWY %in% subset(PWY_Freq, Freq==3)$PWY)$TFid)])
ReplaceName(CommonTFs_in_Top5_pwy)
list_hm = CoExp_heatmap + PDI_heatmap + teQTL_heatmap
# size 4x20
draw(list_hm, heatmap_legend_side = "bottom")  

##################################################
########       heatmap based on TFs       ########
##################################################

# 1. Identification of Common TFs between  
CoExp_PWY_Common_TFs <- CoExp_PWY[CoExp_PWY$TFid %in% PDI_PWY$TFid, ]
PDI_PWY_Common_TFs <- PDI_PWY[PDI_PWY$TFid %in% CoExp_PWY$TFid, ]
#
CoExp_PWY_Common_TFs[,"PWYname"] <- ReplaceNamePWY(CoExp_PWY_Common_TFs$PWY)
PDI_PWY_Common_TFs[,"PWYname"] <- ReplaceNamePWY(PDI_PWY_Common_TFs$PWY)

write.table(CoExp_PWY_Common_TFs, 'CoExp_Commom.39TFs.txt', sep = '\t', quote = F, row.names = F)
write.table(PDI_PWY_Common_TFs, 'PDI_Commom.39TFs.txt', sep = '\t', quote = F, row.names = F)

# 2. heatmap with common TFs in top 5 PWYs
CoExp_heatmapTopTFs <- MakeHeatmap(subset(CoExpDFTop, TFid %in% CommonTFs_in_Top5_pwy))
PDI_heatmapTopTFs <- MakeHeatmap(subset(PDIDFTop, TFid %in% CommonTFs_in_Top5_pwy))

draw(CoExp_heatmapTopTFs+PDI_heatmapTopTFs, heatmap_legend_side = "bottom")  


Top5_pwy_Candiates <- rbind(tibble(layer="PDI",TFid=subset(PDI_PWY, PWY %in% subset(PWY_Freq, Freq==3)$PWY)$TFid),
                   tibble(layer="teQTL",TFid=subset(teQTL_PWY, PWY %in% subset(PWY_Freq, Freq==3)$PWY)$TFid),
                   tibble(layer="CoExp",TFid=subset(CoExp_PWY, PWY %in% subset(PWY_Freq, Freq==3)$PWY)$TFid))


write.table(Top5_pwy_Candiates, "Top5_pwy_Candiates.txt", row.names = F, quote = F, sep = "\t")

subset(CoExp_PWY, PWY %in% subset(PWY_Freq, Freq==3)$PWY)$TFid

TFs_by_ClassesM <- reshape2::dcast(TFs_by_Classes_summary, TF ~ class)
row.names(TFs_by_ClassesM) <- TFs_by_ClassesM$TF
TFs_by_ClassesM <- TFs_by_ClassesM[,-c(1)]

################################################
#### Heatmap TF vs pathway
################################################


Fullnet.table <- gather(Fullnet, key, value, -TF, -Target, -GeneName, -Source)
colnames(Fullnet.table)[6] <- 'Enriched'

Fullnet.table$Enriched[Fullnet.table$Enriched==0] <- "FDR > 0.1"
Fullnet.table$Enriched[Fullnet.table$Enriched==1] <- "FDR < 0.1"

length(unique(subset(Fullnet.table, TF=="MYB31" & Source >1)$Target))
length(unique(subset(Fullnet.table, TF=="BHLH91" & Source >1)$Target))
length(unique(subset(Fullnet.table, TF=="BHLH136" & Source >1)$Target))

Plot.MYB31 <- ggplot(subset(Fullnet.table, TF=="MYB31" & Source >1),
                     aes(y=reorder(GeneName, Source),  x=key, fill=Enriched))+
  geom_tile()  + #scale_x_continuous(breaks = seq(0,100,20), expand = c(0,0)) + 
  #scale_fill_viridis(option = "D", direction = -1, discrete = T) + #
  scale_fill_manual(values=c('darkorchid1', 'gainsboro')) +
  theme_bw() +
  theme(legend.title = element_text(vjust = 0.8), 
        panel.border = element_blank(), 
        panel.grid.major = element_blank())


Plot.bHLH91 <- ggplot(subset(Fullnet.table, TF=="BHLH91" & Source >1),
                     aes(y=reorder(GeneName, Source),  x=key, fill=Enriched))+
  geom_tile()  + #scale_x_continuous(breaks = seq(0,100,20), expand = c(0,0)) + 
  #scale_fill_viridis(option = "D", direction = -1, discrete = T) + #
  scale_fill_manual(values=c('darkorchid1', 'gainsboro')) +
  theme_bw() +
  theme(legend.title = element_text(vjust = 0.8), 
        panel.border = element_blank(), 
        panel.grid.major = element_blank())

Plot.bHLH136 <- ggplot(subset(Fullnet.table, TF=="BHLH136" & Source >1),
                      aes(y=reorder(GeneName, Source),  x=key, fill=Enriched))+
  geom_tile()  + #scale_x_continuous(breaks = seq(0,100,20), expand = c(0,0)) + 
  #scale_fill_viridis(option = "D", direction = -1, discrete = T) + #
  scale_fill_manual(values=c('darkorchid1', 'gainsboro')) +
  theme_bw() +
  theme(legend.title = element_text(vjust = 0.8), 
        panel.border = element_blank(), 
        panel.grid.major = element_blank())


Plot.MYB31 <- ggpar(Plot.MYB31, ylab = "Target", xlab = "Prediction method", font.ytickslab = 5, x.text.angle = 45)
Plot.bHLH91 <- ggpar(Plot.bHLH91, ylab = "Target", xlab = "Prediction method", font.ytickslab = 5, x.text.angle = 45)
Plot.bHLH136 <- ggpar(Plot.bHLH136, ylab = "Target", xlab = "Prediction method", font.ytickslab = 5, x.text.angle = 45)

ggarrange(Plot.MYB31, Plot.bHLH91, Plot.bHLH136,
          common.legend = T, ncol = 3, legend = 'bottom')

###### description of cis
library(reshape2)
Fullnet.table[,"Index"] <- paste(Fullnet.table$TF, Fullnet.table$Target, sep = "_")

top.top.cis <- subset(Fullnet.table, Source >1 & key =='ceQTL')$Index
top.top.top.cis <- subset(Fullnet.table, Source > 2 & key =='ceQTL')$Index


Top.ceQTL[, "Index"] <- paste(Top.ceQTL$TF.Name, Top.ceQTL$target, sep = "_")

Top.ceQTL.source <- subset(Top.ceQTL,  Index %in% top.top.cis)[,c(2,4,1)]
Top.Top.ceQTL.source <- subset(Top.ceQTL,  Index %in% top.top.top.cis)[,c(2,4,1)]

Top.ceQTL.source  <- left_join(Top.ceQTL.source, PheGenes[,1:2], by=c("target"="GeneID"))
Top.Top.ceQTL.source <- left_join(Top.Top.ceQTL.source, PheGenes[,1:2], by=c("target"="GeneID"))

Top.ceQTL.source[,'Mask'] <- 'Yes'
Top.Top.ceQTL.source[,'Mask'] <- 'Yes'

unique(subset(Top.Top.ceQTL.source, TF.Name=="MYB31")$GeneName)


Plot.top.top.ceqtl.MYB31 <- ggplot(subset(Top.Top.ceQTL.source, TF.Name=="MYB31"),
                       aes(y=GeneName,  x=source, fill=Mask))+
  geom_tile()  + #scale_x_continuous(breaks = seq(0,100,20), expand = c(0,0)) + 
  #scale_fill_viridis(option = "D", direction = -1, discrete = T) + #
  scale_fill_manual(values=c('darkorange', 'gainsboro')) +
  theme_bw() +
  theme(legend.title = element_text(vjust = 0.8), 
        panel.border = element_blank(), 
        panel.grid.major = element_blank())

Plot.top.top.ceqtl.MYB31 <- ggpar(Plot.top.top.ceqtl.MYB31,
                                  ylab = "Target", xlab = "eQTL source", 
                                  legend = 'None')

Plot.top.top.ceqtl.bHLH91 <- ggplot(subset(Top.Top.ceQTL.source, TF.Name=="BHLH91"),
                                   aes(y=GeneName,  x=source, fill=Mask))+
  geom_tile()  + #scale_x_continuous(breaks = seq(0,100,20), expand = c(0,0)) + 
  #scale_fill_viridis(option = "D", direction = -1, discrete = T) + #
  scale_fill_manual(values=c('darkorange', 'gainsboro')) +
  theme_bw() +
  theme(legend.title = element_text(vjust = 0.8), 
        panel.border = element_blank(), 
        panel.grid.major = element_blank())

Plot.top.top.ceqtl.bHLH91 <- ggpar(Plot.top.top.ceqtl.bHLH91,
                                  ylab = "Target", xlab = "eQTL source", 
                                  legend = 'None')
Plot.top.top.ceqtl.bHLH91


Plot.top.top.ceqtl.bHLH136 <- ggplot(subset(Top.Top.ceQTL.source, TF.Name=="BHLH136"),
                                    aes(y=GeneName,  x=source, fill=Mask))+
  geom_tile()  + #scale_x_continuous(breaks = seq(0,100,20), expand = c(0,0)) + 
  #scale_fill_viridis(option = "D", direction = -1, discrete = T) + #
  scale_fill_manual(values=c('darkorange', 'gainsboro')) +
  theme_bw() +
  theme(legend.title = element_text(vjust = 0.8), 
        panel.border = element_blank(), 
        panel.grid.major = element_blank())

Plot.top.top.ceqtl.bHLH136 <- ggpar(Plot.top.top.ceqtl.bHLH136,
                                   ylab = "Target", xlab = "eQTL source", 
                                   legend = 'None')
Plot.top.top.ceqtl.bHLH136



###### description of CoExp
Top.CoExp[, "Index"] <- paste(ReplaceName(Top.CoExp$reg.gid), Top.CoExp$tgt.gid, sep = "_")
top.top.CoExp <- subset(Fullnet.table, Source >1 & key =='CoExp')$Index
top.top.top.CoExp <- subset(Fullnet.table, Source >2 & key =='CoExp')$Index


Top.CoExp.source <- subset(Top.CoExp,  Index %in% top.top.CoExp) [,c(2,3,5,1)]
Top.Top.CoExp.source <- subset(Top.CoExp,  Index %in% top.top.top.CoExp) [,c(2,3,5,1)]

Top.CoExp.source    <- left_join(Top.CoExp.source, PheGenes[,1:2], by=c("tgt.gid"="GeneID"))
Top.Top.CoExp.source <- left_join(Top.Top.CoExp.source, PheGenes[,1:2], by=c("tgt.gid"="GeneID"))

Top.CoExp.source$reg.gid <- ReplaceName(Top.CoExp.source$reg.gid)
Top.Top.CoExp.source$reg.gid <- ReplaceName(Top.Top.CoExp.source$reg.gid)

Top.CoExp.source <- left_join(Top.CoExp.source, x, by =c('.id'='nid'))
Top.Top.CoExp.source <- left_join(Top.Top.CoExp.source, x, by =c('.id'='nid'))

#
Plot.top.top.CoExp.MYB31 <- ggplot(subset(Top.Top.CoExp.source, reg.gid=="MYB31"),
                                   aes(y=GeneName,  x=nname, fill=pcc))+
  geom_tile()  + #scale_x_continuous(breaks = seq(0,100,20), expand = c(0,0)) + 
  scale_fill_viridis(option = "D", direction = -1) + #
  #scale_fill_manual(values=c('darkorange', 'gainsboro')) +
  theme_bw() +
  theme(legend.title = element_text(vjust = 0.8), 
        panel.border = element_blank(), 
        panel.grid.major = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

Plot.top.top.CoExp.MYB31 <- ggpar(Plot.top.top.CoExp.MYB31,
                                  ylab = "Target", xlab = "CoExp source", 
                                  legend = 'bottom')

Plot.top.top.CoExp.MYB31

Plot.top.top.CoExp.bHLH91 <- ggplot(subset(Top.Top.CoExp.source, reg.gid=="BHLH91"),
                                   aes(y=GeneName,  x=nname, fill=pcc))+
  geom_tile()  + #scale_x_continuous(breaks = seq(0,100,20), expand = c(0,0)) + 
  scale_fill_viridis(option = "D", direction = -1) + #
  #scale_fill_manual(values=c('darkorange', 'gainsboro')) +
  theme_bw() +
  theme(legend.title = element_text(vjust = 0.8), 
        panel.border = element_blank(), 
        panel.grid.major = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

Plot.top.top.CoExp.bHLH91 <- ggpar(Plot.top.top.CoExp.bHLH91,
                                  ylab = "Target", xlab = "CoExp source", 
                                  legend = 'bottom')
Plot.top.top.CoExp.bHLH91

Plot.top.top.CoExp.bHLH136 <- ggplot(subset(Top.Top.CoExp.source, reg.gid=="BHLH136"),
                                    aes(y=GeneName,  x=nname, fill=pcc))+
  geom_tile()  + #scale_x_continuous(breaks = seq(0,100,20), expand = c(0,0)) + 
  scale_fill_viridis(option = "D", direction = -1) + #
  #scale_fill_manual(values=c('darkorange', 'gainsboro')) +
  theme_bw() +
  theme(legend.title = element_text(vjust = 0.8), 
        panel.border = element_blank(), 
        panel.grid.major = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

Plot.top.top.CoExp.bHLH136 <- ggpar(Plot.top.top.CoExp.bHLH136,
                                   ylab = "Target", xlab = "CoExp source", 
                                   legend = 'bottom')
Plot.top.top.CoExp.bHLH136

################################################################################################
## Integration og PDI candidates with co-expression




