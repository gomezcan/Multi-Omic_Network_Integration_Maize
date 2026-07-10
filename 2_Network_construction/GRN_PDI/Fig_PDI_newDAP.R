library(tidyverse)
library(ggrepel)
library(ggpubr)
library(data.table)
library(GeneOverlap)
library(viridis)
library(ComplexHeatmap)
library(fgsea)
library(reshape2)
library(scales)

##################################################
##########        Annotations       ##############
##################################################
# saf <- as_tibble(read.table("Data/eQTL_data/Zea_mays.B73_RefGen_v4.46.saf", stringsAsFactors = F))
# saf1 <- subset(saf, V5=="+")[,c(1,2,3)]
# saf2 <- subset(saf, V5=="-")[,c(1,2,4)]
# colnames(saf1) <- c("GeneID", "chrAnn", "TSS")
# colnames(saf2) <- c("GeneID", "chrAnn", "TSS")
#
# saf <- rbind(saf1, saf2)

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

#table(subset(PDI.class, TF %in% Y1H$TF.v4)$source)

## syntenic genes 
Syntenic <- as_tibble(read.table("Data/Annotations/Zm.v4.synteny.genes.txt", h=T, stringsAsFactors = F, sep = '\t'))
Syntenic <- subset(Syntenic, gene_synteny=='syntenic')$gene_id

# CoExp
CoExp <- as_tibble(unique(fread("../Fig_Coexpression/CoExp_NetworkFinal.10_11_2021.txt")))
colnames(CoExp)[2] <- "Source"
CoExp <- unique(CoExp[,2:3])

# Peak numbers
#Peak <- as_tibble(read.table("Table_S4.txt", h=T, stringsAsFactors = F))

# PDI network
PDInewDAP <- as_tibble(read.table("../Fig_PDI/NewDAPseq_All_Peals_02.05.2020.txt", h=T, stringsAsFactors = F))

PDI_Up5_Down02_kb <- subset(subset(PDInewDAP, Z > -0.5 & Dis <=0.5), Dis >= -3)
PDI_Up5_Down02_kb[,"TF"] <- gsub("newDAP.", "", PDI_Up5_Down02_kb$TFsample)
PDI_Up5_Down02_kb[,"TF"] <- sapply(strsplit(PDI_Up5_Down02_kb$TF, split='_', fixed=TRUE), `[`, 1) # add methods label
#
PDI_Up5_Down05_kb <- subset(subset(PDInewDAP, Z > -0.5 & Dis <=0.5), Dis >= -5)
PDI_Up5_Down05_kb[,"TF"] <- gsub("newDAP.", "", PDI_Up5_Down05_kb$TFsample)
PDI_Up5_Down05_kb[,"TF"] <- sapply(strsplit(PDI_Up5_Down05_kb$TF, split='_', fixed=TRUE), `[`, 1) # add methods label
PDI_Up5_Down05_kb[,"TFsample"] <- sapply(strsplit(PDI_Up5_Down05_kb$TFsample, split='_', fixed=TRUE), `[`, 2) # add methods label
#
PDI_Up5_Down05_kb <- PDI_Up5_Down05_kb[,c(9,1:8)]

TFsInNewDAP <- unique(PDI_Up5_Down05_kb$TF)

PDI_Up5_Down05_kb["CoExp"] <- "No"

TFsInNewDAP[1]
for (i in TFsInNewDAP) {
  tem_tar <-  unique(subset(CoExp, Source ==i)$Target)
  #tem_tar2 <- unique(subset(PDI_Up5_Down05_kb, TF ==i)$Target)
  
  PDI_Up5_Down05_kb[PDI_Up5_Down05_kb$TF == i,]$CoExp[PDI_Up5_Down05_kb[PDI_Up5_Down05_kb$TF == i,]$Target %in% tem_tar] <- "Yes"
  
  
  #[ (PDI_Up5_Down05_kb$Target %in% tem_tar)] <- "Yes"
}


write.table(PDI_Up5_Down05_kb[,c(9,1:5,7, 10)], "NewDAP_Up5_Down05_kb.txt", row.names = F, quote = F, sep = "\t")

################################################
# Counts common targets by Y1H and DAP-seq
################################################

TFs2test <- unique(PDI_Up5_Down02_kb$TF)

CountCommonY1H <- function(tf){
  
  tY1H <- unique(subset(Y1H, TF.v4==tf)$Target.v4)
  tDAP <- unique(subset(PDI_Up5_Down02_kb, TF==tf)$Target)
  
  common <- intersect(tY1H, tDAP)
  out <- tibble(Y1H=length(tY1H), DAP=length(tDAP), Common=length(common))
  return(out)
}

CountsY1H_DAP <-  lapply(TFs2test, CountCommonY1H)
names(CountsY1H_DAP) <- ReplaceName(TFs2test)
CountsY1H_DAP <- as_tibble(rbindlist(CountsY1H_DAP, idcol = T))

CountsY1H_DAP <- CountsY1H_DAP[CountsY1H_DAP$Common >0,]
CountsY1H_DAP <- gather(CountsY1H_DAP, key , value, -.id)
CountsY1H_DAP$.id <- as.character(CountsY1H_DAP$.id)
CountsY1H_DAP$key <- gsub("Common", "DAP_Y1H", CountsY1H_DAP$key)

Phe_Y1H_DAP_Plot <- ggplot(subset(CountsY1H_DAP, key!="DAP" & value>0), aes(y=.id, x=value, fill=key))+
  geom_bar(stat="identity", position=position_dodge()) +
  geom_text(aes(label=value, x= value+0.2 ), position = position_dodge(0.9))+
  ylab("TF") + xlab("Target")

Phe_Y1H_DAP_Plot <- ggpar(Phe_Y1H_DAP_Plot, font.tickslab=14, font.x = 14, font.y = 14)

###
# Files with targets and coords
###

CountCommonY1H_Net <- function(list_tf){
  
  fileout <- list()
  for (tf in list_tf){
    tY1H <- unique(subset(Y1H, TF.v4==tf)$Target.v4)
    tDAP <- unique(subset(PDI_Up5_Down02_kb, TF==tf)$Target)
    
    common <- intersect(tY1H, tDAP)
    out <- subset(PDI_Up5_Down02_kb, TF==tf & Target %in% common)
    fileout[[tf]] <- out
  }
  fileout <- rbindlist(fileout, idcol = F)
  return(as_tibble(fileout))
}

Y1H_DAP_Net <- CountCommonY1H_Net(TFs2test)
Y1H_DAP_Net["TFname"] <- ReplaceName(Y1H_DAP_Net$TF)
Y1H_DAP_Net <- Y1H_DAP_Net[,c("TF","TFname","TFsample", "Target", "Summit", "Dis", "OCR")]

write.table(Y1H_DAP_Net, "Y1H_DAP_Net.txt", sep = "\t", quote = F, row.names = F)


sum(subset(CountsY1H_DAP, key=="DAP_Y1H" & value>0)$value)

################################################

################################################
# Counts common targets by Y1H and DAP-seq
################################################


T2_CountCommonY1H <- function(tf){
  
  tY1H <- unique(Y1H$TF.v4)
  tDAP <- unique(subset(PDI_Up5_Down02_kb, TF==tf)$Target)
  
  # Common TFs in Y1H
  common1 <- intersect(tDAP, tY1H)
  # Common TFs in DAP
  common2 <- intersect(tDAP, TFs2test)
  
  out <- tibble(TF=ReplaceName(tf), tTFs_Y1H=length(common1), tTFs_DAP=length(common2))
  return(out)
}

TFs_CountsY1H_DAP <-  lapply(TFs2test, T2_CountCommonY1H)
TFs_CountsY1H_DAP <- as_tibble(rbindlist(TFs_CountsY1H_DAP, idcol = F))

TFs_CountsY1H_DAP <- subset(TFs_CountsY1H_DAP, tTFs_Y1H>0)

TFs_CountsY1H_DAP <- gather(TFs_CountsY1H_DAP, key , value, -TF)
TFs_CountsY1H_DAP$TF <- as.character(TFs_CountsY1H_DAP$TF)
TFs_CountsY1H_DAP$key <- gsub("tTFs_Y1H", "TFs_Y1H", TFs_CountsY1H_DAP$key)
TFs_CountsY1H_DAP$key <- gsub("tTFs_DAP", "TFs_DAP", TFs_CountsY1H_DAP$key)


TFs_Y1H_DAP_Plot <- ggplot(TFs_CountsY1H_DAP, aes(y=TF, x=value, fill=key))+
  geom_bar(stat="identity", position=position_dodge()) +
  geom_text(aes(label=value, x= value+0.2 ), position = position_dodge(0.9))+
  ylab("TF") + xlab("Target")

TFs_Y1H_DAP_Plot <- ggpar(TFs_Y1H_DAP_Plot,
                          font.tickslab=14, font.x = 14, font.y = 14)

###
# Files with targets and coords
#

Count_TFs_Y1H_DAP_Net <- function(list_tf){
  
  fileout <- list()
  for (tf in list_tf){
    
    tY1H <- unique(Y1H$TF.v4)
    tDAP <- unique(subset(PDI_Up5_Down02_kb, TF==tf)$Target)
    
    # Common TFs in Y1H
    common1 <- intersect(tDAP, tY1H)
    # Common TFs in DAP
    common2 <- intersect(tDAP, TFs2test)
    
    #TFs of interest
    common <- unique(c(common1, common2))
    
    out <- subset(PDI_Up5_Down02_kb, TF==tf & Target %in% common)
    out[,"TFinDAP"] <- (out$Target %in% TFs2test)
    
    out$TFsample <- gsub("newDAP.", "", out$TFsample)
    out$TFinDAP <- gsub("FALSE", "No", out$TFinDAP)
    out$TFinDAP <- gsub("TRUE", "Yes", out$TFinDAP)
    
    fileout[[tf]] <- out
  }
  fileout <- rbindlist(fileout, idcol = F)
  return(as_tibble(fileout))
}



TFs_Y1H_DAP_Net <- Count_TFs_Y1H_DAP_Net(TFs2test)
TFs_Y1H_DAP_Net["TFname"] <- ReplaceName(TFs_Y1H_DAP_Net$TF)

TFs_Y1H_DAP_Net <- TFs_Y1H_DAP_Net[,c("TF","TFname","TFsample", "Target", "Summit", "Dis", "OCR", "TFinDAP")]
colnames(TFs_Y1H_DAP_Net)[8] <- "TargWithDAP"

write.table(TFs_Y1H_DAP_Net, "TFs_Y1H_DAP_Net.txt", sep="\t", quote = F, row.names = F)



MakePromoterPlot <- function(target){
  
  # regulators
  df <- subset(TFs_Y1H_DAP_Net, Target==target)
  print(dim(df))
  name <- paste0(ReplaceName(target)," promoter")
  plot <- ggplot(df, aes(x=Dis, y=0)) +
    #geom_segment(aes(x = Dis, y = 0, xend = Dis, yend = 0.1)) +
    geom_point(color = "red") +
    geom_text_repel(aes(label=TFname, x=Dis, y = 0.00), 
                    
                    direction= "y",
                    point.padding = 0.5, 
                    size=3, 
                    force = 0.5, 
                    max.overlaps = Inf) +
    scale_y_continuous(expand = c(0,0)) +
    labs(title = name) +
    theme_pubclean() +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.text.y =  element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size=14)) +
    xlab("Promoter region (Kbs)  [0 = TSS] ") +
    ylab("")
    
    
      
    #plot <- ggpar(plot, xlab ="")
return(plot)
}

ReplaceName(TFs2test[9])
ReplaceName(TFs2test)

# size 3x8
MakePromoterPlot(TFs2test[9])
# size 3x8
MakePromoterPlot(TFs2test[8])
MakePromoterPlot(TFs2test[16])
                        


 #View(unique(subset(PDI_public, Method=="DAP")[,c(1,8)]))

#table(unique(PDI_public[,c(1,8)])$Method)

#subset(TFdic, TF.v4=="Zm00001d032923")
#subset(PDI, TF=="Zm00001d032923" & Target %in% PheGenes$GeneID)

################################################

