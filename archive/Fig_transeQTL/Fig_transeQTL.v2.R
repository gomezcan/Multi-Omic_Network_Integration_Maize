library(tidygraph)
library(tidyverse)
library(ggrepel)
library(ggpubr)
library(data.table)
library(GeneOverlap)
library(viridis)
library(ComplexHeatmap)
library(fgsea)
library(reshape)
library(scales)


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

# Mediators
Mediator <- as_tibble(read.table("Data/Annotations/Mediators.txt", h=F, stringsAsFactors = F))
colnames(Mediator) <- c("GeneID", "Class")

# TF and CoReg
TF_CoR <- as_tibble(read.table("Data/Annotations/TF_CoR_Mazie.txt", h=T, stringsAsFactors = F)) 
TF_CoR[, "Class"] <- "TF"

CoRegs <- c("BSD", "DDT", "FHA", "GNAT", "HMG", "IWS1,SPN1", "JUMONJI", "LEUNIG", 
            "LIM", "MBF1", "MED", "p15", "PC4", "Sub1", "PHD", "RBR", "SNF2", "SWI", "SNF", 
            "BAF60, SWI", "SNF", "SWI3", "TAZ", "TRAF", "Ultrapetala", "WD40")

# Label CoRegs
for (i in CoRegs){
  TF_CoR$Class[grepl(i, TF_CoR$Family)] <- "CoReg"
}

TF_CoR <- TF_CoR[!(TF_CoR$GeneID %in% Mediator$GeneID),]

## Y1H network
Y1H <- as_tibble(read.table("Data/Annotations/Y1H.data.txt", h=T, stringsAsFactors = F))[,4:6]

#table(subset(PDI.class, TF %in% Y1H$TF.v4)$source)


#TF	target	weight	Random.mean, Z, Index

## syntenic genes 
Syntenic <- as_tibble(read.table("Data/Annotations/Zm.v4.synteny.genes.txt", h=T, stringsAsFactors = F, sep = '\t'))
Syntenic <- subset(Syntenic, gene_synteny=='syntenic')$gene_id

# CornCYC
CornCYC <- as_tibble(read.table("Data/Annotations/corn_pathways.0210325.reduced.txt", h=T, stringsAsFactors = F))
CornCYC <- CornCYC[!(CornCYC$Gene.id %in% "unknown"),]
CornCYC <- subset(CornCYC, Gene.id %in% Syntenic)
#CornCYC$Gene.id[grepl("Zm", CornCYC$Gene.id) == FALSE]
CornCYC[,"Class"] <- "Enzyme"
colnames(CornCYC)[3] <- "GeneID"
CornCYC$Pathway.id <- gsub("-", "_", CornCYC$Pathway.id)

# kinases
kinases <- as_tibble(read.table("Data/Annotations/kinases_maize_AGPv4.txt", h=T, stringsAsFactors = F))
kinases <- subset(kinases, GeneID %in% Syntenic)
kinases[,"Class"] <- "kinases"

Classes <- rbind(unique(TF_CoR[,c("GeneID", "Class")]), 
                 unique(Mediator[,c("GeneID", "Class")]), 
                 unique(kinases[,c("GeneID", "Class")]),
                 unique(CornCYC[,c("GeneID", "Class")]))

Classes <- Classes[(Classes$GeneID %in% Syntenic),]


#as.data.frame(table(as.data.frame(table(CornCYC$Pathway.id))$Freq))
#as.data.frame(table(as.data.frame(table(kinases$Family))$Freq))

unique(CornCYC$Pathway.name)

# define query list
# filter 
List_query <- rbind(tibble(Target="Phe", GeneID=unique(PheGenes$GeneID)),
                    tibble(Target="TFs", GeneID=unique(TF_CoR$GeneID)), 
                    tibble(Target="TF.phe", GeneID=unique(Y1H$TF.v4)[is.na(unique(Y1H$TF.v4)) == FALSE]))
#
# List_query['Syntenic'] <- List_query$GeneID %in% Syntenic
List_query <- subset(List_query, GeneID %in% Syntenic) # reduce analysis to only Syntenic genes

List_query <- split(List_query$GeneID, List_query$Target)
CornCYC.list <- split(CornCYC$Gene.id, CornCYC$Pathway.id)
length(unique(CornCYC$Pathway.id))
#
table(as.data.frame(table(CornCYC$Pathway.id))$Freq)

List_query <- c(List_query, CornCYC.list)


######################################################################################
#########       Read different eQTL data sets and define summary tables      #########
######################################################################################

# Total cis- & trans-eQTLs
#eQTLs <- as_tibble(read.table("Clean_eQTL_Syntenic.txt", stringsAsFactors = F, h=T))

# cis-eQTLt
cis.eQTLt <- as_tibble(read.table("Clean_cis.eQTLt.txt", stringsAsFactors = F, h=T))

# cis-eQTL
#cis.eQTL <- as_tibble(read.table("Clean_cis.eQTL.txt", stringsAsFactors = F, h=T))

# trans-eQTL
trans.eQTL <- as_tibble(read.table("Clean_trans.eQTL.txt", stringsAsFactors = F, h=T))
# trans-eQTLp
trans.eQTLp <- as_tibble(read.table("Clean_trans.eQTLp.txt", stringsAsFactors = F, h=T))


# Confirm unique values
trans.eQTL <- subset(trans.eQTL, !(Index %in% cis.eQTLt$Index))
cis.eQTL   <- subset(cis.eQTL, !(Index %in% cis.eQTLt$Index))

table(cis.eQTL$Index %in% cis.eQTLt$Index)   # cis in cis target?
table(cis.eQTL$Index %in% trans.eQTLp$Index) # cis in trans promoter?
table(cis.eQTLt$Index %in% trans.eQTL$Index) # cis target in trans?
table(trans.eQTLp$Index %in% cis.eQTL$Index) # trans promoter in cis?

# Count eQTL not assigned
eQTLs <- subset(eQTLs, !(Index %in% trans.eQTL$Index))
eQTLs <- subset(eQTLs, !(Index %in% trans.eQTLp$Index))
eQTLs <- subset(eQTLs, !(Index %in% cis.eQTL$Index))
eQTLs <- subset(eQTLs, !(Index %in% cis.eQTLt$Index))

length(unique(eQTLs$Index))
length(unique(trans.eQTL$Index))
length(unique(trans.eQTLp$Index))
length(unique(cis.eQTL$Index))
length(unique(cis.eQTLt$Index))

# trans-eQTL network
teQTLs <- subset(teQTLs, Target %in% Syntenic)

######################################################################################

######################################################################################
################## QC plots       
######################################################################################


net_trans_eQTL <- unique(trans.eQTL[,c("source", "Target")])
net_trans_eQTLp <- unique(trans.eQTLp[,c("source", "Target")])

# Reduce to only syntenic genes
net_trans_eQTL <- subset(net_trans_eQTL, source %in% Syntenic & Target %in% Syntenic)


# histogram with Peaks by TF
Degree_teQTL <- as_tibble(as.data.frame(table(net_trans_eQTL$source), stringsAsFactors = F))
colnames(Degree_teQTL) <- c("Source", "Targets")
Degree_teQTL[,"isTF"] <- Degree_teQTL$Source %in% TF_CoR$GeneID


Degree_teQTLp <- as_tibble(as.data.frame(table(net_trans_eQTLp$source), stringsAsFactors = F))
colnames(Degree_teQTLp) <- c("Source", "Targets")
Degree_teQTLp[,"isTF"] <- Degree_teQTLp$Source %in% TF_CoR$GeneID



####################################
#### Figure S3b : teQTL degree
####################################

Hisplot_degree.teQTL <- ggplot(Degree_teQTL, aes(x=Targets, fill=isTF)) +
  geom_histogram(alpha= 0.5, bins = 50) +
  #geom_vline(xintercept = median(Degree$Peaks/1000)) +
  theme_pubclean() +
  scale_y_continuous(label=comma, expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  xlab("Targets") +
  ylab("Count") +
  scale_fill_discrete(limits = c("FALSE", "TRUE"), 
                    labels = c("Other", "TFs")) +
  guides(fill=guide_legend(title="Gene source"))



## Figure S3b part2: bar plot with fraction of data
temFre_tQTL <- as.data.frame(table(Degree_teQTL$isTF))

Hisplot_degree.teQTL_bar <- ggplot(temFre_tQTL, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(position="fill", stat="identity", alpha= 0.5) +
  geom_text(aes(label=scales::comma(Freq)), position = position_fill(vjust = 0.5)) +
  theme_pubclean()  +
  scale_y_continuous(expand = c(0,0)) +
  xlab("Gene\nsource") +
  ylab("trans-eQTL fraction") + 
  scale_fill_discrete(limits = c("FALSE", "TRUE"), labels = c("Other", "TFs")) +
  guides(fill=guide_legend(title="Gene source"))



####################################
#### Figure S3c: teQTLp
####################################

temFre_tQTLp <- as.data.frame(table(Degree_teQTLp$isTF))

subset(Degree_teQTL, Source %in% Syntenic)

Hisplot_degree.teQTLp <- ggplot(Degree_teQTLp, aes(x=Targets, fill=isTF)) + # isTF
  geom_histogram(alpha= 0.5, bins = 50) +
  #geom_vline(xintercept = median(Degree$Peaks/1000)) +
  theme_pubclean() +
  scale_y_continuous(label=comma, expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  xlab("Targets") +
  ylab("Count") +
  scale_fill_discrete(limits = c("FALSE", "TRUE"), labels = c("Other", "TFs")) +
  guides(fill=guide_legend(title="Gene source"))


## Figure S3c part2: bar plot with fraction of data
temFre_tQTLp <- as.data.frame(table(Degree_teQTLp$isTF))

Hisplot_degree.teQTLp_bar <- ggplot(temFre_tQTLp, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(position="fill", stat="identity", alpha= 0.5) +
  geom_text(aes(label=scales::comma(Freq)), position = position_fill(vjust = 0.5)) +
  theme_pubclean()  +
  scale_x_continuous(expand = c(0,0)) +
  xlab("Gene\nsource") +
  ylab("trans-eQTL fraction") + 
  scale_fill_discrete(limits = c("FALSE", "TRUE"), 
                      labels = c("Other", "TFs")) +
  guides(fill=guide_legend(title="Gene source"))

##### merge figures

Hisplot_degree.teQTL <- ggpar(Hisplot_degree.teQTL, font.tickslab=14, font.x = 14, font.y = 14)
Hisplot_degree.teQTL_bar <- ggpar(Hisplot_degree.teQTL_bar, font.tickslab=14, font.x = 14, font.y = 14)
Hisplot_degree.teQTLp      <- ggpar(Hisplot_degree.teQTLp, font.tickslab=14, font.x = 14, font.y = 14)
Hisplot_degree.teQTLp_bar <- ggpar(Hisplot_degree.teQTLp_bar, font.tickslab=14, font.x = 14, font.y = 14)

Fig_S3b <- ggarrange(Hisplot_degree.teQTL, Hisplot_degree.teQTL_bar, 
                     ncol = 2, widths = c(3,1), common.legend = T, align = 'h')

Fig_S3c <- ggarrange(Hisplot_degree.teQTLp, Hisplot_degree.teQTLp_bar, 
                     ncol = 2, widths = c(3,1), common.legend = T, align = 'h')

Fig_S3bc <- ggarrange(Fig_S3b, Fig_S3c, nrow = 1, common.legend = T, align = 'h')

# size 4x12
Fig_S3bc

## 
TF_CoR

hist(net_trans.eQTL_freq$Freq, 100)

subset(net_trans_eQTL, source %in% TF_CoR$GeneID)
subset(net_trans_eQTLp, source %in% TF_CoR$GeneID)

subset(Degree_teQTL, Targets >=40) %>%
  arrange(Targets)

net_trans_eQTL

temInterestlist <- c("Zm00001d046170", 
                     subset(net_trans_eQTL, source=="Zm00001d046170")$source,
                     subset(net_trans_eQTL, source=="Zm00001d046170")$Target)


subset(net_trans_eQTL, Target=="Zm00001d047671")
subset(net_trans_eQTLp, Target=="Zm00001d047671")



Degree_teQTL %>% 
  summarise(Degree=mean(Targets))

Degree_teQTL

Get_target_degree <- function(df_net, degreedf){
  
  tfs_t <- unique(subset(df_net, source %in% TF_CoR$GeneID)$Target) # TFs targets
  
  targ_tfs <- unique(subset(df_net, source %in% TF_CoR$GeneID)$Target)
  targ_Others <- unique(subset(df_net, !(source %in% TF_CoR$GeneID))$Target)
  
  print(table(targ_tfs %in% TF_CoR$GeneID))
  print(table(targ_Others %in% TF_CoR$GeneID))
    
  tfs_t <- data.frame(Class="TF",     Value=subset(degreedf,  Source %in% tfs_t)$Targets)
  other_t<- data.frame(Class="Other", Value=subset(degreedf, !(Source %in% tfs_t))$Targets)
  
  
  return(as_tibble(rbind(tfs_t, other_t)))
}

L2_degree_teQTL <- Get_target_degree(net_trans_eQTL, Degree_teQTL)

net_trans_eQTL

ggplot(L2_degree_teQTL, aes(x=Class, y=Value))+
  geom_violin()



######################################################################################
######### Count number of interactions for each class genes set of interest  #########
######################################################################################

Enrichmet.targets <- function(network){
  ## Count TF targets in network
  # Count Total
  #
  network <- unique(network[,c('TF','Target')])
  # Total 
  Total_targtes <- as_tibble(as.data.frame(table(network$TF)))
  colnames(Total_targtes) <- c('TF', 'targets') 
  Total_targtes$TF <- as.character(Total_targtes$TF)
  
  # list input: network
  network.list <- split(network$Target, network$TF)
  
  ## Compare list of predicted targets vs annoated genes in query.vector 
  # length(Syntenic) : 
  go.obj <- newGOM(network.list, List_query, genome.size=24622) # annotated genes in Genome v4
  
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

# Identify candidates
TF.Candidates <- Enrichmet.targets(teQTLs)



########  matrix input to heatmap
write.table(TF.Candidates, "PDI_candidates.06-09-21.txt", sep = '\t', row.names = F, quote = F)

TF.Candidates <- as_tibble(read.table("~/Google_Drive/My Drive/Fabio_Mac/MSU_projects/2020/MaizeENCODE/Project_Folder_eQTL_PDI_CoExp/Data/Candidates_files/PDI_candidates.06-09-21.txt", h=T, stringsAsFactors = F))
#length(unique(subset(TF.Candidates_table.Summary, class == 'TF.phe' & Freq>0)$TF))


######################################################################################
######### Summary and plots
######################################################################################

################################################
#### Pathway enriched by TF
################################################
#### Selection of TFs predicted as regulators of Phe and other
# TF.Source.Freq <- subset(as_tibble(as.data.frame(table(TF.Candidates[,c(1,6)]), stringsAsFactors = F)), Freq>0)
# #
# Source.Freq <- as.data.frame(table(TF.Source.Freq[,c(2,3)]))
# colnames(Source.Freq) <- c("Method", "Pathway", 'Freq')
# 
# Pathway.Dis.Method.plot <- ggplot(Source.Freq, aes(y=Freq, x=Pathway, fill=Method))+
#   geom_bar(stat="identity", position=position_dodge())  +
#   theme_pubclean()
# Pathway.Dis.Method.plot   <- ggpar(Pathway.Dis.Method.plot, xlab = "Pathways Enriched by TF")

################################################  

################################################
#### Heatmap TF vs pathway
################################################

## names
CornCYC.DB <- unique(CornCYC[,1:2])
CornCYC.DB <-  rbind(CornCYC.DB, tibble(Pathway.id=c("Phe","TF.phe", "TFs"), 
                                        Pathway.name=c("Phe","TF.phe", "TFs")))
# size
List_query.DB.size <- data.frame(t(as.data.frame(lapply(List_query, length))))
List_query.DB.size[,"Pathway.id"] <- row.names(List_query.DB.size)
colnames(List_query.DB.size)[1] <- "Genes"
List_query.DB.size <- as_tibble(List_query.DB.size)
List_query.DB.size$Pathway.id <- gsub("\\.", "-", List_query.DB.size$Pathway.id)
List_query.DB.size$Pathway.id[2] <- "TF.phe"

CornCYC.DB <- left_join(CornCYC.DB, List_query.DB.size, by='Pathway.id')
CornCYC.DB[,"Pathway.name.V2"] <- paste(CornCYC.DB$Pathway.name, " [", CornCYC.DB$Genes, "]", sep = "")

#TF.Candidates <- TF.Candidates[,-c(7)]
TF.Candidates <- left_join(TF.Candidates, CornCYC.DB[,-c(2)], by= c('class'="Pathway.id"))
TF.Candidates[,"TFname"] <- ReplaceName(TF.Candidates$TF)
TF.Candidates$TF

Plot.TF_pathway <- ggplot(TF.Candidates, aes(x=reorder(TFname, n.targ), 
                                             y=reorder(Pathway.name.V2, Genes), fill=-log10(padj)))+
  geom_tile()  +
  #scale_x_continuous(breaks = seq(0,100,20), expand = c(0,0)) + 
  scale_fill_viridis(option = "D", direction = 1) + #
  theme_bw() +
  theme(legend.title = element_text(vjust = 0.8), 
        panel.border = element_blank(), 
        panel.grid.major = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

Plot.TF_pathway <- ggpar(Plot.TF_pathway, ylab = "Pathway", xlab = "TF", x.text.angle = 90, font.xtickslab = 10)
Plot.TF_pathway

######


Phe.candidates <- subset(TF.Candidates, class=="Phe")$TF

Top45[Top45$V1 %in% Phe.candidates[Phe.candidates %in% Top45$V1],]
Y1H[Y1H$TF.v4 %in% Phe.candidates[Phe.candidates %in% Y1H$TF.v4],]

TFdic[TFdic$TF.v4 %in% Phe.candidates[Phe.candidates %in% Y1H$TF.v4],]

subset(TF.Candidates, TF %in% c("Zm00001d048208") & class=='Phe')

unique(subset(TF.Candidates, TF %in% Phe.candidates)$TFname)


Plot.TF_pathway.Phe <- ggplot(subset(TF.Candidates, TF %in% Phe.candidates), 
                              aes(x=reorder(TFname, n.targ), 
                                  y=reorder(Pathway.name.V2, Genes), fill=-log10(padj)))+
  geom_tile()  +
  #scale_x_continuous(breaks = seq(0,100,20), expand = c(0,0)) + 
  scale_fill_viridis(option = "D", direction = 1) + #
  theme_bw() +
  theme(legend.title = element_text(vjust = 0.8), 
        panel.border = element_blank(), 
        panel.grid.major = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

Plot.TF_pathway.Phe <- ggpar(Plot.TF_pathway.Phe, ylab = "Pathway", xlab = "TF", x.text.angle = 90, font.xtickslab = 6)
Plot.TF_pathway.Phe


################################################################################################
## Integration og PDI candidates with co-expression

CoexpressionCandidates <- as_tibble(read.table("../Fig_Coexpression/CoExpression_candidates.fullReport.06-10-21.txt", h=T)) 
CoexpressionCandidates <- subset(CoexpressionCandidates, TF %in% Phe.candidates)
CoexpressionCandidates.Reduced <- as_tibble(as.data.frame(table(CoexpressionCandidates[,c(2:3)]), stringsAsFactors = F))

# total set of TF-pathways
PDI.CoExp.phe <- unique(rbind(subset(CoexpressionCandidates.Reduced, Freq >0)[,1:2],
                              subset(TF.Candidates, TF %in% Phe.candidates)[,1:2]))

PDI.CoExp.phe[,"PDI"] <- (paste(PDI.CoExp.phe$TF, PDI.CoExp.phe$class, sep = "_") %in%
                            paste(subset(TF.Candidates, TF %in% Phe.candidates)[,1:2]$TF, subset(TF.Candidates, TF %in% Phe.candidates)[,1:2]$class, sep = "_"))*1

PDI.CoExp.phe[,"CoExp"] <- (paste(PDI.CoExp.phe$TF, PDI.CoExp.phe$class, sep = "_") %in%
                              paste(subset(CoexpressionCandidates.Reduced, Freq >0)[,1:2]$TF, subset(CoexpressionCandidates.Reduced, Freq >0)[,1:2]$class, sep = "_"))*2
PDI.CoExp.phe[,"Support"] <- PDI.CoExp.phe$PDI + PDI.CoExp.phe$CoExp

# add info
PDI.CoExp.phe <- left_join(PDI.CoExp.phe, CornCYC.DB[,-c(2)], by= c('class'="Pathway.id"))

PDI.CoExp.phe[, "TFname"] <- ReplaceName(PDI.CoExp.phe$TF)

PDI.CoExp.phe


# PLOT
Plot.TF_pathway.and.Coexp.Phe <- ggplot(PDI.CoExp.phe, 
                                        aes(x=TFname, y=Pathway.name.V2, 
                                            fill=as.character(Support)))+
  geom_tile()  +
  #scale_x_continuous(breaks = seq(0,100,20), expand = c(0,0)) + 
  scale_fill_viridis(option = "D", direction = 1, discrete = T) + #
  theme_bw() +
  theme(legend.title = element_text(vjust = 0.8), 
        panel.border = element_blank(), 
        panel.grid.major = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

Plot.TF_pathway.and.Coexp.Phe

write.table(unique(subset(PDI, TF=='Zm00001d047017' & Target %in% unique(PheGenes$GeneID))[,2:3]),
            "bHLH91.Phe.targets.pdi.txt", row.names = F, sep = '\t', quote = F)

write.table(unique(subset(PDI, TF=='Zm00001d006236' & Target %in% unique(PheGenes$GeneID))[,2:3]),
            "MYB31.Phe.targets.pdi.txt", row.names = F, sep = '\t', quote = F)



