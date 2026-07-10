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

# CornCYC
CornCYC <- as_tibble(read.table("Data/Annotations/corn_pathways.0210325.reduced.txt", h=T, sep = '\t'))
CornCYC$Gene.id <- gsub('ZM', 'Zm', CornCYC$Gene.id)
CornCYC$Gene.id <- gsub('D', 'd', CornCYC$Gene.id)
CornCYC <- subset(CornCYC, Gene.id %in% Syntenic)

# Peak numbers
Peak <- as_tibble(read.table("Table_S4.txt", h=T, stringsAsFactors = F))

# PDI network
PDI <- as_tibble(read.table("Data/PDI_data/scATAC.Z.Full.Net.Dis2TSS.txt", h=T, stringsAsFactors = F))


#
PDI$Dis <- PDI$Dis/1000 # distance in kbs
PDI[, "Method"] <- sapply(strsplit(PDI$TFsample, split='.', fixed=TRUE), `[`, 1) # add methods label
table(PDI$Method)

PDI_public <- subset(PDI, Method != "newDAP")
PDInewDAP <- subset(PDI, Method == "newDAP")

#View(unique(subset(PDI_public, Method=="DAP")[,c(1,8)]))

#table(unique(PDI_public[,c(1,8)])$Method)

#subset(TFdic, TF.v4=="Zm00001d032923")
#subset(PDI, TF=="Zm00001d032923" & Target %in% PheGenes$GeneID)

# PDI class 
PDI.class  <- as_tibble(read.table("PDI_Class.txt", h=T, stringsAsFactors = F))

PDI.class
# Mediators
Mediator <- as_tibble(read.table("Data/Annotations/Mediators.txt", h=F, stringsAsFactors = F))
colnames(Mediator) <- c("GeneID", "Class")

# Define query list
# Filter 
List_query <- rbind(tibble(Target="Phe", GeneID=unique(PheGenes$GeneID)),
                    tibble(Target="TFs", GeneID=unique(TF_CoR$GeneID)), 
                    tibble(Target="TF.phe", GeneID=unique(Y1H$TF.v4)[is.na(unique(Y1H$TF.v4)) == FALSE]))

#
# List_query['Syntenic'] <- List_query$GeneID %in% Syntenic
List_query <- subset(List_query, GeneID %in% Syntenic) # reduce analysis to only Synthenic genes

List_query <- split(List_query$GeneID, List_query$Target)
CornCYC.list <- split(CornCYC$Gene.id, CornCYC$Pathway.id)
length(unique(CornCYC$Pathway.id))
#
List_query <- c(List_query, CornCYC.list)

#
####################################################################################
########################            Functions               ########################
####################################################################################

ReplaceName <- function(ids){
  
  for (i in 1:nrow(Top45)){
    ids <- gsub(Top45$V1[i], Top45$V2[i], ids)
  }
  
  for (i in 1:nrow(TFdic)){
    ids <- gsub(TFdic$TF.v4[i], TFdic$TF.Name[i], ids)
  }
  return(ids)
}


LabelPeaksClasses <- function(PDIdf){
  
  # class 1 No eQTL support and <3 kb from TSS
  Index_1 <- subset(PDIdf, Class==0 & DisClass==1)$IndexSummit
  
  # class 2 eQTL support and <3 kb from TSS and same eQTL and PDI target
  Index_2 <- subset(PDIdf, Class==2 & DisClass==1)$IndexSummit
  
  # class 3 eQTL support and <3 kb from TSS and different eQTL and PDI target
  Index_3 <- subset(PDIdf, Class==1 & DisClass==1)$IndexSummit
  
  Index_3 <- Index_3[!(Index_3 %in% Index_2)]
  
  # class 4 eQTL support and 3kb < eQTL_tss >=50 kb
  Index_4 <- subset(PDIdf, Class==2 & DisClass==0 & abs(eQTL_tss) <=50)$IndexSummit
  
  # class 5 all the other not in classes 1-4
  Selectec <- c(Index_1, Index_2, Index_3, Index_4)
  
  Index_5 <- subset(PDIdf, !(IndexSummit %in% Selectec))$IndexSummit
  
  PDIdf[,"Class2"] <- "NA"
  
  #
  PDIdf$Class2[PDIdf$IndexSummit %in% Index_5] <- 5
  PDIdf$Class2[PDIdf$IndexSummit %in% Index_4] <- 4
  PDIdf$Class2[PDIdf$IndexSummit %in% Index_3] <- 3
  PDIdf$Class2[PDIdf$IndexSummit %in% Index_2] <- 2
  PDIdf$Class2[PDIdf$IndexSummit %in% Index_1] <- 1
  
  return(PDIdf)
}

######################################################################################

######################################################################################
##################                      QC plots                    ##################
######################################################################################

########################################
####### ACR nd Z-score  analysis  
########################################

# histogram with Peaks by method
Degree %>%
  group_by(Method) %>%
  summarise(total=mean(Peaks))
  
# histogram with Peaks by TF
Degree <- as_tibble(as.data.frame(table(PDI_public$TFsample), stringsAsFactors = F))
colnames(Degree)[2] <- "Peaks"
Degree[,"Method"] <- sapply(strsplit(Degree$Var1, split='.', fixed=TRUE), `[`, 1) # add methods label


## Figure S2a
Hisplot_degree <- ggplot(Degree, aes(x=Peaks/1000, fill=Method)) +
  geom_density(alpha= 0.5) +
  #geom_vline(xintercept = median(Degree$Peaks/1000)) +
  theme_pubclean() +
  scale_x_continuous(label=comma) +
  xlab("Peaks (k)") +
  ylab("Density")

Hisplot_degree <- ggpar(Hisplot_degree, font.tickslab=14, font.x = 14, font.y = 14)


## Figure S2b
### Count peaks mapping to ACR
ACRs_method <- as_tibble(as.data.frame(table(PDI_public[,c("OCR", "Method")]), stringsAsFactors = F))
ACRs_method$Method <- factor(ACRs_method$Method, levels = c("DAP", "ChIP", "pChIP"))

PeaksInACR_plot <- ggplot(ACRs_method, aes(x=Method, y=Freq, fill=OCR)) +
  geom_bar(position="fill", stat="identity") +
  geom_text(aes(label=scales::comma(Freq)), position = position_fill(vjust = 0.5)) +
  theme_pubclean()  +
  xlab("") +
  ylab("Peak percentage") + 
  #scale_fill_discrete() +
  scale_fill_brewer(limits = c("0", "1"), labels = c("Not in ACRs", "In ACRs"), palette="BuPu") +
  guides(fill=guide_legend(title="Peaks"))

PeaksInACR_plot <- ggpar(PeaksInACR_plot, font.tickslab=14, font.x = 14, font.y = 14)

## Figure S2C
### Count peaks with low coverage

Zvalue_method <- PDI_public %>% 
  select(Method, Z) %>% 
  group_by(Method) %>% 
  count( Zf = Z > -0.5)

Zvalue_method$Method <- factor(Zvalue_method$Method, levels = c("DAP", "ChIP", "pChIP"))

PeaksInZvalue_plot <- ggplot(Zvalue_method, aes(x=Method, y=n, fill=Zf)) +
  geom_bar(position="fill", stat="identity") +
  geom_text(aes(label=scales::comma(n)), position = position_fill(vjust = 0.5)) +
  theme_pubclean()  +
  xlab("") +
  ylab("Peak percentage") + 
  #scale_fill_discrete() +
  scale_fill_brewer(limits = c("FALSE", "TRUE"), labels = c(expression("Z " <=" -0.5"), "Z > -0.5 "), palette="BuPu") +
  guides(fill=guide_legend(title="Peaks"))

PeaksInZvalue_plot <- ggpar(PeaksInZvalue_plot, font.tickslab=14, font.x = 14, font.y = 14)
PeaksInZvalue_plot



## Figure S2D
# reduce summit - TSS distance to bins of 10 kbps
breaks <- seq(min(-900), max(900), by = 10)
breaks <- cut(PDI_public$Dis, breaks = breaks)
breaks <- gsub("]", "", breaks)
breaks <- gsub("\\(", "", breaks)

###
PDI_public[, "Filter"] <- (((PDI_public$Z > -0.5)*1 + PDI_public$OCR) > 1)*1
PDI_public[,"DisBin"] <- as.numeric(sapply(strsplit(breaks, split=',', fixed=TRUE), `[`, 2)) # add bin label


# to reduce area of plot
#subset(PDI_public, abs(DisBin) <=100)
PDI_public$Method <- factor(PDI_public$Method, levels = c("DAP", "ChIP", "pChIP"))

table(subset(PDI_public, abs(DisBin) <=50 & Filter ==1)$DisBin)

Z_vs_Dis_plot <- ggplot(subset(PDI_public, abs(DisBin) <=200), 
                        aes(x=DisBin, y=Z,  color=as.factor(OCR))) +
  geom_smooth() +
  theme_pubclean() + 
  ylab("Z-score (Peak)") +
  xlab("Summit - TSS distance (kbp)") +
  scale_fill_discrete(limits = c("0", "1"), labels = c("Not in ACRs", "In ACRs"))+
  scale_color_discrete(limits = c("0", "1"), labels = c("Not in ACRs", "In ACRs")) +
  guides(fill=guide_legend(title="Peaks"), color=guide_legend(title="Peaks")) +
  facet_grid(~ Method ~ .) 

Z_vs_Dis_plot <- ggpar(Z_vs_Dis_plot, font.tickslab=14, font.x = 14, font.y = 14)

Fig_S2abcd <- ggarrange(ggarrange(Hisplot_degree, PeaksInACR_plot, PeaksInZvalue_plot, ncol = 1),
                        Z_vs_Dis_plot, ncol = 2)

Fig_S2abcd

###############################################
####### Target selection and cis-eQTL analyses
###############################################

# define final set of peaks to further experiments
PDIclean <- subset(PDI_public, Filter==1)
PDIclean[, "TF"] <- sapply(strsplit(PDIclean$TFsample, split='.', fixed=TRUE), `[`, 2) # add methods label
PDIclean


# Define targets based on distance to eQTLs
PDIclean_bed <- tibble(Chr=sapply(strsplit(PDIclean$Summit, split=':', fixed=TRUE), `[`, 1),
                       s=sapply(strsplit(PDIclean$Summit, split=':', fixed=TRUE), `[`, 2),
                       e=sapply(strsplit(PDIclean$Summit, split=':', fixed=TRUE), `[`, 2),
                       tf=PDIclean$TFsample,
                       CloseTarget=PDIclean$Target)

PDIclean_bed$Chr <- gsub('chr', "", PDIclean_bed$Chr)
write.table(PDIclean_bed, "../Data/PDI_data/Net.PDIclean.Full.bed", sep = "\t", row.names = F, col.names = F, quote = F)

# read results from PDI and cis-eQTL overlaping
PDIclean_ciseQTL <- as_tibble(read.table("Data/PDI_data/Net.PDIclean.Full_cis.eQTL.txt", h=F, stringsAsFactors = F))
PDIclean_ciseQTL[,"Index"] <- paste(PDIclean_ciseQTL$V1, PDIclean_ciseQTL$V2, PDIclean_ciseQTL$V4, sep = ":")
PDIclean_ciseQTL[,"IndexSummit"] <- paste(PDIclean_ciseQTL$V10, 
                                          PDIclean_ciseQTL$V11, 
                                          PDIclean_ciseQTL$V13,
                                          sep = ":")

PDIclean_ciseQTL <- subset(PDIclean_ciseQTL, V1==V6)
PDIclean_ciseQTL <- PDIclean_ciseQTL[,-c(3,6,7,8,12)]

# add annotated info from saf file for eQTl target
PDIclean_ciseQTL <- left_join(PDIclean_ciseQTL, saf, by=c("V5"="GeneID"))
PDIclean_ciseQTL[,"eQTL_TSS"] <- PDIclean_ciseQTL$TSS - PDIclean_ciseQTL$V2

# remove not informative columns
PDIclean_ciseQTL <- PDIclean_ciseQTL[,-c(1:3,5,6,7,13,14)]
colnames(PDIclean_ciseQTL) <- c("eQTL.Target", "Sample",	"PDI.Target", "eQTL_summit",	"Index",	"IndexSummit", "eQTL_tss")

## add summit index
PDIclean[,"IndexSummit"] <- paste(gsub("chr", "", PDIclean$Summit), PDIclean$TFsample, sep = ":")

# reduce not-informative columns
PDIclean <- PDIclean[,c("TFsample", "Target", "Dis", "Method", "DisBin", "TF", "IndexSummit")]

# add cis-eQTL info
PDIclean[,"Class"] <- (PDIclean$IndexSummit %in% PDIclean_ciseQTL$IndexSummit)*1

PDIclean <- left_join(PDIclean, PDIclean_ciseQTL, by="IndexSummit")
PDIclean$eQTL.Target[is.na(PDIclean$eQTL.Target)] <- 0

# label pdi and eQTL supporting/predicting the same target gene
PDIclean$Class <- PDIclean$Class + (PDIclean$Target == PDIclean$eQTL.Target)*1

# label peaks in close proximity to TSSs
PDIclean[,"DisClass"] <- (abs(PDIclean$Dis) <= 3)*1
PDIclean <- unique(PDIclean)

PDIclean$eQTL_tss <- PDIclean$eQTL_tss/1000

# ADD CLASSES LABEL FROM eQTL AND PDI OVERLAPING
PDIclean <- LabelPeaksClasses(PDIclean)

# Summary tables to describe: Total values
PeaksClean_classesTotal <- PDIclean %>% 
  #select(Class2) %>% 
  group_by(Class2) %>% 
  count( Class2 = Class2)
PeaksClean_classesTotal

# Summary tables to describe: values by method
PeaksClean_classes <- PDIclean %>% 
  select(Method, Class2) %>% 
  group_by(Class2) %>% 
  count( Method = Method)

# Plots total and total by method
PeaksClean_classesTotal$Class2 <- factor(PeaksClean_classesTotal$Class2, levels = c('1',"2", "3", "4", "5"))

Peaks_classesTotal_bar <- ggplot(PeaksClean_classesTotal, aes(x='Total', y=n, fill=forcats::fct_rev(Class2))) +
  geom_bar(position="fill", stat="identity") +
  geom_label_repel(aes(label=scales::comma(n)), 
            position = position_fill(vjust = 0.5), 
            label.size = 0.5,
            box.padding = unit(0.5, "lines"),
            show.legend = F) +
  theme_pubclean()  + xlab("") + ylab("Peak fraction") + 
  scale_fill_manual(limits = c("1", "2", "3", "4", "5"), 
                      labels = c(expression(" "<=" 3 kb"), 
                                 expression(" "<=" 3 kb & cis-eQTL"),
                                 expression(" "<=" 3 kb & cis-eQTL New Targ."), 
                                 expression(" "<=" 50 kb & cis-eQTL"),
                                 expression(" ">=" 3 kb & no-eQTL")),
                    values=c("mediumorchid1", "mediumturquoise", "goldenrod1", "deepskyblue1", "gray84")) +
  guides(fill=guide_legend(title="Peaks", ncol=1, title.hjust = 1,
                           label.position = 'left'))+
  theme(legend.position="left") +
  coord_flip()


PeaksClean_classes$Method <- factor(PeaksClean_classes$Method, levels = c("DAP", "ChIP", "pChIP"))

Peaks_classes_bar <- ggplot(PeaksClean_classes, aes(x=Method, y=n, fill=forcats::fct_rev(Class2))) +
  geom_bar(position="fill", stat="identity") +
  geom_label_repel(aes(label=scales::comma(n)), 
                   position = position_fill(vjust = 0.5), 
                   label.size = 0.5,
                   box.padding = unit(0.5, "lines"),
                   show.legend = F) +
  theme_pubclean()  + xlab("") + ylab("Peak fraction") + 
  scale_fill_manual(limits = c("1", "2", "3", "4", "5"), 
                    labels = c(expression(" "<=" 3 kb"), 
                               expression(" "<=" 3 kb & cis-eQTL"),
                               expression(" "<=" 3 kb & cis-eQTL New Targ."), 
                               expression(" "<=" 50 kb & cis-eQTL"),
                               expression(" ">=" 3 kb & no-eQTL")),
                    values=c("mediumorchid1", "mediumturquoise", "goldenrod1", "deepskyblue1", "gray84")) +
  guides(fill=guide_legend(title="Peaks", ncol=1, title.hjust = 1,
                           label.position = 'left'))+
  theme(legend.position="left") +
  coord_flip()


Peaks_classesTotal_bar <- ggpar(Peaks_classesTotal_bar, font.tickslab=14, font.x = 14, font.y = 14, font.legend = 12)
Peaks_classes_bar <- ggpar(Peaks_classes_bar, font.tickslab=14, font.x = 14, font.y = 14,  font.legend = 12)

########################################################################

########################################################################
#### Main supplementary Figure describing PDI network data source. #####
########################################################################

Fig_S2abcd <- ggarrange(ggarrange(Hisplot_degree, PeaksInACR_plot, PeaksInZvalue_plot, ncol = 1),
                        Z_vs_Dis_plot, ncol = 2)
#3x9
Fig_S2e <- ggarrange(Peaks_classesTotal_bar, Peaks_classes_bar, 
                     ncol = 1, common.legend = T, legend = 'left',
                     font.label = 12, align = "v")

Fig_S2abcd
########################################################################

########################################################################
#### Defining final network 
########################################################################

# Replace with cis-eQTL target
tem <- subset(PDIclean, Class2 %in% c(3))[,c(1:4,6,7,9,16)]
tem$Target <- tem$eQTL.Target

network <- subset(PDIclean, Class2 %in% c(1,2,4))[,c(1:4,6,7,9,16)]
network <- unique(rbind(network, tem)[,c("TF","Target", "Method")])

network <- subset(network, Target %in% Syntenic)

# Target summary table by method
TFdic


OutDegree <- as_tibble(as.data.frame(table(network[,c(1,3)]), stringsAsFactors = F))
OutDegree <- subset(OutDegree, Freq >0)
OutDegree_low <- subset(OutDegree, Freq <= 50)
OutDegree_high <- subset(OutDegree, Freq >= 20000)
hist(OutDegree$Freq, 100) 
OutDegree_low

# Filter TFs with low number of interactions: targets < 50
network <- network[!(paste(network$TF,network$Method, sep = ":") %in% paste(OutDegree_low$TF, OutDegree_low$Method,sep = ":")),]
network

dim(unique(network[,c(1:2)]))
length(unique(network$TF))
length(unique(network$Target))  

## Figure S2a
Peaks_Edegree <- ggplot(Degree, aes(x=Peaks/1000, fill=Method)) +
  geom_density(alpha= 0.5) +
  #geom_vline(xintercept = median(Degree$Peaks/1000)) +
  theme_pubclean() +
  scale_x_continuous(label=comma) +
  xlab("Peaks (k)") +
  ylab("Density")


##############################################################################
##########

PDIclean$TFsample

CornCYC.list <- split(CornCYC$GeneID, CornCYC$Pathway.id)

Enrichmet_classes <- function(network){
  ## Count TF targets in network
  # Count Total
  network <- subset(network, tgt.gid %in% Syntenic & reg.gid %in% Syntenic) 
  #
  Total_targtes <- as_tibble(as.data.frame(table(unique(network[,1:2])$reg.gid)))
  colnames(Total_targtes) <- c('TF', 'targets') 
  Total_targtes$TF <- as.character(Total_targtes$TF)
  
  
  # list input: network
  network.list <- unique(network[,1:2])
  network.list <- split(network.list$tgt.gid, network.list$reg.gid)
  
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
  #Pval_table <- tibble(TF=c("test", "a", "b"), TF2="test2")
  
  write.table(Pval_table, name, sep = "\t", row.names = F, quote = F)
  print(paste("... Done: ", name, " ...", sep = ""))
  #return(Pval_table)
}



######################################################################################
######### Count number of interactions for each class genes set of interest  #########
######################################################################################

Enrichmet.targets <- function(network){
  ## Count TF targets in network
  # Count Total
  #
  network <- unique(network[,c(2,3)])
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
  
  # Pval$Phe <- p.adjust(Pval$Phe, method = 'fdr')
  # Pval$TF.phe <- p.adjust(Pval$TF.phe, method = 'fdr')
  # Pval$TFs <- p.adjust(Pval$TFs, method = 'fdr')
  
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
TF.Candidates <- Enrichmet.targets(PDI)

#TF.Candidates
#subset(TFdic, TF.v4=="Zm00001d032923")
#subset(PDI, TF=="Zm00001d032923" & Target %in% PheGenes$GeneID)
#subset(TF.Candidates, TF=="Zm00001d032923")

# add class PDI  
TF.Candidates <- left_join(TF.Candidates, PDI.class, by='TF')

# PDI <- subset(PDI, Target %in% Syntenic)


########  matrix input to heatmap
write.table(TF.Candidates, "PDI_candidates.06-09-21.txt", sep = '\t', row.names = F, quote = F)

TF.Candidates <- as_tibble(read.table("~/Google_Drive/My Drive/Fabio_Mac/MSU_projects/2020/MaizeENCODE/Project_Folder_eQTL_PDI_CoExp/Data/Candidates_files/PDI_candidates.06-09-21.txt", h=T, stringsAsFactors = F))
#length(unique(subset(TF.Candidates_table.Summary, class == 'TF.phe' & Freq>0)$TF))

rm(TF.Candidates)

################################################  

################################################
#### Heatmap TF vs pathway
################################################
TF.Candidates

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

Plot.TF_pathway <- ggplot(TF.Candidates, aes(x=reorder(TF, n.targ), 
                                             y=reorder(Pathway.name.V2, Genes), fill=-log10(padj)))+
  geom_tile()  +
  #scale_x_continuous(breaks = seq(0,100,20), expand = c(0,0)) + 
  scale_fill_viridis(option = "D", direction = 1) + #
  theme_bw() +
  theme(legend.title = element_text(vjust = 0.8), 
        panel.border = element_blank(), 
        panel.grid.major = element_blank())

Plot.TF_pathway <- ggpar(Plot.TF_pathway, ylab = "Pathway", xlab = "TF", x.text.angle = 90, font.xtickslab = 2)

typeof(as.data.frame(TF.Candidates)[,3])
colnames()

######
TF.Candidates[,"TFname"] <- TF.Candidates$TF

for (i in 1:nrow(Top45)){
  TF.Candidates$TFname <- gsub(Top45$V1[i], Top45$V2[i], TF.Candidates$TFname)
}

for (i in 1:nrow(TFdic)){
  TF.Candidates$TFname <- gsub(TFdic$TF.v4[i], TFdic$TF.Name[i], TF.Candidates$TFname)
}


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

write.table(unique(subset(PDI, TF=='Zm00001d021019' & Target %in% unique(PheGenes$GeneID))[,2:3]),
            "bHLH136.Phe.targets.pdi.txt", row.names = F, sep = '\t', quote = F)

top.phe <- c("Zm00001d047017", "Zm00001d006236", "Zm00001d021019")
write.table(unique(subset(PDI, TF %in% top.phe & Target %in% unique(PheGenes$GeneID))[,2:3]),
            "Top.Phe.targets.pdi.txt", row.names = F, sep = '\t', quote = F)




####################################################################################################
# CoExp.phe.cand <- as.character(subset(CoexpressionCandidates, class=='Phe' & Freq>0)$TF)
# CoExp.TFphe.cand <- as.character(subset(CoexpressionCandidates, class=='TF.phe' & Freq>0)$TF)
# #
# TF.Candidates_table.Summary.Phe <- subset(TF.Candidates, TF %in% CoExp.phe.cand & class=='Phe')[,c(2,6)]
# TF.Candidates_table.Summary.TFPhe <- subset(TF.Candidates, TF %in% CoExp.TFphe.cand & class=='TF.phe')[,c(2,6)]
# #
# TF.Candidates_table.Summary.Phe <- as_tibble(as.data.frame(table(TF.Candidates_table.Summary.Phe)))
# TF.Candidates_table.Summary.TFPhe <- as_tibble(as.data.frame(table(TF.Candidates_table.Summary.TFPhe)))
# #
# TF.Candidates_table.Summary.Phe$class <- as.character(TF.Candidates_table.Summary.Phe$class)
# TF.Candidates_table.Summary.TFPhe$class <- as.character(TF.Candidates_table.Summary.TFPhe$class)
# 
# TF.Candidates_table.Summary.Phe <- subset(TF.Candidates_table.Summary.Phe, Freq>0)
# TF.Candidates_table.Summary.TFPhe <- subset(TF.Candidates_table.Summary.TFPhe, Freq>0)
# #
# Plot.summary.phe.coexp <- ggplot(TF.Candidates_table.Summary.Phe, 
#                        aes(y=Freq, x=source))+
#   geom_bar(stat="identity", position=position_dodge(), fill='indianred1')  +
#   #scale_x_continuous(breaks = seq(0,100,20), expand = c(0,0)) + 
#   #scale_fill_viridis(option = "D", direction = -1) + #
#   theme_pubclean() 
# 
# #
# Plot.summary.TFPhe.coexp <- ggplot(TF.Candidates_table.Summary.TFPhe, 
#                                  aes(y=Freq, x=source))+
#   geom_bar(stat="identity", position=position_dodge(), fill='darkturquoise')  +
#   #scale_x_continuous(breaks = seq(0,100,20), expand = c(0,0)) + 
#   #scale_fill_viridis(option = "D", direction = -1) + #
#   theme_pubclean() 
# 
# 
# Plot.summary.phe.coexp <- ggpar(Plot.summary.phe.coexp, ylab='TFs')
# Plot.summary.TFPhe.coexp <- ggpar(Plot.summary.TFPhe.coexp, ylab='TFs')
# 
# ggarrange(Plot.summary.phe.coexp, Plot.summary.TFPhe.coexp, nrow = 1, 
#           align = 'h', labels = c('Phe', 'TF.phe'))
####################################################################################################



