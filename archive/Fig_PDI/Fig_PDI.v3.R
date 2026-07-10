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
# CornCYC
CornCYC <- as_tibble(read.table("Data/Annotations/corn_pathways.0210325.reduced.txt", h=T, stringsAsFactors = F))
CornCYC <- CornCYC[!(CornCYC$Gene.id %in% "unknown"),]
CornCYC <- subset(CornCYC, Gene.id %in% Syntenic)
CornCYC[,"Class"] <- "Enzyme"
colnames(CornCYC)[3] <- "GeneID"
CornCYC$Pathway.id <- gsub("-", "_", CornCYC$Pathway.id)

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

write.table(PDInewDAP, "NewDAPseq_All_Peals_02.05.2020.txt", sep = "\t", row.names = F, quote = F)

PDInewDAP_1_kb <- subset(subset(PDInewDAP, Z > -0.5 & OCR ==1 & Dis <=0), Dis >= -1)

write.table(unique(PDInewDAP_1_kb[,c(1,3)]), 
            "TF_Target_NewDAP_1kb_promoter.txt", sep="\t", quote = F, row.names = F)

#View(unique(subset(PDI_public, Method=="DAP")[,c(1,8)]))
#table(unique(PDI_public[,c(1,8)])$Method)

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
  
  for (i in 1:nrow(TFdic)){
    ids <- gsub(TFdic$TF.v4[i], TFdic$TF.Name[i], ids)
  }
  return(ids)
}

ReplaceNamePWY <- function(ids){
  
  for (i in 1:nrow(CornCYC)){
    w <- paste0('\\<', CornCYC$Pathway.id[i], '\\>')
    ids <- gsub(w, CornCYC$Pathway.name[i], ids)
    #ids <- gsub("_", " ", ids)
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


## chop a string by a separator and return specified field
chop=function(myStr,mySep,myField){
  
  choppedString=sapply(strsplit(myStr,mySep),"[",myField)
  if(length(myField)>1){
    choppedString=apply(choppedString,2,function(x){paste0(x[!is.na(x)],collapse=mySep)})
  }
  return(choppedString)
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
PDI_public[, "TF"] <- sapply(strsplit(PDI_public$TFsample, split='.', fixed=TRUE), `[`, 2) # add methods label

write.table(subset(PDI_public[,c("TFsample","TF", "Target", "Method","Dis")], abs(Dis) <=2), 
            "GRASSIUS_ChIP_DAP.txt", row.names = F, quote = F)


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
LabelPeaksClasses
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

# read results from PDI and cis-eQTL overlapping
PDIclean_ciseQTL <- as_tibble(read.table("Data/PDI_data/Net.PDIclean.Full_cis.eQTL.txt", h=F, stringsAsFactors = F))
PDIclean_ciseQTL[,"Index"] <- paste(PDIclean_ciseQTL$V1, PDIclean_ciseQTL$V2, PDIclean_ciseQTL$V4, sep = ":")
PDIclean_ciseQTL[,"IndexSummit"] <- paste(PDIclean_ciseQTL$V10, 
                                          PDIclean_ciseQTL$V11, 
                                          PDIclean_ciseQTL$V13,
                                          sep = ":")

PDIclean_ciseQTL <- subset(PDIclean_ciseQTL, V1==V6)   # Same chromosome
PDIclean_ciseQTL <- PDIclean_ciseQTL[,-c(3,6,7,8,12)]  

# add annotated info from saf file for eQTl target
PDIclean_ciseQTL <- left_join(PDIclean_ciseQTL, saf, by=c("V5"="GeneID"))
PDIclean_ciseQTL[,"eQTL_TSS"] <- PDIclean_ciseQTL$TSS - PDIclean_ciseQTL$V2


# Remove not informative columns
PDIclean_ciseQTL <- PDIclean_ciseQTL[,-c(1:3,5,6,7,13,14)]
colnames(PDIclean_ciseQTL) <- c("eQTL.Target", "Sample",	"PDI.Target", "eQTL_summit",	"Index",	"IndexSummit", "eQTL_tss")

## Add summit index
PDIclean[,"IndexSummit"] <- paste(gsub("chr", "", PDIclean$Summit), PDIclean$TFsample, sep = ":")

# Reduce not-informative columns
PDIclean <- PDIclean[,c("TFsample", "Target", "Dis", "Method", "DisBin", "TF", "IndexSummit")]


# add cis-eQTL info
PDIclean[,"Class"] <- (PDIclean$IndexSummit %in% PDIclean_ciseQTL$IndexSummit)*1

PDIclean <- left_join(PDIclean, PDIclean_ciseQTL, by="IndexSummit")
PDIclean$eQTL.Target[is.na(PDIclean$eQTL.Target)] <- 0

# label PDI and eQTL supporting/predicting the same target gene
PDIclean$Class <- PDIclean$Class + (PDIclean$Target == PDIclean$eQTL.Target)*1
#table(PDIclean$Class)

# label peaks in close proximity to TSSs
PDIclean[,"DisClass"] <- (abs(PDIclean$Dis) <= 3)*1
PDIclean <- unique(PDIclean)

PDIclean$eQTL_tss <- PDIclean$eQTL_tss/1000

# ADD CLASSES LABEL FROM eQTL AND PDI OVERLAPING
PDIclean <- LabelPeaksClasses(PDIclean)
table(PDIclean$Class2)/1000

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

NetworkFinal <- subset(PDIclean, Class2 %in% c(1,2,4))[,c(1:4,6,7,9,16)]
PDIclean

# Combine network after correction by cis-eQTL
NetworkFinal <- unique(rbind(NetworkFinal, tem)[,c("TFsample","TF","Target", "Method")])

NetworkFinal <- subset(NetworkFinal, Target %in% Syntenic)
colnames(NetworkFinal)[1:2] <- c("TF", "TFid")

# Target summary table by method
OutDegree <- as_tibble(as.data.frame(table(NetworkFinal[,c(1,4)]), stringsAsFactors = F))
OutDegree <- subset(OutDegree, Freq > 0 )
#
OutDegree_low <- subset(OutDegree, Freq <= 50)
OutDegree_high <- subset(OutDegree, Freq >= 20000)


# Filter TFs with low number of interactions: targets < 50
NetworkFinal <- NetworkFinal[!(paste(NetworkFinal$TF,NetworkFinal$Method, sep = ":") %in% paste(OutDegree_low$TF, OutDegree_low$Method,sep = ":")),]


dim(unique(NetworkFinal[,c(1,3)]))
dim(unique(NetworkFinal))
length(unique(NetworkFinal$TF))
length(unique(NetworkFinal$Target)) 

OutDegree <- subset(OutDegree, Freq >50)
mean(OutDegree$Freq)
sd(OutDegree$Freq)/1000



## Figure S2a
Peaks_MethodDegree <- ggplot(OutDegree, aes(x=Freq, fill=Method)) +
  geom_density(alpha= 0.5) +
  #geom_vline(xintercept = median(Degree$Peaks/1000)) +
  theme_pubclean() +
  scale_x_continuous(label=comma) +
  xlab("Target genes") +
  ylab("Density")

Peaks_MethodDegree




##########
# Write results: Two types of networks, with and without cis-eQTL
##########

write.table(NetworkFinal, "PDI_NetworkFinal.10_11_2021.txt", sep = "\t", row.names = F, quote = F) 
write.table(NetworkEnrichment, "PDI_NetworkFinal.CornCYC.04_18_2021.txt", sep = "\t", row.names = F, quote = F) 

# Tagets only PDIs
NetworkFinal_NoCiseQTL <- subset(PDIclean, Class2 %in% c(1))#[,c(1:4,6,7,9,16)]

# Tagets only PDIs + cis-eQTL
NetworkFinal_YesCiseQTL <- subset(PDIclean, Class2 %in% c(2,4))[,c(1:4,6,7,9,16)]

# Replace with cis-eQTL target
tem <- subset(PDIclean, Class2 %in% c(3))[,c(1:4,6,7,9,16)]
tem$Target <- tem$eQTL.Target

# Combine network after correction by cis-eQTL
NetworkFinal_YesCiseQTL <- unique(rbind(NetworkFinal_YesCiseQTL, tem)[,c("TFsample","TF","Target", "Method")])


# Reduce network to target synthetic 
NetworkFinal_YesCiseQTL <- subset(NetworkFinal_YesCiseQTL, Target %in% Syntenic)
colnames(NetworkFinal_YesCiseQTL)[1:2] <- c("TF", "TFid")

NetworkFinal_NoCiseQTL <- NetworkFinal_NoCiseQTL[,c("TFsample","TF","Target", "Method")]
NetworkFinal_NoCiseQTL  <- subset(NetworkFinal_NoCiseQTL, Target %in% Syntenic)
colnames(NetworkFinal_NoCiseQTL)[1:2] <- c("TF", "TFid")

# Filter TFs with low number of interactions: targets < 50
mask <- paste(OutDegree_low$TF, OutDegree_low$Method,sep = ":")

tfsinyescis <- paste(NetworkFinal_YesCiseQTL$TF, NetworkFinal_YesCiseQTL$Method, sep = ":")
tfsinnocis <- paste(NetworkFinal_NoCiseQTL$TF, NetworkFinal_NoCiseQTL$Method, sep = ":")

NetworkFinal_YesCiseQTL <- unique(NetworkFinal_YesCiseQTL[!(tfsinyescis %in% mask),])
NetworkFinal_NoCiseQTL  <- unique(NetworkFinal_NoCiseQTL[!(tfsinnocis %in% mask),])

write.table(NetworkFinal_NoCiseQTL,  "Only_PDI_NetworkFinal.10_14_2022.txt", sep = "\t", row.names = F, quote = F)
write.table(NetworkFinal_YesCiseQTL, "CisE_PDI_NetworkFinal.10_14_2022.txt", sep = "\t", row.names = F, quote = F)

## Classes summary
# c("1", "2", "3", "4", "5")
# c("3 kb") 
# c("3 kb  & cis-eQTL")
# c("3 kb  & cis-eQTL New Targ.")
# c("50 kb & cis-eQTL")
# c("3 kb  & no-eQTL")

###
length(unique(NetworkEnrichment$TF))
length(unique(NetworkEnrichment$PWY))

####################################
## Heatmap with enrichment test   ##
####################################

NetEnr_Pval <- reshape2::dcast(NetworkEnrichment[,1:3], TF ~ PWY)
row.names(NetEnr_Pval) <- NetEnr_Pval$TF
NetEnr_Pval <- NetEnr_Pval[,-c(1)]
NetEnr_Pval[is.na(NetEnr_Pval)] <- 1

row.names(NetEnr_Pval) <- ReplaceName(row.names(NetEnr_Pval))
colnames(NetEnr_Pval) <- ReplaceNamePWY(colnames(NetEnr_Pval))

Hetmap_Matrix <- Heatmap(t(-log10(NetEnr_Pval)), name="-log10 (P-value)",
                         # column_km = 3,
                         column_names_rot = 45,
                         # row_names_rot = 45,
                         cluster_rows = TRUE, cluster_columns = TRUE,
                         show_column_dend = FALSE, show_row_dend = FALSE, 
                         #clustering_method_columns = "ward.D2",
                         #clustering_method_rows = "ward.D2",
                         col=viridis(4, direction = -1, option = "B"),
                         column_names_gp = gpar(fontsize = 5),
                         row_names_gp = gpar(fontsize = 6),
                         show_heatmap_legend = T,
                         # heatmap_ = unit(10),
                         heatmap_height = unit(13, 'cm'),
                         heatmap_width  = unit(26, 'cm'),
                         heatmap_legend_param = list(title_gp = gpar(fontsize = 10), 
                                                     labels_gp = gpar(fontsize = 10), 
                                                     direction = "vertical"))

#5x10
Hetmap_Matrix
####################################



######################################################################################
######### Count number of interactions for each class genes set of interest  #########
######################################################################################


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

################################################################################################
## Integration og PDI candidates with co-expression


# add info


