library(tidygraph)
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
library(patchwork)


ReplaceNamePWY <- function(ids){
  
  for (i in 1:nrow(CornCYC)){
    w <- paste0('\\<', CornCYC$Pathway.id[i], '\\>')
    ids <- gsub(w, CornCYC$Pathway.name[i], ids)
    #ids <- gsub("_", " ", ids)
  }
  return(ids)
}

ReplaceName <- function(ids){
  
  for (i in 1:nrow(TFdic)){
    ids <- gsub(TFdic$TF.v4[i], TFdic$TF.Name[i], ids)
  }
  return(ids)
}

##################################################
##########        Annotations       ##############
##################################################

saf <- as_tibble(read.table("Data/eQTL_data/Zea_mays.B73_RefGen_v4.46.saf", stringsAsFactors = F))

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

##
TF_CoR <- TF_CoR[!(TF_CoR$GeneID %in% Mediator$GeneID),]

## Y1H network
Y1H <- as_tibble(read.table("Data/Annotations/Y1H.data.txt", h=T, stringsAsFactors = F))[,4:6]

## TF	target	weight	Random.mean, Z, Index

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

table(data.table(table(CornCYC$GeneID))$N)

# kinases
kinases <- as_tibble(read.table("Data/Annotations/kinases_maize_AGPv4.txt", h=T, stringsAsFactors = F))
kinases <- subset(kinases, GeneID %in% Syntenic)
kinases[,"Class"] <- "kinases"

Classes <- rbind(unique(TF_CoR[,c("GeneID", "Class")]), 
                 unique(Mediator[,c("GeneID", "Class")]), 
                 unique(kinases[,c("GeneID", "Class")]),
                 unique(CornCYC[,c("GeneID", "Class")]))

Classes <- Classes[(Classes$GeneID %in% Syntenic),]
Classes_Freq <- as.data.table(table(Classes$Class))


# define query list
# filter 
List_query <- rbind(tibble(Target="Phe", GeneID=unique(PheGenes$GeneID)),
                    tibble(Target="TFs", GeneID=unique(TF_CoR$GeneID)), 
                    tibble(Target="TF.phe", GeneID=unique(Y1H$TF.v4)[is.na(unique(Y1H$TF.v4)) == FALSE]))

# reduce analysis to only Syntenic genes
List_query <- subset(List_query, GeneID %in% Syntenic) 
List_query <- split(List_query$GeneID, List_query$Target)
CornCYC.list <- split(CornCYC$Gene.id, CornCYC$Pathway.id)

#
List_query <- c(List_query, CornCYC.list)
lapply(List_query, length)

######################################################################################
#########       Read different eQTL data sets and define summary tables      #########
######################################################################################

# trans-eQTL: "V2 included trans-eqTL from second part files"
trans.eQTL <- fread("Data/eQTL_data/Clean_trans.eQTL.v2.txt", stringsAsFactors = F, h=T)

# trans-eQTLp
trans.eQTLp <- fread("Data/eQTL_data/Clean_trans.eQTLp.v2.txt", stringsAsFactors = F, h=T)

# Combine trans-eQTL and trans-eQTLp
trans_eQTL <- unique(rbind(trans.eQTL[,c("Target", "source", "Index")], trans.eQTLp[,c("Target", "source", "Index")]))

# Reduce interactions to source-target interactions
Net_teQTL <- as.data.table(table(trans_eQTL[,c('source', 'Target')]))
Net_teQTL <- subset(Net_teQTL, N > 0)
Net_teQTL <- subset(Net_teQTL, Target %in% Syntenic )
Net_teQTL <- subset(Net_teQTL, source %in% Syntenic)
Net_teQTL

# 
Net_teQTL <- left_join(Net_teQTL, saf[,1:2], by=c('source'="V1")) 
Net_teQTL <- left_join(Net_teQTL, saf[,1:2], by=c('Target'="V1")) 
colnames(Net_teQTL)[4:5] <-  c("ChrSource", "ChrTarget")

#table(Net_teQTL$ChrSource == Net_teQTL$ChrTarget)

# Count SNP per target gene
SNPs_per_target <- as.data.table(table(unique(trans_eQTL[,c(1,3)])$Target))
SNPs_per_target <- dplyr::left_join(SNPs_per_target, Classes, by=c("V1"="GeneID"))
SNPs_per_target$Class[is.na(SNPs_per_target$Class)]   <- "Other"

# Count SNP per Source gene
SNPs_per_source <- as.data.table(table(unique(trans_eQTL[,c(2,3)])$source))
SNPs_per_source <- dplyr::left_join(SNPs_per_source, Classes, by=c("V1"="GeneID"))
SNPs_per_source$Class[is.na(SNPs_per_source$Class)]   <- "Other"

######################################################################################



######################################################################################
##################                      QC plots                    ##################
######################################################################################
Net_teQTL
length(unique(Net_teQTL$source))
length(unique(Net_teQTL$Target))

# Add gene category 
Net_teQTL <- left_join(Net_teQTL, Classes, by=c("source"="GeneID"))
Net_teQTL <- left_join(Net_teQTL, Classes, by=c("Target"="GeneID"))
colnames(Net_teQTL)[c(6,7)] <- c("sourceClass", "targetClass")

Net_teQTL$sourceClass[is.na(Net_teQTL$sourceClass)]   <- "Other"
Net_teQTL$targetClass[is.na(Net_teQTL$targetClass)] <- "Other"
  
Total_other <-  unique(c(subset(Net_teQTL, sourceClass=='Other')$source,
                  subset(Net_teQTL, targetClass=='Other')$Target))



# Class interactions frequency
ClassInteractionFreq <- as.data.table(table(Net_teQTL[,6:7])) %>%
  dplyr::arrange(N) %>%
  dplyr::filter(N>0)


# histogram with Peaks by TF
Degree_teQTL <- as_tibble(as.data.frame(table(Net_teQTL$source), stringsAsFactors = F))
colnames(Degree_teQTL) <- c("Source", "Targets")

# histogram with Peaks by TF
InDegree_teQTL <- as_tibble(as.data.frame(table(Net_teQTL$Target), stringsAsFactors = F))
colnames(InDegree_teQTL) <- c("Target", "Regulators")

Degree_teQTL <- left_join(Degree_teQTL, Classes, by=c("Source"="GeneID"))
InDegree_teQTL <- left_join(InDegree_teQTL, Classes, by=c("Target"="GeneID"))

Degree_teQTL$Class[is.na(Degree_teQTL$Class)] <- "Other"
InDegree_teQTL$Class[is.na(InDegree_teQTL$Class)] <- "Other"

Degree_teQTL



#### Figure S3b : teQTL degree
Degree_teQTL$Class <- factor(Degree_teQTL$Class, levels = c("TF", "CoReg", "Mediator", "kinase", "Enzyme", "Other"))

Degree_teQTL

Hisplot_degree.teQTL <- ggplot(Degree_teQTL, aes(y=Targets, x=Class, fill=Class)) +
  geom_boxplot(alpha= 0.5, notch = T) + theme_pubclean() +
  #geom_jitter(size=0.2, alpha=0.5, color='grey', width = 0.15) +
  scale_y_continuous(label=comma, expand = c(0,0),  trans = 'log2') +
  xlab("Source gene") + ylab("Number of targets") +
  scale_fill_viridis(option = "D", direction = 1, discrete = T) + 
  guides(fill=guide_legend(title="Gene source")) +
  stat_compare_means(method = "anova", label.y = 6) +
  stat_compare_means(label = "p.signif", method = "wilcox", ref.group = "TF", label.y = 7) 



InDegree_teQTL$Class <- factor(InDegree_teQTL$Class, levels = c("TF", "CoReg", "Mediator", "kinase", "Enzyme", "Other"))

Hisplot_Indegree.teQTL <- ggplot(InDegree_teQTL, aes(y=Regulators, x=Class, fill=Class)) +
  geom_boxplot(alpha= 0.5,  notch = T) + theme_pubclean() +
  #geom_jitter(size=0.2, alpha=0.5, color='grey', width = 0.15) +
  scale_y_continuous(label=comma, expand = c(0,0),  trans = 'log2') +
  xlab("Target gene") + ylab("Number of regulators") +
  scale_fill_viridis(option = "D", direction = 1, discrete = T) + 
  guides(fill=guide_legend(title="Gene target")) +
  stat_compare_means(method = "anova", label.y = 6) +
  stat_compare_means(label = "p.signif", method = "wilcox", ref.group = "TF", label.y = 7) 

Hisplot_degree.teQTL
Hisplot_Indegree.teQTL

## Figure S3b part2: bar plot with fraction of data
temFre_tQTL <- as.data.frame(table(Degree_teQTL$Class))
## Figure S3b part2: bar plot with fraction of data
IntemFre_tQTL <- as.data.frame(table(InDegree_teQTL$Class))

Hisplot_degree.teQTL_bar <- ggplot(temFre_tQTL, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(position="fill", stat="identity", alpha= 0.5) +
  geom_label_repel(aes(label=scales::comma(Freq)), 
                   position = position_fill(vjust = 0.5), 
                   alpha=0.5,
                   size=1) +
  theme_pubclean()  +
  scale_y_continuous(expand = c(0,0)) +
  xlab("Source\ngene") +
  ylab("trans-eQTL fraction") + 
  scale_fill_viridis(option = "D", direction = 1, discrete = T) + 
  guides(fill=guide_legend(title="Gene source"))

Hisplot_Indegree.teQTL_bar <- ggplot(IntemFre_tQTL, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(position="fill", stat="identity", alpha= 0.5) +
  geom_label_repel(aes(label=scales::comma(Freq)), 
                   position = position_fill(vjust = 0.5), 
                   alpha=0.5, 
                   size=1) +
  theme_pubclean()  +
  scale_y_continuous(expand = c(0,0)) +
  xlab("Target\ngene") +
  ylab("trans-eQTL fraction") + 
  scale_fill_viridis(option = "D", direction = 1, discrete = T) + 
  guides(fill=guide_legend(title="Gene Target"))

Hisplot_degree.teQTL_bar
Hisplot_Indegree.teQTL_bar



######################################################################################


# Figure S3D
ClassInteractionFreq %>%
  mutate(Index=paste0(sourceClass, " -> ", targetClass)) -> ClassInteractionFreq

ClassInteractionFreq$Index <- factor(ClassInteractionFreq$Index, levels = rev(ClassInteractionFreq$Index))

ggplot(ClassInteractionFreq, aes(x=Index, y=N)) +
  geom_bar(stat = "identity") +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  annotation_logticks(sides = 'l') +
  theme_pubclean() +
  theme(axis.text.x = element_text(size=10, angle = 50, vjust = 1, hjust = 1),
        axis.text=element_text(size=10),
        text = element_text(size=10)) +
  ylab(bquote(Log[10] ~ "Counts")) +
  xlab("Interaction classes") -> Fig_3d
  
  

##### merge figures ##### 
Hisplot_degree.teQTL <- ggpar(Hisplot_degree.teQTL, 
                              font.tickslab=10, 
                              font.x = 10, 
                              font.y = 10, 
                              xtickslab.rt = 45)

Hisplot_degree.teQTL_bar <- ggpar(Hisplot_degree.teQTL_bar, 
                                  font.tickslab=10, 
                                  font.x = 10, 
                                  font.y = 10)
#
Hisplot_Indegree.teQTL <- ggpar(Hisplot_Indegree.teQTL, 
                                font.tickslab=10,
                                font.x = 10, font.y = 10, xtickslab.rt = 45)

Hisplot_Indegree.teQTL_bar <- ggpar(Hisplot_Indegree.teQTL_bar, 
                                    font.tickslab=10, 
                                    font.x = 10, font.y = 10)

Fig_S3b <- ggarrange(Hisplot_degree.teQTL, Hisplot_degree.teQTL_bar, ncol = 2, widths = c(3,1.5), legend = 'none', align = 'h')
Fig_S3c <- ggarrange(Hisplot_Indegree.teQTL, Hisplot_Indegree.teQTL_bar, ncol = 2, widths = c(3,1.4), legend = 'none', align = 'h')

# size 4x12
Fig_S3bc <- ggarrange(Fig_S3b, Fig_S3c, ncol = 2, nrow = 1)

Fig_S3bcd <- ggarrange(Fig_S3bc, Fig_3d, ncol = 1, legend = 'none', align = 'v')

pdf("Plots/Fig_S3bcd.pdf", width=8, height=7)
print(Fig_S3bcd)
dev.off()

tiff("Plots/Fig_S3d.tiff", units="in", width=10.5, height=3, res=300)
print(Fig_3d)
dev.off()


######################################################################################

######################################################################################
######## Reviewer responses: Different Number of SNPs by Gene class             #######
######################################################################################


# Total SNPs per class
SNPs_per_class <- rbind(SNPs_per_target, SNPs_per_source)

SNPs_per_class %>%
  dplyr::group_by(Class) %>%
  dplyr::summarise(Total_SNPS=sum(N)) -> SNPs_per_class_summary

# Total Genes per class
Classes_Freq <- rbind(Classes_Freq, dplyr::tibble(V1='Other', N=length(Total_other)))

SNPs_per_class_summary <- dplyr::left_join(SNPs_per_class_summary, Classes_Freq, by=c('Class'='V1'))
colnames(SNPs_per_class_summary)[3] <- "Total_genes"

# Total target per class
total_targets <- Degree_teQTL %>%
  dplyr::group_by(Class) %>%
  dplyr::summarise(Total_Targets=sum(Targets))

# Add total targets to summary
SNPs_per_class_summary <- dplyr::left_join(SNPs_per_class_summary, total_targets, by=c("Class"))

# Calculate normalized targets per SNP-gene pair
SNPs_per_class_summary <- SNPs_per_class_summary %>%
  dplyr::mutate(
    Normalized_Targets = Total_Targets / (as.numeric(Total_SNPS) * as.numeric(Total_genes)),
  )

ggplot(SNPs_per_class_summary, aes(x = reorder(Class, -Normalized_Targets), y = Normalized_Targets, 
                 size = Total_Targets, color = Class)) +
  geom_point(alpha = 0.7) +
  scale_size_continuous(range = c(3, 15)) +  # Adjust bubble size
  labs(
    title = "Trans-eQTL Target Density",
    x = "Class", 
    y = "Normalized Targets (per SNP×Gene)",
    size = "Total Targets (raw)"
  ) +
  theme_minimal()


# add Total SNPs and genes per class
library(dplyr)
Degree_teQTL_corrected <- Degree_teQTL %>%
  dplyr::left_join(
    SNPs_per_class_summary %>%
      select(Class, Normalized_Targets, Total_SNPS, Total_genes), by='Class') %>%
  # Calculate normalized target score per gene (adjust for class-specific SNP/gene counts)
  mutate(
    # Same as class-level normalization
    Normalized_Target_Score = Targets / (as.numeric(Total_SNPS) * as.numeric(Total_genes))
  )
  
Degree_teQTL_corrected <- Degree_teQTL_corrected %>%
  mutate(Normalized_Target_Score_scaled = Normalized_Target_Score * 1e9)  # Scale to "per billion"

# Update plot (change y-axis label accordingly)
Hisplot_degree.teQTL_response + 
  scale_y_log10(labels = scales::scientific, breaks = 10^(-10:0)) +
  labs(y = "Normalized Targets (per billion SNP-gene pairs)")

# 
Degree_teQTL_corrected$Class <- factor(Degree_teQTL_corrected$Class, 
                                       levels = c("TF", "CoReg", "Mediator", "kinase", "Enzyme", "Other"))

Hisplot_degree.teQTL_response <- ggplot(Degree_teQTL_corrected, 
                                        aes(y=Normalized_Target_Score, x=Class, fill=Class)) +
  geom_jitter(size=0.2, alpha=0.5, color='grey', width = 0.15) +
  geom_boxplot(outlier.shape = NA, alpha= 0.5, notch = T) + 
  theme_pubclean() +
  scale_y_log10(labels = scales::scientific, breaks = 10^(-10:0)) +  # Log scale for small normalized values
  labs(
    title = "Normalized Trans-eQTL Targets per Gene by Class",
    x = "Class",
    y = "Normalized Targets\n(per billion SNP-gene pairs)",
    caption = "Normalized by: Targets / (Total_SNPS × Total_genes) for each class"
  ) + 
  scale_fill_viridis(option = "D", direction = 1, discrete = T) + 
  stat_compare_means(method = "anova") +
  stat_compare_means(label = "p.signif", method = "wilcox", ref.group = "TF") 
  

Hisplot_degree.teQTL_response
Hisplot_degree.teQTL_response

Hisplot_degree.teQTL_response
Hisplot_Indegree.teQTL_reponse <- ggplot(InDegree_teQTL, aes(y=Regulators, x=Class, fill=Class)) +
  geom_boxplot(alpha= 0.5,  notch = T) + theme_pubclean() +
  #geom_jitter(size=0.2, alpha=0.5, color='grey', width = 0.15) +
  scale_y_continuous(label=comma, expand = c(0,0),  trans = 'log2') +
  xlab("Target gene") + ylab("Number of regulators") +
  scale_fill_viridis(option = "D", direction = 1, discrete = T) + 
  guides(fill=guide_legend(title="Gene target")) +
  stat_compare_means(method = "anova", label.y = 6) +
  stat_compare_means(label = "p.signif", method = "wilcox", ref.group = "TF", label.y = 7) 


SNPs_per_target_plot
######################################################################################


######################################################################################
#########                         Summary values                             #########
######################################################################################

mean(Degree_teQTL$Targets)

Degree_teQTL %>%
  group_by(Class) %>% summarise(Value = mean(Targets), sd=sd(Targets)) %>%
  dplyr::arrange(Value)

InDegree_teQTL %>%
  group_by(Class) %>% summarise(Value = mean(Regulators), sd=sd(Regulators)) %>%
  dplyr::arrange(Value)


Degree_teQTL %>%
  group_by(Class) %>% summarise(Value = max(Targets))

Degree_teQTL %>%
  group_by(Class) %>% 
  arrange(desc(Targets)) %>% 
  filter(Class == 'TF') %>%
  slice(1:12)


Degree_teQTL %>%
  summarise(Degree=mean(Targets))

######################################################################################

InDegree_teQTL %>%
  dplyr::filter(Class=='TF') %>%
  dplyr::arrange(-Regulators)

ReplaceName("Zm00001d050781")
ARF14_teQTL <- subset(Net_teQTL, source=='Zm00001d050781' | Target == "Zm00001d050781")


########################################################
########       HSF20 network exploration        ########
########################################################

write.table(Classes, "Table_S4.txt", sep = '\t', quote = F, row.names = F)


HSF20_teQTL <- subset(Net_teQTL, Target == "Zm00001d026094")
HSF20_teQTL <- subset(Net_teQTL, source=='Zm00001d026094' | Target == "Zm00001d026094")

write.table(HSF20_teQTL, "HSF20_teQTL.txt", sep = "\t", row.names = F, quote = F)


HeatResponses <- as_tibble(read.table("Data/DEGs/Zhou2021_DEG.tsv", sep = "\t", stringsAsFactors = F, h=T))
HeatResponses <- subset(HeatResponses[grepl("Heat", HeatResponses$cond),], Genotype=="B73")[,c(2,3,5)]


HeatResponses <- rbind(tibble(Treatment="Heat_1h_up", DEGs=unlist(strsplit(HeatResponses[1,]$gid, ","))),
                       tibble(Treatment="Heat_35h_up", DEGs=unlist(strsplit(HeatResponses[2,]$gid, ","))),
                       tibble(Treatment="Heat_2h_down", DEGs=unlist(strsplit(HeatResponses[3,]$gid, ","))),
                       tibble(Treatment="Heat_35h_down", DEGs=unlist(strsplit(HeatResponses[4,]$gid, ","))))
       
HeatResponses[, "teQTL"] <- HeatResponses$DEGs %in% subset(HSF20_teQTL, Class=='TF')$Target
HeatResponses[, "teQTL"] <- HeatResponses$DEGs %in% HSF20_teQTL$Target
table(HeatResponses[,c(1,3)])



temInterestlist <- c("Zm00001d046170", 
                     subset(net_trans_eQTL, source=="Zm00001d046170")$source,
                     subset(net_trans_eQTL, source=="Zm00001d046170")$Target)


subset(net_trans_eQTL, Target=="Zm00001d047671")
subset(net_trans_eQTLp, Target=="Zm00001d047671")
########################################################

######################################################################################
#########                         Enrichment analysis                        #########
######################################################################################

# Make CornCyc list and remove small PWY == 1 
CornCYC  <- subset(CornCYC, GeneID %in% Syntenic)
CornCYC.list <- split(CornCYC$GeneID, CornCYC$Pathway.id)

CornCYC_size <-  as.data.frame(t(as.data.frame(lapply(CornCYC.list, length))))
colnames(CornCYC_size) <- "Freq"

#
CornCYCred  <- subset(CornCYC, !(Pathway.id %in% row.names(subset(CornCYC_size, Freq == 1))))
CornCYC.list <- split(CornCYCred$GeneID, CornCYCred$Pathway.id)


Enrichmet_classes <- function(network){
  # set default names
  colnames(network)[1:2] <- c('TF',"Target")
  network <- network[,c('TF',"Target")]
  
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
  
  # Select significant TFs 
  Pval_table <- subset(Pval_table, Pval <= 0.05)
  #Pval_table <- tibble(TF=c("test", "a", "b"), TF2="test2")
  
  #write.table(Pval_table, name, sep = "\t", row.names = F, quote = F)
  print(paste("... Done ...", sep = ""))
  return(Pval_table)
}


# Identify candidates
TF_Candidates <- Enrichmet_classes(Net_teQTL)
TF_Candidates$PWY <- as.character(TF_Candidates$PWY)


##########################################################################################
#####################                Write results                    ####################
##########################################################################################

write.table(TF_Candidates, "teQTL_NetworkFinal.CornCYC.04_18_2021.txt", sep = '\t', row.names = F, quote = F)

write.table(Net_teQTL[,1:3], "teQTL_NetworkFinal.10_11_2022.txt", sep = '\t', row.names = F, quote = F)



######################################################
# Get_L2_degree <- function(df_net, nodeClass){
#   colnames(nodeClass)[1:2] <- c('GeneID', "Degree")
#   
#   Classes <- unique(nodeClass$Class)
#   
#   #out
#   total <- as_tibble(as.data.frame(matrix(0, nrow = 0, ncol = 3)))
#   colnames(total) <- c("Class", "Degree", "L2_class")
#   
#   for (c in Classes){
#     
#     # Nodes by class
#     L2_nodes <- unique(subset(nodeClass, Class==c)$GeneID)
#     
#     
#     # Subset targets of nodes in class 'c'
#     L2genes <- unique(subset(df_net, source %in% L2_nodes)$Target)
#     
#     
#     # get degree associated with L2 nodes 
#     subnodedegree <- subset(nodeClass, GeneID %in% L2genes)
#     print(subnodedegree)
#     
#     # out <- subnodedegree %>%
#     #   group_by(Class) %>% summarise(Value = mean(Degree))
#     # 
#     # out[,"L2_class"] <- c
#     subnodedegree[,"L2_class"] <- c
#     
#     total <- rbind(total, subnodedegree[,c("Class", "Degree", "L2_class")])
#   } 
#   total <- total[order(total$Class),]
#   #
#   #colnames(total)[2] <- "Degree"
#   
#   return(total)
# }
# 
# L2_OutDegree <- Get_L2_degree(Net_teQTL_full, Degree_teQTL)
# 
# 
# ggplot(subset(L2_OutDegree, Degree<=100), aes(x=Class, y=Degree, fill=L2_class))+
#   geom_boxplot() +
#   #geom_jitter()+
#   #scale_y_continuous(label=comma, expand = c(0,0),  trans = 'log2') +
#   theme_pubclean()
######################################################




################################################################################################################

TFdic <- as_tibble(read.table("Data/Annotations/TF_Id_Name.txt", h=T, stringsAsFactors = F))

TF_Candidates[,"TFname"] <- ReplaceName(TF_Candidates$TF)
TF_Candidates[,"PWYname"] <- ReplaceNamePWY(TF_Candidates$PWY)

subset(TF_Candidates, TF=='Zm00001d028842')
subset(Net_teQTL_full, source=="Zm00001d028842")

subset(PheGenes, GeneID %in% subset(Net_teQTL_full, source=="Zm00001d028842")$Target)




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



