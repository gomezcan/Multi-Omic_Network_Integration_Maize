library(scales)
library(tidyverse)
library(ggrepel)
library(ggpubr)
library(data.table)
library(GeneOverlap)
library(viridis)
library(ggrepel)
library(ComplexHeatmap)



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

ReplaceNamePWY <- function(ids){
  
  for (i in 1:nrow(CornCYC)){
    ids <- gsub(CornCYC$Pathway.id[i], CornCYC$Pathway.name[i], ids)
    #ids <- gsub("_", " ", ids)
  }
  return(ids)
}


Enrichmet_classes <- function(network, name){
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
  #Pval[,1:ncol(Pval)] <- apply(Pval[,1:ncol(Pval)], 2, p.adjust)
  
  Pval_table   <- as_tibble(reshape2::melt(as.matrix(Pval)))
  colnames(Pval_table) <- c('TF', 'class', 'Pval')
  
  #
  Common_table <- as_tibble(reshape2::melt(as.matrix(Common)))
  colnames(Common_table) <- c('TF', 'class', 'n.targ')
  
  # Add predicted target in class by TF
  Pval_table <- left_join(Pval_table, Common_table , by=c('TF', 'class'))
  
  # Add total predicted targets
  Pval_table <- left_join(Pval_table, Total_targtes, by="TF")
  
  # Select significant TFs 
  Pval_table <- subset(Pval_table, Pval <= 0.05)
  #Pval_table <- tibble(TF=c("test", "a", "b"), TF2="test2")
  
  write.table(Pval_table, name, sep = "\t", row.names = F, quote = F)
  print(paste("... Done: ", name, " ...", sep = ""))
  #return(Pval_table)
}

plot_cluster <- function(data, var_cluster, palette, labels) {
  #
  Genelabel <- subset(data, GeneId %in% GenesTarget)
  clusters_target <- unique(as.character(subset(data, GeneId %in% GenesTarget)$km))
  #
  print(Genelabel)
  #
  ggplot(data, aes_string(x="V1", y="V2", color=var_cluster)) +
    geom_point(size=0.05) +
    geom_point(data=subset(data, km %in% clusters_target), size=0.1, aes_string(x="V1", y="V2", color=var_cluster)) + 
    geom_label_repel(data=Genelabel, aes(x=V1, y=V2, label=Name), fill='black', box.padding = 0.5) +
    
    guides(colour=guide_legend(override.aes=list(size=3))) +
    xlab("tsne 1") + ylab("tsne 2") +
    ggtitle("") +
    theme_light(base_size=20) +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          legend.direction = "horizontal", 
          legend.position = "bottom",
          legend.box = "horizontal") + 
    scale_color_manual(values = rainbow(palette)) +
    geom_label_repel(data=labels, aes(label = km)) +
    guides(colour = 'none') 
}

outdegree.counter <- function(network){
  ## Count TF targets in network
  # Count Total
  network <- subset(network, tgt.gid %in% Syntenic & reg.gid %in% Syntenic) 
  #
  Total_targtes <- as_tibble(as.data.frame(table(unique(network[,1:2])$reg.gid)))
  colnames(Total_targtes) <- c('TF', 'targets') 
  Total_targtes$TF <- as.character(Total_targtes$TF)
  
  
  return(Total_targtes)
}



##################################################
##########        Annotations       ##############
##################################################

## syntenic genes 
Syntenic <- as_tibble(read.table("Data/Annotations/Zm.v4.synteny.genes.txt", h=T, stringsAsFactors = F, sep = '\t'))
Syntenic <- subset(Syntenic, gene_synteny=='syntenic')$gene_id

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

# CornCYC
CornCYC <- as_tibble(read.table("Data/Annotations/corn_pathways.0210325.reduced.txt", h=T, stringsAsFactors = F))
CornCYC <- CornCYC[!(CornCYC$Gene.id %in% "unknown"),]
CornCYC <- subset(CornCYC, Gene.id %in% Syntenic)
CornCYC[,"Class"] <- "Enzyme"
colnames(CornCYC)[3] <- "GeneID"
CornCYC$Pathway.id <- gsub("-", "_", CornCYC$Pathway.id)


# kinases
kinases <- as_tibble(read.table("Data/Annotations/kinases_maize_AGPv4.txt", h=T, stringsAsFactors = F))
kinases <- subset(kinases, GeneID %in% Syntenic)
kinases[,"Class"] <- "kinase"


# Coexpression data
# CoexDB <- as_tibble(read.table("ExpressionData/widiv304.grn.tsv.gz", h=T))
# CoexDB <- CoexDB[order(-CoexDB$score),][1:100000,]

CoexDB = readRDS("Data/Coexpression/rf.100k.rds")
# add network names
names(CoexDB$tn) <- CoexDB$nid
names(CoexDB)
CoexDB <- CoexDB$tn

# Coexpression results from 304 lines used on eQTL analysis
CoexDB_304 <- as_tibble(read.table("Data/Coexpression/widiv304.grn.tsv.gz", h=T))
CoexDB_304 <- CoexDB_304[order(-CoexDB_304$score),]
CoexDB_304 <- CoexDB_304[1:100000,] # select top predictions

# add 
CoexDB[['n304']] <- CoexDB_304

# Collapse network
CoexDB_reduce <- unique(as_tibble(rbindlist(CoexDB, idcol = FALSE))[,1:2])

# 
CoexDB_reduce <- subset(CoexDB_reduce, tgt.gid %in% Syntenic & reg.gid %in% Syntenic)
CoexDB_reduce[,1:2] <- apply(CoexDB_reduce[,1:2], 2, as.character)

dim(CoexDB_reduce)
length(unique(CoexDB_reduce$tgt.gid))
length(unique(CoexDB_reduce$tgt.gid))

# # define query list
# List_query <- rbind(tibble(Target="Phe", GeneID=unique(PheGenes$GeneID)),
#                     tibble(Target="TFs", GeneID=unique(TF_CoR$GeneID)), 
#                     tibble(Target="TF.phe", GeneID=unique(Y1H$TF.v4)[is.na(unique(Y1H$TF.v4)) == FALSE]))
# 
# List_query <- split(List_query$GeneID, List_query$Target)
# 
# #
# List_query <- c(List_query, CornCYC.list)


######################################################################################
#########         Count targets by TF, and general description               #########
######################################################################################

Total.OutDegree <- as_tibble(as.data.frame(table(CoexDB_reduce$reg.gid), stringsAsFactors = F))
colnames(Total.OutDegree) <- c("TF","Total.targets")
mean(Total.OutDegree$Total.targets)
sd(Total.OutDegree$Total.targets)


# Count targets by TF and by network
OutDegree <- lapply(CoexDB, outdegree.counter)
OutDegree <- unique(as_tibble(rbindlist(OutDegree, idcol = T)))
colnames(OutDegree)[1] <- "Network"
OutDegree
# Number of TF-target association by net
Net.size <- OutDegree %>%
  group_by(Network) %>%
  summarise(Size =  sum(targets)) %>%
  arrange(-Size)
Net.size

# Sub top 2 TF by network
Top.OutDegree <- OutDegree %>%
  arrange(desc(targets)) %>%
  group_by(Network) %>% slice(1:2)

Top.OutDegree
  
# average targets by network
Mean.OutDegree <- OutDegree %>%
  group_by(Network) %>% 
  summarise(Average =  mean(targets)) 


# Number of TFs by network
TFs.by.network <- as_tibble(as.data.frame(table(OutDegree$Network), stringsAsFactors = F))
colnames(TFs.by.network) <- c("Network", "TFs")


# Add total targets
OutDegree <- left_join(OutDegree, Total.OutDegree, by='TF')

# add network order
OutDegree$Network <- factor(OutDegree$Network, levels = Net.size$Network)
TFs.by.network$Network <- factor(OutDegree$Network, levels = Net.size$Network)

# add name
Top.OutDegree[,"Name"] <- Top.OutDegree$TF
Top.OutDegree$Name <- ReplaceName(Top.OutDegree$Name)

# Plot S1a TFs by network
mean(TFs.by.network$TFs)

Plot.TFs.by.net <- ggplot(TFs.by.network, aes(x=TFs)) +
  geom_histogram(binwidth = 50, fill='black', alpha=0.8) + 
  geom_vline(aes(xintercept=mean(TFs)), color='grey', linetype='dashed', size=1.5) +
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(label=comma) + 
  theme_bw() 


Plot.TFs.by.net  <- ggpar(Plot.TFs.by.net, font.tickslab=14, font.x = 14, font.y = 14,
                          xlab = 'TFs by Network', ylab="Frequence")

mean(TFs.by.network$TFs)

# Notes from plot S1b
mean(Mean.OutDegree$Average) # to describe targets by TF in all netwtoks

Top.TF.hubs <- as_tibble(as.data.frame(table(Top.OutDegree$TF)))
Top.TF.hubs <- Top.TF.hubs[order(-Top.TF.hubs$Freq),]
Top.TF.hubs[,"Name"] <- ReplaceName(Top.TF.hubs$Var1)


# Plot S1b, targets by TF by network
Plot.OutDegree <- ggplot(OutDegree, aes(x=targets, y=Network)) +
  geom_boxplot(alpha=0.5) +
  geom_text_repel(data = subset(Top.OutDegree, TF %in% as.character(Top.TF.hubs$Var1[1:7])),
                  aes(x=targets, Network, label=Name),
                  nudge_x = .15,
                  box.padding = 0.5,
                  nudge_y = 1,
                  segment.curvature = -0.1,
                  segment.ncp = 3,
                  segment.angle = 20,
                  size=2, 
                  color="Orange") +
  theme_bw() +
  scale_x_continuous(label=comma) + 
  theme(legend.title = element_text(vjust = 0.8), 
        panel.border = element_blank(), 
        panel.grid.major = element_blank())

Plot.OutDegree

# Plot S1c total TF-target degree
Total.OutDegree <- Total.OutDegree[order(-Total.OutDegree$Total.targets),]
mean(Total.OutDegree$Total.targets)

Total.OutDegree
length(unique(Total.OutDegree$TF))




Plot.Total.OutDegree <- ggplot(Total.OutDegree, aes(x=Total.targets)) +
  geom_histogram(binwidth = 50, fill='orange1', alpha=0.5) + 
  geom_vline(aes(xintercept=mean(Total.targets)), color='black', linetype='dashed', size=1) +
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0), label=comma) + 
  theme_bw() 


Plot.Total.OutDegree  <- ggpar(Plot.Total.OutDegree, 
                               font.tickslab=14, font.x = 14, font.y = 14,
                               xlab = 'Total targets by TF', ylab="Frequence")


# Figure S1
Plot.OutDegree <- ggpar(Plot.OutDegree, 
                        font.xtickslab=14, 
                        font.ytickslab=10,
                        font.x = 14, font.y = 14)
Plot.OutDegree
# size: 7x7
Fig_S1 <- ggarrange(ggarrange(Plot.TFs.by.net,Plot.Total.OutDegree, ncol = 1),
                    Plot.OutDegree, ncol = 2, widths = c(0.7, 1))
Fig_S1

pdf('Plots/Fig_S1.pdf', width = 8, height = 8)
print(Fig_S1)
dev.off()

Total.OutDegree
CoexDB_reduce
mean(Total.OutDegree$Total.targets)
length(unique(CoexDB_reduce$reg.gid))
length(unique(CoexDB_reduce$tgt.gid))
######################################################################################

######################################################################################
#########                          Enrichment test                           #########
######################################################################################

length(unique(CoexDB_reduce$tgt.gid))

GenesInNet <- tibble(GeneID=unique(CoexDB_reduce$tgt.gid), Class="Other")

Classes <- rbind(unique(TF_CoR[,c("GeneID", "Class")]), 
                 unique(Mediator[,c("GeneID", "Class")]), 
                 unique(kinases[,c("GeneID", "Class")]),
                 unique(CornCYC[,c("GeneID", "Class")]))

Classes <- Classes[(Classes$GeneID %in% Syntenic),]

# add all genes 
Classes <- rbind(Classes, GenesInNet[!(GenesInNet$GeneID %in% Classes$GeneID),])
Classes_list <- split(Classes$GeneID, Classes$Class)

# General numbers
CornCYC
length(unique(CornCYC$Pathway.id))
length(unique(CornCYC$GeneID))

# make list input to "Enrichmet_classes" function
CornCYC.list <- split(CornCYC$GeneID, CornCYC$Pathway.id)

# test CornCYC pathways
mapply(Enrichmet_classes, CoexDB[-c(1:12)], paste0("Candidates_files/CornCYC.", names(CoexDB)[-c(1:12)], ".txt"))




## set parameters to reads enrichment test results
Enrichfile <- list.files("Candidates_files", pattern="*.txt", full.names=TRUE)
Enrichfilename <- gsub("Candidates_files/CornCYC.", "", Enrichfile)
Enrichfilename <- gsub(".txt", "", Enrichfilename)

# read resutls
TFs_by_CornCYC <- lapply(Enrichfile, fread)
names(TFs_by_CornCYC) <- Enrichfilename

TFs_by_CornCYC <- as_tibble(rbindlist(TFs_by_CornCYC, idcol = TRUE))


# Count significant TF-pathway associations 
TFs_by_Classes_summary <- as.data.frame(table(subset(TFs_by_CornCYC, Pval <= 0.05)[,2:3]), stringsAsFactors = F)


CytoNet_TF_CornCYC_Freq <- subset(TFs_by_Classes_summary, Freq>0)
CytoNet_TF_CornCYC_Freq <- left_join(CytoNet_TF_CornCYC_Freq, unique(CornCYC[,1:2]), by=c("class"="Pathway.id"))
dim(CytoNet_TF_CornCYC_Freq)

length(unique(CytoNet_TF_CornCYC_Freq$TF))
length(unique(CytoNet_TF_CornCYC_Freq$Pathway.name))

####################################
###### write main results 
####################################
FinalNet <- unique(as_tibble(rbindlist(CoexDB, idcol = TRUE))[,1:3])
colnames(FinalNet) <- c("Net", "TF", "Target")
write.table(FinalNet, "CoExp_NetworkFinal.10_11_2021.txt", sep = "\t", row.names = F, quote = F) 
head(CytoNet_TF_CornCYC_Freq)
write.table(CytoNet_TF_CornCYC_Freq, "CoExp_NetworkFinal.CornCYC.04_18_2022.txt", sep = "\t", row.names = F, quote = F) 
write.table(TFs_by_CornCYC, "CoExp_NetworkFinal.Full.CornCYC.04_18_2022.txt", sep = "\t", row.names = F, quote = F)

subset(CytoNet_TF_CornCYC_Freq, TF=="Zm00001d028850")PWY_7897

subset(TFs_by_CornCYC, class=="PWY_7897" & n.targ >= 8)$TF
ReplaceName(subset(TFs_by_CornCYC, class=="PWY_7897" & n.targ >= 8)$TF)

ReplaceNamePWY("PWY_7897")

length(subset(CornCYC, Pathway.id=="PWY_7897")$GeneID)

CytoNet_TF_CornCYC_Freq[1:5,]

####################################


# from summary to matrix
TFs_by_ClassesM <- reshape2::dcast(TFs_by_Classes_summary, TF ~ class)
row.names(TFs_by_ClassesM) <- TFs_by_ClassesM$TF
TFs_by_ClassesM <- TFs_by_ClassesM[,-c(1)]
#TFs_by_ClassesM["Zm00001d028854",]
# 
table(colSums(TFs_by_ClassesM_r) > 1)
table(rowSums(TFs_by_ClassesM) > 1)


##########
## Rtsne
##########
library(Rtsne)
library(factoextra) # clustering algorithms & visualization
dim(TFs_by_ClassesM)
TFs_by_ClassesM[1:5,1:5]
tsne_table2 = Rtsne(as.matrix(TFs_by_ClassesM), check_duplicates=FALSE, pca=TRUE, perplexity=30, theta=0.5, dims=2)
d_tsne_2 = as.data.frame(tsne_table2$Y)
dim(d_tsne_2)
row.names(d_tsne_2) <- row.names(TFs_by_ClassesM)

fviz_nbclust(d_tsne_2, kmeans,nboot = 50, k.max = 60) # best opt: 58 k

## Creating k-means clustering model, and assigning the result to the data used to create the tsne

fit_cluster_kmeans2 = kmeans(d_tsne_2, 33) # 
d_tsne_2$km <-  factor(fit_cluster_kmeans2$cluster)
d_tsne_2[,"GeneId"] <- row.names(d_tsne_2)
d_tsne_2 <- as_tibble(d_tsne_2)
d_tsne_2[,"Name"] <- ReplaceName(d_tsne_2$GeneId)

View(d_tsne_2[d_tsne_2$GeneId %in% TFdic$TF.v4,])


#ReplaceName(d_tsne_2$GeneId) <- ReplaceName(d_tsne_2$GeneId)

## Exploratory plot

# labels
clusterLabels2 <- d_tsne_2 %>% 
  group_by(km) %>% select(V1, V2) %>% summarize_all(mean)


GenesTarget = c("Zm00001d028854","Zm00001d047671", "Zm00001d043036")
table(d_tsne_2$GeneId %in% GenesTarget)

##########
d_tsne_2
plot_tsne <- plot_cluster(d_tsne_2, "km", 33, clusterLabels2) 

#########
# mask for  MYB154: cellulose_biosynthesis 
#########
unique(subset(CornCYC, Pathway.id=="PWY_7120")[,1:2])
unique(subset(CornCYC, Pathway.id=="PWY_5080")[,1:2])

#unique(subset(TFs_by_Classes_summary, class=="PWY_1001" & Freq>0)$TF)
TFs_by_ClassesM_r <- TFs_by_ClassesM[subset(d_tsne_2, km==9)$GeneId, ]
#TFs_by_ClassesM_r <- TFs_by_ClassesM_r[rowSums(TFs_by_ClassesM_r) > 5, ]
TFs_by_ClassesM_r <- TFs_by_ClassesM_r[,colSums(TFs_by_ClassesM_r) > 0]
dim(TFs_by_ClassesM_r)

row.names(TFs_by_ClassesM_r) <- ReplaceName(row.names(TFs_by_ClassesM_r))
colnames(TFs_by_ClassesM_r) <- ReplaceNamePWY(colnames(TFs_by_ClassesM_r))

Hetmap_Matrix <- Heatmap(t(TFs_by_ClassesM_r), name="Freq",
                         # column_km = 3,
                         column_names_rot = 45,
                         # row_names_rot = 45,
                         cluster_rows = TRUE, cluster_columns = TRUE,
                         show_column_dend = FALSE, show_row_dend = FALSE, 
                         #clustering_method_columns = "ward.D2",
                         #clustering_method_rows = "ward.D2",
                         col=viridis(4, direction = -1, option = "B"),
                         column_names_gp = gpar(fontsize = 7),
                         row_names_gp = gpar(fontsize = 6),
                         show_heatmap_legend = T,
                         # heatmap_ = unit(10),
                         heatmap_height = unit(12, 'cm'),
                         heatmap_width  = unit(12, 'cm'),
                         heatmap_legend_param = list(title_gp = gpar(fontsize = 10), 
                                                     labels_gp = gpar(fontsize = 10), 
                                                     direction = "vertical"))

#5x6
Hetmap_Matrix
#########

#########
Enrichmet.targets <- function(network){
  ## Count TF targets in network
  # Count Total
  network <- subset(network, tgt.gid %in% Syntenic & reg.gid %in% Syntenic) 
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
  
  ## Compare list of predicted targets vs annotated genes in query.vector 
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

#####################################################

#############################################################
#######     Define Unique cis-QTLs association     ########
#############################################################

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



