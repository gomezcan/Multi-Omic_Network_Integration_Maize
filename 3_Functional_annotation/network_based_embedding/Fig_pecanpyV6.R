library(parallel)
suppressMessages(library(Rgraphviz))
library(tm)
library(SnowballC)
library(wordcloud)
library(RColorBrewer)
library(topGO)
library(GeneOverlap)
library(tidyverse)
library(ggrepel)
library(ggpubr)
library(data.table)
library(viridis)
library(ComplexHeatmap)
library(fgsea)
library(reshape2)
library(circlize)
library(ggVennDiagram)
library(scales)
library(purrr)
library(gplots)
library(ggplot2)
library(PCAtools)
library(factoextra)
#library(networkD3)
library(dplyr)
library(patchwork)



##################################################
##########          Functions        #############
##################################################

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

GetGO <- function(degs, mutant, netname){
  
  # Use a list of DEGs and the name of the mutant (string)
  # to identify GOs enriched. Required to have a background predefined
  # Define background based on genes in network
  
  background_tem <- background[names(background) %in% unique(c(Net$Source, Net$Target))]
  background_IDs_tem <- as.character(unique(names(background_tem)))
  
  # mutant: description of a category of source of the data
  # degs: set of genes to test enrichment
  GeneList <- factor(as.integer(background_IDs_tem %in% degs))
  names(GeneList) <- background_IDs_tem
  
  
  GOdata_BP <- new("topGOdata", ontology = "BP", allGenes = GeneList, 
                   annot = annFUN.gene2GO, gene2GO = background_tem)
  
  #GOdata_MF <- new("topGOdata", ontology = "MF", allGenes = GeneList, 
  # annot = annFUN.gene2GO, gene2GO = background)
  #GOdata_CC <- new("topGOdata", ontology = "CC", allGenes = GeneList, 
  # annot = annFUN.gene2GO, gene2GO = background)
  
  #### Define test ####
  test.stat <- new("classicCount", testStatistic = GOFisherTest, 
                   name = "Fisher test", nodeSize = 10)
  
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
  
  DF_GO_Genes <- unique(DF_GO_Genes)
  
  filename <- 
    write.table(DF_GO_Genes,
                paste("BP_results_targets/Genes_GOBP_", netname, "_", mutant, ".txt", sep = ""),
                sep = '\t', quote = F,
                row.names = F)
  
  #return(list(GOs_DF=Res_DF_BP, Genes=DF_GO_Genes))# return list of GOs-Stats and GeneID-GOs
  return(Res_DF_BP) # Return list of GOs-Stats and GeneID-GOs
}

SuperGO_Modules_Targ <- function(tf, netname){
  ## used TF/Module targets/genes to test GO terms enrichment
  # 1. Select tf/module's Targets
  # 2. Make list file: degs
  # 3. Test enrichment
  
  # Get network by TF
  network <- subset(Net, Source==tf)
  network <- subset(network, Target %in% Syntenic)
  
  Total_targets <- as_tibble(as.data.frame(table(unique(network)$Source), stringsAsFactors = F))
  colnames(Total_targets) <- c('TF', 'nTF') 
  
  # Genes input Net
  degs <- unique(network$Target)
  
  ## Compare list of predicted targets vs annotated genes in query.vector 
  print(". Pre-GO.")
  out <- GetGO(degs, tf, netname)
  #
  colnames(out)[7] <- "TF"
  
  # save GOs
  out <- left_join(out, Total_targets, by="TF")  
  write.table(out, paste0("BP_results_targets/",netname, "_", tf, "_GOs.txt"), sep = '\t', row.names = F, quote = F)
  # Genes_GOBP_
  return(print('.Done.'))
  
}

toSpace <- content_transformer(function (x , pattern ) gsub(pattern, " ", x))

Get_WC <- function(GO_list) {
  
  # Read the text vecto: if file to reads
  # text <- readLines(GO_list)
  
  ## Load the data as a corpus
  # df <- data.frame(doc_id = 'GOs', text = BP_PDIs$Term, stringsAsFactors = FALSE)
  # docs <- Corpus(DataframeSource(df))
  
  docs <- Corpus(VectorSource(GO_list)) # vector soruce
  
  
  ## text transformation
  docs <- tm_map(docs, toSpace, "/")
  docs <- tm_map(docs, toSpace, "@")
  docs <- tm_map(docs, toSpace, "\\|")
  docs <- tm_map(docs, toSpace, "\\.")
  docs <- tm_map(docs, toSpace, "process")
  docs <- tm_map(docs, toSpace, "response")
  
  ## Cleaning the text
  # Convert the text to lower case
  docs <- tm_map(docs, content_transformer(tolower))
  # Remove numbers
  docs <- tm_map(docs, removeNumbers)
  # Remove english common stopwords
  docs <- tm_map(docs, removeWords, stopwords("english"))
  
  # Remove your own stop word
  # specify your stopwords as a character vector
  #docs <- tm_map(docs, removeWords, c("blabla1", "blabla2")) 
  # Remove punctuations
  docs <- tm_map(docs, removePunctuation)
  # Eliminate extra white spaces
  docs <- tm_map(docs, stripWhitespace)
  #inspect(docs)
  
  # Text stemming
  # docs <- tm_map(docs, stemDocument)
  
  ## Build a term-document matrix
  dtm <- TermDocumentMatrix(docs)
  m <- as.matrix(dtm)
  v <- sort(rowSums(m),decreasing=TRUE)
  d <- data.frame(word = names(v),freq=v)
  d <- subset(d, freq >=1)
  
  # plot
  set.seed(1234)
  wordcloud(words = d$word, freq = d$freq, min.freq = 2,
            max.words=400, random.order=FALSE, rot.per=0.40, 
            colors=brewer.pal(10, "Paired"), fixed.asp=T) 
}

chop=function(myStr,mySep,myField){
  ## chop a string by a separator and return specified field
  choppedString=sapply(strsplit(myStr,mySep),"[",myField)
  if(length(myField)>1){
    choppedString=apply(choppedString,2,function(x){paste0(x[!is.na(x)],collapse=mySep)})
  }
  return(choppedString)
}

##################################################
##########        Annotations       ##############
##################################################
## Syntenic genes 
Syntenic <- as_tibble(read.table("Data/Annotations/Zm.v4.synteny.genes.txt", h=T, stringsAsFactors = F, sep = '\t'))
Syntenic <- subset(Syntenic, gene_synteny=='syntenic')$gene_id
length(Syntenic)
#### top 45 ###
#Top45 <- as_tibble(read.table("Data/Annotations/Top45.txt", h=F, stringsAsFactors = F))

# TF names
TFdic <- as_tibble(read.table("Data/Annotations/TF_Id_Name.txt", h=T, stringsAsFactors = F))

# Phenolic related genes
# PheGenes <- as_tibble(read.table("Data/Annotations/LinaPheGenes2020.txt", h=T, sep = "\t", quote="", stringsAsFactors = F))

# TF and CoReg
TF_CoR <- as_tibble(read.table("Data/Annotations/TF_CoR_Mazie.txt", h=T, stringsAsFactors = F))

## Y1H network
#Y1H <- as_tibble(read.table("Data/Annotations/Y1H.data.txt", h=T, stringsAsFactors = F))[,4:6]
#ReplaceName(Y1H$TF.v4)

# CornCYC
CornCYC <- as_tibble(read.table("Data/Annotations/corn_pathways.0210325.reduced.txt", h=T, stringsAsFactors = F))
CornCYC <- CornCYC[!(CornCYC$Gene.id %in% "unknown"),]
CornCYC <- subset(CornCYC, Gene.id %in% Syntenic)
#CornCYC[,"Class"] <- "Enzyme"
colnames(CornCYC)[3] <- "GeneID"
CornCYC$Pathway.id <- gsub("-", "_", CornCYC$Pathway.id)

# make CornCYC list and remove small PWY == 1 
CornCYC  <- subset(CornCYC, GeneID %in% Syntenic)
CornCYC.list <- split(CornCYC$GeneID, CornCYC$Pathway.id)
#
CornCYC_size <-  as.data.frame(t(as.data.frame(lapply(CornCYC.list, length))))
colnames(CornCYC_size) <- "Freq"
#
#CornCYCred  <- subset(CornCYC, !(Pathway.id %in% row.names(subset(CornCYC_size, Freq == 1))))
#CornCYC.list <- split(CornCYCred$GeneID, CornCYCred$Pathway.id)

#CornCYCred_size <- as_tibble(as.data.frame(table(CornCYCred$Pathway.id), stringsAsFactors = F))
#colnames(CornCYCred_size) <- c("PWY", "nPWY")

# GOs term annotations
background <- readMappings("Data/Annotations/synteny.ID_TopGO_V4_GRAMER.txt")
background_IDs <- as.character(unique(names(background)))


##################################################
######        Setup pecanpy entry          #######
##################################################

# PDI
PDI <- unique(fread("../Fig_PDI/Only_PDI_NetworkFinal.10_14_2022.txt")[,c(2,3)])
colnames(PDI)[1] <- "Source"
PDI <- unique(PDI$Source)

PDIeQTL <- unique(fread("../Fig_PDI/CisE_PDI_NetworkFinal.10_14_2022.txt")[,c(2,3)])
colnames(PDIeQTL)[1] <- "Source"
PDIeQTL <- unique(PDIeQTL$Source)

# CoExp
CoExp <- unique(fread("../Fig_Coexpression/CoExp_NetworkFinal.10_11_2021.txt"))
colnames(CoExp)[2] <- "Source"
CoExp <- unique(CoExp$Source)

# Defined total tfs
All_TFs <- unique(c(PDI, PDIeQTL, CoExp, TF_CoR$GeneID))
length(All_TFs)
length(unique(TF_CoR$GeneID))
fwrite(data.table(TF=All_TFs), "All_TFs.txt", row.names = F, col.names = F, quote = F)

# teQTL
teQTL <- unique(fread("../Fig_transeQTL/teQTL_NetworkFinal.10_11_2021.txt"))
teQTL <- unique(teQTL$source)
teQTLtfs <-  teQTL[teQTL %in% All_TFs]

# Full network
FullNet <- fread("uniqFullNets_weighted.txt")
colnames(FullNet) <- c('Source', 'Target', 'Weight')
FullNet <- unique(subset(FullNet, Target %in% Syntenic))

length(unique(c(FullNet$Source, FullNet$Target)))
#FullNet <- unique(subset(FullNet, Source %in% All_TFs))
#FullNet <- unique(subset(FullNet, Source %in% Syntenic))


###################################################
# Pecanpy and MR clusters: V2                     #
# Alpha 0.05                                      #
###################################################

ClustersFull <- fread('DistanceCalculation/InputClusterONE_Dim50_WL80_nW10_0.005_syntenic.txt')
ClustersFull  <- subset(ClustersFull, V2 %in% Syntenic)
ClustersFull  <- subset(ClustersFull, V1 %in% Syntenic)

# Subset of cluster to TFs data
Clusters  <- subset(ClustersFull, V1 %in% All_TFs)

##########################################################
####            Compare Clusters summary              ####
##########################################################


# Genes by cluster
Clusters_Size <- as.data.table(table(Clusters$V1))
colnames(Clusters_Size) <- c("TF", "Genes")

#table(Clusters_Size$Class)
#mean(Clusters_Size$Genes)

# set input object
Net <- unique(Clusters[,1:2])
colnames(Net) <- c("Source", "Target")


test2 <- subset(Clusters, V1=='Zm00001d051520')$V2
test3 <- subset(FullNet, Source=='Zm00001d051520')$Target

venn(list(test2=test2,
          test3=test3))

table(test1 %in% test2)
table(test2 %in% test1)
subset(Clusters_Size, TF=="Zm00001d006236") # MYB31
subset(Clusters_Size, TF=="Zm00001d005016") # WRI1
subset(Clusters_Size, TF=="Zm00001d033859") # KN1
subset(Clusters_Size, TF=="Zm00001d051520") # MYB19

###############################################################

###############################################################
#######              Cluster annotation                 #######
###########          CornCYC Enrichment               #########
###############################################################

targList=CornCYC.list
network=Net
Enrichmet_classesV2 <- function(network, targList){
  
  ## Count TF targets in network
  # Count Total
  network <- unique(network[,c("Target", "Source")])
  #network <- subset(network, Target %in% Syntenic)
  
  # Genes in module
  Total_targets <- as.data.table(table(network$Source))
  colnames(Total_targets) <- c('TF', 'nTF') 
  
  
  # Genes by TF
  totalTarg.byTF <- as.data.frame(mapply(length, targList))
  totalTarg.byTF[,"PWY"] <- row.names(totalTarg.byTF)
  colnames(totalTarg.byTF)[1] <- "nPWY"
  #head(totalTarg.byTF)
  
  # list input: network
  network.list <- split(network$Target, network$Source)
  
  
  ## Compare list of predicted targets vs annotated genes in query.vector 
  #
  print(". Pre-newGOM .")
  
  Tem_length <- Syntenic[Syntenic %in% unique(c(network$Target, network$Source))]
  Tem_length <- length(Tem_length)
  
  go.obj <- newGOM(network.list, targList, genome.size=Tem_length) # annotated genes in Genome v4
  
  Pval <- getMatrix(go.obj, name="pval")
  Common <- getMatrix(go.obj, name="intersection")
  
  print(". Post-newGOM .")
  
  ### Summary tables
  ## adjust p value
  Pval <- as.data.frame(Pval)
  
  Pval_table   <- as_tibble(reshape2::melt(as.matrix(Pval)))
  colnames(Pval_table) <- c('TF', 'PWY', 'Pval')
  Pval_table$TF <- as.character(Pval_table$TF)
  Pval_table$PWY <- as.character(Pval_table$PWY)
  
  #
  Common_table <- as_tibble(reshape2::melt(as.matrix(Common)))
  colnames(Common_table) <- c('TF', 'PWY', 'n.targ')
  Common_table$TF  <- as.character(Common_table$TF)
  Common_table$PWY <- as.character(Common_table$PWY)
  
  # Add predicted target in class by Module
  Pval_table <- left_join(Pval_table, Common_table , by=c('TF', 'PWY'))
  
  # Add total genes in module targets
  Pval_table <- left_join(Pval_table, Total_targets, by="TF")
  
  # add total target by TF
  Pval_table <- left_join(Pval_table, totalTarg.byTF, by="PWY")
  
  #write.table(Pval_table, name, sep = "\t", row.names = F, quote = F)
  print(paste("... Done ...", sep = ""))
  return(Pval_table)
}

# set input object
# PWY enrichment
PWY_Clusters <- Enrichmet_classesV2(ClustersNet, CornCYC.list) 


# Add PWY name
PWY_Clusters  <- dplyr::left_join(PWY_Clusters, unique(CornCYC[,1:2]), by=c("PWY"="Pathway.id"))

write.table(PWY_Clusters, "PWY_GO_results/NetworkBased_PWY_Clusters_enrichment.txt", 
            sep = '\t', quote = F, row.names = F)

#PWY_Clusters <- fread("PWY_GO_results/NetworkBased_PWY_Clusters_enrichment.txt")

# Add FDR by TF
PWY_Clusters %>% 
  group_by(TF) %>%
  mutate("FDR" = p.adjust(Pval, method = 'fdr')) %>%
  filter(Pval <= 0.05) %>%
  dplyr::arrange(TF) -> PWY_Clusters

subset(PWY_Clusters, Pval <=0.01)

subset(PWY_Clusters, FDR <=0.1)

hist(PWY_Clusters$n.targ)
length()
table(test$Module)


# Examples MYB31
subset(PWY_Clusters, TF=="Zm00001d006236") # MYB31
subset(PWY_Clusters, TF=="Zm00001d005016") # WRI1
subset(PWY_Clusters, TF=="Zm00001d051520") # MYB19
subset(PWY_Clusters, TF=="Zm00001d033859") # KN1

table(PDI %in% Clusters$V1)
table(teQTLtfs %in% Clusters$V1)
table(CoExp %in% Clusters$V1)

# TFs not in networks
All_TFsNew <- All_TFs[!(All_TFs %in% c(PDI, teQTLtfs, CoExp))]
All_TFsNew <- All_TFs[!(All_TFs %in% c(PDI, CoExp))]
table(All_TFsNew %in% Clusters$V1)
table(All_TFsNew %in% PDI)


###############################################################

###############################################################
#######               Cluster annotation               #######
###########            GO enrichment                  #########
###############################################################

# Defined number of clusters
GenesList <- unique(Net$Source) 
Lgenes <- length(GenesList)
#Lgenes <- 100

# Size of range to test
w = 2

print(".. Ready to start ..")
for (i in seq(0, Lgenes, w)){
  max=Lgenes
  Start=i+1
  end=i+w
  # 
  if (end <= max){
    listtotest <- GenesList[Start:end]
    print(length(listtotest))
    mclapply(listtotest, SuperGO_Modules_Targ, mc.cores=2)
  }
  else if (end > max){ break}
  else {
    r <- paste(Start, max, sep = "-")
    listtotest <- GenesList[Start:max]
    print(length(listtotest))
    
    w = length(listtotest)+1
    mclapply(listtotest, SuperGO_Modules_Targ, mc.cores=2)
  }
}

GOs_Clusters <- list.files(path = 'BP_results/', pattern = "^MRnet_*")
# GOs_Clusters <- list.files(path = 'BP_results/', pattern = "^MI0.5_*")

# secong round for cluster not tested
ClusterDone <- gsub("Module_","", GOs_Clusters)
ClusterDone <- gsub("_GOs.txt","",ClusterDone)
GenesList <- GenesList[!(GenesList %in% ClusterDone)]
GenesList <- unique(GenesList)
Lgenes <- length(GenesList)

for (i in GenesList){
  SuperGO_Modules_Targ(i)
  
  
}

## Remove empty results
GOs_Clusters <- list.files(path = 'BP_results/', pattern = "^MRnet_*")
GOs_Clusters <- lapply(GOs_Clusters, function(x) fread(paste0("BP_results/",x)))

mask <- unlist(lapply(GOs_Clusters, function(x) is.data.frame(x)))
GOs_Clusters <- GOs_Clusters[mask]
length(GOs_Clusters)

## Combine DF results 
GOs_Clusters <- rbindlist(GOs_Clusters, idcol = F)

write.table(GOs_Clusters, 
            "PWY_GO_results/NetworkBased_GO_Clusters_enrichment.txt", 
            sep = '\t', quote = F, row.names = F)

GOs_ClustersFiltered <- GOs_Clusters %>% 
  group_by(TF) %>%
  mutate("FDR" = p.adjust(classic, method = 'fdr')) %>%
  filter(classic <= 0.05) 

length(table(subset(GOs_ClustersFiltered, FDR <=0.1)$TF))
length(table(GOs_ClustersFiltered$TF))



###############################################################

###############################################################
#######              Cluster Description                #######
###########             TFs-based                     #########
###############################################################

colnames(Clusters_Size)[3] <- 'Total.TFs'

# Total new TFs by cluster
TFs_notTested <- All_TFs[!(All_TFs %in% c(PDI, PDIeQTL, CoExp, teQTLtfs))]


# Total Genes by TF
Clusters %>% 
  dplyr::group_by(V1) %>%
  dplyr::summarise(Genes=n()) %>%
  dplyr::mutate(Class='TFs') -> Clusters_Size

# Total Genes by in other genes
ClustersFull %>%
  dplyr::filter( !(V1 %in% All_TFs)) %>%
  dplyr::group_by(V1) %>%
  dplyr::summarise(Genes=n()) %>%
  dplyr::mutate(Class='Other genes') -> ClustersFull_Size

#Clusters_Size <-  rbind(Clusters_Size, ClustersFull_Size)
  
#### tem analysis
subset(Clusters, V1=="Zm00001d044107") # CPP8
subset(Clusters, V1=="Zm00001d010713") # MYB53
subset(Clusters, V1=="Zm00001d032923") # HSF24
subset(Clusters, V1=="Zm00001d035651") # DOF3
subset(Clusters, V1=="Zm00001d006236") # MYB31
subset(Clusters, V1 =="Zm00001d005016") # WRI1

###############################################################



###############################################################
########                      Plots                    ########
###############################################################


## 1
# Plot cluster size
Clusters_Size %>%
  ggplot(aes(x=Genes)) +
  geom_histogram(aes(y=..count..), fill="#E69F00", bins = 30, position="identity", alpha=0.5) +
  #
  #geom_density(alpha=0.6)+
  #geom_vline(aes(xintercept=mean(Genes)), linetype="dashed") +
  labs(x="Genes", y = "Counts") +
  scale_x_continuous(labels = comma, expand = c(0,0)) + 
  theme_pubclean() +
  theme(strip.text.x = element_text(size = 10), 
        axis.text=element_text(size=10), 
        legend.position = 'bottom',
        text = element_text(size=10, family="Times")) -> Plot_1

Plot_1
length(unique(Clusters_Size$V1))
mean(Clusters_Size$Genes)


pdf("Plots/Plot_FigS6b.pdf", width=4, height=2)
print(Plot_1)
dev.off()

####################################################################################

## Combine PWYs and GOs results by Module: Plots 3 & 4
DF_3_4 <- cbind(PWY_Clusters[,c(1,4,7,8)], data.frame('Class'='PWY'))
colnames(DF_3_4)[c(2,3)] <- c("Significant", "Term") 
DF_3_4 <- rbind(DF_3_4, cbind(GOs_ClustersFiltered[,c(7,4,2,9)], data.frame('Class'='GO')))
DF_3_4$Class <- factor(DF_3_4$Class, levels = c("PWY", 'GO'))

DF_3_4 %>% 
  ggplot(aes(y=paste0("Cluster ", Module), 
             x=-log10(FDR), label=Term))+
  geom_point(size=1, alpha=0.5) + 
  geom_text_repel(data= DF_3_4 %>% group_by(Module) %>% top_n(2, -log10(FDR)),
                  size=1,
                  show.legend = F,
                  max.overlaps= Inf,
                  force = 20,
                  min.segment.length = unit(0.2, "lines"),
                  segment.size = 0.2) +
  ylab("") + xlab(bquote(-Log[10] ~ "FDR")) +
  theme_pubclean() +
  theme(strip.text.x = element_text(size = 12), 
        axis.text=element_text(size=12), 
        legend.position = 'none',
        text = element_text(size=12, family="Helvetica")) +
  facet_wrap(~ Class, scales="free_x") -> Plot_PWY_GOs_Clusters


Plot_FigS6abc <- {Plot_ClusterSize + Plot_TFsInCluster + plot_layout(widths = c(1, 0.5))}/Plot_PWY_GOs_Clusters +
  plot_layout(heights = c(0.5, 1.5))

Plot_FigS6abc

tiff("Plots/Plot_ClusterAnnotation.tiff", units="in", width=8, height=10, res=300)
print(Plot_FigS6abc)
dev.off()

## 3
# Plot number of PWYs enriched on clusters

# Plot number of PWYs enriched on clusters
DF_3_4 %>%
  ggplot(aes(y=paste0("Cluster ", Module), 
             x=-log10(FDR), label=gsub('_', ' ', Pathway.name)))+
  geom_point(size=2,
             alpha=0.5) + 
  geom_text_repel(data=subset(PWY_Clusters, n.targ >= 2),
                  size=1,
                  show.legend = F,
                  max.overlaps= Inf,
                  force = 20,
                  min.segment.length = unit(0.2, "lines"),
                  segment.size = 0.2) +
  ylab("") + xlab(bquote(-Log[10] ~ "FDR")) +
  theme_pubclean() +
  theme(strip.text.x = element_text(size = 12), 
        axis.text=element_text(size=12), 
        legend.position = 'none',
        text = element_text(size=12, family="Helvetica")) -> Plot_PWY_Clusters



## 4
# Plot number of PWYs enriched on clusters

####################################################################################
