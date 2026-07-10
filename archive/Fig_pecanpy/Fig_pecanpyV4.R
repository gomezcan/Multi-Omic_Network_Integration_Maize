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

GetGO <- function(degs, mutant){
  
  # Use a list of DEGs and the name of the mutant (string)
  # to identify GOs enriched. Required to have a background predefined
  # Define bacground based on genes in network
  
  background_tem <- background[names(background) %in% unique(ClustersNet$GeneID)]
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
                paste("BP_results/Genes_GOBP_Cluster_", mutant, ".txt", sep = ""),
                sep = '\t', quote = F,
                row.names = F)
  
  #return(list(GOs_DF=Res_DF_BP, Genes=DF_GO_Genes))# return list of GOs-Stats and GeneID-GOs
  return(Res_DF_BP) # Return list of GOs-Stats and GeneID-GOs
}

SuperGO_Modules_Targ <- function(tf){
  ## used TF/Module targets/genes to test GO terms enrichment
  # 1. Select tf/module's Targets
  # 2. Make list file: degs
  # 3. Test enrichment
  
  # Get network by TF
  network <- subset(ClustersNet, M==tf)
  #network <- subset(Res_PecAll, M=="0_4")
  
  network <- subset(network, GeneID %in% Syntenic)
  
  #
  
  Total_targets <- as_tibble(as.data.frame(table(unique(network)$M), stringsAsFactors = F))
  colnames(Total_targets) <- c('Module', 'nModule') 
  
  # Genes input
  degs <- unique(network$GeneID)
  
  
  ## Compare list of predicted targets vs annotated genes in query.vector 
  print(". Pre-GO.")
  out <- GetGO(degs, tf)
  #out <- GetGO(degs, "0_4")
  
  #
  colnames(out)[7] <- "Module"
  
  # save GOs
  out <- left_join(out, Total_targets, by="Module")  
  write.table(out, paste0("BP_results/Module_",tf, "_GOs.txt"), sep = '\t', row.names = F, quote = F)
  
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
#PheGenes <- as_tibble(read.table("Data/Annotations/LinaPheGenes2020.txt", h=T, sep = "\t", quote="", stringsAsFactors = F))

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

# ClusterONE and Pecanpy clusters
# nW10 S5
# nW10, S5, seed all nodes
Clusters <- readLines('DistanceCalculation/MR_Clusters_Dim50_WL80_nW10.s5.nodes_syntenic.csv')

# As list
Clusters <- lapply(Clusters, function (x) strsplit(x, '\t'))
names(Clusters) <- seq(length(Clusters))

# as a DF and remove non-syntenic genes
Clusters <- rbindlist(Clusters, idcol = T)
Clusters  <- subset(Clusters, V1 %in% Syntenic)

## network used on ClusterONE
InputCONE <- fread('DistanceCalculation/InputClusterONE_Dim50_WL80_nW10_syntenic.txt', header = F)

InputCONE_out <- as.data.table(table(InputCONE$V1))
InputCONE_in <- as.data.table(table(InputCONE$V2))
#rm(InputCONE)

hist(InputCONE_out$N, 100)
hist(InputCONE_in$N, 100)


##########################################################
####            Compare Clusters summary              ####
##########################################################

# Label TFs
Clusters[,"isTF"] <-  Clusters$V1 %in% All_TFs

# Genes by cluster
Clusters_Size <- as.data.table(table(Clusters$.id))
Clusters_SizeTF <- as.data.table(table(subset(Clusters, isTF==TRUE)$.id))

#
Clusters_Size <- left_join(Clusters_Size, Clusters_SizeTF, by='V1')
Clusters_Size$N.y[is.na(Clusters_Size$N.y)] <- 0
colnames(Clusters_Size) <- c("Cluster", "Genes", "TFs")
#
plot(Clusters_Size$Genes, Clusters_Size$TFs/Clusters_Size$Genes)

# filter cluster with high content of TFs 
ClustersPass <- Clusters_Size[(Clusters_Size$TFs/Clusters_Size$Genes)<=0.5,]
ClustersPass <- Clusters[Clusters$.id %in% ClustersPass$Cluster, ]

table(PDI %in% ClustersPass$V1)
table(teQTL %in% ClustersPass$V1)
table(CoExp %in% ClustersPass$V1)

# TFs annotation
TFs_inCluster <- unique(ClustersNet[ClustersNet$GeneID %in% All_TFs,])
write.table(TFs_inCluster, "../Fig_MethodsComparison/Network_TFs_Incluster_annotation.txt", 
            row.names = F, quote = F, sep = "\t")



###############################################################

###############################################################
#######              Cluster annotation                 #######
###########          CornCYC Enrichment               #########
###############################################################

Enrichmet_classesV2 <- function(network, targList){
  
  ## Count TF targets in network
  # Count Total
  network <- unique(network[,c("GeneID", "M")])
  #network <- subset(network, Target %in% Syntenic)
  
  # Genes in module
  Total_targets <- as.data.table(table(network$M))
  colnames(Total_targets) <- c('Module', 'GenesModule') 
  Total_targets$Module <- as.character(Total_targets$Module)
  
  # Genes by TF
  totalTarg.byTF <- as.data.frame(mapply(length, targList))
  totalTarg.byTF[,"TF"] <- row.names(totalTarg.byTF)
  colnames(totalTarg.byTF)[1] <- "TF_targets"
  
  
  # list input: network
  network.list <- split(network$GeneID, network$M)
  
  ## Compare list of predicted targets vs annotated genes in query.vector 
  #
  print(". Pre-newGOM .")
  go.obj <- newGOM(network.list, targList, genome.size=length(Syntenic)) # annotated genes in Genome v4
  
  Pval <- getMatrix(go.obj, name="pval")
  Common <- getMatrix(go.obj, name="intersection")
  
  print(". Post-newGOM .")
  
  ### Summary tables
  ## adjust p value
  Pval <- as.data.frame(Pval)
  #Pval[,1:ncol(Pval)] <- apply(Pval[,1:ncol(Pval)], 2, p.adjust)
  
  Pval_table   <- as_tibble(reshape2::melt(as.matrix(Pval)))
  colnames(Pval_table) <- c('Module', 'TF', 'Pval')
  Pval_table$TF <- as.character(Pval_table$TF)
  Pval_table$Module <- as.character(Pval_table$Module)
  
  #
  Common_table <- as_tibble(reshape2::melt(as.matrix(Common)))
  colnames(Common_table) <- c('Module', 'TF', 'n.targ')
  Common_table$TF <- as.character(Common_table$TF)
  Common_table$Module <- as.character(Common_table$Module)
  
  # Add predicted target in class by Module
  Pval_table <- left_join(Pval_table, Common_table , by=c('Module', 'TF'))
  
  # Add total genes in module targets
  Pval_table <- left_join(Pval_table, Total_targets, by="Module")
  
  # add total target by TF
  Pval_table <- left_join(Pval_table, totalTarg.byTF, by="TF")
  
  #write.table(Pval_table, name, sep = "\t", row.names = F, quote = F)
  print(paste("... Done ...", sep = ""))
  return(Pval_table)
}

# set input object
ClustersNet <- ClustersPass[,1:2]
colnames(ClustersNet) <- c("M", "GeneID")
ClustersNet <- ClustersNet[,c("GeneID", "M")]

# PWY enrichment
PWY_Clusters <- Enrichmet_classesV2(ClustersNet, CornCYC.list) 
colnames(PWY_Clusters)[c(2,6)] <- c("PWY",'n.PWY')


# Add PWY name
PWY_Clusters  <- dplyr::left_join(PWY_Clusters, unique(CornCYC[,1:2]), by=c("PWY"="Pathway.id"))

write.table(PWY_Clusters, 
            "PWY_GO_results/NetworkBased_PWY_Clusters_enrichment.txt", 
            sep = '\t', quote = F, row.names = F)

#PWY_Clusters <- fread("PWY_results/PWY_Clusters_enrichment.txt")

# Add FDR by TF
PWY_Clusters %>% 
  #filter(n.targ >= 1) %>% # if not TF's targets in module, do not even try
  group_by(Module) %>%
  mutate("FDR" = p.adjust(Pval, method = 'fdr')) %>%
  filter(Pval <= 0.05) %>%
  dplyr::arrange(Module) -> PWY_Clusters

subset(PWY_Clusters, Pval <=0.01)
View(subset(PWY_Clusters, FDR <=0.1))
length(table(test$Module))
table(test$Module)

table(subset(ClustersFilter, .id=='Cluster.672')$V1 %in% PheGenes$GeneID)

# Examples MYB31
subset(Clusters, V1=="Zm00001d006236") # MYB31
subset(Clusters, V1=="Zm00001d033859") # KN1
subset(Clusters, V1=="Zm00001d005016") # WRI1
subset(Clusters, V1=="Zm00001d051520") # MYB19

subset(Clusters,  .id=="134") # MYB31
subset(Clusters,  .id=="178")  # KN1
subset(Clusters,  .id=="463") # WRI1

table(PDI %in% Clusters$V1)
table(teQTLtfs %in% Clusters$V1)
table(CoExp %in% Clusters$V1)

# TFs not in networks
All_TFsNew <- All_TFs[!(All_TFs %in% c(PDI, teQTLtfs, CoExp))]
All_TFsNew <- All_TFs[!(All_TFs %in% c(PDI, CoExp))]
table(All_TFsNew %in% Clusters$V1)


###############################################################

###############################################################
#######               Cluster annotation               #######
###########            GO enrichment                  #########
###############################################################

# Defined number of clusters
GenesList <- unique(ClustersNet$M)
Lgenes <- length(GenesList)
#Lgenes <- 100

# Size of range to test
w = 5

print(".. Ready to start ..")
for (i in seq(0, Lgenes, w)){
  max=Lgenes
  Start=i+1
  end=i+w
  # 
  if (end<max){
    listtotest <- GenesList[Start:end]
    print(length(listtotest))
    mclapply(listtotest, SuperGO_Modules_Targ, mc.cores=w)
  }
  else{
    r <- paste(Start, max, sep = "-")
    listtotest <- GenesList[Start:max]
    print(length(listtotest))
    w=length(listtotest)
    mclapply(listtotest, SuperGO_Modules_Targ, mc.cores=w)
  }
}

GOs_Clusters <- list.files(path = 'BP_results/', pattern = "^Module_*")

# secong round for cluster not tested
ClusterDone <- gsub("Module_","", GOs_Clusters)
ClusterDone <- gsub("_GOs.txt","",ClusterDone)
GenesList <- GenesList[!(GenesList %in% ClusterDone)]
GenesList <- unique(GenesList)
Lgenes <- length(GenesList)

for (i in GenesList){
  
  GOs_Clusters <- c(GOs_Clusters, SuperGO_Modules_Targ(i))
}

## Remove empty results
GOs_Clusters <- list.files(path = 'BP_results/', pattern = "^Module_*")
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
  #filter(Significant >= 1) %>% # Remove GOs reported without targets
  group_by(Module) %>%
  mutate("FDR" = p.adjust(classic, method = 'fdr')) %>%
  filter(classic <= 0.05) 

subset(GOs_ClustersFiltered, FDR <= 0.01)

GOs_ClustersFiltered
dim(table(GOs_ClustersFiltered$Module))
dim(GOs_ClustersFiltered)

###############################################################

###############################################################
#######              Cluster Description                #######
###########             TFs-based                     #########
###############################################################

colnames(Clusters_Size)[3] <- 'Total.TFs'

# Total new TFs by cluster
TFs_notTested <- All_TFs[!(All_TFs %in% c(PDI, PDIeQTL, CoExp, teQTLtfs))]

# Total new-TFs by cluster
Clusters %>% 
  dplyr::filter(V1 %in%  TFs_notTested) %>%
  dplyr::group_by(.id) %>%
  dplyr::summarise(New.TFs=n()) -> NoTestedTFs_inCLuster

# Total GRN-TFs by cluster
Clusters %>% 
  dplyr::filter(V1 %in%  PDI) %>%
  dplyr::group_by(.id) %>%
  dplyr::summarise(GRN.TFs=n()) -> GRNTFs_inCLuster

# Total eGRN-TFs by cluster
Clusters %>% 
  dplyr::filter(V1 %in%  PDIeQTL) %>%
  dplyr::group_by(.id) %>%
  dplyr::summarise(eGRN.TFs=n()) -> eGRNTFs_inCLuster

# Total CEN-TFs by cluster
Clusters %>% 
  dplyr::filter(V1 %in%  CoExp) %>%
  dplyr::group_by(.id) %>%
  dplyr::summarise(CEN.TFs=n()) -> CENTFs_inCLuster

# Total GAN-TFs by cluster
Clusters %>% 
  dplyr::filter(V1 %in%  teQTLtfs) %>%
  dplyr::group_by(.id) %>%
  dplyr::summarise(GAN.TFs=n()) -> GANTFs_inCLuster


# add info to DF with cluster size
Clusters_Size <- 
  left_join(Clusters_Size, NoTestedTFs_inCLuster, by=c('Cluster'='.id')) %>% 
  left_join(GRNTFs_inCLuster, by=c('Cluster'='.id')) %>% 
  left_join(eGRNTFs_inCLuster, by=c('Cluster'='.id')) %>% 
  left_join(CENTFs_inCLuster, by=c('Cluster'='.id')) %>% 
  left_join(GANTFs_inCLuster, by=c('Cluster'='.id'))

Clusters_Size[is.na(Clusters_Size)] <- 0

Clusters_Size[,'P.Total.TFs'] <- round((Clusters_Size$Total.TFs/Clusters_Size$Genes)*100,2)

#### tem analysis
subset(Clusters, V1=="Zm00001d044107") # CPP8
subset(Clusters, V1=="Zm00001d010713") # MYB53
subset(Clusters, V1=="Zm00001d032923") # HSF24
subset(Clusters, V1=="Zm00001d035651") # DOF3
subset(Clusters, V1=="Zm00001d006236") # MYB31

###############################################################
#######              Cluster annotation                 #######
###########           target Enrichment               #########
###############################################################

# Combined network
Full_net <- subset(Full_net, Target %in% Syntenic)
length(unique(Full_net$Source))
length(unique(Full_net$Target))


Full_net.list <- split(Full_net$Target, Full_net$Source)

Targ_In_Clusters <- Enrichmet_classesV2(ClustersNet, Full_net.list)

write.table(Targ_In_Clusters, "Targ_In_ClustersEnrichment.02012023.txt",
            sep = "\t", row.names = F, quote = F)

Targ_In_Clusters %>% 
  filter(n.targ >= 1) %>% # if not TF's targets in module, do not even try
  group_by(Module) %>%
  mutate("FDR" = p.adjust(Pval, method = 'fdr')) %>%
  filter(FDR <= 0.1) %>% 
  dplyr::arrange(Module) -> Targ_In_Clusters

Targ_In_Clusters

# Count TFs with 
###############################################################

###############################################################
########                      Plots                    ########
###############################################################


## 1
# Plot cluster size
Clusters_Size %>%
  ggplot(aes(x=Genes)) +
  geom_histogram(aes(y=..count..), fill="#E69F00", 
                 bins = 30, position="identity", alpha=0.5)+
  #geom_density(alpha=0.6)+
  #geom_vline(aes(xintercept=mean(Genes)), linetype="dashed") +
  labs(x="Genes in clusters", y = "Number of clusters") +
  scale_x_continuous(labels = comma, expand = c(0,0)) + 
  theme_pubclean() +
  theme(strip.text.x = element_text(size = 10), 
        axis.text=element_text(size=10), 
        legend.position = 'none',
        text = element_text(size=10, family="Times")) -> Plot_1

length(unique(ClusterSizeNewK$NewK))
mean(ClusterSizeNewK$N)
mean(ClusterSizeNewK$TFs_Per)

## 2.1 (v1)
######
# Plot TFs percentage in Cluster
Clusters_Size %>%
  ggplot(aes(x='', y=P.Total.TFs, label=comma(Genes))) +
  geom_boxplot(fill="#9370DB", alpha=0.5, outlier.shape = NA, notch = T) +
  #geom_jitter(aes(size=Genes), alpha=0.5, position = position_jitter(seed = 1)) +
  #geom_text_repel(position = position_jitter(seed = 1), size=3) + 
  labs(y="% TFs in cluster", x = "Clusters") +
  #scale_y_continuous(labels = comma, limits = c(0, 51)) + 
  theme_pubclean() +
  theme(strip.text.x = element_text(size = 10), 
        axis.text=element_text(size=10), 
        legend.position = 'none',
        text = element_text(size=10, family="Times")) -> Plot_2
######

Plot_TFsInCluster

## 2.1 (v2)
# Plot TFs percentage in Cluster
Clusters_Size %>%
  ggplot(aes(x=Genes, y=P.Total.TFs)) +
  geom_point(alpha=0.5, size=1) +
  labs(y="% TFs in cluster", x = "Genes in cluster") +
  scale_x_continuous(expand = c(0,0)) + 
  geom_hline(yintercept=50, linetype="dashed") +
  theme_pubclean() +
  theme(strip.text.x = element_text(size = 10), 
        axis.text=element_text(size=10), 
        legend.position = 'none',
        text = element_text(size=10, family="Times")) -> Plot_2v2

plot_bc <- Plot_1/Plot_2v2
pdf("Plots/plot_bc.pdf", height = 4, width = 4)
print(plot_bc)
dev.off()


Clusters_Size  %>%
  dplyr::filter(P.Total.TFs <= 50) %>%
  dplyr::select(Cluster, Genes, Total.TFs, New.TFs, GRN.TFs, eGRN.TFs, CEN.TFs, GAN.TFs) %>%
  gather(key, value, -Cluster, -Genes) %>% as.data.table -> Clusters_Size_long


ClusterSizeNewK[ClusterSizeNewK$TFs_Per>10, ]

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

###############################################################

#########################################################################
###################         Exploratory analysis      ###################
#########################################################################


vennfunc2 <- function(list_x){
  venn <- Venn(list_x)
  data <- process_data(venn)
  colorGroups <- c(CEN = 'yellow1', Cluster= "orange2", GRN='steelblue1', GAN='darkorchid1')
  colfunc <- colorRampPalette(colorGroups)
  col <- colfunc(4)
  
  colorGroups <- c(CEN="gray98", Cluster= "gray99", GRN="gray99", GAN="gray98")
  colfunc2 <- colorRampPalette(colorGroups)
  col2 <- colfunc2(15)
  
  ggplot() +
    geom_sf(aes(fill=name), data = venn_region(data), show.legend = F) +
    geom_sf(aes(color=name), size = 1.5, data = venn_setedge(data), show.legend = F) +
    #
    geom_sf_text(aes(label = name), size=6, data = venn_setlabel(data)) +
    geom_sf_text(aes(label= scales::comma(count, accuracy = 1)), size=5, data = venn_region(data)) +
    #
    scale_fill_manual(values = col2) + # 
    scale_color_manual(values = alpha(col, .5), ) +
    #
    theme_void() +
    theme(plot.margin = unit(c(0.5, 1, 1.5, 0.1), "cm")) #+
  #xlim(-150,1000)
}

SummaryTF_cluster <- function(tf){
  # Genes in cluster/Module
  c <- subset(Res_PecAll, GeneID=="Zm00001d006236")$M
  targ <- subset(Res_PecAll, M==c)[,c("GeneID", "M")]
  return(targ)
}


# tem with MYB31
tem <- SummaryTF_cluster("Zm00001d006236")

MYB31_PecanAll <- list(Cluster=tem$GeneID, 
                       CEN=subset(CoExp, Source=="Zm00001d006236")$Target,
                       GRN=subset(PDI, Source=="Zm00001d006236")$Target,
                       GAN=subset(teQTL, Source=="Zm00001d006236")$Target)

vennfunc2(MYB31_PecanAll)

subset(NE_All_CornC_Pecap, Module=="0_4") # MYB31 present in module



Module_MYB31_GOs_Pep <- GetGO(tem$GeneID, "Zm00001d006236", "Pecanpy_All")


NewDAP <- fread("../Fig_PDI/NewDAPseq_All_Peals_02.05.2020.txt", h=T)
NewDAP[,"TF"] <- NewDAP$TFsample <- gsub("newDAP.", "", NewDAP$TFsample)
NewDAP[,"TF"] <- sapply(strsplit(NewDAP$TF, split='_', fixed=TRUE), `[`, 1) # add methods label

NewDAP <- subset(NewDAP, TF %in% Top45[Top45$V1 %in% tem$GeneID, ]$V1)

NewDAP <- subset(NewDAP, OCR==1)
NewDAP <- subset(NewDAP, Z > -0.5)
NewDAP <- subset(NewDAP, abs(Dis) <= 3)
NewDAP <- as_tibble(unique(NewDAP[,c(3,9)]))

ReplaceName(unique(NewDAP$TF))

DOF23_PecanAll <- list(Cluster=tem$GeneID, 
                       CEN=subset(CoExp, Source=="Zm00001d026096")$Target,
                       GRN=subset(NewDAP, TF=="Zm00001d026096")$Target,
                       GAN=subset(teQTL, Source=="Zm00001d026096")$Target)

MYB19_PecanAll <- list(Cluster=tem$GeneID, 
                       CEN=subset(CoExp, Source=="Zm00001d051520")$Target,
                       GRN=subset(NewDAP, TF=="Zm00001d051520")$Target,
                       GAN=subset(teQTL, Source=="Zm00001d051520")$Target)

lapply(DOF23_PecanAll, length)
lapply(MYB19_PecanAll, length)

vennfunc2 <- function(list_x){
  venn <- Venn(list_x)
  data <- process_data(venn)
  colorGroups <- c(CEN = 'yellow1', Cluster= "orange2", GRN='steelblue1', GAN='darkorchid1')
  colfunc <- colorRampPalette(colorGroups)
  col <- colfunc(4)
  
  colorGroups <- c(CEN="gray98", Cluster= "gray99", GRN="gray99", GAN="gray98")
  colfunc2 <- colorRampPalette(colorGroups)
  col2 <- colfunc2(15)
  
  ggplot() +
    geom_sf(aes(fill=name), data = venn_region(data), show.legend = F) +
    geom_sf(aes(color=name), size = 1.5, data = venn_setedge(data), show.legend = F) +
    #
    geom_sf_text(aes(label = name), size=6, data = venn_setlabel(data)) +
    geom_sf_text(aes(label= scales::comma(count, accuracy = 1)), size=5, data = venn_region(data)) +
    #
    scale_fill_manual(values = col2) + # 
    scale_color_manual(values = alpha(col, .5), ) +
    #
    theme_void() +
    theme(plot.margin = unit(c(1, 1, 0.5, 0.1), "cm")) +
    xlim(-0.1,1)
}

ggarrange((vennfunc2(DOF23_PecanAll) + labs(title="DOF23")),
          (vennfunc2(MYB19_PecanAll) + labs(title="MYB19")))

library(ComplexHeatmap)

Enrichmet_List <- function(targList){
  
  ## Count TF targets in network
  # Count Total
  
  # Genes by TF
  totalTarg.byTF <- as.data.frame(mapply(length, targList))
  totalTarg.byTF[,"TF"] <- row.names(totalTarg.byTF)
  colnames(totalTarg.byTF)[1] <- "TF_targets"
  totalTarg.byTF <- as_tibble(totalTarg.byTF)
  
  # list input: network
  network.list <- split(network$GeneID, network$M)
  
  ## Compare list of predicted targets vs annotated genes in query.vector 
  #
  print(". Pre-newGOM .")
  targList = MYB19_PecanAll
  targList= DOF23_PecanAll
  go.obj <- newGOM(targList, targList, genome.size=length(Syntenic)) # annotated genes in Genome v4
  
  MA <- getMatrix(go.obj, name="pval")
  Common <- getMatrix(go.obj, name="intersection")
  
  #Pval <- -log10(Pval)
  diag(Pval) <- 1
  Pval <- (Pval <=0.05)*1
  Heatmap(Pval, 
          cluster_rows = F,
          cluster_columns = F 
  )
  
  
  
  
  print(". Post-newGOM .")
  
  ### Summary tables
  ## adjust p value
  Pval <- as.data.frame(Pval)
  #Pval[,1:ncol(Pval)] <- apply(Pval[,1:ncol(Pval)], 2, p.adjust)
  
  
  Pval_table   <- as_tibble(reshape2::melt(as.matrix(Pval)))
  colnames(Pval_table) <- c('Module', 'TF', 'Pval')
  Pval_table[,1:2] <- apply(Pval_table[,1:2], 2, as.character)
  
  #print(Pval_table)
  #
  
  Common_table <- as_tibble(reshape2::melt(as.matrix(Common)))
  colnames(Common_table) <- c('Module', 'TF', 'n.targ')
  Common_table[,1:2] <- apply(Common_table[,1:2], 2, as.character)
  #print(Common_table)
  
  # Add predicted target in class by Module
  
  Pval_table <- left_join(Pval_table, Common_table , by=c('Module', 'TF'))
  
  # Add total genes in module targets
  Pval_table <- left_join(Pval_table, Total_targets, by="Module")
  
  # add total target by TF
  Pval_table <- left_join(Pval_table, totalTarg.byTF, by="TF")
  
  # Select significant TFs 
  Pval_table <- subset(Pval_table, Pval <= 0.05)
  #Pval_table <- tibble(TF=c("test", "a", "b"), TF2="test2")
  
  #write.table(Pval_table, name, sep = "\t", row.names = F, quote = F)
  print(paste("... Done ...", sep = ""))
  return(Pval_table)
}


#########################################################################
#####               GRN and & CEN only networks                     #####
#########################################################################

NetworkEnrichment_GRN_CEN <- Enrichmet_classes(Res_PecGRN_CEN, nrow(Res_PecAll)) 
NetworkEnrichment_tsne_GRN_CEN <- Enrichmet_classes(Res_tsneGRN_CEN, nrow(Res_tsneGRN_CEN)) 

NetworkEnrichment_GRN_CEN  <- left_join(NetworkEnrichment_GRN_CEN, unique(CornCYC[,1:2]), by=c("PWY"="Pathway.id"))
NetworkEnrichment_tsne_GRN_CEN <- left_join(NetworkEnrichment_tsne_GRN_CEN, unique(CornCYC[,1:2]), by=c("PWY"="Pathway.id"))
NetworkEnrichment_GRN_CEN <- left_join(NetworkEnrichment_GRN_CEN, CornCYCred_size, by="PWY")
NetworkEnrichment_tsne_GRN_CEN <- left_join(NetworkEnrichment_tsne_GRN_CEN, CornCYCred_size, by="PWY")

NE_PDI_PecaP_GRNCEN   <- Enrichmet_classesV2(Res_PecGRN_CEN, nrow(Res_PecGRN_CEN), PDI.list) 
NE_CoExp_PecaP_GRNCEN <- Enrichmet_classesV2(Res_PecGRN_CEN, nrow(Res_PecGRN_CEN), CoExp.list)
NE_teQTL_PecaP_GRNCEN <- Enrichmet_classesV2(Res_PecGRN_CEN, nrow(Res_PecGRN_CEN), teQTL.list)  

NE_PDI_tsne_GRNCEN   <- Enrichmet_classesV2(Res_tsneGRN_CEN, nrow(Res_tsneGRN_CEN), PDI.list) 
NE_CoExp_tsne_GRNCEN <- Enrichmet_classesV2(Res_tsneGRN_CEN, nrow(Res_tsneGRN_CEN), CoExp.list)
NE_teQTL_tsne_GRNCEN <- Enrichmet_classesV2(Res_tsneGRN_CEN, nrow(Res_tsneGRN_CEN), teQTL.list)

#########################################################################









