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
  #
  # mutant: description of a category of source of the data
  # degs: set of genes to test enrichment
  
  GeneList <- factor(as.integer(background_IDs %in% degs))
  names(GeneList) <- background_IDs
  GOdata_BP <- new("topGOdata", ontology = "BP", allGenes = GeneList, annot = annFUN.gene2GO, gene2GO = background)
  
  #GOdata_MF <- new("topGOdata", ontology = "MF", allGenes = GeneList, annot = annFUN.gene2GO, gene2GO = background)
  #GOdata_CC <- new("topGOdata", ontology = "CC", allGenes = GeneList, annot = annFUN.gene2GO, gene2GO = background)
  
  #### Define test ####
  test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test", nodeSize = 10)
  
  ### test enrichment 
  results_BP <- getSigGroups(GOdata_BP, test.stat)
  
  ### save pdf Graph
  #namepdf=paste("GOs_Plots/GO.BP_",mutant, "", sep = "")
  #printGraph(GOdata_BP, results_BP, firstSigNodes=20,  fn.prefix = namepdf, useInfo = "def", pdfSW = TRUE) #
  
  ######## Get Significant GOs ########  
  Res_DF_BP <- as_tibble(GenTable(GOdata_BP, classic = results_BP, topNodes = 500, orderBy='Fis')) # save as dataframe
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
  
  return(list(GOs_DF=Res_DF_BP, Genes=DF_GO_Genes))# return list of GOs-Stats and GeneID-GOs
}

SuperGO_Moludes_Targ <- function(tf){
  ## used TF/Module targets/genes to test GO terms enrichment
  # 1. Select tf/module's Targets
  # 2. Make list file: degs
  # 3. Test enrichment
  
  # Get network by TF
  network <- subset(Res_PecAll, M==tf)
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
  colnames(out$GOs_DF)[7] <- "Module"
  colnames(out$Genes)[5] <- "Module"
  
  # save GOs
  out$GOs_DF <- subset(out$GOs_DF, classic <= 0.05)
  out$GOs_DF <- left_join(out$GOs_DF, Total_targets, by="Module")  
  
  
  write.table(out$GOs_DF, paste0("BP_results/Module_",tf, "_GOs.txt"), sep = '\t', row.names = F, quote = F)
  
  # save gene annotation GOs  
  write.table(out$Genes, paste0("BP_results/Module_",tf, "_Annotation.txt"), sep = '\t', row.names = F, quote = F)
  
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
#TF_CoR <- as_tibble(read.table("Data/Annotations/TF_CoR_Mazie.txt", h=T, stringsAsFactors = F))

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
CornCYCred  <- subset(CornCYC, !(Pathway.id %in% row.names(subset(CornCYC_size, Freq == 1))))
CornCYC.list <- split(CornCYCred$GeneID, CornCYCred$Pathway.id)

CornCYCred_size <- as_tibble(as.data.frame(table(CornCYCred$Pathway.id), stringsAsFactors = F))
colnames(CornCYCred_size) <- c("PWY", "nPWY")

##################################################
######        Setup pecanpy entry          #######
##################################################

# PDI
PDI <- as_tibble(unique(fread("../Fig_PDI/Only_PDI_NetworkFinal.10_14_2022.txt")[,c(2,3)]))
colnames(PDI)[1] <- "Source"

PDIeQTL <- as_tibble(unique(fread("../Fig_PDI/CisE_PDI_NetworkFinal.10_14_2022.txt")[,c(2,3)]))
colnames(PDIeQTL)[1] <- "Source"

# CoExp
CoExp <- as_tibble(unique(fread("../Fig_Coexpression/CoExp_NetworkFinal.10_11_2021.txt")))
colnames(CoExp)[2] <- "Source"
CoExp <- unique(CoExp[,2:3])

# teQTL
teQTL <- as_tibble(unique(fread("../Fig_transeQTL/teQTL_NetworkFinal.10_11_2021.txt")))
colnames(teQTL)[1] <- "Source"

# Pecanpy and k-means results
ClustersPecan <- fread('Two_Kmean_to_Pecanpy_uniqFullNets_weighted.txt')

##################################################
##########        Pencanty          ##############
##################################################
#
teQTLtf <- subset(teQTL, Source %in% unique(c(TF_CoR$GeneID, PDI$Source, CoExp$Source))) 
All_TFs <- unique(c(PDI$Source, CoExp$Source, TF_CoR$GeneID))
#

peGAN <- as.data.frame(fread("pecanpy_GAN.txt" , skip=1, header=F))
peGRN <- as.data.frame(fread("pecanpy_GRN.txt" , skip=1, header=F))
peCEN <- as.data.frame(fread("pecanpy_CEN.txt" , skip=1, header=F))

row.names(peGAN) <- peGAN$V1
row.names(peGRN) <- peGRN$V1
row.names(peCEN) <- peCEN$V1

# peGAN <- peGAN[,-c(1)]
# peGRN <- peGRN[,-c(1)]
# peCEN <- peCEN[,-c(1)]

peGAN[1:5,1:5]
peGRN[1:5,1:5]
peCEN[1:5,1:5]


##################################################
######     dimension reduction: PCA       ########
##################################################

pcaGAN <- pca(peGAN[,-c(1)]) # removeVar = 0.01
pcaGRN <- pca(peGRN[,-c(1)])
pcaCEN <- pca(peCEN[,-c(1)])

screeplot(pcaGAN, axisLabSize = 18, titleLabSize = 22, components = 1:20)
screeplot(pcaGRN, axisLabSize = 18, titleLabSize = 22, components = 1:20)
screeplot(pcaCEN, axisLabSize = 18, titleLabSize = 22, components = 1:20)

##################################################

##################################################
######      Combine pecanpy       ########
##################################################

length(row.names(peGAN))
length(row.names(peGRN))
length(row.names(peCEN))

# scale by gene
peGAN[,-c(1)] <- apply(peGAN[,-c(1)], 2, scale)
peGRN[,-c(1)] <- apply(peGRN[,-c(1)], 2, scale)
peCEN[,-c(1)] <- apply(peCEN[,-c(1)], 2, scale)

colnames(peGAN)[-c(1)] <- paste0("GAN_", colnames(peGAN)[-c(1)])
colnames(peGRN)[-c(1)] <- paste0("GRN_", colnames(peGRN)[-c(1)])
colnames(peCEN)[-c(1)] <- paste0("CEN_", colnames(peCEN)[-c(1)])

peAll <- as_tibble(data.frame(gid=Reduce(intersect, list(row.names(peGAN), row.names(peGRN), row.names(peCEN)))))
peGRN_CEN <- as_tibble(data.frame(gid=Reduce(intersect, list(row.names(peGAN), row.names(peGRN), row.names(peCEN)))))

peAll <- left_join(peAll, peGRN, by=c("gid"="V1"))
peAll <- left_join(peAll, peCEN, by=c("gid"="V1"))
peAll <- left_join(peAll, peGAN, by=c("gid"="V1"))

peGRN_CEN <- left_join(peAll, peGRN, by=c("gid"="V1"))
peGRN_CEN <- left_join(peAll, peCEN, by=c("gid"="V1"))

write.table(peAll, "All_pecanpy_layer.txt", row.names = F, quote = F, sep = "\t")
write.table(peGRN_CEN, "GRN_CEN_pecanpy_layer.txt", row.names = F, quote = F, sep = "\t")

#dim(peAll)
#heatmap.2(as.matrix(peGRN_CEN[1:100,-c(1)]),trace="none")


#################################################
###  Read Modules from "OptimumK_pepanpy.ipynb"
#################################################

Res_PecAll <- as_tibble(read.table("Cluster_Round2_all_Pecapy.txt", h=T))[,2:4]
Res_PecGRN_CEN <- as_tibble(read.table("Cluster_Round2_GRN_CEN_Pecapy.txt", h=T))[,2:4]
colnames(Res_PecAll)[3] <- "M"
colnames(Res_PecGRN_CEN)[3] <- "M"


####################################################
########     Test with known regulators     ########
########     CornCYC Enrichment test        ########
####################################################


Enrichmet_classesV2 <- function(network, targList){
  
  ## Count TF targets in network
  # Count Total
  network <- unique(network[,c("GeneID", "M")])
  #network <- subset(network, Target %in% Syntenic)
  
  # Genes in module
  Total_targets <- as_tibble(as.data.frame(table(network$M)))
  colnames(Total_targets) <- c('Module', 'GenesModule') 
  Total_targets$Module <- as.character(Total_targets$Module)
  
  
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

#####################################################
##### Metabolic genes from Corn Cyc.            
#### testing enrichment of modules in CornCyc PWYs 
#####################################################

# Pecanpy
NE_All_CornC_Pecap <- Enrichmet_classesV2(Res_PecAll, CornCYC.list) 
colnames(NE_All_CornC_Pecap)[2] <- "PWY"


# add PWY name
NE_All_CornC_Pecap  <- left_join(NE_All_CornC_Pecap, unique(CornCYC[,1:2]), by=c("PWY"="Pathway.id"))
write.table(NE_All_CornC_Pecap, "ModulesPecapy.CornCYC.04_19_2021.txt", sep = "\t", row.names = F, quote = F)



#### testing enrichment of modules in Targets genes
PDI.list <- split(PDI$Target, PDI$Source)
CoExp.list <- split(CoExp$Target, CoExp$Source)
teQTL.list <- split(teQTL$Target, teQTL$Source)

# combined network
Full_net <- rbind(PDI, CoExp, teQTL)
Full_net <- unique(subset(Full_net, Source %in% All_TFs))
Full_net.list <- split(Full_net$Target, Full_net$Source)


# Second test with TF's targets
NE_PDI_PecaP   <- Enrichmet_classesV2(Res_PecAll, PDI.list) 
NE_CoExp_PecaP <- Enrichmet_classesV2(Res_PecAll, CoExp.list)
NE_teQTL_PecaP <- Enrichmet_classesV2(Res_PecAll, teQTL.list)
#
NE_FullNet_PecaP <- Enrichmet_classesV2(Res_PecAll, Full_net.list)
#

write.table(NE_PDI_PecaP, "ModulesPecapy.PDI.04_19_2021.txt", sep = "\t", row.names = F, quote = F)
write.table(NE_CoExp_PecaP, "ModulesPecapy.CoExp.04_19_2021.txt", sep = "\t", row.names = F, quote = F)
write.table(NE_teQTL_PecaP, "ModulesPecapy.teQTL.04_19_2021.txt", sep = "\t", row.names = F, quote = F)

dim(subset(as.data.frame(table(NE_PDI_PecaP$TF)), Freq>0))
dim(subset(as.data.frame(table(NE_CoExp_PecaP$TF)), Freq>0))
dim(subset(as.data.frame(table(NE_teQTL_PecaP$TF)), Freq>0))




#####################################################
#### Modules GO enrichment
#####################################################


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

# GOs term annotations
background <- readMappings("../Fig_GOs/synteny.ID_TopGO_V4_GRAMER.txt")
background_IDs <- as.character(unique(names(background)))


Lgenes <- length(unique(Res_PecAll$M))
GenesList <- unique(Res_PecAll$M)

w=40 # Size of range to test
print(".. Ready to start ..")

for (i in seq(0, Lgenes, w)){
  max=Lgenes
  Start=i+1
  end=i+w
  if (end<max){
    listtotest <- GenesList[Start:end]
    print(length(listtotest))
    mclapply(listtotest, SuperGO_Moludes_Targ, mc.cores=w)
    
  }
  else{
    r <- paste(Start, max, sep = "-")
    listtotest <- GenesList[Start:max]
    print(length(listtotest))
    mclapply(listtotest, SuperGO_Moludes_Targ, mc.cores=w)
  }
}

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

Get_WC(subset(Module_MYB31_GOs_Pep$BP, Annotated>5 & classic <= 0.05)$Term)

table(tem$GeneID %in% TF_CoR$GeneID)
table(tem$GeneID %in% Y1H$TF.v4)
table(tem$GeneID %in% Top45$V1)

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









