library(factoextra)
library(circlize)
library(ComplexHeatmap)
library(hrbrthemes)
library(scales)
library(tidyverse)
library(data.table)
library(ggVennDiagram)
library(GeneOverlap)
library(topGO)
library(purrr)
library(gplots)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(viridis)
library(patchwork)
library(reshape2)
library(rrvgo)
library(org.Zmays.eg.db)
library(UpSetR)
library(fgsea)
library(GOSemSim)
library(ggpmisc)
library(ggpointdensity)
library(edgeR)

#packageVersion('edgeR')

##################################################
##########          Functions        #############
##################################################

#
ReplaceName <- function(ids){
  # for (i in 1:nrow(Top45)){
  #   ids <- gsub(Top45$V1[i], Top45$V2[i], ids)
  # }
  
  for (i in 1:nrow(TFdic)){
    ids <- gsub(TFdic$V2[i], TFdic$V1[i], ids)
  }
  return(ids)
}

ReplaceGO <- function(ids){
  # for (i in 1:nrow(Top45)){
  #   ids <- gsub(Top45$V1[i], Top45$V2[i], ids)
  # }
  
  for (i in 1:nrow(GOdic)){
    ids <- gsub(GOdic$parent[i], GOdic$parentTerm[i], ids)
  }
  return(ids)
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

MapGOparent <- function(GO_vector){
  # library(rrvgo)
  
  if (length(GO_vector) >1) {
    # Semantic similarity
    simMatrix <- calculateSimMatrix(GO_vector,  
                                    orgdb=org.Zmays.eg.db,  ont="BP", 
                                    semdata=Zm.GOSemSim.BP,
                                    method="Wang")
    # Reduce term
    reducedTerms <- reduceSimMatrix(simMatrix, 
                                    keytype="GENENAME",
                                    threshold=0.7, 
                                    orgdb=org.Zmays.eg.db)
    
    # treemapPlot(reducedTerms)
    #
    return(list(DF=reducedTerms, simM=simMatrix))
  }
  
}

library(dendextend)
##################################################
##########        Annotations       ##############
##################################################
## Syntenic genes 
Syntenic <- as_tibble(read.table("../Data/Annotations/Zm.v4.synteny.genes.txt", h=T, stringsAsFactors = F, sep = '\t'))
Syntenic <- subset(Syntenic, gene_synteny=='syntenic')$gene_id

#### top 45 ###
#Top45 <- as_tibble(read.table("Data/Annotations/Top45.txt", h=F, stringsAsFactors = F))

# TF names
TFdic <- as_tibble(read.table("Data/Annotations/TF_Id_Name.txt", h=F, stringsAsFactors = F))

# TF names
GOdic <- unique(fread("../Fig_MethodsComparison/ReduceGOterms_All_methods.txt",  header =T)[,1:2])

# Phenolic related genes
#PheGenes <- as_tibble(read.table("Data/Annotations/LinaPheGenes2020.txt", h=T, sep = "\t", quote="", stringsAsFactors = F))

# TF and CoReg
TF_CoR <- as_tibble(read.table("../Data/Annotations/TF_CoR_Mazie.txt", h=T, stringsAsFactors = F))
TF_CoR$Family <- gsub('Others', 'Other', TF_CoR$Family)

# All TFs
All_TFs <- fread("../Fig_pecanpy/All_TFs.txt", header = F)$V1

## Y1H network
Y1H <- as_tibble(read.table("Data/Annotations/Y1H.data.txt", h=T, stringsAsFactors = F))[,4:6]
#ReplaceName(Y1H$TF.v4)

# GOs term annotations
background <- readMappings("Data/Annotations/synteny.ID_TopGO_V4_GRAMER.txt")
background_IDs <- as.character(unique(names(background)))

# GO list 
GO.list <- rbindlist(lapply(background, as_tibble), idcol = T)

# Map GOs to parents 
GO.list.parent <- MapGOparent(unique(GO.list$value))
# scatterPlot(GO.list.parent$simM, GO.list.parent$DF)

# Discard GO not annotated as BP 
mask  <- unique(GO.list$value)[(unique(GO.list$value) %in% c(GO.list.parent$DF$go, GO.list.parent$DF$parent))]
length(mask)

# Add parents to gene table
GO.list <- left_join(GO.list, GO.list.parent$DF[,c("go","parent", "size")], by=c("value"="go"))

# Full network before MR_mi calculation
FullNet <- fread("../Fig_CommonTarg/Full_Final_network.11022022.txt")

TF_out_Fullnet <- as.data.table(table(FullNet$Source))


# TF target net from Network-based
MRMI_full <- fread('../Fig_pecanpy/DistanceCalculation/InputClusterONE_Dim50_WL80_nW10_0.005_syntenic.txt')
MRMI_full  <- subset(MRMI_full, V2 %in% Syntenic)
MRMI_full  <- subset(MRMI_full, V1 %in% Syntenic)


# Subset of cluster to TFs data
MRMI_full  <- unique(subset(MRMI_full, V1 %in% All_TFs)[,1:2])
colnames(MRMI_full) <- c("TF", "Target")
TF_out_MRMI <- as.data.table(table(MRMI_full$TF))

MRMI_TFnet  <- fread('../Fig_pecanpy/DistanceCalculation/InputClusterONE_Dim50_WL80_nW10_0.005_syntenic.txt')
MRMI_TFnet  <- subset(MRMI_TFnet, V2 %in% Syntenic)
MRMI_TFnet  <- subset(MRMI_TFnet, V1 %in% Syntenic)
MRMI_TFnet  <- subset(MRMI_TFnet, V2 %in% All_TFs)
MRMI_TFnet  <- unique(subset(MRMI_TFnet, V1 %in% All_TFs)[,1:2])
colnames(MRMI_TFnet) <- c("TF1", "TF2")

# Files from GSS calculation
# ## 1. Calculate semantic similarity background
suppressMessages(library(GOSemSim))

# 
# # Pre-calculate semantic similarity
Zm.GOSemSim.BP <- GOSemSim::godata(org.Zmays.eg.db, keytype="GENENAME", ont = 'BP')

################################################


################################################
##          Remove GOs To generals            ##
################################################

colnames(GO.list) <- c('GeneID', "GO.ID", "pGO.ID", 'pGO.IDsize')

# Add GO size
GO.list <- left_join(GO.list, 
                     as.data.table(table(GO.list$GO.ID)),
                     by=c('GO.ID'="V1"))

colnames(GO.list)[5] <- 'GO.IDsize'

# keep GOs with BP parent
GO.listV2 <- GO.list[!is.na(GO.list$pGO.ID),] 

# GO list background with parent GO mapped
Genes_GOv1 <- unique(GO.list$GeneID)
Genes_GOv2 <- unique(GO.listV2$GeneID)

length(Genes_GOv1)
length(Genes_GOv2)

# Size: genes in parent term
GO.listV2 <- subset(GO.listV2, size <= 1500)
Genes_in_V2 <- unique(GO.listV2$.id)

saveRDS(GO.listV2, 'MaizeSyntenicGenes_GOparent.rds')
#treemapPlot(GO.list.parent)

################################################


################################################
##   Description of GOs from network: Filters ##
################################################

# Used reduce GO terms data
GO_Network <- fread("../Fig_MethodsComparison/Total_NetworkBased_predictions.txt") %>%
  filter(Annotation == "GO") %>%
  dplyr::select(TF, Ann.ID) %>% unique()

length(unique(GO_Network$Ann.ID))

## Filter 1: Keep GO in BP categories
# Number of genes by GO terms
GO_Network <- left_join(GO_Network, unique(GO.listV2[,2:5]), by=c('Ann.ID'='GO.ID'))

# Remove GO terms without BP parent 
GO_Network <- GO_Network[!is.na(GO_Network$pGO.ID),]
length(unique(GO_Network$Ann.ID))
length(unique(GO_Network$pGO.ID))

## Filter 2: Discard GO to general
## 
# mask1: TFs with GOs without size threshold
mask1 <- subset(GO_Network, pGO.IDsize <= 500 | GO.IDsize <= 500)
mask1 <- mask1$TF

# mask2: TFs with GOs out of size threshold
mask2 <- subset(GO_Network, pGO.IDsize > 500 & GO.IDsize > 500)
mask2 <- mask2$TF

# Discard only TFGOs indexes present in both mask1 and mask2. So, discard
# GO if there  any lower hierarchy GO (defined by genes in GO). 
mask3 <- mask2[mask2 %in% mask1]
length(mask3)

# final TF-GO netowrk: defined as GO in treshold and TF in Mask
mask <- subset(GO_Network, pGO.IDsize > 500 & GO.IDsize > 500 & TF %in% mask3)
mask <- paste0(mask$TF, mask$Ann.ID, mask$pGO.ID)

# GO_Network v2
GO_Networkv2 <- !(paste0(GO_Network$TF, GO_Network$Ann.ID, GO_Network$pGO.ID) %in% mask)
GO_Networkv2 <- GO_Network[GO_Networkv2,]

length(unique(GO_Networkv2$TF))
GO_Networkv2 <- GO_Networkv2[GO_Networkv2$GO.IDsize <= 500,]
GO_Networkv2 <- GO_Networkv2[GO_Networkv2$GO.IDsize <= 500,]
length(unique(test$TF))

#length(unique(GO_Networkv2$TF))
## Filter 3: Map GO to specific to their corresponding parents
# Count total TF-GO associations and its association 
GO_Networkv2[,"Index"] <- seq(1, nrow(GO_Networkv2))

# with GO size
mask1 <- subset(GO_Networkv2, GO.IDsize <= 50)
mask2 <- subset(GO_Networkv2, GO.IDsize > 50)

mask1$Ann.ID <- mask1$pGO.ID
mask1$GO.IDsize <- mask1$pGO.IDsize

GO_Networkv3 <- rbind(mask1, mask2)

length(unique(GO_Networkv2$Ann.ID))
length(unique(GO_Networkv3$Ann.ID))

length(unique(GO_Networkv2$TF))
length(unique(GO_Networkv3$TF))


################################################

################################################
##  Description of GOs from network: Degree   ##
################################################

# Count Tf per GO term
InDegree <- as.data.table(table(GO_Networkv3$Ann.ID)) %>%
  dplyr::rename('Ann.ID'='V1')

InDegree  <- left_join(InDegree, unique(GO_Networkv3[,2:5]), by="Ann.ID")

# Count GOs per TF
OutDegree <- as.data.table(table(GO_Networkv3$TF)) %>%
  dplyr::rename('TF'='V1')

OutDegree <- left_join(OutDegree, TF_out_MRMI, by=c('TF'='V1'))
colnames(OutDegree) <- c("TF", 'GOs', 'Targ.MR')

OutDegree <- left_join(OutDegree, TF_out_Fullnet, by=c('TF'='V1'))
colnames(OutDegree) <- c("TF", 'GOs', 'Targ.MR', 'Targ.Full')

# Degree Freq
InDegree_Table  <- as.data.table(table(InDegree$N)) %>%
  dplyr::rename('Indegree'='V1') %>%
  dplyr::mutate(Indegree= as.numeric(Indegree))

#
OutDegree_Table <- as.data.table(table(OutDegree$GOs)) %>%
  dplyr::rename('Outdegree'='V1') %>%
  dplyr::mutate(Outdegree= as.numeric(Outdegree))


########
# Define terms Freq by degree
########

## map to parent and get Dis matrix
Red_objv2 <- MapGOparent(unique(GO_Networkv2$Ann.ID))
Red_objv3 <- MapGOparent(unique(GO_Networkv3$Ann.ID))

# Red_obj_for top GOs more highly connected
TopGOs_RedObj <-  subset(Red_objv2$DF, go %in% subset(InDegree, N>=10)$Ann.ID)
TailGOs_RedObj <- subset(Red_objv2$DF, go %in% subset(InDegree, N < 10)$Ann.ID)

GOs_RedObj1 <- subset(Red_objv2$DF, go %in% subset(InDegree, N==1)$Ann.ID)
GOs_RedObj2 <- subset(Red_objv2$DF, go %in% subset(InDegree, N==2)$Ann.ID)

 

InDegree[,'DegreeClass'] <- NA

InDegree$DegreeClass[InDegree$N == 1] <- 1
InDegree$DegreeClass[InDegree$N == 2] <- 2
InDegree$DegreeClass[InDegree$N == 3] <- 3
InDegree$DegreeClass[InDegree$N == 4] <- 4
InDegree$DegreeClass[InDegree$N == 5] <- 5
InDegree$DegreeClass[InDegree$N >= 6 &  InDegree$N < 11 ] <- 10
InDegree$DegreeClass[InDegree$N >= 11 &  InDegree$N <= 100] <- 100
InDegree$DegreeClass[InDegree$N > 100] <- 101

wordcloudDF <-  function(reducedTerms, onlyParents = TRUE) {
  if (!all(sapply(c("wordcloud", "tm", "slam"), requireNamespace, 
                  quietly = TRUE))) {
    stop("Package wordcloud and/or its dependencies (tm, slam) not available. ", 
         "Consider installing it before using this function.", 
         call. = FALSE)
  }
  if (onlyParents) {
    x <- tm::Corpus(tm::VectorSource(reducedTerms$term[reducedTerms$parent == 
                                                         reducedTerms$go]))
  }
  else {
    x <- tm::Corpus(tm::VectorSource(reducedTerms$term))
  }
  tdm <- tm::TermDocumentMatrix(x, control = list(removePunctuation = TRUE, 
                                                  stopwords = TRUE))
  m <- as.matrix(tdm)
  v <- sort(rowSums(m), decreasing = TRUE)
  d <- data.frame(word = names(v), freq = v)
  
  out <- data.table(Word = d$word, Freq = d$freq)
  return(out)
}

GOterm_DB <- tibble(matrix(nrow = 0, ncol = 3))
colnames(GOterm_DB) <- c("Word", "Freq","Degree")

for (i in sort(as.numeric(unique(InDegree$DegreeClass)))){
  # DF with GO reduced
  tem_RedObj <-  subset(Red_objv2$DF, go %in% subset(InDegree, DegreeClass==i)$Ann.ID)
  
  # Counts words
  tem_wcObj <- wordcloudDF(tem_RedObj, onlyParents=T)
  #
  tem_wcObj[,'Degree'] <- i
  GOterm_DB <- rbind(GOterm_DB, tem_wcObj)
  #print(i)
}


GOterm_DB %>% 
  group_by(Word) %>%
  mutate(Total=sum(Freq)) %>%
  mutate(FreqN=Freq/Total) -> GOterm_DB
  
GOterm_DB$Degree <- factor(GOterm_DB$Degree,
                           as.character(sort(as.numeric(unique(InDegree$DegreeClass)))))


M_GOterm_DB <- dcast(GOterm_DB[,c("Word", 'Degree', 'FreqN')], Word ~ Degree,  value.var = "FreqN")
M_GOterm_DB[is.na(M_GOterm_DB)] <- 0
row.names(M_GOterm_DB) <- M_GOterm_DB[,1]
M_GOterm_DB <- M_GOterm_DB[,-c(1)]
M_GOterm_DB <- as.matrix(M_GOterm_DB)
########

################################################

####################################################
##          CLR Enrichment for full TF-GO net     ##
####################################################

## Define list of targets and genes by GO input
#
GOList1 <- split(GO.listV2$GeneID, GO.listV2$GO.ID)
GOList2 <- split(GO.listV2$GeneID, GO.listV2$pGO.ID)

GOList <- c(GOList1, GOList2)
GOList <- GOList[unique(names(GOList))]

length(GOList)
#GOListp <- split(GO.listV2$GeneID, GO.listV2$pGO.ID)

length(unique(GO_Networkv3$TF))
length(unique(GO_Networkv3$Ann.ID))

#
FullNet_list <- split(FullNet$Target, FullNet$Source)
MRMI_full <- split(MRMI_full$Target, MRMI_full$TF)

# 
GetCLR_enrichment <- function(tf, TargList, GO_vector , GOList) {
  # Tf targets
  targets <- TargList[[tf]]
  
  # target: p
  p <- length(targets)
  
  # Genes in GO within network
  golist <- unique(GO_vector)
  
  # output DF
  DFout <- tibble("TF"=tf, GO=golist, "p"=p, "c"=NA, "t"=NA, E=NA, Log2E=NA)
  
  for (go in golist) {
    # go = golist[1]
    # genes in go
    g_go <- unique(GOList[[go]])
    
    # Total genes: u
    u <- length(unique(c(g_go, targets)))
    
    # Intersection: c
    c <- length(intersect(targets, g_go))
    
    # Total GOs: t
    t <- length(g_go)
    
    # Enrichment CLR
    E <- (c/t)/(p/u)
    
    # Save results
    DFout[DFout$GO == go,]$E <- E
    DFout[DFout$GO == go,]$c <- c
    DFout[DFout$GO == go,]$t <- t
  }
  
  DFout$Log2E <- log2(DFout$E)
  
  # rm(targets)
  # rm(tf)
  # rm(golist)
  # rm(p)
  # rm(DFout)
  # rm(go)
  # rm(g_go)
  # rm(u)
  # rm(c)
  # rm(t)
  # rm(E)
  # rm(DFout)
  
  return(DFout)

  }


# Enrichment with GOs
Test_DF <- GetCLR_enrichment(unique(GO_Networkv3$TF)[1], 
                                  MRMI_full,
                                  unique(GO_Networkv3$Ann.ID),
                                  GOList)

# Test enrichment with target from MI net
EnrichmentDF <-  lapply(unique(GO_Networkv3$TF), 
                        function(x) GetCLR_enrichment(x, 
                                                      MRMI_full,
                                                      unique(GO_Networkv3$Ann.ID),
                                                      GOList))
# Test enrichment with target from full net
EnrichmentDFt <-  lapply(unique(GO_Networkv3$TF), 
                        function(x) GetCLR_enrichment(x, 
                                                      FullNet_list,
                                                      unique(GO_Networkv3$Ann.ID),
                                                      GOList))


# With full list of target genes
EnrichmentDF <- rbindlist(EnrichmentDF, idcol = F)
EnrichmentDFt <- rbindlist(EnrichmentDFt, idcol = F)

# Keep only interactions already predicted
EnrichmentDFv2  <- EnrichmentDF[paste0(EnrichmentDF$TF, EnrichmentDF$GO) %in% 
                               paste0(GO_Networkv3$TF, GO_Networkv3$Ann.ID), ]

EnrichmentDFtv2 <- EnrichmentDFt[paste0(EnrichmentDFt$TF, EnrichmentDFt$GO) %in% 
                                paste0(GO_Networkv3$TF, GO_Networkv3$Ann.ID), ]


hist(EnrichmentDFv2$Log2E, 100)
hist(EnrichmentDFtv2$Log2E, 100)

EnrichmentDFv2$Log2E[EnrichmentDFv2$Log2E <= -9] <- 9
EnrichmentDFtv2$Log2E[EnrichmentDFtv2$Log2E <= -9] <- 9

# Z-score by TF
EnrichmentDFv2 %>%
  filter(c>0) %>%
  group_by(TF) %>%
  mutate( TFz = (Log2E - mean(Log2E))/sd(Log2E)) -> EnrichmentDFv2

EnrichmentDFtv2 %>%
  filter(c>0) %>%
  group_by(TF) %>%
  mutate( TFz = (Log2E - mean(Log2E))/sd(Log2E)) -> EnrichmentDFtv2

# Z-score by GO
EnrichmentDFv2 %>%
  filter(c>0) %>%
  group_by(GO) %>%
  mutate(GOz = (Log2E - mean(Log2E))/sd(Log2E)) -> EnrichmentDFv2

EnrichmentDFtv2 %>%
  filter(c>0) %>%
  group_by(GO) %>%
  mutate(GOz = (Log2E - mean(Log2E))/sd(Log2E)) -> EnrichmentDFtv2

#table(as.data.table(table(EnrichmentDFv2[is.na(EnrichmentDFv2$GOz),]$GO))$N)
#table(as.data.table(table(EnrichmentDFtv2[is.na(EnrichmentDFtv2$GOz),]$GO))$N)
#table(as.data.table(table(EnrichmentDFtv2[is.na(EnrichmentDFtv2$TFz),]$TF))$N)

EnrichmentDFv2$GOz[is.na(EnrichmentDFv2$GOz)] <- 1

EnrichmentDFtv2$GOz[is.na(EnrichmentDFtv2$GOz)] <- 1
EnrichmentDFtv2$TFz[is.na(EnrichmentDFtv2$TFz)] <- 1

# Remove GOs with large number of genes annotated
EnrichmentDFv2 <- left_join(EnrichmentDFv2,  GO.list.parent$DF[, c('go', "term")], by=c("GO"='go')) 

EnrichmentDFtv2 <- left_join(EnrichmentDFtv2,  GO.list.parent$DF[, c('go', "term")], by=c("GO"='go')) 

# Manual curatoin
GOsTooGeneral <- unique(subset(EnrichmentDFv2, t >= 500)[,c("GO",'t','term')])
GOsTooGeneral <- c('GO:0006351', 'GO:0006355', "GO:0032501", 'GO:0050794', 'GO:0016043',
                   "GO:0007275", "GO:0065007", "GO:0051704", "GO:0009987", "GO:0010467",
                   'GO:0044249', 'GO:0006412', 'GO:0006508', 'GO:0006810', 'GO:0023052', 
                   'GO:0016310',' GO:0009058', 'GO:0006612')

View(GOsTooGeneral)
EnrichmentDFv2 <- subset(EnrichmentDFv2, !(GO %in% GOsTooGeneral))
#EnrichmentDFtv2 <- subset(EnrichmentDFtv2, t <= 500)
table(subset(EnrichmentDFv2, t >= 500)$term)

length(unique(EnrichmentDFv2$TF))
length(unique(EnrichmentDFv2$GO))

#plot(EnrichmentDF$TFz ~ EnrichmentDF$GOz)
EnrichmentDFv2[,"rZ"] <- sqrt(sapply(EnrichmentDFv2$TFz, function(x) max(0, x))^2 + 
                              sapply(EnrichmentDFv2$GOz, function(x) max(0, x))^2)

# Add p-value
EnrichmentDFv2[,'Pval'] <- sapply(EnrichmentDFv2$rZ, function (x) pnorm(x, lower.tail=FALSE))

# Add FDR by TF 
EnrichmentDFv2 %>%
  group_by(TF) %>%
  mutate(FDR = p.adjust(Pval, method = 'fdr')) -> EnrichmentDFv2

# scale rZ values for ploting
#EnrichmentDFv2[,'rZs'] <- scale(EnrichmentDFv2$rZ)

write.table(EnrichmentDFv2, 'TFGO_4337_net.txt', row.names = F, quote = F, sep = '\t')
saveRDS(EnrichmentDFv2, 'TFGO_4337_net.rds')



#####################
# Enrichment matrix
#####################
Enr_M <- dcast(EnrichmentDFv2[,c("TF","GO", "rZ")], TF ~ GO,  value.var = 'rZ')
row.names(Enr_M) <- Enr_M$TF

#
Enr_M <- as.matrix(Enr_M[,-c(1)])
Enr_M[is.na(Enr_M)] <- 0
dim(Enr_M)

# GO mask first: GO with at least 2 TFs
dim(Enr_M)
GOsMask <- colSums((Enr_M > 0)*1) >= 4
table(GOsMask)

Enr_M_mask <- Enr_M[,GOsMask]
dim(Enr_M_mask)

# Mask TF: kepp TFs with at least 2 GOs
GenesMask <- rowSums((Enr_M_mask > 0)*1) >= 4
table(GenesMask)

Enr_M_mask <- Enr_M_mask[GenesMask,]

Enr_M_mask[1:10,1:10]
## Test cluster by gos

# cluster by genes
fviz_nbclust(Enr_M_mask, kmeans, method = "silhouette") + 
  labs(subtitle = "Silhouette method")

# cluster by GOs
fviz_nbclust(t(Enr_M_mask), kmeans, method = "silhouette") + 
  labs(subtitle = "Silhouette method")


#my_palette <- colorRamp2(c(0, 0.0001, 15), c("white", "mediumseagreen", "darkorchid1"))

# Replace tf id by Tf name
#colnames(Matrix.GO$TargGO) <- ReplaceNameVector(colnames(rZ))

# TF names
TFdic <- as_tibble(read.table("Data/Annotations/TF_Id_Name.txt", h=F, stringsAsFactors = F))

# cluster TFs and color dendogram 
col_dend = as.dendrogram(hclust(dist(t(as.matrix(rZ)))))
col_dend = color_branches(col_dend, k = 3, col=c("blueviolet", "darkorange1", "forestgreen")) # `color_branches()` returns a dendrogram object

rownames(rZ_mask) <- ReplaceName(rownames(Enr_M_mask))
colnames(rZ_mask) <- ReplaceGO(colnames(Enr_M_mask))


my_palette <- colorRamp2(seq(min(Enr_M_mask), 
                             max(Enr_M_mask), 
                             0.2), 
                         viridis(length(
                           seq(min(Enr_M_mask), 
                               max(Enr_M_mask), 
                               0.2)), 
                           direction = 1, option = "D"))

Hetmap_Z <- Heatmap(Enr_M_mask, name="rZ",
                       column_km = 3,
                       row_km =  3,
                       # column_names_rot = 90,
                       # row_names_rot = 45,
                       cluster_rows = TRUE, 
                       cluster_columns = TRUE,
                       show_column_dend = F, show_row_dend = F, 
                       clustering_method_columns = "ward.D2",
                       col=my_palette,
                       column_names_gp = gpar(fontsize = 1),
                       row_names_gp = gpar(fontsize = 1),
                       show_heatmap_legend = T,
                       # heatmap_ = unit(10),
                       heatmap_height = unit(15, 'cm'),
                       heatmap_width  = unit(12, 'cm'),
                       heatmap_legend_param = list(title_gp = gpar(fontsize = 10), 
                                                   labels_gp = gpar(fontsize = 10), 
                                                   direction = "horizontal"))

#Hetmap_Z
draw(Hetmap_Z, heatmap_legend_side = "bottom")
#####################


####################################################
##            Examples and predictions            ##
####################################################

# TF names
TFdic <- as_tibble(read.table("Data/Annotations/TF_Id_Name.txt", h=F, stringsAsFactors = F))

as.data.table(table(EnrichmentDFv2[,c("GO", "term")])) %>%
  filter(N>0) %>%
  arrange(-N) -> GO_Candidates
write.table(GO_Candidates, 'GOs_InDegree.txt', row.names = F, quote = F, sep = '\t')

Chop=function(myStr,mySep,myField){
  ## chop a string by a separator and return specified field
  choppedString=sapply(strsplit(myStr,mySep),"[",myField)
  if(length(myField)>1){
    choppedString=apply(choppedString,2,function(x){paste0(x[!is.na(x)],collapse=mySep)})
  }
  return(choppedString)
}

Get_Net_relatedToGO <- function(term){
  #term = "fatty"
  #tfs <- subset(EnrichmentDFv2, GO == go)$TF
  #unique(EnrichmentDFv2$term[grepl(term, EnrichmentDFv2$term)])
  go <- unique(EnrichmentDFv2$GO[grepl(term, EnrichmentDFv2$term)] )
  
  subtieble <- paste0(Chop(term, '[|]', 1) , '-related') # Used only first/main descriptior
  
  net <- subset(EnrichmentDFv2, GO %in% go)
  net[, 'TFname'] <- ReplaceName(net$TF)
  
  # Make plot
  # Top predictions
  ggplot(net, aes(x=TFz, y=GOz, color=rZ)) +
    geom_point() +
    geom_text_repel(data=subset(net, rZ>=0.5),
                    aes(x=TFz, y=GOz, label=TFname, color=rZ),
                    size=0.8, segment.size=0.2, 
                    point.padding = 0, 
                    box.padding = 0.5,
                    max.overlaps = Inf) +
    ylab(bquote(GO[Z~";TFs per GO"])) + 
    xlab(bquote(TF[Z~";GOs per TF"])) +
    scale_color_viridis(labels=comma, option = 'B', direction = -1)  + 
    theme_pubclean() + 
    labs(subtitle = subtieble) +
    theme(strip.text.x = element_text(size = 10), 
          axis.text=element_text(size=10),
          legend.position = 'bottom',
          legend.direction='horizontal', 
          legend.text = element_text(angle = 90, hjust=0.1, vjust = 0.5),
          legend.key.size = unit(0.3, 'cm'),
          plot.subtitle=element_text(size=6, hjust=0, color="black"),
          text = element_text(size=10)) -> Plot_out
  
  
  return(Plot_out)
  
}


# length(unique(EnrichmentDFv2$TF))
# length(unique(EnrichmentDFv2$GO))

(nrow(subset(EnrichmentDFv2, TFz >= 1 & GOz >=1))/nrow(EnrichmentDFv2))*100

table(subset(EnrichmentDFv2, TFz >= 1 & GOz >=1)$GO)

####################################################

####################################################
##        TF OutDegree vs TF-TF nets              ##
####################################################

Get_TFnet_relatedToGO <- function(term){
  #
  go <- unique(EnrichmentDFv2$GO[grepl(term, EnrichmentDFv2$term)] )
  
  tfs <- subset(EnrichmentDFv2, GO %in% go)$TF
  net <- subset(MRMI_TFnet, TF1 %in% tfs & TF2 %in% tfs)
  net[, 'TF1name'] <- ReplaceName(net$TF1)
  net[, 'TF2name'] <- ReplaceName(net$TF2)
  
  net[,'Egde'] <- 1
  
  degree1 <- as.data.table(table(net$TF1name)) %>%
    arrange(-N) %>% top_n(5)
  
  
  degree2 <- as.data.table(table(net$TF2name)) %>%
    arrange(-N) 
  
  
  net$TF1name <- factor(net$TF1name, levels = rev(degree1$V1))
  net$TF2name <- factor(net$TF2name, levels = rev(degree2$V1))
  
  Top <- as.data.table(table(net$TF1name)) %>%
    arrange(-N) %>% top_n(5)
  
  
  net <- subset(net, TF1name %in% Top$V1)
  # Make plot
  
  ggplot(net, aes(y=TF1name, x=TF2name, fill=Egde)) +
    geom_tile(fill="#FF69B4") +
    theme_bw() + 
    scale_y_discrete(position = "right") + 
    ylab('') + xlab("") +
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(), 
          axis.text = element_text(size = 10), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=4),
          plot.title = element_text(size = 10),
          legend.position = 'none',
          text = element_text(size=10, family="Times")) -> Plot_out
  Plot_out
  
  return(Plot_out)
  
}

Get_TFnet_relatedToGOv2 <- function(term){
  # term='brassinosteroid'
  # GOs id for the corresponding with term
  go <- unique(EnrichmentDFv2$GO[grepl(term, EnrichmentDFv2$term)] )
  
  #  TFs link to the tested GOs
  tfs <- subset(EnrichmentDFv2, GO %in% go)$TF
  
  # net of TF2s associated with TF1s:
  net <- subset(FullNet, Source %in% tfs & Target %in% tfs)
  colnames(net)[1:2] <- c('TF1', 'TF2')
  
  if (nrow(net) >0) {
    # TFs targeted by TF
    degree1 <- as.data.table(table(net$TF1)) %>%
      arrange(-N) 
    colnames(degree1) <- c('TF', 'TFprocess')
    
    # Total TF target by TF
    degree1 <- left_join(degree1, TF_TF_outDegree, by=c('TF'='V1'))
    colnames(degree1)[3] <- 'TFtotal'
    
    # Master regulator score
    degree1[,'MasterTFscore'] <- (((degree1$TFprocess/degree1$TFtotal)*(degree1$TFprocess/nrow(degree1)))/nrow(degree1))*100
    
    #Master regulator rank
    degree1[,'RankMTFS'] <- rank(-degree1$MasterTFscore)
    
    # Add TF name
    degree1[, 'TFname'] <- ReplaceName(degree1$TF)
    
    degree1 %>% 
      arrange(RankMTFS) -> degree1
    
    degree1[,'Term'] <- paste0(Chop(term, '[|]', 1) , '-related')
    
    return(degree1)
    
  } else {
    print(paste0('.. Salado ..: ', term,' '))
  }
    
  
}

# Calculate degree in TF-Tf network
TF_TF_outDegree <- as.data.table(table(subset(FullNet, Target %in% All_TFs)$Target))

# Calculate master-regulator score
MTFprosess <- list(Get_TFnet_relatedToGOv2("abscisic"),
                   Get_TFnet_relatedToGOv2("auxin"),
                   Get_TFnet_relatedToGOv2("brassinosteroid"),
                   Get_TFnet_relatedToGOv2("cytokinin"),
                   Get_TFnet_relatedToGOv2("salicylic"),
                   Get_TFnet_relatedToGOv2("jasmonic"),
                   Get_TFnet_relatedToGOv2("ethylene"),
                   Get_TFnet_relatedToGOv2("gibberellin"),
                   Get_TFnet_relatedToGOv2("lipid|fatty|tricarboxylic|oxylipin|galactolipid|tricarboxylic"),
                   Get_TFnet_relatedToGOv2('phenylpropanoid|flavonoid'),
                   Get_TFnet_relatedToGOv2('nitrogen|nitrate'),
                   Get_TFnet_relatedToGOv2('cell wall|cellulose|xylan|pectin'),
                   Get_TFnet_relatedToGOv2('pentose|trehalose'),
                   Get_TFnet_relatedToGOv2('sugar|glucose|carbon|pentose|trehalose|oxalate|starch|carbohydrate|polysaccharide'),
                   Get_TFnet_relatedToGOv2('carotenoid|carotene|sterol'),
                   Get_TFnet_relatedToGOv2('leaf'),
                   Get_TFnet_relatedToGOv2('shoot|meristem|cotyledon'),
                   Get_TFnet_relatedToGOv2('pollen|anther|carpel'),
                   Get_TFnet_relatedToGOv2('flower'),
                   Get_TFnet_relatedToGOv2('seed|embryo|embryonic'))

# Combined predictions
MTFprosess <- rbindlist(MTFprosess[unlist(lapply(MTFprosess, function(x) is.data.frame(x)))], idcol = T)


# Create term index
MTFprosess[,'id2'] <- paste0('[',MTFprosess$.id,'] ', MTFprosess$Term)

# defined factors
MTFprosess$id2 <- factor(MTFprosess$id2, 
                          levels = c('[1] abscisic-related',
                                     '[2] auxin-related',
                                     '[3] brassinosteroid-related',
                                     '[4] cytokinin-related',
                                     '[5] salicylic-related',
                                     '[6] jasmonic-related',
                                     '[7] ethylene-related',
                                     '[8] gibberellin-related',
                                     '[9] lipid-related',
                                     '[10] phenylpropanoid-related',
                                     '[11] nitrogen-related',
                                     '[12] cell wall-related',
                                     '[13] pentose-related',
                                     '[14] sugar-related',
                                     '[15] carotenoid-related',
                                     '[16] leaf-related',
                                     '[17] shoot-related',
                                     '[18] pollen-related',
                                     '[19] flower-related',
                                     '[20] seed-related'))

# Create label DF with TFs within the first ranks
MTFprosess_Label <- MTFprosess %>%
  dplyr::group_by(Term) %>% 
  dplyr::filter(RankMTFS <= 2)

# Create label DF with TFs within the first ranks
Get_TFnet_Class_ToGOv2 <- function(term, rankfilter){
  # term='brassinosteroid'
  # GOs id for the corresponding with term
  go <- unique(EnrichmentDFv2$GO[grepl(term, EnrichmentDFv2$term)] )
  
  #  TFs link to the tested GOs
  tfs <- subset(EnrichmentDFv2, GO %in% go)$TF
  
  # net of TF2s associated with TF1s:
  net <- subset(FullNet, Source %in% tfs & Target %in% tfs)
  colnames(net)[1:2] <- c('TF1', 'TF2')
  
  if (nrow(net) >0) {
    # TFs targeted by TF
    degree1 <- as.data.table(table(net$TF1)) %>%
      arrange(-N) 
    colnames(degree1) <- c('TF', 'TFprocess')
    
    # Total TF target by TF
    degree1 <- left_join(degree1, TF_TF_outDegree, by=c('TF'='V1'))
    colnames(degree1)[3] <- 'TFtotal'
    
    # Master regulator score
    degree1[,'MasterTFscore'] <- (((degree1$TFprocess/degree1$TFtotal)*(degree1$TFprocess/nrow(degree1)))/nrow(degree1))*100
    
    #Master regulator rank
    degree1[,'RankMTFS'] <- rank(-degree1$MasterTFscore)
    
    # Add TF name
    degree1[, 'TFname'] <- ReplaceName(degree1$TF)
    degree1 %>% 
      arrange(RankMTFS) -> degree1
  
    degree1[,'Term'] <- paste0(Chop(term, '[|]', 1) , '-related')
    # select top TFs based on rank filter
    tftop <- subset(degree1, RankMTFS <= rankfilter)$TF
    
    # target for specifc TF
    net <- net[net$TF1 %in% tftop,]
    net[,'TF1name'] <- ReplaceName(net$TF1)
    net[,'TF2name'] <- ReplaceName(net$TF2)
    
    # add term label
    net[,'Term'] <- paste0(Chop(term, '[|]', 1) , '-related')
    
    net$TF1name <- factor(net$TF1name, levels = rev(ReplaceName(tftop)))

      
    return(net)
    
  } else {
    print(paste0('.. Salado ..: ', term,' '))
  }
}

# TF names
TFdic <- as_tibble(read.table("Data/Annotations/TF_Id_Name.txt", h=F, stringsAsFactors = F))

MTFNetTop <- list(Get_TFnet_Class_ToGOv2("abscisic", 2),
                  Get_TFnet_Class_ToGOv2("lipid|fatty|tricarboxylic|oxylipin|galactolipid|tricarboxylic", 2),
                  Get_TFnet_Class_ToGOv2('phenylpropanoid|flavonoid', 2),
                  Get_TFnet_Class_ToGOv2('leaf',2))

MTFNetTop <- rbindlist(MTFNetTop)
#subset(FullNet, Source =='Zm00001d015815' & Target %in% MTFNetTop$Target)


as.data.table(table(MTFNetTop[,c(1,6)])) %>%
  filter(N>0) %>%
  mutate(Name=ReplaceName(TF1))


########
OutDegree <- left_join(OutDegree, TF_TF_outDegree, by=c('TF'="V1"))
colnames(OutDegree)[5] <- 'TFsMR'

plot(OutDegree$GOs, OutDegree$N)
OutDegree

# Out degree vs tag per TF
ggplot(OutDegree, aes(y=GOs, x=TFsMR)) + 
  geom_pointdensity() +
  scale_color_viridis(labels=comma) +
  scale_y_log10() +
  scale_x_log10() +
  annotation_logticks(sides = "lb", color = 'black') +
  theme_pubclean() +
  stat_poly_line() +
  stat_poly_eq(color='#FF4500') +
  ylab('GOs per TF') + 
  xlab('TF Targets per TF') +
  theme(axis.text=element_text(size=10), 
        legend.position = 'bottom', 
        legend.text = element_text(angle = 90, hjust=0.1, vjust = 0.5),
        legend.key.size = unit(0.3, 'cm'),
        text = element_text(size=10, family="Times"))
########


####################################################
##              Y1H and predictions               ##
####################################################
View(Top_EnrichmentDF)
View(subset(Top_EnrichmentDF, TF %in% Y1H$TF.v4))

Y1H[Y1H$TF.v4=="Zm00001d024353", ]
subset(Y1H, TF.v4 %in% "Zm00001d007436")
subset(Y1H, TF.v4 %in% "Zm00001d011496")

EnrichmentDFv2_Y1H <- subset(EnrichmentDFv2, TF %in% Y1H$TF.v4)
EnrichmentDFv2_Y1H <- subset(EnrichmentDFtv2, TF %in% Y1H$TF.v4)

EnrichmentDFv2_Y1H[,"TFname"] <- ReplaceName(EnrichmentDFv2_Y1H$TF)
EnrichmentDFv2_Y1H <- left_join(EnrichmentDFv2_Y1H,  GO.list.parent$DF[, c('go', "term")], by=c("GO"='go')) 
View(EnrichmentDFv2_Y1H)

EnrichmentDFv2_Y1H

################################################


##################################################
######              Plots               ##########
##################################################

######
# TF-GO net description: Fig 3a-d
######


# Out degree: GOs per TF
ggplot(OutDegree_Table, aes(x=Outdegree, y=N)) + 
  geom_point(size=1) + 
  scale_y_log10(label=comma, limit=c(1, 1000)) +
  scale_x_log10() +
  annotation_logticks(sides = "lb", color = 'black') +
  theme_pubclean() +
  stat_poly_line() +
  stat_poly_eq() +
  xlab('Out degree\n[GOs per TF]') + 
  ylab('Counts') +
  theme(axis.text=element_text(size=10), 
        legend.position = 'none',
        text = element_text(size=10, family="Times")) -> Plot_3a

# In degree: TFs per GO
ggplot(InDegree_Table, aes(x=Indegree, y=N)) + 
  geom_point(size=1) + 
  scale_y_log10(label=comma, limit=c(1, 1000)) +
  scale_x_log10() +
  annotation_logticks(sides = "lb", color = 'black') +
  theme_pubclean() +
  stat_poly_line() +
  stat_poly_eq() +
  xlab('In degree\n[TFs per GO]') + 
  ylab('Counts') +
  theme(axis.text=element_text(size=10), 
        legend.position = 'none',
        text = element_text(size=10, family="Times")) -> Plot_3b


# In degree vs genes in GO
ggplot(InDegree, aes(y=N, x=GO.IDsize)) + 
  geom_pointdensity() +
  scale_color_viridis(labels=comma) +
  scale_y_log10() +
  scale_x_log10() +
  annotation_logticks(sides = "lb", color = 'black') +
  theme_pubclean() +
  stat_poly_line() +
  stat_poly_eq(color='#FF4500') +
  ylab('TFs per GO') + 
  xlab('Genes per GO') +
  theme(axis.text=element_text(size=10), 
        legend.position = 'bottom', 
        legend.direction='horizontal', 
        legend.text = element_text(angle = 90, hjust=0.1, vjust = 0.5),
        legend.key.size = unit(0.4, 'cm'),
        text = element_text(size=10, family="Times")) -> Plot_3c

# Out degree vs tag per TF
ggplot(OutDegree, aes(y=GOs, x=Targ.MR)) + 
  geom_pointdensity() +
  scale_color_viridis(labels=comma) +
  scale_y_log10() +
  scale_x_log10() +
  annotation_logticks(sides = "lb", color = 'black') +
  theme_pubclean() +
  stat_poly_line() +
  stat_poly_eq(color='#FF4500') +
  ylab('GOs per TF') + 
  xlab('Targets per TF') +
  theme(axis.text=element_text(size=10), 
        legend.position = 'bottom', 
        legend.text = element_text(angle = 90, hjust=0.1, vjust = 0.5),
        legend.key.size = unit(0.4, 'cm'),
        text = element_text(size=10, family="Times")) -> Plot_3d

######

######
# GO word freq by degree 
######

M_GOterm_DB

Hetmap_GOsDegree <- Heatmap(M_GOterm_DB, name="Fraction",
                            # column_km = 3,
                            # column_names_rot = 90,
                            # row_names_rot = 45,
                            cluster_rows = T, 
                            cluster_columns = F,
                            show_column_dend = F, show_row_dend = T, 
                            clustering_method_columns = "ward.D2",
                            col=my_palette,
                            column_names_gp = gpar(fontsize = 6),
                            row_names_gp = gpar(fontsize = 1.5),
                            show_heatmap_legend = T,
                            # heatmap_ = unit(10),
                            heatmap_height = unit(20, 'cm'),
                            heatmap_width  = unit(5, 'cm'),
                            heatmap_legend_param = list(title_gp = gpar(fontsize = 10), 
                                                        labels_gp = gpar(fontsize = 10), 
                                                        direction = "horizontal"))

pdf("Plots/Fig_S12_TermFreq.pdf", width=5, height=10)
draw(Hetmap_GOsDegree, heatmap_legend_side = "bottom")
dev.off()


ggplot(GOterm_DB, aes(x=Degree, y=FreqN, label=Word)) +
  geom_text_repel()




######

######
# Z-score of TFs vs GOs plot: Fig 3e
######

ggplot(EnrichmentDFv2,  aes(y=GOz, x=TFz)) +
  geom_pointdensity(size=0.3) +
  #scale_color_viridis(labels=comma, option = "A")  + 
  scale_color_gradient2(low = "white", high = 'black') + 
  # geom_text_repel(data=Top_EnrichmentDF, 
  #                 aes(y=GOz, x=TFz, label=Label), 
  #                 box.padding = 1, 
  #                 max.overlaps = Inf,
  #                 size=1, segment.size=0.2) + 
  geom_hline(yintercept = 1, colour='#D2691E', linetype='dotted') +
  geom_vline(xintercept = 1, colour='#D2691E', linetype='dotted') +
  theme_pubclean() + 
  ylab(bquote(GO[Z~";TFs per GO"])) + 
  xlab(bquote(TF[Z~";GOs per TF"])) +
  theme(strip.text.x = element_text(size = 10), 
        axis.text=element_text(size=10),
        legend.position = 'bottom',
        legend.direction='horizontal', 
        legend.text = element_text(angle = 90, hjust=0.1, vjust = 0.5),
        legend.key.size = unit(0.3, 'cm'),
        text = element_text(size=10)) -> Plot_3e

Plot_3e

######

######
# Plot 3f 
######
TFdic <- as_tibble(read.table("Data/Annotations/TF_Id_Name.txt", h=F, stringsAsFactors = F))

Plot_3f1 <- Get_Net_relatedToGO("abscisic")
Plot_3f2 <- Get_Net_relatedToGO('lipid|fatty|tricarboxylic|oxylipin|galactolipid')
Plot_3f3 <- Get_Net_relatedToGO('phenylpropanoid|flavonoid')
#Plot_3f4 <- Get_Net_relatedToGO("glucosinolate")
Plot_3f5 <- Get_Net_relatedToGO("leaf")


# plot_layout(guides = "collect")
# Plot_3f <- Plot_3f & theme(legend.position = "bottom")  
# Plot_3f + plot_layout(guides = "collect")

######

######
# Plot 3g
######

# Set labels
TFdic <- as_tibble(read.table("Data/Annotations/TF_Id_Name.txt", h=F, stringsAsFactors = F))
MTFprosess_Label$TFname <- ReplaceName(MTFprosess_Label$TF)
MTFprosess_Label[,'idTF'] <- paste0('[',MTFprosess_Label$.id,'] ', MTFprosess_Label$TFname)

library(RColorBrewer)
getPalette=colorRampPalette(brewer.pal(12, "Paired"))


# MTFprosess$MasterTFscore[MTFprosess$MasterTFscore > 0.2] <- 0.2

# set max rank
MTFprosess$RankMTFS[MTFprosess$RankMTFS >= 10] <- 10

ggplot(MTFprosess, aes(y=MasterTFscore, x=RankMTFS, color=id2)) + 
  geom_point(size=1) +
  scale_x_reverse(limits=c(10, -5)) +
  scale_colour_manual(values=getPalette(20)) + 
  theme_pubclean() +
  xlab('TF Rank') + 
  ylab('Regulator score per process') +
  geom_text_repel(data=MTFprosess_Label,
                  aes(y=MasterTFscore, x=RankMTFS, color=id2, label=idTF),
                  size=1.5, 
                  force = 0.5,
                  segment.size=0.2, 
                  point.padding = 0, 
                  box.padding = 0.5,
                  nudge_x = 0.1,
                  max.overlaps = Inf,
                  direction = "y",
                  xlim=c(1, NA),
                  hjust= 0,
                  show.legend = FALSE) +
  theme(axis.text=element_text(size=10), 
        legend.position = 'bottom', legend.text.align = 0,
        legend.key.size = unit(0.2, 'cm'),
        text = element_text(size=10, family="Times")) -> Plot_3g

Plot_3g

pdf("Plots/Fig_3g.pdf", width = 4, height = 5)
print(Plot_3g)
dev.off()

######

#####
# Fig 3h
#####
MTFNetTop$Term <- factor(MTFNetTop$Term, 
                         levels = c('abscisic-related', 
                                    'lipid-related',
                                    'phenylpropanoid-related',
                                    'leaf-related'))
# Plot
ggplot(MTFNetTop, aes(y=TF1name, x=TF2name, fill=Edge)) +
  geom_tile() +
  theme_bw() + 
  scale_y_discrete(expand = c(0,0)) + 
  scale_fill_brewer(palette="Set3") +
  ylab('') + xlab("") +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text = element_text(size = 10), 
        axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1, size=10),
        plot.title = element_text(size = 10),
        legend.position = 'bottom',
        text = element_text(size=10, family="Times")) +
  facet_wrap( . ~ Term, nrow = 4, ncol = 1, scales = 'free') -> Plot_3h

Plot_3h

pdf("Plots/Fig_3h.pdf", width = 4.5, height = 6)
print(Plot_3h)
dev.off()

#####

########
## Combine all Plots plots
########

# Plots a to d
Plot_3ad <- { Plot_3a + Plot_3b + Plot_3c + Plot_3d + plot_layout(nrow = 1, ncol = 4)} 

pdf("Plots/Fig_3ad.pdf", width = 8, height = 3)
print(Plot_3ad)
dev.off()

# Plot_3f <- {Plot_3e|Plot_3f1|Plot_3f2|Plot_3f3|Plot_3f5}

Plot_3ef <- Plot_3e + Plot_3f1 + Plot_3f2 + Plot_3f3 + Plot_3f5 + plot_layout(nrow = 1, ncol = 5)
pdf("Plots/Fig_3ef.pdf", width = 8, height = 2.5)
print(Plot_3ef)
dev.off()

#
Plot_3g <- Plot_3g1+Plot_3g2+Plot_3g3+Plot_3g4+Plot_3g5 + plot_layout(nrow = 1, ncol = 5)
pdf("Plots/Fig_3g.pdf", width = 9, height = 1.6)
print(Plot_3g)
dev.off()

########

########
#### Supplementary figures 
########

## hormones
Plot_12a1 <- Get_Net_relatedToGO("auxin")
Plot_12a2 <- Get_Net_relatedToGO("brassinosteroid")
Plot_12a3 <- Get_Net_relatedToGO("cytokinin")
Plot_12a4 <- Get_Net_relatedToGO("salicylic")
Plot_12a5 <- Get_Net_relatedToGO("jasmonic")
Plot_12a6 <- Get_Net_relatedToGO("ethylene")
Plot_12a7 <- Get_Net_relatedToGO("gibberelli")

Plot_12a <- (Plot_12a1|Plot_12a2|Plot_12a3|Plot_12a4|Plot_12a5|Plot_12a6|Plot_12a7)  + plot_layout(ncol = 4, nrow = 2)

## Metabolism
Plot_12b1 <- Get_Net_relatedToGO("nitrogen|nitrate")
Plot_12b2 <- Get_Net_relatedToGO("cell wall|cellulose|xylan|pectin")
Plot_12b3 <- Get_Net_relatedToGO("pentose|trehalose")
Plot_12b4 <- Get_Net_relatedToGO("sugar|glucose|carbon|pentose|trehalose|oxalate|starch|carbohydrate|polysaccharide|tricarboxylic")
Plot_12b5 <- Get_Net_relatedToGO("carotenoid|carotene|sterol")

Plot_12b <- (Plot_12b1|Plot_12b2|Plot_12b3|Plot_12b4|Plot_12b5)  + plot_layout(ncol = 4, nrow = 2)

## Development
Plot_12c1 <- Get_Net_relatedToGO("meristem")
Plot_12c2 <- Get_Net_relatedToGO("pollen|anther|carpel")
Plot_12c3 <- Get_Net_relatedToGO("flower")
Plot_12c4 <- Get_Net_relatedToGO("seed|embryo|embryonic")
Plot_12c5 <- Get_Net_relatedToGO("shoot|cotyledon")

Plot_12c <- (Plot_12c1|Plot_12c2|Plot_12c3|Plot_12c4|Plot_12c5)  + plot_layout(ncol = 4, nrow = 2)

Plot_s12 <- Plot_12a/Plot_12b/Plot_12c

######
pdf("Plots/Fig_S12a.pdf", width = 8, height = 5)
print(Plot_12a)
dev.off()

pdf("Plots/Fig_S12b.pdf", width = 8, height = 5)
print(Plot_12b)
dev.off()

pdf("Plots/Fig_S12c.pdf", width = 8, height = 5)
print(Plot_12c)
dev.off()
########

#############


#############
# tem files with df summary
#############

Net_relatedToGO <- function(term){
  
  #term'
  go <- unique(EnrichmentDFv2$GO[grepl(term, EnrichmentDFv2$term)] )
  #
  net <- subset(EnrichmentDFv2, GO %in% go)
  net[, 'TFname'] <- ReplaceName(net$TF)
  
  return(net[,-c(13)])
}

tem1 <- Net_relatedToGO("abscisic")
tem2 <- Net_relatedToGO('lipid|fatty|tricarboxylic|oxylipin|galactolipid')
tem3 <- Net_relatedToGO('phenylpropanoid|flavonoid')
tem4 <- Net_relatedToGO("glucosinolate")
tem5 <- Net_relatedToGO("leaf")


length(unique(tem1$TF))
length(unique(tem2$TF))
length(unique(tem3$TF))
length(unique(tem4$TF))
length(unique(tem5$TF))

tem1

unique(ReplaceName(subset(tem1, rZ>=1)$TF[subset(tem1, rZ>=1)$TF %in% subset(tem2, rZ>=1)$TF]))
unique(ReplaceName(subset(tem1, rZ>=1)$TF[subset(tem1, rZ>=1)$TF %in% subset(tem3, rZ>=1)$TF]))
unique(ReplaceName(subset(tem1, rZ>=1)$TF[subset(tem1, rZ>=1)$TF %in% subset(tem4, rZ>=1)$TF]))
unique(ReplaceName(subset(tem1, rZ>=1)$TF[subset(tem1, rZ>=1)$TF %in% subset(tem5, rZ>=1)$TF]))

unique(ReplaceName(subset(tem2, rZ>=1)$TF[subset(tem2, rZ>=1)$TF %in% subset(tem3, rZ>=1)$TF]))
unique(ReplaceName(subset(tem2, rZ>=1)$TF[subset(tem2, rZ>=1)$TF %in% subset(tem4, rZ>=1)$TF]))
unique(ReplaceName(subset(tem2, rZ>=1)$TF[subset(tem2, rZ>=1)$TF %in% subset(tem5, rZ>=1)$TF]))


length(unique(subset(tem1, rZ>=0.5)$TF))
length(unique(subset(tem2, rZ>=0.5)$TF))
length(unique(subset(tem3, rZ>=0.5)$TF))
length(unique(subset(tem4, rZ>=0.5)$TF))
length(unique(subset(tem5, rZ>=0.5)$TF))


ReplaceName(unique(subset(tem1, rZ>=0.5)$TF))
ReplaceName(unique(subset(tem2, rZ>=0.5)$TF))
ReplaceName(unique(subset(tem3, rZ>=0.5)$TF))
ReplaceName(unique(subset(tem4, rZ>=0.5)$TF))
ReplaceName(unique(subset(tem5, rZ>=0.5)$TF))

ReplaceName(tem3$TF[tem3$TF %in% Y1H$TF.v4])
subset(tem3, TF %in% Y1H$TF.v4)

unique(subset(tem3, rZ>=0.5)$TF) %in% Y1H$TF.v4

tem5 %>% 
  arrange(-rZ)

subset(EnrichmentDFv2, TF=="Zm00001d020430") %>% 
  arrange(-rZ)

subset(EnrichmentDFtv2, TF=="Zm00001d051891") %>% 
  arrange(-E)
#############

#############
########

########
Top_EnrichmentDF <- subset(EnrichmentDFv2, t >=10 & Pval <= 0.01)[,-c(6)] %>%
  dplyr::group_by(TF) %>% 
  dplyr::arrange(TF, -rZ)

Top_EnrichmentDF[,"TFname"] <- ReplaceName(Top_EnrichmentDF$TF)
Top_EnrichmentDF[,'Label'] <- paste0(Top_EnrichmentDF$TFname, ":", Top_EnrichmentDF$term)
########

######