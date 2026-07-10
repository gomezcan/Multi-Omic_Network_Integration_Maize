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
#library(DECIPHER)
# install.packages("UpSetR")


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
    ids <- gsub(TFdic$V2[i], TFdic$V1[i], ids)
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

ReduceGOs <- function(GOsDB, tf){
  #library(rrvgo)
  GOsDB <- unique(subset(GOsDB, TF==tf))
  scores <- setNames(-log10(GOsDB$FDR), GOsDB$GO.ID) 
  GO_vector = GOsDB$GO.ID
  
  if (length(GO_vector) >1) {
    # Semantic similarity
    cat(paste0(tf,'\n'))
    simMatrix <- calculateSimMatrix(GO_vector,  orgdb=org.Zmays.eg.db,  ont="BP", 
                                    semdata=Zm.GOSemSim.BP,
                                    method="Wang")
    
    # Reduce term
    reducedTerms <- reduceSimMatrix(simMatrix, scores, keytype="GENENAME",
                                    threshold=0.7, orgdb=org.Zmays.eg.db)
    
    # treemapPlot(reducedTerms)
    reducedTerms[,"TF"] <- tf
    return(reducedTerms)
  }
  else{
    cat(paste0('.. salado: ', tf,'\n'))
  }
}

Enrichmet_classes <- function(network){
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
  
  #write.table(Pval_table, name, sep = "\t", row.names = F, quote = F)
  print(paste("... Done ...", sep = ""))
  return(Pval_table)
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
  # Define background based on genes in network
  
  #background_tem <- background[names(background) %in% unique(c(Net$Source, Net$Target))]
  #background_IDs_tem <- as.character(unique(names(background_tem)))
  
  # mutant: description of a category of source of the data
  # degs: set of genes to test enrichment
  GeneList <- factor(as.integer(background_IDs %in% degs))
  names(GeneList) <- background_IDs
  
  
  GOdata_BP <- new("topGOdata", ontology = "BP", allGenes = GeneList, 
                   annot = annFUN.gene2GO, gene2GO = background)
  
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
  
  # ##### get all GOs and their genes from the topGO result #####
  # gs <- genesInTerm(GOdata_BP) # list genes by GO
  # # 
  # ANOTATION = lapply(gs,function(x) x[x %in% degs]) ## Get only my Differential expressed genes
  # 
  # ### Get only the GO's located in the result of topGO in Res_DF
  # DF_GO_Genes <- ANOTATION[Res_DF_BP$GO.ID] # list
  # 
  # ## Transform it to a data frame.
  # DF_GO_Genes = list_to_DF(DF_GO_Genes)
  # DF_GO_Genes <- left_join(DF_GO_Genes, Res_DF_BP[,c(1,2,6)], by='GO.ID')   # left join to add GO info
  # DF_GO_Genes["Mutant"] <- mutant # add Mutant column name
  # 
  # DF_GO_Genes <- unique(DF_GO_Genes)
  # 
  # filename <- 
  #   write.table(DF_GO_Genes,
  #               paste("BP_results_targets/Genes_GOBP_", netname, "_", mutant, ".txt", sep = ""),
  #               sep = '\t', quote = F,
  #               row.names = F)
  # 
  return(Res_DF_BP) # Return list of GOs-Stats and GeneID-GOs
}

SuperGO_DEGs <- function(tf){
  ## used TF/Module targets/genes to test GO terms enrichment
  # 1. Select tf/module's Targets
  # 2. Make list file: degs
  # 3. Test enrichment
  
  # Get network by TF
  network <- subset(DEGs_1, TF==tf)
  network <- subset(network, Target %in% Syntenic)
  
  Total_targets <- as_tibble(as.data.frame(table(network$TF), stringsAsFactors = F))
  colnames(Total_targets) <- c('TF', 'nTF') 
  
  # Genes input Net
  degs <- unique(network$Target)
  
  ## Compare list of predicted targets vs annotated genes in query.vector 
  print(". Pre-GO.")
  out <- GetGO(degs, tf)
  #
  colnames(out)[7] <- "TF"
  
  # save GOs
  out <- left_join(out, Total_targets, by="TF")  
  
  # Genes_GOBP_
  return(out)
  
}

Get_PWY_methogd_DEGs <- function(tf){
  # pwy predicted
  pred_pwy <- subset(PWYs, TF==tf)
  #pred_pwy <- subset(PWYs, TF==DEGs_Peng[6])
  pred_pwy <- split(pred_pwy$Annotation, pred_pwy$Type)
  
  # pwy obsedverd in degs
  obs_pwy <- subset(DEGs_1_PWY, TF==tf)
  
  #obs_pwy <- subset(DEGs_1_PWY, TF==DEGs_Peng[6])
  obs_pwy <- split(obs_pwy$PWY, obs_pwy$Treatment)
  
  if(length(pred_pwy) > 1 ){
    ## total pwys in lists
    
    # predicted
    pred_pwy_total <- as.data.frame(unlist(lapply(pred_pwy, length))) 
    colnames(pred_pwy_total)  <- 'n.Method'
    pred_pwy_total[,"Method"] <- row.names(pred_pwy_total)
    pred_pwy_total$Method <- as.character(pred_pwy_total$Method)
    
    # observed
    obs_pwy_total <- as.data.frame(unlist(lapply(obs_pwy, length))) 
    colnames(obs_pwy_total)  <- 'n.DEG'
    obs_pwy_total[,"DEG"] <- row.names(obs_pwy_total)
    obs_pwy_total$DEG <- as.character(obs_pwy_total$DEG)
    
    #
    #print(". Pre-newGOM .")
    pwy.obj <- newGOM(pred_pwy, obs_pwy, genome.size=length(CornCYC.list))
    
    Pval <- getMatrix(pwy.obj, name="pval")
    Common <- getMatrix(pwy.obj, name="intersection")
    
    ### Summary tables
    ## adjust p value
    Pval <- as.data.frame(Pval)
    #Pval[,1:ncol(Pval)] <- apply(Pval[,1:ncol(Pval)], 2, p.adjust)
    Pval_table   <- as_tibble(reshape2::melt(as.matrix(Pval)))
    colnames(Pval_table) <- c('Method', 'DEG', 'Pval')
    
    #
    Common_table <- as_tibble(reshape2::melt(as.matrix(Common)))
    colnames(Common_table) <- c('Method', 'DEG', 'n.Common')
    
    # Add predicted target in class by TF
    Pval_table <- left_join(Pval_table, Common_table , by=c('Method', 'DEG'))
    Pval_table$Method <- as.character(Pval_table$Method)
    Pval_table$DEG <- as.character(Pval_table$DEG)
    
    ## Add total predicted values
    # predicted
    Pval_table <- left_join(Pval_table, pred_pwy_total, by="Method")
    # observed
    Pval_table <- left_join(Pval_table, obs_pwy_total, by="DEG")
    
    #write.table(Pval_table, name, sep = "\t", row.names = F, quote = F)
    cat(paste("... Done ...", tf, '\n', sep = ""))
    return(Pval_table)
  }
  
  if(length(pred_pwy) == 1 ){
    
    # Create empty DF with nrows equal to number of experiments in DEGs: length(obs_pwy)
    out <- as.data.frame(matrix(0, nrow = length(obs_pwy), ncol = 4))
    colnames(out) <- c("Method", "DEG", "Pval", "n.Common")
    
    # Method equal to names in pred_pwys
    out$Method <- as.character(names(pred_pwy))
    
    # Method equal to names in obs_pwy
    out$DEG <- as.character(names(obs_pwy))
    
    ## total pwys in lists
    # predicted
    pred_pwy_total <- as.data.frame(unlist(lapply(pred_pwy, length))) 
    colnames(pred_pwy_total)  <- 'n.Method'
    pred_pwy_total[,"Method"] <- row.names(pred_pwy_total)
    pred_pwy_total$Method <- as.character(pred_pwy_total$Method)
    
    # observed
    obs_pwy_total <- as.data.frame(unlist(lapply(obs_pwy, length))) 
    colnames(obs_pwy_total)  <- 'n.DEG'
    obs_pwy_total[,"DEG"] <- row.names(obs_pwy_total)
    obs_pwy_total$DEG <- as.character(obs_pwy_total$DEG)
    
    #
    #print(". Pre-newGOM .")
    pwy.obj <- newGOM(pred_pwy, obs_pwy, genome.size=length(CornCYC.list))
    
    Pval <- getMatrix(pwy.obj, name="pval")
    Common <- getMatrix(pwy.obj, name="intersection")
    
    
    ### Summary tables
    ## P-value
    Pval <- as.numeric(Pval)
    # Common PWYs
    Common <- as.numeric(Common)
    
    # values to out table
    out$Pval <- Pval
    out$n.Common <- Common
    
    ## Add total predicted values
    # predicted
    out <- left_join(out, pred_pwy_total, by="Method")
    # observed
    out <- left_join(out, obs_pwy_total, by="DEG")
    
    #write.table(Pval_table, name, sep = "\t", row.names = F, quote = F)
    cat(paste("... Done ...", tf, '\n', sep = ""))
    return(out)
  }
  else(return(cat(paste("... Salado ...", tf, '\n', sep = ""))))
  
}

Get_GSS_method_DEGs <- function(tf){
  
  out <- as.data.table(matrix(0, nrow = 0, ncol = 4)) 
  colnames(out) <- c("GO_obs", "GO_pred", "GSS", 'Class')
  
  # make df to save output
  # DEGs
  go_obs <- subset(DEGs_1_GOs, TF == tf)$GO.ID
  
  # define tf id
  tfid <- chop(tf, '[:]', 2)
  
  # ComTarget
  go_ComTarget <- subset(GO_CommTarg, TF == tfid)$GO.ID
  
  # ComFunct
  go_ComFunt <- subset(GO_CommFunt, TF == tfid)$GO.ID
  
  # Network
  go_Network <- subset(GO_Network, TF == tfid)$GO.ID
  
  print(cat(paste0(ReplaceName(tfid), ' ', 
             length(go_ComTarget), ' ',
             length(go_ComFunt), ' ',
             length(go_Network), ' ')))
  
  if(length(go_ComTarget)>1){
    # Matrix of GOs
    M_ss <- mgoSim(go_obs, 
                   go_ComTarget, 
                   semData=Zm.GOSemSim.BP, measure="Wang", combine=NULL)
    M_ss <- as_tibble(melt(M_ss))
    colnames(M_ss) <- c("GO_obs", "GO_pred", "GSS")
    M_ss[,'Class'] <- 'Comm.Target'
    out <- rbind(out, M_ss)
  }
  
  if(length(go_ComFunt)>1){
    # Matrix of GOs
    M_ss2 <- mgoSim(go_obs, 
                    go_ComFunt, 
                    semData=Zm.GOSemSim.BP, measure="Wang", combine=NULL)
    M_ss2 <- as_tibble(melt(M_ss2))
    colnames(M_ss2) <- c("GO_obs", "GO_pred", "GSS")
    M_ss2[,'Class'] <- 'Comm.Funct'
    
    out <- rbind(out, M_ss2)
    
  }
  
  if(length(go_Network)>1){
    # Matrix of GOs
    M_ss3 <- mgoSim(go_obs, 
                    go_Network, 
                    semData=Zm.GOSemSim.BP, measure="Wang", combine=NULL)
    M_ss3 <- as_tibble(melt(M_ss3))
    colnames(M_ss3) <- c("GO_obs", "GO_pred", "GSS")
    M_ss3[,'Class'] <- 'Network-based'
    
    out <- rbind(out, M_ss3)
  }
  
  out[,'TF'] <- tf
  
  #BP_GRN_test <- subset(GOsDB_GRN, Mutant==tf & GO.ID %in% M_ss$GO_GRN)
  #BP_CEN_test <- subset(GOsDB_CEN, Mutant==tf & GO.ID %in% M_ss$GO_CEN)
  
  #M_ss <- left_join(M_ss, BP_GRN_test[,c(1,2,4,9)], by=c("GO_GRN"="GO.ID"))
  #M_ss <- left_join(M_ss, BP_CEN_test[,c(1,2,4,9)], by=c("GO_CEN"="GO.ID"))
  #colnames(M_ss) <- c("GO.GRN", "GO.CEN", "GSS", 
  #                    "Term.GRN",  "Sig.Terms.GRN", "FDR.GRN",
  #                    "Term.CEN", "Sig.Terms.CEN", "FDR.CEN")
  
  #namef <- paste0("GO_SS_data/GSS_",tf,".GRN_CEN.txt")
  #cat(namef)
  #cat('\n')
  #write.table(M_ss, namef, row.names = F, sep = '\t', quote = F)
  #print(cat(" .. Done TFs ..\n"))
  #c = c+1
  
  return(out) 
}


##################################################
##########        Annotations       ##############
##################################################
## Syntenic genes 
Syntenic <- as_tibble(read.table("../Data/Annotations/Zm.v4.synteny.genes.txt", h=T, stringsAsFactors = F, sep = '\t'))
Syntenic <- subset(Syntenic, gene_synteny=='syntenic')$gene_id
length(Syntenic)
#### top 45 ###
#Top45 <- as_tibble(read.table("Data/Annotations/Top45.txt", h=F, stringsAsFactors = F))

# TF names
TFdic <- as_tibble(read.table("../Data/Annotations/TF_Id_Name.txt", h=F, stringsAsFactors = F))

# Phenolic related genes
#PheGenes <- as_tibble(read.table("Data/Annotations/LinaPheGenes2020.txt", h=T, sep = "\t", quote="", stringsAsFactors = F))

# TF and CoReg
TF_CoR <- as_tibble(read.table("../Data/Annotations/TF_CoR_Mazie.txt", h=T, stringsAsFactors = F))
TF_CoR$Family <- gsub('Others', 'Other', TF_CoR$Family)

table(TF_CoR$Family)

## Y1H network
#Y1H <- as_tibble(read.table("Data/Annotations/Y1H.data.txt", h=T, stringsAsFactors = F))[,4:6]
#ReplaceName(Y1H$TF.v4)

# CornCYC
CornCYC <- as_tibble(read.table("../Data/Annotations/corn_pathways.0210325.reduced.txt", h=T, stringsAsFactors = F))
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

# Pre-calculate semantic similarity
Zm.GOSemSim.BP <- GOSemSim::godata(org.Zmays.eg.db, keytype="GENENAME", ont = 'BP')
typeof(org.Zmays.eg.db)

# GOs term annotations
background <- readMappings("Data/Annotations/synteny.ID_TopGO_V4_GRAMER.txt")
background_IDs <- as.character(unique(names(background)))

# GO list 
GO.list <- rbindlist(lapply(background, as_tibble), idcol = T) 
GO.list <- split(GO.list$.id, GO.list$value)

# Files from GSS calculation
## 1. Calculate semantic similarity background
suppressMessages(library(GOSemSim))
suppressMessages(library(AnnotationHub))


# Pre-calculate semantic similarity
Zm.GOSemSim.BP <- GOSemSim::godata(org.Zmays.eg.db, keytype="GENENAME", ont = 'BP')
typeof(org.Zmays.eg.db)
################################################

################################################
##        Read enrichment results: PWYs       ##
################################################

# Common target
PWY_CommTarg <- fread("CommonTarg_PWY_enrichment.txt")
PWY_CommTarg <- subset(PWY_CommTarg, Pval <=0.05)

# Common function
PWY_CommFunt <- fread("CommonFunction_PWY_enrichment.txt")
PWY_CommFunt <- subset(PWY_CommFunt, Pval <=0.05)

# Network-based
PWY_Network <- fread("NetworkBased_PWY_Clusters_enrichment.txt")
PWY_Network <- subset(PWY_Network, Pval <=0.05)

#TF_Network_annotation <- fread("Network_TFs_Incluster_annotation.txt")
#

################################################

################################################
##        Read enrichment results: GOs        ##
################################################

GO_CommTarg <- fread("CommonTarg_GO_enrichment.txt") %>%
  group_by(TF) %>%
  mutate(FDR=p.adjust(classic, method = 'fdr')) %>%
  filter(FDR <= 0.1)

# Common function based on GOs
GO_CommFunt <- fread("CommonFunction_GO_enrichment.txt") %>%
  filter(FDR.1 <= 0.1 & FDR.2 <=0.1 & GSS >= 0.6)  

# Extract each pair of GOs that pass GSS filter
## Part 1
GO_CommFunt %>%
  dplyr::select(GO1, FDR.1, TF) -> GO_CommFunt1
colnames(GO_CommFunt1) <- c('GO.ID',"FDR", "TF")

## Part 2
GO_CommFunt %>%
  dplyr::select(GO2, FDR.2, TF) -> GO_CommFunt2 
colnames(GO_CommFunt2) <- c('GO.ID',"FDR", "TF")

# Combined GOs 
GO_CommFunt <- unique(rbind(GO_CommFunt1, GO_CommFunt2))

#
GO_Network <- fread("NetworkBased_GO_Clusters_enrichment.txt") %>%
  group_by(TF) %>%
  mutate(FDR=p.adjust(classic, method = 'fdr')) %>%
  filter(FDR <= 0.1)
GO_Network

####
## Map GOs to parents
####
GO_CommTarg_TFs <- unique(GO_CommTarg$TF)
GO_CommTarg_Red <- lapply(GO_CommTarg_TFs, function(x) ReduceGOs(GO_CommTarg, x)) 
GO_CommTarg_Red <- rbindlist(GO_CommTarg_Red, idcol = F)
#
table(GO_CommFunt$TF)
GO_CommFunt_TFs <- unique(GO_CommFunt$TF)
GO_CommFunt_Red <- lapply(GO_CommFunt_TFs, function(x) ReduceGOs(GO_CommFunt, x))
GO_CommFunt_Red <- rbindlist(GO_CommFunt_Red, idcol = F)


GO_Network_TFs <- as.character(unique(GO_Network$TF))
GO_Network$TF  <- as.character(GO_Network$TF)

GO_Network_Red <- lapply(GO_Network_TFs, function(x) ReduceGOs(GO_Network, x))
mask <- unlist(lapply(GO_Network_Red, function(x) is.data.frame(x)))
GO_Network_Red <- GO_Network_Red[mask]
GO_Network_Red <- rbindlist(GO_Network_Red, idcol = F)

length(unique(unique(GO_Network_Red[,c("parent", "TF")])$TF))

GONAMES_DB <- unique(rbind(GO_CommTarg_Red[,c('go', "term")] %>% dplyr::rename(GO.ID = go, GO.name=term),
                    GO_CommTarg_Red[,c('parent', "parentTerm")] %>% dplyr::rename(GO.ID = parent, GO.name=parentTerm),
                    GO_CommFunt_Red[,c('go', "term")] %>% dplyr::rename(GO.ID = go, GO.name=term),
                    GO_CommFunt_Red[,c('parent', "parentTerm")] %>% dplyr::rename(GO.ID = parent, GO.name=parentTerm),
                    GO_Network_Red[,c('go', "term")] %>% dplyr::rename(GO.ID = go, GO.name=term),
                    GO_Network_Red[,c('parent', "parentTerm")] %>% dplyr::rename(GO.ID = parent, GO.name=parentTerm)))

GONAMES_DB

# save Reduce GO terms by TF and method
library(tidyverse)



write.table(rbind(GO_CommTarg_Red %>% dplyr::mutate(Method='Com.Target'),
                  GO_CommFunt_Red %>% dplyr::mutate(Method='Com.Function'), 
                  GO_Network_Red %>% dplyr::mutate(Method='Network-base')) %>%
                  dplyr::select(parent, parentTerm, TF, Method) %>%
                  unique(),
            'ReduceGOterms_All_methods.txt', 
            sep = '\t', row.names = F, col.names = T, quote = F)



################################################

########################################################
##   number of genes in network-approach comparison   ##
########################################################

## network used on ClusterONE
#Clusters <- fread('../Fig_pecanpy/DistanceCalculation/MR_Clusters_Dim50_WL80_nW10.s5.nodes_syntenic.csv')
########################################################

#########################################################
#####       Combine results and reduce GO terms     #####
#########################################################

##  PWYs
PWY_CommTarg <- PWY_CommTarg[,1:2]
#
PWY_CommFunt1 <- PWY_CommFunt[,c(1,2)]
colnames(PWY_CommFunt1)[2] <- "PWY"
PWY_CommFunt2 <- PWY_CommFunt[,c(1,3)]
colnames(PWY_CommFunt2)[2] <- "PWY"
PWY_CommFunt <- rbind(PWY_CommFunt1, PWY_CommFunt2)
#
PWY_Network <- unique(PWY_Network[,c("TF","PWY")])

PWY_CommTarg[,"Type"] <- "Comm.Target"
PWY_CommFunt[,"Type"] <- "Comm.Function"
PWY_Network[,"Type"] <- "Network.base"

###
PWYs <- rbind(PWY_CommTarg,
              PWY_CommFunt,
              PWY_Network) %>% unique()

##  GOs
GOs <- rbind(GO_CommFunt_Red[,c("TF","parent")] %>% dplyr::mutate(Type='Comm.Target'),
             GO_CommTarg_Red[,c("TF","parent")] %>% dplyr::mutate(Type='Comm.Function'),
             GO_Network_Red[,c("TF","parent")] %>% dplyr::mutate(Type='Network.base')) %>% 
  unique()  %>% dplyr::rename(GO.ID=parent)

# Add GOs from cluster with sinlge GO that could not be mapped
# to a parental GO
GOs <- rbind(GOs,
             GO_Network[!(GO_Network$TF %in% GO_Network_Red$TF),c('TF', 'GO.ID')] %>%
               unique() %>%
               dplyr::mutate(Type='Network.base'))
GOs <- unique(GOs)


## 1. Total TFs annotated by Method
#
GOs[,"Class"] <- paste0('GO.', GOs$Type)
PWYs[,"Class"] <- paste0('PWY.', PWYs$Type)

table(GOs$Class)
#
colnames(PWYs)[2] <- 'Annotation'
colnames(GOs)[2] <- 'Annotation'

## 1. TFs commonly annotated among methods
CombinedAnnotation <- rbind(GOs, PWYs)
head(CombinedAnnotation)

# write table as Sup. table S8
write.table(CombinedAnnotation[,c(1,2,4)], "Table_S8.txt",row.names = F, sep = '\t', quote = F)


## 2. PWYs/GOs by TF by Method
CombinedAnnotation[,c("TF","Class")]
multi_TFs <- as.data.table(table(CombinedAnnotation[,c("TF","Class")]))
multi_TFs[,"ClassType"] <- chop(multi_TFs$Class, "[.]", 1)
multi_TFs[,"TFname"] <- ReplaceName(multi_TFs$TF)
multi_TFs[,"ClassMethod"] <- paste0(chop(multi_TFs$Class, "[.]", 2),
                                    ". ",chop(multi_TFs$Class, "[.]", 3))
multi_TFs$ClassMethod <- gsub('Network. base', 'Network-based', multi_TFs$ClassMethod)
multi_TFs <- multi_TFs[multi_TFs$N > 0,]

sum(multi_TFs[multi_TFs$Class=='PWY.Network.base',]$N)
sum(multi_TFs[multi_TFs$Class=='GO.Network.base',]$N)

mean(multi_TFs[multi_TFs$Class=='PWY.Network.base',]$N)
mean(multi_TFs[multi_TFs$Class=='GO.Network.base',]$N)

multi_TFs %>% 
  dplyr::group_by(Class) %>%
  dplyr::summarise(mean=mean(N))

length(multi_TFs[multi_TFs$Class=='PWY.Network.base',]$N)
length(multi_TFs[multi_TFs$Class=='GO.Network.base',]$N)

# 
write.table(multi_TFs, 'Summary.Total.Annotation.txt', 
            col.names =T, row.names = F, quote = F, sep = '\t')

#########################################################


############################################################################
#####        PWY_Network and GO_Network_Red description Fig_S6         #####
############################################################################

PWY_Network
GO_Network_Red

# Combine GOs with reduce GO terms and unreduce 
Network_summary <- unique(GO_Network_Red[, c('TF', 'parent')]) %>% 
  rename("GO.ID"="parent") %>%
  rbind(GO_Network[!(GO_Network$TF %in% GO_Network_Red$TF), c('TF', 'GO.ID')]) %>%
  mutate(Annotation='GO') %>%
  rename("Ann.ID"="GO.ID")

fwrite(Network_summary, 
       'Total_NetworkBased_predictions.txt', row.names = F, col.names = T, sep = '\t')

# Combined GOs and PWY annotation into a sinlge DF
unique(PWY_Network[,c(1,2)]) %>%
  rename("Ann.ID"="PWY") %>%
  mutate(Annotation='PWY') %>%
  rbind(Network_summary) -> Network_summary
  
# Total TFs annotated 
list_NetSummary <- split(unique(Network_summary[,c(1,3)])$TF, unique(Network_summary[,c(1,3)])$Annotation)  
list_NetSummary <- lapply(list_NetSummary, unique)
#
list_NetSummary <- venn(list_NetSummary)
list_NetSummary <- as.list(attr(list_NetSummary, "intersections"))
list_NetSummary <- as.data.table(plyr::ldply(list_NetSummary, data.table)) 
list_NetSummary <- as.data.table(table(list_NetSummary$.id))
#
list_NetSummary <- plyr::ddply(list_NetSummary, "V1", transform, label_ypos=cumsum(N))
list_NetSummary$V1 <- factor(list_NetSummary$V1, levels = rev(list_NetSummary$V1))

# Total TFs annotated by Family
list_NetSummaryFamily <- left_join(Network_summary, TF_CoR[,c(1,3)], by=c('TF'="GeneID"))


## list
list_NetSummaryFamily <- unique(list_NetSummaryFamily[,c(1,3,4)])
list_NetSummaryFamily <- split(list_NetSummaryFamily[,c(1,2)]$TF, 
                               list_NetSummaryFamily[,c(1,2)]$Annotation)  

list_NetSummaryFamily <- lapply(list_NetSummaryFamily, unique)

## Make DFs with classes
list_NetSummaryFamily <- venn(list_NetSummaryFamily)
list_NetSummaryFamily <- as.list(attr(list_NetSummaryFamily, "intersections"))
list_NetSummaryFamily <- as.data.table(plyr::ldply(list_NetSummaryFamily, data.table)) 
list_NetSummaryFamily <- unique(list_NetSummaryFamily)

## Add families
list_NetSummaryFamily <- left_join(list_NetSummaryFamily, TF_CoR[,c(1,3)], by=c('V1'="GeneID"))

## Total TFs annotated per family by PWY and GO terms
FamilyAnnotate <- as.data.table(table(list_NetSummaryFamily[list_NetSummaryFamily$.id == 'GO:PWY',]$Family)) %>%
  left_join(as.data.table(table(TF_CoR$Family)), by='V1')
colnames(FamilyAnnotate) <- c('Family', 'Annotated', 'Total')

FamilyAnnotate[,'PerAnnotate'] <- (FamilyAnnotate$Annotated/FamilyAnnotate$Total)*100
mean(FamilyAnnotate$PerAnnotate)


Plot_FigS6c <- ggplot(list_NetSummary, aes(y='', x=N, fill=V1)) +
  geom_bar(stat="identity") +
  geom_text(aes(x=label_ypos + 5, label=scales::comma(N)), vjust=1,  color="black", size=3.5) +
  scale_fill_brewer(palette="Paired", name='') +
  theme_pubclean() +
  scale_y_discrete(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0), labels = comma) + 
  ylab("Annotation")+
  xlab("TFs") +
  theme(strip.text = element_text(size = 10), 
        legend.position = 'right',
        #axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1),
        #axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1),
        axis.text=element_text(size=10), 
        text = element_text(size=10, family="Times"))

Plot_FigS6c

Plot_FigS6d <- ggplot(FamilyAnnotate, aes(y=reorder(Family, Total), x=PerAnnotate)) +
  geom_bar(stat="identity", alpha=0.5) +
  geom_text(aes(x=PerAnnotate + 0.4, label=scales::comma(Annotated)), color="black", size=1.5) +
  #scale_fill_brewer(palette="Paired", name='') +
  theme_pubclean() +
  scale_y_discrete(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0), limits = c(0,101), labels = comma) + 
  ylab("TF Familiy")+
  xlab("Percentage of TF Familiy annotated") +
  theme(strip.text = element_text(size = 10), 
        legend.position = 'none',
        # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
        axis.text=element_text(size=10), 
        text = element_text(size=10, family="Times"))

Plot_S6cd <- Plot_FigS6c/Plot_FigS6d + plot_layout(heights = c(0.1, 1.5))
Plot_S6cd


pdf('../Fig_pecanpy/Plots/Plot_S6cd.pdf', width = 5, height = 9)
print(Plot_S6cd)
dev.off()

# Total PWYs/GOs by clusters
#list_NetSummary2 <- split(Network_summary[,c(1,3)]$TF, Network_summary[,c(1,3)]$Annotation)
# list_NetSummary2 <- as.data.table(table(Network_summary[,c(1,3)]))
# list_NetSummary2 <- subset(list_NetSummary2, N > 0)
# 
# mean(list_NetSummary2[list_NetSummary2$Annotation=="GO",]$N)
# mean(list_NetSummary2[list_NetSummary2$Annotation=="PWY",]$N)
# 
# 
# Plot_FigS6e <- ggplot(list_NetSummary2, aes(y=Annotation, x=N, fill=Annotation)) +
#   geom_boxplot(notch = T,  outlier.shape = NA) +
#   theme_pubclean() +
#   scale_x_continuous(expand = c(0,0), limits = c(0, 30), labels = comma) + 
#   ylab("")+
#   xlab("PWYs/GOs by TF") +
#   theme(strip.text = element_text(size = 10), 
#         # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
#         axis.text=element_text(size=10), 
#         legend.position = 'none',
#         text = element_text(size=10, family="Times")) 
#   
# Plot_FigS6de <- Plot_FigS6d/Plot_FigS6e
# Plot_FigS6de

############################################################################

#############################################################
#####          DEGs comparison: overlapping             #####
#############################################################

# Read DEGs from Peng et al 2022
DEGs <- readRDS("degs.rds")
unique(DEGs[,2:3])

DEGs <- subset(DEGs, reg.gid %in% DEGs_Peng) 

DEGs_names <- paste(DEGs$tf, DEGs$reg.gid, DEGs$tissue, sep = ':')
DEGs <- DEGs$ds
names(DEGs) <- DEGs_names

DEGs <- rbindlist(DEGs, idcol = T)

# Read DEGs from Erika
DEGs2 <- unique(fread('DEGs.part2.tsv')[,-c(3,6,10,11)])
DEGs2_names <- fread('DEGs.part2.tf.gid.tsv')
colnames(DEGs2_names)[2] <- 'TFid'

DEGs2 <- left_join(DEGs2, DEGs2_names[,c("tf","TFid")], by='tf')[,-c(1)]
DEGs2[,'.id'] <- paste0(toupper(chop(DEGs2$Genotype, '[_]',1)), '_', chop(DEGs2$Genotype, '[_]',2),
                        ":", DEGs2$TFid, ":", DEGs2$Tissue)

DEGs2 <- DEGs2[,c(".id", "gid","padj", "log2fc")]
DEGs


## Defined DEGs for test by overlapping of PWYs and GO terms
#DEGs_1 <- subset(DEGs, padj <= 0.05 & abs(log2fc) >=1)
DEGs_1 <- subset(DEGs, padj <= 0.05)
colnames(DEGs_1)[1:2] <- c('TF', "Target")
DEGs_1[,'TFid'] <- chop(DEGs_1$TF, '[:]',2)

DEGs_2 <- subset(DEGs2, padj <= 0.05)
colnames(DEGs_2)[1:2] <- c('TF', "Target")
DEGs_2[,'TFid'] <- chop(DEGs_2$TF, '[:]',2)

## 
DEGs_1 <- unique(rbind(DEGs_1, DEGs_2))
#j
DEGs_1 <- subset(DEGs_1, padj <= 0.05)
table(DEGs_1$TF)

unique(DEGs_1$TF)

head(DEGs_1)

DEGs2$gid

Tale_S9 <- DEGs_1

Tale_S9[,"Reference"] <- Tale_S9$TF %in% unique(DEGs2$.id)
Tale_S9$Reference[Tale_S9$Reference == TRUE] <- "Ellison et al., 2023"
Tale_S9$Reference[Tale_S9$Reference == FALSE] <- "Zhou et al., 2022"

write.table(Tale_S9, "Table_S9.txt", row.names = F, sep = '\t', quote = F)

#### DEGs PYW Enrichment 
#
DEGs_1_PWY <- Enrichmet_classes(DEGs_1)

DEGs_1_PWY <- subset(DEGs_1_PWY, Pval <= 0.05)
colnames(DEGs_1_PWY)[1] <- 'Index'

DEGs_1_PWY[,'TF'] <- chop(DEGs_1_PWY$Index, '[:]', 2)
DEGs_1_PWY[,'Treatment'] <- paste(chop(DEGs_1_PWY$Index, '[:]', 1),
                                  chop(DEGs_1_PWY$Index, '[:]', 3), sep = ":")
DEGs_1_PWY$PWY <- as.character(DEGs_1_PWY$PWY)
####

## Obs vs pred. PWYs 
DEGs_Peng <- c("Zm00001d033859", "Zm00001d018971", "Zm00001d020430", "Zm00001d033673", 
               "Zm00001d021191", "Zm00001d039694", "Zm00001d037317", "Zm00001d002654", 
               "Zm00001d028129")

# with out NKD1, and GT1
DEGs_Peng <- c("Zm00001d033859", "Zm00001d018971", "Zm00001d020430", "Zm00001d033673", 
               "Zm00001d021191", "Zm00001d037317")
# with Erika's TFs
DEGs_Peng <- unique(c(DEGs_Peng, DEGs2_names$TFid))

# TF names
TFdic <- as_tibble(read.table("../Data/Annotations/TF_Id_Name.txt", h=F, stringsAsFactors = F))


DEGs_Peng_in_PWY <- unique(DEGs_Peng[DEGs_Peng %in% PWYs$TF])
DEGs_Peng_in_PWY <- unique(DEGs_Peng[DEGs_Peng_in_PWY %in% unique(DEGs_1_PWY$TF)]) 
# ReplaceName(DEGs_Peng_in_PWY[17])

Pred_and_obs_PWYs <- lapply(DEGs_Peng_in_PWY, Get_PWY_methogd_DEGs)
names(Pred_and_obs_PWYs) <- DEGs_Peng_in_PWY

# Discard empty DFs
mask <- unlist(lapply(Pred_and_obs_PWYs, function(x) is.data.frame(x)))
Pred_and_obs_PWYs <- Pred_and_obs_PWYs[mask]

# FINAL DF
Pred_and_obs_PWYs <- rbindlist(Pred_and_obs_PWYs, idcol = T)

########
## DEG GO terms enrichment
########
DEGs_1$TF <- as.character(DEGs_1$TF)

DEGs_1_GOs <- lapply(unique(DEGs_1$TF), SuperGO_DEGs)
DEGs_1_GOs <- rbindlist(DEGs_1_GOs, idcol = F)
DEGs_1_GOs[,'TFid'] <- chop(DEGs_1_GOs$TF, '[:]', 2)

# Mask Na with significant numbers
mask <- DEGs_1_GOs$Significant >  DEGs_1_GOs$Expected & is.na(DEGs_1_GOs$classic)
DEGs_1_GOs$classic[mask] <- min(DEGs_1_GOs$classic, na.rm = T)

#temDEGs <- DEGs_1_GOs

DEGs_1_GOs <- DEGs_1_GOs %>%
  group_by(TF) %>%
  mutate(FDR=p.adjust(classic, method = 'fdr')) %>%
  filter(FDR <= 0.1)

idstestables <- unique(DEGs_1_GOs$TFid)[unique(DEGs_1_GOs$TFid) %in% unique(c(GO_CommTarg$TF, GO_CommFunt$TF, GO_Network$TF))]

unique(DEGs_1_GOs$TF[DEGs_1_GOs$TFid %in% idstestables])

DEGs_1_GOs

########

########
## Calculate GSS for TFs by DEG dataset
########

DEGs_1_GOs_TFs <- unique(DEGs_1_GOs$TF)

write.table(DEGs_1_GOs, 'DEGs_GO.db.txt', row.names = F, 
            sep = '\t', quote = T)


# GSS between predicted vs obs GO terms
Pred_and_obs_GOs <- lapply(DEGs_1_GOs_TFs, Get_GSS_method_DEGs)
mask <- unlist(lapply(Pred_and_obs_GOs, function(x) nrow(x)>0))
table(mask)

Pred_and_obs_GOs <- Pred_and_obs_GOs[mask]
Pred_and_obs_GOs <- rbindlist(Pred_and_obs_GOs)
Pred_and_obs_GOs[,'TFid'] <- chop(Pred_and_obs_GOs$TF, '[:]', 2)

Pred_and_obs_GOs[,'DEG'] <-  paste(chop(Pred_and_obs_GOs$TF, '[:]', 1),
                                   ":",
                                   chop(Pred_and_obs_GOs$TF, '[:]', 3), sep = '')

########

########
## Calculate GSS for DEG and random GOs
########

# Get 
set.seed(123)

GetGSS_DEGs_Random_CTarg <- function(tf){
  # make df to save output
  # DEGs
  go_obs <- subset(DEGs_1_GOs, TF == tf)$GO.ID
  
  # define tf id
  tfid <- chop(tf, '[:]', 2)
  
  # ComTarget
  n=length(subset(GO_CommTarg, TF == tfid)$GO.ID)
  
  print(paste0(ReplaceName(tfid), ' ', n))
  
  
  if( n > 1 ){
    # get real GOs
    ReadGOs <- subset(GO_CommTarg, TF == tfid)$GO.ID
    
    # get random GOs
    dbCTarg <- lapply(seq(1:100), 
                      function(x) sample(GO_CommTarg$GO.ID, size = n, replace = T))
    names(dbCTarg) <- paste0('R', seq(1:100))
    
    # GSS with random GOs
    RandomGSS <- lapply(dbCTarg, function(x)
                  mgoSim(go_obs, x, semData=Zm.GOSemSim.BP, 
                         measure="Wang", combine='BMA'))
    
    # GSS with reads set of GOs
    RandomGSS <- as_tibble(unlist(RandomGSS))
    colnames(RandomGSS) <- c("rGSS")
    
    # real obs value
    GSS <- mgoSim(go_obs, ReadGOs, semData=Zm.GOSemSim.BP, measure="Wang", combine='BMA')
    
    # add "realGSS"
    RandomGSS[,'RealGSS'] <- GSS
    
    RandomGSS[,'Class'] <- 'Comm.Target'
    RandomGSS[,'TF'] <- tf
    return(RandomGSS) 
  }
  
  
}

GetGSS_DEGs_Random_CFunct <- function(tf){
  # make df to save output
  # DEGs
  go_obs <- subset(DEGs_1_GOs, TF == tf)$GO.ID
  
  # define tf id
  tfid <- chop(tf, '[:]', 2)
  
  # ComTarget
  n=length(subset(GO_CommFunt, TF == tfid)$GO.ID)
  
  print(paste0(ReplaceName(tfid), ' ', n))
  
  
  if( n > 1 ){
    # get real GOs
    ReadGOs <- subset(GO_CommFunt, TF == tfid)$GO.ID
    
    # get random GOs
    dbCFunct <- lapply(seq(1:100), 
                      function(x) sample(GO_CommFunt$GO.ID, size = n, replace = T))
    names(dbCFunct) <- paste0('R', seq(1:100))
    
    # GSS with random GOs
    RandomGSS <- lapply(dbCFunct, function(x)
      mgoSim(go_obs, x, semData=Zm.GOSemSim.BP, 
             measure="Wang", combine='BMA'))
    
    # GSS with reads set of GOs
    RandomGSS <- as_tibble(unlist(RandomGSS))
    colnames(RandomGSS) <- c("rGSS")
    
    # real obs value
    GSS <- mgoSim(go_obs, ReadGOs, semData=Zm.GOSemSim.BP, measure="Wang", combine='BMA')
    
    # add "realGSS"
    RandomGSS[,'RealGSS'] <- GSS
    
    RandomGSS[,'Class'] <- 'Comm.Funct'
    RandomGSS[,'TF'] <- tf
    return(RandomGSS) 
  }
  
  
}

GetGSS_DEGs_Random_nbase <- function(tf){
  # make df to save output
  # DEGs
  go_obs <- subset(DEGs_1_GOs, TF == tf)$GO.ID
  
  # define tf id
  tfid <- chop(tf, '[:]', 2)
  
  # ComTarget
  n=length(subset(GO_Network, TF == tfid)$GO.ID)
  
  print(paste0(ReplaceName(tfid), ' ', n))
  
  
  if( n > 1 ){
    # get real GOs
    ReadGOs <- subset(GO_Network, TF == tfid)$GO.ID
    
    # get random GOs
    dbnbase <- lapply(seq(1:100), 
                       function(x) sample(GO_Network$GO.ID, size = n, replace = T))
    names(dbnbase) <- paste0('R', seq(1:100))
    
    # GSS with random GOs
    RandomGSS <- lapply(dbnbase, function(x)
      mgoSim(go_obs, x, semData=Zm.GOSemSim.BP, 
             measure="Wang", combine='BMA'))
    
    # GSS with reads set of GOs
    RandomGSS <- as_tibble(unlist(RandomGSS))
    colnames(RandomGSS) <- c("rGSS")
    
    # real obs value
    GSS <- mgoSim(go_obs, ReadGOs, semData=Zm.GOSemSim.BP, measure="Wang", combine='BMA')
    
    # add "realGSS"
    RandomGSS[,'RealGSS'] <- GSS
    
    RandomGSS[,'Class'] <- 'Comm.nBase'
    RandomGSS[,'TF'] <- tf
    return(RandomGSS) 
  }
  
  
}

RandomGSS_ctarg <- list()

RandomGSS_ctarg <- lapply(DEGs_1_GOs_TFs[1:10], function(x) GetGSS_DEGs_Random_CTarg(x))
RandomGSS_ctargv2 <- lapply(DEGs_1_GOs_TFs[11:20], function(x) GetGSS_DEGs_Random_CTarg(x))
RandomGSS_ctargv3 <- lapply(DEGs_1_GOs_TFs[21:length(DEGs_1_GOs_TFs)], function(x) GetGSS_DEGs_Random_CTarg(x))

mask <- unlist(lapply(RandomGSS_ctarg, function(x) is.data.frame(x)))
RandomGSS_ctarg <- rbindlist(RandomGSS_ctarg[mask])

mask <- unlist(lapply(RandomGSS_ctargv2, function(x) is.data.frame(x)))
RandomGSS_ctargv2 <- rbindlist(RandomGSS_ctargv2[mask])

mask <- unlist(lapply(RandomGSS_ctargv3, function(x) is.data.frame(x)))
RandomGSS_ctargv3 <- rbindlist(RandomGSS_ctargv3[mask])

RandomGSS_ctarg <- rbind(RandomGSS_ctarg, RandomGSS_ctargv2, RandomGSS_ctargv3)

#
RandomGSS_cFunct <- lapply(DEGs_1_GOs_TFs[1:2], function(x) GetGSS_DEGs_Random_CFunct(x))
mask <- unlist(lapply(RandomGSS_cFunct, function(x) is.data.frame(x)))
RandomGSS_cFunct <- rbindlist(RandomGSS_cFunct[mask])

RandomGSS_cFunctv2 <- lapply(DEGs_1_GOs_TFs[3:10], function(x) GetGSS_DEGs_Random_CFunct(x))
mask <- unlist(lapply(RandomGSS_cFunctv2, function(x) is.data.frame(x)))
RandomGSS_cFunctv2 <- rbindlist(RandomGSS_cFunctv2[mask])

RandomGSS_cFunctv3 <- lapply(DEGs_1_GOs_TFs[11:25], function(x) GetGSS_DEGs_Random_CFunct(x))
mask <- unlist(lapply(RandomGSS_cFunctv3, function(x) is.data.frame(x)))
RandomGSS_cFunctv3 <- RandomGSS_cFunctv3[mask]

RandomGSS_cFunctv4 <- lapply(DEGs_1_GOs_TFs[30:length(DEGs_1_GOs_TFs)], function(x) GetGSS_DEGs_Random_CFunct(x))
mask <- unlist(lapply(RandomGSS_cFunctv4, function(x) is.data.frame(x)))
RandomGSS_cFunctv4 <- RandomGSS_cFunctv4[mask]

RandomGSS_cFunctv <- rbind(RandomGSS_cFunctv, RandomGSS_cFunctv2,
                           RandomGSS_cFunctv3, RandomGSS_cFunctv4)


## Filter to comparable TFs
#Pred_and_obs_PWYs <- subset(Pred_and_obs_PWYs, .id %in% DEGs_Peng)
#Pred_and_obs_GOs  <- subset(Pred_and_obs_GOs, TFid %in% DEGs_Peng)

# Reduce version of GO term in DEGs
#DEGs_1_GOs_Red <- lapply(DEGs_1_GOs_TFs, function(x) ReduceGOs(DEGs_1_GOs, x)) 
#DEGs_1_GOs_Red <- rbindlist(DEGs_1_GOs_Red, idcol = F)


@@@
#############################################################

#############################################################
#####          DEGs comparison: fgsea                   #####
#############################################################

# Combine all DEGs
AllDEGs <- rbind(DEGs, DEGs2)

DEG.fgsea_PWY <- function(tf){
  
  degsvector <- subset(AllDEGs, .id == tf)
  id_degsvector <- degsvector$gid
  degsvector <- degsvector$log2fc
  names(degsvector) <- id_degsvector
  
  # define tf id
  tfid <- chop(tf, '[:]', 2)
  
  # PWYs in Comm. targets
  commtarg <- unique(subset(PWY_CommTarg, TF==tfid)$PWY)
  
  # PWYs in Comm. function
  commfunct <- unique(subset(PWY_CommFunt, TF==tfid)$PWY)
  
  # PWYs in network
  netbased <- unique(subset(PWY_Network, TF==tfid)$PWY)
  
  # list of target genes by PWY/GO
  commtarg <- CornCYC.list[commtarg]
  
  # list of target genes by PWY/GO
  commfunct <- CornCYC.list[commfunct]
  
  # list of target genes by PWY/GO
  netbased <- CornCYC.list[netbased]
  
  fgseaRes <- list()
  
  if(length(commtarg)>0){
    
    fgseaRes_1 <- fgseaMultilevel(pathways = commtarg, stats = degsvector, 
                                  minSize  = 2, maxSize  = 500, eps=0)
    
    fgseaRes_1 <- as_tibble(fgseaRes_1)
    fgseaRes_1 <- fgseaRes_1[order(fgseaRes_1$pval),]
    fgseaRes_1[,'TF'] <- tf
    fgseaRes_1[,'Class'] <- 'Comm.Target'
    
    fgseaRes[["Comm.Target"]] <- fgseaRes_1
    
  }
  if(length(commfunct)>0){
    
    fgseaRes_2 <- fgseaMultilevel(pathways = commfunct, stats = degsvector, 
                                  minSize  = 2, maxSize  = 500, eps=0)
    
    fgseaRes_2 <- as_tibble(fgseaRes_2)
    fgseaRes_2 <- fgseaRes_2[order(fgseaRes_2$pval),]
    fgseaRes_2[,'TF'] <- tf
    fgseaRes_2[,'Class'] <- 'Comm.Function'
    
    fgseaRes[["Comm.Function"]] <- fgseaRes_2
  }
  if(length(netbased)>0){
    
    fgseaRes_3 <- fgseaMultilevel(pathways = netbased, stats = degsvector, 
                                  minSize  = 2, maxSize  = 500, eps=0)
    
    fgseaRes_3 <- as_tibble(fgseaRes_3)
    fgseaRes_3 <- fgseaRes_3[order(fgseaRes_3$pval),]
    fgseaRes_3[,'TF'] <- tf
    fgseaRes_3[,'Class'] <- 'Network-based'
    
    fgseaRes[["Network.based"]] <- fgseaRes_3
  }
  
  fgseaRes <- as_tibble(rbindlist(fgseaRes, idcol = F))
  
  return(fgseaRes)
  
}

DEG.fgsea_GO <- function(tf){
  
  degsvector <- subset(AllDEGs, .id == tf)
  id_degsvector <- degsvector$gid
  degsvector <- degsvector$log2fc
  names(degsvector) <- id_degsvector
  
  # define tf id
  tfid <- chop(tf, '[:]', 2)
  
  # GOs in Comm. targets
  commtarg <- unique(subset(GO_CommTarg_Red, TF==tfid)$parent)
  
  # GOs in Comm. function
  commfunct <- unique(subset(GO_CommFunt_Red, TF==tfid)$parent)
  
  # GOs in network
  netbased <- unique(subset(GO_Network_Red, TF==tfid)$parent)
  
  # list of target genes by PWY/GO
  commtarg <- GO.list[commtarg]
  
  # list of target genes by PWY/GO
  commfunct <- GO.list[commfunct]
  
  # list of target genes by PWY/GO
  netbased <- GO.list[netbased]
  
  # GSEA by GO and by PCC network
  fgseaRes <- list()
  
  if(length(commtarg)>0){
    
    fgseaRes_1 <- fgseaMultilevel(pathways = commtarg, stats = degsvector, 
                                  minSize  = 5, maxSize  = 1000, eps=0)
    
    fgseaRes_1 <- as_tibble(fgseaRes_1)
    fgseaRes_1 <- fgseaRes_1[order(fgseaRes_1$pval),]
    fgseaRes_1[,'TF'] <- tf
    fgseaRes_1[,'Class'] <- 'Comm.Target'
    
    fgseaRes[["Comm.Target"]] <- fgseaRes_1
    
  }
  if(length(commfunct)>0){
    
    fgseaRes_2 <- fgseaMultilevel(pathways = commfunct, stats = degsvector, 
                                  minSize  = 5, maxSize  = 1000, eps=0)
    
    fgseaRes_2 <- as_tibble(fgseaRes_2)
    fgseaRes_2 <- fgseaRes_2[order(fgseaRes_2$pval),]
    fgseaRes_2[,'TF'] <- tf
    fgseaRes_2[,'Class'] <- 'Comm.Function'
    
    fgseaRes[["Comm.Function"]] <- fgseaRes_2
  }
  if(length(netbased)>0){
    
    fgseaRes_3 <- fgseaMultilevel(pathways = netbased, stats = degsvector, 
                                  minSize  = 5, maxSize  = 1000, eps=0)
    
    fgseaRes_3 <- as_tibble(fgseaRes_3)
    fgseaRes_3 <- fgseaRes_3[order(fgseaRes_3$pval),]
    fgseaRes_3[,'TF'] <- tf
    fgseaRes_3[,'Class'] <- 'Network-based'
    
    fgseaRes[["Network.based"]] <- fgseaRes_3
  }
  
  fgseaRes <- as_tibble(rbindlist(fgseaRes, idcol = F))
  
  return(fgseaRes)  
  
}

# PWYs
GSEA_PWYs <- lapply(DEGs_1_GOs_TFs, DEG.fgsea_PWY)
mask <- unlist(lapply(GSEA_PWYs, function(x) nrow(x) >0 ))
GSEA_PWYs <- GSEA_PWYs[mask]
GSEA_PWYs <- rbindlist(GSEA_PWYs)

# PWYs P-val class
GSEA_PWYs[,'PvalClass'] <- GSEA_PWYs$pval <= 0.05

GSEA_PWYs_Freq <- as.data.table(table(GSEA_PWYs[,c("TF", "Class", "PvalClass")]))
GSEA_PWYs_Freq <- GSEA_PWYs_Freq[GSEA_PWYs_Freq$N > 0,]

GSEA_PWYs_Freq %>%
  group_by(TF, Class) %>% 
  mutate(Total=sum(N), Fraction=N/sum(N)) %>%
  arrange(TF) -> GSEA_PWYs_Freq

GSEA_PWYs_Freq[,'DEG'] <- paste0(chop(GSEA_PWYs_Freq$TF, '[:]',1),
                                 ":",
                                 chop(GSEA_PWYs_Freq$TF, '[:]',3))

# GOs
GSEA_GOs <- lapply(DEGs_1_GOs_TFs, DEG.fgsea_GO)
mask <- unlist(lapply(GSEA_GOs, function(x) nrow(x)>0))
GSEA_GOs <- GSEA_GOs[mask]
GSEA_GOs <- rbindlist(GSEA_GOs)

# PWYs P-val class
GSEA_GOs[,'PvalClass'] <- GSEA_GOs$pval <= 0.05
GSEA_GOs_Freq <- as.data.table(table(GSEA_GOs[,c("TF", "Class", "PvalClass")]))
GSEA_GOs_Freq <- GSEA_GOs_Freq[GSEA_GOs_Freq$N > 0,]

GSEA_GOs_Freq %>%
  group_by(TF, Class) %>% 
  mutate(Total=sum(N), Fraction=N/sum(N)) %>%
  arrange(TF) -> GSEA_GOs_Freq

GSEA_GOs_Freq[,'DEG'] <- paste0(chop(GSEA_GOs_Freq$TF, '[:]',1),
                                 ":",
                                 chop(GSEA_GOs_Freq$TF, '[:]',3))

#############################################################


##############################################
#####              Plots                 #####
##############################################

#####
# 1. Total TFs-function interactions
#####
CombinedAnnotation %>% 
  dplyr::group_by(Class) %>%
  dplyr::select(Class) %>%
  dplyr::summarise(n=n()) -> DF1

DF1[,"ClassType"] <- chop(DF1$Class, "[.]", 1)
DF1[,"ClassMethod"] <- paste0(chop(DF1$Class, "[.]", 2),". ",
                              chop(DF1$Class, "[.]", 3))
DF1$ClassMethod <- gsub('Network. base', 'Network-based', DF1$ClassMethod)

DF1$ClassMethod <- factor(DF1$ClassMethod, levels = rev(c("Comm. Target", 
                                                          "Comm. Function",
                                                          "Network-based")))
DF1$ClassType <- factor(DF1$ClassType, levels = c('PWY', 'GO'))

ggplot(DF1, aes(x=n, y=ClassMethod, fill=ClassMethod)) +
  geom_bar(stat="identity") +
  geom_text(aes(x=n-(n*0.4), label=scales::comma(n)), vjust=0.5,  color="black", size=3.5) +
  theme_pubclean() +
  scale_fill_viridis(discrete = T, option = "D", alpha = 0.5, direction = 1) + 
  facet_grid(ClassType ~ . )  + # scales = 'free_x'
  scale_x_log10(expand=c(0,0), lim=c(1, 30000), labels = comma) +
  annotation_logticks(sides = "b", color = 'black') +
  ylab("")+
  xlab("TFs-function ") +
  theme(strip.text = element_text(size = 10), 
        # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
        axis.text=element_text(size=10), 
        legend.position = 'none',
        text = element_text(size=10, family="Times")) -> Plot_1

Plot_1
#####

#####
## 2. GOs/PWYs by TF and by Method
#####

# levels
multi_TFs$ClassMethod <- factor(multi_TFs$ClassMethod, levels = rev(c("Comm. Target", 
                                                                      "Comm. Function",
                                                                      "Network-based")))
multi_TFs$ClassType <- factor(multi_TFs$ClassType, levels = c('PWY', 'GO'))


ggplot(multi_TFs, aes(x=N, y=ClassMethod, fill=ClassMethod))+
  geom_boxplot(notch = T, outlier.size = 0.5)+
  theme_pubclean() +
  scale_fill_viridis(discrete = T, option = "D", alpha = 0.5, direction = 1) + 
  scale_x_continuous(expand=c(0,0), limits = c(0,230)) + 
  facet_grid(ClassType ~ . , scales = 'free_x')  +
  ylab("")+
  xlab("PWYs/GOs by TF") +
  theme(strip.text = element_text(size = 10), 
        # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
        axis.text=element_text(size=10), 
        legend.position = 'none',
        text = element_text(size=10, family="Times")) + 
  stat_compare_means(label = "p.signif", method = "t.test",
                     comparisons=list(c("Comm. Target", "Comm. Function"),
                                      c('Comm. Function', 'Network-based'),
                                      c('Comm. Target', 'Network-based'))) -> Plot_2
  
Plot_2
pdf("Plots/TotalAnnotation_plot.pdf", width=3.5, height=2.5)
print(Plot_2)
dev.off()

Plot_2ab <- Plot_1/Plot_2

pdf("Plots/Plot_2ab.pdf", width=3.5, height=3)
print(Plot_2ab)
dev.off()


multi_TFs %>%
  group_by(Class, ClassType) %>%
  dplyr::summarise(M=mean(N))

#####

#####
# 3. TFs commonly annotated among methods
#####
upsetInput <- as.data.table(table(unique(CombinedAnnotation[,c('TF',"Class")])))
upsetInput <- dcast(upsetInput,  formula = TF~Class, value.var = "N")

Plotupset <- upset(upsetInput, nsets = 6, number.angles = 0, point.size = 2, line.size = 1, 
      order.by = "freq", scale.intersections="identity",
      mb.ratio = c(0.6, 0.4),
      mainbar.y.label = "Annotation\nintersections", 
      sets.x.label = "Total TF Annotated", 
      text.scale = c(1.3, 1.3, 1, 1, 1, 0.7))

pdf("Plots/upset_plot.pdf", width=6, height=3)
print(Plotupset)
dev.off()

#####

#####
# 4. TFs with prediction in all three methods
#####
TFdic <- as_tibble(read.table("../Data/Annotations/TF_Id_Name.txt", h=F, stringsAsFactors = F))

# mask
dim(upsetInput)
To_PWY <- rowSums(upsetInput[,c("PWY.Comm.Function", "PWY.Comm.Target", "PWY.Network.base")]) == 3
To_GO <- rowSums(upsetInput[,c("GO.Comm.Function", "GO.Comm.Target", "GO.Network.base")]) == 3


To_PWY_GO <- rowSums(upsetInput[,c("PWY.Comm.Function", "PWY.Comm.Target", "PWY.Network.base",
                                   "GO.Comm.Function", "GO.Comm.Target", "GO.Network.base")]) == 6


To_PWY <- upsetInput$TF[To_PWY]
To_GO <- upsetInput$TF[To_GO]
To_PWY_GO <- upsetInput$TF[To_PWY_GO]

To_PWY <- To_PWY[!(To_PWY %in% To_PWY_GO)]
To_GO <- To_GO[!(To_GO %in% To_PWY_GO)]

length(To_PWY)
length(To_GO)
length(To_PWY_GO)

View(as.data.frame(ReplaceName(To_PWY)))
ReplaceName(To_GO)
ReplaceName(To_PWY_GO)

write.table(data.frame(TF=c(To_GO, To_PWY_GO)), 
            'TFs_2_test_RNets.txt', sep = '\t', quote = F, 
            row.names = F, 
            col.names = F)


#####

#####
# 5. TFs with DEGs in Peng data: Possible comparisons 
#####
# KN1	Zm00001d033859
# O2	Zm00001d018971
# RA1	Zm00001d020430
# TB1	Zm00001d033673
# bZIP22	Zm00001d021191
# RA2	Zm00001d039694
# FEA4	Zm00001d037317
# NKD1	Zm00001d002654
# GT1	Zm00001d028129



ReplaceName(DEGs_Peng)

# Heatmap description
DF2 <- upsetInput[upsetInput$TF %in% unique(c(DEGs_1$TFid, DEGs_2$TFid)),] %>% 
  gather(key, value, -TF) %>% as_tibble()

DF2[,"TFname"] <- ReplaceName(DF2$TF)

DF2$key <- factor(DF2$key, levels = c("PWY.Comm.Target", "PWY.Comm.Function", "PWY.Network.base",
                                      "GO.Comm.Target", "GO.Comm.Function", "GO.Network.base"))
DF2$value <- DF2$value > 0

ggplot(DF2, aes(x=key, y=TFname, fill=value)) +
  geom_tile() +
  scale_fill_manual(values = c('white',"#9400D3")) +
  theme_bw() + 
  ylab('TF') + xlab("Annotation") +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text = element_text(size = 10), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.title = element_text(size = 10, face = 'bold'), 
        axis.title = element_text(face = 'bold'),
        legend.position = 'none',
        text = element_text(size=10, family="Times")) -> Plot_3
Plot_3

#pdf("Plots/Plot_2d.pdf", width=2.5, height=3)
#print(Plot_3)
#dev.off()
#####

#####
# 6. TFs with DEGs: PWYs overlapping and GSEA
#####

# Overlaping of PWYs
Pred_and_obs_PWYs[,'DEGname'] <- paste0(Pred_and_obs_PWYs$DEG,' [',
                                        Pred_and_obs_PWYs$n.DEG,
                                        ']')

Pred_and_obs_PWYs$Method <- factor(Pred_and_obs_PWYs$Method, 
                               levels = c("Comm.Target", "Comm.Function", "Network.base"))

Pred_and_obs_PWYs$DEGname <- gsub('FEA4:ear', 'FEA4_ear', Pred_and_obs_PWYs$DEGname)
Pred_and_obs_PWYs$DEGname <- gsub('O2:endosperm', 'O2', Pred_and_obs_PWYs$DEGname)
Pred_and_obs_PWYs$DEGname <- gsub('TB1:tiller_buds_12DAP', 'TB1_buds_12DAP', Pred_and_obs_PWYs$DEGname)
Pred_and_obs_PWYs$DEGname <- gsub('TB1:tiller_buds_8DAP', 'TB1_buds_8DAP', Pred_and_obs_PWYs$DEGname)
Pred_and_obs_PWYs$DEGname <- gsub('TB1:tiller_buds_8DAP', 'TB1_buds_8DAP', Pred_and_obs_PWYs$DEGname)
Pred_and_obs_PWYs$DEGname <- gsub('KN1:ear', 'KN1_ear', Pred_and_obs_PWYs$DEGname)
Pred_and_obs_PWYs$DEGname <- gsub('KN1:leaf', 'KN1_leaf', Pred_and_obs_PWYs$DEGname)
Pred_and_obs_PWYs$DEGname <- gsub('KN1:SAM', 'KN1_SAM', Pred_and_obs_PWYs$DEGname)
Pred_and_obs_PWYs$DEGname <- gsub('KN1:tassel', 'KN1_tassel', Pred_and_obs_PWYs$DEGname)
Pred_and_obs_PWYs$DEGname <- gsub('RA1:ear_2mm', 'RA1', Pred_and_obs_PWYs$DEGname)
Pred_and_obs_PWYs$DEGname <- gsub("BAF6021_m1:tassel_stem", 'BAF6021m1_tassel', Pred_and_obs_PWYs$DEGname)
Pred_and_obs_PWYs$DEGname <- gsub("BAF6021_m2:tassel_stem", 'BAF6021m2_tassel', Pred_and_obs_PWYs$DEGname)
Pred_and_obs_PWYs$DEGname <- gsub("bZIP22:kernel", 'bZIP22_kernel', Pred_and_obs_PWYs$DEGname)
Pred_and_obs_PWYs$DEGname <- gsub("BZIP76_m2:leaf", 'BZIP76m2_leaf', Pred_and_obs_PWYs$DEGname)
Pred_and_obs_PWYs$DEGname <- gsub("BZIP76_m3:leaf", 'BZIP76m3_leaf', Pred_and_obs_PWYs$DEGname)
Pred_and_obs_PWYs$DEGname <- gsub("C3H42_m1:tassel_stem", 'C3H42m1_tassel_stem', Pred_and_obs_PWYs$DEGname)
Pred_and_obs_PWYs$DEGname <- gsub("E2F13_m1:coleoptile_tip", 'E2F13m1_coleoptile', Pred_and_obs_PWYs$DEGname)
Pred_and_obs_PWYs$DEGname <- gsub("E2F19_m1:leaf", 'E2F19m1_leaf', Pred_and_obs_PWYs$DEGname)
Pred_and_obs_PWYs$DEGname <- gsub("E2F19_m2:leaf", 'E2F19m2_leaf', Pred_and_obs_PWYs$DEGname)
Pred_and_obs_PWYs$DEGname <- gsub("GRAS52_m1:embryo_imbibed", 'GRAS52m1_embryo', Pred_and_obs_PWYs$DEGname)
Pred_and_obs_PWYs$DEGname <- gsub("GRAS75_m1:embryo_imbibed", 'GRAS75m1_embryo', Pred_and_obs_PWYs$DEGname)
Pred_and_obs_PWYs$DEGname <- gsub("HSF13_m1m2:leaf", 'HSF13m1m2_leaf', Pred_and_obs_PWYs$DEGname)
Pred_and_obs_PWYs$DEGname <- gsub("HSF18_m1:embryo_imbibed", 'HSF18m1_embryo', Pred_and_obs_PWYs$DEGname)
Pred_and_obs_PWYs$DEGname <- gsub("HSF20_m1:embryo_imbibed", 'HSF20m1_embryo', Pred_and_obs_PWYs$DEGname)
Pred_and_obs_PWYs$DEGname <- gsub("HSF24_m3:tassel", 'HSF24m3_tassel', Pred_and_obs_PWYs$DEGname)
Pred_and_obs_PWYs$DEGname <- gsub("HSF24_m4:tassel", 'HSF24m4_tassel', Pred_and_obs_PWYs$DEGname)
Pred_and_obs_PWYs$DEGname <- gsub("HSF29_m1:embryo_imbibed", 'HSF29m1_embryo', Pred_and_obs_PWYs$DEGname)
Pred_and_obs_PWYs$DEGname <- gsub("HSF29_m2:embryo_imbibed", 'HSF29m2_embryo', Pred_and_obs_PWYs$DEGname)
Pred_and_obs_PWYs$DEGname <- gsub("HSF6_m1:embryo_imbibed", 'HSF6m1_embryo', Pred_and_obs_PWYs$DEGname)
Pred_and_obs_PWYs$DEGname <- gsub("HSF6_m2:embryo_imbibed", 'HSF6m2_embryo', Pred_and_obs_PWYs$DEGname)
Pred_and_obs_PWYs$DEGname <- gsub("JMJ13_m4:tassel_stem", 'JMJ13m4_tassel', Pred_and_obs_PWYs$DEGname)
Pred_and_obs_PWYs$DEGname <- gsub("MYBR21_m1:embryo_imbibed", 'MYBR21m1_embryo', Pred_and_obs_PWYs$DEGname)
Pred_and_obs_PWYs$DEGname <- gsub("MYBR32_m1:leaf", 'MYBR32m1_leaf', Pred_and_obs_PWYs$DEGname)
Pred_and_obs_PWYs$DEGname <- gsub("ORPHAN249_m2:embryo_imbibed", 'ORPHAN249m2_embryo', Pred_and_obs_PWYs$DEGname)
Pred_and_obs_PWYs$DEGname <- gsub("SBP20_m2:embryo_imbibed", 'SBP20m2_embryo', Pred_and_obs_PWYs$DEGname)
Pred_and_obs_PWYs$DEGname <- gsub("SBP20_m3:embryo_imbibed", 'SBP20m3_embryo', Pred_and_obs_PWYs$DEGname)
Pred_and_obs_PWYs$DEGname <- gsub("WRKY2_m2:coleoptile_tip", 'WRKY2m2_coleoptile', Pred_and_obs_PWYs$DEGname)
Pred_and_obs_PWYs$DEGname <- gsub("WRKY8_m1:embryo_imbibed", 'WRKY8m1_embryo', Pred_and_obs_PWYs$DEGname)
Pred_and_obs_PWYs$DEGname <- gsub("WRKY8_m2:embryo_imbibed", 'WRKY8m2_embryo', Pred_and_obs_PWYs$DEGname)
Pred_and_obs_PWYs$DEGname <- gsub("WRKY82_m1:embryo_imbibed", 'WRKY82m1_embryo', Pred_and_obs_PWYs$DEGname)
Pred_and_obs_PWYs$DEGname <- gsub("WRKY87_m1:embryo_imbibed", 'WRKY87m1_embryo', Pred_and_obs_PWYs$DEGname)
Pred_and_obs_PWYs$DEGname <- gsub("WRKY87_m2:embryo_imbibed", 'WRKY87m2_embryo', Pred_and_obs_PWYs$DEGname)

ggplot(Pred_and_obs_PWYs, 
       aes(x=Method, y=DEGname, fill=Pval <= 0.05, label=n.Common)) +
  geom_tile() +
  geom_text(size=2.5) +
  scale_fill_manual(values = c('#E0E0E0',"#FF66FF")) +
  theme_bw() + 
  ylab('TF [PWYs in DEGs]') + xlab("") +
  scale_y_discrete(limits=rev) + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text = element_text(size = 10), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.title = element_text(size = 10, face = 'bold'), 
        axis.title = element_text(face = 'bold'),
        legend.position = 'none',       
        text = element_text(size=10, family="Times")) -> Plot_4

#  GSEA plot
# keep TFs with data in at least two methods
# GSEA_PWYs_Freq <- GSEA_PWYs_Freq[!(GSEA_PWYs_Freq$DEG %in% c("GT1:tiller_buds_12DAP",
#                                                              "GT1:tiller_buds_8DAP", 
#                                                              "NKD1:aleurone",
#                                                              "NKD1:endosperm")),]

unique(GSEA_PWYs_Freq$DEG)
GSEA_PWYs_Freq$DEG <- gsub('FEA4:ear', 'FEA4_ear', GSEA_PWYs_Freq$DEG)
GSEA_PWYs_Freq$DEG <- gsub('O2:endosperm', 'O2', GSEA_PWYs_Freq$DEG)
GSEA_PWYs_Freq$DEG <- gsub('TB1:tiller_buds_12DAP', 'TB1_buds_12DAP', GSEA_PWYs_Freq$DEG)
GSEA_PWYs_Freq$DEG <- gsub('TB1:tiller_buds_8DAP', 'TB1_buds_8DAP', GSEA_PWYs_Freq$DEG)
GSEA_PWYs_Freq$DEG <- gsub('TB1:tiller_buds_8DAP', 'TB1_buds_8DAP', GSEA_PWYs_Freq$DEG)
GSEA_PWYs_Freq$DEG <- gsub('KN1:ear', 'KN1_ear', GSEA_PWYs_Freq$DEG)
GSEA_PWYs_Freq$DEG <- gsub('KN1:leaf', 'KN1_leaf', GSEA_PWYs_Freq$DEG)
GSEA_PWYs_Freq$DEG <- gsub('KN1:SAM', 'KN1_SAM', GSEA_PWYs_Freq$DEG)
GSEA_PWYs_Freq$DEG <- gsub('KN1:tassel', 'KN1_tassel', GSEA_PWYs_Freq$DEG)
GSEA_PWYs_Freq$DEG <- gsub('RA1:ear_2mm', 'RA1', GSEA_PWYs_Freq$DEG)
GSEA_PWYs_Freq$DEG <- gsub("BAF6021_m1:tassel_stem", 'BAF6021m1_tassel', GSEA_PWYs_Freq$DEG)
GSEA_PWYs_Freq$DEG <- gsub("BAF6021_m2:tassel_stem", 'BAF6021m2_tassel', GSEA_PWYs_Freq$DEG)
GSEA_PWYs_Freq$DEG <- gsub("bZIP22:kernel", 'bZIP22_kernel', GSEA_PWYs_Freq$DEG)
GSEA_PWYs_Freq$DEG <- gsub("BZIP76_m2:leaf", 'BZIP76m2_leaf', GSEA_PWYs_Freq$DEG)
GSEA_PWYs_Freq$DEG <- gsub("BZIP76_m3:leaf", 'BZIP76m3_leaf', GSEA_PWYs_Freq$DEG)
GSEA_PWYs_Freq$DEG <- gsub("C3H42_m1:tassel_stem", 'C3H42m1_tassel_stem', GSEA_PWYs_Freq$DEG)
GSEA_PWYs_Freq$DEG <- gsub("E2F13_m1:coleoptile_tip", 'E2F13m1_coleoptile', GSEA_PWYs_Freq$DEG)
GSEA_PWYs_Freq$DEG <- gsub("E2F19_m1:leaf", 'E2F19m1_leaf', GSEA_PWYs_Freq$DEG)
GSEA_PWYs_Freq$DEG <- gsub("E2F19_m2:leaf", 'E2F19m2_leaf', GSEA_PWYs_Freq$DEG)
GSEA_PWYs_Freq$DEG <- gsub("GRAS52_m1:embryo_imbibed", 'GRAS52m1_embryo', GSEA_PWYs_Freq$DEG)
GSEA_PWYs_Freq$DEG <- gsub("GRAS75_m1:embryo_imbibed", 'GRAS75m1_embryo', GSEA_PWYs_Freq$DEG)
GSEA_PWYs_Freq$DEG <- gsub("HSF13_m1m2:leaf", 'HSF13m1m2_leaf', GSEA_PWYs_Freq$DEG)
GSEA_PWYs_Freq$DEG <- gsub("HSF18_m1:embryo_imbibed", 'HSF18m1_embryo', GSEA_PWYs_Freq$DEG)
GSEA_PWYs_Freq$DEG <- gsub("HSF20_m1:embryo_imbibed", 'HSF20m1_embryo', GSEA_PWYs_Freq$DEG)
GSEA_PWYs_Freq$DEG <- gsub("HSF24_m3:tassel", 'HSF24m3_tassel', GSEA_PWYs_Freq$DEG)
GSEA_PWYs_Freq$DEG <- gsub("HSF24_m4:tassel", 'HSF24m4_tassel', GSEA_PWYs_Freq$DEG)
GSEA_PWYs_Freq$DEG <- gsub("HSF29_m1:embryo_imbibed", 'HSF29m1_embryo', GSEA_PWYs_Freq$DEG)
GSEA_PWYs_Freq$DEG <- gsub("HSF29_m2:embryo_imbibed", 'HSF29m2_embryo', GSEA_PWYs_Freq$DEG)
GSEA_PWYs_Freq$DEG <- gsub("HSF6_m1:embryo_imbibed", 'HSF6m1_embryo', GSEA_PWYs_Freq$DEG)
GSEA_PWYs_Freq$DEG <- gsub("HSF6_m2:embryo_imbibed", 'HSF6m2_embryo', GSEA_PWYs_Freq$DEG)
GSEA_PWYs_Freq$DEG <- gsub("JMJ13_m4:tassel_stem", 'JMJ13m4_tassel', GSEA_PWYs_Freq$DEG)
GSEA_PWYs_Freq$DEG <- gsub("MYBR21_m1:embryo_imbibed", 'MYBR21m1_embryo', GSEA_PWYs_Freq$DEG)
GSEA_PWYs_Freq$DEG <- gsub("MYBR32_m1:leaf", 'MYBR32m1_leaf', GSEA_PWYs_Freq$DEG)
GSEA_PWYs_Freq$DEG <- gsub("ORPHAN249_m2:embryo_imbibed", 'ORPHAN249m2_embryo', GSEA_PWYs_Freq$DEG)
GSEA_PWYs_Freq$DEG <- gsub("SBP20_m2:embryo_imbibed", 'SBP20m2_embryo', GSEA_PWYs_Freq$DEG)
GSEA_PWYs_Freq$DEG <- gsub("SBP20_m3:embryo_imbibed", 'SBP20m3_embryo', GSEA_PWYs_Freq$DEG)
GSEA_PWYs_Freq$DEG <- gsub("WRKY2_m2:coleoptile_tip", 'WRKY2m2_coleoptile', GSEA_PWYs_Freq$DEG)
GSEA_PWYs_Freq$DEG <- gsub("WRKY8_m1:embryo_imbibed", 'WRKY8m1_embryo', GSEA_PWYs_Freq$DEG)
GSEA_PWYs_Freq$DEG <- gsub("WRKY8_m2:embryo_imbibed", 'WRKY8m2_embryo', GSEA_PWYs_Freq$DEG)
GSEA_PWYs_Freq$DEG <- gsub("WRKY82_m1:embryo_imbibed", 'WRKY82m1_embryo', GSEA_PWYs_Freq$DEG)
GSEA_PWYs_Freq$DEG <- gsub("WRKY87_m1:embryo_imbibed", 'WRKY87m1_embryo', GSEA_PWYs_Freq$DEG)
GSEA_PWYs_Freq$DEG <- gsub("WRKY87_m2:embryo_imbibed", 'WRKY87m2_embryo', GSEA_PWYs_Freq$DEG)


# method levels
GSEA_PWYs_Freq$Class <- factor(GSEA_PWYs_Freq$Class, 
                               levels = c("Comm.Target", "Comm.Function", "Network-based"))


ggplot(GSEA_PWYs_Freq, aes(y=Class, x=Fraction, fill=PvalClass)) +
  geom_bar(stat="identity") +
  #geom_text(aes(x=label_ypos + 2, label=scales::comma(N)), vjust=1,  color="black", size=3.5) +
  scale_fill_brewer(palette="Paired", name='P-Value', 
                    label=c(expression("">=0.05), expression(""<=0.05))) +
  #bquote(Log[10] ~ "TMM")
  theme_pubclean() +
  #scale_y_discrete(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0)) + 
  ylab("Method") +
  xlab("PWYs Fraction") + 
  facet_grid(DEG ~., scales = "free_y") +
  theme(strip.text.y = element_text(size = 5, angle = 0), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text=element_text(size=10), 
        legend.position = 'bottom',
        text = element_text(size=10, family="Times")) -> Plot_6

Plot_S7 <- Plot_4 + Plot_6
pdf("Plots/Plot_S7.pdf", width=8, height=10)
print(Plot_S7)
dev.off()
# PWYs
# keep TFs with data in at least two methods
# GSEA_PWYs_Freq <- GSEA_PWYs_Freq[!(GSEA_PWYs_Freq$DEG %in% c("GT1:tiller_buds_12DAP",
#                                                              "GT1:tiller_buds_8DAP", 
#                                                              "NKD1:aleurone",
#                                                              "NKD1:endosperm")),]


#####

#####
# 7 TFs with DEGs: GOs overlapping and GSEA
#####
## GSS with GO from DEGs

# GSS GOs 
unique(Pred_and_obs_GOs$DEG)
Pred_and_obs_GOs$Class <- factor(Pred_and_obs_GOs$Class,
                                 levels = c("Comm.Target", "Comm.Funct", "Network-based"))
unique(Pred_and_obs_GOs$DEG)

Pred_and_obs_GOs %>%
  group_by(DEG, GO_obs) %>%
  arrange(-GSS, .by_group = TRUE) %>%
  top_n(20) -> Pred_and_obs_GOs_max

unique(Pred_and_obs_GOs_max$DEG)
Pred_and_obs_GOs_max$DEG <- gsub('RA1:ear_2mm', 'RA1_ear', Pred_and_obs_GOs_max$DEG)
Pred_and_obs_GOs_max$DEG <- gsub('FEA4:ear', 'FEA4_ear', Pred_and_obs_GOs_max$DEG)
Pred_and_obs_GOs_max$DEG <- gsub('O2:endosperm', 'O2', Pred_and_obs_GOs_max$DEG)
Pred_and_obs_GOs_max$DEG <- gsub('TB1:tiller_buds_12DAP', 'TB1\nbuds_12DAP', Pred_and_obs_GOs_max$DEG)
Pred_and_obs_GOs_max$DEG <- gsub('TB1:tiller_buds_8DAP', 'TB1\nbuds_8DAP', Pred_and_obs_GOs_max$DEG)
Pred_and_obs_GOs_max$DEG <- gsub('KN1:ear', 'KN1\near', Pred_and_obs_GOs_max$DEG)
Pred_and_obs_GOs_max$DEG <- gsub('KN1:leaf', 'KN1\nleaf', Pred_and_obs_GOs_max$DEG)
Pred_and_obs_GOs_max$DEG <- gsub('KN1:SAM', 'KN1\nSAM', Pred_and_obs_GOs_max$DEG)
Pred_and_obs_GOs_max$DEG <- gsub('KN1:tassel', 'KN1\ntassel', Pred_and_obs_GOs_max$DEG)
Pred_and_obs_GOs_max$DEG <- gsub('HSF13_m1m2:leaf', 'HSF13m1m2\nleaf', Pred_and_obs_GOs_max$DEG)
Pred_and_obs_GOs_max$DEG <- gsub('HSF18_m1:embryo_imbibed', 'HSF18m1\nembryo', Pred_and_obs_GOs_max$DEG)
Pred_and_obs_GOs_max$DEG <- gsub('HSF20_m1:embryo_imbibed', 'HSF20m1\nembryo', Pred_and_obs_GOs_max$DEG)
Pred_and_obs_GOs_max$DEG <- gsub('HSF24_m3:tassel', 'HSF24_m3\ntassel', Pred_and_obs_GOs_max$DEG)
Pred_and_obs_GOs_max$DEG <- gsub('HSF24_m4:tassel', 'HSF24_m4\ntassel', Pred_and_obs_GOs_max$DEG)
Pred_and_obs_GOs_max$DEG <- gsub('HSF29_m1:embryo_imbibed', 'HSF29m1\nembryo', Pred_and_obs_GOs_max$DEG)
Pred_and_obs_GOs_max$DEG <- gsub('HSF29_m2:embryo_imbibed', 'HSF29m2\nembryo', Pred_and_obs_GOs_max$DEG)
Pred_and_obs_GOs_max$DEG <- gsub('MYBR32_m1:leaf', 'MYBR32m1\nleaf', Pred_and_obs_GOs_max$DEG)
Pred_and_obs_GOs_max$DEG <- gsub('MYBR32_m1:leaf', 'MYBR32m1\nleaf', Pred_and_obs_GOs_max$DEG)
Pred_and_obs_GOs_max$DEG <- gsub('WRKY2_m2:coleoptile_tip', 'WRKY2m2\ncoleoptile', Pred_and_obs_GOs_max$DEG)
Pred_and_obs_GOs_max$DEG <- gsub('WRKY8_m1:embryo_imbibed', 'WRKY8m1\nembryo', Pred_and_obs_GOs_max$DEG)
Pred_and_obs_GOs_max$DEG <- gsub('WRKY8_m2:embryo_imbibed', 'WRKY8m2\nembryo', Pred_and_obs_GOs_max$DEG)
Pred_and_obs_GOs_max$DEG <- gsub('WRKY82_m1:embryo_imbibed', 'WRKY82m1\nembryo', Pred_and_obs_GOs_max$DEG)
unique(Pred_and_obs_GOs_max$DEG)

DF_3.1 <- subset(Pred_and_obs_GOs_max, (DEG %in% c("TB1\nbuds_12DAP",
                                                   "TB1\nbuds_8DAP")))

#DF_3.2 <- subset(Pred_and_obs_GOs_max, DEG=='FEA4\tear')
DF_3.2 <- Pred_and_obs_GOs_max[grepl('FEA4', Pred_and_obs_GOs_max$DEG),]

DF_3.3 <- subset(Pred_and_obs_GOs_max, DEG %in% c("KN1\near",
                                                  "KN1\nleaf",
                                                  "KN1\nSAM",
                                                  "KN1\ntassel",
                                                  "MYBR32m1\nleaf",
                                                  "O2",
                                                  "WRKY82m1\nembryo"))

DF_3.4 <- subset(Pred_and_obs_GOs_max, DEG %in% c('RA1_ear',
                                                  'HSF13m1m2\nleaf',
                                                  "HSF18m1\nembryo",
                                                  "HSF20m1\nembryo",
                                                  "HSF24m3\ntassel",
                                                  "HSF24m4\ntassel",
                                                  "HSF29m1\nembryo",
                                                  "HSF29m2\nembryo",
                                                  "WRKY2m2\ncoleoptile",
                                                  'WRKY8m1\nembryo',
                                                  "WRKY8m2\nembryo"))



#  viridis(n=3, option = "D", alpha = 0.5, direction = -1)
ggplot(DF_3.1, aes(x=Class, y=GSS, fill=Class)) +
  geom_boxplot(notch = T, outlier.shape = NA) +
  scale_fill_viridis(discrete = T, option = "D", alpha = 0.5, direction = -1,
                     name='Method') + 
  scale_y_continuous(expand=c(0,0), limits = c(0, 0.6)) + 
  theme_pubclean() +
  ylab("GO semantic similarity")+
  xlab("") +
  theme(strip.text = element_text(size = 6, angle = 90), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text=element_text(size=10), 
        legend.position = 'none',
        text = element_text(size=10, family="Times")) + 
  facet_grid(. ~DEG, scales = "free_x") +
  stat_compare_means(label = "p.signif",
                     method = "wilcox.test",
                     label.y = c(0.35, 0.4, 0.45),
                     comparisons=list(c("Comm.Target", "Comm.Funct"), 
                                      c('Comm.Target', 'Network-based'),
                                      c('Comm.Funct', 'Network-based'))) -> Plot_5.1
#Plot_5.1
ggplot(DF_3.2, aes(x=Class, y=GSS, fill=Class)) +
  geom_boxplot(notch = T, outlier.shape = NA) +
  scale_fill_manual(values = c("#21908C80", "#44015480"), name='Method') + 
  scale_y_continuous(expand=c(0,0), limits = c(0, 0.6)) + 
  theme_pubclean() +
  ylab("")+
  xlab("") +
  theme(strip.text = element_text(size = 6, angle = 90), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text=element_text(size=10), 
        legend.position = 'none',
        text = element_text(size=10, family="Times")) + 
  facet_grid(. ~DEG, scales = "free_x") +
  stat_compare_means(label = "p.signif", 
                     label.y = 0.5,
                     method = "wilcox.test") -> Plot_5.2


ggplot(DF_3.3, aes(x=Class, y=GSS, fill=Class)) +
  geom_boxplot(notch = T, outlier.shape = NA) +
  #scale_fill_viridis(discrete = T, option = "D", alpha = 0.5, direction = -1, name='Method') + 
  scale_fill_manual(values = c("#FDE72580", "#21908C80"), name='Method') + 
  scale_y_continuous(expand=c(0,0), limits = c(0, 0.6)) + 
  theme_pubclean() +
  ylab("")+
  xlab("") +
  theme(strip.text = element_text(size = 6, angle = 90), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text=element_text(size=10), 
        legend.position = 'none',
        text = element_text(size=10, family="Times")) + 
  facet_grid(. ~DEG, scales = "free_x") +
  stat_compare_means(label = "p.signif", 
                     label.y = 0.5,
                     method = "wilcox.test")  -> Plot_5.3

ggplot(DF_3.4, aes(x=Class, y=GSS, fill=Class)) +
  geom_boxplot(notch = T, outlier.shape = NA) +
  #scale_fill_viridis(discrete = T, option = "D", alpha = 0.5, direction = -1, name='Method') + 
  scale_fill_manual(values = c("#FDE72580", "#44015480"), name='Method') + 
  scale_y_continuous(expand=c(0,0), limits = c(0, 0.6)) + 
  theme_pubclean() +
  ylab("")+
  xlab("") +
  theme(strip.text = element_text(size = 6, angle = 90), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text=element_text(size=10), 
        legend.position = 'none',
        text = element_text(size=10, family="Times")) + 
  facet_grid(. ~DEG, scales = "free_x")  -> Plot_5.4
#


## GSEA with GO predicted and DEGS logFCs
unique(GSEA_GOs_Freq$DEG)
GSEA_GOs_Freq$DEG <- gsub('FEA4:ear', 'FEA4\near', GSEA_GOs_Freq$DEG)
GSEA_GOs_Freq$DEG <- gsub("HSF13_m1m2:leaf", 'HSF13m1m2\nleaf', GSEA_GOs_Freq$DEG)
GSEA_GOs_Freq$DEG <- gsub("HSF18_m1:embryo_imbibed", 'HSF18m1\nembryo', GSEA_GOs_Freq$DEG)
GSEA_GOs_Freq$DEG <- gsub("HSF20_m1:embryo_imbibed", 'HSF20m1\nembryo', GSEA_GOs_Freq$DEG)
GSEA_GOs_Freq$DEG <- gsub("HSF29_m1:embryo_imbibed", 'HSF29m1\nembryo', GSEA_GOs_Freq$DEG)
GSEA_GOs_Freq$DEG <- gsub("HSF29_m2:embryo_imbibed", 'HSF29m2\nembryo', GSEA_GOs_Freq$DEG)
GSEA_GOs_Freq$DEG <- gsub('KN1:ear', 'KN1\near', GSEA_GOs_Freq$DEG)
GSEA_GOs_Freq$DEG <- gsub('KN1:leaf', 'KN1\nleaf', GSEA_GOs_Freq$DEG)
GSEA_GOs_Freq$DEG <- gsub('KN1:SAM', 'KN1\nSAM', GSEA_GOs_Freq$DEG)
GSEA_GOs_Freq$DEG <- gsub('KN1:tassel', 'KN1\ntassel', GSEA_GOs_Freq$DEG)
GSEA_GOs_Freq$DEG <- gsub("MYBR32_m1:leaf", 'MYBR32m1\nleaf', GSEA_GOs_Freq$DEG)
GSEA_GOs_Freq$DEG <- gsub('O2:endosperm', 'O2', GSEA_GOs_Freq$DEG)
GSEA_GOs_Freq$DEG <- gsub('RA1:ear_2mm', 'RA1\near', GSEA_GOs_Freq$DEG)
GSEA_GOs_Freq$DEG <- gsub('TB1:tiller_buds_12DAP', 'TB1\nbuds_12DAP', GSEA_GOs_Freq$DEG)
GSEA_GOs_Freq$DEG <- gsub('TB1:tiller_buds_8DAP', 'TB1\nbuds_8DAP', GSEA_GOs_Freq$DEG)
GSEA_GOs_Freq$DEG <- gsub("WRKY2_m2:coleoptile_tip", 'WRKY2m2\ncoleoptile', GSEA_GOs_Freq$DEG)
GSEA_GOs_Freq$DEG <- gsub("WRKY8_m1:embryo_imbibed", 'WRKY8m1\nembryo', GSEA_GOs_Freq$DEG)
GSEA_GOs_Freq$DEG <- gsub("WRKY8_m2:embryo_imbibed", 'WRKY8m2\nembryo', GSEA_GOs_Freq$DEG)
GSEA_GOs_Freq$DEG <- gsub("WRKY82_m1:embryo_imbibed", 'WRKY82m1\nembryo', GSEA_GOs_Freq$DEG)

## No pred GOs
# GSEA_GOs_Freq$DEG <- gsub("BAF6021_m1:tassel_stem", 'BAF6021m1\ntassel', GSEA_GOs_Freq$DEG)
# GSEA_GOs_Freq$DEG <- gsub("BAF6021_m2:tassel_stem", 'BAF6021m2\ntassel', GSEA_GOs_Freq$DEG)
# GSEA_GOs_Freq$DEG <- gsub("bZIP22:kernel", 'bZIP22\nkernel', GSEA_GOs_Freq$DEG)
# GSEA_GOs_Freq$DEG <- gsub("BZIP76_m2:leaf", 'BZIP76m2\nleaf', GSEA_GOs_Freq$DEG)
# GSEA_GOs_Freq$DEG <- gsub("BZIP76_m3:leaf", 'BZIP76m3\nleaf', GSEA_GOs_Freq$DEG)
# GSEA_GOs_Freq$DEG <- gsub("C3H42_m1:tassel_stem", 'C3H42m1\ntassel_stem', GSEA_GOs_Freq$DEG)
# GSEA_GOs_Freq$DEG <- gsub("E2F13_m1:coleoptile_tip", 'E2F13m1\ncoleoptile', GSEA_GOs_Freq$DEG)
# GSEA_GOs_Freq$DEG <- gsub("E2F19_m1:leaf", 'E2F19m1\nleaf', GSEA_GOs_Freq$DEG)
# GSEA_GOs_Freq$DEG <- gsub("E2F19_m2:leaf", 'E2F19m2\nleaf', GSEA_GOs_Freq$DEG)
# GSEA_GOs_Freq$DEG <- gsub("GRAS52_m1:embryo_imbibed", 'GRAS52m1\nembryo', GSEA_GOs_Freq$DEG)
# GSEA_GOs_Freq$DEG <- gsub("GRAS75_m1:embryo_imbibed", 'GRAS75m1\nembryo', GSEA_GOs_Freq$DEG)
# GSEA_GOs_Freq$DEG <- gsub("HSF24_m3:tassel", 'HSF24m3\ntassel', GSEA_GOs_Freq$DEG)
# GSEA_GOs_Freq$DEG <- gsub("HSF24_m4:tassel", 'HSF24m4\ntassel', GSEA_GOs_Freq$DEG)
# GSEA_GOs_Freq$DEG <- gsub("HSF6_m1:embryo_imbibed", 'HSF6m1\nembryo', GSEA_GOs_Freq$DEG)
# GSEA_GOs_Freq$DEG <- gsub("HSF6_m2:embryo_imbibed", 'HSF6m2\nembryo', GSEA_GOs_Freq$DEG)
# GSEA_GOs_Freq$DEG <- gsub("JMJ13_m4:tassel_stem", 'JMJ13m4\ntassel', GSEA_GOs_Freq$DEG)
# GSEA_GOs_Freq$DEG <- gsub("MYBR21_m1:embryo_imbibed", 'MYBR21m1\nembryo', GSEA_GOs_Freq$DEG)
# GSEA_GOs_Freq$DEG <- gsub("ORPHAN249_m2:embryo_imbibed", 'ORPHAN249m2\nembryo', GSEA_GOs_Freq$DEG)
# GSEA_GOs_Freq$DEG <- gsub("SBP20_m2:embryo_imbibed", 'SBP20m2\nembryo', GSEA_GOs_Freq$DEG)
# GSEA_GOs_Freq$DEG <- gsub("SBP20_m3:embryo_imbibed", 'SBP20m3\nembryo', GSEA_GOs_Freq$DEG)
# GSEA_GOs_Freq$DEG <- gsub("WRKY87_m1:embryo_imbibed", 'WRKY87m1\nembryo', GSEA_GOs_Freq$DEG)
# GSEA_GOs_Freq$DEG <- gsub("WRKY87_m2:embryo_imbibed", 'WRKY87m2\nembryo', GSEA_GOs_Freq$DEG)

GSEA_GOs_Freq[,'TFid'] <- chop(GSEA_GOs_Freq$TF, '[:]', 2)
unique(GSEA_GOs_Freq$DEG)

# method levels
GSEA_GOs_Freq$Class <- factor(GSEA_GOs_Freq$Class, 
                              levels = c("Comm.Target", "Comm.Function", "Network-based"))

GSEA_GOs_Freq$DEG <- factor(GSEA_GOs_Freq$DEG,
                            levels = c("TB1\nbuds_12DAP",'TB1\nbuds_8DAP', 'FEA4\near',
                                       'KN1\near', 'KN1\nleaf', 'KN1\nSAM', 'KN1\ntassel',
                                       'MYBR32m1\nleaf', 'O2', 'WRKY82m1\nembryo',
                                       'HSF13m1m2\nleaf', 'HSF18m1\nembryo',
                                       'HSF20m1\nembryo', 'HSF29m1\nembryo',
                                       'HSF29m2\nembryo', 'RA1\near',
                                       'WRKY2m2\ncoleoptile', 'WRKY8m1\nembryo',
                                       'WRKY8m2\nembryo'))

ggplot(GSEA_GOs_Freq, aes(x=Class, y=Fraction, fill=PvalClass)) +
  geom_bar(stat="identity") +
  scale_fill_brewer(palette="Paired", name='P-Value', 
                    label=c(expression("">=0.05), expression(""<=0.05))) +
  theme_pubclean() +
  scale_y_continuous(expand = c(0,0)) + 
  ylab("Method") +
  xlab("GOs Fraction") + 
  facet_grid(.~ DEG, scales = "free_x") +
  theme(strip.text = element_text(size = 6, angle=90), 
        legend.position = 'bottom',
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text=element_text(size=10), 
        text = element_text(size=10, family="Times")) -> Plot_7
#axis.text.y=element_blank(),
#axis.ticks.y=element_blank()


# Plot_layout(widths = c(0.4, 4))
Plot_5cd <- {Plot_5.1 + Plot_5.2 + Plot_5.3 + Plot_5.4 + 
    plot_layout(nrow = 1, widths = c(0.45, 0.15, 1.2, 1.1))}  / 
    Plot_7 +  plot_layout(heights = c(1, 1))


pdf("Plots/Plot_2cd.pdf", width=8, height=6)
print(Plot_5cd)
dev.off()

# Plot_defg <- {Plot_4 + 
#              {Plot_5.1 + Plot_5.2 + Plot_5.3 + plot_layout(nrow = 1, widths = c(0.4, 1.6, 1))} + 
#               plot_layout(widths = c(0.4, 4))} / 
#              {Plot_6 + Plot_7 + plot_layout(widths = c(1, 0.9))} + 
#               plot_layout(nrow = 2, heights = c(1, 0.8), )
# Plot_defg




##########################################################################################
# multi_TFs[multi_TFs$TF %in% DEGs_Peng,][multi_TFs[multi_TFs$TF %in% To_GO,]$ClassType == "GO",]
# Example_1 <- subset(CombinedAnnotation[grepl('GO', CombinedAnnotation$Class),], TF=="Zm00001d042777")
# Example_1 <- left_join(Example_1, GONAMES_DB, by=c('Annotation'="GO.ID"))
# View(Example_1)
# 
# 
# ## 4 
# Example_PWY_CommTarg <- subset(PWY_CommTarg, TF=="Zm00001d025770")
# Example_PWY_CommFunt <- subset(PWY_CommFunt, TF=="Zm00001d025770")
# Example_PWY_Network  <- subset(PWY_Network, TF=="Zm00001d025770")
# 
# Example_PWY_CommTarg[,"PWYname"] <- ReplaceNamePWY(Example_PWY_CommTarg$PWY)
# Example_PWY_CommFunt[,"PWYname"] <- ReplaceNamePWY(Example_PWY_CommFunt$PWY)
# Example_PWY_Network[,"PWYname"] <- ReplaceNamePWY(Example_PWY_Network$PWY)
# Example_PWY_CommTarg[,4]
# Example_PWY_CommFunt[,4]
# 
# PWYs_Freq <- as.data.table(table(PWYs[,c(1,3)]))
# PWYs_Freq <- subset(PWYs_Freq, N>0)
# 
# ggplot(PWYs_Freq , aes(x=Type, y=N))+
#   geom_jitter() +
#   theme_pubclean() +
#   ylab("PWY by TF")

##########################################################################################


