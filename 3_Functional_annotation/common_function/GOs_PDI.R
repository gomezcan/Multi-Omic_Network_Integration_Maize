suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(scales))
suppressMessages(library(topGO))
suppressMessages(library(Rgraphviz))



###################################
#####        Functions        #####
###################################

ReplaceName <- function(ids){
  
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

GetGO <- function(degs, mutant, class){
  
  # Use a list of DEGs and the name of the mutant (string)
  # to identify GOs enriched. Required to have a background predefined
  #
  # mutant: description of a category of source of the data
  # degs: set of genes to test enrichment
  print(length(degs))
  GeneList <- factor(as.integer(background_IDs %in% degs))
  names(GeneList) <- background_IDs
  
  GOdata_BP <- new("topGOdata", ontology = "BP", allGenes = GeneList, annot = annFUN.gene2GO, gene2GO = background)
  GOdata_MF <- new("topGOdata", ontology = "MF", allGenes = GeneList, annot = annFUN.gene2GO, gene2GO = background)
  #GOdata_CC <- new("topGOdata", ontology = "CC", allGenes = GeneList, annot = annFUN.gene2GO, gene2GO = background)
  
  #### Define test ####
  test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test", nodeSize = 10)
  
  ### test enrichment 
  results_BP <- getSigGroups(GOdata_BP, test.stat)
  results_MF <- getSigGroups(GOdata_MF, test.stat)
  
  ### save pdf Graph
  #namepdf=paste("GOs_Plots/GO.BP_",mutant, "", sep = "")
  #printGraph(GOdata_BP, results_BP, firstSigNodes=20,  fn.prefix = namepdf, useInfo = "def", pdfSW = TRUE) #
  
  ######## Get Significant GOs ########  
  # save as dataframe
  Res_DF_BP <- as_tibble(GenTable(GOdata_BP, classic = results_BP, topNodes = 500, orderBy='Fis')) 
  Res_DF_BP["Mutant"] <- mutant # add Mutant column name
  Res_DF_BP$classic <- as.numeric(Res_DF_BP$classic)
  # save as dataframe
  Res_DF_MF <- as_tibble(GenTable(GOdata_MF, classic = results_MF, topNodes = 500, orderBy='Fis')) 
  Res_DF_MF["Mutant"] <- mutant # add Mutant column name
  Res_DF_MF$classic <- as.numeric(Res_DF_MF$classic)
  
  
  ##### get all GOs and their genes from the topGO result #####
  gs <- genesInTerm(GOdata_BP) # list genes by GO
  
  ##### get all GOs and their genes from the topGO result #####
  gs_mf <- genesInTerm(GOdata_MF) # list genes by GO
  # 
  ## Get only my Differential expressed genes
  ANOTATION = lapply(gs,function(x) x[x %in% degs]) 
  
  ## Get only my Differential expressed genes
  ANOTATION_mf = lapply(gs_mf, function(x) x[x %in% degs]) 
  
  ### Get only the GO's located in the result of topGO in Res_DF
  DF_GO_Genes <- ANOTATION[Res_DF_BP$GO.ID] # list
  
  ### Get only the GO's located in the result of topGO in Res_DF
  DF_GO_Genes_mf <- ANOTATION_mf[Res_DF_MF$GO.ID] # list
  
  ## Transform list to a dataframe.
  DF_GO_Genes = list_to_DF(DF_GO_Genes)
  DF_GO_Genes <- left_join(DF_GO_Genes, Res_DF_BP[,c(1,2,6)], by='GO.ID')   # left join to add GO info
  DF_GO_Genes["Mutant"] <- mutant # add Mutant column name
  
  ## Transform list to a dataframe.
  DF_GO_Genes_MF = list_to_DF(DF_GO_Genes_mf)
  DF_GO_Genes_MF <- left_join(DF_GO_Genes_MF, Res_DF_MF[,c(1,2,6)], by='GO.ID')   # left join to add GO info
  DF_GO_Genes_MF["Mutant"] <- mutant # add Mutant column name
  
  Genesname <- paste0(class,'.Genes.', mutant, ".txt")
  GOsname <-   paste0(class,'.GOs.', mutant, ".txt")
  
  # GOs results
  write.table(Res_DF_BP, paste0("BP_results/BP_", GOsname), quote = F, row.names=F, sep = '\t')
  write.table(Res_DF_MF, paste0("MF_results/MF_", GOsname), quote = F, row.names=F, sep = '\t')
  
  write.table(DF_GO_Genes, paste0("BP_results/BP_", Genesname), quote = F, row.names=F, sep = '\t')
  write.table(DF_GO_Genes_MF, paste0("MF_results/MF_", Genesname), quote = F, row.names=F, sep = '\t')
  
  
}

###################################

# GOs term annotations
background <- readMappings("synteny.ID_TopGO_V4_GRAMER.txt")
background_IDs <- as.character(unique(names(background)))

#print(background_IDs)

# PDI
PDI <- as_tibble(unique(fread("../Fig_PDI/PDI_NetworkFinal.10_11_2021.txt")[,c(2,3)]))
colnames(PDI)[1] <- "Source"


################################################################
########   Count common TFs and number of associations  ######## 
################################################################

#print(PDI[1:10,])

TFsPDI <- unique(PDI$Source)
TFsPDI <- TFsPDI[!(TFsPDI %in% c("Zm00001d002025", "Zm00001d005737"))] # TFs with low number of targets

print(" .. TFs to test ..")
print(length(TFsPDI))


for (t in TFsPDI){
  tgs <- as.character(subset(PDI, Source==t)$Target)
  #
  GetGO(tgs, t, "PDI")
}

