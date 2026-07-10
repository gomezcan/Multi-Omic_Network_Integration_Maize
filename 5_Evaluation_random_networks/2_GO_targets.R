suppressMessages(library(parallel))
suppressMessages(library(topGO))
suppressMessages(library(GeneOverlap))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(reshape2))
suppressMessages(library(dplyr))

args = commandArgs(trailingOnly=TRUE)

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


## Input file

## Syntenic genes 
Syntenic <- as_tibble(read.table("Zm.v4.synteny.genes.txt", h=T, stringsAsFactors = F, sep = '\t'))
Syntenic <- subset(Syntenic, gene_synteny=='syntenic')$gene_id

# GOs term annotations
background <- readMappings("synteny.ID_TopGO_V4_GRAMER.txt")
background_IDs <- as.character(unique(names(background)))

# List of TFs to test enrichment
AllTFs <- unique(fread("All_TFs.txt", h=F)$V1)

# Net
filenet = args[1]   # filenet = 'GRN.txt'
# filenet = 'MI0.5_Net.txt'
# Net name
filename = gsub('.txt', '', filenet)

# Read net input
Net <- unique(fread(filenet, header = F)[,1:2])
colnames(Net) <- c("Source", "Target")

# Defined number of TFs/modules
GenesList <- unique(Net$Source)
GenesList <- GenesList[(GenesList %in% AllTFs)]

GenesList <- GenesList[1:2]
Lgenes <- length(GenesList)

## Check already tested TFs/modules
GOs_Done <- list.files(path = 'BP_results_targets/', pattern = "*_GOs.txt")
GOs_Done <- GOs_Done[grepl(filename, GOs_Done)]

# # secong round for cluster not tested
GOs_Done <- gsub(paste0(filename,'_'), "", GOs_Done)
GOs_Done <- gsub("_GOs.txt","", GOs_Done)
# 
GenesList <- GenesList[!(GenesList %in% GOs_Done)]
GenesList <- unique(GenesList)
Lgenes <- length(GenesList)


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
    
    #SuperGO_Modules_Targ("Zm00001d006236", filename)
    mclapply(listtotest, function(x) SuperGO_Modules_Targ(x, filename), mc.cores=w)
  }
  else if (end > max){ break}
  else {
    r <- paste(Start, max, sep = "-")
    listtotest <- GenesList[Start:max]
    print(length(listtotest))
    
    w = length(listtotest)+1
    mclapply(listtotest, function(x) SuperGO_Modules_Targ(x, filename), mc.cores=w)
    
  }
}
