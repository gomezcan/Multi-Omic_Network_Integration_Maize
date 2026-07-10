suppressMessages(library(tidyverse))
suppressMessages(library(viridis))
suppressMessages(library(igraph))
suppressMessages(library(scales))
suppressMessages(library(parallel))
suppressMessages(library(data.table))
library(parallel)
suppressMessages(library(topGO))

#####################################################
#############      Functions            #############
#####################################################

RewireNet <- function(net_igraph) {
  ####  
  # This function take as df_net-like to convert it into an igraph 
  # object. rewire it, and return a random network.
  ####
  # Rewire network with similar degree
  igraph_R <- rewire(net_igraph, with = keeping_degseq(loops = FALSE, niter = vcount(net_igraph)*10000))
  # Get  igraph DF
  out <- as_data_frame(igraph_R, what = "edges")
  colnames(out) <- c("Source", "Target")
  return(out)
}


GetGO <- function(degs, mutant, netindex){
  
  # Use a list of DEGs and the name of the mutant (string)
  # to identify GOs enriched. Required to have a background predefined
  #
  # mutant: description of a category of source of the data
  # degs: set of genes to test enrichment
  GeneList <- factor(as.integer(background_IDs %in% degs))
  names(GeneList) <- background_IDs
  
  GOdata_BP <- new("topGOdata", ontology = "BP", allGenes = GeneList, annot = annFUN.gene2GO, gene2GO = background)
  #GOdata_MF <- new("topGOdata", ontology = "MF", allGenes = GeneList, annot = annFUN.gene2GO, gene2GO = background)
  
  #### Define test ####
  test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test", nodeSize = 10)
  
  ### test enrichment 
  results_BP <- getSigGroups(GOdata_BP, test.stat)
  #results_MF <- getSigGroups(GOdata_MF, test.stat)
  
  ### save pdf Graph
  #namepdf=paste("GOs_Plots/GO.BP_",mutant, "", sep = "")
  #printGraph(GOdata_BP, results_BP, firstSigNodes=20,  fn.prefix = namepdf, useInfo = "def", pdfSW = TRUE) #
  
  ######## Get Significant GOs ########  
  # save as dataframe
  Res_DF_BP <- as_tibble(GenTable(GOdata_BP, classic = results_BP, topNodes = 500, orderBy='Fis')) 
  Res_DF_BP["Mutant"] <- mutant # add Mutant column name
  Res_DF_BP$classic <- as.numeric(Res_DF_BP$classic)
  
  
  GOsname <-   paste0("BP_rNet.teQTL.", netindex,'.GOs.', mutant, ".txt")
  
  # GOs results
  write.table(Res_DF_BP, paste0("RandomNetsData/teQTL_GOsV2/", GOsname), quote = F, row.names=F, sep = '\t')
  
}

#####################################################
#############     Data input            #############
#####################################################

# TF and CoReg
TF_CoR <- as_tibble(read.table("Data/Annotations/TF_CoR_Mazie.txt", h=T, stringsAsFactors = F)) 

# GOs term annotations
background <- readMappings("synteny.ID_TopGO_V4_GRAMER.txt")
background_IDs <- as.character(unique(names(background)))

#print(background_IDs)

# teQTL
teQTL <- as_tibble(unique(fread("../Fig_transeQTL/teQTL_NetworkFinal.10_11_2021.txt")))
colnames(teQTL)[1] <- "Source"

# TFs to be tested
#TFslist <- as.character(read.table("TFs2Test_PDI.txt", h=T)$TFs)
TFslist <- unique(as.character(teQTL$Source[teQTL$Source %in% TF_CoR$GeneID]))

TFslist <- c("Zm00001d001945", "Zm00001d009599", "Zm00001d010634", "Zm00001d011953", "Zm00001d012916", 
             "Zm00001d014995", "Zm00001d015407", "Zm00001d015421", "Zm00001d016838", "Zm00001d017726", 
             "Zm00001d020267", "Zm00001d024324", "Zm00001d031182", "Zm00001d033719", "Zm00001d042267",
             "Zm00001d050195", "Zm00001d050781")
#TFslist <- unique(as.character(teQTL$Source[teQTL$Source %in% unique(PDI$Source)]))


print(" .. TFs to test ..")
print(length(TFslist))

# make igraph object from true network
NetTFs_igraph <- graph_from_data_frame(teQTL, directed = T) 

set.seed(123)


RamdonSample <- function(times){
  
  # Get random network
  Rnet <- RewireNet(NetTFs_igraph)
  write.table(Rnet, paste0("RandomNetsData/teQTL_NetsV2/Random_teQTL.",times,".txt"), 
              sep = '\t', row.names = F, quote = F)
  
  for (t in TFslist){
    tgs <- as.character(subset(Rnet, Source==t)$Target)
    #
    print(length(tgs))
    GetGO(tgs, t, times) # genes to test; TF; net_index=times
  }
  
}

w=50 # Size of range to test

print(".. Ready to start ..")
Samples = 1000

for (i in seq(0, Samples, w)){
  max=Samples
  Start=i+1
  end=i+w
  if (end<max){
    #
    r.nets <- seq(Start,end,1)
    #print(r.nets)
    mclapply(r.nets, RamdonSample, mc.cores=w)
    r <- paste(Start, end, sep = "-")
    print(paste0(" .. Done ", r, " Samples .."))
  }
  else{
    r <- paste(Start, max, sep = "-")
    r.nets <- seq(Start, max, 1)
    #print(r.nets)
    mclapply(r.nets, RamdonSample, mc.cores=w)
    print(paste0(" .. Done ", r, " Samples .."))
  }
}
