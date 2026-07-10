library(data.table)
suppressMessages(library(tidyverse))


Get_spearman_id1_id2 <- function(id1, DB){
  # gene id from maize copy 2
  id2 <- DB[DB$GeneID1 == id1,]$GeneID2
  # MR_MI.pecanpy.Zm00001d036331.txt
  file1 <- paste0('MR_edgesDB_syntenic/MR_MI.pecanpy.',id1,'.txt')
  file2 <- paste0('MR_edgesDB_syntenic/MR_MI.pecanpy.',id2,'.txt')
  #
  mr.mi1 <- fread(file1)[,1:4]
  geneorder <- mr.mi1$GeneID2
  
  #
  mr.mi1 <- unlist(split(mr.mi1$MR, mr.mi1$GeneID2))
  
  # Read DFs from paralog paris: id2
  mr.mi2 <- lapply(file2, function(x) fread(x)[,1:4])
  names(mr.mi2) <- id2
  
  # keep id2 as vector with geneid and MR values
  mr.mi2 <- lapply(mr.mi2, function(x) unlist(split(x$MR, x$GeneID2)))
  names(mr.mi2) <- id2
  
  # Define vector 1: id1
  v1 <- mr.mi1[geneorder]
  
  # Define empty list to save results
  spearmanResults <- list()
  
  for(id in id2){
    # To keep gene order
    v2 <- mr.mi2[[id]][geneorder] 
    
    # spearman cor of MR profile
    val <- cor(v1, v2, method = 'spearman')
    spearmanResults <- c(spearmanResults, val)
  }
  names(spearmanResults) <- id2
  
  t(as.data.frame(spearmanResults)) %>%
    as.data.frame() %>%
    rownames_to_column("GeneID2") -> spearmanResults
  
  spearmanResults[,'GeneID1'] <- id1
  spearmanResults <- spearmanResults[,c("GeneID1", "GeneID2", "V1")]
  colnames(spearmanResults)[3] <- "SCC"
  
  return(spearmanResults)
}



# Read random pairs file index
args = commandArgs(trailingOnly=TRUE)
index <- args[1]
# index <- 1

# file name
file <- paste0("ParalogsRandom/RandomParalogsPairs_",index, ".txt")

DFrandom <- fread(file, sep = ',')

# calculate SCC for each pair of paralogs
DFrandom_SCC <- lapply(unique(DFrandom$GeneID1), function(x) Get_spearman_id1_id2(x, DFrandom))

# Combined df from list
DFrandom_SCC <- as.data.table(rbindlist(DFrandom_SCC, idcol = F))

# Save results by random set
fwrite(DFrandom_SCC,
       paste0('SCC_PatalogsRandom/SSC_RandomParalogsPairs.', index, '.txt'),
       quote = F, row.names = F)


cat(paste0(" .. donde: ", index, " ..\n"))






