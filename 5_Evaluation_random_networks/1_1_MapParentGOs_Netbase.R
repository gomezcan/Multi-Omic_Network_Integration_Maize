suppressMessages(library(parallel))
suppressMessages(library(data.table))
suppressMessages(library(purrr))
suppressMessages(library(rrvgo))
suppressMessages(library(reshape2))
suppressMessages(library(GOSemSim))
suppressMessages(library(dplyr))
suppressMessages(library(org.Zmays.eg.db))

args = commandArgs(trailingOnly=TRUE)

## Reduce files to targeted TFs
tf <- args[1]
#tf <- "Zm00001d018571"

#system("conda activate /maindisk/fabio/miniconda3")

# Pre-calculated semantic similarity for from org.Zmays.eg.db for BP term
Zm.GOSemSim.BP <- readRDS("Zm.GOSemSim.BP.rds")

###############################################
###               Functions                 ###
###############################################
chop=function(myStr,mySep,myField){
  ## chop a string by a separator and return specified field
  choppedString=sapply(strsplit(myStr,mySep),"[",myField)
  if(length(myField)>1){
    choppedString=apply(choppedString,2,function(x){paste0(x[!is.na(x)],collapse=mySep)})
  }
  return(choppedString)
}

#GOsDB = List_ComTargRam[[1]]
ReduceGOs <- function(GOsDB, tf){
  #library(rrvgo)
  GOsDB <- unique(subset(GOsDB, TF==tf))
  scores <- setNames(-log10(GOsDB$FDR), GOsDB$GO.ID) 
  GO_vector = GOsDB$GO.ID
  netId <- unique(GOsDB$.id)
  
  if(length(GO_vector) >2) {
    # Semantic similarity
    simMatrix <- calculateSimMatrix(GO_vector,  
                                    orgdb=org.Zmays.eg.db,  
                                    ont="BP", 
                                    semdata=Zm.GOSemSim.BP,
                                    method="Wang")
    
    ##
    # check for TFs with few GOs terms
    ##
    # Reduce term
    reducedTerms <- reduceSimMatrix(simMatrix, scores, keytype="GENENAME",
                                    threshold=0.7, orgdb=org.Zmays.eg.db)
    
    # treemapPlot(reducedTerms)
    reducedTerms[,"TF"] <- tf
    reducedTerms[,".id"] <- netId
    #
    reducedTerms <- unique(reducedTerms[,c("TF",".id", 'parent')])
    colnames(reducedTerms)[c(3)] <- c("GO.ID")
    file <- paste0("BP_RandomGOs_NetBase_Parents/Net.", netId,".RedRandom_NetBased_", tf, ".txt")
    # BP_results_NetBase_random
    write.table(reducedTerms, file, row.names = F, quote = F, sep = '\t')
    return(cat(paste0(tf,': ',netId,'\n')))
  }
  if(length(GO_vector) <= 2 & length(GO_vector) > 0) {
    reducedTerms <- GOsDB[,c('TF', ".id", 'GO.ID')]
    file <- paste0("BP_RandomGOs_NetBase_Parents/Net.", netId,".RedRandom_NetBased_", tf, ".txt")
    write.table(reducedTerms, file, row.names = F, quote = F, sep = '\t')
    return(cat(paste0(tf,': ',netId,'\n')))
  }
  else {cat(paste0('.. Salado: ', tf,'\n'))}
}

# list of files
List_NetBase <- list.files(path = 'BP_results_NetBase_random/', 
                              pattern = "^Random.FullMR.*")

## Reduce files list to only TF tested
# keep GO files associates with TF "tf"
List_NetBase <- List_NetBase[grepl(tf, List_NetBase)]
# length(List_NetBase)
# List_NetBase <- List_NetBase[1:100]

if(length(List_NetBase) >=1){
  ## list by GOs by netid to test
  netids <- gsub('Random.FullMR.','', chop(List_NetBase, '[_]', 1))

  # Check for combination of nets and TFs already tested
  netsDone <- list.files(path = 'BP_RandomGOs_NetBase_Parents/', pattern = "^Net.*")
  # length(netsDone) 
  
  # Filter for method type
  netsDone <- netsDone[grepl('NetBased', netsDone)]
  # length(netsDone)
  
  # Filter for TF
  netsDone <- netsDone[grepl(tf, netsDone)] 
  # length(netsDone)
  
  # netids of nets already tested
  netsDoneID <- chop(netsDone, '[.]', 2)
  
  # Make mask with net_comparison already tested
  mask <- !(netids %in% netsDoneID)
  # length(List_ComFuncRam)
  List_NetBase <- List_NetBase[mask]
  
  netids <- netids[mask]
  
  # Read files
  List_NetBase <- lapply(List_NetBase, function(x) fread(paste0("BP_results_NetBase_random/", x)))
  names(List_NetBase) <- netids
  #length(List_NetBase)
  
  # filter out GOs terms which do not pass FDR filters
  List_NetBase <- rbindlist(List_NetBase, idcol = T)
  
  # Combine random results and keep GO Sig. enriched after FDR
  List_NetBase %>% 
    dplyr::group_by(.id) %>% # group by network
    dplyr::mutate("FDR" = p.adjust(classic, method = 'fdr')) %>%
    dplyr::filter(FDR <= 0.1)  -> List_NetBase
  
  ## list by GOs by netid
  netids <- unique(List_NetBase$.id)
  
  List_NetBase <- lapply(netids, function(x) subset(List_NetBase, .id==x))
  names(List_NetBase) <- netids
  
  # Reduce GO terms by random netid
  lapply(List_NetBase, function(x) ReduceGOs(x, tf))
  # List_ComTargRam_red <- lapply(List_ComTargRam, function(x) ReduceGOs(x, tf))
  # mask <- unlist(lapply(List_ComTargRam_red, function(x) is.data.frame(x)))
  # List_ComTargRam_red <- List_ComTargRam_red[mask]
  # List_ComTargRam_red <- rbindlist(List_ComTargRam_red, idcol = F)
  
  # write.table(List_ComTargRam_red, paste0("BP_RandomGOs_Parents/RedRandom_Commtarg_", tf, ".txt"),
  #             row.names = F, quote = F, sep = '\t')
  
}
