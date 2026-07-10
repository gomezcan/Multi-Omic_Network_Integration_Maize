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
#tf <- "Zm00001d005578"

#conda activate /maindisk/fabio/miniconda3
# Pre-calculated semantic similarity for from org.Zmays.eg.db for BP term
Zm.GOSemSim.BP <- readRDS("Zm.GOSemSim.BP.rds")

# rm(list=ls())
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
    cat(paste0(tf,' :',netId,'\n'))
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
    #head(reducedTerms)
    #
    reducedTerms <- unique(reducedTerms[,c("TF",".id", 'parent','size')])
    colnames(reducedTerms)[c(3,4)] <- c("GO.ID",'Size')
    file <- paste0("BP_RandomGOs_CommTarg_Parents/Net.", netId,".RedRandom_Commtarg_", tf, ".txt")
    write.table(reducedTerms, file, row.names = F, quote = F, sep = '\t')
    
    return(reducedTerms)
  }
  if(length(GO_vector) <= 2 & length(GO_vector) > 0) {
    reducedTerms <- GOsDB[,c('TF', ".id", 'GO.ID', 'Annotated')]
    colnames(reducedTerms)[4] <- 'Size'
    file <- paste0("BP_RandomGOs_CommTarg_Parents/Net.", netId,".RedRandom_Commtarg_", tf, ".txt")
    write.table(reducedTerms, file, row.names = F, quote = F, sep = '\t')
    return(reducedTerms)
  }
    else {cat(paste0('.. Salado: ', tf,'\n'))}
}


# list of files
List_ComTargRam <- list.files(path = 'BP_results_CommonTargets_random/', 
                              pattern = "^Random.CommonTarget.*")

#Full_List_ComTargRam[1:2]
# reduce files list to only TF tested
# keep GO files associates with TF "tf"
List_ComTargRam <- List_ComTargRam[grepl(tf, List_ComTargRam)]
# length(List_ComTargRam)
# List_ComTargRam <- List_ComTargRam[1:100]
#List_ComTargRam[1]
if(length(List_ComTargRam) >=1){
  
  ## list by GOs by netid
  netids <- gsub("Random.CommonTarget.", "", chop(List_ComTargRam, '[_]', 1))
  
  # Check for combination of nets and TFs already tested
  netsDone <- list.files(path = 'BP_RandomGOs_CommTarg_Parents/', pattern = "^Net.*")
  #length(netsDone)
  netsDone <- netsDone[grepl(tf, netsDone)]
  #length(netsDone)
  netsDone <- chop(netsDone, '[.]',2)
  
  #
  mask <- !(netids %in% netsDone)
  # selected untested samples
  List_ComTargRam <- List_ComTargRam[mask]
  #netids <- netids[mask]
  
  netids <- gsub("Random.CommonTarget.", "", chop(List_ComTargRam, '[_]', 1))
  # Read files
  List_ComTargRam <- lapply(List_ComTargRam, function(x) fread(paste0("BP_results_CommonTargets_random/", x)))
  names(List_ComTargRam) <- netids
  List_ComTargRam <- rbindlist(List_ComTargRam, idcol = T)
  colnames(List_ComTargRam)[8] <- 'TF'
  
  # Combine random results and keep GO sig. enrichd after FDR
  List_ComTargRam %>% 
    dplyr::group_by(.id) %>% # group by network
    dplyr::mutate("FDR" = p.adjust(classic, method = 'fdr')) %>%
    dplyr::filter(FDR <= 0.1)  -> List_ComTargRam
  #
  netids <- unique(List_ComTargRam$.id)
  
  List_ComTargRam <- lapply(netids, function(x) subset(List_ComTargRam, .id==x))
  names(List_ComTargRam) <- netids
  
  
  # Reduce GO terms by random netid
  lapply(List_ComTargRam, function(x) ReduceGOs(x, tf))
  # List_ComTargRam_red <- lapply(List_ComTargRam, function(x) ReduceGOs(x, tf))
  # mask <- unlist(lapply(List_ComTargRam_red, function(x) is.data.frame(x)))
  # List_ComTargRam_red <- List_ComTargRam_red[mask]
  # List_ComTargRam_red <- rbindlist(List_ComTargRam_red, idcol = F)
  
  # write.table(List_ComTargRam_red, paste0("BP_RandomGOs_Parents/RedRandom_Commtarg_", tf, ".txt"),
  #             row.names = F, quote = F, sep = '\t')
  
}