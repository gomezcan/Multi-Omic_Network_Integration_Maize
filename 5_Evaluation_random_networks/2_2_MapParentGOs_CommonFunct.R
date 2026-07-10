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

#conda activate /maindisk/fabio/miniconda3
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


ReduceGOs <- function(GOsDB){
  #library(rrvgo)
  tf <- unique(GOsDB$Index)
  
  GOsDB <- unique(subset(GOsDB, Index==tf))
  scores <- setNames(-log10(GOsDB$FDR), GOsDB$GO.ID) 
  GO_vector = GOsDB$GO.ID
  
  if(length(GO_vector) >2) {
    # Semantic similarity
    cat(paste0(tf,'\n'))
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

    #head(reducedTerms)
    #
    reducedTerms <- unique(reducedTerms[,c("TF",'parent')])
    colnames(reducedTerms)[2] <- c("GO.ID")
    file <- paste0("BP_RandomGOs_CommFunct_Parents/RedRandom_CommFunct_", tf, ".txt")
    # RedRandom_CommFunct_Net_100_Zm00001d033719_GRN_CEN.txt
    write.table(reducedTerms, file, row.names = F, quote = F, sep = '\t')
    
    return(reducedTerms)
  }
  if(length(GO_vector) <= 2 & length(GO_vector) > 0) {
    
    reducedTerms <- GOsDB[,c("Index", 'GO.ID')]
    colnames(reducedTerms)[1] <- 'TF'
    file <- paste0("BP_RandomGOs_CommFunct_Parents/RedRandom_CommFunct_", tf, ".txt")
    write.table(reducedTerms, file, row.names = F, quote = F, sep = '\t')
    return(reducedTerms)
  }
  else {cat(paste0('.. Salado: ', tf,'\n'))}
}

# list of files
List_ComFuncRam <- list.files(path = 'BP_GSS_CommonFunction_random/', 
                              pattern = "^GSS.*")
#length(List_ComFuncRam)

# Reduce files list to only TF tested
# keep GO files associates with TF "tf"
# tf="Zm00001d005892"
List_ComFuncRam <- List_ComFuncRam[grepl(tf, List_ComFuncRam)]
#List_ComFuncRam <- List_ComFuncRam[1:5]
# length(List_ComFuncRam)

if(length(List_ComFuncRam) >=1){
  
  ## list by GOs by netid to test
  netids <- paste0(chop(List_ComFuncRam, '[.]', 2),'.',
                   chop(List_ComFuncRam, '[.]', 4))
  netids <- gsub('_','.',netids)
  
  # Check for combination of nets and TFs already tested
  netsDone <- list.files(path = 'BP_RandomGOs_CommFunct_Parents/', pattern = "^RedRandom*")
  # length(netsDone) 
  
  # Filter for method type
  netsDone <- netsDone[grepl('CommFunct', netsDone)]
  # length(netsDone)
  
  # Filter for TF
  netsDone <- netsDone[grepl(tf, netsDone)] 
  # length(netsDone)
  
  # netids of nets already tested
  netsDoneID <- paste0(chop(netsDone, '[_]', 4),'.',
                       chop(netsDone, '[_]', 6),'.',
                       chop(netsDone, '[_]', 7))
  
  netsDoneID <- gsub('.txt', '', netsDoneID)
  
  # Make mask with net_comparison already tested
  mask <- !(netids %in% netsDoneID)
  # length(List_ComFuncRam)
  List_ComFuncRam <- List_ComFuncRam[mask]
  # length(List_ComFuncRam)
  
  # New list of nets
  netids <- paste0(chop(List_ComFuncRam, '[.]', 2),'.',
                   chop(List_ComFuncRam, '[.]', 4))
  netids <- gsub('_','.',netids)
  
  # Read files
  List_ComFuncRam <- lapply(List_ComFuncRam, function(x) fread(paste0("BP_GSS_CommonFunction_random/", x)))
  names(List_ComFuncRam) <- as.character(netids)
  
  # Filter out GOs with low GSS
  List_ComFuncRam <- lapply(List_ComFuncRam, function(x) subset(x, GSS >=  0.6))
  mask <- unlist(lapply(List_ComFuncRam, function(x) nrow(x) >0))
  #table(mask)
  # Filter out empty DFs
  List_ComFuncRam <- List_ComFuncRam[mask]
  # length(List_ComFuncRam)
  
  # set FDR.net wwith generic names
  List_ComFuncRam <- lapply(List_ComFuncRam, function(x) {
    colnames(x)[c(6,9)] <- c('FDR1', "FDR2") 
    return(x)}
    )
  
  # filter out GOs terms which do not pass FDR filters
  List_ComFuncRam <- lapply(List_ComFuncRam, function(x) {
    subset(x, FDR1 <= 0.1 & FDR2 <= 0.1)
    }
  )
  mask <- unlist(lapply(List_ComFuncRam, function(x) nrow(x) >0))
  
  # Filter out empty DFs
  List_ComFuncRam <- List_ComFuncRam[mask]
  
  # Create name index
  List_ComFuncRam <- lapply(List_ComFuncRam, function(x) {
    x[,'Index'] <- gsub('GO.', '', paste0(colnames(x)[1:2], collapse = '_'))
    return(x)
  }
  )
  
  # Set
  List_ComFuncRam <- lapply(List_ComFuncRam, setNames, c('GO1', "GO2", "GSS",
                                                         'Term1', 'nGO1', 'FDR1',
                                                         'Term2', 'nGO2', 'FDR2', 
                                                         'Index'))
  #
  List_ComFuncRam <- rbindlist(List_ComFuncRam, idcol = T)
  List_ComFuncRam <- List_ComFuncRam[,-c('Term1', 'Term2')]
  List_ComFuncRam[,'TF'] <- tf
  #  Net.2145.RedRandom_CommFunct_Zm00001d046925.txt
  List_ComFuncRam[,'Index'] <- paste('Net',
                                     chop(List_ComFuncRam$.id, '[.]',1),
                                     List_ComFuncRam$TF,
                                     List_ComFuncRam$Index, sep='_')
  # RedRandom_CommFunct_Net_100_Zm00001d033719_GRN_CEN.txt
  #ist_ComFuncRam

  # Extract each pair of GOs that pass GSS filter
  ## Part 1
  List_ComFuncRam %>%
    dplyr::select(GO1, FDR1, TF, Index) -> DF1
  colnames(DF1) <- c('GO.ID',"FDR", "TF", "Index")
  
  ## Part 2
  List_ComFuncRam %>%
    dplyr::select(GO2, FDR2, TF, Index) -> DF2
  colnames(DF2) <- c('GO.ID',"FDR", "TF", "Index")
  
  # Combine random results and keep GO sig. enriched after FDR
  List_ComFuncRam <- unique(rbind(DF1, DF2))
  
  netids <- unique(List_ComFuncRam$Index)
  ## keep net testable
  List_ComFuncRam <- lapply(netids, function(x) subset(List_ComFuncRam, Index==x))
  names(List_ComFuncRam) <- netids
  # check empty DFs
  mask <- unlist(lapply(List_ComFuncRam, function(x) nrow(x) >0))
  #table(mask)
  # Filter out empty DFs
  List_ComFuncRam <- List_ComFuncRam[mask]
  netids <- names(List_ComFuncRam)
  # length(List_ComFuncRam)
  
  ## Reduce GO terms by random netid
  lapply(List_ComFuncRam, function(x) ReduceGOs(x))
  # List_ComTargRam_red <- lapply(List_ComTargRam, function(x) ReduceGOs(x, tf))
  # mask <- unlist(lapply(List_ComTargRam_red, function(x) is.data.frame(x)))
  # List_ComTargRam_red <- List_ComTargRam_red[mask]
  # List_ComTargRam_red <- rbindlist(List_ComTargRam_red, idcol = F)
  
  # write.table(List_ComTargRam_red, paste0("BP_RandomGOs_Parents/RedRandom_Commtarg_", tf, ".txt"),
  #             row.names = F, quote = F, sep = '\t')
  
}
