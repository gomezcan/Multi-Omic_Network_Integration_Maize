suppressMessages(library(parallel))
suppressMessages(library(topGO))
suppressMessages(library(GeneOverlap))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(reshape2))
suppressMessages(library(dplyr))
suppressMessages(library(gplots))

### Functions
## chop a string by a separator and return specified field
chop=function(myStr,mySep,myField){
  
  choppedString=sapply(strsplit(myStr,mySep),"[",myField)
  if(length(myField)>1){
    choppedString=apply(choppedString,2,function(x){paste0(x[!is.na(x)],collapse=mySep)})
  }
  return(choppedString)
}

GetGO <- function(degs, mutant){
  
  # Use a list of DEGs and the name of the mutant (string)
  # to identify GOs enriched. Required to have a background predefined
  # Define background based on genes in network

  # mutant: description of a category of source of the data
  # degs: set of genes to test enrichment
  GeneList <- factor(as.integer(background_IDs %in% degs))
  names(GeneList) <- background_IDs
  
  GOdata_BP <- new("topGOdata", ontology = "BP", 
                   allGenes = GeneList, 
                   annot = annFUN.gene2GO,
                   gene2GO = background)
  
  
  #### Define test ####
  test.stat <- new("classicCount", testStatistic = GOFisherTest, 
                   name = "Fisher test", nodeSize = 10)
  
  ### test enrichment 
  results_BP <- getSigGroups(GOdata_BP, test.stat)
  
  ######## Get Significant GOs ########  
  Res_DF_BP <- as_tibble(GenTable(GOdata_BP, classic = results_BP, topNodes = 1000, orderBy='Fis')) # save as dataframe
  Res_DF_BP["Mutant"] <- mutant # add Mutant column name
  Res_DF_BP$classic <- as.numeric(Res_DF_BP$classic)

  return(Res_DF_BP) # Return list of GOs-Stats and GeneID-GOs
}



SuperGO_CommomTarg <- function(tf, netname){
  ## used TF and its targets to test GO terms enrichment
  # 1. Select tf's Targets
  # 2. Make list file: degs
  # 3. Test enrichment
  # 4. Save top 1000 GOs table
  
  # Get network by TF
  network <- subset(Reduced_InteractionDB, Source==tf)
  
  Total_targets <- as_tibble(as.data.frame(table(unique(network)$Source), stringsAsFactors = F))
  colnames(Total_targets) <- c('TF', 'nTF') 
  
  # Genes input Net
  degs <- unique(network$Target)
  
  if(length(degs) > 10){
    # If targets largest than 
    # Include 3_later targets into paired comparison if required 
    out <- GetGO(degs, tf)
    
    # Save gene results GOs  
    out <- left_join(out, Total_targets, by=c("Mutant"="TF"))  
    write.table(out, paste0("BP_results_CommonTargets/Random.CommonTarget.", netname, "_", tf, "_GOs.txt"),
                sep = '\t', row.names = F, quote = F)
    
  }
  else { print(paste0(" Salado: ", tf, " .."))}
  
  # Genes_GOBP_
  return(print('.Done.'))
  
}


#########
args = commandArgs(trailingOnly=TRUE)

Netid <- args[1] # Random net id
#Netid <- 1

# PDI
GRN <- fread(paste0('RandomNets/Random.GRN.',Netid, '.txt'))

# PDI eQTL
eGRN <- fread(paste0('RandomNets/Random.eGRN.',Netid, '.txt'))

# CoExp
CEN <- fread(paste0('RandomNets/Random.CEN.',Netid, '.txt'))

# teQTL
GAN <- fread(paste0('RandomNets/Random.GAN.',Netid, '.txt'))

## Syntenic genes 
Syntenic <- as_tibble(read.table("Zm.v4.synteny.genes.txt", h=T, stringsAsFactors = F, sep = '\t'))
Syntenic <- subset(Syntenic, gene_synteny=='syntenic')$gene_id

# GOs term annotations
background <- readMappings("synteny.ID_TopGO_V4_GRAMER.txt")
background_IDs <- as.character(unique(names(background)))

# TFs targeted in analysis 
# List of TFs to test enrichment
AllTFs <- unique(fread("TFs_2_test_RNets.txt", h=F)$V1)

# Reduce nets to TFs of interest
GRN <- GRN[GRN$Source %in% AllTFs,]
eGRN <- eGRN[eGRN$Source %in% AllTFs,]
CEN <- CEN[CEN$Source %in% AllTFs,]
GAN <- GAN[GAN$Source %in% AllTFs,]

# Reduce nets to targets in syntenic
GRN <- GRN[GRN$Target %in% Syntenic,]
eGRN <- eGRN[eGRN$Target %in% Syntenic,]
CEN <- CEN[CEN$Target %in% Syntenic,]
GAN <- GAN[GAN$Target %in% Syntenic,]

# Identify common interactions

# Define list
InteractionDB <- list()

InteractionDB[["GRN"]] <- paste0(GRN$Source, "_", GRN$Target)
InteractionDB[["eGRN"]] <- paste0(eGRN$Source, "_", eGRN$Target)
InteractionDB[["CEN"]] <- paste0(CEN$Source, "_", CEN$Target)
InteractionDB[["GAN"]] <- paste0(GAN$Source, "_", GAN$Target)

# keep unique interactions
InteractionDB <- lapply(InteractionDB, unique)

# Interaction groups
InteractionDB_venn <- venn(InteractionDB)
InteractionDB_venn <- as.list(attr(InteractionDB_venn, "intersections"))
InteractionDB_venn <- plyr::ldply(InteractionDB_venn, data.table)
InteractionDB_venn <- as.data.table(InteractionDB_venn)

# Defined interactions with at least two lines of evidence 
Reduced_InteractionDB <- subset(InteractionDB_venn, !(.id %in% c("CEN", "GRN", "GAN", 'eGRN')))

# Table(Reduced_InteractionDB$.id)
# Split edges into Tf and targets
Reduced_InteractionDB[,"Source"] <- chop(Reduced_InteractionDB$V1, '[_]', 1)
Reduced_InteractionDB[,"Target"] <- chop(Reduced_InteractionDB$V1, '[_]',2)

# BP_results_CommonTargets/Random.CommonTarget
GOs_Done <- list.files(path = 'BP_results_CommonTargets/', pattern = "*_GOs.txt")

# Get tf/modules from the same 'filename' network
#Netid = 1049
filname=paste0("Random.CommonTarget.", Netid,"_")

GOs_Done <- GOs_Done[grepl(filname, GOs_Done)]
GOs_Done <- gsub("_GOs.txt","", GOs_Done)
GOs_Done <- gsub(filname,"", GOs_Done)

AllTFs <- AllTFs[!(AllTFs %in% GOs_Done)]

lapply(AllTFs, function(x) SuperGO_CommomTarg(x, Netid))





