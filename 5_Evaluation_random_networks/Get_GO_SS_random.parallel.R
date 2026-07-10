suppressMessages(library(topGO))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(reshape2))
suppressMessages(library(scales))
suppressMessages(library(AnnotationHub))
suppressMessages(library(GOSemSim))
library(parallel)
# Note: this version used 



###################################
#####        Functions        #####
###################################

ReadGOs <- function(list_DEGs) {
  # Uses a list of files to read
  names <- sapply(strsplit(list_DEGs, split='/', fixed=TRUE), `[`, 2) # add methods label
  names <- gsub('GOs.', "", names)
  names <- gsub('.txt', "", names)
  
  
  go_df <- lapply(list_DEGs, fread)
  names(go_df) <- names
  
  go_df <- as_tibble(rbindlist(go_df, idcol = T))
  
  #
  go_df$classic[is.na(go_df$classic)] <- 0
  go_df <- go_df[go_df$classic <= 0.05, ]
  
  return(go_df)
  
}

Get_SS_PDI_Random_CoExp_TF <- function(tf){
  
  # Number of Sig GOs in PDI data
  Lgo.pdi <- length(subset(BP_PDIs, Mutant == tf)$GO.ID)
  
  # Number of Sig GOs in CoExp data
  Lgo.coe <- length(subset(BP_CoExp, Mutant == tf)$GO.ID)
  
  # Empty vector to storage results
  Random <- c()
  
  # Read GOs from random PDIs
  # PDI
  TF_R_pdi_gos <- R_pdi_gos[grepl(tf, R_pdi_gos)]
  TF_R_pdi_gos <- paste0("RandomNetsData/PDI_GOs/", TF_R_pdi_gos)
  # CoExp
  TF_R_coe_gos <- R_coe_gos[grepl(tf, R_coe_gos)]
  TF_R_coe_gos <- paste0("RandomNetsData/CoExp_GOs/", TF_R_coe_gos)
  
  # Maxinum numner of comparison equal to min lenght
  totalsamples <- min(c(length(TF_R_pdi_gos), length(TF_R_coe_gos)))
  print(totalsamples)
  
  # Read files
  rdb_pdi <- lapply(TF_R_pdi_gos[1:totalsamples], fread)
  rdb_coe <- lapply(TF_R_coe_gos[1:totalsamples], fread)
  
  for (s in 1:totalsamples){
    # Select top GOs with lenght eaul to actual observation
    # pdi
    gos.pdi <- as.character(as.data.frame(rdb_pdi[[s]])[,1])[1:Lgo.pdi]
    # coe
    gos.coe <- as.character(as.data.frame(rdb_coe[[s]])[,1])[1:Lgo.coe]
    
    # Sample pool of GOs significantly enrich observed in CoExp
    #r.gos <- sample(background_list, Lgo2_leng, replace = T)
    r.ss <- mgoSim(gos.pdi, gos.coe, semData=zmGO_BP, measure="Wang", combine="BMA")
    Random <- c(Random, r.ss)
    print(s)
  }

  df = tibble(ID=tf, r.ss=Random)
  #print(df)
  write.table(df, paste0("GO_SS_RandomBackground/GO_SS_random_PDI_CoExp.", tf, ".txt"), sep = '\t', quote = F, row.names = F)
  
  #
  
}


##############################################################################
##################         Read data input           #########################
##############################################################################

# PDI
PDI <- as_tibble(unique(fread("../Fig_PDI/PDI_NetworkFinal.10_11_2021.txt")[,c(2,3)]))
colnames(PDI)[1] <- "Source"
# teQTL
teQTL <- as_tibble(unique(fread("../Fig_transeQTL/teQTL_NetworkFinal.10_11_2021.txt")))
colnames(teQTL)[1] <- "Source"

# CoExp
CoExp <- as_tibble(unique(fread("../Fig_Coexpression/CoExp_NetworkFinal.10_11_2021.txt")))
colnames(CoExp)[2] <- "Source"

# GOs term annotations
background <- readMappings("synteny.ID_TopGO_V4_GRAMER.txt")
background_IDs <- as.character(unique(names(background)))
background_list <- unique(as.character(unlist(background)))

# list of GOs from random networks
R_pdi_gos <- list.files(path="RandomNetsData/PDI_GOs/", pattern = 'BP_rNet.*')
R_coe_gos <- list.files(path="RandomNetsData/CoExp_GOs/", pattern = 'BP_rNet.*')


################################################################
########   Count common TFs and number of associations  ######## 
################################################################

## Count common TFs and number of associations
#
PDI_counts <- as_tibble(as.data.frame(table(PDI$Source), stringsAsFactors = F))
colnames(PDI_counts) <- c("Source", "Targets")
#
CoExp_counts <- as_tibble(as.data.frame(table(unique(CoExp[,2:3])$Source), stringsAsFactors = F))
colnames(CoExp_counts) <- c("Source", "Targets")
#
teQTL_counts <- as_tibble(as.data.frame(table(teQTL$Source), stringsAsFactors = F))
colnames(teQTL_counts) <- c("Source", "Targets")


TFs_common <- tibble(Source=unique(c(PDI_counts$Source, CoExp_counts$Source, teQTL_counts$Source)))
#
TFs_common <- left_join(TFs_common, PDI_counts, by='Source')
TFs_common <- left_join(TFs_common, CoExp_counts, by='Source')
TFs_common <- left_join(TFs_common, teQTL_counts, by='Source')
colnames(TFs_common) <- c("Source", "PDI", "CoExp", "teQTL")

TFs_common[is.na(TFs_common)] <- 0

mask <- ((TFs_common$PDI>50)*1 + (TFs_common$CoExp>50)*1 ) >=2

Common_TFs_to_test <- TFs_common[mask,]$Source # based on PDI and PDI data

#################################################################
#############           Read GOs terms results      #############
#################################################################

#################################################################
##  Make file names of Sig. GOs to reads
#################################################################

# GO BP PDI
GO_bp_pdi <- list.files(pattern = "BP_PDI.GOs.*", path = 'BP_results/')
GO_bp_pdi <- paste0('BP_results/', GO_bp_pdi)

# GO BP CoExp
GO_bp_coexp <- list.files(pattern = "BP_CoExp.GOs.*", path = 'BP_results/')
GO_bp_coexp <- paste0('BP_results/', GO_bp_coexp)

# GO BP teQTL: 
GO_bp_teQTL <- list.files(pattern = "BP_teQTL.GOs.*", path = 'BP_results/')
GO_bp_teQTL <- paste0('BP_results/', GO_bp_teQTL)

####################################################

################################################################
## Read sig. GOs for corresponding TF 
################################################################
# BP
BP_PDIs <- ReadGOs(GO_bp_pdi)      # Read all GOs files for PDI data
BP_CoExp <- ReadGOs(GO_bp_coexp)   # Read all GOs files for CoExp data

# Reduce datasets to TFs with PDI and CoExp data
BP_PDIs <- subset(BP_PDIs, Mutant %in% Common_TFs_to_test)
BP_CoExp <- subset(BP_CoExp, Mutant %in% Common_TFs_to_test)
################################################################


########################################################
#####    Calculate semantic similarities        ########
########################################################

####################################################
## 1. Calculate semantic similarity background
ah = AnnotationHub()
zm <- query(ah, c("Zea mays"))
zm <- zm["AH85440"]

zm_hub <- ah[[zm$ah_id]] # Extract info for specific maize id

# Create Semantic similarity calculations: in this cases based on BP
zmGO_BP <- godata(zm_hub, ont ="BP")

## 2. Calculate semantic similarity measurement for PDI:GOs vs random CoExp:GOs
Lgenes <- length(Common_TFs_to_test)

w=40 # Size of range to test
print(".. Ready to start ..")

for (i in seq(0, Lgenes, w)){
  max=Lgenes
  Start=i+1
  end=i+w
  if (end<max){
    listtotest <- Common_TFs_to_test[Start:end]
    print(length(listtotest))
    
    mclapply(listtotest, Get_SS_PDI_Random_CoExp_TF, mc.cores=w)
    
  }
  else{
    r <- paste(Start, max, sep = "-")
    listtotest <- Common_TFs_to_test[Start:max]
    print(length(listtotest))
    mclapply(listtotest, Get_SS_PDI_Random_CoExp_TF, mc.cores=w)
  }
}


#Random_pdi_coExp_ss <- Get_SS_PDI_Random_CoExp(Common_TFs_to_test)
#df_random <- plyr::ldply(Random_pdi_coExp_ss_1, data.frame)
#colnames(df_random) <- c("ID", "r.ss")
#write.table(df_random, "SS_GOs_PDI_CoExp_random.txt", sep = '\t', quote = F, row.names = F)
