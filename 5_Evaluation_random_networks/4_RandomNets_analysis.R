suppressMessages(library(parallel))
suppressMessages(library(data.table))
suppressMessages(library(scales))
suppressMessages(library(tidyverse))
suppressMessages(library(ggVennDiagram))
suppressMessages(library(purrr))
suppressMessages(library(gplots))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(ggpubr))
suppressMessages(library(viridis))
suppressMessages(library(patchwork))
suppressMessages(library(rrvgo))
suppressMessages(library(reshape2))
suppressMessages(library(GOSemSim))
suppressMessages(library(enrichplot))
suppressMessages(library(GeneOverlap))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(circlize))
library(ggridges)
library(patchwork)


#library(org.Zmays.eg.db)
# To install org.Zmays.eg.db: ("/maindisk/fabio/Projects/Genomes/Maize/B73v4/FunAnntation_GRAMER/org.Zmays.eg.db")

rm(list=ls())

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

ReplaceName <- function(ids){
  # for (i in 1:nrow(Top45)){
  #   ids <- gsub(Top45$V1[i], Top45$V2[i], ids)
  # }
  
  for (i in 1:nrow(TFdic)){
    ids <- gsub(TFdic$V2[i], TFdic$V1[i], ids)
  }
  return(ids)
}

CountGO_nbase <- function(netid_TF){
  netid <- chop(netid_TF, '[_]', 1)
  tf    <- chop(netid_TF, '[_]', 2)
  
  # GO mapped to parent
  file1 <- paste0('BP_RandomGOs_NetBase_Parents/Net.', netid,".RedRandom_NetBased_",tf,'.txt')
  
  # Original GO enrichment test to get Pvalue
  file2 <- paste0('BP_results_NetBase_random/Random.FullMR.', netid,"_",tf,'_GOs.txt')
  
  # BP_results_NetBase_random/Random.FullMR.1000_Zm00001d005578_GOs.txt
  
  # BP_RandomGOs_NetBase_Parents/Net.1000.RedRandom_NetBased_Zm00001d005578.txt
  df1  <- fread(file1) 
  df2  <- fread(file2)[,c('GO.ID','classic')] 
  df2[,'FDR'] <- p.adjust(df2$classic, method = 'fdr')
  df1 <- left_join(df1, df2, 'GO.ID')
  df1 <- as_tibble(df1[,-c(4)])
  
  return(df1)
  
}

CountGO_cfunct <- function(netid_TF){
  tf    <- chop(netid_TF, '[_]', 2)
  netid <- chop(chop(netid_TF, '[_]', 1),'[.]', 1)
  typeNet <- gsub(':','_', chop(chop(netid_TF, '[_]', 1),'[.]', 2))
  # GO mapped to parent
  # BP_RandomGOs_CommFunct_Parents/RedRandom_CommFunct_Net_1000_Zm00001d005892_GRN_eGRN.txt
  file1 <- paste0('BP_RandomGOs_CommFunct_Parents/RedRandom_CommFunct_Net_', netid,"_",tf,"_",typeNet,'.txt')
  
  # Original GO enrichment test to get Pvalue
  #BP_GSS_CommonFunction_random/GSS.1000.Zm00001d005578.GRN_CEN.txt
  file2 <- paste0('BP_GSS_CommonFunction_random/GSS.', netid,".",tf,'.',typeNet,'.txt')
  
  # Read files after reduction
  df1  <- fread(file1)[,2]
  
  # Read files before reduction
  df2  <- fread(file2)[,c(1,2,6,9)]
  colnames(df2) <- c('GO1', "GO2", "FDR1", "FDR2")
  df2a <- unique(df2[,c(1,3)])
  df2b <- unique(df2[,c(2,4)])
  colnames(df2a) <- c('GO.ID', "FDR")
  colnames(df2b) <- c('GO.ID', "FDR")
  df2 <- unique(rbind(df2a, df2b))
  
  #
  df1 <- left_join(df1, df2, "GO.ID")
  
  
  # group by GOs and keep min FDR
  df1 %>% group_by(GO.ID) %>%
    summarise(FDR = min(FDR)) -> df1
  
  df1[,'.id'] <- netid
  df1[,'TF'] <- tf
  df1 <- df1[,c("TF", ".id", "GO.ID","FDR")]
  
  return(df1)
  
}

###############################################



###############################################
###            Input files                  ###
###############################################

## Syntenic genes 
Syntenic <- as_tibble(read.table("../Data/Annotations/Zm.v4.synteny.genes.txt", h=T, stringsAsFactors = F, sep = '\t'))
Syntenic <- subset(Syntenic, gene_synteny=='syntenic')$gene_id

# TF names
TFdic <- as_tibble(read.table("../Data/Annotations/TF_Id_Name.txt", h=F, stringsAsFactors = F))

# TFs of interest
AllTFs <- unique(fread("TFs_2_test_RNets.txt", h=F)$V1)

######
# Total random nets tested by method
######
## Network-based
NetBased <- list.files(path = 'BP_results_NetBase_random/', pattern = "^Random.FullMR.*")
NetBased <- tibble(TF=chop(NetBased, '[_]', 2), 
                   Rnet=gsub("Random.FullMR.", "", chop(NetBased, '[_]', 1))) %>%
  dplyr::arrange(TF)
# NetBased
## Common function
CommFunct <- list.files(path = 'BP_GSS_CommonFunction_random/', pattern = "^GSS.*")
CommFunct <- tibble(TF=chop(CommFunct, '[.]', 3), 
                    Rnet=gsub("_", ":", paste0(chop(CommFunct, '[.]', 2),'.', 
                                               chop(CommFunct, '[.]', 4)))) %>%
  dplyr::arrange(TF)

#CommFunct
## Common targets
#BP_results_CommonTargets_random/Random.CommonTarget.1000_Zm00001d005578_GOs.txt
CommTargt <- list.files(path = 'BP_results_CommonTargets_random/', pattern = "^Random.CommonTarget.*")

CommTargt <- tibble(TF =chop(CommTargt, '[_]', 2), 
                    Rnet=gsub("Random.CommonTarget.", "", chop(CommTargt, '[_]', 1))) %>%
  dplyr::arrange(TF)

# Pre-calculate semantic similarity
Zm.GOSemSim.BP <- GOSemSim::godata(org.Zmays.eg.db, keytype="GENENAME", ont = 'BP')

# Save a single object to a file
#saveRDS(Zm.GOSemSim.BP, "mtcars.rds")

# Pre-calculated semantic similarity for from org.Zmays.eg.db for BP term
#Zm.GOSemSim.BP <- readRDS("Zm.GOSemSim.BP.rds")
######

######
# Total random nets after FDR filtering and GO reductions
######

## Network-based
Red_NetBased <- list.files(path = 'BP_RandomGOs_NetBase_Parents/', pattern = "^Net.*")
Red_NetBased <- tibble(TF=gsub(".txt", "", chop(Red_NetBased, '[_]', 3)), 
                       Rnet=chop(Red_NetBased, '[.]', 2)) %>%
  dplyr::arrange(TF)
# NetBased

## Common function
Red_CommFunct <- list.files(path = 'BP_RandomGOs_CommFunct_Parents/', pattern = "^RedRandom*")
Red_CommFunct <- tibble(TF=gsub(".txt", "", chop(Red_CommFunct, '[_]', 5)), 
                        Rnet=gsub(".txt", "", paste0(chop(Red_CommFunct, '[_]', 4),'.', 
                                                     chop(Red_CommFunct, '[_]', 6),':',
                                                     chop(Red_CommFunct, '[_]', 7)))) %>%
  dplyr::arrange(TF)
#Red_CommFunct
## Common targets
#
Red_CommTargt <- list.files(path = 'BP_RandomGOs_CommTarg_Parents/', pattern = "^Net.*")
length(Red_CommTargt)

#Red_CommTargt <- Red_CommTargt[grepl('Commtarg', Red_CommTargt)]
Red_CommTargt <- tibble(TF =gsub(".txt", "",  chop(Red_CommTargt, '[_]', 3)), 
                        Rnet=gsub(".txt", "", chop(Red_CommTargt, '[.]', 2))) %>%
  dplyr::arrange(TF)



######

######
# Real GO term results
######
TrueGOs <- fread("../Fig_MethodsComparison/ReduceGOterms_All_methods.txt")
######

###############################################

###############################################
###    Comparison tested vs keep  nets      ###
###############################################

# Add class to count overlapping
Red_NetBased[,'tested']  <- 1
Red_CommFunct[,'tested'] <- 1
Red_CommTargt[,'tested'] <- 1

# Add 'tested' class
NetBased  <- left_join(NetBased, Red_NetBased, by=c('TF', "Rnet"))
CommFunct <- left_join(CommFunct, Red_CommFunct, by=c('TF', "Rnet"))
CommTargt <- left_join(CommTargt, Red_CommTargt, by=c('TF', "Rnet"))

NetBased$tested[is.na(NetBased$tested)] <- 0
CommFunct$tested[is.na(CommFunct$tested)] <- 0
CommTargt$tested[is.na(CommTargt$tested)] <- 0

#
Total_Nets <- rbind(NetBased %>%
                      dplyr::group_by(TF) %>%
                      dplyr::summarise(Total=n(), Tested=sum(tested), Method='Network-based'),
                    CommFunct %>%
                      dplyr::group_by(TF) %>%
                      dplyr::summarise(Total=n(), Tested=sum(tested), Method='Comm.Function'),
                    CommTargt %>%
                      dplyr::group_by(TF) %>%
                      dplyr::summarise(Total=n(), Tested=sum(tested), Method='Comm.Target'))


Total_Nets$Tested <- Total_Nets$Tested/Total_Nets$Total



###############################################


###############################################
###           Counting total GOs            ###
###############################################

# Net-based: id to reads GO files
NetBased_ids  <- paste0(Red_NetBased$Rnet, '_', Red_NetBased$TF)

# Comm-Function: id to reads GO files
ComFunct_ids <- paste0(Red_CommFunct$Rnet, '_', Red_CommFunct$TF)

#
ComTarg_ids    <- paste0(Red_CommTargt$Rnet, '_', Red_CommTargt$TF)
RawComTarg_ids <- paste0(CommTargt$Rnet, '_', CommTargt$TF)
length(ComTarg_ids)
length(RawComTarg_ids)

ComTarg_ids <- ComTarg_ids[ComTarg_ids %in% RawComTarg_ids]
length(ComTarg_ids)

CountGO_ctarg <- function(netid_TF){
  netid <- chop(netid_TF, '[_]', 1)
  tf    <- chop(netid_TF, '[_]', 2)
  
  # GO mapped to parent
  #BP_RandomGOs_CommTarg_Parents/Net.1000.RedRandom_Commtarg_Zm00001d005892.txt
  file1 <- paste0('BP_RandomGOs_CommTarg_Parents/Net.', netid,".RedRandom_Commtarg_",tf,'.txt')
  
  # Original GO enrichment test to get Pvalue
  # BP_results_CommonTargets_random/Random.CommonTarget.1000_Zm00001d005578_GOs.txt
  file2 <- paste0('BP_results_CommonTargets_random/Random.CommonTarget.', netid,"_",tf,'_GOs.txt')
  
  df1  <- fread(file1)[,c('TF','.id','GO.ID')] 
  df2  <- fread(file2)[,c('GO.ID','classic')] 
  df2[,'FDR'] <- p.adjust(df2$classic, method = 'fdr')
  df1 <- left_join(df1, df2, 'GO.ID')
  df1 <- as_tibble(df1[,-c(4)])
  # input FDR from the mean of FDR observed
  # into parents not presented in initial top 1000 GO tested
  df1$FDR[(is.na(df1$FDR))] <- mean(df1$FDR[!(is.na(df1$FDR))])
  
  return(df1)
  
}

# Counts GOs: Net. based
DFGO_nbase <- lapply(NetBased_ids, CountGO_nbase)
DFGO_nbase <- rbindlist(DFGO_nbase, idcol = F)

write.table(DFGO_ctarg, 
            'Obs_vs_Random_GSS_results/Input_Random_nbase.txt', 
            sep='\t', row.names = F, quote = F)

# Counts GOs: Com. funct
DFGO_cfunct <- lapply(ComFunct_ids, CountGO_cfunct)
DFGO_cfunct <- rbindlist(DFGO_cfunct, idcol = F)

write.table(DFGO_cfunct, 
            'Obs_vs_Random_GSS_results/Input_Random_cFunction.txt', 
            sep='\t', row.names = F, quote = F)

# Counts GOs: Com. targ
DFGO_ctarg <- lapply(ComTarg_ids, CountGO_ctarg)
DFGO_ctarg <- rbindlist(DFGO_ctarg, idcol = F)

write.table(DFGO_ctarg, 
            'Obs_vs_Random_GSS_results/Input_Random_cTarget.txt', 
            sep='\t', row.names = F, quote = F)

# Count GOs in R.net per TF
## net-based
Table_DFGO_nbase <- as.data.table(table(DFGO_nbase[,c(1,2)])) %>%
  dplyr::filter(N>0) %>%
  dplyr::mutate(Method='Network-based')

## com-function
Table_DFGO_cfunct <- as.data.table(table(DFGO_cfunct[,c(1,2)])) %>%
  dplyr::filter(N>0) %>%
  dplyr::mutate(Method='Comm.Function')

## c-target
Table_DFGO_ctarg <- as.data.table(table(DFGO_ctarg[,c(1,2)])) %>%
  dplyr::filter(N>0) %>%
  dplyr::mutate(Method='Comm.Target')

## average GOs per TF in all Random nets
Table_Counts <- rbind(Table_DFGO_nbase, Table_DFGO_cfunct, Table_DFGO_ctarg) %>%
  dplyr::group_by(TF, Method) %>%
  dplyr::summarise(MeanN=mean(N))

## 


###############################################

###############################################
###              Summary FDR                ###
###############################################
 
Table_FDR_ctarg <- DFGO_ctarg %>%
  dplyr::group_by(TF, .id) %>%
  dplyr::summarise(MeanFDR=mean(FDR)) %>% 
  dplyr::mutate(Method='Comm.Target')

Table_FDR_cfunct <- DFGO_cfunct %>%
  dplyr::group_by(TF, .id) %>%
  dplyr::summarise(MeanFDR=mean(FDR)) %>% 
  dplyr::mutate(Method='Comm.Function')
Table_FDR_cfunct$.id <- as.numeric(Table_FDR_cfunct$.id)

Table_FDR_nbase <- DFGO_nbase %>%
  dplyr::group_by(TF, .id) %>%
  dplyr::summarise(MeanFDR=mean(FDR)) %>% 
  dplyr::mutate(Method='Network-based')

# summary by TFs, Method, and Network
Table_FDR <- rbind(Table_FDR_ctarg, Table_FDR_cfunct, Table_FDR_nbase) %>%
  dplyr::mutate(log10FDR=-log10(MeanFDR))

# summary by TFs and Method
Table_FDR_all <- Table_FDR %>%
  dplyr::group_by(TF, Method) %>%
  dplyr::summarise(Meanlog10FDR=mean(log10FDR))
###############################################

###############################################
###           GSS observed vs random        ###
###############################################

# results from:
## 4_1_GSS_Obs_vs_random_nbase.R
## 4_1_GSS_Obs_vs_random_cFunc.R
## 4_1_GSS_Obs_vs_random_cTarg.R

# net-base
GSSr_nbase <- list.files(path = '4_1_Obs_vs_Random_GSS_results/GSSrDB/', pattern = "^GSSr_nbase_*")
length(GSSr_nbase)
# Com. Funct
GSSr_cFunc <- list.files(path = '4_1_Obs_vs_Random_GSS_results/GSSrDB/', pattern = "^GSSr_cFunctionn_*")
length(GSSr_cFunc)
# Com. Target
GSSr_cTarg <- list.files(path = '4_1_Obs_vs_Random_GSS_results/GSSrDB/', pattern = "^GSSr_cTarget_*")
length(GSSr_cTarg)

#
GSSr_nbase <- lapply(GSSr_nbase, function(x) fread(paste0("4_1_Obs_vs_Random_GSS_results/GSSrDB/", x)))
GSSr_cFunc <- lapply(GSSr_cFunc, function(x) fread(paste0("4_1_Obs_vs_Random_GSS_results/GSSrDB/", x)))
GSSr_cTarg <- lapply(GSSr_cTarg, function(x) fread(paste0("4_1_Obs_vs_Random_GSS_results/GSSrDB/", x)))

# 
GSSr_nbase <- rbindlist(GSSr_nbase, idcol = F)
GSSr_cFunc <- rbindlist(GSSr_cFunc, idcol = F)
GSSr_cTarg <- rbindlist(GSSr_cTarg, idcol = F)
GSSr <- rbind(GSSr_nbase, GSSr_cFunc, GSSr_cTarg)
# add names
GSSr <- left_join(GSSr, TFdic, by=c("TF"="V2"))
colnames(GSSr)[5] <- 'Name'


# summary by TFs and Method
GSSr_all <- GSSr %>%
  dplyr::group_by(TF, Method) %>%
  dplyr::summarise(MeangSS=mean(GSSr, na.rm = T))


###############################################

##############################################################
###   Z-score using random values all null: network based  ###
##############################################################

#DFGO_cfunct
#DFGO_ctarg

DFGO_nbase

# number og GOs: Defined background models
n_GOs_background <- split(Table_DFGO_nbase$N, Table_DFGO_nbase$TF)

# Average FDR by network: Defined background models
## fdr_GOs_background <- split(Table_FDR_nbase$MeanFDR, Table_DFGO_nbase$TF)

# Get Zscore from observed values               
Zscore_nGOs <- function(tf){
  obs <- subset(TrueGOs, TF==tf)
  obs <- length(obs$parent)
  val <- (obs - mean(n_GOs_background[[tf]])) / sd(n_GOs_background[[tf]])
  return(val)
}

# # Get Zscore from observed values               
# Zscore_fdrGOs <- function(tf){
#   obs <- subset(TrueGOs, TF==tf)
#   obs <- length(obs$parent)
#   val <- (obs - mean(n_GOs_background[[tf]])) / sd(n_GOs_background[[tf]])
#   return(val)
# }

# Defined out DF
DF_Zscore_nGOs <- tibble(TF=unique(Table_DFGO_nbase$TF))

# Z score from random background
DF_Zscore_nGOs[,'Z_nGOs'] <- unlist(lapply(DF_Zscore_nGOs$TF, Zscore_nGOs)) 

# P values from Z score from random background
#DF_Zscore_nGOs[,"Pval_neg"] <- unlist(lapply(DF_Zscore_nGOs$Z_nGOs, function (x) pnorm(x, lower.tail=TRUE)))
DF_Zscore_nGOs[,"Pval"] <- unlist(lapply(DF_Zscore_nGOs$Z_nGOs, function (x) pnorm(x, lower.tail=FALSE)))

# Add Zvelues to DF to make plot
Table_DFGO_nbase %>%
  dplyr::group_by(TF) %>%
  dplyr::mutate( Z = (N-mean(N))/sd(N)) -> Table_DFGO_nbase
                

# DF_Zscore_nGOs[,'Name'] <- paste0(ReplaceName(DF_Zscore_nGOs$TF),
#                                   ' (',
#                                   formatC(DF_Zscore_nGOs$Pval, format = "e", digits = 1),
#                                   ")")

###############################################

###############################################

GONAMES_DB <- unique(rbind(GO_CommTarg_Red[,c('go', "term")] %>% dplyr::rename(GO.ID = go, GO.name=term),
                           GO_CommTarg_Red[,c('parent', "parentTerm")] %>% dplyr::rename(GO.ID = parent, GO.name=parentTerm),
                           GO_CommFunt_Red[,c('go', "term")] %>% dplyr::rename(GO.ID = go, GO.name=term),
                           GO_CommFunt_Red[,c('parent', "parentTerm")] %>% dplyr::rename(GO.ID = parent, GO.name=parentTerm),
                           GO_Network_Red[,c('go', "term")] %>% dplyr::rename(GO.ID = go, GO.name=term),
                           GO_Network_Red[,c('parent', "parentTerm")] %>% dplyr::rename(GO.ID = parent, GO.name=parentTerm)))
######

###############################################

###############################################
###                Plots                    ###
###############################################

#####
# Plot2f Total nets tested
#####
Total_Nets$Method <- factor(Total_Nets$Method, 
                            levels = c("Comm.Target",
                                       "Comm.Function",
                                       "Network-based"))

ggplot(Total_Nets, aes(x=Tested, y=Method, fill=Method)) +
  geom_jitter(alpha=0.4, size=0.7, height = 0.1) +
  geom_violin() +
  scale_fill_viridis(discrete = T, option = "D", alpha = 0.5, direction = -1) + 
  theme_pubclean() +
  scale_x_continuous(expand = c(0,0), limits = c(0, 1.4), breaks = c(0, 0.5, 1)) + 
  scale_y_discrete(expand = c(0,0)) + 
  xlab("Fraction of Random networks\nwith at least a GO term") +
  ylab("") +
  theme(strip.text.x = element_text(size = 10), 
        legend.position = 'none',
        #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text=element_text(size=10), 
        text = element_text(size=10, family="Times")) +
  stat_compare_means(method = 'wilcox.test', size=2, paired = F, label="p.signif",
                      comparisons=list(c("Comm.Target", "Comm.Function"),
                                       c('Comm.Function', 'Network-based'),
                                       c('Comm.Target', 'Network-based'))) -> Plot_2f

Plot_2f

Total_Nets %>% 
  dplyr::group_by(Method) %>%
  dplyr::summarise(Total=mean(Total), Tested=mean(Tested))

#####

#####
# Plot2g Counts nets tested
#####
Table_Counts$Method <- factor(Table_Counts$Method, 
                            levels = c("Comm.Target",
                                       "Comm.Function",
                                       "Network-based"))

ggplot(Table_Counts, aes(x=MeanN, y=Method, fill = Method)) +
  geom_boxplot(alpha=0.5, notch = T) +
  geom_jitter(alpha=0.4, height = 0.1, size=1) +
  scale_fill_viridis(discrete = T, option = "D", alpha = 0.5, direction = -1) + 
  theme_pubclean() +
  scale_x_continuous(expand = c(0,0), limits = c(0, 300), breaks = c(0,100,200,300)) + 
  xlab("Average GOs by\nrandom net. per TF") +
  ylab("") +
  theme(strip.text.x = element_text(size = 10), 
        legend.position = 'none',
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text=element_text(size=10), 
        text = element_text(size=10, family="Times")) +
  stat_compare_means(method = 'wilcox.test', size=2, paired = F, 
                     label="p.signif", label.y = c(200, 230, 260),
                      comparisons=list(c("Comm.Target", "Comm.Function"),
                                       c('Comm.Function', 'Network-based'),
                                       c('Comm.Target', 'Network-based'))) -> Plot_2g
Plot_2g


Table_Counts %>% 
  dplyr::group_by(Method) %>%
  dplyr::summarise(N=mean(MeanN))

#####

#####
# Plot2h -log10(FDR) by net
#####

Table_FDR_all$Method <- factor(Table_FDR_all$Method, levels = c("Comm.Target", "Comm.Function", "Network-based"))

ggplot(Table_FDR_all, aes(y=Method, x=Meanlog10FDR, fill = Method)) +
  geom_boxplot(alpha=0.5, notch = T) +
  geom_jitter(alpha=0.4, height = 0.1, size=1) +
  scale_fill_viridis(discrete = T, option = "D", alpha = 0.5, direction = -1) + 
  theme_pubclean() +
  scale_x_continuous(expand = c(0,0), limits = c(1, 3)) + 
  ylab("") +
  xlab(bquote("Average ("~-Log[10] ~ "FDR)")) +
  theme(strip.text.x = element_text(size = 10), 
        legend.position = 'none',
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text=element_text(size=10), 
        text = element_text(size=10, family="Times")) +
  stat_compare_means(method = 'wilcox.test', size=2, paired = F, label="p.signif",
                      comparisons=list(c("Comm.Target", "Comm.Function"),
                                       c('Comm.Function', 'Network-based'),
                                       c('Comm.Target', 'Network-based'))) -> Plot_2h

#####

#####
# Plot3i GSS by Method
#####

ggplot(GSSr_all, aes(y=Method, x=MeangSS, fill = Method)) +
  geom_boxplot(alpha=0.5, notch = T) +
  geom_jitter(alpha=0.4, height = 0.1, size=1) +
  scale_fill_viridis(discrete = T, option = "D", alpha = 0.5, direction = -1) + 
  theme_pubclean() +
  scale_x_continuous(expand = c(0,0), limits = c(0, 1)) + 
  ylab("") +
  xlab("Average (GSS)") +
  theme(strip.text.x = element_text(size = 10), 
        legend.position = 'none',
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text=element_text(size=10), 
        text = element_text(size=10, family="Times")) +
  stat_compare_means(method = 'wilcox.test', size=2, paired = F, label="p.signif",
                     comparisons=list(c("Com.Target", "Com.Function"),
                                      c('Com.Function', 'Network-base'),
                                      c('Com.Target', 'Network-base'))) -> Plot_2i


Plot_2fghi <- Plot_2f|Plot_2g|Plot_2h|Plot_2i

pdf("Plots/Plot_2fghi.pdf", width=8, height=2.5)
print(Plot_2fghi)
dev.off()
#####

#####
# Plot S10a: -log10(FDR) by net
#####

Table_FDR[,'Name'] <- ReplaceName(Table_FDR$TF)

Table_FDR$Method <- factor(Table_FDR$Method, levels = c("Comm.Target", "Comm.Function", "Network-based"))

ggplot(Table_FDR, aes(x = log10FDR, y = Name, fill = Method)) +
  geom_density_ridges(bandwidth=0.05) +
  theme_ridges() + 
  scale_fill_viridis(discrete = T, option = "D", alpha = 0.6, direction = -1) + 
  scale_x_continuous(expand = c(0,0), limits = c(0.8, 3)) + 
  scale_y_discrete(expand = c(0,0)) + 
  xlab(bquote(-Log[10] ~ "FDR")) +
  ylab("TF") +
  theme(strip.text.x = element_text(size = 10), 
        legend.position = 'none',
        axis.text=element_text(size=10), 
        text = element_text(size=10, family="Times"),
        legend.key.size = unit(0.2, 'cm')) -> Plot_S10a

#####

#####
# Plot s10b: GSSr by TF and methods
#####

GSSr$Method <- factor(GSSr$Method, levels = c("Com.Target", "Com.Function", "Network-base"))
unique(GSSr$Method)

ggplot(GSSr, aes(x = GSSr, y = Name, fill = Method)) +
  geom_density_ridges(bandwidth=0.05) +
  theme_ridges() + 
  scale_fill_viridis(discrete = T, option = "D", alpha = 0.6, direction = -1) + 
  scale_x_continuous(expand = c(0,0)) +  # limits = c(0, 1)
  scale_y_discrete(expand = c(0,0)) + 
  xlab("GO semantic similarity (GSS)") +
  ylab("TF") +
  theme(strip.text.x = element_text(size = 10), 
        legend.position = 'bottom',
        #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text=element_text(size=10), 
        text = element_text(size=10, family="Times"),
        legend.key.size = unit(0.2, 'cm')) -> Plot_S10b

##
Plot_S10 <- (Plot_S10a + Plot_S10b) + plot_layout(guides = "collect") & theme(legend.position = "bottom")

pdf("Plots/Plot_S10.pdf", width=8, height=8)
print(Plot_S10)
dev.off()

#####

#####
# 
#####
DF_Zscore_nGOs[,'Name'] <- ReplaceName(DF_Zscore_nGOs$TF)

DF_Zscore_nGOs$Z_nGOs[DF_Zscore_nGOs$Z_nGOs >= 10] <- 10

ggplot(Table_DFGO_nbase, aes(x=Z)) + 
  geom_density(alpha=.1, colour="#C0C0C0", fill="#DCDCDC") +
  geom_segment(data=DF_Zscore_nGOs,
               aes(y=0, yend=1.5, x=Z_nGOs, xend=Z_nGOs),
               linetype="dashed", linewidth=1.2, color = "#FFA500") +
  scale_x_continuous(expand = c(0,0), limits = c(-1,11)) + 
  geom_text_repel(data=DF_Zscore_nGOs,
                  aes(x=Z_nGOs, y=2, label=formatC(Pval, format = "e", digits = 1)),
                  direction    = "y",
                  xlim = c(NA, 9),
                  ylim = c(1.5, NA),
                  vjust = 1, 
                  segment.size = 0.3,
                  max.iter = 1e4, 
                  box.padding = 0.1,
                  max.time = 1,
                  max.overlaps = Inf,
                  force_pull = 1,
                  size=3) + 
  facet_wrap( .~ Name, ncol = 4) +
  theme_pubclean() +
  ylab('Density') + 
  xlab('Z-score (GOs in random networks)') +
  #labs(subtitle='Tau null distribution for n=14') + 
  theme(
    #strip.text.y = element_text(size = 5, angle = 0), 
    #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.text=element_text(size=10), 
    #legend.position = 'bottom',
    text = element_text(size=10, family="Times")) -> Plot_Z_nGOs_nbase

Plot_Z_nGOs_nbase

pdf("Plots/Plot_S11.pdf", width=7, height=10)
print(Plot_Z_nGOs_nbase)
dev.off()


#####
