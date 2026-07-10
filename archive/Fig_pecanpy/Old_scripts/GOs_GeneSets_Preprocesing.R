library(lsa)
library(parallel)
suppressMessages(library(Rgraphviz))
library(RColorBrewer)
library(topGO)
library(GeneOverlap)
library(tidyverse)
library(ggrepel)
library(ggpubr)
library(data.table)
library(viridis)
library(circlize)
library(ggVennDiagram)
library(scales)
library(purrr)
library(gplots)
library(ggplot2)
library(rrvgo)
library(ggdark)
library(org.Zmays.eg.db)
library(glmnet)

##################################################
##########          Functions        #############
##################################################

WideAndScale <- function(long_df){
  df.wide <- as.data.frame(pivot_wider(long_df, names_from = TF, values_from = Cos))
  row.names(df.wide) <- df.wide$GO
  df.wide <- df.wide[,-c(1)]
  
  
  df.wide1 <- apply(df.wide, 2, Zscale)
  df.wide1[df.wide1 < 0] <- 0
  
  df.wide2 <- t(apply(t(df.wide), 2, Zscale))
  df.wide2[df.wide2 < 0] <- 0
  
  df.wide <- sqrt(df.wide1^2 + df.wide2^2)
  
  # Filter to TFs  with at least X beta values largers than zero
  #df.wide <- df.wide[,(colSums(df.wide) > 11)]
  return(df.wide)
}

Zscale <- function(vals){
  mval <- mean(vals)
  sval <- sd(vals)
  
  z <- unlist(lapply(vals, function(x) (x - mval)/sval))
  z <- unlist(lapply(z, function(x) max(0, x)))
  #
  z[is.na(z)] <- 0
  
  return(z)
}

chop=function(myStr,mySep,myField){
  
  choppedString=sapply(strsplit(myStr,mySep),"[",myField)
  if(length(myField)>1){
    choppedString=apply(choppedString,2,function(x){paste0(x[!is.na(x)],collapse=mySep)})
  }
  return(choppedString)
}

lmg_go <- function(go){  
  
  # Define Genes associated with GO target
  GeneTARGET <- unique(GOzm.DF_4$GeneId[GOzm.DF_4$GO.ID == go])
  
  # Define labels
  GeneClass <- rep(1, nrow(GOs_Pecan))
  
  # Define Positive labels 
  GeneClass[row.names(GOs_Pecan) %in% GeneTARGET] <- 2
  
  # 
  cvbeta <- cv.glmnet(GOs_Pecan, GeneClass,  family = "binomial", type.measure = "deviance")
  # 
  out <- as.matrix(coef(cvbeta, s='lambda.min'))
  out <- data.table(Embedding=row.names(out)[-c(1)], B=out[2:nrow(out),])
  
  colnames(out)[2] <- go
  
  return(out)
}

ReplaceName <- function(ids){
  
  for (i in 1:nrow(TFdic)){
    ids <- gsub(TFdic$TF.v4[i], TFdic$TF.Name[i], ids)
  }
  return(ids)
}

##################################################

##################################################
########           data Input        #############
##################################################

# TF names
TFdic <- as_tibble(read.table("Data/Annotations/TF_Id_Name.txt", h=T, stringsAsFactors = F))

# TF and CoReg
TF_CoR <- as_tibble(read.table("Data/Annotations/TF_CoR_Mazie.txt", h=T, stringsAsFactors = F)) 

# PDI
GRN <- as_tibble(unique(fread("../Fig_PDI/Only_PDI_NetworkFinal.10_14_2022.txt")[,c(2,3)]))
colnames(GRN)[1] <- "Source"

eGRN <- as_tibble(unique(fread("../Fig_PDI/CisE_PDI_NetworkFinal.10_14_2022.txt")[,c(2,3)]))
colnames(eGRN)[1] <- "Source"

# add eGRN to GRN to count it twice: two evidences 
GRN <- rbind(GRN, eGRN)

# CoExp
CEN <- as_tibble(unique(fread("../Fig_Coexpression/CoExp_NetworkFinal.10_11_2021.txt")))
colnames(CEN)[2] <- "Source"
CEN <- unique(CEN[,2:3])

# teQTL
GAN <- as_tibble(unique(fread("../Fig_transeQTL/teQTL_NetworkFinal.10_11_2021.txt")))
colnames(GAN)[1] <- "Source"

# teQTL associated with TFs
GAN_tf <- subset(GAN, Source %in% unique(c(TF_CoR$GeneID, PDI$Source, CoExp$Source))) 

# Selection of TF-target with at least two evidences
Net_HighConfident <- c(paste(GRN$Source, GRN$Target, sep = "_"),
                       paste(CEN$Source, CEN$Target, sep = "_"),
                       paste(GAN$Source, GAN$Target, sep = "_"),
                       paste(eGRN$Source, eGRN$Target, sep = "_"))

Net_HighConfident <- data.table(table(Net_HighConfident)) 


Net_HighConfident %>%
  filter(N >= 2) -> Net_HighConfident

Net_HighConfident[,"TF"] <- chop(Net_HighConfident$Net_HighConfident, '[_]', 1)
Net_HighConfident[,"Targ"] <- chop(Net_HighConfident$Net_HighConfident, '[_]', 2)

# Count Targets by TF
Net_HighConfident %>%
  select(TF) %>%
  table() %>%
  as_tibble() %>%
  filter(n > 20) -> TFs_HighConfident

Net_HighConfident
TFs_HighConfident
##################################################

##############################################################
#########      Process GO terms: Specificity     #############
##############################################################

## Syntenic genes 
Syntenic <- as_tibble(read.table("Data/Annotations/Zm.v4.synteny.genes.txt", h=T, stringsAsFactors = F, sep = '\t'))
Syntenic <- subset(Syntenic, gene_synteny=='syntenic')$gene_id

# Maize GOs term annotations
GOzm <- readMappings("synteny.ID_TopGO_V4_GRAMER.txt")
GOzm_IDs <- as.character(unique(names(GOzm)))

# Unlist as reshape as DF
GOzm.DF <- rbindlist(lapply(GOzm, function(x) data.table(unlist(x))), idcol = T)
colnames(GOzm.DF) <- c("GeneId", 'GO.ID')
GOzm.DF <- unique(GOzm.DF)


# Count freq. as formatted it as DF
GOzmCount_1 <- data.table(table(GOzm.DF$GO.ID)) %>%
  arrange(N) 
colnames(GOzmCount_1)[1] <- "GO.ID"

table(GOzmCount_1$N)


rbindlist(lapply(GOzmCount_1[GOzmCount_1$N >= 500,]$GO.ID, GetGODescription), fill = T)

# GO vectors all GOs in maize
GOzm.vector <- unique(GOzm.DF$GO.ID)
GOzm.vector <- GOzm.vector[!(GOzm.vector %in% GOzmCount_1[GOzmCount_1$N >= 400,]$GO.ID)]
length(GOzm.vector)

######################## 
# Filter 1: Redundancy
########################

# Pre-calculate sematic similarity
Zm.GOSemSim.BP <- GOSemSim::godata(org.Zmays.eg.db, keytype="GENENAME", ont = 'BP')

# Semantic Similarity matrix
simMatrix <- calculateSimMatrix(GOzm.vector, 
                                orgdb=org.Zmays.eg.db, 
                                ont="BP", 
                                semdata=Zm.GOSemSim.BP,
                                method="Wang")

# Reduce term
reducedTerms <- reduceSimMatrix(simMatrix,
                                keytype="GENENAME",
                                threshold=0.4,
                                orgdb=org.Zmays.eg.db)


# Parents
GOzmCount_2 <- as.data.table(unique(reducedTerms[,c("parent", "parentTerm")]))

GOzm.DF_2 <- GOzm.DF[GOzm.DF$GO.ID %in% GOzmCount_2$parent,]

GOzmCount_2 <- data.table(table(GOzm.DF_2$GO.ID)) %>%
  arrange(N) 
colnames(GOzmCount_2)[1] <- "GO.ID"
GOzmCount_2
#Descarted <- rbindlist(lapply(GOzmCount_2[GOzmCount_2$N < 10,]$GO.ID, GetGODescription), fill = T)

# Final set OGs that pass 
#GOzmCount_2 <- subset(GOzmCount_2,  N >= 20)


######################## 

######################## 
# Filter 2: Redundancy
########################

GOzm.DF_2 <- subset(GOzm.DF, GO.ID %in% GOzmCount_2$GO.ID)
TotalGenesGenome <- length(Syntenic)


# define list to test redundancy
GOdf2_list <- split(GOzm.DF_2$GeneId, GOzm.DF_2$GO.ID)

namesGOdf2_list <- names(GOdf2_list)

CompareGeneSets <- function(goid){
  
  ## Compare list of predicted targets vs annotated genes in query.vector 
  #
  print(". Pre-newGOM .")
  go.obj <- newGOM(GOdf2_list[goid], GOdf2_list, genome.size=TotalGenesGenome) # annotated genes in Genome v4
  #
  Pval <- getMatrix(go.obj, name="pval")
  Pval <- data.table("GO"=chop(names(Pval), '[.]', 1), "Pval"=Pval)
  
  #
  Common <- getMatrix(go.obj, name="intersection")
  Common <- data.table("GO"=chop(names(Common), '[.]', 1), "CommonGenes" = Common)
  
  # Combine Pvals and Common genes
  out <- dplyr::left_join(Pval, Common, by="GO")
  
  # FDR corrections
  out[,"FDR"] <- p.adjust(out$Pval, method = 'fdr')
  out <- out[out$FDR <= 0.01,]
  
  out[,"GOs"] <- goid
  
  return(out)
}

****************************
** aca voy
****************************
# Define GOs to test
max=length(GOdf2_list)
w=50 # Size of range to test

# Empty list
GO_GO_net <- list()

## parallel Geneset overlap testing 
for (i in seq(0, max, w)){
  # define rank
  Start=i+1
  end=i+w
  
  cat('\n')
  if (end < max){
    print("path 1")
    cat('Start..end:', Start,':',end)
    cat('\n')
    GO_GO_net <- c(GO_GO_net, mclapply(namesGOdf2_list[Start:end], CompareGeneSets, mc.cores=w))
    
  } else{
    print("path 2")
    cat('Start..end:', Start,':',max)
    cat('\n')
    w= max-Start 
    GO_GO_net <- c(GO_GO_net, mclapply(namesGOdf2_list[Start:max], CompareGeneSets, mc.cores=w))
  }
}

GO_GO_net <- rbindlist(GO_GO_net)
colnames(GO_GO_net)[c(1,5)] <- c("GO2","GO1")

# Add genes in GO
GO_GO_net <- left_join(GO_GO_net, GOzmCount_2, by=c('GO1'='GO.ID'))
GO_GO_net <- left_join(GO_GO_net, GOzmCount_2, by=c('GO2'='GO.ID'))

GO_GO_net <- subset(GO_GO_net, GO2 != GO1)

# smallest groups for each comparison
minvals <- apply(GO_GO_net[,c("N.x", "N.y")], 1, min)

GO_GO_net[,"OverlapIndex"] <- GO_GO_net$CommonGenes/minvals

DiscardOverlapIndex  <- GO_GO_net[GO_GO_net$OverlapIndex > 0.7,]
DiscardOverlapIndex <- unique(c(DiscardOverlapIndex$GO2, DiscardOverlapIndex$GO1))

# Select GO-GO to filter by sig and OverlapIndex
GO_GO_net <- GO_GO_net[GO_GO_net$OverlapIndex > 0.7,]

GO_GO_net <- GO_GO_net[,c("GO1", "GO2", "FDR", "CommonGenes", "N.x","N.y", "OverlapIndex")]

# Keep largest GOs by group of interactions
GO_GO_net %>%
  group_by(GO1) %>%
  mutate("MaxSize"= max(c(N.x, N.y))) %>%
  filter( N.x == MaxSize| N.y == MaxSize) -> PassRedundancy

PassRedundancy <- unique(c(PassRedundancy$GO1,PassRedundancy$GO2))

# GOs not DiscardOverlapIndex but present in PassRedundancy: Discarted
DiscardOverlapIndex <- PassRedundancy[DiscardOverlapIndex %in% PassRedundancy]

# Gene-GO that pass second filter
GOzm.DF_3 <- GOzm.DF_2[!(GOzm.DF_2$GO.ID %in% DiscardOverlapIndex),]

# Count Freq after filter
GOzmCount_3 <- data.table(table(GOzm.DF_3$GO.ID)) %>%
  arrange(N) 
colnames(GOzmCount_3)[1] <- "GO.ID"

# GO's description 
GOzmDFdesc <- rbindlist(lapply(GOzmCount_3$GO.ID, function(x) data.table(unlist(x))), idcol = T)

######################## 

######################## 
# Filter 3: Multi-task
########################

# Counts Genes frequency
GenesMultitask <- data.table(table(GOzm.DF_3$GeneId)) %>%
  arrange(N) %>%
  magrittr::set_colnames(c("GeneID", "N"))

rbindlist(lapply(GOBPzmDF_MulF$GO.ID[GOBPzmDF_MulF$N == 22], GetGODescription), fill = T)

## Remove genes with more than 10 functions
# Genes with less than 10 GOs
GOzm.DF_4 <- subset(GOzm.DF_3, !(GeneId %in% GenesMultitask[GenesMultitask$N > 10,]$GeneID))

# Counts Genes frequency
GOzmCount_4 <- data.table(table(GOzm.DF_4$GO.ID)) %>%
  arrange(N) %>%
  magrittr::set_colnames(c("GO.ID", "N"))

####
# Plot Multi task
#######

GenesMultitask %>%
  select(N) %>%
  count(N) %>% as_tibble () %>%
  mutate(Fraction = (n/sum(n)), CumFraction = cumsum((n/sum(n)))) 

GenesMultitask %>%
  select(N) %>%
  count(N) %>% as_tibble () %>%
  mutate(Fraction = (n/sum(n)), CumFraction = cumsum((n/sum(n)))) %>%
  ggplot(aes(x=N, y=n)) +
  geom_point() +
  geom_line(aes(x=N, y=CumFraction*max(n)), color='#FF99FF') +
  scale_y_continuous(labels = scales::comma,
                     # Features of the first axis
                     name = "Gene counts",
                     # Add a second axis and specify its features
                     sec.axis = sec_axis(~. /max(.), name="Cumulative gene fraction")) +
  dark_theme_linedraw(base_line_size = 0.1) +
  geom_text(label="99 % genes", y=3000, x=11, size=5, fontface='plain', color='white') +
  xlab("GOs") +
  geom_vline(xintercept = 10, linetype="dashed", color = "white", size=1) + 
  theme(strip.text.x = element_text(size = 14), 
        axis.text=element_text(size=14),
        text = element_text(size=14), 
        legend.position="none") -> Plot_MultiFunctionally


tiff("Plots/Plot_MultiFunctionally.tiff", units="in", width=5, height=3, res=300)
print(Plot_MultiFunctionally)
dev.off()
#######
######################## 

################################
### Summary Plot GOs filtered
################################

# Summary
length(unique(GOzm.DF$GeneId))
length(unique(GOzm.DF_2$GeneId))
length(unique(GOzm.DF_3$GeneId))
length(unique(GOzm.DF_4$GeneId))

## Compare before and after filters
GOzmCount_1[,"Class"] <- "Raw"
GOzmCount_2[,"Class"] <- "Redundancy"
GOzmCount_3[,"Class"] <- "Specificity"
GOzmCount_4[,"Class"] <- "Multi-functionally"

Input_GOToUse <- rbind(GOzmCount_1, GOzmCount_2, GOzmCount_3, GOzmCount_4[,c("GO.ID", "N","Class")])

Input_GOToUse$Class <- factor(Input_GOToUse$Class, levels = c("Raw", "Redundancy", "Specificity", "Multi-functionally"))

Plot_GOToUse <- ggplot(Input_GOToUse, aes(y=N,  fill=Class)) +
  geom_histogram(bins = 50) +
  scale_fill_manual(values = c("Raw"="gray83", 
                               "Redundancy"="tomato", 
                               "Specificity"="#FFFF66",
                               "Multi-functionally"="#66FF66")) +
  dark_theme_linedraw() + 
  scale_x_continuous(expand = c(0,0), labels = comma) +
  scale_y_log10(labels = comma) +
  ylab("Genes in GO term") + 
  xlab("Counts") +
  annotation_logticks(sides = "l", color = 'white') +
  theme(strip.text.x = element_text(size = 14), 
        axis.text=element_text(size=14, angle = 45, vjust = 0.5, hjust = 1),
        text = element_text(size=14), 
        legend.position="none") +
  facet_grid(. ~ Class, scales = "free_x")

tiff("Plots/Plot_GOToUseFiter_Specificity.tiff", units="in", width=7, height=3, res=300)
print(Plot_GOToUse)
dev.off()

Plot_GOToUse
################################

##############################################################


##############################################################
#########   Training GO 
##############################################################


##############
dfpecanpy <- fread("weighted_pecanpy_ALL.txt", skip = 1, header = F)

# Filter genes with GO annoated
GOs_Pecan <- subset(dfpecanpy, V1 %in% unique(GOzm.DF_4$GeneId))
Geneids <- GOs_Pecan$V1
GOs_Pecan <- as.matrix(GOs_Pecan[,-c(1)])
row.names(GOs_Pecan) <- Geneids


GetGODescription(TARGET_GOs[1])

go = TARGET_GOs[1]
go = 'GO:0044281'



lmg_tf <- function(tfid){  
  #
  ## Discard GOs similar to target GOid
  
  ## Genes not to use
  
  
  # keep genes from GO of interest
  GeneTARGET <- Net_HighConfident$Targ[Net_HighConfident$TF == tfid]
  
  #GeneIN <- row.names(GOs_Pecan)[!(row.names(GOs_Pecan) %in% GeneTARGET)]
  
  
  
  # GO pos lables
  GeneClass <- rep(0, nrow(GOs_Pecan))
  GeneClass[row.names(GOs_Pecan) %in% GeneTARGET] <- 1
  
  cvbeta <- cv.glmnet(GOs_Pecan, GeneClass,  family = "binomial", type.measure = "deviance")
  #plot(cvbeta)
  # 
  out <- as.matrix(coef(cvbeta, s='lambda.min'))
  out <- data.table(Embedding=row.names(out)[-c(1)], B=out[2:nrow(out),])
  colnames(out)[2] <- tfid
  return(out)
}


############################################################
###########            GOs models                 ##########
############################################################
TARGET_GOs <- unique(GOzmCount_4$GO.ID)

# Define GOs to test
max=length(TARGET_GOs)
w=50 # Size of range to test

# Empty list
ModelGO_DB <- list()

## Parallel GeneSet overlap testing 
for (i in seq(0, max, w)){
  # define rank
  Start=i+1
  end=i+w
  
  cat('\n')
  if (end < max){
    print("path 1")
    cat('Start..end:', Start,':',end)
    cat('\n')
    ModelGO_DB <- c(ModelGO_DB, mclapply(TARGET_GOs[Start:end], lmg_go, mc.cores=w))
    
  } else{
    print("path 2")
    cat('Start..end:', Start,':',max)
    cat('\n')
    w= max-Start 
    ModelGO_DB <- c(ModelGO_DB, mclapply(TARGET_GOs[Start:max], lmg_go, mc.cores=w))
  }
}

# Combine results in a wide matrix
ModelGO_DB_wide <- data.table(Embedding=ModelGO_DB[[1]]$Embedding)

for (i in 1:length(ModelGO_DB)){
  ModelGO_DB_wide <- left_join(ModelGO_DB_wide, ModelGO_DB[[i]], by='Embedding')
  
}

# GOs as rows and embedding as cols
ModelGO_DB_wide <- t(as.matrix(ModelGO_DB_wide[,-c(1)]))


GOBP_zma_Description <- rbindlist(lapply(TARGET_GOs, GetGODescription), fill = T)

############################################################
###########            TFs models                 ##########
############################################################

TARGET_tfs <- unique(TFs_HighConfident$.)

# Define GOs to test
max=length(TARGET_tfs)
w=40 # Size of range to test

# Empty list
ModelTF_DB <- list()

## parallel Geneset overlap testing 
for (i in seq(0, max, w)){
  # define rank
  Start=i+1
  end=i+w
  
  cat('\n')
  if (end < max){
    print("path 1")
    cat('Start..end:', Start,':',end)
    cat('\n')
    ModelTF_DB <- c(ModelTF_DB, mclapply(TARGET_tfs[Start:end], lmg_tf, mc.cores=w))
    
  } else{
    print("path 2")
    cat('Start..end:', Start,':',max)
    cat('\n')
    w= max-Start 
    ModelTF_DB <- c(ModelTF_DB, mclapply(TARGET_tfs[Start:max], lmg_tf, mc.cores=w))
  }
}

names(ModelTF_DB) <- TARGET_tfs

#Model_P1 <- lmg_tf("Zm00001d028854")
#ModelTF_DB[["Zm00001d028854"]] <- Model_P1

ReplaceName(TARGET_tfs)

############################################################

############################################################
###########        Cosine distance                 #########
############################################################

# Calculate cos distance 
CosDis <- function(TFid){
  
  tfvasl <-as.numeric(unlist(ModelTF_DB[[TFid]][,2]))
  
  
  cosVecor <- apply(ModelGO_DB_wide, 1, function(x) cosine(tfvasl, x)) 
  
  cosVecor = data.table(GO=names(cosVecor), Cos=cosVecor)
  cosVecor$Cos[is.na(cosVecor$Cos)] <- 0
  
  cosVecor[,"TF"] <- TFid
  
  return(cosVecor)
}

#GO_TF_CosDis_DB <- CosDis("Zm00001d006236")

GO_TF_CosDis_DB <- lapply(TARGET_tfs, CosDis)
GO_TF_CosDis_DB <- rbindlist(GO_TF_CosDis_DB)

GO_TF_CosDis_DB_wide <- WideAndScale(GO_TF_CosDis_DB)

GO_TF_CosDis_DBz <- as_tibble(melt(as.matrix(GO_TF_CosDis_DB_wide)))

#GO_TF_CosDis_DBz <- subset(GO_TF_CosDis_DBz, !(Var1 %in% c("GO:0090558", "GO:0044451", "GO:0001085")))

GO_TF_CosDis_DBz %>% 
  group_by(Var2) %>%
  dplyr::slice_max(value, n=10) -> TopGO_TFs

TopGO_TFs[,1:2] <- apply(TopGO_TFs[,1:2], 2, as.character)

# Add TFs names
TopGO_TFs[, "TFname"] <- ReplaceName(TopGO_TFs$Var2)

# add GO descriptions
TopGO_TFs <- left_join(TopGO_TFs, GOBP_zma_Description, by=c("Var1"="GO.ID"))

TopGO_TFs[(TopGO_TFs$TFname=='KN1'), ]

tfcandiates <- c("MYB31", "KN1", "COL18")


ggplot(subset(TopGO_TFs, TFname %in% tfcandiates),
       aes(x=TFname, y=Term, fill=value)) +
  geom_tile() +
  scale_fill_viridis("rZ", option = 'C', discrete = F) +
  scale_x_discrete(expand = c(0,0)) +
  ylab("GO term") + xlab("") +
  dark_theme_minimal() + 
  theme(strip.text.y = element_text(size = 14), 
        axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 1),
        text = element_text(size=14), legend.title = element_text(colour = "white")) -> Plot_TopExample

Plot_TopExample

tiff("Plots/Plot_ML_heatmap_example.tiff", units="in", width=8, height=5, res=300)
print(Plot_TopExample)
dev.off()


# "Flavonoid biosynthetic process"

GO_TF_CosDis_DBz 

GO_TF_CosDis_DB %>% 
  dplyr::arrange(Cos) %>%
  dplyr::slice_min(Cos, n=10) -> TailGO_MYB31

TopGO_MYB31
TailGO_MYB31$GO

rbindlist(lapply(TailGO_MYB31$GO, GetGODescription), fill = T)

pheatmap::pheatmap(GO_TF_CosDis_DB_wide)

