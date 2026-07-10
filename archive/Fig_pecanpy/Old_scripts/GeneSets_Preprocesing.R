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

#library(org.At.tair.db)
# file:///Users/fabiogomezcano/Downloads/mdgsa_vignette.pdf
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

GetGODescription <- function(id) {
  
  # If go in GO.db, then get info
  if( id %in% names(GOID(GOTERM)) ) {
    out <- data.table("GO.ID"=GOID(GOTERM[[id]]),
                      "Term" = Term(GOTERM[[id]]),
                      "Ontology"= Ontology(GOTERM[[id]]))
    return(out)
  }
  
}

# Make ChildDB
ChildDB <- as.list(GOBPOFFSPRING)
head(ChildDB)

# Define names
ParentsNames <- names(ChildDB)

# replace paretn without child with own name
ChildDB[is.na(ChildDB)] <- ParentsNames[is.na(ChildDB)]

# as data table
ChildDB <- rbindlist(lapply(ChildDB, function(x) data.table(unlist(x))), idcol = T)
colnames(ChildDB) <- c("Parent", "Child")


GetGODescription("GO:0032042")


# PDI
GRN <- as_tibble(unique(fread("../Fig_PDI/PDI_NetworkFinal.10_11_2021.txt")[,c(2,3)]))
colnames(GRN)[1] <- "Source"

# CoExp
CEN <- as_tibble(unique(fread("../Fig_Coexpression/CoExp_NetworkFinal.10_11_2021.txt")))
colnames(CEN)[2] <- "Source"
CEN <- unique(CEN[,2:3])

# teQTL
GAN <- as_tibble(unique(fread("../Fig_transeQTL/teQTL_NetworkFinal.10_11_2021.txt")))
colnames(GAN)[1] <- "Source"

# teQTL associated with TFs
teQTLtf <- subset(teQTL, Source %in% unique(c(TF_CoR$GeneID, PDI$Source, CoExp$Source))) 

# https://ssayols.github.io/rrvgo/articles/rrvgo.html
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
GOzmList <- rbindlist(lapply(GOzm, function(x) data.table(unlist(x))), idcol = T)
colnames(GOzmList) <- c("GeneId", 'GO.ID')
GOzmList <- unique(GOzmList)

subset(GOzmList, GO.ID=="GO:0008152")
table(subset(GOzmList, GO.ID=="GO:0016051")$GeneId %in% subset(GOzmList, GO.ID=="GO:0005975")$GeneId)


# Count freq. as formatted it as DF
GOzmDF <- data.table(table(GOzmList$GO.ID)) %>%
  arrange(N) 
colnames(GOzmDF)[1] <- "GO.ID"
  
# Get description and ontology for GOs in Zm genome
GOBP_zma_Description <- rbindlist(lapply(GOzmDF$GO.ID, GetGODescription))

# Combined GO size with descriptions
GOzmDF <- left_join(GOzmDF, GOBP_zma_Description, by="GO.ID")
GOzmDF$Ontology[is.na(GOzmDF$Ontology)] <- "NN"

# Filter to only BP
GOBPzmDF <- subset(GOzmDF, Ontology %in% c("BP", "NN")) 

# Define GO too specific and too general
GOBPzmDF_Small <- subset(GOBPzmDF, N<=20)
#GOBPzmDF_Large <- subset(GOBPzmDF, N >= 400)
GOBPzmDF_ready <- subset(GOBPzmDF, N <= 300 & N > 20)


# Filter Genes:GO DF to only GOS in BP and NN
GOzmList <- subset(GOzmList, GO.ID %in% unique(GOBPzmDF$GO.ID))


## Child to parent parent with large GO sizes
GOzmList_small <- subset(GOzmList, GO.ID %in% GOBPzmDF_Small$GO.ID)
GOzmList_ready <- subset(GOzmList, GO.ID %in% GOBPzmDF_ready$GO.ID)


## Add Childs to large GOs
GOzmList_small <- left_join(GOzmList_small, ChildDB, by=c("GO.ID"="Child"))
GOzmList_small <- unique(GOzmList_small)

# keep only Syntenic genes
GOzmList_small <- subset(GOzmList_small, GeneId %in% Syntenic)


# Count freq. of parent terms
small_2_keep <- data.table(table(unique(GOzmList_small[,c(1,3)])$Parent)) %>%
  arrange(N) %>%
  magrittr::set_colnames(c("Parent", "Np")) %>%
  dplyr::filter(Np > 20 & Np < 300)

# Reduce to parent with size expected
GOzmList_small <- GOzmList_small[GOzmList_small$Parent %in% small_2_keep$Parent,]

# Replace GO child wit its corresponding parent
GOzmList_small$GO.ID <- GOzmList_small$Parent
GOzmList_small <- GOzmList_small[,-c(3)]
GOzmList_small <- unique(GOzmList_small)

## Combine ready with new small mapped to its corresponding parents
GOzmList_new <- unique(rbind(GOzmList_ready, GOzmList_small))

# define gene sets after cleaning
GOzmList_new <- subset(GOzmList_new, GeneId %in% Syntenic)

# Counts Genes by GOs
GOBPzmDF_new <- data.table(table(GOzmList_new$GO.ID)) %>%
  arrange(N) %>%
  magrittr::set_colnames(c("GO.ID", "N"))


##############################################################

##############################################################
#########      Process GO terms: Redundancy      #############
##############################################################
chop=function(myStr,mySep,myField){
  
  choppedString=sapply(strsplit(myStr,mySep),"[",myField)
  if(length(myField)>1){
    choppedString=apply(choppedString,2,function(x){paste0(x[!is.na(x)],collapse=mySep)})
  }
  return(choppedString)
}

TotalGenesGenome <- length(Syntenic)

GO_list <- split(GOzmList_new$GeneId, GOzmList_new$GO.ID)

GO_round_1 <- names(GO_list)

CompareGeneSets <- function(goid){
  
  ## Compare list of predicted targets vs annotated genes in query.vector 
  #
  print(". Pre-newGOM .")
  go.obj <- newGOM(GO_list[goid], GO_list, genome.size=TotalGenesGenome) # annotated genes in Genome v4
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

# Define GOs to test
max=length(GO_round_1)
w=50 # Size of range to test

# Empty list
GO_GO_net_1 <- list()

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
    GO_GO_net_1 <- c(GO_GO_net_1, mclapply(GO_round_1[Start:end], CompareGeneSets, mc.cores=w))
    
  } else{
    print("path 2")
    cat('Start..end:', Start,':',max)
    cat('\n')
    w= max-Start 
    GO_GO_net_1 <- c(GO_GO_net_1, mclapply(GO_round_1[Start:max], CompareGeneSets, mc.cores=w))
  }
}

GO_GO_net_1 <- rbindlist(GO_GO_net_1)
colnames(GO_GO_net_1)[c(1,5)] <- c("GO2","GO1")

# GO_GO_net_1 <- GO_GO_net_1[GO_GO_net_1$GO1 %in% GO_round_1,]
# GO_GO_net_1 <- GO_GO_net_1[GO_GO_net_1$GO2 %in% GO_round_1,]

# Add genes in GO
GO_GO_net_1 <- left_join(GO_GO_net_1, GOBPzmDF_new, by=c('GO1'='GO.ID'))
GO_GO_net_1 <- left_join(GO_GO_net_1, GOBPzmDF_new, by=c('GO2'='GO.ID'))

GO_GO_net_1 <- subset(GO_GO_net_1, GO2 != GO1)

minvals <- apply(GO_GO_net_1[,c("N.x", "N.y")], 1, min)

GO_GO_net_1[,"OverlapIndex"] <- GO_GO_net_1$CommonGenes/minvals

# filter by 
GO_GO_net_1 <- GO_GO_net_1[GO_GO_net_1$OverlapIndex > 0.7,]
GO_GO_net_1 <- GO_GO_net_1[,c("GO1", "GO2", "FDR", "CommonGenes", "N.x","N.y", "OverlapIndex")]

# Keep largest GOs by group of interactions
GO_GO_net_1 %>%
  group_by(GO1) %>%
  mutate("MaxSize"= max(c(N.x, N.y))) %>%
  filter( N.x == MaxSize| N.y == MaxSize) -> GOs_2_keep


# Select GOs equal to max size from GO1 and GO2
GOs_2_keep <- unique(c(GOs_2_keep$GO1[GOs_2_keep$N.x == GOs_2_keep$MaxSize], 
                       GOs_2_keep$GO2[GOs_2_keep$N.y == GOs_2_keep$MaxSize]))
length(GOs_2_keep)

# Define GOs to keep after first round
GO_Not_InGO_GO <- GO_round_1[!(GO_round_1 %in% unique(c(GO_GO_net_1$GO1, GO_GO_net_1$GO2)))]


GO_After_1 <- unique(c(GO_Not_InGO_GO, GOs_2_keep))
GOzmList_After_1 <- GO_list[GO_After_1]

# redefine as DF
GOzmList_After_1 <- rbindlist(lapply(GOzmList_After_1, function(x) data.table(unlist(x))), idcol = T)
colnames(GOzmList_After_1) <- c("GO.ID", "GeneID")


GOBPzmDF_After1 <- data.table(table(GOzmList_After_1$GO.ID)) %>%
  arrange(N) %>%
  magrittr::set_colnames(c("GO.ID", "N"))

##############################################################
#########      Process GO terms: Redundancy      #############
##############################################################

# Counts Genes frequency
GenesMultitask <- data.table(table(GOzmList_After_1$GeneID)) %>%
  arrange(N) %>%
  magrittr::set_colnames(c("GeneID", "N"))


GenesMultitask %>%
  select(N) %>%
  count(N) %>% as_tibble () %>%
  mutate(Fraction = (n/sum(n)), CumFraction = cumsum((n/sum(n)))) %>% View()

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
  geom_text(label="92 % genes", y=3000, x=25, size=5, fontface='plain', color='white') +
  xlab("GOs") +
  geom_vline(xintercept = 15, linetype="dashed", color = "white", size=1) + 
  theme(strip.text.x = element_text(size = 14), 
        axis.text=element_text(size=14),
        text = element_text(size=14), 
        legend.position="none") -> Plot_MultiFunctionally
  

tiff("Plots/Plot_MultiFunctionally.tiff", units="in", width=5, height=3, res=300)
print(Plot_MultiFunctionally)
dev.off()

#rbindlist(lapply(GOzmList_After_1[GOzmList_After_1$GeneID == "Zm00001d001780", ]$GO.ID, GetGODescription), fill = T)

rbindlist(lapply(GOBPzmDF_MulF$GO.ID[GOBPzmDF_MulF$N == 22], GetGODescription), fill = T)

GetGODescription("GO:0000122")
## Remove genes with more than 15 functions


# genes with less than 20 GOs
GOzmList_MulF <- subset(GOzmList_After_1, (GeneID %in% GenesMultitask[GenesMultitask$N <= 20,]$GeneID))

GOBPzmDF_MulF <- data.table(table(GOzmList_MulF$GO.ID)) %>%
  arrange(N) %>%
  magrittr::set_colnames(c("GO.ID", "N")) %>%
  filter(N>=20)

GOzmList_MulF <- subset(GOzmList_MulF, GO.ID %in% unique(GOBPzmDF_MulF$GO.ID))

length(unique(GOzmList_new$GeneId))
length(unique(GOzmList_After_1$GeneID))
length(unique(GOzmList_MulF$GeneID))


## Compare before and after filters
GOBPzmDF_MulF[,"Class"] <- "Multi-functionally"
GOBPzmDF_After1[,"Class"] <- "Specificity"
GOBPzmDF_new[,"Class"] <- "Redundancy"
GOBPzmDF[,"Class"] <- "Raw"

GOBPzmDF_MulF

length(unique(GOzmList$GeneId))

go="GO:0016192"
GetGODescription(go)

GOBPzmDF_MulF[GOBPzmDF_MulF$GO.ID == go]
GOBPzmDF_MulF$GO.ID[GOBPzmDF_MulF$GO.ID %in% GO_GO_net_1$GO1]

###########
Input_GOToUse <- rbind(GOBPzmDF_new, GOBPzmDF_After1, GOBPzmDF_MulF, GOBPzmDF[,c("GO.ID", "N","Class")])

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

###########

#################################################################################################
### Redefine group of GOs with overlapping number of genes  to remove them from the regression ##
#################################################################################################

TotalGenesGenome <- length(Syntenic)

GO_list_Clean <- split(GOzmList_MulF$GeneID, GOzmList_MulF$GO.ID)

GO_round_Clean <- names(GO_list_Clean)

TotalGenes = length(unique(GOzmList_MulF$GeneID))

CompareGeneSets2 <- function(goid){
  
  ## Compare list of predicted targets vs annotated genes in query.vector 
  #
  
  print(". Pre-newGOM .")
  go.obj <- newGOM(GO_list_Clean[goid], GO_list_Clean, genome.size=TotalGenes) # annotated genes in Genome v4
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
  out <- out[out$FDR <= 0.05,]
  
  out[,"GOs"] <- goid
  
  return(out)
}

# Define GOs to test
max=length(GO_round_Clean)
w=50 # Size of range to test

# Empty list
GO_GO_net_Clean <- list()

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
    GO_GO_net_Clean <- c(GO_GO_net_Clean, mclapply(GO_round_Clean[Start:end], CompareGeneSets2, mc.cores=w))
    
  } else{
    print("path 2")
    cat('Start..end:', Start,':',max)
    cat('\n')
    w= max-Start 
    GO_GO_net_Clean <- c(GO_GO_net_Clean, mclapply(GO_round_Clean[Start:max], CompareGeneSets2, mc.cores=w))
  }
}

GO_GO_net_Clean <- rbindlist(GO_GO_net_Clean)
colnames(GO_GO_net_Clean)[c(1,5)] <- c("GO2","GO1")

#################################################################################################

##############
dfpecanpy <- fread("weighted_pecanpy_ALL.txt", skip = 1, header = F)

# Filter genes with GO annoated
GOs_Pecan <- subset(dfpecanpy, V1 %in% GOzmList_MulF$GeneID)
Geneids <- GOs_Pecan$V1
GOs_Pecan <- as.matrix(GOs_Pecan[,-c(1)])
row.names(GOs_Pecan) <- Geneids
library(glmnet)


lmg_go <- function(go){  
  #
  ## Discard GOs similar to target GOid
  GOsNotouse <- unique(c(subset(GO_GO_net_Clean, GO1 == go)$GO2,
                  subset(GO_GO_net_Clean, GO2 == go)$GO1))
  
  ## Genes not to use
  GeneOUT <- unique(GOzmList_MulF$GeneID[(GOzmList_MulF$GO.ID %in% GOsNotouse)])
  
  # keep genes from GO of interest
  GeneTARGET <- unique(GOzmList_MulF$GeneID[(GOzmList_MulF$GO.ID %in% go)])
  GeneOUT <- GeneOUT[!(GeneOUT %in% GeneTARGET)]
  
  #
  Input <- GOs_Pecan[ !(row.names(GOs_Pecan) %in% GeneOUT),]
  
  # GO pos lables
  GeneClass <- rep(0, nrow(Input))
  GeneClass[row.names(Input) %in% GeneTARGET] <- 1
  
  cvbeta <- cv.glmnet(Input, GeneClass,  family = "binomial", type.measure = "deviance")
  # 
  out <- as.matrix(coef(cvbeta, s='lambda.min'))
  out <- data.table(Embedding=row.names(out)[-c(1)], B=out[2:nrow(out),])
  colnames(out)[2] <- go
  return(out)
}

lmg_tf <- function(tfid){  
  #
  ## Discard GOs similar to target GOid
  
  ## Genes not to use
  
  
  # keep genes from GO of interest
  GeneTARGET <- GRN$Target[GRN$Source == tfid]
  
  #GeneIN <- row.names(GOs_Pecan)[!(row.names(GOs_Pecan) %in% GeneTARGET)]
  
  #
  Input <- GOs_Pecan# [GeneIN,]
  
  # GO pos lables
  GeneClass <- rep(0, nrow(Input))
  GeneClass[row.names(Input) %in% GeneTARGET] <- 1
  
  cvbeta <- cv.glmnet(Input, GeneClass,  family = "binomial", type.measure = "deviance")
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
TARGET_GOs <- unique(GOBPzmDF_MulF$GO.ID)

# Define GOs to test
max=length(TARGET_GOs)
w=50 # Size of range to test

# Empty list
ModelGO_DB <- list()

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
    ModelGO_DB <- c(ModelGO_DB, mclapply(TARGET_GOs[Start:end], lmg_go, mc.cores=w))
    
  } else{
    print("path 2")
    cat('Start..end:', Start,':',max)
    cat('\n')
    w= max-Start 
    ModelGO_DB <- c(ModelGO_DB, mclapply(TARGET_GOs[Start:max], lmg_go, mc.cores=w))
  }
}

# combine resuls in a wide matrix
ModelGO_DB_wide <- data.table(Embedding=ModelGO_DB[[1]]$Embedding)


for (i in 1:length(ModelGO_DB)){
  ModelGO_DB_wide <- left_join(ModelGO_DB_wide, ModelGO_DB[[i]], by='Embedding')
  
}

# GOs as rows and embedding as cols
ModelGO_DB_wide <- t(as.matrix(ModelGO_DB_wide[,-c(1)]))



############################################################
###########            TFs models                 ##########
############################################################

TARGET_tfs <- unique(GRN$Source)
TARGET_tfs <- TARGET_tfs[TARGET_tfs %in% unique(CEN$Source)] 
TARGET_tfs <- TARGET_tfs[TARGET_tfs %in% unique(GAN$Source)] 

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

Model_P1 <- lmg_tf("Zm00001d028854")

ModelTF_DB[["Zm00001d028854"]] <- Model_P1


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
Model_P1 <- 

GO_TF_CosDis_DB <- CosDis("Zm00001d006236")

GO_TF_CosDis_DB <- lapply(c(TARGET_tfs, "Zm00001d028854"), CosDis)

GO_TF_CosDis_DB <- rbindlist(GO_TF_CosDis_DB)

GO_TF_CosDis_DB_wide <- WideAndScale(GO_TF_CosDis_DB)

GO_TF_CosDis_DBz <- as_tibble(melt(as.matrix(GO_TF_CosDis_DB_wide)))



GO_TF_CosDis_DBz <- subset(GO_TF_CosDis_DBz, !(Var1 %in% c("GO:0090558", "GO:0044451", "GO:0001085")))

GO_TF_CosDis_DBz %>% 
  group_by(Var2) %>%
  dplyr::slice_max(value, n=5) -> TopGO_TFs

TopGO_TFs[, "TFname"] <- ReplaceName(TopGO_TFs$Var2)

TopGO_TFs$Var1 <- as.character(TopGO_TFs$Var1)
TopGO_TFs$Var2 <- as.character(TopGO_TFs$Var2)

TopGO_TFs <- left_join(TopGO_TFs, GOBP_zma_Description, by=c("Var1"="GO.ID"))

TopGO_TFs[(TopGO_TFs$TFname=='KN1'), ][3,5] <- "Regulation of reproductive process"
TopGO_TFs[(TopGO_TFs$TFname=='P1'), ][4,5] <- "Organic acid transport"

tfcandiates <- c("MYB31", "KN1", "COL18", "P1")



ggplot(subset(TopGO_TFs, TFname %in% tfcandiates),
       aes(x=TFname, y=Term, fill=value)) +
  geom_tile() +
  scale_fill_viridis("rZ", option = 'B', discrete = F) +
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


"flavonoid biosynthetic process"

GO_TF_CosDis_DBz 

GO_TF_CosDis_DB %>% 
  dplyr::arrange(Cos) %>%
  dplyr::slice_min(Cos, n=10) -> TailGO_MYB31

TopGO_MYB31
TailGO_MYB31$GO
  
rbindlist(lapply(TailGO_MYB31$GO, GetGODescription), fill = T)

pheatmap::pheatmap(GO_TF_CosDis_DB_wide)

