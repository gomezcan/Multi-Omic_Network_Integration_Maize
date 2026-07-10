#library(tidygraph)
#library(ggrepel)
#library(ggpubr)
#library(GeneOverlap)
#library(viridis)
#library(ComplexHeatmap)
#library(fgsea)
#library(reshape)
library(scales)
library(tidyverse)
library(data.table)


ReplaceName <- function(ids){
  
  for (i in 1:nrow(Top45)){
    ids <- gsub(Top45$V1[i], Top45$V2[i], ids)
  }
  
  for (i in 1:nrow(TFdic)){
    ids <- gsub(TFdic$TF.v4[i], TFdic$TF.Name[i], ids)
  }
  return(ids)
}

trans = trans.eQTL

Get_transcis <- function(trans){
  
  #trans <- trans[1:10,]
  trans[,"tem"] <- paste(trans$Target, trans$source, trans$Index, sep = "_")
  out <- trans
  
  #
  trans <- trans[,c("Target","source","Index", "tem")]
  
  # Add coordenates
  trans <- left_join(trans, saf, by=c("Target"="V1"))
  trans <- left_join(trans, saf, by=c("source"="V1"))
  
  # define trans-eQTL because of same chr
  trans <- subset(trans, V2.x == V2.y) # same chr
  
  # distance in bps
  trans[,"DisS_T"] <- abs(trans$V3.x - trans$V3.y)
  
  # Selects targen::snp_index and disS_T: ditance to Start of termination site
  trans <- unique(trans[,c("tem", "DisS_T")])
  
  
  # keep initial format for defined trans-eQTLs
  out <- subset(out, tem %in% trans$tem)
  
  # Add DisS_T value
  out <- left_join(out, trans, by='tem')
  
  # remove index col
  out <- out[, -c("tem")]
  
  # redefine class from "trasneQTL" to "transciseQTL"
  out[,"class"] <- "transciseQTL"
  
  # trans-cis: defined as those close proximity to target genes
  out <- subset(out, abs(DisS_T) <= 50000) 
  out <- unique(out)
  
  return(out)
}

subset(Net_teQTL, source=='Zm00001d037784')

######################################################################################
#########       Read different eQTL data sets and define summary tables      #########
######################################################################################
## syntenic genes 
Syntenic <- as_tibble(read.table("../Annotations/Zm.v4.synteny.genes.txt", h=T, stringsAsFactors = F, sep = '\t'))
Syntenic <- subset(Syntenic, gene_synteny=='syntenic')$gene_id

#
saf <- as_tibble(read.table("Zea_mays.B73_RefGen_v4.46.saf", stringsAsFactors = F))

# Total cis- & trans-eQTLs
eQTLs <- fread("Final_cis_trans_all_eQTL_10012021.bed", stringsAsFactors = F, header = T)
colnames(eQTLs) <- c("chr",'s', "e","snp", "Target", "Support")

eQTLs[,"Index"] <- paste(eQTLs$chr, eQTLs$s, eQTLs$snp, sep = ":")
eQTLs <- eQTLs[,-c(1,2,3)]
eQTLs <- unique(eQTLs[,c(2,4)])

# trans-eQTL
trans.eQTL <- fread("trans.eQTL_target.txt", stringsAsFactors = F)
colnames(trans.eQTL) <- c("chr",'s','e','snp','Target','Support','chr.source', 'tss.source', 'e.source', 'source',"score","strand","dis","class")
trans.eQTL[,"Index"] <- paste(trans.eQTL$chr, trans.eQTL$s, trans.eQTL$snp, sep = ":")
trans.eQTL <- unique(trans.eQTL[,-c(1,2,3,4,9,11,13)])

# trans-eQTL: Second Part
trans.eQTL_SP <- fread("trans.SecondPart.eQTL_target.txt", stringsAsFactors = F)
colnames(trans.eQTL_SP) <- c("chr",'s','e','snp','Target','Support','chr.source', 
                             'tss.source', 'e.source', 'source',"score","strand","dis","class")

trans.eQTL_SP[,"Index"] <- paste(trans.eQTL_SP$chr, trans.eQTL_SP$s, trans.eQTL_SP$snp, sep = ":")
trans.eQTL_SP <- unique(trans.eQTL_SP[,-c(1,2,3,4,9,11,13)])

# transcis-eQTL definition 
transcis_eQTL <- Get_transcis(trans.eQTL)

# transcis-eQTL definition 
transcis_eQTL_SP <- Get_transcis(trans.eQTL_SP)

# Redefine trans-eQTL
trans.eQTL[,"tem"] <- paste(trans.eQTL$Target, trans.eQTL$source, trans.eQTL$Index, sep = "_")
transcis_eQTL[,"tem"] <- paste(transcis_eQTL$Target, transcis_eQTL$source, transcis_eQTL$Index, sep = "_")

trans.eQTL_SP[,"tem"] <- paste(trans.eQTL_SP$Target, trans.eQTL_SP$source, trans.eQTL_SP$Index, sep = "_")
transcis_eQTL_SP[,"tem"] <- paste(transcis_eQTL_SP$Target, transcis_eQTL_SP$source, transcis_eQTL_SP$Index, sep = "_")

trans.eQTL <- subset(trans.eQTL, !(tem %in% transcis_eQTL$tem))
trans.eQTL_SP <- subset(trans.eQTL_SP, !(tem %in% transcis_eQTL_SP$tem))


# trans-eQTLp
trans.eQTLp <- fread("trans.eQTLp_target.txt", stringsAsFactors = F)
colnames(trans.eQTLp) <- c("chr",'s','e','snp','Target','Support','chr.source', 'tss.source', 'e.source', 'source',"score","strand","dis","class")
trans.eQTLp[,"Index"] <- paste(trans.eQTLp$chr, trans.eQTLp$s, trans.eQTLp$snp, sep = ":")
trans.eQTLp <- unique(trans.eQTLp[,-c(1,2,3,4,9,11,13)])

# trans-eQTLp Second Part
trans.eQTLp_SP <- fread("trans.SecondPart.eQTLp_target.txt", stringsAsFactors = F)
colnames(trans.eQTLp_SP) <- c("chr",'s','e','snp','Target','Support','chr.source', 'tss.source', 'e.source', 'source',"score","strand","dis","class")
trans.eQTLp_SP[,"Index"] <- paste(trans.eQTLp_SP$chr, trans.eQTLp_SP$s, trans.eQTLp_SP$snp, sep = ":")
trans.eQTLp_SP <- unique(trans.eQTLp_SP[,-c(1,2,3,4,9,11,13)])

# define target_source_snp index
trans.eQTLp[,"tem"] <- paste(trans.eQTLp$Target, trans.eQTLp$source, trans.eQTLp$Index, sep = "_")
trans.eQTLp_SP[,"tem"] <- paste(trans.eQTLp_SP$Target, trans.eQTLp_SP$source, trans.eQTLp_SP$Index, sep = "_")


# Cis-eQTLt
cis.eQTLt <- fread("cis.eQTLt_target.txt", stringsAsFactors = F)
colnames(cis.eQTLt) <- c("chr",'s','e','snp','Target','Support','chr.source', 'tss.source', 'e.source', 'source',"score","strand","dis","class")
cis.eQTLt[,"Index"] <- paste(cis.eQTLt$chr, cis.eQTLt$s, cis.eQTLt$snp, sep = ":")
cis.eQTLt <- unique(cis.eQTLt[,-c(1,2,3,4,9,11,13)])
cis.eQTLt[,"tem"] <- paste(cis.eQTLt$Target, cis.eQTLt$source, cis.eQTLt$Index, sep = "_")

# cis-eQTL
cis.eQTL <- fread("ciseQTL_noFiter.txt", stringsAsFactors = F)
colnames(cis.eQTL) <- c("chr",'s','e','snp','Target','Support','chr.source', 'tss.source', 'e.source', 'source',"score","strand","dis","class")
cis.eQTL[,"Index"] <- paste(cis.eQTL$chr, cis.eQTL$s, cis.eQTL$snp, sep = ":")
cis.eQTL <- subset(cis.eQTL, abs(dis) <= 50000)
cis.eQTL <- unique(cis.eQTL[,-c(1,2,3,4,9,11)])
cis.eQTL[,"tem"] <- paste(cis.eQTL$Target, cis.eQTL$source, cis.eQTL$Index, sep = "_")

#
table(trans.eQTL$Support)
table(trans.eQTL_SP$Support)


trans.eQTL <- unique(trans.eQTL[, -c(2)])   # remove Support
trans.eQTLp <- unique(trans.eQTLp[, -c(2)]) # remove Support 
transcis_eQTL <- unique(transcis_eQTL[, -c(2)]) # remove Support 

trans.eQTL_SP <- unique(trans.eQTL_SP[, -c(2)])   # remove Support
trans.eQTLp_SP <- unique(trans.eQTLp_SP[, -c(2)]) # remove Support 
transcis_eQTL_SP <- unique(transcis_eQTL_SP[, -c(2)]) # remove Support 

cis.eQTL <- unique(cis.eQTL[, -c(2)]) # remove Support 
cis.eQTLt <- unique(cis.eQTLt[, -c(2)]) # remove Support 

# Confirm unique values
table(trans.eQTL$tem %in% cis.eQTL$tem)      # trans in cis target?
table(cis.eQTL$tem %in% cis.eQTLt$tem)       # cis in cis target?
table(cis.eQTL$tem %in% trans.eQTLp$tem)     # cis in trans promoter?
table(cis.eQTLt$tem %in% trans.eQTL$tem)     # cis target in trans?

table(trans.eQTL_SP$tem %in% cis.eQTL$tem)    # trans_SP in cis target?
table(cis.eQTL$tem %in% trans.eQTLp_SP$tem)   # cis in trans promoter SP?
table(cis.eQTLt$tem %in% trans.eQTL_SP$tem)   # cis target in trans SP?


table(trans.eQTLp$tem %in% cis.eQTL$tem)      # trans promoter in cis?
table(trans.eQTLp_SP$tem %in% cis.eQTL$tem)   # trans promoter in cis?

table(cis.eQTL$tem %in% cis.eQTLt$tem)        # cis in cis target?
table(cis.eQTLt$tem %in% cis.eQTL$tem)        # cis in cis target?

#############################
## Count eQTL not assigned ##
#############################

eQTLSyn <- subset(eQTLs, Target %in% Syntenic)

#dim(eQTLSyn)
#length(unique(eQTLSyn$Index))
#length(unique(eQTLSyn$Target))

# Define eQTL eQTL-target unassigned
# eQTL in trans-eQTL
eQTLSyn <- eQTLSyn[!(paste(eQTLSyn$Target, eQTLSyn$Index, sep = "_") %in% paste(trans.eQTL$Target, trans.eQTL$Index, sep = "_")),]

# eQTL in trans-eQTL_SP
eQTLSyn <- eQTLSyn[!(paste(eQTLSyn$Target, eQTLSyn$Index, sep = "_") %in% paste(trans.eQTL_SP$Target, trans.eQTL_SP$Index, sep = "_")),]

# eQTL in trans-eQTLp
eQTLSyn <- eQTLSyn[!(paste(eQTLSyn$Target, eQTLSyn$Index, sep = "_") %in% paste(trans.eQTLp$Target, trans.eQTLp$Index, sep = "_")),]
eQTLSyn <- eQTLSyn[!(paste(eQTLSyn$Target, eQTLSyn$Index, sep = "_") %in% paste(trans.eQTLp_SP$Target, trans.eQTLp_SP$Index, sep = "_")),]

# eQTL in trans-eQTLp
eQTLSyn <- eQTLSyn[!(paste(eQTLSyn$Target, eQTLSyn$Index, sep = "_") %in% paste(transcis_eQTL$Target, transcis_eQTL$Index, sep = "_")),]
eQTLSyn <- eQTLSyn[!(paste(eQTLSyn$Target, eQTLSyn$Index, sep = "_") %in% paste(transcis_eQTL_SP$Target, transcis_eQTL_SP$Index, sep = "_")),]

# eQTL in cis-eQTL
eQTLSyn <- eQTLSyn[!(paste(eQTLSyn$Target, eQTLSyn$Index, sep = "_") %in% paste(cis.eQTL$Target, cis.eQTL$Index, sep = "_")),]

# eQTL in cis-eQTLt
eQTLSyn <- eQTLSyn[!(paste(eQTLSyn$Target, eQTLSyn$Index, sep = "_") %in% paste(cis.eQTLt$Target, cis.eQTLt$Index, sep = "_")),]

# Combined first ans second group of trans-eQTLs
trans.eQTL <- rbind(trans.eQTL, trans.eQTL_SP)
trans.eQTLp <- rbind(trans.eQTLp, trans.eQTLp_SP)

fwrite(eQTLSyn,                  "Clean_Unasigned_eQTL.v2.txt", sep = "\t", row.names = F, quote = F)
fwrite(trans.eQTL[,-c("tem")],   "Clean_trans.eQTL.v2.txt", sep = "\t", row.names = F, quote = F)
fwrite(trans.eQTLp[,-c("tem")],  "Clean_trans.eQTLp.v2.txt", sep = "\t", row.names = F, quote = F)
fwrite(transcis_eQTL[,-c("tem")], "Clean_transcis_eQTL.v2.txt", sep = "\t", row.names = F, quote = F)
fwrite(cis.eQTL[,-c("tem")],     "Clean_cis.eQTL.v2.txt", sep = "\t", row.names = F, quote = F)
fwrite(cis.eQTLt[,-c("tem")],    "Clean_cis.eQTLt.v2.txt", sep = "\t", row.names = F, quote = F)

cat("Total eQTL-target associations:",(nrow(eQTLSyn) + nrow(trans.eQTL) + nrow(trans.eQTLp) + nrow(transcis_eQTL) + nrow(cis.eQTL) + nrow(cis.eQTLt))/1e6)
cat("Total SNPs: ", length(unique(c(eQTLSyn$Index, trans.eQTL$Index, trans.eQTLp$Index, transcis_eQTL$Index, cis.eQTL$Index, cis.eQTLt$Index)))/1e6)
cat("Total genes: ", length(unique(c(eQTLSyn$Target, trans.eQTL$Target, trans.eQTLp$Target, transcis_eQTL$Target, cis.eQTL$Target, cis.eQTLt$Target)))/1e6)


# 1180857 Clean_cis.eQTLt.txt
# 3540659 Clean_cis.eQTL.txt
# 10235805 Clean_eQTL_Syntenic.txt
# 1208608 Clean_transcis_eQTL.txt
# 1722780 Clean_trans.eQTLp.txt
# 1319533 Clean_trans.eQTL.txt

#(1180857 + 3540659 + 10235805 + 1208608 + 1722780 + 1319533)/1e6


