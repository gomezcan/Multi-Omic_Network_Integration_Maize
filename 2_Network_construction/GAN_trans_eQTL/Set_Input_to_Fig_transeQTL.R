#library(tidygraph)
library(tidyverse)
#library(ggrepel)
#library(ggpubr)
library(data.table)
#library(GeneOverlap)
#library(viridis)
#library(ComplexHeatmap)
#library(fgsea)
#library(reshape)
library(scales)


ReplaceName <- function(ids){
  
  for (i in 1:nrow(Top45)){
    ids <- gsub(Top45$V1[i], Top45$V2[i], ids)
  }
  
  for (i in 1:nrow(TFdic)){
    ids <- gsub(TFdic$TF.v4[i], TFdic$TF.Name[i], ids)
  }
  return(ids)
}

subset(Net_teQTL, source=='Zm00001d037784')

######################################################################################
#########       Read different eQTL data sets and define summary tables      #########
######################################################################################
## syntenic genes 
Syntenic <- as_tibble(read.table("Data/Annotations/Zm.v4.synteny.genes.txt", h=T, stringsAsFactors = F, sep = '\t'))
Syntenic <- subset(Syntenic, gene_synteny=='syntenic')$gene_id

#
saf <- as_tibble(read.table("Data/eQTL_data/Zea_mays.B73_RefGen_v4.46.saf", stringsAsFactors = F))


# Total cis- & trans-eQTLs
eQTLs <- as_tibble(read.table("Data/eQTL_data/Final_cis_trans_all_eQTL_10012021.bed", stringsAsFactors = F))
colnames(eQTLs) <- c("chr",'s', "e","snp", "Target", "Support")
eQTLs[,"Index"] <- paste(eQTLs$chr, eQTLs$s, eQTLs$snp, sep = ":")
eQTLs <- eQTLs[,-c(1,2,3)]
eQTLs <- unique(eQTLs[,c(2,4)])

# trans-eQTL
trans.eQTL <- as_tibble(read.table("Data/eQTL_data/trans.eQTL_target.txt", stringsAsFactors = F))
colnames(trans.eQTL) <- c("chr",'s','e','snp','Target','Support','chr.source', 'tss.source', 'e.source', 'source',"score","strand","dis","class")
trans.eQTL[,"Index"] <- paste(trans.eQTL$chr, trans.eQTL$s, trans.eQTL$snp, sep = ":")
trans.eQTL <- unique(trans.eQTL[,-c(1,2,3,4,9,11,13)])

# transcis-eQTL definition 
Get_transcis <- function(trans){
  
  #trans <- trans[1:10,]
  trans[,"tem"] <- paste(trans$Target, trans$source, trans$Index, sep = "_")
  out <- trans
  
  #
  trans <- trans[,c("Target","source","Index", "tem")]
  #
  trans <- left_join(trans, saf, by=c("Target"="V1"))
  trans <- left_join(trans, saf, by=c("source"="V1"))
  #
  trans <- subset(trans, V2.x == V2.y) # same chr
  #
  trans[,"DisS_T"] <- abs(trans$V3.x - trans$V3.y)
  trans <- unique(trans[,c("tem", "DisS_T")])
  print(trans)
  #
  out <- subset(out, tem %in% trans$tem)
  
  out <- left_join(out, trans, by='tem')
  out <- out[, (colnames(out) != "tem")]
  out$class <- "transciseQTL"
  out <- unique(out)
  
  return(out)
}

transcis_eQTL <- Get_transcis(trans.eQTL)
transcis_eQTL <- subset(transcis_eQTL, abs(DisS_T) <= 50000) # trans-cis: defined as those close proximity to target genes

# Redefine trans-eQTL
trans.eQTL[,"tem"] <- paste(trans.eQTL$Target, trans.eQTL$source, trans.eQTL$Index, sep = "_")
transcis_eQTL[,"tem"] <- paste(transcis_eQTL$Target, transcis_eQTL$source, transcis_eQTL$Index, sep = "_")
trans.eQTL <- subset(trans.eQTL, !(tem %in% transcis_eQTL$tem))

# trans-eQTLp
trans.eQTLp <- as_tibble(read.table("Data/eQTL_data/trans.eQTLp_target.txt", stringsAsFactors = F))
colnames(trans.eQTLp) <- c("chr",'s','e','snp','Target','Support','chr.source', 'tss.source', 'e.source', 'source',"score","strand","dis","class")
trans.eQTLp[,"Index"] <- paste(trans.eQTLp$chr, trans.eQTLp$s, trans.eQTLp$snp, sep = ":")
trans.eQTLp <- unique(trans.eQTLp[,-c(1,2,3,4,9,11,13)])
# define target_source_snp index
trans.eQTLp[,"tem"] <- paste(trans.eQTLp$Target, trans.eQTLp$source, trans.eQTLp$Index, sep = "_")

# cis-eQTLt
cis.eQTLt <- as_tibble(read.table("Data/eQTL_data/cis.eQTLt_target.txt", stringsAsFactors = F))
colnames(cis.eQTLt) <- c("chr",'s','e','snp','Target','Support','chr.source', 'tss.source', 'e.source', 'source',"score","strand","dis","class")
cis.eQTLt[,"Index"] <- paste(cis.eQTLt$chr, cis.eQTLt$s, cis.eQTLt$snp, sep = ":")
cis.eQTLt <- unique(cis.eQTLt[,-c(1,2,3,4,9,11,13)])
cis.eQTLt[,"tem"] <- paste(cis.eQTLt$Target, cis.eQTLt$source, cis.eQTLt$Index, sep = "_")

# cis-eQTL
cis.eQTL <- as_tibble(read.table("Data/eQTL_data/ciseQTL_noFiter.txt", stringsAsFactors = F))
colnames(cis.eQTL) <- c("chr",'s','e','snp','Target','Support','chr.source', 'tss.source', 'e.source', 'source',"score","strand","dis","class")
cis.eQTL[,"Index"] <- paste(cis.eQTL$chr, cis.eQTL$s, cis.eQTL$snp, sep = ":")
cis.eQTL <- subset(cis.eQTL, abs(dis) <= 50000)
cis.eQTL <- unique(cis.eQTL[,-c(1,2,3,4,9,11)])
cis.eQTL[,"tem"] <- paste(cis.eQTL$Target, cis.eQTL$source, cis.eQTL$Index, sep = "_")


#
trans.eQTL <- unique(trans.eQTL[, -c(2)])   # remove Support
trans.eQTLp <- unique(trans.eQTLp[, -c(2)]) # remove Support 
transcis_eQTL <- unique(transcis_eQTL[, -c(2)]) # remove Support 
cis.eQTL <- unique(cis.eQTL[, -c(2)]) # remove Support 
cis.eQTLt <- unique(cis.eQTLt[, -c(2)]) # remove Support 

# Confirm unique values
table(trans.eQTL$tem %in% cis.eQTL$tem)      # trans in cis target?
table(cis.eQTL$tem %in% cis.eQTLt$tem)       # cis in cis target?
table(cis.eQTL$tem %in% trans.eQTLp$tem)     # cis in trans promoter?
table(cis.eQTLt$tem %in% trans.eQTL$tem)     # cis target in trans?
table(trans.eQTLp$tem %in% cis.eQTL$tem)     # trans promoter in cis?
table(cis.eQTL$tem %in% cis.eQTLt$tem)       # cis in cis target?
table(cis.eQTLt$tem %in% cis.eQTL$tem)       # cis in cis target?

#cis.eQTL <- subset(cis.eQTL, !(tem %in% cis.eQTLt$tem))
#subset(trans.eQTL, Index %in% cis.eQTLt$Index)
#subset(cis.eQTLt, Index %in% trans.eQTL$Index)

#table(trans.eQTL$Index %in% cis.eQTLt$Index)
#table(trans.eQTL$tem %in% cis.eQTLt$tem)

# Count eQTL not assigned
eQTLSyn <- subset(eQTLs, Target %in% Syntenic)

#dim(eQTLSyn)
#length(unique(eQTLSyn$Index))
#length(unique(eQTLSyn$Target))

# Define eQTL eQTL-target unassigned

# eQTL in trans-eQTL
eQTLSyn <- eQTLSyn[!(paste(eQTLSyn$Target, eQTLSyn$Index, sep = "_") 
                     %in% paste(trans.eQTL$Target, trans.eQTL$Index, sep = "_")),]
# eQTL in trans-eQTLp
eQTLSyn <- eQTLSyn[!(paste(eQTLSyn$Target, eQTLSyn$Index, sep = "_") %in% paste(trans.eQTLp$Target, trans.eQTLp$Index, sep = "_")),]
# eQTL in trans-eQTLp
eQTLSyn <- eQTLSyn[!(paste(eQTLSyn$Target, eQTLSyn$Index, sep = "_") %in% paste(transcis_eQTL$Target, transcis_eQTL$Index, sep = "_")),]

# eQTL in cis-eQTL
eQTLSyn <- eQTLSyn[!(paste(eQTLSyn$Target, eQTLSyn$Index, sep = "_") %in% paste(cis.eQTL$Target, cis.eQTL$Index, sep = "_")),]
# eQTL in cis-eQTLt
eQTLSyn <- eQTLSyn[!(paste(eQTLSyn$Target, eQTLSyn$Index, sep = "_") %in% paste(cis.eQTLt$Target, cis.eQTLt$Index, sep = "_")),]

write.table(eQTLSyn, "Clean_eQTL_Syntenic.txt", sep = "\t", row.names = F, quote = F)
write.table(trans.eQTL, "Clean_trans.eQTL.txt", sep = "\t", row.names = F, quote = F)
write.table(trans.eQTLp, "Clean_trans.eQTLp.txt", sep = "\t", row.names = F, quote = F)
write.table(transcis_eQTL, "Clean_transcis_eQTL.txt", sep = "\t", row.names = F, quote = F)
write.table(cis.eQTL, "Clean_cis.eQTL.txt", sep = "\t", row.names = F, quote = F)
write.table(cis.eQTLt, "Clean_cis.eQTLt.txt", sep = "\t", row.names = F, quote = F)



# 1180857 Clean_cis.eQTLt.txt
# 3540659 Clean_cis.eQTL.txt
# 10235805 Clean_eQTL_Syntenic.txt
# 1208608 Clean_transcis_eQTL.txt
# 1722780 Clean_trans.eQTLp.txt
# 1319533 Clean_trans.eQTL.txt

#(1180857 + 3540659 + 10235805 + 1208608 + 1722780 + 1319533)/1e6


