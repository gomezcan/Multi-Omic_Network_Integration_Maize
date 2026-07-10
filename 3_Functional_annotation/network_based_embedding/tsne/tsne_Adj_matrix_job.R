#install.packages("clValid")
library(Rtsne)
library(tidyverse)
library(data.table)

#

##################################################
##########          Functions        #############
##################################################

ReplaceNamePWY <- function(ids){
  
  for (i in 1:nrow(CornCYC)){
    w <- paste0('\\<', CornCYC$Pathway.id[i], '\\>')
    ids <- gsub(w, CornCYC$Pathway.name[i], ids)
    #ids <- gsub("_", " ", ids)
  }
  return(ids)
}
#
ReplaceName <- function(ids){
  # for (i in 1:nrow(Top45)){
  #   ids <- gsub(Top45$V1[i], Top45$V2[i], ids)
  # }
  
  for (i in 1:nrow(TFdic)){
    ids <- gsub(TFdic$TF.v4[i], TFdic$TF.Name[i], ids)
  }
  return(ids)
}

##################################################
##########        Annotations       ##############
##################################################

# PDI
PDI <- as_tibble(unique(fread("../Fig_PDI/PDI_NetworkFinal.10_11_2021.txt")[,c(2,3)]))
colnames(PDI)[1] <- "Source"

# CoExp
CoExp <- as_tibble(unique(fread("../Fig_Coexpression/CoExp_NetworkFinal.10_11_2021.txt")))
colnames(CoExp)[2] <- "Source"
CoExp <- unique(CoExp[,2:3])

# teQTL
teQTL <- as_tibble(unique(fread("../Fig_transeQTL/teQTL_NetworkFinal.10_11_2021.txt")))
colnames(teQTL)[1] <- "Source"


Matrix_freq <- function(TF_target_DF) {
  
  # Count TF-target associations
  TF_target_DF <- as_tibble(as.data.frame(table(TF_target_DF), stringsAsFactors = F))
  
  # filter 
  TF_target_DF <- subset(TF_target_DF, Freq > 0) # remove TF-target with zero freq to reduce speed during dcast
  
  TF_target_DF <- reshape2::dcast(TF_target_DF, Source ~ Target)
  row.names(TF_target_DF) <- TF_target_DF$Source
  TF_target_DF <- TF_target_DF[,-c(1)]
  # 
  TF_target_DF[is.na(TF_target_DF)] <- 0
  
  # scale values
  TF_target_DF_M <- t(apply(TF_target_DF, 1, scale)) # scale by TF, which are the features: rows on this format
  row.names(TF_target_DF_M) <- row.names(TF_target_DF)
  colnames(TF_target_DF_M) <- colnames(TF_target_DF)
  print(dim(TF_target_DF_M) == dim(TF_target_DF))
  
  return(list(Freq=TF_target_DF, Z=TF_target_DF_M))
  
}




# Combing layers
All_adj <- Matrix_freq(rbind(teQTL, CoExp, PDI))

CEN_GRN_adj <- Matrix_freq(rbind(CoExp, PDI))


write.table(t(All_adj$Freq), "Tsne_Files/AdjM_all.txt", quote = F, row.names = F, sep = '\t')
write.table(t(All_adj$Z), "Tsne_Files/AdjM_Scaled_all.txt", quote = F, row.names = F, sep = '\t')

write.table(t(CEN_GRN_adj$Freq), "Tsne_Files/AdjM_CEN_GRN.txt", quote = F, row.names = F, sep = '\t')
write.table(t(CEN_GRN_adj$Z), "Tsne_Files/AdjM_Scaled_CEN_GRN.txt", quote = F, row.names = F, sep = '\t')

# tsene with  all networks
tsne_All = Rtsne(as.matrix(t(All_adj$Z)), check_duplicates=FALSE, pca=TRUE, perplexity=30, theta=0.2, dims=2, num_threads=50)
tsne_All = as_tibble(as.data.frame(tsne_All$Y))
tsne_All["gid"] <- colnames(All_adj$Z)

write.table(tsne_All, "Tsne_Files/tsne_All.txt", quote = F, row.names = F, sep = '\t')

# tsene with  without GRN network
tsne_CEN_GRN = Rtsne(as.matrix(t(CEN_GRN_adj$Z)), check_duplicates=FALSE, pca=TRUE, perplexity=30, theta=0.2, dims=2, num_threads=50)
tsne_CEN_GRN = as_tibble(as.data.frame(tsne_CEN_GRN$Y))
tsne_CEN_GRN["gid"] <- colnames(CEN_GRN_adj$Z)

write.table(tsne_CEN_GRN, "Tsne_Files/tsne_CEN_GRN.txt", quote = F, row.names = F, sep = '\t')

# 
# PDI_adj <- Matrix_freq(PDI)
# CoExp_adj <- Matrix_freq(CoExp)
# teQTL_adj <- Matrix_freq(teQTL)
# 
# 
# tsne_PDI = Rtsne(as.matrix(t(PDI_adj$Z)),
#                  check_duplicates=FALSE, pca=TRUE,
#                  perplexity=30, theta=0.2, dims=2,
#                  num_threads=40)
# 
# tsne_CoExp = Rtsne(as.matrix(t(CoExp_adj$Z)),
#                    check_duplicates=FALSE, pca=TRUE,
#                    perplexity=30, theta=0.2, dims=2,
#                    num_threads=40)
# 
# tsne_teQTL = Rtsne(as.matrix(t(teQTL_adj$Z)),
#                    check_duplicates=FALSE, pca=TRUE,
#                    perplexity=30, theta=0.2, dims=2,
#                    num_threads=40)
# 
# 
# 
