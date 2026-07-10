####################################################
########           adjacency                ########
####################################################


library(Rtsne)

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
  TF_target_DF_M <- apply(TF_target_DF, 2, scale) # scale by targets: Columns on this format
  row.names(TF_target_DF_M) <- row.names(TF_target_DF)
  print(dim(TF_target_DF_M) == dim(TF_target_DF))
  
  return(list(Freq=TF_target_DF, Z=TF_target_DF_M))
  
}

PDI_adj <- Matrix_freq(PDI)
CoExp_adj <- Matrix_freq(CoExp)
teQTL_adj <- Matrix_freq(teQTL)

All_adj <- Matrix_freq(teQTL)


tsne_PDI = Rtsne(as.matrix(t(PDI_adj$Z)), 
                 check_duplicates=FALSE, pca=TRUE, 
                 perplexity=30, theta=0.2, dims=2,
                 num_threads=40)

tsne_CoExp = Rtsne(as.matrix(t(CoExp_adj$Z)), 
                   check_duplicates=FALSE, pca=TRUE, 
                   perplexity=30, theta=0.2, dims=2,
                   num_threads=40)

tsne_teQTL = Rtsne(as.matrix(t(teQTL_adj$Z)), 
                   check_duplicates=FALSE, pca=TRUE, 
                   perplexity=30, theta=0.2, dims=2,
                   num_threads=40)


tsne_PDI = as.data.frame(tsne_PDI$Y)
row.names(tsne_PDI) <- colnames(PDI_adj$Z)

tsne_CoExp = as.data.frame(tsne_CoExp$Y)
row.names(tsne_CoExp) <- colnames(CoExp_adj$Z)

tsne_teQTL = as.data.frame(tsne_teQTL$Y)
row.names(tsne_teQTL) <- colnames(teQTL_adj$Z)

colnames(tsne_PDI) <- paste0("GRN_", colnames(tsne_PDI))
colnames(tsne_CoExp) <- paste0("CEN_", colnames(tsne_CoExp))
colnames(tsne_teQTL) <- paste0("GAN_", colnames(tsne_teQTL))


tsne_All <- as_tibble(data.frame(gid=Reduce(intersect, list(row.names(tsne_PDI), row.names(tsne_CoExp), row.names(tsne_teQTL)))))

tsne_All <- left_join(tsne_All, tsne_PDI %>% mutate(gid = rownames(tsne_PDI)), by="gid")
tsne_All <- left_join(tsne_All, tsne_CoExp %>% mutate(gid = rownames(tsne_CoExp)), by="gid")
tsne_All <- left_join(tsne_All, tsne_teQTL %>% mutate(gid = rownames(tsne_teQTL)), by="gid")

tsne_All[,2:7] <- apply(tsne_All[,2:7], 2, scale)

write.table(tsne_All, "tsne_All.txt", sep = "\t", quote = F, row.names = F)

row.names(peAll) <- peAll$gid
