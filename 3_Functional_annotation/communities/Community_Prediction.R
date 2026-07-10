library(reshape2)
library(scales)
library(tidyverse)  # data manipulation
library(ggpubr)
library(data.table)
library(GeneOverlap)
library(viridis)
library(Rtsne)
library(RColorBrewer)
library(ggrepel)
library(igraph)
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization

### Methods 
# 1. Latent variable modeling (eigenvalue decomposition)
# 2. Spinglass algorithm
# 3. Walktrap algorithm

##################################################
##########          Functions       ##############
##################################################


Get_Eigenvalue <- function(net, class){
  # decast to get adjacency network 
  
  colnames(net) <- c("Source", "Target", "Val")
  colnames(PDI) <- c("Source", "Target", "Val")
  print(dim(net))
  net <- reshape2::dcast(PDI, Target ~ Source, value.var = "Val")
  #net <- reshape2::dcast(net, Target ~ Source, value.var = "Val")
  #
  row.names(net) <- net$Target
  net <- net[,-c(1)]
  #
  net[is.na(net)] <- 0
  print(dim(net))
  
  # cluster by TFs
  TFscormatrix <- cor(net, method = 'spearman')
  
  # cluster by Targets
  Targscormatrix <- cor(t(net), method = 'spearman')
  
  #net[1:5,1:5]
  # Eigenvalue
  egTFs <- tibble(eigenClass="tfs", values=eigen(TFscormatrix)$values)
  
  # Eigenvalue
  egtgs <- tibble(eigenClass="tgs", values=eigen(Targscormatrix)$values)
  
  out <- rbind(egTFs, egtgs)
  
  write.table(out, paste0(class,"_eigenvalue_decomposition.txt"), sep = '\t', row.names = F, quote = F)
}



##################################################
##########          Data       ##############
##################################################

PDI <- as_tibble(fread("../Fig_pecanpy/PDI.pecanpy.input.txt", header = F))
PDI[,"val"] <- 1

PDI$V1 <- gsub("pChIP.", "",PDI$V1)
PDI$V1 <- gsub("DAP.", "",PDI$V1)
PDI$V1 <- gsub("ChIP.", "",PDI$V1)
PDI <- unique(PDI)


write.table(PDI, '../Fig_pecanpy/PDI.pecanpy.input.txt', sep='\t', row.names = F, quote = F)

Get_Eigenvalue(PDI, "PDI")



##
g <- graph_from_data_frame(PDI[,1:2], directed = F)
rw <- cluster_walktrap(g, steps = 10)
length(unique(rw$membership))
length(unique(rw$modularity))
plot(rw, g)

rw$membership

table(rw$membership)
hist(rw$membership)

rw$names



####################################################################################
# library(qgraph)
# library(CliquePercolation)
# 
# W <- qgraph::qgraph(PDI[sample(1:nrow(PDI), 10000), 1:2], DoNotPlot=T)
# W <- qgraph::qgraph(PDI[, 1:2], DoNotPlot=T)
# 
# thresholds <- cpThreshold(W,
#                           method = "unweighted", k.range = c(3,4),
#                           threshold = c("largest.components.ratio","chi"))
#   
# 
# cpk4 <- cpAlgorithm(W, k = 3, method = "unweighted")
# summary(cpk4)
####################################################################################