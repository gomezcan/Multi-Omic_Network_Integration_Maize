library(ggmap)
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
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization

################################################
######            First round             ######
################################################
##

PDI <- as_tibble(fread("PDI.pecanoy.input.txt", header = F))
PDI[,"val"] <- 1

Get_tsne <- function(net){
  # decast to get adjacency network 
  
  colnames(net) <- c("Source", "Target", "Val")
  print(dim(net))
  net <- reshape2::dcast(net, Target ~ Source, value.var = "Val")
  #
  row.names(net) <- net$Target
  net <- net[,-c(1)]
  #
  net[is.na(net)] <- 0
  print(dim(net))
  
  #net[1:5,1:5]
  tsne_table = Rtsne(as.matrix(net), check_duplicates=FALSE, pca=TRUE, 
                     perplexity=30, theta=0.2, dims=3)
  #
  
  tsne_table = as.data.frame(tsne_table$Y)
  print(dim(tsne_table))
  row.names(tsne_table) <- row.names(net)
  
  return(list(tsne=tsne_table, A=net))
}

EgV_tsne_AdjM <- Get_Eigenvalue_tsne(PDI)

PDI_tsne <- EgV_tsne_AdjM$tsne

fviz_nbclust(EgV_tsne_AdjM$tsne[,1:3], kmeans, nboot = 10, k.max = 100) # best opt: 58 k

## Creating k-means clustering model, and assigning the result to the data used to create the tsne
fit_cluster_kmeans = kmeans(PDI_tsne, 139) # 

PDI_tsne$km <-  factor(fit_cluster_kmeans$cluster)

PDI_tsne[,"GeneID"] <- row.names(PDI_tsne)
PDI_tsne <- as_tibble(PDI_tsne)


## Exploratory plot
# labels
clusterLabels <- PDI_tsne %>% 
  group_by(km) %>% select(V1, V2) %>% summarize_all(mean)


# Plot function
plot_cluster <- function(data, var_cluster, palette, labels) {
  #
  #Genelabel <- subset(data, GeneId %in% GenesTarget)
  #clusters_target <- unique(as.character(subset(data, GeneId %in% GenesTarget)$km))
  #
  #print(Genelabel)
  #
  ggplot(data, aes_string(x="V1", y="V2", color="km")) +
    geom_point(size=0.05) +
    #geom_point(data=subset(data, km %in% clusters_target), size=0.1, aes_string(x="V1", y="V2", color="km")) + 
    #geom_label_repel(data=Genelabel, aes_string(x="V1", y="V2", label="GeneId", color="km"),  
    #                 max.overlaps = Inf, box.padding = 0.5) +
    
    guides(colour=guide_legend(override.aes=list(size=3))) +
    xlab("tsne 1") + ylab("tsne 2") +
    ggtitle("") +
    geom_label_repel(data=labels, aes(label = km)) +
    theme_nothing()+
    guides(colour = 'none') 
}

# Plot
plot_k <- plot_cluster(PDI_tsne, "km", 139, clusterLabels)
plot_k

