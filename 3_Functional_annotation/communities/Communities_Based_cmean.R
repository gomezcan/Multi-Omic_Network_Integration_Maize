library(RSpectra)
library(reshape2)
library(scales)
library(tidyverse)  # data manipulation
library(ggpubr)
library(data.table)
library(GeneOverlap)
library(viridis)
library(ComplexHeatmap)
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization
library(RColorBrewer)
library(ggrepel)
library(corrplot)
library(ppclust)
library(cluster)
library(fclust)

set.seed(123)

########################################
########       Functions        ########
########################################

Get_Eigenvalue <- function(squareM){
  
  # cluster by Row
  Targscormatrix <- cor(t(squareM), method = 'spearman')
  print(dim(Targscormatrix))

  # Eigenvalue
  eg_val <- eigs_sym(Targscormatrix, 500, opts = list(retvec = FALSE))$values
  
  #write.table(out, paste0(class,"_eigenvalue_decomposition.txt"), sep = '\t', row.names = F, quote = F)
  return(eg_val)
}

Get_Eigenvalue_cor_abj <- function(net){
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
  
  # cluster by Targets
  Targscormatrix <- cor(t(net), method = 'spearman')
  
  eigs_cor <- eigs_sym(Targscormatrix, 500, opts = list(retvec = FALSE))$values
  
  return(list(eigs=eigs_cor, corM=Targscormatrix, Adj=net))
}

Get_Eigenvalue_tsne <- function(net){
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
  
  print(head(tsne_table))
  # cluster by Targets
  Targscormatrix <- cor(t(tsne_table), method = 'spearman')
  
  eigs_tsne_cor <- eigs_sym(Targscormatrix, 500, opts = list(retvec = FALSE))$values
  
  return(list(tsne=tsne_table, eigs=eigs_tsne_cor))
}

########################################

########################################
########          Data          ########
########################################

# TF and CoReg
TF_CoR <- as_tibble(read.table("Data/Annotations/TF_CoR_Mazie.txt", h=T, stringsAsFactors = F)) 

# Read pecanpy results
PDIpecanpy <- as.data.frame(fread("../Fig_pecanpy/PDI.pecanpy.W20.txt", header = F))
row.names(PDIpecanpy) <- PDIpecanpy$V1
PDIpecanpy <- PDIpecanpy[,-c(1)]

# Read PDIs
PDI <- as_tibble(fread("../Fig_pecanpy/PDI.pecanpy.input.txt", header = F))
PDI <- unique(PDI)

unique(PDI$V1)

########################################


## 1 Eigen values with pecanpy 
EgV_pecanpy <- Get_Eigenvalue(PDIpecanpy)

## 2 Eigen values with cor adjacence matrix
# First 500 Eigen values
EgV_AdjM <- Get_Eigenvalue_cor_abj(PDI)
EgV_AdjM$eigs

EgV_tsne_AdjM <- Get_Eigenvalue_tsne(PDI)
EgV_tsne_AdjM$eigs

EgV_tsne_AdjM$tsne

EgV_results <- tibble(x=seq_along(EgV_pecanpy), 
                      EgVp=EgV_pecanpy, 
                      EgVA=EgV_AdjM$eigs)
EgV_results


# Plot
Plor_EgV_results <- ggplot(EgV_results, aes(x=x)) +
  geom_line(aes(y=EgVp), size=2, color='#69b3a2') +
  geom_line(aes(y=EgVA), size=2) +
  scale_y_continuous(expand = c(0,0), limits = c(0,150)) + 
  scale_x_continuous(expand = c(0,0), limits = c(0,250)) + 
  xlab("Component") + ylab("Eigenvalue") +
  geom_hline(yintercept = 1, color='red', linetype="dashed") + 
  #
  annotate("segment", x=length(EgV_pecanpy[EgV_pecanpy >=1])+1,
           xend=length(EgV_pecanpy[EgV_pecanpy >=1]), y=0, yend=50,color='#69b3a2', arrow=arrow()) +
  annotate("segment", x=length(EgV_AdjM$eigs[EgV_AdjM$eigs >=1])+1,
           xend=length(EgV_AdjM$eig[EgV_AdjM$eigs >=1]), y=0, yend=70,color='black', arrow=arrow()) +
  #
  annotate("text", label=paste0("Pecanpy (egv >= 1)\n", length(EgV_pecanpy[EgV_pecanpy >=1]), " components"), 
           color ='#69b3a2', size=4,
           x = length(EgV_pecanpy[EgV_pecanpy >=1]), y = 60) +
  annotate("text", label=paste0("Abj. M (egv >= 1)\n", length(EgV_AdjM$eigs[EgV_AdjM$eigs >=1]), " components"), 
           color ='black', size=4,
           x = length(EgV_AdjM$eigs[EgV_AdjM$eigs >=1]), y = 80) +
  #
  annotate("text", label="Kaiser criterion (eigenvalue > 1)", color='red', size=5, x = 50, y = 10) +
  theme_pubclean()

Plor_EgV_results   <- ggpar(Plor_EgV_results, font.tickslab = 14, font.xtickslab = 14, font.ytickslab = 14)
Plor_EgV_results  

# Fuzzy c-Means

cm <- cmeans(PDIpecanpy, 127)

hist(cm$membership[,2])

dim(as.data.frame(cm$cluster))


table(cm$membership)

dim(cm$membership)

cm_adj <- cmeans(g, 139)

PDIpecanpy[1:5,1:5]

table(cm$cluster)

cm_membership <- as.data.frame(cm$membership)
head(cm_membership)

hist(cm$membership[,1], 100)

apply(cm$membership, 1, max)



M <- as.matrix(PDIpecanpy[1:100,])
M <- as.matrix(PDIpecanpy[row.names(PDIpecanpy) %in% TF_CoR$GeneID,])
#
fviz_cluster(list(data = M, cluster=cm$cluster), 
             #ellipse.type = "norm",
             #ellipse.level = 0.68,
             palette = "jco",
             ggtheme = theme_minimal(), 
             labelsize = 1
            )


Heatmap(M,
        #km=50,
        show_row_names = F,
        show_column_names = F, 
        heatmap_width  = unit(10, 'cm'))
