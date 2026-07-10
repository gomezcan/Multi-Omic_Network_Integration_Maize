library(scales)
library(tidyverse)
library(ggrepel)
library(ggpubr)
library(data.table)
library(viridis)
library(reshape2)
library(circlize)
library(purrr)
library(gplots)
library(ggplot2)
library(parallel)
library(org.Zmays.eg.db)
library(patchwork)
library(ggthemes)
library(edgeR)

##################################################
#####                 Functions             ######
##################################################

#
ReplaceName <- function(ids){
  # for (i in 1:nrow(Top45)){
  #   ids <- gsub(Top45$V1[i], Top45$V2[i], ids)
  # }
  
  for (i in 1:nrow(TFdic)){
    ids <- gsub(TFdic$V2[i], TFdic$V1[i], ids)
  }
  return(ids)
}

vennfuncInt <- function(list_x){
  venn <- Venn(list_x)
  data <- process_data(venn)
  
  colorGroups <- c(CEN = '#FF9933', GRN='#1E90FF', GAN='#FFD700', eGRN='#FF1493')
  colfunc <- colorRampPalette(colorGroups)
  col <- colfunc(4)
  
  colorGroups <- c(CEN="gray100",GRN="gray99", GAN="gray98", eGRN='gray97')
  colfunc2 <- colorRampPalette(colorGroups)
  col2 <- colfunc2(15)
  
  ggplot() +
    geom_sf(aes(fill=name), data = venn_region(data), show.legend = F) +
    geom_sf(aes(color=name), size = 1.5, data = venn_setedge(data), show.legend = F) +
    #
    geom_sf_text(aes(label = name), size=6, data = venn_setlabel(data)) +
    geom_sf_text(aes(label= scales::comma(count, accuracy = 1)), size=5, data = venn_region(data)) +
    #
    scale_fill_manual(values = col2) + # 
    scale_color_manual(values = alpha(col, .5)) +
    #
    theme_void() +
    theme(plot.margin = unit(c(0.5,0.1, 0.5, 0.1), "cm"),
          text = element_text(family="Helvetica")) # +
  #xlim(-150,1000)
}


## chop a string by a separator and return specified field
chop=function(myStr,mySep,myField){
  
  choppedString=sapply(strsplit(myStr,mySep),"[",myField)
  if(length(myField)>1){
    choppedString=apply(choppedString,2,function(x){paste0(x[!is.na(x)],collapse=mySep)})
  }
  return(choppedString)
}

TMM <- function(Counts){
  ### Calculate the normalization factor for every sample
  factors_norm <- calcNormFactors(Counts, lib_size = library_size, method = 'TMM')
  ### Effective size
  effective_size <- factors_norm*as.numeric(library_size)/10**6
  ### Normalize the counts
  counts_norm <- sweep(Counts, 2, effective_size, FUN = '/')
  return(counts_norm)
}

calc_tau <- function(m,  byRow = TRUE) {
  #' Calculate Tau scores for each gene
  #' 
  #' @param m Matrix of expression values
  #' @param byRow if TRUE, treats genes as row values and samples as columns.
  #' 
  #' @return a vector of Tau scores
  #' 
  if(!byRow){
    m <- t(m)
  }
  
  row_maxes <- apply(m, 1, max)
  
  m <- m / row_maxes
  tau <- Matrix::rowSums(1 - m) / (ncol(m) - 1)
  tau[is.na(tau)] <- 0
  return(tau)
}

##################################################

##################################################
##########        Annotations       ##############
##################################################

#saf <- as_tibble(read.table("Data/eQTL_data/Zea_mays.B73_RefGen_v4.46.saf", stringsAsFactors = F))
#saf1 <- subset(saf, V5=="+")[,c(1,2,3)]
#saf2 <- subset(saf, V5=="-")[,c(1,2,4)]
#colnames(saf1) <- c("GeneID", "chrAnn", "TSS")
#colnames(saf2) <- c("GeneID", "chrAnn", "TSS")
#
#saf <- rbind(saf1, saf2)

# PDI
GRN <- unique(fread("../Fig_PDI/Only_PDI_NetworkFinal.10_14_2022.txt")[,c(2,3)])
colnames(GRN)[1] <- "Source"
dim(PDI)

# PDI eQTL
eGRN <- unique(fread("../Fig_PDI/CisE_PDI_NetworkFinal.10_14_2022.txt")[,c(2,3)])
colnames(eGRN)[1] <- "Source"

GRN <- rbind(GRN, eGRN)

# CoExp
CEN <- unique(fread("../Fig_Coexpression/CoExp_NetworkFinal.10_11_2021.txt"))
colnames(CEN)[2] <- "Source"
CEN <- unique(CEN[,2:3])

# teQTL
GAN <- unique(fread("../Fig_transeQTL/teQTL_NetworkFinal.10_11_2021.txt"))
colnames(GAN)[1] <- "Source"

# Genes in synteny
GenesMaize <- unique(as.character(read.table("Data/Annotations/Zm.v4.synteny.genes.txt", h=T)$gene_id))
length(GenesMaize)

# TF names
TFdic <- as_tibble(read.table("Data/Annotations/TF_Id_Name.txt", h=F, stringsAsFactors = F))

# Total - annotations
annotations <- as.data.table(read.table("../Fig_MethodsComparison/Summary.Total.Annotation.txt", 
                                    h=T, stringsAsFactors = F, sep='\t'))

TF_annotated <- unique(annotations$TF)
length(TF_annotated)
##################################################

##################################################
######        Common interactions       ##########
##################################################

TFsinDEGs <- c("Zm00001d033859","Zm00001d020430","Zm00001d037317",
               "Zm00001d018971","Zm00001d033673","Zm00001d016052",
               "Zm00001d016255","Zm00001d026094","Zm00001d016520",
               "Zm00001d053369","Zm00001d038843","Zm00001d027757",
               "Zm00001d029963","Zm00001d032923")

ReplaceName(TFsinDEGs)


##################################################

##################################################
######        Targets distribution      ##########
##################################################

Total_freq <- as.data.table(table(rbind(CEN, GRN, eGRN, GAN)$Source))
Total_freq[,"Z"] <- scale(Total_freq$N, center = T, scale = T)

CEN_free <- as.data.table(table(CEN$Source)) %>%
  mutate(Z=scale(N, center = T, scale = T), Network='CEN')

GAN_free <- as.data.table(table(GAN$Source)) %>%
  mutate(Z=scale(N, center = T, scale = T), Network='GAN')

GRN_free <- as.data.table(table(GRN$Source)) %>%
  mutate(Z=scale(N, center = T, scale = T), Network='GRN')

eGRN_free <- as.data.table(table(eGRN$Source)) %>%
  mutate(Z=scale(N, center = T, scale = T), Network='eGRN')


zTotal_freq <- Total_freq[Total_freq$V1 %in% TF_annotated,]
hist(Total_freq$N)

zTotal_freq_DEGs <- zTotal_freq[zTotal_freq$V1 %in% TFsinDEGs,]
zTotal_freq_DEGs[,'TFname'] <- ReplaceName(zTotal_freq_DEGs$V1)


SamplingTargets <- function(ntf){
  scores <- mean(sample(zTotal_freq$Z, ntf))
  return(scores)
}

background <- tibble(z=unlist(lapply(c(1:1000), function(x) SamplingTargets(length(TFsinDEGs)))))

##################################################

##################################################
######      Summary TF by networks      ##########
##################################################

MethodsFreq <- unique(rbind(eGRN_free,
                           GRN_free,
                           CEN_free, 
                           GAN_free)[,c(1,4)]) %>%
  dplyr::filter(V1 %in% TFsinDEGs) %>%
  table %>%  as.data.table() %>%
  dplyr::mutate(TFname=ReplaceName(V1))

unique(MethodsFreq$TFname)

##################################################

##################################################
######      TF expression analysis      ##########
##################################################

ExpressionFiles <- list.files(path = 'Expression/', pattern = "*.tsv")
ExpressionFiles <- lapply(ExpressionFiles, function(x) fread(paste0("Expression/", x)))

# Count total CMPs: shout we replace by raw reads
library_size <- lapply(ExpressionFiles, function(x)  colSums(x[,-c(1)]))
library_size <- unique(unlist(library_size))


# Keep expresion for tfs with annotation
ExpressionDB <- ExpressionFiles[[1]]
ExpressionDB <- subset(ExpressionDB, gid %in% TF_annotated)



for (i in 2:length(ExpressionFiles)){
  tem <- subset(ExpressionFiles[[i]], gid %in% TF_annotated)
  ExpressionDB <- left_join(ExpressionDB, tem, by='gid')
}

## Remove duplicated  samples

# lists of TFs in expression
mask <- colnames(ExpressionDB)
mask <- mask[!(grepl('.y', mask))]


ids <- ExpressionDB$gid
ExpressionDB <- ExpressionDB[,-c(1)]
ExpressionDB <- as.matrix(ExpressionDB)[,(colnames(ExpressionDB) %in% mask)]
row.names(ExpressionDB) <- ids
#ExpressionDB[1:5,1:5]
#dim(ExpressionDB)

colnames(ExpressionDB) <- (gsub('.x', '', colnames(ExpressionDB)))

## Calculate TMM
DF_TMMs <- TMM(ExpressionDB)
#DF_TMMs[,"GeneID"] <- ids

DF_Tau <- calc_tau(DF_TMMs, byRow = TRUE) %>%
  as.data.frame() %>%
  tibble::rownames_to_column()
colnames(DF_Tau) <- c('TF', 'Tau')
DF_Tau <- as_tibble(DF_Tau)

# background distribution
SamplingTau <- function(ntf){
  scores <- mean(sample(DF_Tau$Tau, ntf))
  return(scores)
}

Tau_background <- tibble(Background=unlist(lapply(c(1:2000), 
                                                  function(x) SamplingTau(length(TFsinDEGs)))))
# 
# test significance of tau
DF_Tau_tested <- DF_Tau[DF_Tau$TF %in% TFsinDEGs,]

Tau_background

Zcore_background <- function(val){
  val <- (val - mean(Tau_background$Background)) / sd(Tau_background$Background)
  return(val)
}

# Z score from random background
DF_Tau_tested[,"Ztau"] <- Zcore_background(DF_Tau_tested$Tau)

# P values from Z score from random background
DF_Tau_tested[,"Pval_neg"] <- unlist(lapply(DF_Tau_tested$Ztau, function (x) pnorm(x, lower.tail=TRUE)))
DF_Tau_tested[,"Pval_pos"] <- unlist(lapply(DF_Tau_tested$Ztau, function (x) pnorm(x, lower.tail=FALSE)))

# if z negative keep left tail
DF_Tau_tested$Pval[DF_Tau_tested$Ztau < 0]  <- DF_Tau_tested$Pval_neg[DF_Tau_tested$Ztau < 0] 

# if z positive keep right tail
DF_Tau_tested$Pval[DF_Tau_tested$Ztau > 0]  <- DF_Tau_tested$Pval_pos[DF_Tau_tested$Ztau > 0] 

# clean up unused values
DF_Tau_tested <- DF_Tau_tested[,1:4]
DF_Tau_tested[, 'Name'] <- ReplaceName(DF_Tau_tested$TF)
DF_Tau_tested$Name <- paste0(DF_Tau_tested$Name, ' (',formatC(DF_Tau_tested$Pval, format = "e", digits = 1),
                             ')')


##################################################

##################################################
######              Plots               ##########
##################################################

######
# density plots
######

# Histogram with density plot: eGRN
ggplot(eGRN_free, aes(y=Z)) + 
  geom_histogram(aes(x=..density..), colour="black", fill="white") +
  geom_density(alpha=.5, fill="#C0C0C0") +
  geom_segment(data=subset(eGRN_free, V1 %in% TFsinDEGs),
               aes(x=0, xend=0.50, y=Z, yend=Z),
               linetype="dashed", 
               color = "red") +
  scale_x_continuous(expand = c(0,0), limits = c(0,0.8)) + 
  geom_text_repel(data=subset(eGRN_free, V1 %in% TFsinDEGs),
                  aes(x=0.50, y=Z, label=ReplaceName(V1)),
                  direction    = "y",
                  xlim = c(0.51, NA),
                  ylim = c(-1, NA),
                  vjust = 1, 
                  segment.size = 0.2,
                  max.iter = 1e4, 
                  max.time = 1,
                  force_pull   = 0,
                  size=2) + 
  theme_pubclean() +
  ylab('Z score (# targets)') + 
  xlab('Density') + 
  labs(subtitle='eGRN') + 
  theme(strip.text.y = element_text(size = 5, angle = 0), 
        #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text=element_text(size=10), 
        legend.position = 'bottom',
        text = element_text(size=10, family="Times")) -> Plot_eGRN

Plot_eGRN
# Histogram with density plot: GRN
ggplot(GRN_free, aes(y=Z)) + 
  geom_histogram(aes(x=..density..), colour="black", fill="white") +
  geom_density(alpha=.5, fill="#C0C0C0") +
  geom_segment(data=subset(GRN_free, V1 %in% TFsinDEGs),
               aes(x=0, xend=0.50, y=Z, yend=Z),
               linetype="dashed", 
               color = "red") +
  scale_x_continuous(expand = c(0,0), limits = c(0,0.8)) + 
  geom_text_repel(data=subset(GRN_free, V1 %in% TFsinDEGs),
                  aes(x=0.50, y=Z, label=ReplaceName(V1)),
                  direction    = "y",
                  xlim = c(0.51, NA),
                  ylim = c(-1, NA),
                  vjust = 1, 
                  segment.size = 0.2,
                  max.iter = 1e4, 
                  max.time = 1,
                  force_pull   = 0,
                  size=2) + 
  theme_pubclean() +
  ylab('Z score (# targets)') + 
  xlab('Density') + 
  labs(subtitle='GRN') + 
  theme(strip.text.y = element_text(size = 5, angle = 0), 
        #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text=element_text(size=10), 
        legend.position = 'bottom',
        text = element_text(size=10, family="Times")) -> Plot_GRN

# Histogram with density plot: GAN
ggplot(GAN_free, aes(y=Z)) + 
  geom_histogram(aes(x=..density..), colour="black", fill="white") +
  geom_density(alpha=.5, fill="#C0C0C0") +
  geom_segment(data=subset(GAN_free, V1 %in% TFsinDEGs),
               aes(x=0, xend=0.50, y=Z, yend=Z),
               linetype="dashed", 
               color = "red") +
  scale_x_continuous(expand = c(0,0), limits = c(0,0.8)) + 
  scale_y_continuous(expand = c(0,0), limits = c(-2, 3)) + 
  geom_text_repel(data=subset(GAN_free, V1 %in% TFsinDEGs),
                  aes(x=0.50, y=Z, label=ReplaceName(V1)),
                  direction    = "y",
                  xlim = c(0.51, NA),
                  ylim = c(-1, NA),
                  vjust = 1, max.overlaps=Inf,
                  segment.size = 0.2,
                  max.iter = 1e4, 
                  max.time = 1,
                  force_pull   = 0,
                  size=2) + 
  theme_pubclean() +
  ylab('Z score (# targets)') + 
  xlab('Density') + 
  labs(subtitle='GAN') + 
  theme(strip.text.y = element_text(size = 5, angle = 0), 
        #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text=element_text(size=10), 
        legend.position = 'bottom',
        text = element_text(size=10, family="Times")) -> Plot_GAN
Plot_GAN

# Histogram with density plot: GAN
ggplot(CEN_free, aes(y=Z)) + 
  geom_histogram(aes(x=..density..), colour="black", fill="white") +
  geom_density(alpha=.5, fill="#C0C0C0") +
  geom_segment(data=subset(CEN_free, V1 %in% TFsinDEGs),
               aes(x=0, xend=0.50, y=Z, yend=Z),
               linetype="dashed", 
               color = "red") +
  scale_x_continuous(expand = c(0,0), limits = c(0,0.8)) + 
  geom_text_repel(data=subset(CEN_free, V1 %in% TFsinDEGs),
                  aes(x=0.50, y=Z, label=ReplaceName(V1)),
                  direction    = "y",
                  xlim = c(0.51, NA),
                  ylim = c(-1, NA),
                  vjust = 1, max.overlaps=Inf,
                  segment.size = 0.2,
                  max.iter = 1e4, 
                  max.time = 1,
                  force_pull   = 0,
                  size=2) + 
  theme_pubclean() +
  ylab('Z score (# targets)') + 
  xlab('Density') + 
  labs(subtitle='CEN') + 
  theme(strip.text.y = element_text(size = 5, angle = 0), 
        #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text=element_text(size=10), 
        legend.position = 'bottom',
        text = element_text(size=10, family="Times")) -> Plot_CEN
######


######
# Summary method plots
######

Fig_9a <- {Plot_GAN + Plot_GRN + Plot_eGRN +Plot_CEN + plot_layout(nrow = 1)} 


MethodsFreq$TFname <- factor(MethodsFreq$TFname,
                             levels = c("TB1", "FEA4", "KN1", "MYBR32",
                                        "O2", "WRKY82", "HSF13","HSF18",
                                        "HSF20", "HSF29","RA1",
                                        "WRKY2","WRKY8","HSF24"))

MethodsFreq$N <- gsub('0', 'Without targ.', MethodsFreq$N)
MethodsFreq$N <- gsub('1', 'With targ.', MethodsFreq$N)

ggplot(MethodsFreq, aes(x=TFname, y=Network, fill=as.character(N))) +
  geom_tile(color="white", size=0.1) +
  #scale_fill_viridis(discrete = T, alpha = 0.5, name='In network') +
  scale_fill_manual(values = c("#CC99FF", "#E0E0E0"), name='Network')+
  #coord_equal() +
  labs(x=NULL, y=NULL) +
  theme_tufte(base_family="Helvetica") + 
  theme(axis.ticks=element_blank(),
        axis.text=element_text(size=7),
        axis.text.x =element_text(size=7, angle = 90, hjust=1),
        legend.position = 'bottom',) ->  Plot_b
Plot_b


Fig_9ab <- {Plot_GAN + Plot_GRN + Plot_eGRN + Plot_CEN + plot_layout(nrow = 1)} / Plot_b  +
  plot_layout(nrow = 2, heights = c(1,0.7))
  
pdf("Plots/Fig_S9.pdf", width=8, height=6)
print(Fig_9ab)
dev.off()


######

######
# tau distribution
#####

# Histogram with density plot: eGRN
ggplot(DF_Tau, aes(x=Tau)) + 
  geom_histogram(aes(y=..density..), colour="#FF1493", fill="white") +
  geom_density(alpha=.3, fill="#FF69B4") +
  geom_segment(data=subset(DF_Tau, TF %in% TFsinDEGs),
               aes(y=0, yend=2.50, x=Tau, xend=Tau),
               linetype="dashed", 
               color = "#8A2BE2") +
  scale_y_continuous(expand = c(0,0), limits = c(0,3)) + 
  scale_x_continuous(expand = c(0,0)) + 
  geom_text_repel(data=subset(DF_Tau, TF %in% TFsinDEGs),
                  aes(x=Tau, y=2.5, label=ReplaceName(TF)),
                  #direction    = "y",
                  xlim = c(NA, 0.95),
                  ylim = c(2.6, NA),
                  vjust = 1, 
                  segment.size = 0.2,
                  max.iter = 1e4, 
                  max.time = 1,
                  force_pull   = 0,
                  size=2) + 
  theme_pubclean() +
  ylab('Density') + 
  labs(subtitle='Tau distribution for TF annotated') + 
  theme(#strip.text.y = element_text(size = 5, angle = 0), 
    #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.text=element_text(size=10), 
    #legend.position = 'bottom',
    text = element_text(size=10, family="Times")) -> Plot_tau

Plot_tau

## bacground distribution
ggplot(Tau_background, aes(x=Background)) + 
  geom_histogram(aes(y=..density..), colour="#C0C0C0", fill="#DCDCDC") +
  geom_density(alpha=.1, fill="#000000") +
  geom_segment(data=DF_Tau_tested,
               aes(y=0, yend=6, x=Tau, xend=Tau),
               linetype="dashed", 
               color = "#8A2BE2") +
  scale_y_continuous(expand = c(0,0), limits = c(0,8)) + 
  scale_x_continuous(expand = c(0,0)) + 
  geom_text_repel(data=DF_Tau_tested,
                  aes(x=Tau, y=6, label=Name),
                  direction    = "y",
                  xlim = c(NA, 0.95),
                  ylim = c(6.3, NA),
                  vjust = 1, 
                  segment.size = 0.3,
                  max.iter = 1e4, 
                  box.padding = 0.1,
                  max.time = 1,
                  max.overlaps = Inf,
                  force_pull = 1,
                  size=2) + 
  theme_pubclean() +
  ylab('Density') + 
  xlab('Z-score (Density)') +
  labs(subtitle='Tau null distribution for n=14') + 
  theme(#strip.text.y = element_text(size = 5, angle = 0), 
    #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.text=element_text(size=10), 
    #legend.position = 'bottom',
    text = element_text(size=10, family="Times")) -> Plot_backpround

Plot_backpround


DF_Tau_tested[DF_Tau_tested$Pval > 0.05,]

Fig_9 <- (Plot_GAN / Plot_GRN / Plot_eGRN / Plot_CEN + plot_layout(nrow = 4, ncol = 1)) | 
         (Plot_b  / Plot_tau /  Plot_backpround + plot_layout(nrow = 3, ncol = 1, heights = c(0.3, 1, 1))) 
Fig_9

pdf("Plots/Fig_S9.pdf", width=8, height=9)
print(Fig_9)
dev.off()

##################################################


