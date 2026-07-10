library(tm)
library(SnowballC)
library(wordcloud)
library(RColorBrewer)
library(tidyverse)
library(ggrepel)
library(ggpubr)
library(data.table)
library(GeneOverlap)
library(viridis)
library(fgsea)
library(reshape2)
library(scales)
library(topGO)
library(GOSemSim)
library(enrichplot)
library(Rgraphviz)
library(ComplexHeatmap)
library(AnnotationHub)
library(GOSemSim)
library(enrichplot)
library(ggdark)
library(ggVennDiagram)

###################################
#####        Functions        #####
###################################

ReplaceName <- function(ids){
  
  for (i in 1:nrow(TFdic)){
    ids <- gsub(TFdic$TF.v4[i], TFdic$TF.Name[i], ids)
  }
  return(ids)
}

ReplaceNamePWY <- function(ids){
  
  for (i in 1:nrow(CornCYC)){
    w <- paste0('\\<', CornCYC$Pathway.id[i], '\\>')
    ids <- gsub(w, CornCYC$Pathway.name[i], ids)
    #ids <- gsub("_", " ", ids)
  }
  return(ids)
}

list_to_DF <- function(list){
  # This function get a list of GOs and GeneIds to produce a table
  Net <- as_tibble(as.data.frame(matrix(0, nrow = 0, ncol = 2)))
  colnames(Net) <- c("GO.ID", "GeneID")
  
  GOs <- names(list)
  
  for (n in GOs){
    tem <- tibble(GO.ID=n, GeneID=list[[n]])
    Net <- rbind(Net, tem)
  }
  Net <- unique(Net)
  return(Net)
}

GetGO <- function(degs, mutant){
  
  # Use a list of DEGs and the name of the mutant (string)
  # to identify GOs enriched. Required to have a background predefined
  #
  # mutant: description of a category of source of the data
  # degs: set of genes to test enrichment
  print(length(degs))
  GeneList <- factor(as.integer(background_IDs %in% degs))
  names(GeneList) <- background_IDs
  GOdata_BP <- new("topGOdata", ontology = "BP", allGenes = GeneList, annot = annFUN.gene2GO, gene2GO = background)
  #GOdata_MF <- new("topGOdata", ontology = "MF", allGenes = GeneList, annot = annFUN.gene2GO, gene2GO = background)
  #GOdata_CC <- new("topGOdata", ontology = "CC", allGenes = GeneList, annot = annFUN.gene2GO, gene2GO = background)
  
  #### Define test ####
  test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test", nodeSize = 10)
  
  ### test enrichment 
  results_BP <- getSigGroups(GOdata_BP, test.stat)
  
  ### save pdf Graph
  #namepdf=paste("GOs_Plots/GO.BP_",mutant, "", sep = "")
  #printGraph(GOdata_BP, results_BP, firstSigNodes=20,  fn.prefix = namepdf, useInfo = "def", pdfSW = TRUE) #
  
  ######## Get Significant GOs ########  
  Res_DF_BP <- as_tibble(GenTable(GOdata_BP, classic = results_BP, topNodes = 500, orderBy='Fis')) # save as dataframe
  Res_DF_BP["Mutant"] <- mutant # add Mutant column name
  Res_DF_BP$classic <- as.numeric(Res_DF_BP$classic)
  
  ##### get all GOs and their genes from the topGO result #####
  gs <- genesInTerm(GOdata_BP) # list genes by GO
  # 
  ANOTATION = lapply(gs,function(x) x[x %in% degs]) ## Get only my Differential expressed genes
  
  ### Get only the GO's located in the result of topGO in Res_DF
  DF_GO_Genes <- ANOTATION[Res_DF_BP$GO.ID] # list
  
  ## Transform it to a data frame.
  DF_GO_Genes = list_to_DF(DF_GO_Genes)
  DF_GO_Genes <- left_join(DF_GO_Genes, Res_DF_BP[,c(1,2,6)], by='GO.ID')   # left join to add GO info
  DF_GO_Genes["Mutant"] <- mutant # add Mutant column name
  
  return(list(GOs_DF=Res_DF_BP, Genes=DF_GO_Genes))# return list of GOs-Stats and GeneID-GOs
}

toSpace <- content_transformer(function (x , pattern ) gsub(pattern, " ", x))

Get_WC <- function(GO_list) {
  
  # Read the text vecto: if file to reads
  # text <- readLines(GO_list)
  
  ## Load the data as a corpus
  # df <- data.frame(doc_id = 'GOs', text = BP_PDIs$Term, stringsAsFactors = FALSE)
  # docs <- Corpus(DataframeSource(df))
  
  docs <- Corpus(VectorSource(GO_list)) # vector soruce
  
  
  ## text transformation
  docs <- tm_map(docs, toSpace, "/")
  docs <- tm_map(docs, toSpace, "@")
  docs <- tm_map(docs, toSpace, "\\|")
  docs <- tm_map(docs, toSpace, "\\.")
  docs <- tm_map(docs, toSpace, "process")
  docs <- tm_map(docs, toSpace, "response")
  
  ## Cleaning the text
  # Convert the text to lower case
  docs <- tm_map(docs, content_transformer(tolower))
  # Remove numbers
  docs <- tm_map(docs, removeNumbers)
  # Remove english common stopwords
  docs <- tm_map(docs, removeWords, stopwords("english"))
  
  # Remove your own stop word
  # specify your stopwords as a character vector
  #docs <- tm_map(docs, removeWords, c("blabla1", "blabla2")) 
  # Remove punctuations
  docs <- tm_map(docs, removePunctuation)
  # Eliminate extra white spaces
  docs <- tm_map(docs, stripWhitespace)
  #inspect(docs)
  
  # Text stemming
  # docs <- tm_map(docs, stemDocument)
  
  ## Build a term-document matrix
  dtm <- TermDocumentMatrix(docs)
  m <- as.matrix(dtm)
  v <- sort(rowSums(m),decreasing=TRUE)
  d <- data.frame(word = names(v),freq=v)
  d <- subset(d, freq >=3)
  
  # plot
  set.seed(1234)
  wordcloud(words = d$word, freq = d$freq, min.freq = 3,
            max.words=400, random.order=FALSE, rot.per=0.35, 
            colors=brewer.pal(10, "Paired"), fixed.asp=T) 
}

ReadGOs <- function(list_DEGs) {
  # Uses a list of files to read
  names <- sapply(strsplit(list_DEGs, split='/', fixed=TRUE), `[`, 2) # add methods label
  names <- gsub('GOs.', "", names)
  names <- gsub('.txt', "", names)
  
  
  go_df <- lapply(list_DEGs, fread)
  names(go_df) <- names
  
  go_df <- as_tibble(rbindlist(go_df, idcol = T))
  
  #
  go_df$classic[is.na(go_df$classic)] <- 0
  go_df <- go_df[go_df$classic <= 0.05, ]
  
  return(go_df)
  
}

Get_SS_PDI_CoExp <- function(tfslist){
  # make df to save output
  GOSIM_DF <- tibble(Source=tfslist, PDI_CoExp_SS=0)
  
  c= 1
  for (tf in tfslist){
    # PDI
    go1 <- subset(BP_PDIs, Mutant == tf)$GO.ID
    # CoExp
    go2 <- subset(BP_CoExp, Mutant == tf)$GO.ID
    
    ss <- mgoSim(go1, go2, semData=zmGO_BP, measure="Wang", combine="BMA")
    
    GOSIM_DF$PDI_CoExp_SS[GOSIM_DF$Source==tf]  <- ss
    
    print(paste0(" .. Done ", c, " TFs .."))
    c = c+1
    
  }
  return(GOSIM_DF)
  
}

Get_SS_QTL_PDI <- function(tfslist){
  # make df to save output
  GOSIM_DF <- tibble(Source=tfslist, QTL_PDI_SS=0)
  
  c= 1
  for (tf in tfslist){
    
    # QTL
    go1 <- subset(BP_teQTL, Mutant == tf)$GO.ID
    # CoExp
    go2 <- subset(BP_PDIs, Mutant == tf)$GO.ID
    
    ss <- mgoSim(go1, go2, semData=zmGO_BP, measure="Wang", combine="BMA")
    
    GOSIM_DF$QTL_PDI_SS[GOSIM_DF$Source==tf]  <- ss
    
    print(paste0(" .. Done ", c, " TFs .."))
    c = c+1
    
  }
  return(GOSIM_DF)
  
}

Get_SS_QTL_CoExp <- function(tfslist){
  # make df to save output
  GOSIM_DF <- tibble(Source=tfslist, QTL_CoExp_SS=0)
  
  c= 1
  for (tf in tfslist){
    
    # QTL
    go1 <- subset(BP_teQTL, Mutant == tf)$GO.ID
    # CoExp
    go2 <- subset(BP_CoExp, Mutant == tf)$GO.ID
    
    ss <- mgoSim(go1, go2, semData=zmGO_BP, measure="Wang", combine="BMA")
    
    GOSIM_DF$QTL_CoExp_SS[GOSIM_DF$Source==tf]  <- ss
    
    print(paste0(" .. Done ", c, " TFs .."))
    c = c+1
    
  }
  return(GOSIM_DF)
  
}


Ztest <- function(value, randomlist){
  zv <- (value - mean(randomlist))/sd(randomlist)
  Pval <- pnorm(zv, 0,1, lower.tail = F)         ## zv is larged thant Random values? 
  Zresults <- c(round(zv,3), Pval, value)
  return(Zresults)
}


Get_SS_PDI_CoExp_Matrix <- function(tfslist){
  # make df to save output
  GOSIM_DF <- list()
  
  c= 1
  for (tf in tfslist){
    # PDI
    #c("Zm00001d006236")
    
    go1 <- subset(BP_PDIs, Mutant == tf)$GO.ID
    # CoExp
    go2 <- subset(BP_CoExp, Mutant == tf)$GO.ID
    
    ss <- mgoSim(go1, go2, semData=zmGO_BP, measure="Wang", combine=NULL)
    
    GOSIM_DF[[tf]] <- ss
    
    print(paste0(" .. Done ", c, " TFs .."))
    c = c+1
    
  }
  return(GOSIM_DF)
  
}

vennfunc <- function(list_x){
  venn <- Venn(list_x)
  data <- process_data(venn)
  colorGroups <- c(CEN = 'goldenrod1', GRN='steelblue1', GAN='darkorchid1')
  colfunc <- colorRampPalette(colorGroups)
  col <- colfunc(3)
  
  colorGroups <- c(CEN="gray100",GRN="gray99", GAN="gray98")
  colfunc2 <- colorRampPalette(colorGroups)
  col2 <- colfunc2(7)
  
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
    theme(plot.margin = unit(c(0.5,1,1,0.1), "cm")) +
    xlim(-150,1000)
}

##############################################################################
##################         Read data input           #########################
##############################################################################

#### top 45
Top45 <- as_tibble(read.table("Data/Annotations/Top45.txt", h=F, stringsAsFactors = F))

# TF names
TFdic <- as_tibble(read.table("Data/Annotations/TF_Id_Name.txt", h=T, stringsAsFactors = F))

# Phenolic related genes
PheGenes <- as_tibble(read.table("Data/Annotations/LinaPheGenes2020.txt", h=T, sep = "\t", quote="", stringsAsFactors = F))

# TF and CoReg
TF_CoR <- as_tibble(read.table("Data/Annotations/TF_CoR_Mazie.txt", h=T, stringsAsFactors = F)) 

## Y1H network
Y1H <- as_tibble(read.table("Data/Annotations/Y1H.data.txt", h=T, stringsAsFactors = F))[,4:6]

#table(subset(PDI.class, TF %in% Y1H$TF.v4)$source)

## syntenic genes 
Syntenic <- as_tibble(read.table("Data/Annotations/Zm.v4.synteny.genes.txt", h=T, stringsAsFactors = F, sep = '\t'))
Syntenic <- subset(Syntenic, gene_synteny=='syntenic')$gene_id

# GOs term annotations
background <- readMappings("synteny.ID_TopGO_V4_GRAMER.txt")
background_IDs <- as.character(unique(names(background)))
background_list <- unique(as.character(unlist(background)))

# PDI
PDI <- as_tibble(unique(fread("../Fig_PDI/PDI_NetworkFinal.10_11_2021.txt")[,c(2,3)]))
colnames(PDI)[1] <- "Source"
# teQTL
teQTL <- as_tibble(unique(fread("../Fig_transeQTL/teQTL_NetworkFinal.10_11_2021.txt")))
colnames(teQTL)[1] <- "Source"

# CoExp
CoExp <- as_tibble(unique(fread("../Fig_Coexpression/CoExp_NetworkFinal.10_11_2021.txt")))
colnames(CoExp)[2] <- "Source"

# PDI class 
PDI.class <- as_tibble(fread("Data/PDI_data/scATAC.Z.Full.Net.Dis2TSS.txt"))
PDI.class <- PDI.class[,c(1)]
PDI.class[, "Method"] <- sapply(strsplit(PDI.class$TFsample, split='.', fixed=TRUE), `[`, 1) # add methods label
PDI.class <- unique(PDI)
table(PDI.class$Method)
PDI.class$TFsample <- gsub("ChIP.", "", PDI.class$TFsample)
PDI.class$TFsample <- gsub("newDAP.", "", PDI.class$TFsample)
PDI.class$TFsample <- gsub("_Met", "", PDI.class$TFsample)
PDI.class$TFsample <- gsub("_deM", "", PDI.class$TFsample)
PDI.class$TFsample <- gsub(".tasse", "", PDI.class$TFsample)
PDI.class$TFsample <- gsub("DAP.", "", PDI.class$TFsample)
PDI.class$TFsample <- gsub(".ear", "", PDI.class$TFsample)
PDI.class$TFsample <- gsub("pZm", "Zm", PDI.class$TFsample)
PDI.class$TF <- gsub("pGRMZM", "GRMZM", PDI.class$TF)
PDI.class$TF <- gsub("Zm00001d020430l", "Zm00001d020430", PDI.class$TF)

PDI.class$TF
colnames(PDI.class)[1] <- "TF"
PDI.class <- unique(PDI.class)

################################################################
########   Count common TFs and number of associations  ######## 
################################################################

## Count common TFs and number of associations
#
PDI_counts <- as_tibble(as.data.frame(table(PDI$Source), stringsAsFactors = F))
colnames(PDI_counts) <- c("Source", "Targets")
#
CoExp_counts <- as_tibble(as.data.frame(table(unique(CoExp[,2:3])$Source), stringsAsFactors = F))
colnames(CoExp_counts) <- c("Source", "Targets")
#
teQTL_counts <- as_tibble(as.data.frame(table(teQTL$Source), stringsAsFactors = F))
colnames(teQTL_counts) <- c("Source", "Targets")

TFs_common <- tibble(Source=unique(c(PDI_counts$Source, CoExp_counts$Source, teQTL_counts$Source)))
#
TFs_common <- left_join(TFs_common, PDI_counts, by='Source')
TFs_common <- left_join(TFs_common, CoExp_counts, by='Source')
TFs_common <- left_join(TFs_common, teQTL_counts, by='Source')
colnames(TFs_common) <- c("Source", "PDI", "CoExp", "teQTL")

TFs_common[is.na(TFs_common)] <- 0

mask <- ((TFs_common$PDI>50)*1 + (TFs_common$CoExp>50)*1 ) >=2

TFs_common_20 <- as.data.frame(TFs_common[mask,])
row.names(TFs_common_20) <- TFs_common_20$Source
TFs_common_20 <- TFs_common_20[,-c(1)]

row.names(TFs_common_20) <- ReplaceName(row.names(TFs_common_20))

Hm_TFs_common <- Heatmap(t(as.matrix(log2(TFs_common_20+1))), name="log2 Targets",
                       # column_km = 3,
                       # column_names_rot = 90,
                       # row_names_rot = 45,
                       cluster_rows = TRUE, 
                       cluster_columns = FALSE,
                       show_column_dend = TRUE, show_row_dend = FALSE, 
                       clustering_method_columns = "ward.D2",
                       col=viridis(100, direction = 1, option = "A"),
                       column_names_gp = gpar(fontsize = 5),
                       row_names_gp = gpar(fontsize = 10),
                       show_heatmap_legend = T,
                       # heatmap_ = unit(10),
                       heatmap_height = unit(7, 'cm'),
                       heatmap_width  = unit(20, 'cm'),
                       heatmap_legend_param = list(title_gp = gpar(fontsize = 10), 
                                                   labels_gp = gpar(fontsize = 10), 
                                                   direction = "horizontal"))
Hm_TFs_common
# size=3x10
draw(Hm_TFs_common, heatmap_legend_side = "bottom")


################################################################


#################################################################
#############           Read GOs terms results      #############
#################################################################

Common_TFs_to_test <- TFs_common[mask,]$Source # based on PDI and PDI data

TFs_112 <- subset(TFs_common, PDI>0 & CoExp> 0 & teQTL >0)$Source

length(Common_TFs_to_test[(Common_TFs_to_test %in% TFs_112)])

#################################################################
##  Make file names of Sig. GOs to reads
#################################################################

# GO BP PDI
GO_bp_pdi <- list.files(pattern = "BP_PDI.GOs.*", path = 'BP_results/')
GO_bp_pdi <- paste0('BP_results/', GO_bp_pdi)

# GO BP CoExp
GO_bp_coexp <- list.files(pattern = "BP_CoExp.GOs.*", path = 'BP_results/')
GO_bp_coexp <- paste0('BP_results/', GO_bp_coexp)

# GO BP teQTL: 
GO_bp_teQTL <- list.files(pattern = "BP_teQTL.GOs.*", path = 'BP_results/')
GO_bp_teQTL <- paste0('BP_results/', GO_bp_teQTL)


# GO MF PDI
GO_mf_pdi <- list.files(pattern = "MF_PDI.GOs.*", path = 'MF_results/')
GO_mf_pdi <- paste0('MF_results/', GO_mf_pdi)

# GO MF CoExp
GO_mf_coexp <- list.files(pattern = "MF_CoExp.GOs.*", path = 'MF_results/')
GO_mf_coexp <- paste0('MF_results/', GO_mf_coexp)

#teQTL_counts[teQTL_counts$Source=="Zm00001d033859", ]
#teQTL[teQTL$Source=="Zm00001d033859", ]
#
#################################################################

####################################################
########            Example: KN1 and O2
####################################################
kn1_bp_pdi <- GO_bp_pdi[grepl("Zm00001d033859", GO_bp_pdi)]
kn1_bp_coexp <- GO_bp_coexp[grepl("Zm00001d033859", GO_bp_coexp)]
#
O2_bp_pdi <- GO_bp_pdi[grepl("Zm00001d018971", GO_bp_pdi)]
O2_bp_coexp <- GO_bp_coexp[grepl("Zm00001d018971", GO_bp_coexp)]

#kn1_bp_teQTL <- GO_bp_teQTL[grepl("Zm00001d033859", GO_bp_teQTL)]

#kn1_mf_pdi <- GO_mf_pdi[grepl("Zm00001d033859", GO_mf_pdi)]
#kn1_mf_coexp <- GO_mf_coexp[grepl("Zm00001d033859", GO_mf_coexp)]

#O2_mf_pdi <- GO_mf_pdi[grepl("Zm00001d018971", GO_mf_pdi)]
#O2_mf_coexp <- GO_mf_coexp[grepl("Zm00001d018971", GO_mf_coexp)]

ReadGOs(kn1_bp_pdi)     # kn1_bp_pdi
ReadGOs(kn1_bp_coexp)   # kn1_bp_coex
#
ReadGOs(O2_bp_pdi)      # kn1_bp_pdi
ReadGOs(O2_bp_coexp)    # kn1_bp_coexp

# MF
#MF_PDIs <- ReadGOs(O2_mf_pdi)
#MF_CoExp <- ReadGOs(O2_mf_coexp)

# size 7x6
Get_WC(BP_PDIs$Term)
Get_WC(BP_CoExp$Term)
####################################################

################################################################
## Read sig. GOs for corresponding TF 
################################################################

# BP
BP_PDIs <- ReadGOs(GO_bp_pdi)      # Read all GOs files for PDI data
BP_CoExp <- ReadGOs(GO_bp_coexp)   # Read all GOs files for CoExp data
BP_teQTL <- ReadGOs(GO_bp_teQTL)   # Read all GOs files for CoExp data

# Filter to TFs
BP_PDIs   <- subset(BP_PDIs, Mutant %in% TF_CoR$GeneID)
BP_CoExp  <- subset(BP_CoExp, Mutant %in% TF_CoR$GeneID)
BP_teQTL  <- subset(BP_teQTL, Mutant %in% TF_CoR$GeneID)


BP_TFs <- rbind(tibble(TFs=unique(BP_PDIs$Mutant), Net="GRN"),
                tibble(TFs=unique(BP_CoExp$Mutant), Net="CEN"),
                tibble(TFs=unique(BP_teQTL$Mutant), Net="GAN"))

BP_TFs <- subset(BP_TFs, TFs %in% TFs_112)
BP_TFs <- split(BP_TFs$TFs, BP_TFs$Net)

vennfunc <- function(list_x){
  venn <- Venn(list_x)
  data <- process_data(venn)
  colorGroups <- c(CEN = 'goldenrod1', GRN='steelblue1', GAN='darkorchid1')
  colfunc <- colorRampPalette(colorGroups)
  col <- colfunc(3)
  
  colorGroups <- c(CEN="gray100",GRN="gray99", GAN="gray98")
  colfunc2 <- colorRampPalette(colorGroups)
  col2 <- colfunc2(7)
  
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
    theme(plot.margin = unit(c(0.5,1,1,0.1), "cm")) +
    xlim(-150,1000)
}

vennfunc(BP_TFs)

length(Common_TFs_to_test)

write.table(BP_PDIs, "Results_BP_PDIs.txt", row.names = F, sep = "\t", quote = F)
write.table(BP_CoExp, "Results_BP_CoExp.txt", row.names = F, sep = "\t", quote = F)
write.table(BP_teQTL, "Results_BP_teQTL.txt", row.names = F, sep = "\t", quote = F)

length(table(BP_CoExp$Mutant))

# Reduce data sets to TFs with PDI and CoExp data
BP_PDIs <- subset(BP_PDIs, Mutant %in% Common_TFs_to_test)
BP_CoExp <- subset(BP_CoExp, Mutant %in% Common_TFs_to_test)
BP_teQTL <- subset(BP_teQTL, Mutant %in% Common_TFs_to_test)

# Count GOs by TF
BP_PDIs_freq <- as_tibble(as.data.frame(table(unique(BP_PDIs[,c(2,8)])$Mutant)))
BP_CoExp_freq <- as_tibble(as.data.frame(table(unique(BP_CoExp[,c(2,8)])$Mutant)))
BP_teQTL_freq <- as_tibble(as.data.frame(table(unique(BP_teQTL[,c(2,8)])$Mutant)))

BP_PDIs_freq["Net"] <- "GRN"
BP_CoExp_freq["Net"] <- "CEN"
BP_teQTL_freq["Net"] <- "GAN"

BP_freq <- rbind(BP_PDIs_freq, BP_CoExp_freq, BP_teQTL_freq)
table(BP_freq$Net)

length(unique(teQTL[teQTL$Source %in% TF_CoR$GeneID,]$Source))
length(unique(PDI[PDI$Source %in% TF_CoR$GeneID,]$Source))
length(unique(CoExp[CoExp$Source %in% TF_CoR$GeneID,]$Source))

################################################################
## Figure 
################################################################

Plot_GO_freq <- ggplot(BP_freq, aes(x=Net, y=Freq))+
  #geom_point() +
  geom_jitter(size=2, alpha=0.5, width = 0.2, color="black") + 
  theme_pubclean()

Plot_GO_freq

################################################################


########################################################
#####    TFs with at least a GO term enriched.  ########
########################################################

TFs_with_GOs_list = list(CEN=unique(BP_CoExp$Mutant), 
                         GRN=unique(BP_PDIs$Mutant),
                         GAN=unique(BP_teQTL$Mutant))

vennfunc(TFs_with_GOs_list) + labs(title="TFs with GOs Enriched")
########################################################


########################################################
#####    Calculate semantic similarities        ########
########################################################

####################################################
## 1. Calculate semantic similarity background

ah = AnnotationHub()
zm <- query(ah, c("Zea mays"))
zm <- zm["AH85440"]

zm_hub <- ah[[zm$ah_id]] # Extract info for specific maize id

# Create Semantic similarity calculations: in this cases based on BP
zmGO_BP <- godata(zm_hub, ont="BP")

## 2. Calculate semantic similarity measurement for pair of PDI and CoExp GOs
DF_SS_PDI_CoExp <- Get_SS_PDI_CoExp(Common_TFs_to_test)

## 3. Calculate semantic similarity measurement for pair of teQTL and PDI GOs
Common_TFs_to_test_qtl <- c("Zm00001d001945", "Zm00001d009599", "Zm00001d010634", 
                            "Zm00001d011953", "Zm00001d012916", "Zm00001d014995", 
                            "Zm00001d015407", "Zm00001d015421", "Zm00001d016838", 
                            "Zm00001d017726", "Zm00001d020267", "Zm00001d024324",
                            "Zm00001d031182", "Zm00001d033719", "Zm00001d042267",
                            "Zm00001d050195", "Zm00001d050781")


DF_SS_QTL_PDI <- Get_SS_QTL_PDI(Common_TFs_to_test_qtl)

## 4. Calculate semantic similarity measurement for pair of teQTL and CoExp GOs
DF_SS_QTL_CoExp <- Get_SS_QTL_CoExp(Common_TFs_to_test_qtl)

####################################################
###   Read GO_ss PDI vs random CoExp files       ###
####################################################

rGOss_bp_pdi <- list.files(pattern = "GO_SS_random_PDI_CoExp.*", path = 'GO_SS_RandomBackground/')

# TFs list in random data
TFs_rGOss_bp_pdi <- gsub(".txt", "", gsub("GO_SS_random_PDI_CoExp.", "", rGOss_bp_pdi))
# set files to read
rGOss_bp_pdi <- paste0('GO_SS_RandomBackground/', rGOss_bp_pdi)

# Read random GO SS files
Random_pdi_coExp_ss <- lapply(rGOss_bp_pdi, fread)
Random_pdi_coExp_ss <- as_tibble(rbindlist(Random_pdi_coExp_ss, idcol = F))
Random_pdi_coExp_ss <- split(Random_pdi_coExp_ss$r.ss, Random_pdi_coExp_ss$ID)


##############################################################
###   Read GO_ss teQTL vs random PDI and CoExp files       ###
##############################################################

rGOss_bp_tqtls <- list.files(pattern = "GO_SS_random_teQTL_PDI_CoExp.*", path = 'GO_SS_RandomBackground/')

# TFs list in random data
TFs_rGOss_bp_tqtls <- gsub(".txt", "", gsub("GO_SS_random_teQTL_PDI_CoExp.", "", rGOss_bp_tqtls))
# set files to read
rGOss_bp_tqtls <- paste0('GO_SS_RandomBackground/', rGOss_bp_tqtls)

# Read random GO SS files
Random_teqtl_pdi_coExp_ss <- lapply(rGOss_bp_tqtls, fread)
Random_teqtl_pdi_coExp_ss <- as_tibble(rbindlist(Random_teqtl_pdi_coExp_ss, idcol = F))

Random_teqtl_pdi_ss  <- split(Random_teqtl_pdi_coExp_ss$r.ss.pdi, Random_teqtl_pdi_coExp_ss$ID)
Random_teqtl_coExp_ss<- split(Random_teqtl_pdi_coExp_ss$r.ss.coe, Random_teqtl_pdi_coExp_ss$ID)


####################################################
######           Z test for all TFs           ######
####################################################

# Make observed GSS values a list
# PDI and CoEXP
DF_SS_PDI_CoExp_list <- split(DF_SS_PDI_CoExp$PDI_CoExp_SS, DF_SS_PDI_CoExp$Source)
DF_SS_PDI_CoExp_list
# etQTL and PDI
DF_SS_QTL_PDI_list <- split(DF_SS_QTL_PDI$QTL_PDI_SS, DF_SS_QTL_PDI$Source)
# etQTL and CoExp
DF_SS_QTL_CoExp_list <- split(DF_SS_QTL_CoExp$QTL_CoExp_SS, DF_SS_QTL_CoExp$Source)


# make Ztest using random SS: PDI vs CoExp
Table_rGss_Reasult  <- as.data.frame(t(mapply(Ztest, DF_SS_PDI_CoExp_list[TFs_rGOss_bp_pdi], 
                                              Random_pdi_coExp_ss[TFs_rGOss_bp_pdi])),
                                    stringAsFactor=F)

Table_rGss_Reasult[,"TFs"] <- TFs_rGOss_bp_pdi
Table_rGss_Reasult <- as_tibble(Table_rGss_Reasult)
colnames(Table_rGss_Reasult) <- c("Z", "Pval", "ss","TF")

# make Ztest using random SS: teQTL vs PDI
Table_rGss_Result_teqtl_pdi  <- as.data.frame(t(mapply(Ztest, DF_SS_QTL_PDI_list[TFs_rGOss_bp_tqtls], 
                                                       Random_teqtl_pdi_ss[TFs_rGOss_bp_tqtls])),
                                     stringAsFactor=F)

Table_rGss_Result_teqtl_pdi[,"TFs"] <- TFs_rGOss_bp_tqtls
Table_rGss_Result_teqtl_pdi <- as_tibble(Table_rGss_Result_teqtl_pdi)
colnames(Table_rGss_Result_teqtl_pdi) <- c("Z", "Pval", "ss","TF")

# make Ztest using random SS: teQTL vs CoExp
Table_rGss_Result_teqtl_coe  <- as.data.frame(t(mapply(Ztest, 
                                                       DF_SS_QTL_CoExp_list[TFs_rGOss_bp_tqtls], 
                                                       Random_teqtl_coExp_ss[TFs_rGOss_bp_tqtls])),
                                              stringAsFactor=F)

Table_rGss_Result_teqtl_coe[,"TFs"] <- TFs_rGOss_bp_tqtls
Table_rGss_Result_teqtl_coe <- as_tibble(Table_rGss_Result_teqtl_coe)
colnames(Table_rGss_Result_teqtl_coe) <- c("Z", "Pval", "ss","TF")


# Add FDr value and TF name
Table_rGss_Reasult[,"FDR"] <- p.adjust(Table_rGss_Reasult$Pval, method = 'fdr')
Table_rGss_Result_teqtl_pdi[,"FDR"] <- p.adjust(Table_rGss_Result_teqtl_pdi$Pval, method = 'fdr')
Table_rGss_Result_teqtl_coe[,"FDR"] <- p.adjust(Table_rGss_Result_teqtl_coe$Pval, method = 'fdr')

Table_rGss_Reasult["Net"] <- "CEN-GRN"
Table_rGss_Result_teqtl_pdi["Net"] <- "GAN-GRN"
Table_rGss_Result_teqtl_coe["Net"] <- "GAN-CEN"

Table_rGss_Reasult["Rank"] <- rank(-Table_rGss_Reasult$ss)
Table_rGss_Result_teqtl_pdi["Rank"] <- rank(-Table_rGss_Result_teqtl_pdi$ss)
Table_rGss_Result_teqtl_coe["Rank"] <- rank(-Table_rGss_Result_teqtl_coe$ss)

# Combine all GSS results
All_rGss_Results <- rbind(Table_rGss_Reasult, 
                          Table_rGss_Result_teqtl_pdi,
                          Table_rGss_Result_teqtl_coe)

All_rGss_Results["TFname"] <- ReplaceName(All_rGss_Results$TF)

# Add PDI method source
#Table_rGss_Reasult <- left_join(Table_rGss_Reasult, PDI.class, by="TF") 
#Table_rGss_Reasult$Method <- factor(Table_rGss_Reasult$Method, levels = c("ChIP", "DAP", "pChIP"))

# add Plot Figure 3a
subset(All_rGss_Results, ss>=0.85)


table(subset(All_rGss_Results, TF %in% TFs_112 & FDR <=0.1)$Net)


Plot_rGss_Results <- ggplot(All_rGss_Results, aes(x=Net, y=ss, color=(FDR<=0.1), label=TFname))+
  geom_jitter(size=2, alpha=0.4, width = 0.02) +
  geom_text_repel(data = subset(All_rGss_Results, TFname %in% c("bHLH145", "COL18") & Net=="CEN-GRN"), show.legend = F)+
  scale_color_manual(values=c("TRUE" = "#F8766D", "FALSE" = "#00BFC4"), labels=c(expression("FDR" <= 0.1), "FDR > 0.1")) +
  theme_pubclean() +
  ylab("GO terms Semantic\nsimilarity (GSS)") + #ylab(expression(-Log[10]~(FDR))) + #xlab("GO terms Semantic similarity (GSS)") +
  xlab("Networks")  + 
  coord_flip()


Plot_rGss_Results <- ggpar(Plot_rGss_Results, font.tickslab = 14, font.x = 14, font.y = 14, 
                           legend = "bottom", legend.title = "")
# Size: 5x4
Plot_rGss_Results

All_rGss_Results <- All_rGss_Results[order(All_rGss_Results$Rank),]

####################################################

####################################################
## Common GOs that pass the GSS filter: CEN-GRN   ##           
####################################################

GetGOGroups <- function(tf){
  
  GO_matrix <- GSS_PASS[[tf]]
  # cor matrix
  test <- melt(GO_matrix)
  colnames(test) <- c("GO_GRN", "GO_CEN", "GSS")
  
  BP_PDIs_test <- subset(BP_PDIs, Mutant==tf & GO.ID %in% test$GO_GRN)
  BP_CoExp_test <- subset(BP_CoExp, Mutant==tf & GO.ID %in% test$GO_CEN)
  
  test <- left_join(test, BP_PDIs_test[,c(2,3,5,7)], by=c("GO_GRN"="GO.ID"))
  test <- left_join(test, BP_CoExp_test[,c(2,3,5,7)], by=c("GO_CEN"="GO.ID"))
  colnames(test) <- c("GO.GRN", "GO.CEN", "GSS", 
                      "Term.GRN",  "Sig.Terms.GRN", "Pval.GRN",
                      "Term.CEN", "Sig.Terms.CEN", "Pval.CEN")
  # GSS filter based on overall CEN and GRN Gss
  val <- subset(Table_rGss_Reasult, TF==tf)$ss
  test <- subset(test, GSS >= val)
  
  
  return(as_tibble(test))
  
}

GetGOGroupsV2 <- function(tf){
  
  GO_matrix <- GSS_PASS[[tf]]
  # cor matrix
  test <- melt(GO_matrix)
  colnames(test) <- c("GO_GRN", "GO_CEN", "GSS")
  
  BP_PDIs_test <- subset(BP_PDIs, Mutant==tf & GO.ID %in% test$GO_GRN)
  BP_CoExp_test <- subset(BP_CoExp, Mutant==tf & GO.ID %in% test$GO_CEN)
  
  test <- left_join(test, BP_PDIs_test[,c(2,3,5,7)], by=c("GO_GRN"="GO.ID"))
  test <- left_join(test, BP_CoExp_test[,c(2,3,5,7)], by=c("GO_CEN"="GO.ID"))
  colnames(test) <- c("GO.GRN", "GO.CEN", "GSS", 
                      "Term.GRN",  "Sig.Terms.GRN", "Pval.GRN",
                      "Term.CEN", "Sig.Terms.CEN", "Pval.CEN")
  
  # GSS filter based on overall CEN and GRN Gss
  TablePvals  <- as.data.frame(t(sapply(test$GSS, Ztest, Random_pdi_coExp_ss[[tf]])))
  TablePvals <- as_tibble(TablePvals)[,1:2]
  colnames(TablePvals) <- c("Z", "Pval")
  
  test <- cbind(test, TablePvals)
  test["FDR"] <- p.adjust(test$Pval, method = "fdr")
  
  test <- subset(test, FDR <= 0.1)
  
  
  return(as_tibble(test))
  
}


# TFs with Significant GSS: Coexp and PDI
SigTFs_GSS <- All_rGss_Results[All_rGss_Results$FDR <=0.1,]$TF
ReplaceName(SigTFs_GSS)

# list of GSS score among all PDI and Coexp GO terms
GSS_PASS <- Get_SS_PDI_CoExp_Matrix(SigTFs_GSS)

# list of individual GSS that pass global score: PDI and Coexp GO terms
GO_table  <- lapply(SigTFs_GSS, GetGOGroups) 
names(GO_table) <- SigTFs_GSS

# list of individual GSS that pass Z test: random PDI and Coexp GO terms
GO_table_v2  <- lapply(SigTFs_GSS, GetGOGroupsV2) 
names(GO_table_v2) <- SigTFs_GSS


# Save
ReplaceName('Zm00001d033673')

write.table(as_tibble(rbindlist(GO_table, idcol = T)), "Sig_GOs_CEN_GRN.txt", row.names = F, sep = "\t", quote = F)
write.table(as_tibble(rbindlist(GO_table_v2, idcol = T)), "Sig_GOs_CEN_GRN.v2.txt", row.names = F, sep = "\t", quote = F)
 
Get_WC(GO_table$Zm00001d033673$Term.GRN)
Get_WC(GO_table$Zm00001d033673$Term.CEN)

Sig_GOs_Summary <- function(DFpass) {
  
  g1 <- DFpass[,1]
  colnames(g1) <- "GO"
  g2 <- DFpass[,2]
  colnames(g2) <- "GO"
  
  g <- unique(rbind(g1, g2))
  
  return(length(unique(g1$GO)))
  
}

# Count number of significant GOs by TF combined both list
TableSig_GOs_Summary <- as_tibble(t(as.data.frame(lapply(GO_table_v2, Sig_GOs_Summary))))
TableSig_GOs_Summary["TF"] <- names(GO_table_v2)
colnames(TableSig_GOs_Summary)[1] <- "GOs"

TableSig_GOs_Summary["TFname"] <- ReplaceName(TableSig_GOs_Summary$TF)
TableSig_GOs_Summary

summary(TableSig_GOs_Summary$GOs)

FigureSigCommonGOs <- ggplot(TableSig_GOs_Summary, aes(x="TFs", y=GOs, label=TFname))+
  geom_jitter(position = position_jitter(seed = 1, width = 0.15)) +
  geom_text_repel(position = position_jitter(seed = 1, width = 0.15))+
  theme_pubclean() +
  ylab("Number of GOs significanly similars") +
  xlab("") +  coord_flip()

# Size: 3x10
FigureSigCommonGOs <- ggpar(FigureSigCommonGOs, 
                            font.tickslab = 14, font.x = 14, font.y = 14)


View(TableSig_GOs_Summary)
write.table(TableSig_GOs_Summary[,2], 
            "/maindisk/fabio/Projects/MaizeENCODE/Data_45_net/wPCC_net_only_TFs/TF_list.txt",
            sep="\t", quote = F, row.names = F)

#unique(All_rGss_Results[,4])
#write.table(TableSig_GOs_Summary[,2], "/maindisk/fabio/Projects/MaizeENCODE/Data_45_net/wPCC_net_only_TFs/TF_list.txt", sep="\t", quote = F, row.names = F)

Get_WC_Freq <- function(GO_list) {
  
  # Read the text vecto: if file to reads
  # text <- readLines(GO_list)
  
  ## Load the data as a corpus
  df <- data.frame(doc_id = 1, text = GO_list, stringsAsFactors = FALSE)
  docs <- Corpus(DataframeSource(df))
  
  #docs <- Corpus(VectorSource(GO_list)) # vector soruce
  
  ## text transformation
  # Remove english common stopwords
  docs <- tm_map(docs, removeWords, stopwords("english"))
  # Remove punctuations
  docs <- tm_map(docs, removePunctuation)
  # Eliminate extra white spaces
  docs <- tm_map(docs, stripWhitespace)
  #docs <- tm_map(docs, toSpace, "/")
  #docs <- tm_map(docs, toSpace, "@")
  #docs <- tm_map(docs, toSpace, "\\|")
  #docs <- tm_map(docs, toSpace, "\\.")
  docs <- tm_map(docs, removeWords, c("process", "response", "cell", "metabolic", "cellular", "transport"))
  
  
  ## Cleaning the text
  # Convert the text to lower case
  docs <- tm_map(docs, content_transformer(tolower))
  # Remove numbers
  docs <- tm_map(docs, removeNumbers)
  # Remove english common stopwords
  docs <- tm_map(docs, removeWords, stopwords("english"))
  
  # Remove your own stop word
  # specify your stopwords as a character vector
  #docs <- tm_map(docs, removeWords, c("blabla1", "blabla2")) 
  # Remove punctuations
  docs <- tm_map(docs, removePunctuation)
  # Eliminate extra white spaces
  docs <- tm_map(docs, stripWhitespace)
  #inspect(docs)
  
  # Text stemming
  # docs <- tm_map(docs, stemDocument)
  
  ## Build a term-document matrix
  dtm <- TermDocumentMatrix(docs)
  m <- as.matrix(dtm)
  v <- sort(rowSums(m),decreasing=TRUE)
  d <- data.frame(word = names(v),freq=v)
  d <- subset(d, freq >=3)
  
  # plot
  set.seed(1234)
  wordcloud(words = d$word, freq = d$freq, min.freq = 3, max.words=300, random.order=FALSE, rot.per=0,  
            colors=brewer.pal(10, "Paired"), fixed.asp=T) 
  return(d)
}

WC_freq_COL18 <- Get_WC_Freq(unique(c(GO_table_v2$Zm00001d015468$Term.GRN, GO_table_v2$Zm00001d015468$Term.CEN)))
WC_freq_MYB31 <- Get_WC_Freq(unique(c(GO_table_v2$Zm00001d006236$Term.GRN, GO_table_v2$Zm00001d006236$Term.CEN)))
WC_freq_KN1 <- Get_WC_Freq(unique(c(GO_table_v2$Zm00001d033859$Term.GRN, GO_table_v2$Zm00001d033859$Term.CEN)))



####################################################


####################################################
## Identification of top GOs Terms: Heatmap       ##           
## COL18 and bHLH145 Examples                     ##
####################################################

ReplaceName("Zm00001d015468")

length(unique(GO_table$Zm00001d015468$GO.GRN))
length(unique(GO_table$Zm00001d015468$GO.CEN))

subset(Table_rGss_Reasult, -log10(FDR)>=3)
length(unique(subset(BP_PDIs, Mutant == "Zm00001d015468")$GO.ID))
length(unique(subset(BP_CoExp, Mutant == "Zm00001d015468")$GO.ID))

length(unique(subset(BP_PDIs, Mutant == "Zm00001d031717")$GO.ID))
length(unique(subset(BP_CoExp, Mutant == "Zm00001d031717")$GO.ID))

# Get example
COL18_GO_matrix <- GSS_PASS$Zm00001d015468
bHLH145_GO_matrix <- GSS_PASS$Zm00001d031717


#gap_stat_rows <- clusGap(bHLH91_GO_matrix, FUN = kmeans, nstart = 25, K.max = 40, B = 5)
#gap_stat_row_plot <- fviz_gap_stat(gap_stat_rows)
#km_rows <- kmeans(bHLH91_GO_matrix, 14)

#mask <- colnames(bHLH91_GO_matrix)[colnames(bHLH91_GO_matrix) %in% row.names(bHLH91_GO_matrix)]

Hm_SGO_1 <- Heatmap(COL18_GO_matrix, name="GSS",
                  #column_km = 10, #row_km = 10, # column_names_rot = 90, # row_names_rot = 45,
                  column_title = "GOs CEN",
                  row_title = "GOs GRN",
                  cluster_rows = TRUE, cluster_columns = TRUE,
                  #row_split = km_rows$cluster, #column_split = km_cols$cluster,
                  show_column_dend = F, show_row_dend = F, 
                  show_row_names = F, show_column_names = F,
                  clustering_method_columns = "ward.D2",
                  col=viridis(50, direction = -1, option = "A"),
                  column_names_gp = gpar(fontsize = 2),
                  row_names_gp = gpar(fontsize = 2),
                  show_heatmap_legend = T,
                  # heatmap_ = unit(10),
                  heatmap_height = unit(8, 'cm'),
                  heatmap_width  = unit(10, 'cm'),
                  heatmap_legend_param = list(title_gp = gpar(fontsize = 10), 
                                              labels_gp = gpar(fontsize = 10),
                                              direction = "horizontal")) #

Hm_SGO_2 <- Heatmap(bHLH145_GO_matrix, name="GSS",
                    #column_km = 10, #row_km = 10, # column_names_rot = 90, # row_names_rot = 45,
                    column_title = "GOs CEN",
                    row_title = "GOs GRN",
                    cluster_rows = TRUE, cluster_columns = TRUE,
                    #row_split = km_rows$cluster, #column_split = km_cols$cluster,
                    show_column_dend = F, show_row_dend = F, 
                    show_row_names = F, show_column_names = F,
                    clustering_method_columns = "ward.D2",
                    col=viridis(50, direction = -1, option = "A"),
                    column_names_gp = gpar(fontsize = 2),
                    row_names_gp = gpar(fontsize = 2),
                    show_heatmap_legend = T,
                    # heatmap_ = unit(10),
                    heatmap_height = unit(8, 'cm'),
                    heatmap_width  = unit(10, 'cm'),
                    heatmap_legend_param = list(title_gp = gpar(fontsize = 10), 
                                                labels_gp = gpar(fontsize = 10),
                                                direction = "horizontal")) #

Hm_SGO_1
# size=4x4
draw(Hm_SGO_1, heatmap_legend_side = "bottom")
draw(Hm_SGO_2, heatmap_legend_side = "bottom")


# Figure WC for example
GO_table_bHLH91 <- GO_table$Zm00001d047017

Get_WC((c(GO_table_bHLH91$Term.GRN, GO_table_bHLH91$Term.CEN)))


####################################################


####################################################
######               Old plot                #######
####################################################

PDI.class <- subset(PDI.class, Method!="newDAP")
PDI.class[PDI.class$Method=='ChIP',]
ReplaceName(PDI.class[PDI.class$Method=='ChIP',]$TF)
ReplaceName(Table_rGss_Reasult[Table_rGss_Reasult$Method=='ChIP',]$TF)

DF_SS_PDI_CoExp <- left_join(DF_SS_PDI_CoExp, PDI.class, by=c("Source"="TF"))
DF_SS_PDI_CoExp$Method <- factor(DF_SS_PDI_CoExp$Method, levels = c("pChIP", "ChIP", "DAP"))


Plot_DF_SS_PDI_CoExp <- ggplot(DF_SS_PDI_CoExp, aes(x=Method, y=PDI_CoExp_SS))+
  geom_jitter(size=2, alpha=0.9, width = 0.2) + 
  dark_theme_bw()
Plot_DF_SS_PDI_CoExp <- ggpar(Plot_DF_SS_PDI_CoExp, font.tickslab = 14, 
                           font.x = 14, font.y = 14)

Plot_DF_SS_PDI_CoExp
subset(Table_rGss_Reasult, Pval <=0.05)
subset(Table_rGss_Reasult, FDR <= 0.1)
####################################################


####################################################
## Observed vs random SS: Example plot MYB31  ######
####################################################

## Number of GOs by BP data type
#
BP_PDIs[BP_PDIs$Mutant=="Zm00001d006236",]
#
BP_CoExp[BP_CoExp$Mutant=="Zm00001d006236",]

# observed value
DF_SS_PDI_CoExp[DF_SS_PDI_CoExp$Source=='Zm00001d006236',]
#DF_SS_PDI_CoExp[DF_SS_PDI_CoExp$Source=='Zm00001d000184',]

# Random values
MYB31_r.gss <- Random_pdi_coExp_ss$Zm00001d006236


R.GSS_MYB31_example <- ggplot(tibble(gss=MYB31_r.gss), aes(x=gss))+
  geom_histogram(fill='grey', color = 'grey', binwidth = 0.01)  +
  coord_flip() +
  scale_y_continuous(expand = c(0,0)) +
  geom_vline(xintercept=0.736, color="chocolate1", size=1.5, linetype='dashed') +
  xlab("GSS PDI vs Random GOs") +
  ylab("Counts") +
  theme_pubclean()
#dark_theme_bw()

R.GSS_MYB31_example <- ggpar(R.GSS_MYB31_example, font.tickslab = 14, 
                             font.x = 14, font.y = 14)
# size 4x5
R.GSS_MYB31_example

Ztest(0.736, MYB31_r.gss)

values_obs <- c(0.75 , 0.622)

####################################################

####################################################
## Identification of top GOs Terms: Heatmap       ##           
## MYB31 Example                                  ##
####################################################

subset(Table_rGss_Reasult, -log10(FDR)>=3)

MYB31_GO_matrix <- Get_SS_PDI_CoExp_Matrix(c("Zm00001d006236"))
dim(MYB31_GO_matrix)
nrow(MYB31_GO_matrix)

mask <- colnames(MYB31_GO_matrix)[colnames(MYB31_GO_matrix) %in% row.names(MYB31_GO_matrix)]
length(mask)

dim(MYB31_GO_matrix)


Hm_SGO <- Heatmap(MYB31_GO_matrix[mask,mask], name="GSS",
                         # column_km = 3,
                         # column_names_rot = 90,
                         # row_names_rot = 45,
                         cluster_rows = TRUE, 
                         cluster_columns = TRUE,
                         show_column_dend = FALSE, 
                          show_row_dend = FALSE, 
                         clustering_method_columns = "ward.D2",
                         col=viridis(100, direction = 1, option = "A"),
                         column_names_gp = gpar(fontsize = 5),
                         row_names_gp = gpar(fontsize = 5),
                         show_heatmap_legend = T,
                         # heatmap_ = unit(10),
                         heatmap_height = unit(30, 'cm'),
                         heatmap_width  = unit(30, 'cm'),
                         heatmap_legend_param = list(title_gp = gpar(fontsize = 10), 
                                                     labels_gp = gpar(fontsize = 10), 
                                                     direction = "horizontal"))
Hm_SGO
# size=3x10
draw(Hm_SGO, heatmap_legend_side = "bottom")




subset(BP_PDIs, 
       Mutant=="Zm00001d006236" & GO.ID %in% c("GO:0015711", "GO:0006865", "GO:0046942", "GO:0006835", 
                                                        "GO:0015800", "GO:0043090", "GO:0015804", "GO:0015807", 
                                                        "GO:0015812", "GO:0015824", "GO:0015912")) [,-c(1,8)]


ReplaceName(names(Random_pdi_coExp_ss))
# GO:0033993 response to lipid # PDI
# GO:0055088 lipid homeostasis # CoExp
# GO:0006554 lysine catabolic process # CoExp

TFdic <- as_tibble(read.table("Data/Annotations/TF_Id_Name.txt", h=T, stringsAsFactors = F))



