ReplaceNamePWY <- function(ids){
  # 
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

chop=function(myStr,mySep,myField){
  ## chop a string by a separator and return specified field
  choppedString=sapply(strsplit(myStr,mySep),"[",myField)
  if(length(myField)>1){
    choppedString=apply(choppedString,2,function(x){paste0(x[!is.na(x)],collapse=mySep)})
  }
  return(choppedString)
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
    geom_sf(aes(color=name), size = 1.5, 
            data = venn_setedge(data), show.legend = F) +
    #
    geom_sf_text(aes(label = name), size=6, data = venn_setlabel(data)) +
    geom_sf_text(aes(label= scales::comma(count, accuracy = 1)), size=5,
                 data = venn_region(data)) +
    #
    scale_fill_manual(values = col2) + # 
    scale_color_manual(values = alpha(col, .5)) +
    #
    theme_void() +
    theme(plot.margin = unit(c(0.5,1,1,0.1), "cm"),
          text = element_text(family="Helvetica")) +
    xlim(-150,1000)
}

vennfuncInt <- function(list_x){
  venn <- Venn(list_x)
  data <- process_data(venn)
  
  colorGroups <- c(CEN = '#FFD700', GRN='#6A5ACD', GAN='#1E90FF', eGRN='#FF1493')
  colfunc <- colorRampPalette(colorGroups)
  col <- colfunc(4)
  
  colorGroups <- c(CEN="gray100",GRN="gray99", GAN="gray98", eGRN='gray97')
  colfunc2 <- colorRampPalette(colorGroups)
  col2 <- colfunc2(15)
  
  ggplot() +
    geom_sf(aes(fill=name), data = venn_region(data), show.legend = F) +
    geom_sf(aes(color=name), size = 1.5,
            data = venn_setedge(data), show.legend = F) +
    #
    geom_sf_text(aes(label = name), size=4, data = venn_setlabel(data)) +
    geom_sf_text(aes(label= scales::comma(count, accuracy = 1)), size=3, data = venn_region(data)) +
    #
    scale_fill_manual(values = col2) + # 
    scale_color_manual(values = alpha(col, .5)) +
    #
    theme_void() +
    theme(plot.margin = unit(c(0.5,1,1,0.1), "cm"), 
          text=element_text(family="Helvetica")) #+
  #xlim(-150,1000)
  #text = element_text(family="Helvetica")
}

## PYW Enrichment 
Enrichmet_classes <- function(network){
  ## Count TF targets in network
  # Count Total
  network <- subset(network, Target %in% Syntenic)
  #
  Total_targtes <- as_tibble(as.data.frame(table(unique(network[,c("TF", "Target")])$TF)))
  colnames(Total_targtes) <- c('TF', 'targets') 
  Total_targtes$TF <- as.character(Total_targtes$TF)
  
  # list input: network
  network.list <- unique(network[,c("TF", "Target")])
  network.list <- split(network.list$Target, network.list$TF)
  
  ## Compare list of predicted targets vs annotated genes in query.vector 
  #
  print(". Pre-newGOM .")
  go.obj <- newGOM(network.list, CornCYC.list, genome.size=length(Syntenic)) # annotated genes in Genome v4
  
  Pval <- getMatrix(go.obj, name="pval")
  Common <- getMatrix(go.obj, name="intersection")
  
  print(". Post-newGOM .")
  
  ### Summary tables
  ## adjust p value
  Pval <- as.data.frame(Pval)
  #Pval[,1:ncol(Pval)] <- apply(Pval[,1:ncol(Pval)], 2, p.adjust)
  
  Pval_table   <- as_tibble(reshape2::melt(as.matrix(Pval)))
  colnames(Pval_table) <- c('TF', 'PWY', 'Pval')
  
  #
  Common_table <- as_tibble(reshape2::melt(as.matrix(Common)))
  colnames(Common_table) <- c('TF', 'PWY', 'n.targ')
  
  # Add predicted target in class by TF
  Pval_table <- left_join(Pval_table, Common_table , by=c('TF', 'PWY'))
  
  # Add total predicted targets
  Pval_table <- left_join(Pval_table, Total_targtes, by="TF")
  
  #write.table(Pval_table, name, sep = "\t", row.names = F, quote = F)
  print(paste("... Done ...", sep = ""))
  return(Pval_table)
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

GetGO <- function(degs, mutant, net){
  
  # Use a list of DEGs and the name of the mutant (string)
  # to identify GOs enriched. Required to have a background predefined
  #
  # mutant: description of a category of source of the data
  # degs: set of genes to test enrichment
  print(length(degs))
  GeneList <- factor(as.integer(background_IDs %in% degs))
  names(GeneList) <- background_IDs
  
  GOdata_BP <- new("topGOdata", ontology = "BP", 
                   allGenes = GeneList, annot = annFUN.gene2GO, 
                   gene2GO = background)
  
  #GOdata_MF <- new("topGOdata", ontology = "MF", 
  # allGenes = GeneList, annot = annFUN.gene2GO, gene2GO = background)
  #GOdata_CC <- new("topGOdata", ontology = "CC", 
    # allGenes = GeneList, annot = annFUN.gene2GO, gene2GO = background)
  
  #### Define test ####
  test.stat <- new("classicCount", testStatistic = GOFisherTest, 
                   name = "Fisher test", nodeSize = 10)
  
  ### test enrichment 
  results_BP <- getSigGroups(GOdata_BP, test.stat)
  
  ### save pdf Graph
  #namepdf=paste("GOs_Plots/GO.BP_",mutant, "", sep = "")
  #printGraph(GOdata_BP, results_BP, firstSigNodes=20,  fn.prefix = namepdf, useInfo = "def", pdfSW = TRUE) #
  
  ######## Get Significant GOs ########  
  Res_DF_BP <- as_tibble(GenTable(GOdata_BP, classic = results_BP, topNodes = 1000, orderBy='Fis')) # save as dataframe
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
  
  DF_GO_Genes <- unique(DF_GO_Genes)
  filename <- 
    write.table(DF_GO_Genes,
                paste("BP_results/Genes_GOBP_", net, "_", mutant, ".txt", sep = ""),
                sep = '\t', quote = F,
                row.names = F)
  
  #
  Res_DF_BP["Type"] <- net
  #return(list(GOs_DF=Res_DF_BP, Genes=DF_GO_Genes)) # Return list of GOs-Stats and GeneID-GOs
  return(Res_DF_BP) # Return list of GOs-Stats and GeneID-GOs
}

SuperGO <- function(tf, network){
  ## used TF's targets to test GO terms enrichment
  # 1. Select tf's Targets
  # 2. Filter out Targets predicted by only a layers
  # 3. Add up targets predicted by at least two layers
  # 4. test enrichment
  
  # Get network by TF
  if (network=='GRN') {
    cat(' .. path grn ..')
    targets <- unique(PDI[PDI$TF==tf,]$Target)
    targets <- targets[targets %in% Syntenic]
    #
    Total_targets <- length(targets)
  }
  
  else if (network=='eGRN') {
    cat(' .. path grn ..')
    targets <- unique(PDIeQTL[PDIeQTL$TF==tf,]$Target)
    targets <- targets[targets %in% Syntenic]
    #
    Total_targets <- length(targets)
  }
  
  else if (network=='CEN') {
    cat(' .. path grn ..')
    targets <- unique(CoExp[CoExp$TF==tf,]$Target)
    targets <- targets[targets %in% Syntenic]
    #
    Total_targets <- length(targets)
  }
  else if (network=='GAN') {
    cat(' .. path grn ..')
    targets <- unique(teQTL[teQTL$TF==tf,]$Target)
    targets <- targets[targets %in% Syntenic]
    #
    Total_targets <- length(targets)
  }
  else{ print('.. salado .. ') }
  
  if(Total_targets > 5){
    # If targets largert than 
    # Include 3_later targets into paired comparison if required 
    out <- GetGO(targets, tf, network)
    # Save gene results GOs  
    write.table(out, paste0("BP_results/GOsBP_", network, "_",tf, ".txt"), 
                sep = '\t', row.names = F, quote = F)
  }
  else { print(paste0(" Salado: ", tf, " ..")) }
  return(out)
  
}

Get_SS_GRN_CEN <- function(tfslist){
  # make df to save output
  GOSIM_DF <- tibble(Source=tfslist, SS=0, Net="GRN:CEN")
  
  c= 1
  for (tf in tfslist){
    # PDI
    go1 <- subset(GOsDB_GRN, Mutant == tf)$GO.ID
    # CoExp
    go2 <- subset(GOsDB_CEN, Mutant == tf)$GO.ID
    #
    ss <- mgoSim(go1, go2, semData=Zm.GOSemSim.BP, measure="Wang", combine="BMA")
    GOSIM_DF$SS[GOSIM_DF$Source==tf]  <- ss
    
    # Matrix of GOs
    M_ss <- mgoSim(go1, go2, semData=Zm.GOSemSim.BP, measure="Wang", combine=NULL)
    M_ss <- melt(M_ss)
    colnames(M_ss) <- c("GO_GRN", "GO_CEN", "GSS") 
    
    BP_GRN_test <- subset(GOsDB_GRN, Mutant==tf & GO.ID %in% M_ss$GO_GRN)
    BP_CEN_test <- subset(GOsDB_CEN, Mutant==tf & GO.ID %in% M_ss$GO_CEN)
    
    M_ss <- left_join(M_ss, BP_GRN_test[,c(1,2,4,9)], by=c("GO_GRN"="GO.ID"))
    M_ss <- left_join(M_ss, BP_CEN_test[,c(1,2,4,9)], by=c("GO_CEN"="GO.ID"))
    colnames(M_ss) <- c("GO.GRN", "GO.CEN", "GSS", 
                        "Term.GRN",  "Sig.Terms.GRN", "FDR.GRN",
                        "Term.CEN", "Sig.Terms.CEN", "FDR.CEN")
    
    namef <- paste0("GO_SS_data/GSS_",tf,".GRN_CEN.txt")
    cat(namef)
    cat('\n')
    write.table(M_ss, namef, row.names = F, sep = '\t', quote = F)
    print(paste0(" .. Done ", c, " TFs .."))
    c = c+1
  }
  return(GOSIM_DF) 
}

Get_SS_GRN_GAN <- function(tfslist){
  # make df to save output
  GOSIM_DF <- tibble(Source=tfslist, SS=0, Net="GRN:GAN")
  
  c= 1
  for (tf in tfslist){
    # PDI
    go1 <- subset(GOsDB_GRN, Mutant == tf)$GO.ID
    # CoExp
    go2 <- subset(GOsDB_GAN, Mutant == tf)$GO.ID
    #
    ss <- mgoSim(go1, go2, semData=Zm.GOSemSim.BP, measure="Wang", combine="BMA")
    GOSIM_DF$SS[GOSIM_DF$Source==tf]  <- ss
    
    # Matrix of GOs
    M_ss <- mgoSim(go1, go2, semData=Zm.GOSemSim.BP, measure="Wang", combine=NULL)
    M_ss <- melt(M_ss)
    colnames(M_ss) <- c("GO_GRN", "GO_GAN", "GSS") 
    
    BP_GRN_test <- subset(GOsDB_GRN, Mutant==tf & GO.ID %in% M_ss$GO_GRN)
    BP_GAN_test <- subset(GOsDB_GAN, Mutant==tf & GO.ID %in% M_ss$GO_GAN)
    
    M_ss <- left_join(M_ss, BP_GRN_test[,c(1,2,4,9)], by=c("GO_GRN"="GO.ID"))
    M_ss <- left_join(M_ss, BP_GAN_test[,c(1,2,4,9)], by=c("GO_GAN"="GO.ID"))
    colnames(M_ss) <- c("GO.GRN", "GO.GAN", "GSS", 
                        "Term.GRN", "Sig.Terms.GRN", "FDR.GRN",
                        "Term.GAN", "Sig.Terms.GAN", "FDR.GAN")
    
    namef <- paste0("GO_SS_data/GSS_",tf,".GRN_GAN.txt")
    cat(namef)
    cat('\n')
    write.table(M_ss, namef, row.names = F, sep = '\t', quote = F)
    print(paste0(" .. Done ", c, " TFs .."))
    c = c+1
  }
  return(GOSIM_DF) 
}
##
Get_SS_GRN_eGRN <- function(tfslist){
  # make df to save output
  GOSIM_DF <- tibble(Source=tfslist, SS=0, Net="GRN:eGRN")
  
  c= 1
  for (tf in tfslist){
    # PDI
    go1 <- subset(GOsDB_GRN, Mutant == tf)$GO.ID
    # CoExp
    go2 <- subset(GOsDB_eGRN, Mutant == tf)$GO.ID
    
    ss <- mgoSim(go1, go2, semData=Zm.GOSemSim.BP, measure="Wang", combine="BMA")    
    GOSIM_DF$SS[GOSIM_DF$Source==tf]  <- ss
    
    # Matrix of GOs
    M_ss <- mgoSim(go1, go2, semData=Zm.GOSemSim.BP, measure="Wang", combine=NULL)
    M_ss <- melt(M_ss)
    colnames(M_ss) <- c("GO_GRN", "GO_eGRN", "GSS")
    
    BP_GRN_test <- subset(GOsDB_GRN, Mutant==tf & GO.ID %in% M_ss$GO_GRN)
    BP_eGRN_test <- subset(GOsDB_eGRN, Mutant==tf & GO.ID %in% M_ss$GO_eGRN)
    
    M_ss <- left_join(M_ss, BP_GRN_test[,c(1,2,4,9)],  by=c("GO_GRN"="GO.ID"))
    M_ss <- left_join(M_ss, BP_eGRN_test[,c(1,2,4,9)], by=c("GO_eGRN"="GO.ID"))
    colnames(M_ss) <- c("GO.GRN", "GO.eGRN", "GSS", 
                        "Term.GRN",  "Sig.Terms.GRN", "FDR.GRN",
                        "Term.eGRN", "Sig.Terms.eGRN", "FDR.eGRN")
    
    namef <- paste0("GO_SS_data/GSS_",tf,".GRN_eGRN.txt")
    cat(namef)
    cat('\n')
    write.table(M_ss, namef, row.names = F, sep = '\t', quote = F)
    
    print(paste0(" .. Done ", c, " TFs .."))
    c = c+1
    
  }
  return(GOSIM_DF)
  
}

Get_SS_CEN_GAN <- function(tfslist){
  # make df to save output
  GOSIM_DF <- tibble(Source=tfslist, SS=0, Net="CEN:GAN")
  
  c= 1
  for (tf in tfslist){
    # PDI
    go1 <- subset(GOsDB_CEN, Mutant == tf)$GO.ID
    
    # CoExp
    go2 <- subset(GOsDB_GAN, Mutant == tf)$GO.ID
    
    ss <- mgoSim(go1, go2, semData=Zm.GOSemSim.BP, measure="Wang", combine="BMA")
    GOSIM_DF$SS[GOSIM_DF$Source==tf]  <- ss
    
    # Matrix of GOs
    M_ss <- mgoSim(go1, go2, semData=Zm.GOSemSim.BP, measure="Wang", combine=NULL)
    M_ss <- melt(M_ss)
    colnames(M_ss) <- c("GO_CEN", "GO_GAN", "GSS")
    
    BP_CEN_test <- subset(GOsDB_CEN, Mutant==tf & GO.ID %in% M_ss$GO_GRN)
    BP_GAN_test <- subset(GOsDB_GAN, Mutant==tf & GO.ID %in% M_ss$GO_GAN)
    
    M_ss <- left_join(M_ss, BP_CEN_test[,c(1,2,4,9)], by=c("GO_CEN"="GO.ID"))
    M_ss <- left_join(M_ss, BP_GAN_test[,c(1,2,4,9)], by=c("GO_GAN"="GO.ID"))
    colnames(M_ss) <- c("GO.CEN", "GO.GAN", "GSS", 
                        "Term.CEN", "Sig.Terms.CEN", "FDR.CEN",
                        "Term.GAN", "Sig.Terms.GAN", "FDR.GAN")
    
    namef <- paste0("GO_SS_data/GSS_",tf,".CEN_GAN.txt")
    cat(namef)
    cat('\n')
    write.table(M_ss, namef, row.names = F, sep = '\t', quote = F)  
      
    print(paste0(" .. Done ", c, " TFs .."))
    c = c+1
    
  }
  return(GOSIM_DF)
  
}

Get_SS_eGRN_CEN <- function(tfslist){
  # make df to save output
  GOSIM_DF <- tibble(Source=tfslist, SS=0, Net="eGRN:CEN")
  
  c= 1
  for (tf in tfslist){
    # PDI
    go1 <- subset(GOsDB_eGRN, Mutant == tf)$GO.ID
    # CoExp
    go2 <- subset(GOsDB_CEN, Mutant == tf)$GO.ID
    
    ss <- mgoSim(go1, go2, semData=Zm.GOSemSim.BP, measure="Wang", combine="BMA")    
    GOSIM_DF$SS[GOSIM_DF$Source==tf]  <- ss
    
    # Matrix of GOs
    M_ss <- mgoSim(go1, go2, semData=Zm.GOSemSim.BP, measure="Wang", combine=NULL)
    M_ss <- melt(M_ss)
    colnames(M_ss) <- c("GO_eGRN", "GO_CEN", "GSS")
    
    BP_eGRN_test <- subset(GOsDB_eGRN, Mutant==tf & GO.ID %in% M_ss$GO_eGRN)
    BP_CEN_test <- subset(GOsDB_CEN, Mutant==tf & GO.ID %in% M_ss$GO_CEN)
    
    M_ss <- left_join(M_ss, BP_eGRN_test[,c(1,2,4,9)], by=c("GO_eGRN"="GO.ID"))
    M_ss <- left_join(M_ss, BP_CEN_test[,c(1,2,4,9)], by=c("GO_CEN"="GO.ID"))
    colnames(M_ss) <- c("GO.eGRN", "GO.CEN", "GSS", 
                        "Term.eGRN", "Sig.Terms.eGRN", "FDR.eGRN",
                        "Term.CEN", "Sig.Terms.CEN", "FDR.CEN")
    
    namef <- paste0("GO_SS_data/GSS_",tf,".eGRN_CEN.txt")
    cat(namef)
    cat('\n')
    write.table(M_ss, namef, row.names = F, sep = '\t', quote = F)
    print(paste0(" .. Done ", c, " TFs .."))
    c = c+1
    
  }
  return(GOSIM_DF)
  
}

Get_SS_GAN_eGRN <- function(tfslist){
  # make df to save output
  GOSIM_DF <- tibble(Source=tfslist, SS=0, Net="GAN:eGRN")
  
  c= 1
  for (tf in tfslist){
    # PDI
    go1 <- subset(GOsDB_GAN, Mutant == tf)$GO.ID
    # CoExp
    go2 <- subset(GOsDB_eGRN, Mutant == tf)$GO.ID
    
    ss <- mgoSim(go1, go2, semData=Zm.GOSemSim.BP, measure="Wang", combine="BMA")    
    GOSIM_DF$SS[GOSIM_DF$Source==tf]  <- ss
    
    # Matrix of GOs
    M_ss <- mgoSim(go1, go2, semData=Zm.GOSemSim.BP, measure="Wang", combine=NULL)
    M_ss <- melt(M_ss)
    colnames(M_ss) <- c("GO_GAN", "GO_eGRN", "GSS")
    
    BP_GAN_test <- subset(GOsDB_GAN, Mutant==tf & GO.ID %in% M_ss$GO_GAN)
    BP_eGRN_test <- subset(GOsDB_eGRN, Mutant==tf & GO.ID %in% M_ss$GO_eGRN)
    
    M_ss <- left_join(M_ss, BP_GAN_test[,c(1,2,4,9)], by=c("GO_GAN"="GO.ID"))
    M_ss <- left_join(M_ss, BP_eGRN_test[,c(1,2,4,9)], by=c("GO_eGRN"="GO.ID"))
    colnames(M_ss) <- c("GO.GAN", "GO.eGRN", "GSS", 
                        "Term.GAN", "Sig.Terms.GAN", "FDR.GAN",
                        "Term.eGRN", "Sig.Terms.eGRN", "FDR.eGRN")
    
    namef <- paste0("GO_SS_data/GSS_",tf,".GAN_eGRN.txt")
    cat(namef)
    cat('\n')
    write.table(M_ss, namef, row.names = F, sep = '\t', quote = F)
    print(paste0(" .. Done ", c, " TFs .."))
    c = c+1
    
  }
  return(GOSIM_DF)
  
}
                     
#
Read_GO_ssTF <- function(file){
  # Read GO_SS file from GO_SS_data
  name <- gsub('.txt', '', gsub('GSS_', '', file))
  namedf <- chop(name, '[.]',1)
  namenet <- gsub('_',':', chop(name, '[.]',2))
  
  file <- fread(paste0("GO_SS_data/", file))
  colnames(file) <- c("GO1","GO2", "GSS", "Term.1", "Sig.Term.1", "FDR.1",
                      "Term.2", "Sig.Term.2", "FDR.2")
  file[,'TF'] <- namedf
  file[,'Nets'] <- namenet
  
  return(file)
  
}

# compare TFs by network pairs
tf="Zm00001d005160"

Func.CEN_GAN_PWY <- function (tf){
  # pwys from cen
  pwy1 <- CEN_PWY[CEN_PWY$TF == tf,]$PWY
  pwy2 <- GAN_PWY[GAN_PWY$TF == tf,]$PWY
  
  # Get vals from PWY similarity matrix
  SIMm <- subset(PWY_similarity, PWY1 %in% pwy1 & PWY2 %in% pwy2)
  colnames(SIMm)[c(1,2, 4)] <- c("PWYcen", "PWYgan", "CommonGenesPWYs")
  
  # make amalysis if at least a common gene between PWYs
  SIMm <- subset(SIMm, CommonGenesPWYs >= 1)
  
  if (nrow(SIMm) >0) {
    # Count targets in original CYC list 
    cen.target <- unique(CoExp[CoExp$TF == tf,]$Target)
    gan.target <- unique(teQTL[teQTL$TF == tf,]$Target)
    
    # PWYs targets
    CornCYC_target <- subset(CornCYC, Pathway.id %in% c(pwy1, pwy2))
    CornCYC_target <- subset(CornCYC_target, GeneID %in% Syntenic)
    
    # PWY size
    CornCYC_target_size <- as.data.table(table(CornCYC_target$Pathway.id))
    
    # Add PWY size to Sim matrix
    SIMm <- left_join(SIMm, CornCYC_target_size, by=c("PWYcen"="V1"))
    SIMm <- left_join(SIMm, CornCYC_target_size, by=c("PWYgan"="V1"))
    colnames(SIMm)[c(5,6)] <- c("PWYcenSize", "PWYganSize")
    
    # targets in PWYs
    CornCYC_target_1 <- unique(subset(CornCYC_target, GeneID %in% cen.target))
    CornCYC_target_2 <- unique(subset(CornCYC_target, GeneID %in% gan.target))
    
    # Counts targets in PWY
    CornCYC_target_1 <- as.data.table(table(CornCYC_target_1$Pathway.id))
    colnames(CornCYC_target_1) <- c("PWYcen", "Targ.PWYcen")
    
    CornCYC_target_2 <- as.data.table(table(CornCYC_target_2$Pathway.id))
    colnames(CornCYC_target_2) <- c("PWYgan", "Targ.PWYgan")
    
    # Add targets information to similarity matrix
    SIMm <- left_join(SIMm, CornCYC_target_1, by='PWYcen')
    SIMm <- left_join(SIMm, CornCYC_target_2, by='PWYgan')
    
    SIMm[,'PWYcen.name'] <- ReplaceNamePWY(SIMm$PWYcen)
    SIMm[,'PWYgan.name'] <- ReplaceNamePWY(SIMm$PWYgan)
    
    return(SIMm)
  }
  else{
    cat('.. Salado ..')
  }
}

Func.CEN_eGRN_PWY <- function (tf){
  # pwys from cen
  pwy1 <- CEN_PWY[CEN_PWY$TF == tf,]$PWY
  pwy2 <- eGRN_PWY[eGRN_PWY$TF == tf,]$PWY
  
  # Get vals from PWY similarity matrix
  SIMm <- subset(PWY_similarity, PWY1 %in% pwy1 & PWY2 %in% pwy2)
  colnames(SIMm)[c(1,2, 4)] <- c("PWYcen", "PWYegrn", "CommonGenesPWYs")
  
  # make amalysis if at least a common gene between PWYs
  SIMm <- subset(SIMm, CommonGenesPWYs >= 1)
  
  if (nrow(SIMm) >0) {
    # Count targets in original CYC list 
    cen.target <- unique(CoExp[CoExp$TF == tf,]$Target)
    egrn.target <- unique(PDIeQTL[PDIeQTL$TF == tf,]$Target)
    
    # PWYs targets
    CornCYC_target <- subset(CornCYC, Pathway.id %in% c(pwy1, pwy2))
    CornCYC_target <- subset(CornCYC_target, GeneID %in% Syntenic)
    
    # PWY size
    CornCYC_target_size <- as.data.table(table(CornCYC_target$Pathway.id))
    
    # Add PWY size to Sim matrix
    SIMm <- left_join(SIMm, CornCYC_target_size, by=c("PWYcen"="V1"))
    SIMm <- left_join(SIMm, CornCYC_target_size, by=c("PWYegrn"="V1"))
    colnames(SIMm)[c(5,6)] <- c("PWYcenSize", "PWYegrnSize")
    
    # targets in PWYs
    CornCYC_target_1 <- unique(subset(CornCYC_target, GeneID %in% cen.target))
    CornCYC_target_2 <- unique(subset(CornCYC_target, GeneID %in% egrn.target))
    
    # Counts targets in PWY
    CornCYC_target_1 <- as.data.table(table(CornCYC_target_1$Pathway.id))
    colnames(CornCYC_target_1) <- c("PWYcen", "Targ.PWYcen")
    
    CornCYC_target_2 <- as.data.table(table(CornCYC_target_2$Pathway.id))
    colnames(CornCYC_target_2) <- c("PWYegrn", "Targ.PWYegrn")
    
    # Add targets information to similarity matrix
    SIMm <- left_join(SIMm, CornCYC_target_1, by='PWYcen')
    SIMm <- left_join(SIMm, CornCYC_target_2, by='PWYegrn')
    
    SIMm[,'PWYcen.name'] <- ReplaceNamePWY(SIMm$PWYcen)
    SIMm[,'PWYegrn.name'] <- ReplaceNamePWY(SIMm$PWYegrn)
    
    return(SIMm)
  }
  else{
    cat('.. Salado ..')
  }
}

Func.GAN_eGRN_PWY <- function (tf){
  # pwys from cen
  pwy1 <- GAN_PWY[GAN_PWY$TF == tf,]$PWY
  pwy2 <- eGRN_PWY[eGRN_PWY$TF == tf,]$PWY
  
  # Get vals from PWY similarity matrix
  SIMm <- subset(PWY_similarity, PWY1 %in% pwy1 & PWY2 %in% pwy2)
  colnames(SIMm)[c(1,2, 4)] <- c("PWYgan", "PWYegrn", "CommonGenesPWYs")
  
  # make amalysis if at least a common gene between PWYs
  SIMm <- subset(SIMm, CommonGenesPWYs >= 1)
  
  if (nrow(SIMm) >0) {
    # Count targets in original CYC list 
    gan.target <- unique(teQTL[teQTL$TF == tf,]$Target)
    egrn.target <- unique(PDIeQTL[PDIeQTL$TF == tf,]$Target)
    
    # PWYs targets
    CornCYC_target <- subset(CornCYC, Pathway.id %in% c(pwy1, pwy2))
    CornCYC_target <- subset(CornCYC_target, GeneID %in% Syntenic)
    
    # PWY size
    CornCYC_target_size <- as.data.table(table(CornCYC_target$Pathway.id))
    
    # Add PWY size to Sim matrix
    SIMm <- left_join(SIMm, CornCYC_target_size, by=c("PWYgan"="V1"))
    SIMm <- left_join(SIMm, CornCYC_target_size, by=c("PWYegrn"="V1"))
    colnames(SIMm)[c(5,6)] <- c("PWYganSize", "PWYegrnSize")
    
    # targets in PWYs
    CornCYC_target_1 <- unique(subset(CornCYC_target, GeneID %in% gan.target))
    CornCYC_target_2 <- unique(subset(CornCYC_target, GeneID %in% egrn.target))
    
    # Counts targets in PWY
    CornCYC_target_1 <- as.data.table(table(CornCYC_target_1$Pathway.id))
    colnames(CornCYC_target_1) <- c("PWYgan", "Targ.PWYgan")
    
    CornCYC_target_2 <- as.data.table(table(CornCYC_target_2$Pathway.id))
    colnames(CornCYC_target_2) <- c("PWYegrn", "Targ.PWYegrn")
    
    # Add targets information to similarity matrix
    SIMm <- left_join(SIMm, CornCYC_target_1, by='PWYgan')
    SIMm <- left_join(SIMm, CornCYC_target_2, by='PWYegrn')
    
    SIMm[,'PWYgan.name'] <- ReplaceNamePWY(SIMm$PWYgan)
    SIMm[,'PWYegrn.name'] <- ReplaceNamePWY(SIMm$PWYegrn)
    
    return(SIMm)
  }
  else{
    cat('.. Salado ..')
  }
}

Func.GAN_GRN_PWY <- function (tf){
  # pwys from cen
  pwy1 <- GAN_PWY[GAN_PWY$TF == tf,]$PWY
  pwy2 <- GRN_PWY[GRN_PWY$TF == tf,]$PWY
  
  # Get vals from PWY similarity matrix
  SIMm <- subset(PWY_similarity, PWY1 %in% pwy1 & PWY2 %in% pwy2)
  colnames(SIMm)[c(1,2, 4)] <- c("PWYgan", "PWYgrn", "CommonGenesPWYs")
  
  # make amalysis if at least a common gene between PWYs
  SIMm <- subset(SIMm, CommonGenesPWYs >= 1)
  
  if (nrow(SIMm) >0) {
    # Count targets in original CYC list 
    gan.target <- unique(teQTL[teQTL$TF == tf,]$Target)
    grn.target <- unique(PDI[PDI$TF == tf,]$Target)
    
    # PWYs targets
    CornCYC_target <- subset(CornCYC, Pathway.id %in% c(pwy1, pwy2))
    CornCYC_target <- subset(CornCYC_target, GeneID %in% Syntenic)
    
    # PWY size
    CornCYC_target_size <- as.data.table(table(CornCYC_target$Pathway.id))
    
    # Add PWY size to Sim matrix
    SIMm <- left_join(SIMm, CornCYC_target_size, by=c("PWYgan"="V1"))
    SIMm <- left_join(SIMm, CornCYC_target_size, by=c("PWYgrn"="V1"))
    colnames(SIMm)[c(5,6)] <- c("PWYganSize", "PWYgrnSize")
    
    # targets in PWYs
    CornCYC_target_1 <- unique(subset(CornCYC_target, GeneID %in% gan.target))
    CornCYC_target_2 <- unique(subset(CornCYC_target, GeneID %in% grn.target))
    
    # Counts targets in PWY
    CornCYC_target_1 <- as.data.table(table(CornCYC_target_1$Pathway.id))
    colnames(CornCYC_target_1) <- c("PWYgan", "Targ.PWYgan")
    
    CornCYC_target_2 <- as.data.table(table(CornCYC_target_2$Pathway.id))
    colnames(CornCYC_target_2) <- c("PWYgrn", "Targ.PWYgrn")
    
    # Add targets information to similarity matrix
    SIMm <- left_join(SIMm, CornCYC_target_1, by='PWYgan')
    SIMm <- left_join(SIMm, CornCYC_target_2, by='PWYgrn')
    
    SIMm[,'PWYgan.name'] <- ReplaceNamePWY(SIMm$PWYgan)
    SIMm[,'PWYgrn.name'] <- ReplaceNamePWY(SIMm$PWYgrn)
    
    return(SIMm)
  }
  else{
    cat('.. Salado ..')
  }
}

Func.GRN_eGRN_PWY <- function (tf){
  # pwys from cen
  pwy1 <- GRN_PWY[GRN_PWY$TF == tf,]$PWY
  pwy2 <- eGRN_PWY[eGRN_PWY$TF == tf,]$PWY
  
  # Get vals from PWY similarity matrix
  SIMm <- subset(PWY_similarity, PWY1 %in% pwy1 & PWY2 %in% pwy2)
  colnames(SIMm)[c(1,2, 4)] <- c("PWYgrn", "PWYegrn", "CommonGenesPWYs")
  
  # make amalysis if at least a common gene between PWYs
  SIMm <- subset(SIMm, CommonGenesPWYs >= 1)
  
  if (nrow(SIMm) >0) {
    # Count targets in original CYC list 
    grn.target <- unique(PDI[PDI$TF == tf,]$Target)
    egrn.target <- unique(PDIeQTL[PDIeQTL$TF == tf,]$Target)
    
    # PWYs targets
    CornCYC_target <- subset(CornCYC, Pathway.id %in% c(pwy1, pwy2))
    CornCYC_target <- subset(CornCYC_target, GeneID %in% Syntenic)
    
    # PWY size
    CornCYC_target_size <- as.data.table(table(CornCYC_target$Pathway.id))
    
    # Add PWY size to Sim matrix
    SIMm <- left_join(SIMm, CornCYC_target_size, by=c("PWYgrn"="V1"))
    SIMm <- left_join(SIMm, CornCYC_target_size, by=c("PWYegrn"="V1"))
    colnames(SIMm)[c(5,6)] <- c("PWYgrnSize", "PWYegrnSize")
    
    # targets in PWYs
    CornCYC_target_1 <- unique(subset(CornCYC_target, GeneID %in% grn.target))
    CornCYC_target_2 <- unique(subset(CornCYC_target, GeneID %in% egrn.target))
    
    # Counts targets in PWY
    CornCYC_target_1 <- as.data.table(table(CornCYC_target_1$Pathway.id))
    colnames(CornCYC_target_1) <- c("PWYgrn", "Targ.PWYgrn")
    
    CornCYC_target_2 <- as.data.table(table(CornCYC_target_2$Pathway.id))
    colnames(CornCYC_target_2) <- c("PWYegrn", "Targ.PWYegrn")
    
    # Add targets information to similarity matrix
    SIMm <- left_join(SIMm, CornCYC_target_1, by='PWYgrn')
    SIMm <- left_join(SIMm, CornCYC_target_2, by='PWYegrn')
    
    SIMm[,'PWYgrn.name'] <- ReplaceNamePWY(SIMm$PWYgrn)
    SIMm[,'PWYegrn.name'] <- ReplaceNamePWY(SIMm$PWYegrn)
    
    return(SIMm)
  }
  else{
    cat('.. Salado ..')
  }
}
