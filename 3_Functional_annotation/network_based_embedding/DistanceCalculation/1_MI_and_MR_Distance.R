# suppressMessages(library(propagate)): bitcor function
suppressMessages(library(reshape2))
library(parmigene)
suppressMessages(library(dplyr))
library(parallel)
suppressMessages(library(data.table))

#system('export OMP_NUM_THREADS=15') 
#OMP_NUM_THREADS=15

Syntenic <- as_tibble(read.table("../Data/Annotations/Zm.v4.synteny.genes.txt", h=T, stringsAsFactors = F, sep = '\t'))
Syntenic <- subset(Syntenic, gene_synteny=='syntenic')$gene_id

# Read pecanpy files
pecanpyM <- fread("Pecanpy_uFNetsW.Dim50_WL80_nW100.txt", h=F, skip = 1)
#pecanpyM <- fread("Pecanpy_uFNetsW.Dim100_WL100_nW50.txt", h=F, skip = 1)
#pecanpyM <- pecanpyM[1:100,]

# Distance for only syntenic genes
#pecanpyM <- subset(pecanpyM, V1 %in% Syntenic)

GeneIDs <- pecanpyM$V1
pecanpyM <- pecanpyM[,-c(1)]

print(paste0('.. Genes tested: ', length(GeneIDs)))

pecanpyM <- as.matrix(pecanpyM)


# add geneid as row names 
row.names(pecanpyM) <- GeneIDs 
#pecanpyM[1:5,1:5]


######## Using as input the MAtrix ##########
MI <- round(knnmi.all(pecanpyM, k=3, noise=1e-12), 4)
MI_Rank <- apply(-MI, 1, rank)            # Rank from smaller to bigger: large MI values (good) with "-" will be rank close to 1.
MR_MI <- sqrt(MI_Rank*t(MI_Rank))         # Mutual rank calcuation
MR_MI <- round(MR_MI, 2)

########################################


###### Save data #######
########################
#count=1
#for (i in GeneIDs){

SaveMI <- function(gid){

  # combined MI and MR 
  TemDF <- as.data.frame(cbind(MI=MI[,gid], MR=MR_MI[,gid]))
  
  # Order by MR value
  TemDF <- TemDF[order(TemDF$MR),]
  
  # Add gene ID
  TemDF <- cbind(GeneID1=gid, GeneID2=row.names(TemDF), TemDF)
  TemDF[,'edgeW25'] <- round(exp(-( (TemDF$MR-1)/25) ),4)
  TemDF[,'edgeW50'] <- round(exp(-( (TemDF$MR-1)/50) ),4)
  
  # MR_edgesDB_Dim50_WL80_nW10_syntenic
  #write.table(TemDF, paste("MR_edgesDB_Dim100_WL100_nW50_syntenic/MR_MI.pecanpy.",gid,".txt",sep=""), row.names = F, sep='\t', quote = F)
  #
  write.table(TemDF, paste("MR_edgesDB_Dim50_WL80_nW100/MR_MI.pecanpy.",gid,".txt",sep=""), row.names = F, sep='\t', quote = F)
  print(paste(".......Printed Gene ", gid, " .......",sep = ""))
  
}

GenesList <- unique(GeneIDs)
Lgenes <- length(GenesList)

w=50  # Size of range to test
print(".. Ready to start ..")

for (i in seq(0, Lgenes, w)){
  max=Lgenes
  Start=i+1
  end=i+w
  
  if (end<max){
    listtotest <- GenesList[Start:end]
    cat(" working on:", Start:end, "\n")
    mclapply(listtotest, function(x) SaveMI(x), mc.cores=w)
  }
  else{
    r <- paste(Start, max, sep = "-")
    listtotest <- GenesList[Start:max]
    cat(" working on:", Start:max, "\n")
    w <- max - (Start-1)
    mclapply(listtotest, function(x) SaveMI(x), mc.cores=w)
  }
}
