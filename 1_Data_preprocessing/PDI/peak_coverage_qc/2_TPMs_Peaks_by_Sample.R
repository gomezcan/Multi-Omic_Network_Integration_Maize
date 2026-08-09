# libraries requared
library(Rsubread)
suppressMessages(library(tidyverse))


args = commandArgs(trailingOnly=TRUE)


###
tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate/sum(rate) * 1e6
}


#

make_df_tpms <- function(counts){
  
  DF_Counts <- data.frame(counts$annotation[,c("GeneID","Strand","Length")], counts$counts, stringsAsFactors=FALSE)
  colnames(DF_Counts) <- gsub(".bam", "", colnames(DF_Counts))
  #
  DF_Description <- DF_Counts[,c(1,2,3)] # remove Geneid, strand and gene length from counts table
  DF_Counts <- DF_Counts[,4] # remove Geneid, strand and gene length from counts table
  #
  # tmp
  DF_Description$Length <- DF_Description$Length/1000
  DF_tpms <- tpm(DF_Counts, DF_Description$Length)
  #
  DF_Description$cpm <- round(DF_tpms, 3)
  DF_Description <- DF_Description[,-c(2,3)]
 
  colnames(DF_Description)[1] <- "Summit"
  #
  #DF_tpms$Summit <- paste("chr", DF_tpms$Summit, sep='')
  
  return(DF_Description)
}

# Lybrary layout
LibrariesLayout <- read.table("LibrariesLayout.txt", header = F, stringsAsFactors = F)


Single <- subset(LibrariesLayout, V2=="Single")$V1
Paired <- subset(LibrariesLayout, V2=="Paired")$V1

# get list of bam files singles from part of single reads
myfiles = list.files(pattern = "*.bam") 

myfiles.singles =  myfiles[myfiles %in% Single]
myfiles.paired =  myfiles[myfiles %in% Paired]

# myfiles.singles = c(myfiles.singles, c("Halo_deM_R1.bam", "Halo_deM_R2.bam")) 

print(".. ready to count ..")

############ count alingments single end ######################## 
#print(myfiles.singles)

for (i in 1:length(myfiles.singles)){
	
	name <- myfiles.singles[i]
	name <- gsub(".bam", "", name)
	Counts_SE <- featureCounts(files=myfiles.singles[i], isGTFAnnotationFile = FALSE, annot.ext="All.peaks.saf",
				countMultiMappingReads=FALSE, isPairedEnd=FALSE, useMetaFeatures=FALSE, nthreads=30, 
				allowMultiOverlap=TRUE, minOverlap=25)
	
	DF_Counts_SE <- make_df_tpms(Counts_SE)

	colnames(DF_Counts_SE)[2] <- name

	write.table(DF_Counts_SE, paste("CPMs.",name,".txt", sep=""), row.names = F, quote = F, sep = '\t')


}

for (i in 1:length(myfiles.paired)){
	
	name <- myfiles.paired[i]
	name <- gsub(".bam", "", name)

	Counts_PE <- featureCounts(files=myfiles.paired[i], isGTFAnnotationFile=FALSE, annot.ext="All.peaks.saf", 
				countMultiMappingReads=FALSE, isPairedEnd=TRUE, 
				useMetaFeatures=FALSE, 	nthreads=40, allowMultiOverlap=TRUE, minOverlap=25)

	DF_Counts_PE <- make_df_tpms(Counts_PE)
	colnames(DF_Counts_PE)[2] <- name
	write.table(DF_Counts_PE, paste("CPMs.",name,".txt", sep=""), row.names = F, quote = F, sep = '\t')

}

############ Get DF  with counts data  ######################## 
#

#DF_Counts_SE <- make_df_tpms(Counts_SE)
#DF_Counts_PE <- make_df_tpms(Counts_PE)
#

#write.table(DF_Counts_SE, "CPMs_SE.txt", row.names = F, quote = F, sep = '\t')
#write.table(DF_Counts_PE, "CPMs_PE.txt", row.names = F, quote = F, sep = '\t')
#

#DF_Counts <- left_join(DF_Counts_SE, DF_Counts_PE, by="Summit")
#write.table(DF_Counts, "CPMs.all.txt", row.names = F, quote = F, sep = '\t')
