
#args = commandArgs(trailingOnly=TRUE)
library(Rsubread)


##### Feacture counts single end bam files ######

myfiles.singles = list.files(pattern = "*.bam") # get list of bam files singles from part of paired-end reads


############ count alingments single end ######################## 
test_single <-featureCounts(files=myfiles.singles, isGTFAnnotationFile=FALSE, 
                            annot.ext="ZmV4.46.gene.saf",
                            countMultiMappingReads=FALSE,
                            allowMultiOverlap=FALSE,
                            isPairedEnd=FALSE, useMetaFeatures=TRUE, nthreads=25)

featureCounts
############ Get DF  with counts data  ######################## 

DF_testsingle <- data.frame(test_single$annotation[,c("GeneID","Strand","Length")], 
                           test_single$counts, stringsAsFactors=FALSE)

colnames(DF_testsingle) <- gsub("protoplast_","", colnames(DF_testsingle))
colnames(DF_testsingle) <- gsub(".bam", "", colnames(DF_testsingle))

write.table(DF_testsingle, "Counts_SE.txt", row.names = F, quote = F, sep = '\t')
