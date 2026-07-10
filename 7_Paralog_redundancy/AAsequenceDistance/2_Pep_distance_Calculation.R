suppressMessages(library(DECIPHER))
suppressMessages(library(data.table))

chop=function(myStr,mySep,myField){
  ## chop a string by a separator and return specified field
  choppedString=sapply(strsplit(myStr,mySep),"[",myField)
  if(length(myField)>1){
    choppedString=apply(choppedString,2,function(x){paste0(x[!is.na(x)],collapse=mySep)})
  }
  return(choppedString)
}

flattenMatrix <- function(cormat) {
  ut <- upper.tri(cormat)
  data.frame(
    GeneID1 = rownames(cormat)[row(cormat)[ut]],
    GeneID2 = rownames(cormat)[col(cormat)[ut]],
    Dis  =(cormat)[ut]
  )
}

#####

args = commandArgs(trailingOnly=TRUE)

# Import name
fastafile <- args[1] 
# fastafile <- "Zm00001d001824_Zm00001d026628.fa"

fas_pp <- list.files(path = '.', pattern = fastafile)

# Read pep sequences
fas_pp <- readAAStringSet(fas_pp)

# reduce pep names
fas_pp@ranges@NAMES <- chop(fas_pp@ranges@NAMES, '[ ]', 1)

# Noters:
## The uncorrected(correction="none") distance 
## matrix represents the hamming distance between
## each of these quences in myXStringSet.
## includeTerminalGaps=FALSE: terminal gaps("-"or"."characters) are not included in sequence length


disM <- DistanceMatrix(fas_pp, type="dist",
                       includeTerminalGaps=TRUE,
                       penalizeGapLetterMatches=TRUE,
                       correction="none", processors=1, verbose=F)


# Get flatten Matrix
disM <- data.table::as.data.table(flattenMatrix(as.matrix(disM)))


fwrite(disM, paste0('Dis.', gsub('.fa', '.txt', fastafile)),
       row.names = F, col.names = T, quote = F)



