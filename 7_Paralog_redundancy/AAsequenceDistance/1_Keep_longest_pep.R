
chop <- function(myStr,mySep,myField){
  ## chop a string by a separator and return specified field
  choppedString=sapply(strsplit(myStr,mySep),"[",myField)
  
  if(length(myField)>1){
    
    choppedString=apply(choppedString,2,function(x){paste0(x[!is.na(x)],collapse=mySep)})
  }
  
  return(choppedString)
}

###################################################################################
#######                             Libraries                               #######
###################################################################################

library(tidyverse)
library(data.table)

###################################################################################

## Read pep lengths to selected longest version of react protein

Peplength <- fread("All_Paralog.pep.IDs.len.txt", header = F) 
Peplength[,'id'] <- chop(Peplength$V1, '[_]', 1)

Peplength %>% 
  arrange(desc(V2)) %>%
  group_by(id) %>%
  slice(1) -> Peplength_filtered

fwrite(Peplength_filtered[,1], 
       "All_Paralog.pep.IDs.Filterlen.txt", row.names=F, col.names = F, quote = F)

#Peplength[grepl("Zm00001d001879", Peplength$V1),] %>%
#  arrange(desc(V2))
