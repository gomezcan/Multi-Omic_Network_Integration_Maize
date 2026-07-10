library(data.table)
library(scales)
library(purrr)
library(stringr)


##################################################
######        Setup pecanpy entry          #######
##################################################

# PDI
PDI <- unique(fread("../Fig_PDI/Only_PDI_NetworkFinal.10_14_2022.txt")[,c(2,3)])
colnames(PDI)[1] <- "Source"

PDIeQTL <- unique(fread("../Fig_PDI/CisE_PDI_NetworkFinal.10_14_2022.txt")[,c(2,3)])
colnames(PDIeQTL)[1] <- "Source"

# CoExp
CoExp <- unique(fread("../Fig_Coexpression/CoExp_NetworkFinal.10_11_2021.txt"))
colnames(CoExp)[2] <- "Source"
CoExp <- unique(CoExp[,2:3])

# teQTL
teQTL <- unique(fread("../Fig_transeQTL/teQTL_NetworkFinal.10_11_2021.txt"))
colnames(teQTL)[1] <- "Source"

####
PDI_index <- paste0(PDI$Source, "_", PDI$Target)
PDIeQTL_index <- paste0(PDIeQTL$Source, "_", PDIeQTL$Target)

CoExp_index <- paste0(CoExp$Source, "_", CoExp$Target)
teQTL_index <- paste0(teQTL$Source, "_", teQTL$Target)

uniNetwork <- data.table(table(c(PDI_index, PDIeQTL_index, CoExp_index, teQTL_index)))

weight <- 0.5 + (0.5/4)*uniNetwork$N

uniNetwork <-  data.table(str_split_fixed(uniNetwork$V1, pattern = "_", n = 2))
uniNetwork[,"w"] <- weight

fwrite(uniNetwork, "uniqFullNets_weighted.txt", sep = '\t', quote = F, col.names = F)