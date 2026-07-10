suppressMessages(library(data.table))
suppressMessages(library(scales))
suppressMessages(library(purrr))
suppressMessages(library(stringr))


##################################################
######        Setup pecanpy entry          #######
##################################################

#########
args = commandArgs(trailingOnly=TRUE)

Netid <- args[1] # Random net id
#Netid <- 1

# GRN
GRN <- fread(paste0('RandomNets/Random.GRN.',Netid, '.txt'))

# eGRN
eGRN <- fread(paste0('RandomNets/Random.eGRN.',Netid, '.txt'))

# CEN
CEN <- fread(paste0('RandomNets/Random.CEN.',Netid, '.txt'))

# GAN
GAN <- fread(paste0('RandomNets/Random.GAN.',Netid, '.txt'))

## Syntenic genes 
Syntenic <- as_tibble(read.table("Zm.v4.synteny.genes.txt", h=T, stringsAsFactors = F, sep = '\t'))
Syntenic <- subset(Syntenic, gene_synteny=='syntenic')$gene_id

# Reduce nets to targets in syntenic
GRN <- GRN[GRN$Target %in% Syntenic,]
eGRN <- eGRN[eGRN$Target %in% Syntenic,]
CEN <- CEN[CEN$Target %in% Syntenic,]
GAN <- GAN[GAN$Target %in% Syntenic,]

####
GRN <-  paste0(GRN$Source, "_",  GRN$Target)
eGRN <- paste0(eGRN$Source, "_", eGRN$Target)
CEN <-  paste0(CEN$Source, "_", CEN$Target)
GAN <-  paste0(GAN$Source, "_", GAN$Target)

uniNetwork <- data.table(table(c(GRN, eGRN, CEN, GAN)))
weight <- 0.5 + (0.5/4)*uniNetwork$N

uniNetwork <-  data.table(str_split_fixed(uniNetwork$V1, pattern = "_", n = 2))
uniNetwork[,"w"] <- weight

fwrite(uniNetwork, paste0('Random_wNets/Random.',Netid, '.FullwNets.txt'), sep = '\t', quote = F, col.names = F)





