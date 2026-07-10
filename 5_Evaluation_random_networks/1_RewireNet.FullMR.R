suppressMessages(library(igraph))
suppressMessages(library(parallel))
suppressMessages(library(data.table))


args = commandArgs(trailingOnly=TRUE)

#####################################################
#############      Functions            #############
#####################################################

RewireNet2 <- function(net_igraph) {
  ####  
  # This function take as df_net-like to convert it into an igraph 
  # object. rewire it, and return a random network.
  ####
  # Rewire network with similar degree
  igraph_R <- rewire(net_igraph, with = keeping_degseq(loops = FALSE,
                                           niter = vcount(net_igraph)*10000))
  # Get  igraph DF
  out <- as.data.table(as_data_frame(igraph_R, what = "edges"))
  
  colnames(out) <- c("Source", "Target")
  return(out)
}

RamdonSample2 <- function(times){
  # Get random network
  Rnet <- RewireNet2(NetTFs_igraph)
  fwrite(Rnet, paste0("RandomNets/Random.FullMR.",times,".txt"), sep = '\t', row.names = F, quote = F)
}

#####################################################
#############     Data input            #############
#####################################################
set.seed(123)

# wNEt from MR after pecanpy analysis
filenet = 'MRnet_syntenic.txt'

# filenet = 'GRN.txt'
filename = args[1]
# filename = '1'

# read net input
Net <- unique(fread(filenet, header = F))[,1:2]
colnames(Net) <- c("Source", "Target")

# from DF to igraph obj.
NetTFs_igraph <- graph_from_data_frame(Net, directed = T) 
# is.weighted(NetTFs_igraph)

#####################################################
#############     Data input            #############
#####################################################

RamdonSample2(filename)


