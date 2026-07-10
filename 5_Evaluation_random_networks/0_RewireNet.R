suppressMessages(library(igraph))
suppressMessages(library(parallel))
suppressMessages(library(data.table))


args = commandArgs(trailingOnly=TRUE)

#####################################################
#############      Functions            #############
#####################################################


RewireNet <- function(net_igraph) {
  ####  
  # This function take as df_net-like to convert it into an igraph 
  # object. rewire it, and return a random network.
  ####
  # Rewire network with similar degree
  igraph_R <- rewire(net_igraph, with = keeping_degseq(loops = FALSE, niter = vcount(net_igraph)*10000))
  # Get  igraph DF
  out <- as.data.table(as_data_frame(igraph_R, what = "edges"))
  
  colnames(out) <- c("Source", "Target")
  return(out)
}

RamdonSample <- function(times, type_net){
  
  # Get random network
  Rnet <- RewireNet(NetTFs_igraph)
  fwrite(Rnet, paste0("RandomNets/Random.",type_net,".",times,".txt"), sep = '\t', row.names = F, quote = F)
  
  
}

#####################################################
#############     Data input            #############
#####################################################

# PDI
filenet = args[1]
# filenet = 'GRN.txt'
filename = args[2]
# filename = 'GRN'

# read net input
Net <- unique(fread(filenet, header = F))
colnames(Net) <- c("Source", "Target")

# from DF to igraph obj.
NetTFs_igraph <- graph_from_data_frame(Net, directed = T) 

set.seed(123)

#####################################################
#############     Data input            #############
#####################################################

RamdonSample(1, filename)

w=100 # Size of range to test
print(".. Ready to start ..")
Samples = 5000

for (i in seq(0, Samples, w)){
  max=Samples
  Start=i+1
  end=i+w
  
  if (end <= max ){
    #
    r.nets <- seq(Start, end, 1)
    #print(r.nets)
    
    mclapply(r.nets, function(x) RamdonSample(x, filename), mc.cores=w)
    
    r <- paste(Start, end, sep = "-")
    print(paste0(" .. Done ", r, " Samples .."))
  }
  
  else if (end > max){ break}
  else{
    r <- paste(Start, max, sep = "-")
    r.nets <- seq(Start, max, 1)
    #print(r.nets)
    w = length(r.nets)+1
    mclapply(r.nets, function(x) RamdonSample(x, filename), mc.cores=w)
    
    print(paste0(" .. Done ", r, " Samples .."))
  }
}


