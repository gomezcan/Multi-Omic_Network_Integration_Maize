suppressMessages(library(data.table))
S <- Sys.getenv("SENS_DIR")
# syntenic set
syn <- fread(file.path(S,"Zm.v4.synteny.genes.txt"), sep="\t")
syn <- syn[gene_synteny=="syntenic", gene_id]
tfs <- fread(file.path(S,"All_TFs.txt"), header=FALSE)$V1
# edge file
E <- fread(file.path(S,"InputClusterONE_Dim50_WL80_nW10_0.005_syntenic.txt"),
           col.names=c("V1","V2","w"))
E <- E[V1 %chin% syn & V2 %chin% syn]
cat("edges after syntenic filter:", nrow(E), "\n")
Etf <- unique(E[V1 %chin% tfs, .(V1,V2)])
sizes <- Etf[, .N, by=V1]
cat("TFs with >=1 edge:", nrow(sizes), "\n")
# baseline nTF
B <- fread(file.path(S,"NetworkBased_PWY_Clusters_enrichment.txt"),
           select=c("TF","nTF"))
B <- unique(B)
cat("baseline TFs:", nrow(B), "\n")
M <- merge(sizes, B, by.x="V1", by.y="TF", all=TRUE)
cat("TFs matching size exactly:", M[N==nTF, .N], "of", nrow(M), "\n")
print(head(M[N!=nTF | is.na(N) | is.na(nTF)], 10))
# genome size used by newGOM
gs <- length(intersect(syn, unique(c(Etf$V1, Etf$V2))))
cat("genome.size (syntenic in TF-network):", gs, "\n")
