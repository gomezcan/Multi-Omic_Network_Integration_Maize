suppressMessages(library(data.table))
S <- Sys.getenv("SENS_DIR")
syn <- fread(file.path(S,"Zm.v4.synteny.genes.txt"), sep="\t")
syn <- syn[gene_synteny=="syntenic", gene_id]
tfs <- fread(file.path(S,"All_TFs.txt"), header=FALSE)$V1
# CornCYC exactly as V6
cc <- fread(file.path(S,"corn_pathways.0210325.reduced.txt"), sep="\t")
setnames(cc, c("Pathway.id","Pathway.name","GeneID"))
cc <- cc[GeneID!="unknown" & GeneID %chin% syn]
cc[, Pathway.id := gsub("-","_",Pathway.id)]
ccl <- split(cc$GeneID, cc$Pathway.id)             # with duplicates (nPWY)
ccl_u <- lapply(ccl, unique)                       # unique sets (test)
# edges → neighborhoods
E <- fread(file.path(S,"InputClusterONE_Dim50_WL80_nW10_0.005_syntenic.txt"),
           col.names=c("V1","V2","w"))
E <- E[V1 %chin% syn & V2 %chin% syn]
Etf <- unique(E[V1 %chin% tfs, .(V1,V2)])
N <- 24597
B <- fread(file.path(S,"NetworkBased_PWY_Clusters_enrichment.txt"))
set.seed(1)
sampleTFs <- sample(unique(B[n.targ>0, TF]), 8)
res <- list()
for (tf in sampleTFs) {
  A <- unique(Etf[V1==tf, V2])
  for (pwy in names(ccl_u)) {
    Bset <- ccl_u[[pwy]]
    k <- length(intersect(A, Bset))
    p <- phyper(k-1, length(Bset), N-length(Bset), length(A), lower.tail=FALSE)
    res[[length(res)+1]] <- data.table(TF=tf, PWY=pwy, myP=p, myK=k)
  }
}
R <- rbindlist(res)
M <- merge(R, B[TF %in% sampleTFs, .(TF,PWY,Pval,n.targ)], by=c("TF","PWY"))
M[, dP := abs(myP-Pval)]
cat("pairs compared:", nrow(M), "\n")
cat("k matches:", M[myK==n.targ,.N], "/", nrow(M), "\n")
cat("p matches (<1e-8):", M[dP<1e-8,.N], "/", nrow(M), "\n")
print(M[order(-dP)][1:5])
