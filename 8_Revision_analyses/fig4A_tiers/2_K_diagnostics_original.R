suppressMessages({library(data.table); library(cluster); library(ggplot2)})
S <- Sys.getenv("SENS_DIR"); OUT <- file.path(S,"results")
L <- readRDS(file.path(S,"fig4A_originals.rds"))
M <- L$M
scaled <- !is.null(attr(M,"scaled:center"))
cat("original M scaled?", scaled, "\n")
Ms <- if (scaled) M else scale(M)
Mraw <- if (scaled) sweep(sweep(M,2,attr(M,"scaled:scale"),"*"),2,attr(M,"scaled:center"),"+") else M
cl <- as.data.table(L$ClusterDF)
cat("== original cluster sizes ==\n"); print(cl[, .N, by=Cluster])
mr <- data.table(TF=rownames(Mraw), meanRatio=rowMeans(Mraw))
mm <- merge(cl, mr, by="TF")
cat("== original cluster mean %-enriched ==\n")
print(mm[, .(nTFs=.N, mean_pct=round(mean(meanRatio),1)), by=Cluster][order(Cluster)])
## K diagnostics on the original scaled matrix
set.seed(42); ks <- 1:10
wss <- sapply(ks, function(k) kmeans(Ms,k,nstart=25,iter.max=50)$tot.withinss)
sil <- c(NA, sapply(2:10, function(k){km<-kmeans(Ms,k,nstart=25,iter.max=50)
  mean(silhouette(km$cluster, dist(Ms))[,"sil_width"])}))
gap <- clusGap(Ms, FUN=function(x,k) kmeans(x,k,nstart=25,iter.max=50), K.max=10, B=100, verbose=FALSE)
DT <- data.table(k=ks, WSS=wss, avg_silhouette=sil, gap=gap$Tab[,"gap"], gap_SE=gap$Tab[,"SE.sim"])
fwrite(DT, file.path(OUT,"fig4A_K_diagnostics_ORIGINAL.txt"), sep="\t")
print(DT)
cat("gap firstSEmax k:", maxSE(gap$Tab[,"gap"], gap$Tab[,"SE.sim"], "firstSEmax"),
    "| silhouette argmax k:", which.max(sil), "\n")
## does seeded kmeans k=3 reproduce the original partition?
km3 <- kmeans(Ms, 3, nstart=25, iter.max=50)
tab <- table(orig=mm[match(rownames(Ms), TF), Cluster], new=km3$cluster)
print(tab)
cat("agreement (best mapping):", round(sum(apply(tab,1,max))/sum(tab),3), "\n")
## compare my rebuilt matrix to original
G <- rbindlist(lapply(list.files(file.path(S,"gsea/GSEA_results"), full.names=TRUE), fread))
tot <- unique(G[,.(.id,pathway,TF)])[, .(TotalGOs=.N), by=.(Net=.id,TF)]
sig <- unique(G[padj<=0.1,.(.id,pathway,TF)])[, .(SigGOs=.N), by=.(Net=.id,TF)]
R <- merge(tot, sig, by=c("Net","TF"), all.x=TRUE)
R[is.na(SigGOs), SigGOs:=0][, Ratio := round(SigGOs/TotalGOs*100,2)]
M2 <- dcast(R, TF ~ Net, value.var="Ratio")
tf2 <- M2$TF; M2 <- as.matrix(M2[,-1]); rownames(M2) <- tf2; M2[is.na(M2)] <- 0
common_tf <- intersect(rownames(Mraw), rownames(M2))
common_ds <- intersect(colnames(Mraw), colnames(M2))
d <- Mraw[common_tf, common_ds] - M2[common_tf, common_ds]
cat("rebuild vs original:", length(common_tf),"TFs x",length(common_ds),"ds; max|diff|:",
    max(abs(d)), "; frac cells differing >0.01:", round(mean(abs(d)>0.01),4), "\n")
