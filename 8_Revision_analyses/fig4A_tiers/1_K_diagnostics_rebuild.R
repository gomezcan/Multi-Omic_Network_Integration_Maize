## K-selection diagnostics for the Fig 4A TF clustering (R2-14) +
## reproduction of the 3 clusters (sizes/means; feeds R2-13 check).
suppressMessages({library(data.table); library(cluster); library(ggplot2)})
S <- Sys.getenv("SENS_DIR"); OUT <- file.path(S,"results")
fs <- list.files(file.path(S,"gsea/GSEA_results"), full.names=TRUE)
G <- rbindlist(lapply(fs, fread))
cat("GSEA rows:", nrow(G), " TFs:", uniqueN(G$TF), " datasets:", uniqueN(G$.id), "\n")
tot <- unique(G[,.(.id,pathway,TF)])[, .(TotalGOs=.N), by=.(Net=.id,TF)]
sig <- unique(G[padj<=0.1, .(.id,pathway,TF)])[, .(SigGOs=.N), by=.(Net=.id,TF)]
R <- merge(tot, sig, by=c("Net","TF"), all.x=TRUE)
R[is.na(SigGOs), SigGOs:=0][, Ratio := round(SigGOs/TotalGOs*100, 2)]
M <- dcast(R, TF ~ Net, value.var="Ratio")
tfnames <- M$TF; M <- as.matrix(M[,-1]); rownames(M) <- tfnames
M[is.na(M)] <- 0
Ms <- scale(M)
cat("matrix:", nrow(Ms), "TFs x", ncol(Ms), "datasets\n")
set.seed(42)
ks <- 1:10
wss <- sapply(ks, function(k) kmeans(Ms, k, nstart=25, iter.max=50)$tot.withinss)
sil <- c(NA, sapply(2:10, function(k){
  km <- kmeans(Ms, k, nstart=25, iter.max=50)
  mean(silhouette(km$cluster, dist(Ms))[, "sil_width"])}))
gap <- clusGap(Ms, FUN=function(x,k) kmeans(x,k,nstart=25,iter.max=50), K.max=10, B=100, verbose=FALSE)
DT <- data.table(k=ks, WSS=wss, avg_silhouette=sil,
                 gap=gap$Tab[,"gap"], gap_SE=gap$Tab[,"SE.sim"])
fwrite(DT, file.path(OUT,"fig4A_K_diagnostics.txt"), sep="\t")
print(DT)
cat("gap-stat firstSEmax k:", maxSE(gap$Tab[,"gap"], gap$Tab[,"SE.sim"], "firstSEmax"), "\n")
cat("silhouette best k:", which.max(sil), "\n")
## reproduce k=3 clusters
km3 <- kmeans(Ms, 3, nstart=25, iter.max=50)
cl <- data.table(TF=rownames(Ms), cluster=km3$cluster, meanRatio=rowMeans(M))
smry <- cl[, .(nTFs=.N, mean_pct=round(mean(meanRatio),1)), by=cluster][order(-mean_pct)]
print(smry)
fwrite(cl, file.path(OUT,"fig4A_cluster_assignments.txt"), sep="\t")
## figure: 3 diagnostics panels
DL <- melt(DT, id.vars="k", measure.vars=c("WSS","avg_silhouette","gap"))
lab <- c(WSS="Elbow (total within-SS)", avg_silhouette="Average silhouette width", gap="Gap statistic")
DL[, panel := lab[as.character(variable)]]
p <- ggplot(DL, aes(k, value)) + geom_line(linewidth=.5) + geom_point(size=1.8) +
  geom_vline(xintercept=3, linetype=2, colour="grey40") +
  facet_wrap(~panel, scales="free_y") +
  scale_x_continuous(breaks=1:10) + theme_bw(11) +
  labs(x="Number of clusters (k)", y=NULL,
       title="K-selection diagnostics for TF clustering by condition-enrichment profile (Fig 4A)")
ggsave(file.path(OUT,"SFig_fig4A_K_diagnostics.pdf"), p, width=9, height=3.2)
cat("saved SFig_fig4A_K_diagnostics.pdf\n")
