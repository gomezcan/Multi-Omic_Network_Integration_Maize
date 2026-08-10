## Sensitivity of network-based TF annotation to the MR-MI similarity threshold
## (reviewer points R1-24 / R2-14).
## Replicates Fig_pecanpyV6.R's neighborhood construction + GeneOverlap Fisher
## (validated: neighborhood sizes and p-values match the published baseline
##  to <1e-12) and re-runs CornCyc PWY enrichment at stricter decay cutoffs
## D = exp(-(MRMI-1)/50); edge file stores w = D^2 (4 dp).
suppressMessages(library(data.table))

S <- Sys.getenv("SENS_DIR")
OUT <- file.path(S, "results")

syn <- fread(file.path(S,"Zm.v4.synteny.genes.txt"), sep="\t")
syn <- syn[gene_synteny=="syntenic", gene_id]
tfs <- fread(file.path(S,"All_TFs.txt"), header=FALSE)$V1

cc <- fread(file.path(S,"corn_pathways.0210325.reduced.txt"), sep="\t")
setnames(cc, c("Pathway.id","Pathway.name","GeneID"))
cc <- cc[GeneID!="unknown" & GeneID %chin% syn]
cc[, Pathway.id := gsub("-","_",Pathway.id)]
ccu <- unique(cc[, .(Pathway.id, GeneID)])          # unique gene sets
pwyN <- ccu[, .(nB=.N), by=Pathway.id]

E <- fread(file.path(S,"InputClusterONE_Dim50_WL80_nW10_0.005_syntenic.txt"),
           col.names=c("V1","V2","w"))
E <- E[V1 %chin% syn & V2 %chin% syn]

grid <- data.table(D=c(0.005,0.02,0.05,0.10,0.20),
                   wmin=c(0, 0.02^2, 0.05^2, 0.10^2, 0.20^2))
grid[, MRmax := round(1 - 50*log(D))]

summary_rows <- list()
for (i in seq_len(nrow(grid))) {
  D <- grid$D[i]; wmin <- grid$wmin[i]
  Ei <- if (wmin==0) E else E[w >= wmin]
  Etf <- unique(Ei[V1 %chin% tfs, .(V1,V2)])
  netGenes <- unique(c(Etf$V1, Etf$V2))
  N <- sum(syn %chin% netGenes)                      # genome.size per threshold
  sizes <- Etf[, .(nA=.N), by=V1]
  ## overlap counts per TF x PWY via join on neighbor gene
  K <- merge(Etf, ccu, by.x="V2", by.y="GeneID", allow.cartesian=TRUE)[
        , .(k=.N), by=.(V1, Pathway.id)]
  ## full grid TF x PWY (dense, incl. zero overlaps for faithful FDR)
  full <- CJ(V1=sizes$V1, Pathway.id=pwyN$Pathway.id)
  full <- merge(full, K, by=c("V1","Pathway.id"), all.x=TRUE)
  full[is.na(k), k := 0L]
  full <- merge(full, sizes, by="V1")
  full <- merge(full, pwyN, by="Pathway.id")
  full[, Pval := phyper(k-1, nB, N-nB, nA, lower.tail=FALSE)]
  full[, FDR := p.adjust(Pval, method="fdr"), by=V1]
  setnames(full, c("V1","k"), c("TF","n.targ"))
  fwrite(full[, .(TF, PWY=Pathway.id, Pval, FDR, n.targ, nTF=nA, nPWY=nB)],
         file.path(OUT, sprintf("PWY_enrichment_D%g.txt", D)), sep="\t")
  summary_rows[[i]] <- data.table(
    D=D, MRmax=grid$MRmax[i], edges=nrow(Ei), TF_edges=nrow(Etf),
    genome.size=N, nTFs=nrow(sizes),
    median_nb=as.numeric(median(sizes$nA)), mean_nb=mean(sizes$nA),
    TFs_sig_FDR05=full[FDR<=0.05, uniqueN(TF)],
    TFs_sig_FDR10=full[FDR<=0.10, uniqueN(TF)],
    pairs_sig_FDR05=full[FDR<=0.05, .N],
    pairs_sig_FDR10=full[FDR<=0.10, .N])
  cat(sprintf("done D=%g (MR<=%d): %d TFs, %d sig pairs FDR<=0.05\n",
      D, grid$MRmax[i], nrow(sizes), summary_rows[[i]]$pairs_sig_FDR05))
}
SUM <- rbindlist(summary_rows)
fwrite(SUM, file.path(OUT,"sweep_summary.txt"), sep="\t")
print(SUM)
