suppressMessages(library(data.table))
S <- Sys.getenv("SENS_DIR"); OUT <- file.path(S,"results")
Ds <- c(0.005,0.02,0.05,0.10,0.20)
tabs <- lapply(Ds, function(d) fread(file.path(OUT, sprintf("PWY_enrichment_D%g.txt", d))))
names(tabs) <- as.character(Ds)
base <- tabs[["0.005"]]

cat("== all-pairs Spearman of -log10 Pval vs baseline ==\n")
for (d in as.character(Ds)[-1]) {
  m <- merge(base[,.(TF,PWY,Pb=Pval)], tabs[[d]][,.(TF,PWY,Pd=Pval)], by=c("TF","PWY"))
  cat(sprintf("D=%-5s rho=%.3f  (n=%d pairs)\n", d,
      cor(-log10(m$Pb+1e-300), -log10(m$Pd+1e-300), method="spearman"), nrow(m)))
}

cat("\n== retention by baseline strength (pairs FDR<=0.10, tiered by baseline Pval) ==\n")
bs <- base[FDR<=0.10]
bs[, tier := cut(Pval, c(0,1e-5,1e-4,1e-3,1), labels=c("<1e-5","1e-5..1e-4","1e-4..1e-3",">1e-3"))]
for (d in as.character(Ds)[-1]) {
  s <- tabs[[d]][FDR<=0.10, paste(TF,PWY)]
  ret <- bs[, .(retained=mean(paste(TF,PWY) %in% s), n=.N), by=tier][order(tier)]
  cat(sprintf("D=%-5s  %s\n", d,
      paste(sprintf("%s: %.2f (n=%d)", ret$tier, ret$retained, ret$n), collapse="  ")))
}

cat("\n== why baseline-significant pairs drop out (at each D): k->0 vs still-overlapping ==\n")
for (d in as.character(Ds)[-1]) {
  m <- merge(bs[,.(TF,PWY)], tabs[[d]][,.(TF,PWY,FDR,n.targ)], by=c("TF","PWY"))
  lost <- m[FDR>0.10]
  cat(sprintf("D=%-5s lost=%3d  of which k=0: %3d (%.0f%%)\n", d, nrow(lost),
      lost[n.targ==0,.N], 100*lost[n.targ==0,.N]/max(1,nrow(lost))))
}

cat("\n== robust core: pairs significant (FDR<=0.10) in >=3 of 5 thresholds ==\n")
allsig <- rbindlist(lapply(as.character(Ds), function(d)
  tabs[[d]][FDR<=0.10, .(TF,PWY,D=d)]))
core <- allsig[, .N, by=.(TF,PWY)]
cat(sprintf("pairs sig >=1x: %d; >=3x: %d; 5x: %d; TFs in >=3x core: %d\n",
    nrow(core), core[N>=3,.N], core[N==5,.N], core[N>=3, uniqueN(TF)]))

cat("\n== exemplar TFs: best PWY raw Pval at baseline ==\n")
for (tf in c("Zm00001d006236","Zm00001d005016","Zm00001d033859","Zm00001d051520")) {
  r <- base[TF==tf][order(Pval)][1:2]
  cat(sprintf("%s: %s (P=%.3g), %s (P=%.3g)\n", tf, r$PWY[1], r$Pval[1], r$PWY[2], r$Pval[2]))
}
