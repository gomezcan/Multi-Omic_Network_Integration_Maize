## Tidy per-arm metrics: yield, retention, jaccard, all-pairs rho + stability scores
suppressMessages(library(data.table))
S <- Sys.getenv("SENS_DIR"); OUT <- file.path(S,"results")
Ds <- c(0.005,0.02,0.05,0.10,0.20)
arm <- Sys.getenv("ARM", "PWY")
f <- function(d) file.path(OUT, sprintf("%s_enrichment_D%g.txt", arm, d))
tabs <- lapply(Ds, function(d) fread(f(d)))
names(tabs) <- as.character(Ds)
if (arm=="GO") for (n in names(tabs)) setnames(tabs[[n]], c("GO.ID"), c("PWY"), skip_absent=TRUE)
base <- tabs[["0.005"]]
rows <- list()
for (d in as.character(Ds)) {
  t <- tabs[[d]]
  for (cut in c(0.05,0.10)) {
    bset <- base[FDR<=cut, paste(TF,PWY)]; sset <- t[FDR<=cut, paste(TF,PWY)]
    m <- merge(base[,.(TF,PWY,Pb=Pval%||%classic)], t[,.(TF,PWY,Pd=Pval%||%classic)],
               by=c("TF","PWY")) |> suppressWarnings()
    rows[[length(rows)+1]] <- data.table(arm=arm, D=as.numeric(d), FDRcut=cut,
      nTFs_sig=t[FDR<=cut, uniqueN(TF)], nPairs_sig=length(sset),
      retention=mean(bset %in% sset), jaccard=length(intersect(bset,sset))/length(union(bset,sset)))
  }
}
`%||%` <- function(a,b) if (!is.null(a)) a else b
M <- rbindlist(rows)
fwrite(M, file.path(OUT, sprintf("metrics_%s.txt", arm)), sep="\t")
## stability score n-of-5 (FDR<=0.1)
allsig <- rbindlist(lapply(as.character(Ds), function(d) tabs[[d]][FDR<=0.10, .(TF,PWY,D=d)]))
stab <- allsig[, .(n_thresholds_sig=.N), by=.(TF,PWY)]
fwrite(stab[order(-n_thresholds_sig)], file.path(OUT, sprintf("stability_%s.txt", arm)), sep="\t")
cat(sprintf("[%s] metrics + stability written; core >=3/5: %d pairs, %d TFs\n",
    arm, stab[n_thresholds_sig>=3,.N], stab[n_thresholds_sig>=3, uniqueN(TF)]))
print(M[FDRcut==0.10])
