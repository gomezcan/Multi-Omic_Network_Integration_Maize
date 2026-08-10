## Process GO sweep: yield/retention/jaccard at FDR<=0.1 (paper rule) +
## S8 parent-level stability via the recovered baseline child->parent map.
suppressMessages(library(data.table))
S <- Sys.getenv("SENS_DIR"); OUT <- file.path(S,"results")
Ds <- c("0.005","0.02","0.05","0.1","0.2")
tabs <- lapply(Ds, function(d) fread(file.path(OUT, sprintf("GO_enrichment_D%s.txt", d))))
names(tabs) <- Ds
for (d in Ds) tabs[[d]][, sig := (!is.na(FDR) & FDR <= 0.1) | is.na(classic)]
cat("NA-classic rows (ultra-significant, counted as sig):",
    sum(sapply(tabs, function(t) t[is.na(classic), .N])), "\n")
base <- tabs[["0.005"]]
cat(sprintf("baseline: %d sig associations, %d TFs (original run: 28,664 / 2,915)\n",
    base[sig==TRUE, .N], base[sig==TRUE, uniqueN(TF)]))
rows <- list()
bset <- base[sig==TRUE, paste(TF, GO.ID)]
for (d in Ds) {
  s <- tabs[[d]][sig==TRUE, paste(TF, GO.ID)]
  rows[[d]] <- data.table(D=d, nTFs=tabs[[d]][sig==TRUE, uniqueN(TF)], nAssoc=length(s),
    retention=round(mean(bset %in% s),3),
    jaccard=round(length(intersect(bset,s))/length(union(bset,s)),3))
}
M <- rbindlist(rows); print(M)
fwrite(M, file.path(OUT,"metrics_GO.txt"), sep="\t")

## S8 parent stability
RED <- readRDS(file.path(S,"GO_Network_Red_baseline.rds"))          # go, parent, TF (+terms)
map <- unique(as.data.table(RED)[, .(TF, parent, go)])
## singles: S8 GO rows not covered by the reduction (parent = the term itself)
TN <- fread(file.path(S,"Total_NetworkBased_predictions.txt"))
s8go <- unique(TN[Annotation=="GO", .(TF, parent=Ann.ID)])
covered <- unique(map[, .(TF, parent)])
singles <- fsetdiff(s8go, covered)
cat("S8 GO rows:", nrow(s8go), "| via reduction:", nrow(covered),
    "| singles (self-mapped):", nrow(singles), "\n")
map <- rbind(map, singles[, .(TF, parent, go=parent)])
## support per threshold
for (d in Ds) {
  sd <- tabs[[d]][sig==TRUE, .(TF, go=GO.ID)]
  sd[, hit := TRUE]
  map <- merge(map, sd, by=c("TF","go"), all.x=TRUE)
  map[is.na(hit), hit := FALSE]
  setnames(map, "hit", paste0("t", d))
}
stab <- map[, lapply(.SD, any), by=.(TF, parent), .SDcols=paste0("t", Ds)]
stab[, n_thresholds_supported := rowSums(.SD), .SDcols=paste0("t", Ds)]
res <- merge(s8go, stab[, .(TF, parent, n_thresholds_supported)], by=c("TF","parent"), all.x=TRUE)
res[is.na(n_thresholds_supported), n_thresholds_supported := 0]
cat("== S8 GO rows by stability (n of 5 thresholds supported) ==\n")
print(table(res$n_thresholds_supported))
cat(sprintf(">=3/5: %d of %d (%.1f%%)\n", res[n_thresholds_supported>=3,.N], nrow(res),
    100*res[n_thresholds_supported>=3,.N]/nrow(res)))
fwrite(res[order(-n_thresholds_supported)], file.path(OUT,"stability_GO_S8rule.txt"), sep="\t")
cat("saved metrics_GO.txt + stability_GO_S8rule.txt\n")
