## Deterministic tercile tiers for Fig 4A (replaces the original k-means clusters).
## Input:  GOsByTF — long table with one row per (TF, dataset) and a Ratio column
##         (% of the TF's GO terms significantly enriched in that dataset), as
##         produced in Fig_GSEAv2.R ("GOsByTF" object).
## Output: per-TF tier assignments (Tier1_high / Tier2_mid / Tier3_low) by
##         terciles of the mean Ratio over tested datasets.
## Published cutoffs (829 TFs x 39 datasets): 45.5% and 60.9%;
## tiers of 277 / 276 / 276 TFs with means 72.5 / 53.3 / 34.0%.
suppressMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)   # 1: GOsByTF table (tsv), 2: output path
GB  <- fread(args[1])
out <- if (length(args) > 1) args[2] else "fig4A_tercile_assignments.txt"

pm <- GB[, .(meanRatio = mean(Ratio)), by = TF]
q  <- quantile(pm$meanRatio, c(1/3, 2/3))
pm[, Tier := fifelse(meanRatio >= q[2], "Tier1_high",
             fifelse(meanRatio >= q[1], "Tier2_mid", "Tier3_low"))]

cat("tercile cutoffs (mean % enriched):", round(q, 1), "\n")
print(pm[, .(nTFs = .N, mean_pct = round(mean(meanRatio), 1)), by = Tier][order(Tier)])
fwrite(pm[order(-meanRatio)], out, sep = "\t")
cat("saved", out, "\n")
