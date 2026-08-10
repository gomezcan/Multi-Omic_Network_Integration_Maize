## R2-17: random/null control for the paralog shared-association Jaccard (Fig 5G)
## Observed: J_MR per paralog pair (932), binned I..IX. Null: random TF pairs,
## global + neighborhood-size-matched. Empirical p per bin.
suppressMessages(library(data.table))
set.seed(42)
S <- Sys.getenv("SENS_DIR"); OUT <- file.path(S,"results")
L <- readRDS(file.path(S,"paralog_originals.rds"))
DJ <- L$DFjaccard; MS <- L$MRdb_SCCdb
cat("Index classes: MS=", class(MS$Index), " DJ=", class(DJ$Index), "\n")
MS[, Key := paste0(GeneID1, "_", GeneID2)]
DJ[, Key := if (is.character(Index)) Index else MS$Key[match(seq_len(.N), seq_len(.N))]]
if (!is.character(DJ$Index)) DJ[, Key := MS[match(DJ$ZonesIdex, ZonesIdex), Key]]
tfs <- fread(file.path(S,"All_TFs.txt"), header=FALSE)$V1

E <- fread(file.path(S,"InputClusterONE_Dim50_WL80_nW10_0.005_syntenic.txt"),
           col.names=c("V1","V2","w"))
## neighbor sets: union of out (V1->V2) and in (V2->V1)
outN <- E[, .(nb=list(unique(V2))), by=V1]; setnames(outN, "V1", "g")
inN  <- E[, .(nb=list(unique(V1))), by=V2]; setnames(inN,  "V2", "g")
NB <- merge(outN, inN, by="g", all=TRUE)
NB[, nbAll := mapply(function(a,b) unique(c(a,b)), nb.x, nb.y, SIMPLIFY=FALSE)]
nbs <- setNames(NB$nbAll, NB$g)
sizes <- setNames(lengths(nbs), names(nbs))
jac <- function(a,b){ i <- length(intersect(a,b)); if (i==0) return(0); i/length(union(a,b)) }

## 1) validate replication on the 932 observed pairs
MS[, myJ := mapply(function(g1,g2) jac(nbs[[g1]], nbs[[g2]]), GeneID1, GeneID2)]
V <- if (is.character(DJ$Index)) merge(MS[,.(Key,myJ)], DJ[,.(Key=Index,J_MR)], by="Key") else data.table(myJ=MS$myJ, J_MR=DJ$J_MR)
cat("pairs:", nrow(V), " max|myJ - J_MR|:", max(abs(V$myJ - V$J_MR)), "\n")

## 2) null universes
tfU <- intersect(tfs, names(nbs))
paral <- unique(c(MS$Key, paste0(MS$GeneID2,"_",MS$GeneID1)))
rp <- function(n){
  g1 <- sample(tfU, n, replace=TRUE); g2 <- sample(tfU, n, replace=TRUE)
  k <- g1!=g2 & !(paste0(g1,"_",g2) %in% paral); data.table(g1=g1[k], g2=g2[k]) }
NULLG <- unique(rp(30000))[1:20000]
NULLG[, J := mapply(function(a,b) jac(nbs[[a]], nbs[[b]]), g1, g2)]
cat("global null: n=", nrow(NULLG), " mean J=", round(mean(NULLG$J),4),
    " q95=", round(quantile(NULLG$J,.95),4), "\n")

## 3) size-matched null (quintiles of neighborhood size within TF universe)
qs <- quantile(sizes[tfU], seq(0,1,.2)); qs[1] <- -Inf; qs[length(qs)] <- Inf
qbin <- function(g) findInterval(sizes[[g]], qs, rightmost.closed=TRUE)
tfByQ <- split(tfU, sapply(tfU, qbin))
K <- 20
matched <- MS[, {
  q1 <- as.character(qbin(GeneID1)); q2 <- as.character(qbin(GeneID2))
  a <- sample(tfByQ[[q1]], K, replace=TRUE); b <- sample(tfByQ[[q2]], K, replace=TRUE)
  ok <- a!=b; .(mJ = mapply(function(x,y) jac(nbs[[x]], nbs[[y]]), a[ok], b[ok]))
}, by=.(Key, ZonesIdex)]

## 4) per-bin comparison + empirical p (bootstrap of matched-null bin means)
obsBin <- MS[, .(n=.N, obs_mean=mean(myJ), obs_median=median(myJ)), by=ZonesIdex]
nullBin <- matched[, .(null_mean=mean(mJ)), by=ZonesIdex]
TAB <- merge(obsBin, nullBin, by="ZonesIdex")
TAB[, ratio := round(obs_mean/null_mean,2)]
emp <- matched[, {
  om <- MS[ZonesIdex==.BY$ZonesIdex, mean(myJ)]
  bs <- replicate(2000, mean(sample(mJ, obsBin[ZonesIdex==.BY$ZonesIdex, n], replace=TRUE)))
  .(p_emp = (sum(bs >= om)+1)/2001)
}, by=ZonesIdex]
TAB <- merge(TAB, emp, by="ZonesIdex")
TAB[, gnull_mean := mean(NULLG$J)]
setorder(TAB, ZonesIdex)
print(TAB)
fwrite(TAB, file.path(OUT,"paralog_null_bins.txt"), sep="\t")
fwrite(NULLG, file.path(OUT,"paralog_null_global_pairs.txt"), sep="\t")
fwrite(matched, file.path(OUT,"paralog_null_matched.txt"), sep="\t")
cat("saved paralog_null_{bins,global_pairs,matched}.txt\n")
