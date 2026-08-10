## R2-07: gene-class target/interaction counts normalized by SNP frequency
## Rebuilds the GAN exactly as Fig_transeQTL.v3.R, then adds per-class
## SNP- and gene-count-normalized rates for the S2B-D claims.
suppressMessages(library(data.table))
S <- Sys.getenv("SENS_DIR"); SF <- file.path(S,"snpfreq"); OUT <- file.path(S,"results")

syn <- fread(file.path(S,"Zm.v4.synteny.genes.txt"), sep="\t")
syn <- syn[gene_synteny=="syntenic", gene_id]

## Classes exactly as Fig_transeQTL.v3.R
med <- fread(file.path(SF,"Mediators.txt"), header=FALSE); setnames(med, c("GeneID","Class"))
tfc <- as.data.table(read.table(file.path(SF,"TF_CoR_Mazie.txt"), h=TRUE, stringsAsFactors=FALSE, fill=TRUE))
tfc[, Class := "TF"]
CoRegs <- c("BSD","DDT","FHA","GNAT","HMG","IWS1,SPN1","JUMONJI","LEUNIG",
            "MBF1","MED6","MED7","PHD","SET","SNF2","SOH1","SWI/SNF-BAF60b",
            "SWI/SNF-SWI3","TAZ","TRAF")   # from script list
for (i in CoRegs) tfc[grepl(i, Family, fixed=TRUE), Class := "CoReg"]
tfc <- tfc[!(GeneID %in% med$GeneID)]
cc <- fread(file.path(S,"corn_pathways.0210325.reduced.txt"), sep="\t")
setnames(cc, c("Pathway.id","Pathway.name","GeneID"))
cc <- cc[GeneID!="unknown"][, Class := "Enzyme"]
kin <- fread(file.path(SF,"kinases_maize_AGPv4.txt"))
kin <- kin[GeneID %chin% syn][, Class := "kinases"]
Classes <- unique(rbind(tfc[,.(GeneID,Class)], med[,.(GeneID,Class)],
                        kin[,.(GeneID,Class)], unique(cc[,.(GeneID,Class)])))
Classes <- Classes[GeneID %chin% syn]
cat("class overlap (genes with >1 class):", Classes[, .N, by=GeneID][N>1, .N], "\n")
## priority dedup for clean accounting (TF > CoReg > Mediator > kinases > Enzyme)
pri <- c(TF=1, CoReg=2, Mediator=3, kinases=4, Enzyme=5)
ClassesU <- Classes[order(pri[Class])][, .SD[1], by=GeneID]
classSize <- ClassesU[, .(n_genes=.N), by=Class]

## GAN rebuild
t1 <- fread(file.path(SF,"Clean_trans.eQTL.v2.txt"))
t2 <- fread(file.path(SF,"Clean_trans.eQTLp.v2.txt"))
TE <- unique(rbind(t1[,.(Target,source,Index)], t2[,.(Target,source,Index)]))
cat("combined trans-eQTL rows (Target,source,Index):", nrow(TE), "\n")
TE <- TE[Target %chin% syn & source %chin% syn]
NET <- unique(TE[,.(source,Target)])
cat("GAN edges:", nrow(NET), " sources:", uniqueN(NET$source),
    " targets:", uniqueN(NET$Target), "\n")

addClass <- function(dt, col) {
  dt <- merge(dt, ClassesU, by.x=col, by.y="GeneID", all.x=TRUE)
  dt[is.na(Class), Class := "Other"]; setnames(dt, "Class", paste0(col,"Class")); dt }
NET <- addClass(addClass(NET, "source"), "Target")

## SNPs per class (unique SNP x source-gene pairs; and unique SNPs)
snp_src <- unique(TE[,.(source,Index)])
snp_src <- merge(snp_src, ClassesU, by.x="source", by.y="GeneID", all.x=TRUE)
snp_src[is.na(Class), Class := "Other"]
snpPerClass <- snp_src[, .(n_SNPgene=.N, n_uniqSNP=uniqueN(Index),
                           n_src_genes=uniqueN(source)), by=Class]

## raw claims (S2B-D)
byS <- NET[, .(n_edges=.N, n_targets=uniqueN(Target)), by=sourceClass]
byT <- NET[, .(n_edges=.N, n_target_genes=uniqueN(Target)), by=TargetClass]

## normalized table (source side)
TAB <- merge(merge(byS, snpPerClass, by.x="sourceClass", by.y="Class", all=TRUE),
             classSize, by.x="sourceClass", by.y="Class", all.x=TRUE)
TAB[sourceClass=="Other", n_genes := length(setdiff(syn, ClassesU$GeneID))]
TAB[, `:=`(edges_per_SNP = round(n_edges/n_uniqSNP,3),
           edges_per_gene = round(n_edges/n_genes,2),
           SNPs_per_gene = round(n_uniqSNP/n_genes,2),
           src_rate = round(n_src_genes/n_genes,3))]
setorder(TAB, -n_edges)
cat("\n== SOURCE side: raw vs normalized ==\n"); print(TAB)

TAB2 <- merge(byT, classSize, by.x="TargetClass", by.y="Class", all.x=TRUE)
TAB2[TargetClass=="Other", n_genes := length(setdiff(syn, ClassesU$GeneID))]
TAB2[, `:=`(edges_per_gene = round(n_edges/n_genes,2),
            targets_per_gene = round(n_target_genes/n_genes,3))]
setorder(TAB2, -n_edges)
cat("\n== TARGET side: raw vs normalized ==\n"); print(TAB2)

fwrite(TAB,  file.path(OUT,"snpfreq_source_normalized.txt"), sep="\t")
fwrite(TAB2, file.path(OUT,"snpfreq_target_normalized.txt"), sep="\t")
cat("saved snpfreq_{source,target}_normalized.txt\n")
