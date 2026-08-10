## GO (topGO BP, classic Fisher) threshold-sensitivity sweep — mirrors
## Fig_pecanpyV6.R GetGO/SuperGO_Modules_Targ verbatim, incl. GenTable
## topNodes=1000 and per-TF FDR over that table.
S <- Sys.getenv("SENS_DIR"); lib <- file.path(S,"Rlib")
.libPaths(c(lib, .libPaths()))
suppressMessages({library(data.table); library(topGO); library(parallel)})
OUT <- file.path(S,"results"); dir.create(OUT, showWarnings=FALSE)
VALIDATE <- Sys.getenv("VALIDATE")=="1"

syn <- fread(file.path(S,"Zm.v4.synteny.genes.txt"), sep="\t")
syn <- syn[gene_synteny=="syntenic", gene_id]
tfs <- fread(file.path(S,"All_TFs.txt"), header=FALSE)$V1
background <- readMappings(file.path(S,"synteny.ID_TopGO_V4_GRAMER.txt"))

E <- fread(file.path(S,"InputClusterONE_Dim50_WL80_nW10_0.005_syntenic.txt"),
           col.names=c("V1","V2","w"))
E <- E[V1 %chin% syn & V2 %chin% syn]

Ds <- if (VALIDATE) 0.005 else c(0.005,0.02,0.05,0.10,0.20)
wmins <- (Ds/0.005==1)*0 + (Ds/0.005!=1)*Ds^2   # 0 for baseline else D^2

for (i in seq_along(Ds)) {
  D <- Ds[i]; wmin <- wmins[i]
  Ei <- if (wmin==0) E else E[w >= wmin]
  Etf <- unique(Ei[V1 %chin% tfs, .(V1,V2)])
  netGenes <- unique(c(Etf$V1, Etf$V2))
  background_tem <- background[names(background) %in% netGenes]
  uni <- as.character(unique(names(background_tem)))
  nbrs <- split(Etf$V2, Etf$V1)
  tflist <- if (VALIDATE) c("Zm00001d001824","Zm00001d006236","Zm00001d005016") else names(nbrs)
  cat(sprintf("D=%g: %d TFs, universe %d genes\n", D, length(nbrs), length(uni)))

  ## build topGOdata once (verbatim GetGO calls; dummy gene list)
  gl0 <- factor(as.integer(uni %in% nbrs[[tflist[1]]])); names(gl0) <- uni
  GOdata <- new("topGOdata", ontology="BP", allGenes=gl0,
                annot=annFUN.gene2GO, gene2GO=background_tem)
  test.stat <- new("classicCount", testStatistic=GOFisherTest,
                   name="Fisher test", nodeSize=10)

  oneTF <- function(tf) {
    degs <- unique(nbrs[[tf]])
    gl <- factor(as.integer(uni %in% degs)); names(gl) <- uni
    gd <- updateGenes(GOdata, gl)
    res <- getSigGroups(gd, test.stat)
    tab <- as.data.table(GenTable(gd, classic=res, topNodes=1000, orderBy='Fis'))
    tab[, classic := suppressWarnings(as.numeric(classic))]
    tab[, `:=`(TF=tf, nTF=length(degs))]
    tab[, FDR := p.adjust(classic, method="fdr")]
    tab[classic <= 0.05 | FDR <= 0.2]
  }
  cores <- if (VALIDATE) 3 else 10
  resL <- mclapply(tflist, function(tf)
    tryCatch(oneTF(tf), error=function(e) NULL), mc.cores=cores)
  ok <- !sapply(resL, is.null)
  cat(sprintf("  completed %d/%d TFs\n", sum(ok), length(tflist)))
  RES <- rbindlist(resL[ok])
  fwrite(RES, file.path(OUT, sprintf("GO_enrichment_D%g%s.txt", D,
         ifelse(VALIDATE,"_VALIDATE",""))), sep="\t")
}
cat("ALL DONE\n")
