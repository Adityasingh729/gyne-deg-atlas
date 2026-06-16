#!/usr/bin/env Rscript
# Per-dataset DEG via limma (microarray). Same output schema as DESeq2.
suppressMessages({library(optparse); library(limma); library(data.table)})
opt <- parse_args(OptionParser(option_list=list(
  make_option("--expr"), make_option("--group"), make_option("--case"),
  make_option("--control"), make_option("--out"))))
x <- readRDS(opt$expr); E <- as.matrix(x$counts); meta <- x$meta   # E = log2 intensities
grp <- factor(meta[[opt$group]], levels=c(opt$control, opt$case))
design <- model.matrix(~grp)
fit <- eBayes(lmFit(E, design))
tt <- topTable(fit, coef=2, number=Inf)
out <- data.table(gene=rownames(tt), log2FC=tt$logFC, p=tt$P.Value,
                  padj=tt$adj.P.Val, n=ncol(E))
fwrite(out, opt$out, sep="\t")
gc()
