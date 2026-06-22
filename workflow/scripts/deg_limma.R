#!/usr/bin/env Rscript
# Per-dataset DEG via limma (microarray). Same output schema as DESeq2.
suppressMessages({library(optparse); library(limma); library(data.table)})
opt <- parse_args(OptionParser(option_list=list(
  make_option("--expr"), make_option("--group"), make_option("--case"),
  make_option("--control"), make_option("--out"))))
x <- readRDS(opt$expr); E <- as.matrix(x$counts); meta <- x$meta   # E = log2 intensities
colnames(meta) <- make.names(colnames(meta))
opt$group <- make.names(opt$group)

# Normalize Mu character to standard 'u' to prevent Windows command line encoding issues
clean_label <- function(s) {
  s <- gsub("\u00B5", "u", s, fixed=TRUE)
  s <- gsub("\u03BC", "u", s, fixed=TRUE)
  return(s)
}
meta[[opt$group]] <- clean_label(as.character(meta[[opt$group]]))
opt$case <- clean_label(opt$case)
opt$control <- clean_label(opt$control)

# Subset samples to keep only case and control groups
keep_samples <- !is.na(meta[[opt$group]]) & (meta[[opt$group]] %in% c(opt$control, opt$case))
if (sum(keep_samples) < 2) {
  stop("Fewer than 2 samples found matching the specified case and control labels.")
}
meta <- meta[keep_samples, , drop=FALSE]
E <- E[, keep_samples, drop=FALSE]

grp <- factor(meta[[opt$group]], levels=c(opt$control, opt$case))
design <- model.matrix(~grp)
fit <- eBayes(lmFit(E, design))
tt <- topTable(fit, coef=2, number=Inf)
out <- data.table(gene=rownames(tt), log2FC=tt$logFC, p=tt$P.Value,
                  padj=tt$adj.P.Val, n=ncol(E))
fwrite(out, opt$out, sep="\t")
gc()
