#!/usr/bin/env Rscript
# Per-dataset DEG via DESeq2. Writes the standardized schema (see docs/DEG_output_schema.md).
suppressMessages({library(optparse); library(DESeq2); library(data.table)})
opt <- parse_args(OptionParser(option_list=list(
  make_option("--expr"), make_option("--group"), make_option("--case"),
  make_option("--control"), make_option("--padj",type="double",default=0.05),
  make_option("--lfc",type="double",default=1.0), make_option("--out"))))
x <- readRDS(opt$expr); counts <- round(as.matrix(x$counts)); meta <- x$meta
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
counts <- counts[, keep_samples, drop=FALSE]

meta[[opt$group]] <- factor(meta[[opt$group]], levels=c(opt$control, opt$case))
keep <- rowSums(counts >= 10) >= max(2, floor(0.2*ncol(counts)))   # low-count filter
dds <- DESeqDataSetFromMatrix(counts[keep,], meta, as.formula(paste0("~ ", opt$group)))
dds <- DESeq(dds); res <- as.data.frame(results(dds))
res$gene <- rownames(res)
out <- data.table(gene=res$gene, log2FC=res$log2FoldChange, p=res$pvalue,
                  padj=res$padj, n=ncol(counts))
fwrite(out, opt$out, sep="\t")
gc()
 