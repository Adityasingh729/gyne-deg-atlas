#!/usr/bin/env Rscript
# Pull analysis-ready expression for ONE dataset and save a tidy RDS:
#   list(counts=matrix genes x samples, meta=data.frame samples x covariates)
# RNA-seq  -> recount3 (preferred) or ARCHS4
# Microarray -> GEOquery getGEO() processed matrix
# 8 GB: load ONE dataset at a time; gc() after.
suppressMessages({library(optparse)})
opt <- parse_args(OptionParser(option_list=list(
  make_option("--dataset"), make_option("--technique"),
  make_option("--source"), make_option("--out"))))

# TODO: implement per source.
# if (grepl("RNA-seq", opt$technique)) {
#   library(recount3)
#   rse <- ... create_rse for opt$dataset ...
#   counts <- compute_read_counts(rse); meta <- as.data.frame(colData(rse))
# } else {
#   library(GEOquery)
#   g <- getGEO(opt$dataset, GSEMatrix=TRUE)[[1]]
#   counts <- exprs(g); meta <- pData(g)
# }
# saveRDS(list(counts=counts, meta=meta), opt$out)
stop("fetch_expr.R is a stub - implement recount3/ARCHS4/GEO retrieval for your datasets.")
