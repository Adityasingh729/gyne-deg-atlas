#!/usr/bin/env Rscript
# clusterProfiler enrichment + STRING PPI hub detection on the meta signature.
suppressMessages({library(optparse); library(data.table)})
opt <- parse_args(OptionParser(option_list=list(
  make_option("--meta"), make_option("--out_enrich"), make_option("--out_hubs"))))
# TODO:
#  - sig <- meta genes with meta_padj < 0.05 & |meta_log2FC| > 1
#  - clusterProfiler::enrichGO / enrichKEGG (org.Hs.eg.db)
#  - STRINGdb -> igraph -> degree/betweenness -> top hub genes
stop("enrich_hubs.R is a stub - implement clusterProfiler + STRINGdb.")
