#!/usr/bin/env Rscript
# Combine per-dataset DEGs within a stratum.
# RRA: rank up- and down- lists separately, RobustRankAggreg::aggregateRanks.
# REM: random-effects on log2FC via metafor::rma (needs SE; derive from p + n).
suppressMessages({library(optparse); library(data.table)})
opt <- parse_args(OptionParser(option_list=list(
  make_option("--stratum"), make_option("--method",default="RRA"),
  make_option("--min_datasets",type="integer",default=4),
  make_option("--deg_dir"), make_option("--samplesheet"), make_option("--out"))))
# TODO:
#  1) read samplesheet, select dataset_ids in this stratum
#  2) if (n < min_datasets) write a flagged single/low-N descriptive result and STOP meta
#  3) RRA: library(RobustRankAggreg); up <- aggregateRanks(list_of_gene_vectors_by_lfc)
#  4) REM: library(metafor); per-gene rma(yi=log2FC, sei=SE)
#  5) write: gene, meta_log2FC, meta_p, meta_padj, n_datasets, direction
stop("meta_combine.R is a stub - implement RRA/REM as commented.")
