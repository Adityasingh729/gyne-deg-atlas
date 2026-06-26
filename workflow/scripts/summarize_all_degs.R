#!/usr/bin/env Rscript
# Script to compile a master DEG summary across all datasets in config/samplesheet.csv

suppressMessages({
  library(data.table)
})

samplesheet_path <- "config/samplesheet.csv"
if (!file.exists(samplesheet_path)) {
  stop("Samplesheet not found at: ", samplesheet_path)
}

samples <- fread(samplesheet_path)
results <- list()

for (i in 1:nrow(samples)) {
  ds <- samples$dataset_id[i]
  cond <- samples$condition[i]
  tissue <- samples$tissue_primary[i]
  tech <- samples$technique[i]
  group_col <- samples$group_column[i]
  case_lbl <- samples$case_label[i]
  ctrl_lbl <- samples$control_label[i]
  
  expr_file <- file.path("data", "expr", paste0(ds, ".expr.rds"))
  deg_file <- file.path("results", "deg", paste0(ds, ".deg.tsv"))
  
  n_case <- 0
  n_control <- 0
  n_deg <- 0
  status <- "pass"
  notes <- ""
  
  if (file.exists(expr_file)) {
    x <- readRDS(expr_file)
    meta <- x$meta
    
    if (group_col %in% colnames(meta)) {
      n_case <- sum(meta[[group_col]] == case_lbl, na.rm=TRUE)
      n_control <- sum(meta[[group_col]] == ctrl_lbl, na.rm=TRUE)
    } else {
      notes <- paste0("Group column '", group_col, "' not found. ")
    }
  } else {
    status <- "missing expression"
  }
  
  if (file.exists(deg_file)) {
    deg <- fread(deg_file)
    if ("padj" %in% colnames(deg)) {
      n_deg <- sum(deg$padj < 0.05, na.rm=TRUE)
    } else {
      notes <- paste0(notes, "Missing padj column. ")
    }
  } else {
    status <- "missing deg"
  }
  
  results[[ds]] <- data.table(
    dataset_id = ds,
    condition = cond,
    tissue = tissue,
    technique = tech,
    n_case = n_case,
    n_control = n_control,
    n_DEG_padj05 = n_deg,
    status = status,
    notes = trimws(notes)
  )
}

dt_results <- rbindlist(results)
fwrite(dt_results, "results/all_datasets_deg_summary.tsv", sep="\t")
fwrite(dt_results, "results/all_datasets_deg_summary.csv", sep=",")

# Output markdown table
cat("\n# Master DEG Analysis Summary (All Datasets)\n\n")
cat("| Dataset ID | Condition | Tissue | Technique | Case (N) | Control (N) | Sig DEGs (padj < 0.05) | Status | Notes |\n")
cat("|---|---|---|---|---|---|---|---|---|\n")
for (i in 1:nrow(dt_results)) {
  r <- dt_results[i]
  cat(sprintf("| [%s](file:///c:/Users/adity/OneDrive/Documents/gyne-deg-atlas/gyne-deg-atlas/results/deg/%s.deg.tsv) | %s | %s | %s | %d | %d | %d | %s | %s |\n",
              r$dataset_id, r$dataset_id, r$condition, r$tissue, r$technique, r$n_case, r$n_control, r$n_DEG_padj05, r$status, r$notes))
}
cat("\n")
 