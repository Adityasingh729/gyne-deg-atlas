#!/usr/bin/env Rscript
# QC script to validate schemas, run PCA separation t-tests, and output the tracking sheet for PCOS Granulosa datasets.

suppressMessages({
  library(data.table)
  library(limma)
})

samplesheet_path <- "config/samplesheet.csv"
samples <- fread(samplesheet_path)

# Filter for PCOS Granulosa stratum
stratum_ds <- samples[condition == "PCOS" & grepl("Granulosa", tissue_primary, ignore.case=TRUE)]

results <- list()

for (i in 1:nrow(stratum_ds)) {
  ds <- stratum_ds$dataset_id[i]
  cond <- stratum_ds$condition[i]
  tissue <- stratum_ds$tissue_primary[i]
  tech <- stratum_ds$technique[i]
  group_col <- stratum_ds$group_column[i]
  case_lbl <- stratum_ds$case_label[i]
  ctrl_lbl <- stratum_ds$control_label[i]
  
  expr_file <- file.path("data", "expr", paste0(ds, ".expr.rds"))
  deg_file <- file.path("results", "deg", paste0(ds, ".deg.tsv"))
  
  n_case <- 0
  n_control <- 0
  n_deg <- 0
  pca_sep <- "No"
  schema_ok <- "No"
  status <- "pass"
  notes <- ""
  
  if (!file.exists(expr_file)) {
    status <- "drop"
    notes <- "Missing expression file"
  } else {
    x <- readRDS(expr_file)
    meta <- x$meta
    
    # Calculate group sizes
    if (group_col %in% colnames(meta)) {
      # Make sure to handle NA values and matching labels
      n_case <- sum(meta[[group_col]] == case_lbl, na.rm=TRUE)
      n_control <- sum(meta[[group_col]] == ctrl_lbl, na.rm=TRUE)
    } else {
      notes <- paste0("Group column '", group_col, "' not found in metadata. ")
    }
    
    # PCA Calculation
    E <- x$counts
    v <- apply(E, 1, var)
    E <- E[v > 0, , drop=FALSE]
    
    if (nrow(E) > 10 && n_case >= 2 && n_control >= 2) {
      top_n <- min(2000, nrow(E))
      v <- apply(E, 1, var)
      t_genes <- names(sort(v, decreasing=TRUE))[1:top_n]
      p <- prcomp(t(E[t_genes, , drop=FALSE]), scale.=TRUE)
      pc1 <- p$x[, 1]
      pc2 <- p$x[, 2]
      
      case_idx <- meta[[group_col]] == case_lbl
      ctrl_idx <- meta[[group_col]] == ctrl_lbl
      
      # T-test on PC coordinates to check separation statistically
      t1 <- tryCatch(t.test(pc1[case_idx], pc1[ctrl_idx])$p.value, error = function(e) 1.0)
      t2 <- tryCatch(t.test(pc2[case_idx], pc2[ctrl_idx])$p.value, error = function(e) 1.0)
      
      if (t1 < 0.05) {
        pca_sep <- "Yes"
      } else if (t2 < 0.05) {
        pca_sep <- "Partial"
      } else {
        pca_sep <- "No"
      }
    } else {
      pca_sep <- "N/A"
      if (n_case < 3 || n_control < 3) {
        notes <- paste0(notes, "Group size small: case=", n_case, ", ctrl=", n_control, ". ")
      }
    }
  }
  
  if (file.exists(deg_file)) {
    deg <- fread(deg_file)
    expected_cols <- c("gene", "log2FC", "p", "padj", "n")
    if (all(expected_cols %in% colnames(deg)) && ncol(deg) == 5) {
      schema_ok <- "Yes"
    } else {
      schema_ok <- paste0("No (cols: ", paste(colnames(deg), collapse=","), ")")
    }
    
    n_deg <- sum(deg$padj < 0.05, na.rm=TRUE)
    
    # Check if probe codes instead of symbols
    non_symbol_fraction <- sum(grepl("^[0-9]+(_[a-z])?_at$", deg$gene, ignore.case=TRUE)) / nrow(deg)
    if (non_symbol_fraction > 0.5) {
      notes <- paste0(notes, "Gene column contains probe codes. ")
      status <- "drop"
    }
    
    if (n_deg == 0) {
      notes <- paste0(notes, "0 significant DEGs. ")
    }
  } else {
    schema_ok <- "No (file missing)"
    status <- "drop"
  }
  
  results[[ds]] <- data.table(
    dataset_id = ds,
    condition = cond,
    tissue = tissue,
    technique = tech,
    n_case = n_case,
    n_control = n_control,
    n_DEG_padj05 = n_deg,
    pca_separation = pca_sep,
    schema_valid = schema_ok,
    outliers_removed = "No",
    status = status,
    owner = "Aman",
    date = "2026-06-22",
    notes = trimws(notes)
  )
}

dt_results <- rbindlist(results)
fwrite(dt_results, "results/tracking_sheet.tsv", sep="\t")

# Output a nice markdown table to console
cat("\n# QC & Tracking Sheet Summary\n\n")
cat("| dataset_id | technique | n_case | n_control | n_DEGs | PCA sep | Schema Ok | Status | Notes |\n")
cat("|---|---|---|---|---|---|---|---|---|\n")
for (i in 1:nrow(dt_results)) {
  r <- dt_results[i]
  cat(sprintf("| [%s](file:///c:/Users/adity/OneDrive/Documents/gyne-deg-atlas/gyne-deg-atlas/results/deg/%s.deg.tsv) | %s | %d | %d | %d | %s | %s | %s | %s |\n",
              r$dataset_id, r$dataset_id, r$technique, r$n_case, r$n_control, r$n_DEG_padj05, r$pca_separation, r$schema_valid, r$status, r$notes))
}
cat("\n")
 