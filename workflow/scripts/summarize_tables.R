#!/usr/bin/env Rscript
# Script to compile DEG summary tables for PCOS, Endometriosis, and Adenomyosis master lists

suppressMessages({
  library(data.table)
})

# Read samplesheet to cross-reference mapped details
samplesheet_path <- "config/samplesheet.csv"
samplesheet_exists <- file.exists(samplesheet_path)
samples_map <- list()
if (samplesheet_exists) {
  ss <- fread(samplesheet_path)
  for (i in 1:nrow(ss)) {
    samples_map[[ss$dataset_id[i]]] <- ss[i]
  }
}

parse_disease_csv <- function(file_path, condition_name) {
  if (!file.exists(file_path)) {
    warning("File not found: ", file_path)
    return(NULL)
  }
  
  lines <- readLines(file_path)
  results <- list()
  current_tech <- "RNA-seq" # default starting
  
  for (line in lines) {
    # Detect technique sections
    if (grepl("Microarray", line, ignore.case=TRUE)) {
      current_tech <- "Microarray"
      next
    }
    if (grepl("RNA-seq", line, ignore.case=TRUE)) {
      current_tech <- "RNA-seq"
      next
    }
    
    # Split CSV line (handling quotes properly)
    parts <- tryCatch({
      # Simple CSV parser
      res <- utils::read.csv(text=line, header=FALSE, stringsAsFactors=FALSE)
      if (nrow(res) > 0) as.character(res[1, ]) else NULL
    }, error = function(e) NULL)
    
    if (is.null(parts) || length(parts) < 1) next
    
    accession <- trimws(parts[1])
    if (grepl("^GSE[0-9]+", accession)) {
      title <- if (length(parts) >= 2) trimws(parts[2]) else ""
      samples_est <- if (length(parts) >= 3) trimws(parts[3]) else ""
      tissue <- if (length(parts) >= 4) trimws(parts[4]) else ""
      flag <- if (length(parts) >= 5) trimws(parts[5]) else ""
      
      # Cross-reference with samplesheet
      n_case <- NA
      n_control <- NA
      n_deg <- NA
      status <- "unmapped"
      notes <- flag
      
      # Check if exists in samplesheet map
      if (accession %in% names(samples_map)) {
        s_row <- samples_map[[accession]]
        tissue <- s_row$tissue_primary
        current_tech <- s_row$technique
        
        # Check expression RDS
        expr_file <- file.path("data", "expr", paste0(accession, ".expr.rds"))
        if (file.exists(expr_file)) {
          x <- readRDS(expr_file)
          meta <- x$meta
          group_col <- s_row$group_column
          case_lbl <- s_row$case_label
          ctrl_lbl <- s_row$control_label
          
          if (group_col %in% colnames(meta)) {
            n_case <- sum(meta[[group_col]] == case_lbl, na.rm=TRUE)
            n_control <- sum(meta[[group_col]] == ctrl_lbl, na.rm=TRUE)
          }
        }
        
        # Check DEG table
        deg_file <- file.path("results", "deg", paste0(accession, ".deg.tsv"))
        if (file.exists(deg_file)) {
          deg <- fread(deg_file)
          status <- "analyzed"
          if ("padj" %in% colnames(deg)) {
            n_deg <- sum(deg$padj < 0.05, na.rm=TRUE)
          } else {
            notes <- paste0(notes, " | Missing padj column.")
          }
        } else {
          status <- "missing deg"
        }
      } else {
        # Check if files exist anyway
        expr_file <- file.path("data", "expr", paste0(accession, ".expr.rds"))
        deg_file <- file.path("results", "deg", paste0(accession, ".deg.tsv"))
        if (file.exists(deg_file)) {
          status <- "analyzed (not in samplesheet)"
          deg <- fread(deg_file)
          if ("padj" %in% colnames(deg)) n_deg <- sum(deg$padj < 0.05, na.rm=TRUE)
        } else if (file.exists(expr_file)) {
          status <- "raw data fetched"
        }
      }
      
      results[[accession]] <- data.table(
        dataset_id = accession,
        condition = condition_name,
        tissue = tissue,
        technique = current_tech,
        samples_estimate = samples_est,
        n_case_analyzed = n_case,
        n_control_analyzed = n_control,
        n_DEG_padj05 = n_deg,
        status = status,
        notes = notes
      )
    }
  }
  
  if (length(results) == 0) return(NULL)
  return(rbindlist(results))
}

# Try relative to project root first, then relative to parent (workspace root)
find_file <- function(fn) {
  if (file.exists(fn)) return(fn)
  parent_fn <- file.path("..", fn)
  if (file.exists(parent_fn)) return(parent_fn)
  # Try absolute path based on workspace location
  abs_fn <- file.path("C:", "Users", "adity", "OneDrive", "Documents", "gyne-deg-atlas", fn)
  if (file.exists(abs_fn)) return(abs_fn)
  return(fn)
}

# Generate summaries for each table
pcos_summary <- parse_disease_csv(find_file("PCOS.csv"), "PCOS")
endo_summary <- parse_disease_csv(find_file("Endometriosis.csv"), "Endometriosis")
adeno_summary <- parse_disease_csv(find_file("Adenomyosis.csv"), "Adenomyosis")

# Write output files
if (!is.null(pcos_summary)) {
  fwrite(pcos_summary, "results/PCOS_deg_summary.tsv", sep="\t")
  fwrite(pcos_summary, "results/PCOS_deg_summary.csv", sep=",")
  message("Generated results/PCOS_deg_summary.csv and .tsv")
}
if (!is.null(endo_summary)) {
  fwrite(endo_summary, "results/Endometriosis_deg_summary.tsv", sep="\t")
  fwrite(endo_summary, "results/Endometriosis_deg_summary.csv", sep=",")
  message("Generated results/Endometriosis_deg_summary.csv and .tsv")
}
if (!is.null(adeno_summary)) {
  fwrite(adeno_summary, "results/Adenomyosis_deg_summary.tsv", sep="\t")
  fwrite(adeno_summary, "results/Adenomyosis_deg_summary.csv", sep=",")
  message("Generated results/Adenomyosis_deg_summary.csv and .tsv")
}
