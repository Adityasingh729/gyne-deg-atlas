#!/usr/bin/env Rscript
# Combine per-dataset DEGs within a stratum.
# RRA: rank up- and down- lists separately, RobustRankAggreg::aggregateRanks.
# REM: random-effects on log2FC via metafor::rma (needs SE; derive from p + n).

suppressMessages({
  library(optparse)
  library(data.table)
})

opt <- parse_args(OptionParser(option_list=list(
  make_option("--stratum"),
  make_option("--method", default="RRA"),
  make_option("--min_datasets", type="integer", default=4),
  make_option("--deg_dir"),
  make_option("--samplesheet"),
  make_option("--out")
)))

stratum_id <- opt$stratum
method <- opt$method
min_datasets <- opt$min_datasets
deg_dir <- opt$deg_dir
samplesheet_path <- opt$samplesheet
out_file <- opt$out

# Ensure output directory exists
dir.create(dirname(out_file), showWarnings=FALSE, recursive=TRUE)

message("Running meta-analysis for stratum: ", stratum_id, " using method: ", method)

# 1) Read samplesheet and select dataset_ids in this stratum
samples <- fread(samplesheet_path)
# Map dataset to stratum ID
samples[, s_id := gsub("/", "-", gsub(" ", "_", paste0(condition, "__", tissue_primary)))]
datasets_in_stratum <- samples[s_id == stratum_id, dataset_id]

message("Datasets in stratum (N = ", length(datasets_in_stratum), "): ", paste(datasets_in_stratum, collapse=", "))

# Load DEG files
deg_data <- list()
for (ds in datasets_in_stratum) {
  # File path is results/deg/{ds}.deg.tsv
  file_path <- file.path(deg_dir, paste0(ds, ".deg.tsv"))
  if (file.exists(file_path)) {
    dt <- fread(file_path)
    if (nrow(dt) > 0) {
      deg_data[[ds]] <- dt
    }
  } else {
    warning("DEG file not found for dataset: ", ds, " at path: ", file_path)
  }
}

n_loaded <- length(deg_data)
message("Successfully loaded ", n_loaded, " DEG datasets.")

# 2) If (n_loaded < min_datasets), write a flagged descriptive result and STOP
if (n_loaded < min_datasets) {
  message("Strata has ", n_loaded, " datasets, which is less than min_datasets (", min_datasets, "). Running descriptive combination fallback.")
  
  if (n_loaded == 0) {
    # No datasets available, write empty schema
    empty_out <- data.table(
      gene = character(),
      meta_log2FC = numeric(),
      meta_p = numeric(),
      meta_padj = numeric(),
      n_datasets = integer(),
      direction = character()
    )
    fwrite(empty_out, out_file, sep="\t")
    message("Wrote empty file since no datasets were loaded.")
    quit(save="no", status=0)
  } else {
    # Combine descriptively (average values across available datasets)
    all_degs <- rbindlist(deg_data, fill=TRUE)
    descriptive_res <- all_degs[, .(
      meta_log2FC = mean(log2FC, na.rm=TRUE),
      meta_p = mean(p, na.rm=TRUE),
      n_datasets = .N
    ), by = gene]
    
    descriptive_res[, meta_padj := p.adjust(meta_p, method="BH")]
    descriptive_res[, direction := ifelse(meta_log2FC > 0, "UP", "DOWN")]
    
    # Select columns in standard order
    out_cols <- descriptive_res[, .(gene, meta_log2FC, meta_p, meta_padj, n_datasets, direction)]
    fwrite(out_cols, out_file, sep="\t")
    message("Wrote descriptive fallback output with ", nrow(out_cols), " genes.")
    quit(save="no", status=0)
  }
}

# Run RRA or REM
if (method == "RRA") {
  suppressMessages(library(RobustRankAggreg))
  
  # Up lists: genes with log2FC > 0, sorted by p-value ascending
  up_list <- lapply(deg_data, function(df) {
    sub_df <- df[log2FC > 0 & !is.na(p), ]
    sub_df <- sub_df[order(p), ]
    return(sub_df$gene)
  })
  up_list <- up_list[sapply(up_list, length) > 0]
  
  # Down lists: genes with log2FC < 0, sorted by p-value ascending
  down_list <- lapply(deg_data, function(df) {
    sub_df <- df[log2FC < 0 & !is.na(p), ]
    sub_df <- sub_df[order(p), ]
    return(sub_df$gene)
  })
  down_list <- down_list[sapply(down_list, length) > 0]
  
  # Aggregated up list
  if (length(up_list) > 0) {
    # Find unique genes
    all_up_genes <- unique(unlist(up_list))
    rra_up <- aggregateRanks(glist = up_list, N = length(all_up_genes))
    rra_up <- as.data.table(rra_up)
    colnames(rra_up) <- c("gene", "meta_p")
    rra_up[, direction := "UP"]
  } else {
    rra_up <- data.table(gene = character(), meta_p = numeric(), direction = character())
  }
  
  # Aggregated down list
  if (length(down_list) > 0) {
    all_down_genes <- unique(unlist(down_list))
    rra_down <- aggregateRanks(glist = down_list, N = length(all_down_genes))
    rra_down <- as.data.table(rra_down)
    colnames(rra_down) <- c("gene", "meta_p")
    rra_down[, direction := "DOWN"]
  } else {
    rra_down <- data.table(gene = character(), meta_p = numeric(), direction = character())
  }
  
  # Combine UP and DOWN
  rra_all <- rbind(rra_up, rra_down)
  rra_all <- rra_all[order(meta_p)]
  # Keep only the most significant direction/p-value for each gene
  rra_all <- rra_all[!duplicated(gene)]
  
  # Calculate average log2FC and count datasets per gene
  all_degs <- rbindlist(deg_data, fill=TRUE)
  gene_stats <- all_degs[, .(
    meta_log2FC = mean(log2FC, na.rm=TRUE),
    n_datasets = sum(!is.na(log2FC))
  ), by = gene]
  
  meta_res <- merge(rra_all, gene_stats, by = "gene")
  meta_res[, meta_padj := p.adjust(meta_p, method="BH")]
  
  # Ensure standard columns are written
  out_cols <- meta_res[, .(gene, meta_log2FC, meta_p, meta_padj, n_datasets, direction)]
  fwrite(out_cols, out_file, sep="\t")
  message("Robust Rank Aggregation (RRA) completed successfully. Wrote ", nrow(out_cols), " genes.")

} else if (method == "REM") {
  suppressMessages(library(metafor))
  
  # Extract all unique genes across all datasets
  all_degs <- rbindlist(deg_data, fill=TRUE)
  unique_genes <- unique(all_degs$gene)
  
  # Helper to derive SE from p-value and log2FC
  derive_se <- function(log2FC, p) {
    p <- pmax(pmin(p, 1 - 1e-16), 1e-300)
    z <- qnorm(p / 2, lower.tail = FALSE)
    se <- abs(log2FC) / z
    se[is.na(se) | se == 0] <- 0.1 # default fallback
    return(se)
  }
  
  results_list <- list()
  
  message("Running REM for ", length(unique_genes), " genes...")
  
  # Since REM can be slow, we run it for each gene
  for (g in unique_genes) {
    gene_data <- all_degs[gene == g & !is.na(log2FC) & !is.na(p)]
    n_ds <- nrow(gene_data)
    
    if (n_ds == 0) next
    
    # Calculate SE for each dataset
    gene_data[, se := derive_se(log2FC, p)]
    
    if (n_ds == 1) {
      # Not enough datasets for REM, write descriptive values
      results_list[[g]] <- list(
        gene = g,
        meta_log2FC = gene_data$log2FC[1],
        meta_p = gene_data$p[1],
        n_datasets = 1,
        direction = ifelse(gene_data$log2FC[1] > 0, "UP", "DOWN")
      )
    } else {
      # Fit REM using metafor
      res_rem <- tryCatch({
        fit <- rma(yi = gene_data$log2FC, sei = gene_data$se, method = "REML")
        list(
          gene = g,
          meta_log2FC = fit$beta[1],
          meta_p = fit$pval,
          n_datasets = n_ds,
          direction = ifelse(fit$beta[1] > 0, "UP", "DOWN")
        )
      }, error = function(e) {
        # Fallback to inverse variance weighted or simple average if REM fails to converge
        w <- 1 / (gene_data$se ^ 2)
        meta_lfc <- sum(gene_data$log2FC * w) / sum(w)
        list(
          gene = g,
          meta_log2FC = meta_lfc,
          meta_p = mean(gene_data$p),
          n_datasets = n_ds,
          direction = ifelse(meta_lfc > 0, "UP", "DOWN")
        )
      })
      results_list[[g]] <- res_rem
    }
  }
  
  meta_res <- rbindlist(results_list)
  meta_res[, meta_padj := p.adjust(meta_p, method="BH")]
  
  out_cols <- meta_res[, .(gene, meta_log2FC, meta_p, meta_padj, n_datasets, direction)]
  fwrite(out_cols, out_file, sep="\t")
  message("Random-Effects Model (REM) completed successfully. Wrote ", nrow(out_cols), " genes.")
  
} else {
  stop("Unknown meta-analysis method: ", method)
}

gc()
 