#!/usr/bin/env Rscript
# Standalone R script to run DEG analysis for all 13 PCOS Granulosa datasets
# following the manual GEO-first workflow guidelines.

suppressMessages({
  library(DESeq2)
  library(limma)
  library(GEOquery)
  library(org.Hs.eg.db)
  library(data.table)
})

# Set working directory to project root
setwd("C:/Users/adity/OneDrive/Documents/gyne-deg-atlas/gyne-deg-atlas")

# Output directory
deg_dir <- "results/deg"
dir.create(deg_dir, showWarnings=FALSE, recursive=TRUE)

# Read samplesheet
samples <- fread("config/samplesheet.csv")
ds_list <- c("GSE34526", "GSE168404", "GSE216609", "GSE328499", "GSE294074", 
             "GSE191022", "GSE193123", "GSE98595", "GSE80432", "GSE106724", 
             "GSE137684", "GSE102293", "GSE114419")

# Helper to clean and format output table
write_output <- function(res_df, base_path, symbol_mapping=NULL) {
  # Columns: gene, log2FC, p, padj, n
  # Sorted by p (smallest first)
  res_df <- res_df[!is.na(res_df$p), ]
  res_df <- res_df[order(res_df$p), ]
  
  # Write primary table (Ensembl for RNA-seq, Probe ID for microarray)
  fwrite(res_df[, c("gene", "log2FC", "p", "padj", "n")], paste0(base_path, ".deg.tsv"), sep="\t")
  
  # Map and write symbol table
  if (!is.null(symbol_mapping)) {
    # symbol_mapping is a named vector mapping gene_id to symbol
    res_df$symbol <- symbol_mapping[as.character(res_df$gene)]
    # Filter empty or NA symbols
    res_df_sym <- res_df[!is.na(res_df$symbol) & res_df$symbol != "" & res_df$symbol != "---", ]
    
    # Resolve duplicate symbols by keeping the one with max absolute log2FC
    res_df_sym$abs_logFC <- abs(res_df_sym$log2FC)
    res_df_sym <- res_df_sym[order(res_df_sym$symbol, -res_df_sym$abs_logFC), ]
    res_df_sym <- res_df_sym[!duplicated(res_df_sym$symbol), ]
    
    # Replace gene column with symbol
    res_df_sym$gene <- res_df_sym$symbol
    res_df_sym <- res_df_sym[order(res_df_sym$p), ]
    
    fwrite(res_df_sym[, c("gene", "log2FC", "p", "padj", "n")], paste0(base_path, ".deg.symbols.tsv"), sep="\t")
  }
}

# ----------------- RUN MICROARRAY (limma) -----------------
run_microarray_dataset <- function(ds) {
  message("\n--- Running Microarray: ", ds, " ---")
  r <- samples[dataset_id == ds]
  
  # Load gset using local series matrix if present, otherwise download
  # destdir="data/expr" allows GEOquery to cache/read GPL and matrix files there
  g <- getGEO(ds, GSEMatrix=TRUE, destdir="data/expr")
  gset <- g[[1]]
  
  E <- exprs(gset)
  meta <- as.data.frame(pData(gset))
  
  # Group column, case and control labels
  col <- r$group_column
  case <- r$case_label
  ctrl <- r$control_label
  
  # Handle special case for GSE34526 (which has disease:ch1 or characteristics_ch1.1)
  if (ds == "GSE34526") {
    if (!col %in% colnames(meta)) {
      col_source <- if ("disease:ch1" %in% colnames(meta)) "disease:ch1" else "characteristics_ch1.1"
      meta[[col]] <- ifelse(grepl("PCOS|Polycystic", meta[[col_source]], ignore.case=TRUE), "PCOS", "Control")
    }
  }
  
  # Normalize labels (e.g. Mu to u)
  clean_lbl <- function(s) {
    s <- gsub("\u00B5", "u", s, fixed=TRUE)
    s <- gsub("\u03BC", "u", s, fixed=TRUE)
    return(s)
  }
  meta[[col]] <- clean_lbl(as.character(meta[[col]]))
  case <- clean_lbl(case)
  ctrl <- clean_lbl(ctrl)
  
  # Subset to case and control
  keep <- !is.na(meta[[col]]) & (meta[[col]] %in% c(case, ctrl))
  if (sum(keep) < 2) {
    stop("Fewer than 2 samples matched grouping.")
  }
  meta_sub <- meta[keep, , drop=FALSE]
  E_sub <- E[, keep, drop=FALSE]
  
  # Log2 transform if necessary
  qx <- as.numeric(quantile(E_sub, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=TRUE))
  LogC <- (qx[5] > 100) || (qx[6] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
  if (LogC) {
    E_sub[E_sub <= 0] <- NaN
    E_sub <- log2(E_sub)
  }
  
  # Design and Fit
  grp <- factor(meta_sub[[col]], levels=c(ctrl, case))
  message("Group sizes: ")
  print(table(grp))
  
  design <- model.matrix(~grp)
  fit <- lmFit(E_sub, design)
  fit <- eBayes(fit)
  tt <- topTable(fit, coef=2, number=Inf)
  
  # Create result data frame
  res <- data.frame(
    gene = rownames(tt),
    log2FC = tt$logFC,
    p = tt$P.Value,
    padj = tt$adj.P.Val,
    n = ncol(E_sub),
    stringsAsFactors = FALSE
  )
  
  # Probe-to-Symbol mapping
  fdata <- fData(gset)
  symbol_map <- NULL
  
  if ("gene_assignment" %in% colnames(fdata)) {
    # Extract symbols from gene_assignment (e.g. Affymetrix Gene 1.0 ST)
    symbols <- fdata$gene_assignment
    # Split by /// and take the first mapping
    symbols <- sapply(strsplit(as.character(symbols), " /// ", fixed=TRUE), `[`, 1)
    # Split by // and take the second token (gene symbol)
    symbols <- sapply(strsplit(symbols, " // ", fixed=TRUE), `[`, 2)
    symbols <- trimws(symbols)
    symbols[symbols == "" | symbols == "---" | symbols == "NULL"] <- NA
    symbol_map <- setNames(symbols, rownames(fdata))
  } else {
    symbol_col <- NULL
    possible_cols <- c("Gene Symbol", "gene_symbol", "GENE_SYMBOL", "SYMBOL", "Gene.Symbol", "GENE", "Gene Symbol (RefSeq)", "GeneName", "Gene.Name")
    for (pc in possible_cols) {
      if (pc %in% colnames(fdata)) {
        symbol_col <- pc
        break
      }
    }
    
    if (!is.null(symbol_col)) {
      symbols <- fdata[[symbol_col]]
      symbols <- sapply(strsplit(as.character(symbols), " /// "), `[`, 1)
      symbols <- sapply(strsplit(as.character(symbols), "//"), `[`, 1)
      symbols <- trimws(symbols)
      symbol_map <- setNames(symbols, rownames(fdata))
    } else {
      warning("Could not find gene symbol column in fData for ", ds)
    }
  }
  
  write_output(res, file.path(deg_dir, ds), symbol_map)
  message("Completed Microarray: ", ds)
}

# ----------------- RUN RNA-SEQ (DESeq2) -----------------
run_rnaseq_dataset <- function(ds) {
  message("\n--- Running RNA-seq: ", ds, " ---")
  r <- samples[dataset_id == ds]
  
  # Load metadata from RDS
  x <- readRDS(file.path("data", "expr", paste0(ds, ".expr.rds")))
  meta <- as.data.frame(x$meta)
  
  col <- r$group_column
  case <- r$case_label
  ctrl <- r$control_label
  
  clean_lbl <- function(s) {
    s <- gsub("\u00B5", "u", s, fixed=TRUE)
    s <- gsub("\u03BC", "u", s, fixed=TRUE)
    return(s)
  }
  meta[[col]] <- clean_lbl(as.character(meta[[col]]))
  case <- clean_lbl(case)
  ctrl <- clean_lbl(ctrl)
  
  # Subset metadata to case and control
  keep_meta <- !is.na(meta[[col]]) & (meta[[col]] %in% c(case, ctrl))
  meta_sub <- meta[keep_meta, , drop=FALSE]
  
  # Load raw counts from parent workspace directory
  counts <- NULL
  
  if (ds == "GSE191022") {
    # File: GSE191022_expressed_gene_reads.txt.gz
    cnt <- read.delim(gzfile("GSE191022_expressed_gene_reads.txt.gz"), check.names=FALSE)
    # Extract version-less ENSG IDs
    ensg_ids <- sub("\\..*", "", cnt$Gene)
    
    # Columns are NC_1, NC_2, NC_3, OE_IGF2BP2_1, OE_IGF2BP2_2, OE_IGF2BP2_3
    # Map to GSM IDs using meta titles
    # GSM5737627 = NC 1 -> NC_1, GSM5737630 = OE IGF2BP2 1 -> OE_IGF2BP2_1, etc.
    # NC samples are controls, OE are cases
    map_cols <- c(
      "GSM5737627" = "NC_1", "GSM5737628" = "NC_2", "GSM5737629" = "NC_3",
      "GSM5737630" = "OE_IGF2BP2_1", "GSM5737631" = "OE_IGF2BP2_2", "GSM5737632" = "OE_IGF2BP2_3"
    )
    counts_raw <- as.matrix(cnt[, map_cols])
    
    # Aggregate duplicate Ensembl Gene IDs by summing counts
    counts_dt <- data.table(gene = ensg_ids, counts_raw)
    counts_agg <- counts_dt[, lapply(.SD, sum), by = gene]
    counts <- as.matrix(counts_agg[, -1, with=FALSE])
    rownames(counts) <- counts_agg$gene
    colnames(counts) <- names(map_cols)
    
    # Symbol mapping
    ensg_to_sym <- setNames(cnt$Symbol, ensg_ids)
    symbol_map <- ensg_to_sym[rownames(counts)]
    
  } else if (ds == "GSE193123") {
    # File: GSE193123_gene_count.txt.gz
    cnt <- read.delim(gzfile("GSE193123_gene_count.txt.gz"), check.names=FALSE)
    rownames(cnt) <- cnt$gene_id
    
    # Columns: C1, C2, C3, P1, P2, P3
    # Map to GSM IDs
    # Control/Healthy: GSM5773736 = C1, GSM5773737 = C2, GSM5773738 = C3
    # Case/PCOS: GSM5773739 = P1, GSM5773740 = P2, GSM5773741 = P3
    map_cols <- c(
      "GSM5773736" = "C1", "GSM5773737" = "C2", "GSM5773738" = "C3",
      "GSM5773739" = "P1", "GSM5773740" = "P2", "GSM5773741" = "P3"
    )
    counts <- as.matrix(cnt[, map_cols])
    colnames(counts) <- names(map_cols)
    
    # Symbol mapping
    symbol_map <- setNames(cnt$gene_name, rownames(cnt))
    
  } else if (ds == "GSE216609") {
    # File: GSE216609_compiled_counts.tsv.gz
    cnt <- read.delim(gzfile("GSE216609_compiled_counts.tsv.gz"), check.names=FALSE)
    # gene has transcript IDs starting with ENST
    enst_ids <- sub("\\..*", "", cnt$gene)
    counts <- as.matrix(cnt[, 2:ncol(cnt)])
    
    # Map ENST transcript IDs to Ensembl Gene IDs and Symbols using org.Hs.eg.db
    map_df <- tryCatch({
      select(org.Hs.eg.db, keys=enst_ids, columns=c("ENSEMBL", "SYMBOL"), keytype="ENSEMBLTRANS")
    }, error = function(e) NULL)
    
    if (!is.null(map_df)) {
      # Match transcripts to genes
      ensg_map <- setNames(map_df$ENSEMBL, map_df$ENSEMBLTRANS)
      sym_map <- setNames(map_df$SYMBOL, map_df$ENSEMBLTRANS)
      
      # Assign ENSG as row names, resolving duplicates by summing counts
      ensg_ids <- ensg_map[enst_ids]
      valid_idx <- !is.na(ensg_ids)
      counts_sub <- counts[valid_idx, , drop=FALSE]
      ensg_ids_sub <- ensg_ids[valid_idx]
      
      counts_dt <- data.table(gene = ensg_ids_sub, counts_sub)
      counts_agg <- counts_dt[, lapply(.SD, sum), by = gene]
      counts <- as.matrix(counts_agg[, -1, with=FALSE])
      rownames(counts) <- counts_agg$gene
      
      # Map ENSG to Symbol for the secondary file
      map_df_valid <- map_df[!is.na(map_df$ENSEMBL) & !is.na(map_df$SYMBOL), ]
      symbol_map <- setNames(map_df_valid$SYMBOL, map_df_valid$ENSEMBL)
    } else {
      # Fallback
      rownames(counts) <- enst_ids
      symbol_map <- setNames(enst_ids, enst_ids)
    }
    
  } else if (ds == "GSE328499") {
    # File: GSE328499_all.anno_converted.tsv
    cnt <- read.delim("GSE328499_all.anno_converted.tsv", check.names=FALSE)
    rownames(cnt) <- cnt$Gene
    
    # Columns are Ctrl_1, Ctrl_2, Ctrl_3 (controls), IAA_1, IAA_2, IAA_3 (cases)
    # Map to GSM IDs using meta titles
    # GSM9684461 = Ctrl_1, GSM9684462 = Ctrl_2, GSM9684463 = Ctrl_3
    # GSM9684464 = IAA_1, GSM9684465 = IAA_2, GSM9684466 = IAA_3
    map_cols <- c(
      "GSM9684461" = "Ctrl_1", "GSM9684462" = "Ctrl_2", "GSM9684463" = "Ctrl_3",
      "GSM9684464" = "IAA_1", "GSM9684465" = "IAA_2", "GSM9684466" = "IAA_3"
    )
    counts <- as.matrix(cnt[, map_cols])
    colnames(counts) <- names(map_cols)
    
    # Map Symbols to Ensembl Gene IDs
    symbols <- rownames(counts)
    map_df <- tryCatch({
      select(org.Hs.eg.db, keys=symbols, columns="ENSEMBL", keytype="SYMBOL")
    }, error = function(e) NULL)
    
    if (!is.null(map_df)) {
      ensg_map <- setNames(map_df$ENSEMBL, map_df$SYMBOL)
      ensg_ids <- ensg_map[symbols]
      
      # Keep rows with valid Ensembl IDs
      valid_idx <- !is.na(ensg_ids)
      counts_sub <- counts[valid_idx, , drop=FALSE]
      ensg_ids_sub <- ensg_ids[valid_idx]
      symbols_sub <- symbols[valid_idx]
      
      # Resolve duplicates
      counts_dt <- data.table(gene = ensg_ids_sub, counts_sub)
      counts_agg <- counts_dt[, lapply(.SD, sum), by = gene]
      counts <- as.matrix(counts_agg[, -1, with=FALSE])
      rownames(counts) <- counts_agg$gene
      
      symbol_map <- setNames(symbols_sub, ensg_ids_sub)
    } else {
      symbol_map <- setNames(symbols, symbols)
    }
    
  } else if (ds == "GSE294074") {
    # File: GSE294074_compiled_counts.tsv.gz
    cnt <- read.delim(gzfile("GSE294074_compiled_counts.tsv.gz"), check.names=FALSE)
    rownames(cnt) <- cnt$gene
    counts <- as.matrix(cnt[, 2:ncol(cnt)])
    
    # Read annotation file
    anno <- read.delim(gzfile("GSE294074_annotation.xls.gz"))
    pb_ids <- sub("\\(.*", "", anno$Unigene)
    symbols <- gsub(".*GN=([^ ]+).*", "\\1", anno$Swissprot)
    symbols[!grepl("GN=", anno$Swissprot)] <- NA
    
    pb_to_sym <- setNames(symbols, pb_ids)
    
    # Map PacBio IDs to Symbols
    pb_ids_counts <- rownames(counts)
    mapped_symbols <- pb_to_sym[pb_ids_counts]
    
    # Map Symbols to Ensembl Gene IDs
    valid_idx <- !is.na(mapped_symbols)
    counts_sub <- counts[valid_idx, , drop=FALSE]
    mapped_symbols_sub <- mapped_symbols[valid_idx]
    
    map_df <- tryCatch({
      select(org.Hs.eg.db, keys=unique(mapped_symbols_sub), columns="ENSEMBL", keytype="SYMBOL")
    }, error = function(e) NULL)
    
    if (!is.null(map_df)) {
      ensg_map <- setNames(map_df$ENSEMBL, map_df$SYMBOL)
      ensg_ids <- ensg_map[mapped_symbols_sub]
      
      valid_ensg <- !is.na(ensg_ids)
      counts_final <- counts_sub[valid_ensg, , drop=FALSE]
      ensg_ids_final <- ensg_ids[valid_ensg]
      symbols_final <- mapped_symbols_sub[valid_ensg]
      
      # Resolve duplicate Ensembl IDs by summing counts
      counts_dt <- data.table(gene = ensg_ids_final, counts_final)
      counts_agg <- counts_dt[, lapply(.SD, sum), by = gene]
      counts <- as.matrix(counts_agg[, -1, with=FALSE])
      rownames(counts) <- counts_agg$gene
      
      symbol_map <- setNames(symbols_final, ensg_ids_final)
    } else {
      symbol_map <- setNames(mapped_symbols_sub, mapped_symbols_sub)
    }
    
  } else if (ds == "GSE168404") {
    # File: GSE168404_mRNA_C-vs-P.all.txt.gz (UTF-16)
    con <- gzfile("GSE168404_mRNA_C-vs-P.all.txt.gz", "rt", encoding="UTF-16")
    cnt <- read.delim(con, check.names=FALSE)
    close(con)
    # The second column is the gene symbol, rename to Symbol
    colnames(cnt)[2] <- "Symbol"
    
    # First column 'id' is Entrez ID
    entrez_ids <- as.character(cnt$id)
    counts <- as.matrix(cnt[, 3:ncol(cnt)]) # columns are C1-C5, P1-P5
    
    # Map Entrez IDs to Ensembl Gene IDs
    map_df <- tryCatch({
      select(org.Hs.eg.db, keys=entrez_ids, columns="ENSEMBL", keytype="ENTREZID")
    }, error = function(e) NULL)
    
    # Map columns to GSM IDs
    # C1-C5 -> GSM5138304-GSM5138308
    # P1-P5 -> GSM5138299-GSM5138303
    map_cols <- c(
      "GSM5138304" = "C1", "GSM5138305" = "C2", "GSM5138306" = "C3", "GSM5138307" = "C4", "GSM5138308" = "C5",
      "GSM5138299" = "P1", "GSM5138300" = "P2", "GSM5138301" = "P3", "GSM5138302" = "P4", "GSM5138303" = "P5"
    )
    counts <- counts[, map_cols]
    colnames(counts) <- names(map_cols)
    
    if (!is.null(map_df)) {
      ensg_map <- setNames(map_df$ENSEMBL, map_df$ENTREZID)
      ensg_ids <- ensg_map[entrez_ids]
      
      valid_idx <- !is.na(ensg_ids)
      counts_sub <- counts[valid_idx, , drop=FALSE]
      ensg_ids_sub <- ensg_ids[valid_idx]
      symbols_sub <- cnt$Symbol[valid_idx]
      
      # Resolve duplicate Ensembl IDs
      counts_dt <- data.table(gene = ensg_ids_sub, counts_sub)
      counts_agg <- counts_dt[, lapply(.SD, sum), by = gene]
      counts <- as.matrix(counts_agg[, -1, with=FALSE])
      rownames(counts) <- counts_agg$gene
      
      symbol_map <- setNames(symbols_sub, ensg_ids_sub)
    } else {
      rownames(counts) <- cnt$Symbol
      symbol_map <- setNames(cnt$Symbol, cnt$Symbol)
    }
  }
  
  # Align counts and metadata
  common_samples <- intersect(colnames(counts), rownames(meta_sub))
  counts <- counts[, common_samples, drop=FALSE]
  meta_sub <- meta_sub[common_samples, , drop=FALSE]
  
  # Force counts to integer mode and round
  counts <- round(counts)
  mode(counts) <- "integer"
  
  # DESeq2 setup
  grp <- factor(meta_sub[[col]], levels=c(ctrl, case))
  message("Group sizes: ")
  print(table(grp))
  
  # Low count filtering
  keep <- rowSums(counts >= 10) >= max(2, floor(0.2*ncol(counts)))
  
  # Rename column to safe name to prevent formula/design issues with colons/spaces in col
  meta_sub$group_DESeq2 <- grp
  dds <- DESeqDataSetFromMatrix(counts[keep, ], meta_sub, ~ group_DESeq2)
  dds <- DESeq(dds)
  res_deseq <- as.data.frame(results(dds))
  
  res <- data.frame(
    gene = rownames(res_deseq),
    log2FC = res_deseq$log2FoldChange,
    p = res_deseq$pvalue,
    padj = res_deseq$padj,
    n = ncol(counts),
    stringsAsFactors = FALSE
  )
  
  write_output(res, file.path(deg_dir, ds), symbol_map)
  message("Completed RNA-seq: ", ds)
}

# ----------------- EXECUTE ALL 13 DATASETS -----------------
for (ds in ds_list) {
  tech <- samples[dataset_id == ds, technique]
  tryCatch({
    if (tech == "Microarray") {
      run_microarray_dataset(ds)
    } else if (tech == "RNA-seq") {
      run_rnaseq_dataset(ds)
    }
  }, error = function(e) {
    message("ERROR running dataset ", ds, ": ", e$message)
  })
}

message("\nAll 13 DEG analyses executed successfully!")
