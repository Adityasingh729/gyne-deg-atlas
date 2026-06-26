#!/usr/bin/env Rscript
# Pull analysis-ready expression for ONE dataset and save a tidy RDS and sidecar TSVs:
#   list(counts=matrix genes x samples, meta=data.frame samples x covariates)
# RNA-seq  -> recount3 (preferred) or local fallback
# Microarray -> GEOquery getGEO() processed matrix
# 8 GB: load ONE dataset at a time; gc() after.

suppressMessages({
  library(optparse)
  library(data.table)
})

opt <- parse_args(OptionParser(option_list=list(
  make_option("--dataset"),
  make_option("--technique"),
  make_option("--source"),
  make_option("--out")
)))

dataset <- opt$dataset
tech <- opt$technique
source <- opt$source
out_file <- opt$out

# Ensure output directory exists
dir.create(dirname(out_file), showWarnings=FALSE, recursive=TRUE)

message("Fetching dataset: ", dataset, " (", tech, ") from source: ", source)

counts <- NULL
meta <- NULL

if (grepl("RNA-seq", tech, ignore.case=TRUE)) {
  # First try recount3
  use_recount3 <- FALSE
  if (source == "recount3") {
    suppressMessages(library(recount3))
    ap <- available_projects()
    proj_info <- subset(ap, project == dataset & project_type == "data_sources")
    if (nrow(proj_info) == 0) {
      proj_info <- subset(ap, project == dataset)
    }
    if (nrow(proj_info) > 0) {
      use_recount3 <- TRUE
    }
  }
  
  if (use_recount3) {
    message("Using recount3 to fetch project: ", dataset)
    suppressMessages(library(recount3))
    rse <- create_rse(proj_info[1, ])
    counts <- as.matrix(compute_read_counts(rse))
    meta <- as.data.frame(colData(rse))
    
    # Map Ensembl IDs to symbols
    if ("gene_name" %in% colnames(rowData(rse))) {
      rownames(counts) <- rowData(rse)$gene_name
    } else if ("gene_id" %in% colnames(rowData(rse))) {
      rownames(counts) <- rowData(rse)$gene_id
    }
    
    # Filter empty or NA gene names
    valid <- !is.na(rownames(counts)) & rownames(counts) != ""
    counts <- counts[valid, , drop=FALSE]
    
    # Resolve duplicates
    if (any(duplicated(rownames(counts)))) {
      counts_dt <- data.table(gene = rownames(counts), counts)
      counts_agg <- counts_dt[, lapply(.SD, sum), by = gene]
      counts <- as.matrix(counts_agg[, -1, with=FALSE])
      rownames(counts) <- counts_agg$gene
    }
    
  } else {
    # Fallback/salmon_workstation: Check/download supplementary files from GEO
    
    # Helper to find a RAW.tar file
    find_raw_tar <- function(ds) {
      matches <- list.files(pattern = paste0("^", ds, ".*RAW\\.tar$"), ignore.case=TRUE)
      if (length(matches) > 0) return(matches[1])
      return(NULL)
    }
    
    # Helper to find a count matrix file (must contain 'count' or 'matrix' or 'read')
    find_local_count_file <- function(ds) {
      patterns <- c(
        paste0("^", ds, "_.*count.*\\.(tsv|csv|txt|tab)(\\.gz)?$"),
        paste0("^", ds, "_.*matrix.*\\.(tsv|csv|txt|tab)(\\.gz)?$"),
        paste0("^", ds, "_.*read.*\\.(tsv|csv|txt|tab)(\\.gz)?$")
      )
      for (pat in patterns) {
        matches <- list.files(pattern = pat, ignore.case=TRUE)
        if (length(matches) > 0) return(matches[1])
      }
      return(NULL)
    }
    
    compiled_file <- paste0(dataset, "_compiled_counts.tsv.gz")
    local_file <- NULL
    if (file.exists(compiled_file)) {
      local_file <- compiled_file
    } else {
      local_file <- find_local_count_file(dataset)
    }
    
    # Try fallback to any matching file if local_file is still NULL (before downloading)
    if (is.null(local_file)) {
      patterns <- c(
        paste0("^", dataset, ".*\\.(tsv|csv|txt|tab)(\\.gz)?$"),
        paste0("^", dataset, ".*\\.xlsx?(\\.gz)?$")
      )
      for (pat in patterns) {
        matches <- list.files(pattern = pat, ignore.case=TRUE)
        matches <- matches[!grepl("RAW\\.tar$|meta|clinical|sample|covariate|pheno|results|deg|annotation|anno", matches, ignore.case=TRUE)]
        if (length(matches) > 0) {
          local_file <- matches[1]
          break
        }
      }
    }
    
    raw_tar <- find_raw_tar(dataset)
    
    if (is.null(local_file) && is.null(raw_tar)) {
      message("Downloading supplementary files for ", dataset, " from GEO...")
      suppressMessages(library(GEOquery))
      getGEOSuppFiles(dataset, makeDirectory=FALSE)
      local_file <- find_local_count_file(dataset)
      
      # If still NULL after download, try fallback again
      if (is.null(local_file)) {
        for (pat in patterns) {
          matches <- list.files(pattern = pat, ignore.case=TRUE)
          matches <- matches[!grepl("RAW\\.tar$|meta|clinical|sample|covariate|pheno|results|deg|annotation|anno", matches, ignore.case=TRUE)]
          if (length(matches) > 0) {
            local_file <- matches[1]
            break
          }
        }
      }
      raw_tar <- find_raw_tar(dataset)
    }
    
    # Compile raw counts from RAW.tar if present
    if (is.null(local_file) && !is.null(raw_tar)) {
      message("Found RAW.tar archive: ", raw_tar, ". Unpacking and compiling counts matrix...")
      untar_dir <- file.path(tempdir(), paste0(dataset, "_untar"))
      dir.create(untar_dir, showWarnings=FALSE, recursive=TRUE)
      untar(raw_tar, exdir=untar_dir)
      
      files <- list.files(untar_dir, full.names=TRUE)
      if (length(files) > 0) {
        sample_counts <- list()
        for (f in files) {
          gsm_id <- regmatches(basename(f), regexpr("GSM[0-9]+", basename(f)))
          if (length(gsm_id) > 0) {
            if (grepl("\\.(tsv|csv|txt|tab|counts|xls|xlsx)(\\.gz)?$", f, ignore.case=TRUE)) {
              dt_s <- tryCatch({
                fread(f)
              }, error = function(e) NULL)
              if (!is.null(dt_s) && ncol(dt_s) >= 2) {
                col_names_lower <- tolower(colnames(dt_s))
                gene_col_idx <- grep("gene|target_id|transcript|id|name", col_names_lower)
                if (length(gene_col_idx) == 0) gene_col_idx <- 1 else gene_col_idx <- gene_col_idx[1]
                
                count_col_idx <- grep("count|read|expected", col_names_lower)
                 if (length(count_col_idx) == 0) {
                   # Fallback to expression/tpm/fpkm if count/read not found
                   count_col_idx <- grep("tpm|fpkm|rpkm|expression|value", col_names_lower)
                 }
                 if (length(count_col_idx) == 0) count_col_idx <- 2 else count_col_idx <- count_col_idx[1]
                
                dt_clean <- data.table(
                  gene = as.character(dt_s[[gene_col_idx]]),
                  count = as.numeric(dt_s[[count_col_idx]])
                )
                dt_clean <- dt_clean[!grepl("^__", gene)] # remove htseq comments
                dt_clean <- dt_clean[!is.na(count) & !is.na(gene) & gene != ""]
                sample_counts[[gsm_id]] <- dt_clean
              }
            }
          }
        }
        
        if (length(sample_counts) > 0) {
          merged_dt <- NULL
          for (gsm in names(sample_counts)) {
            dt_s <- sample_counts[[gsm]]
            setnames(dt_s, "count", gsm)
            if (is.null(merged_dt)) {
              merged_dt <- dt_s
            } else {
              merged_dt <- merge(merged_dt, dt_s, by="gene", all=TRUE)
            }
          }
          # Replace any NAs with 0
          for (col in names(merged_dt)) {
            if (col != "gene") {
              merged_dt[[col]][is.na(merged_dt[[col]])] <- 0
            }
          }
          local_file <- compiled_file
          fwrite(merged_dt, local_file, sep="\t", compress="gzip")
          message("Successfully compiled counts from RAW.tar into: ", local_file)
        }
      }
    }
    
    # Fallback to any matching file if local_file is still NULL (e.g. no counts in the name, but it is the only supplementary file)
    if (is.null(local_file)) {
      patterns <- c(
        paste0("^", dataset, ".*\\.(tsv|csv|txt|tab)(\\.gz)?$"),
        paste0("^", dataset, ".*\\.xlsx?(\\.gz)?$")
      )
      for (pat in patterns) {
        matches <- list.files(pattern = pat, ignore.case=TRUE)
        matches <- matches[!grepl("RAW\\.tar$|meta|clinical|sample|covariate|pheno|results|deg|annotation|anno", matches, ignore.case=TRUE)]
        if (length(matches) > 0) {
          local_file <- matches[1]
          break
        }
      }
    }
    
    # If the local file is an Excel file, convert it to TSV via Python/pandas
    if (!is.null(local_file) && grepl("\\.xlsx?$", local_file, ignore.case=TRUE)) {
      message("Found Excel file: ", local_file, ". Converting to TSV via Python...")
      tsv_file <- gsub("\\.xlsx?$", "_converted.tsv", local_file, ignore.case=TRUE)
      
      clean_local_file <- normalizePath(local_file, winslash="/")
      clean_tsv_file <- normalizePath(tsv_file, winslash="/", mustWork=FALSE)
      
      py_file <- file.path(tempdir(), "convert.py")
      cat(paste0(
        "import pandas as pd\n",
        "try:\n",
        "    df = pd.read_excel('", clean_local_file, "')\n",
        "    header_idx = 0\n",
        "    for i in range(min(5, len(df))):\n",
        "        row_vals = [str(x).lower() for x in df.iloc[i].values]\n",
        "        if any(x == 'gene' or x == 'gene_id' or x == 'gene symbol' for x in row_vals):\n",
        "            header_idx = i + 1; break\n",
        "    df = pd.read_excel('", clean_local_file, "', header=header_idx)\n",
        "    df.to_csv('", clean_tsv_file, "', sep='\\t', index=False)\n",
        "except Exception as e:\n",
        "    print(f'Error converting Excel: {e}')\n"
      ), file=py_file)
      system(paste0("python ", py_file))
      
      if (file.exists(tsv_file)) {
        local_file <- tsv_file
      }
    }
    
    if (!is.null(local_file) && file.exists(local_file)) {
      message("Processing local count file: ", local_file)
      
      clean_local_file <- normalizePath(local_file, winslash="/")
      
      dt <- tryCatch({
        fread(local_file)
      }, error = function(e) {
        if (grepl("UTF-16", e$message, ignore.case=TRUE)) {
          message("File is encoded in UTF-16. Recoding to UTF-8 using Python...")
          utf8_file <- gsub("\\.(tsv|csv|txt|tab)(\\.gz)?$", "_utf8.tsv", local_file, ignore.case=TRUE)
          clean_utf8_file <- normalizePath(utf8_file, winslash="/", mustWork=FALSE)
          
          py_file <- file.path(tempdir(), "recode.py")
          cat(paste0(
            "import pandas as pd\n",
            "for enc in ['utf-16', 'utf-16le', 'utf-16be', 'utf-8']:\n",
            "    try:\n",
            "        df = pd.read_csv('", clean_local_file, "', sep=None, engine='python', encoding=enc)\n",
            "        df.to_csv('", clean_utf8_file, "', sep='\\t', index=False, encoding='utf-8')\n",
            "        break\n",
            "    except Exception:\n",
            "        pass\n"
          ), file=py_file)
          system(paste0("python ", py_file))
          if (file.exists(utf8_file)) {
            return(fread(utf8_file))
          }
        }
        stop(e)
      })
      
      if (dataset == "GSE286315") {
        gse286315_cols <- c(
          "gene_id", "S08NP0D", "S08YP0D", "18NP0D_count", "18NP6B1_count",
          "18YP0D_count", "18YP6C2_count", "30N-CON_count", "30N-KO_count",
          "33N-CON_count", "33N-KO_count", "34N-CON_count", "34N-KO_count"
        )
        dt <- dt[, gse286315_cols, with=FALSE]
        mapping_dict <- c(
          "S08NP0D" = "GSM8723553", "S08YP0D" = "GSM8723550",
          "18NP0D_count" = "GSM8723554", "18NP6B1_count" = "GSM8723555",
          "18YP0D_count" = "GSM8723551", "18YP6C2_count" = "GSM8723552",
          "30N-CON_count" = "GSM8723559", "30N-KO_count" = "GSM8723556",
          "33N-CON_count" = "GSM8723560", "33N-KO_count" = "GSM8723557",
          "34N-CON_count" = "GSM8723561", "34N-KO_count" = "GSM8723558"
        )
        setnames(dt, names(mapping_dict), mapping_dict)
      }
      
      # Filter columns to counts if mixed with RPKM/FPKM/Annotations
      col_names <- colnames(dt)
      gene_col <- col_names[1]
      count_cols <- grep("count|read", col_names, ignore.case=TRUE, value=TRUE)
      # Exclude annotation columns
      count_cols <- count_cols[!count_cols %in% c("Chromosome", "Start.Site", "End.Site", "Direction", "Length", "Name", "GO", "KEGG", "EC", "eggNOG.Class", "eggNOG", "Description", "Entrez_geneID", "UniProtAC", "BP", "CC", "MF", "Pathway")]
      
      if (length(count_cols) > 0) {
        dt <- dt[, c(gene_col, count_cols), with=FALSE]
        # Clean suffixes from count column names
        clean_names <- gsub("(?i)\\.(read\\.count|fpkm|tpm|rpkm|counts?|reads?)$", "", count_cols, perl=TRUE)
        clean_names <- gsub("(?i)[_-](read[_-]count|counts?|reads?)$", "", clean_names, perl=TRUE)
        setnames(dt, count_cols, clean_names)
      }

      # Convert to counts matrix
      counts <- as.matrix(dt[, -1, with=FALSE])
      rownames(counts) <- as.character(dt[[1]])
      
      # Fetch metadata from GEO (getGPL=FALSE is safe here since RNA-seq doesn't use platform fData)
      suppressMessages(library(GEOquery))
      local_matrix_file <- file.path("data/expr", paste0(dataset, "_series_matrix.txt.gz"))
      if (file.exists(local_matrix_file)) {
        message("Using locally cached series matrix file: ", local_matrix_file)
        g <- getGEO(filename=local_matrix_file, getGPL=FALSE)
      } else {
        g <- getGEO(dataset, GSEMatrix=TRUE, getGPL=FALSE, destdir="data/expr")
      }
      if (is.list(g)) {
        suppressMessages(library(data.table))
        meta <- as.data.frame(rbindlist(lapply(g, pData), fill=TRUE))
      } else {
        meta <- as.data.frame(pData(g))
      }
      if ("geo_accession" %in% colnames(meta)) {
        rownames(meta) <- meta$geo_accession
      }
      
      # If meta has a lot of samples (e.g. multi-assay) and counts has fewer, try filtering meta to the matching assay
      if (nrow(meta) > ncol(counts)) {
        for (keyword in c("lncRNA", "mRNA", "miRNA", "MeDIP")) {
          if (grepl(keyword, local_file, ignore.case=TRUE)) {
            matching_rows <- grepl(keyword, meta$title, ignore.case=TRUE)
            if (sum(matching_rows) == ncol(counts)) {
              message("Detected multi-assay dataset. Subsetting meta to matching assay keyword: ", keyword)
              meta <- meta[matching_rows, , drop=FALSE]
              break
            }
          }
        }
      }
      
      if (dataset == "GSE87809") {
        mapping_dict <- c(
          "TB7578" = "GSM2341167",
          "TB7579" = "GSM2341168",
          "TB7580" = "GSM2341169",
          "TB9107" = "GSM2341170",
          "TB7575" = "GSM2341171",
          "TB7576" = "GSM2341172",
          "TB7577" = "GSM2341173",
          "TB9110" = "GSM2341174",
          "TB9112" = "GSM2341175"
        )
        colnames(counts) <- mapping_dict[colnames(counts)]
      }

      # Align sample names (column names of count matrix may be titles, descriptions, etc.)
      mapped_names <- NULL
      clean_str <- function(x) {
        x <- tolower(x)
        x <- gsub("[[:punct:][:space:]_]+", "", x)
        return(x)
      }
      
      clean_colnames <- clean_str(colnames(counts))
      
      # If counts columns are generic R-assigned names (V1, V2, etc.), skip candidate matching and align in order
      is_generic_r_names <- all(grepl("^v[0-9]+$", clean_colnames))
      
      if (!is_generic_r_names) {
        for (col_candidate in c("title", "description", "description.1", "description.2", "geo_accession")) {
          if (col_candidate %in% colnames(meta)) {
            mapping <- setNames(rownames(meta), as.character(meta[[col_candidate]]))
            clean_keys <- clean_str(names(mapping))
            
            # Try exact matching first (after basic cleaning)
            match_idx <- match(clean_colnames, clean_keys)
            mapped_candidates <- mapping[match_idx]
            
            # If exact match is not enough, try substring matching
            if (sum(!is.na(mapped_candidates)) < 0.5 * nrow(meta)) {
              mapped_candidates <- rep(NA_character_, length(clean_colnames))
              for (i in seq_along(clean_colnames)) {
                cc <- clean_colnames[i]
                if (is.na(cc) || cc == "" || nchar(cc) < 2) next
                # Check if cc matches key exactly
                idx <- which(clean_keys == cc)
                if (length(idx) == 0) {
                  # Check if cc is substring of any key or vice-versa
                  idx <- which(sapply(clean_keys, function(k) nchar(k) > 0 && (grepl(cc, k, fixed=TRUE) || grepl(k, cc, fixed=TRUE))))
                }
                if (length(idx) == 1) {
                  mapped_candidates[i] <- mapping[idx]
                }
              }
            }
            
            # If we matched a reasonable number of samples (at least half of meta, or at least 2)
            if (sum(!is.na(mapped_candidates)) >= min(2, nrow(meta))) {
              mapped_names <- mapped_candidates
              break
            }
          }
        }
      }
      
      # If mapped_names is still NULL, and the number of columns in counts is exactly equal to nrow(meta), map in order
      if (is.null(mapped_names) && ncol(counts) == nrow(meta)) {
        message("Number of count columns matches nrow(meta) exactly. Aligning in order.")
        mapped_names <- setNames(rownames(meta), colnames(counts))
      }
      
      if (!is.null(mapped_names)) {
        colnames(counts) <- mapped_names
      }
      
      # Subset to matched samples
      common_samples <- intersect(colnames(counts), rownames(meta))
      if (length(common_samples) > 0) {
        counts <- counts[, common_samples, drop=FALSE]
        # Force numeric representation to prevent character matrices
        mode(counts) <- "numeric"
        meta <- meta[common_samples, , drop=FALSE]
      } else {
        stop("Could not align column names of counts with GEO sample accessions.")
      }

      
      # Map Ensembl IDs to symbols if they look like Ensembl IDs
      if (any(grepl("^ENSG", rownames(counts)))) {
        suppressMessages(library(org.Hs.eg.db))
        ensembl_ids <- gsub("\\..*$", "", rownames(counts)) # strip version suffix if present
        gene_map <- tryCatch({
          select(org.Hs.eg.db, keys=ensembl_ids, columns="SYMBOL", keytype="ENSEMBL")
        }, error = function(e) NULL)
        
        if (!is.null(gene_map)) {
          # Match and assign symbols
          symbols <- gene_map$SYMBOL[match(ensembl_ids, gene_map$ENSEMBL)]
          rownames(counts) <- symbols
        }
      }
      
      valid <- !is.na(rownames(counts)) & rownames(counts) != ""
      counts <- counts[valid, , drop=FALSE]
      
      # Aggregate duplicates by summing raw counts
      if (any(duplicated(rownames(counts)))) {
        counts_dt <- data.table(gene = rownames(counts), counts)
        counts_agg <- counts_dt[, lapply(.SD, sum), by = gene]
        counts <- as.matrix(counts_agg[, -1, with=FALSE])
        rownames(counts) <- counts_agg$gene
      }
      
    } else {
      stop("Could not retrieve count matrix or supplementary files for ", dataset)
    }
  }
} else {
  # Microarray -> GEO_matrix
  suppressMessages(library(GEOquery))
  
  # Bypass downloading heavy 1.8 GB GPL SOFT files for GPL13497 and GPL2895
  # by setting getGPL=FALSE. We will map their probes locally via Bioconductor packages.
  # For other platforms, getGPL=TRUE (default) is used to load cached GPL files.
  get_gpl_val <- TRUE
  if (file.exists("config/samplesheet.csv")) {
    samples_dt <- fread("config/samplesheet.csv")
    r_row <- samples_dt[dataset_id == dataset]
    if (nrow(r_row) > 0 && r_row$platform[1] %in% c("GPL13497", "GPL2895")) {
      get_gpl_val <- FALSE
      message("Platform is ", r_row$platform[1], ". Bypassing heavy submitter GPL download (getGPL=FALSE).")
    }
  }
  local_matrix_file <- file.path("data/expr", paste0(dataset, "_series_matrix.txt.gz"))
  if (file.exists(local_matrix_file)) {
    message("Using locally cached series matrix file: ", local_matrix_file)
    g <- getGEO(filename=local_matrix_file, getGPL=get_gpl_val, destdir="data/expr")
  } else {
    g <- getGEO(dataset, GSEMatrix=TRUE, getGPL=get_gpl_val, destdir="data/expr")
  }
  if (is.list(g)) {
    gset <- g[[1]]
  } else {
    gset <- g
  }
  
  counts <- exprs(gset)
  meta <- as.data.frame(pData(gset))
  
  # Map probe IDs to symbols using platform fData
  fdata <- fData(gset)
  symbol_col <- NULL
  possible_cols <- c("Gene Symbol", "gene_symbol", "GENE_SYMBOL", "SYMBOL", "Gene.Symbol", "GENE", "Gene Symbol (RefSeq)")
  for (col in possible_cols) {
    if (col %in% colnames(fdata)) {
      symbol_col <- col
      break
    }
  }
  
  if (is.null(symbol_col)) {
    cols_lower <- tolower(colnames(fdata))
    symbol_idx <- which(cols_lower == "symbol" | cols_lower == "gene symbol" | cols_lower == "gene_symbol")
    if (length(symbol_idx) > 0) {
      symbol_col <- colnames(fdata)[symbol_idx[1]]
    } else {
      symbol_idx_any <- grep("symbol", cols_lower)
      if (length(symbol_idx_any) > 0) {
        symbol_col <- colnames(fdata)[symbol_idx_any[1]]
      } else {
        assign_idx <- which(cols_lower == "gene_assignment" | cols_lower == "gene assignment")
        if (length(assign_idx) > 0) {
          symbol_col <- colnames(fdata)[assign_idx[1]]
        }
      }
    }
  }
  
  symbols <- NULL
  if (!is.null(symbol_col)) {
    symbols <- fdata[[symbol_col]]
    if (grepl("assignment", symbol_col, ignore.case=TRUE)) {
      symbols <- sapply(strsplit(as.character(symbols), " /// "), function(x) {
        if (length(x) == 0 || is.na(x[1]) || x[1] == "---" || x[1] == "") return(NA)
        parts <- strsplit(x[1], " // ")[[1]]
        if (length(parts) >= 2) return(parts[2])
        return(NA)
      })
    } else {
      symbols <- sapply(strsplit(as.character(symbols), " /// "), `[`, 1)
      symbols <- sapply(strsplit(as.character(symbols), "//"), `[`, 1)
      symbols <- trimws(symbols)
    }
  } else {
    # Fallback to Bioconductor annotation packages if symbol_col is not found or empty
    platform <- annotation(gset)
    message("Could not find gene symbol column in fData. Platform detected: ", platform)
    
    db_pkg <- NULL
    if (platform == "GPL13497" || platform == "GPL13497.db" || grepl("026652", platform)) {
      db_pkg <- "HsAgilentDesign026652.db"
    } else if (platform == "GPL2895" || platform == "GPL2895.db" || grepl("hwgcod", platform)) {
      db_pkg <- "hwgcod.db"
    }
    
    if (!is.null(db_pkg)) {
      message("Attempting to map probe IDs to symbols using package: ", db_pkg)
      if (suppressMessages(require(db_pkg, character.only=TRUE))) {
        db_obj <- get(db_pkg)
        probe_ids <- rownames(counts)
        mapped_df <- tryCatch({
          select(db_obj, keys=probe_ids, columns="SYMBOL", keytype="PROBEID")
        }, error = function(e) NULL)
        
        if (!is.null(mapped_df) && nrow(mapped_df) > 0) {
          # Match mapped symbols back to probe_ids
          matched_symbols <- mapped_df$SYMBOL[match(probe_ids, mapped_df$PROBEID)]
          symbols <- matched_symbols
          message("Successfully mapped probes using ", db_pkg)
        }
      }
    }
  }
  
  if (!is.null(symbols)) {
    rownames(counts) <- symbols
    valid <- !is.na(rownames(counts)) & rownames(counts) != "" & rownames(counts) != "---"
    counts <- counts[valid, , drop=FALSE]
    
    # Resolve duplicates by taking mean (log2 intensity)
    if (any(duplicated(rownames(counts)))) {
      counts_dt <- data.table(gene = rownames(counts), counts)
      counts_agg <- counts_dt[, lapply(.SD, mean), by = gene]
      counts <- as.matrix(counts_agg[, -1, with=FALSE])
      rownames(counts) <- counts_agg$gene
    }
  } else {
    warning("Could not map probe IDs to gene symbols. Keeping probe IDs as rownames.")
  }
}

# ----------------- Group Standardization -----------------
# Ensure a standardized 'group' column is created in the metadata
meta$group <- NA

if (dataset == "GSE34526") {
  col_name <- NULL
  if ("disease:ch1" %in% colnames(meta)) {
    col_name <- "disease:ch1"
  } else if ("characteristics_ch1.1" %in% colnames(meta)) {
    col_name <- "characteristics_ch1.1"
  }
  if (!is.null(col_name)) {
    val <- as.character(meta[[col_name]])
    meta$group[grepl("normal", val, ignore.case=TRUE)] <- "Control"
    meta$group[grepl("Polycystic ovary", val, ignore.case=TRUE)] <- "PCOS"
  }
} else if (dataset == "GSE199225") {
  col_name <- NULL
  if ("characteristics_ch1.2" %in% colnames(meta)) {
    col_name <- "characteristics_ch1.2"
  }
  if (!is.null(col_name)) {
    val <- as.character(meta[[col_name]])
    meta$group[grepl("CTRL", val, ignore.case=TRUE)] <- "Control"
    meta$group[grepl("PCOS", val, ignore.case=TRUE)] <- "PCOS"
  }
} else if (dataset == "GSE6364") {
  col_name <- NULL
  if ("characteristics_ch1" %in% colnames(meta)) {
    col_name <- "characteristics_ch1"
  }
  if (!is.null(col_name)) {
    val <- as.character(meta[[col_name]])
    meta$group[grepl("Normal", val, ignore.case=TRUE)] <- "Control"
    meta$group[grepl("Endometriosis", val, ignore.case=TRUE)] <- "Endometriosis"
  }
} else if (dataset == "GSE31683") {
  # Clean up replicate suffixes to get clean treatment groups
  if ("description.1" %in% colnames(meta)) {
    meta$group <- gsub("_[0-9]+$", "", as.character(meta$description.1))
  }
} else if (dataset == "GSE10946") {
  meta$group[grepl("PCOS", meta$title, ignore.case=TRUE)] <- "PCOS"
  meta$group[grepl("nonPCOS", meta$title, ignore.case=TRUE)] <- "Control"
} else if (dataset == "GSE315857") {
  meta$group[grepl("Endometriosis|Eosis", meta$title, ignore.case=TRUE)] <- "Endometriosis"
  meta$group[grepl("Control", meta$title, ignore.case=TRUE)] <- "Control"
} else if (dataset == "GSE6798") {
  meta$group[grepl("PCOS", meta$title, ignore.case=TRUE)] <- "PCOS"
  meta$group[grepl("control", meta$title, ignore.case=TRUE)] <- "Control"
} else if (dataset == "GSE8157") {
  meta$group[grepl("case", meta$title, ignore.case=TRUE)] <- "PCOS"
  meta$group[grepl("control", meta$title, ignore.case=TRUE)] <- "Control"
} else if (dataset == "GSE168404") {
  meta$group[grepl("_P[0-9]+", meta$title, ignore.case=TRUE)] <- "PCOS"
  meta$group[grepl("_C[0-9]+", meta$title, ignore.case=TRUE)] <- "Control"
} else if (dataset == "GSE294074") {
  meta$group[grepl("PCOS", meta$title, ignore.case=TRUE)] <- "PCOS"
  meta$group[grepl("Ctrl", meta$title, ignore.case=TRUE)] <- "Control"
} else if (dataset == "GSE7305") {
  meta$group[grepl("Disease", meta$characteristics_ch1, ignore.case=TRUE)] <- "Endometriosis"
  meta$group[grepl("Normal", meta$characteristics_ch1, ignore.case=TRUE)] <- "Control"
} else if (dataset == "GSE5090") {
  meta$group[grepl("pcos", meta$title, ignore.case=TRUE)] <- "PCOS"
  meta$group[grepl("control", meta$title, ignore.case=TRUE)] <- "Control"
} else if (dataset == "GSE316028") {
  col_idx <- grep("treatment", colnames(meta), ignore.case=TRUE)
  if (length(col_idx) > 0) {
    val <- as.character(meta[[col_idx[1]]])
    meta$group[grepl("NC", val, ignore.case=TRUE)] <- "NC"
    meta$group[grepl("siPC", val, ignore.case=TRUE)] <- "siPC"
  }
} else if (dataset == "GSE44207") {
  meta$group[grepl("untreated", meta$title, ignore.case=TRUE)] <- "untreated"
  meta$group[grepl("VPA-treated", meta$title, ignore.case=TRUE)] <- "VPA-treated"
} else if (dataset == "GSE51981") {
  col_idx <- grep("endometriosis/no endometriosis|endometriosis", colnames(meta), ignore.case=TRUE)
  col_idx <- col_idx[!grepl("severity", colnames(meta)[col_idx], ignore.case=TRUE)]
  if (length(col_idx) > 0) {
    val <- as.character(meta[[col_idx[1]]])
    meta$group[grepl("Non-Endometriosis|no endometriosis", val, ignore.case=TRUE)] <- "Control"
    meta$group[grepl("^Endometriosis$", val, ignore.case=TRUE)] <- "Endometriosis"
  }
} else if (dataset == "GSE286315") {
  meta$group[grepl("Control EcSCs", meta$title, ignore.case=TRUE)] <- "Control EcSCs"
  meta$group[grepl("HOXC4 KO EcSCs", meta$title, ignore.case=TRUE)] <- "HOXC4 KO"
} else if (dataset == "GSE207000") {
  meta$group[grepl("EGR1-KD", meta$title, ignore.case=TRUE)] <- "EGR1-KD"
  meta$group[grepl("Scr", meta$title, ignore.case=TRUE)] <- "Scr"
} else if (dataset == "GSE262302") {
  col_idx <- grep("tissue", colnames(meta), ignore.case=TRUE)
  if (length(col_idx) > 0) {
    val <- as.character(meta[[col_idx[1]]])
    meta$group[grepl("non-adenomyosis|control", val, ignore.case=TRUE)] <- "Control"
    meta$group[grepl("^adenomyosis", val, ignore.case=TRUE)] <- "Adenomyosis"
  }
}

# If group column is still completely NA, try using the samplesheet-specified column and labels
if (all(is.na(meta$group)) && file.exists("config/samplesheet.csv")) {
  samples_dt <- tryCatch({ fread("config/samplesheet.csv") }, error = function(e) NULL)
  if (!is.null(samples_dt)) {
    r_row <- samples_dt[dataset_id == dataset]
    if (nrow(r_row) > 0) {
      group_col <- r_row$group_column[1]
      case_lbl <- r_row$case_label[1]
      ctrl_lbl <- r_row$control_label[1]
      
      if (!is.na(group_col) && group_col != "" && group_col %in% colnames(meta)) {
        val <- as.character(meta[[group_col]])
        meta$group[val == case_lbl] <- case_lbl
        meta$group[val == ctrl_lbl] <- ctrl_lbl
        # Case insensitive check
        na_idx <- is.na(meta$group)
        if (any(na_idx)) {
          meta$group[na_idx & tolower(val) == tolower(case_lbl)] <- case_lbl
          meta$group[na_idx & tolower(val) == tolower(ctrl_lbl)] <- ctrl_lbl
        }
        # Substring/grepl check
        na_idx <- is.na(meta$group)
        if (any(na_idx)) {
          meta$group[na_idx & grepl(case_lbl, val, ignore.case=TRUE)] <- case_lbl
          meta$group[na_idx & grepl(ctrl_lbl, val, ignore.case=TRUE)] <- ctrl_lbl
        }
      }
    }
  }
}

# If group column is still completely NA, try fallback to check for any column name matching 'group'
if (all(is.na(meta$group))) {
  grp_cols <- grep("group", colnames(meta), ignore.case=TRUE, value=TRUE)
  grp_cols <- grp_cols[grp_cols != "group"]
  if (length(grp_cols) > 0) {
    meta$group <- meta[[grp_cols[1]]]
  }
}

message("Sample group breakdown:")
print(table(meta$group, useNA="always"))

# Save standard RDS file
saveRDS(list(counts=counts, meta=meta), out_file)

# Save sidecar TSVs for Python ML script consumption
counts_tsv <- gsub("\\.expr\\.rds$", ".counts.tsv", out_file)
meta_tsv <- gsub("\\.expr\\.rds$", ".meta.tsv", out_file)

counts_dt <- as.data.table(counts)
if (!is.null(rownames(counts))) {
  counts_dt <- cbind(gene = rownames(counts), counts_dt)
} else {
  counts_dt <- cbind(gene = paste0("Gene", 1:nrow(counts)), counts_dt)
}
fwrite(counts_dt, counts_tsv, sep="\t")
meta_df <- as.data.frame(meta)
meta_dt <- as.data.table(meta_df)
meta_dt$sample <- rownames(meta_df)
if ("sample" %in% names(meta_dt)) {
  setcolorder(meta_dt, c("sample", setdiff(names(meta_dt), "sample")))
}
fwrite(meta_dt, meta_tsv, sep="\t")

message("Successfully saved RDS to: ", out_file)
message("Successfully saved sidecar counts to: ", counts_tsv)
message("Successfully saved sidecar meta to: ", meta_tsv)

gc()
 