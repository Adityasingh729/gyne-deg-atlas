# =====================================================================
# inspect_dataset.R  -  GEO dataset fact extractor (for teammates)
# =====================================================================
# WHAT THIS DOES:
#   Reads a GEO series matrix file + a raw count file and PRINTS the
#   facts you need to answer the dataset self-check questions.
#
# WHAT THIS DOES *NOT* DO:
#   It does NOT decide case vs control for you.
#   It does NOT decide include vs exclude.
#   It does NOT check if your grouping was correct.
#   Those are YOUR decisions - use the Decision Document.
#
# ---------------------------------------------------------------------
# IMPORTANT - PUT BOTH FILES IN THE SAME FOLDER, then run with TWO paths:
#
#   source("inspect_dataset.R")
#   inspect_dataset(
#     series_matrix = "PATH/TO/GSExxxxx_series_matrix.txt.gz",
#     count_file    = "PATH/TO/your_raw_counts.xlsx"
#   )
#
# Series matrix: .txt or .txt.gz
# Count file:    .xlsx, .csv, .tsv, or .txt
# ---------------------------------------------------------------------
# LIMITATION - READ THIS:
#   This script is RELIABLE only for WELL-LABELLED GEO datasets.
#   Some datasets hide disease/tissue info in free text, or use odd
#   field names. For those, the output will be INCOMPLETE.
#   >>> ALWAYS cross-check this output against the live GEO web page. <<<
#   This script was tested only on ONE dataset (GSE190580). Its
#   behaviour on other datasets is NOT guaranteed. Verify every time.
# =====================================================================

inspect_dataset <- function(series_matrix, count_file) {

  line <- function(ch="-") cat(paste0(paste(rep(ch,70),collapse=""),"\n"))
  hdr  <- function(t){ line("="); cat(t,"\n"); line("=") }

  # ---------- checks ----------
  if (missing(series_matrix) || missing(count_file))
    stop("You must give BOTH file paths: series_matrix= and count_file=")
  if (!file.exists(series_matrix)) stop("Series matrix file not found: ", series_matrix)
  if (!file.exists(count_file))    stop("Count file not found: ", count_file)

  hdr("GEO DATASET INSPECTION REPORT")
  cat("Series matrix file :", basename(series_matrix), "\n")
  cat("Count file         :", basename(count_file), "\n")
  cat("\n>>> REMINDER: cross-check everything below against the GEO web page. <<<\n\n")

  # =================================================================
  # PART 1 - SERIES MATRIX (the grouping evidence)
  # =================================================================
  con <- if (grepl("\\.gz$", series_matrix, ignore.case=TRUE)) gzfile(series_matrix,"rt") else file(series_matrix,"rt")
  sm <- readLines(con); close(con)

  get1 <- function(tag){
    r <- sm[startsWith(sm, tag)]
    if (length(r)==0) return(NA_character_)
    v <- sub(paste0("^",tag,"\\t?"), "", r[1])
    gsub('"', '', v)
  }
  # multi-value sample rows -> character vector (one per sample)
  getrow <- function(tag){
    r <- sm[startsWith(sm, tag)]
    if (length(r)==0) return(list())
    lapply(r, function(x){
      v <- sub(paste0("^",tag,"\\t"), "", x)
      gsub('"', '', strsplit(v, "\t")[[1]])
    })
  }

  hdr("PART 1: STUDY-LEVEL INFO (from series matrix)")
  cat("TITLE:\n  ", get1("!Series_title"), "\n\n")
  cat("SUMMARY:\n  ", get1("!Series_summary"), "\n\n")
  cat("OVERALL DESIGN:\n  ", get1("!Series_overall_design"), "\n")

  titles  <- getrow("!Sample_title")
  geoacc  <- getrow("!Sample_geo_accession")
  source_ <- getrow("!Sample_source_name_ch1")
  chars   <- getrow("!Sample_characteristics_ch1")   # may be several rows

  nsamp <- if (length(titles)>0) length(titles[[1]]) else
           if (length(geoacc)>0) length(geoacc[[1]]) else 0

  hdr(paste0("PART 2: PER-SAMPLE TABLE  (", nsamp, " samples in series matrix)"))
  if (nsamp==0){
    cat("Could not read sample rows. This series matrix may be non-standard.\n")
    cat(">>> Inspect the GEO page manually. <<<\n")
  } else {
    st <- if (length(titles)>0) titles[[1]] else rep("",nsamp)
    sg <- if (length(geoacc)>0) geoacc[[1]] else rep("",nsamp)
    ss <- if (length(source_)>0) source_[[1]] else rep("",nsamp)
    for (i in seq_len(nsamp)) {
      cat(sprintf("[%02d] %-14s %-12s source: %s\n", i, st[i], sg[i], ss[i]))
      if (length(chars)>0) for (cr in chars) {
        val <- if (length(cr)>=i) cr[i] else ""
        if (nzchar(val)) cat("       -", val, "\n")
      }
    }
    # distinct-value summary of each characteristics row (helps spot arms)
    hdr("PART 3: WHAT VARIES ACROSS SAMPLES  (use this to spot multi-arm designs)")
    cat("source_name unique values:\n")
    if (length(source_)>0) for (u in sort(unique(source_[[1]]))) cat("   -", u, "  (n=", sum(source_[[1]]==u), ")\n")
    if (length(chars)>0) {
      # each chars row: figure out the field key from "key: value"
      for (cr in chars) {
        key <- sub(":.*$","", cr[1])
        cat("characteristic '", key, "' unique values:\n", sep="")
        for (u in sort(unique(cr))) cat("   -", u, "  (n=", sum(cr==u), ")\n")
      }
    }
    cat("\n>>> If more than one tissue / treatment / timepoint appears above,\n")
    cat(">>> this is a MULTI-ARM study. Pick ONE clean comparison. Do NOT pool.\n")
  }

  # =================================================================
  # PART 4 - COUNT FILE (id type + sample count + cross-check)
  # =================================================================
  hdr("PART 4: RAW COUNT FILE")
  ext <- tolower(tools::file_ext(count_file))
  read_counts_head <- function(){
    if (ext=="xlsx"){
      if (!requireNamespace("readxl", quietly=TRUE)) {
        cat("This is an .xlsx file. The 'readxl' package is required but not installed.\n")
        cat("Install it once with:  install.packages(\"readxl\")\n")
        cat("Then run this script again.\n")
        return(NULL)
      }
      as.data.frame(readxl::read_excel(count_file, n_max=5))
    } else {
      sep <- if (ext %in% c("csv")) "," else "\t"
      utils::read.delim(count_file, sep=sep, nrows=5, check.names=FALSE)
    }
  }
  # full column count needs a full header read
  read_ncol <- function(){
    if (ext=="xlsx"){
      if (!requireNamespace("readxl", quietly=TRUE)) return(NA_integer_)
      ncol(as.data.frame(readxl::read_excel(count_file, n_max=1)))
    } else {
      sep <- if (ext=="csv") "," else "\t"
      h <- readLines(count_file, n=1)
      length(strsplit(h, sep)[[1]])
    }
  }

  hd <- tryCatch(read_counts_head(), error=function(e){cat("Could not read count file:", conditionMessage(e), "\n"); NULL})
  if (!is.null(hd)) {
    nc <- read_ncol()
    first_ids <- as.character(hd[[1]])
    # detect ID type from first column
    idtype <- "UNKNOWN"
    if (all(grepl("^ENSG[0-9]+", first_ids)))            idtype <- "Ensembl (ENSG...)"
    else if (all(grepl("^[0-9]+$", first_ids)))          idtype <- "Entrez (numbers)"
    else if (all(grepl("^[A-Z][A-Z0-9-]+$", first_ids))) idtype <- "Gene SYMBOL (e.g. TP53)"
    else if (any(grepl("_at$|^A_[0-9]|ILMN|^[0-9]+_", first_ids))) idtype <- "Probe ID (microarray)"
    else idtype <- "MIXED or unrecognised - inspect manually"

    n_sample_cols <- if (!is.na(nc)) nc-1 else NA   # minus the gene-id column
    cat("First column name :", names(hd)[1], "\n")
    cat("First gene IDs    :", paste(head(first_ids,3),collapse=", "), "\n")
    cat("DETECTED ID TYPE  :", idtype, "\n")
    cat("Sample columns    :", n_sample_cols, "\n\n")

    # cross-check against series matrix
    cat("CROSS-CHECK (samples):\n")
    cat("  series matrix samples :", nsamp, "\n")
    cat("  count file samples    :", n_sample_cols, "\n")
    if (!is.na(n_sample_cols) && nsamp>0) {
      if (n_sample_cols==nsamp) cat("  -> MATCH.\n")
      else cat("  -> MISMATCH. Samples were added or dropped. INVESTIGATE before trusting this file.\n")
    }
  }

  hdr("WHICH ANSWERS COME FROM WHERE")
  cat("FROM THIS SCRIPT : sample list, disease/tissue/treatment labels, how many\n")
  cat("                   arms exist, gene ID type, sample count, cross-check.\n")
  cat("FROM *YOU* (GEO page + Decision Document):\n")
  cat("                   which samples are case vs control, control type,\n")
  cat("                   is it a treatment/knockdown/culture study, include or\n")
  cat("                   exclude, which arm to use. THE SCRIPT DOES NOT DECIDE.\n")
  line("=")
  cat(">>> Cross-check against the GEO web page. This script can be incomplete. <<<\n")
  invisible(NULL)
}
