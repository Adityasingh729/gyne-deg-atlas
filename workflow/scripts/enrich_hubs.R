#!/usr/bin/env Rscript
# clusterProfiler enrichment + STRING PPI hub detection on the meta signature.

suppressMessages({
  library(optparse)
  library(data.table)
})

opt <- parse_args(OptionParser(option_list=list(
  make_option("--meta"),
  make_option("--out_enrich"),
  make_option("--out_hubs")
)))

meta_file <- opt$meta
out_enrich <- opt$out_enrich
out_hubs <- opt$out_hubs

# Ensure output directory exists
dir.create(dirname(out_enrich), showWarnings=FALSE, recursive=TRUE)
dir.create(dirname(out_hubs), showWarnings=FALSE, recursive=TRUE)

message("Running enrichment and hub analysis for file: ", meta_file)

# Read meta signature
meta_dt <- fread(meta_file)
if (nrow(meta_dt) == 0) {
  # Write empty files
  fwrite(data.table(), out_enrich, sep="\t")
  fwrite(data.table(), out_hubs, sep="\t")
  quit(save="no", status=0)
}

# Select significant genes with meta_padj < 0.05 and |meta_log2FC| > 1
sig_df <- meta_dt[meta_padj < 0.05 & abs(meta_log2FC) > 1, ]

# Fallback: if there are too few significant genes, relax thresholds to get a list to analyze
if (nrow(sig_df) < 10) {
  message("Fewer than 10 genes passed the padj < 0.05 & |log2FC| > 1 threshold. Relaxing to top 100 genes by p-value.")
  sig_df <- meta_dt[order(meta_p), ]
  sig_df <- sig_df[!is.na(meta_p), ]
  if (nrow(sig_df) > 100) {
    sig_df <- sig_df[1:100, ]
  }
}

sig_genes <- unique(sig_df$gene)
message("Number of genes selected for enrichment & PPI network: ", length(sig_genes))

if (length(sig_genes) == 0) {
  # Write empty output files if no genes available
  fwrite(data.table(), out_enrich, sep="\t")
  fwrite(data.table(gene=character(), degree=numeric(), betweenness=numeric()), out_hubs, sep="\t")
  quit(save="no", status=0)
}

# 1) clusterProfiler GO and KEGG enrichment
enrich_all <- data.table()

eg_map <- tryCatch({
  suppressMessages({
    library(clusterProfiler)
    library(org.Hs.eg.db)
  })
  bitr(sig_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
}, error = function(e) {
  warning("Mapping gene symbols to Entrez IDs failed: ", e$message)
  NULL
})

if (!is.null(eg_map) && nrow(eg_map) > 0) {
  # GO BP Enrichment
  ego <- tryCatch({
    enrichGO(gene          = eg_map$ENTREZID,
             OrgDb         = org.Hs.eg.db,
             ont           = "BP",
             pAdjustMethod = "BH",
             pvalueCutoff  = 0.05,
             qvalueCutoff  = 0.2,
             readable      = TRUE)
  }, error = function(e) {
    warning("enrichGO failed: ", e$message)
    NULL
  })
  
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    ego_df <- as.data.frame(ego)
    ego_df$source <- "GO:BP"
    enrich_all <- rbindlist(list(enrich_all, as.data.table(ego_df)), fill=TRUE)
  }
  
  # KEGG Pathway Enrichment
  ekegg <- tryCatch({
    enrichKEGG(gene         = eg_map$ENTREZID,
               organism     = "hsa",
               pAdjustMethod = "BH",
               pvalueCutoff  = 0.05,
               qvalueCutoff  = 0.2)
  }, error = function(e) {
    warning("enrichKEGG failed: ", e$message)
    NULL
  })
  
  if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
    ekegg_df <- as.data.frame(ekegg)
    # Check if we can convert KEGG Entrez IDs to gene symbols for compatibility
    ekegg_df <- as.data.table(ekegg_df)
    # We can try to map ENTREZ IDs in geneID back to SYMBOL for KEGG
    # For now, write it as is with source tagged
    ekegg_df[, source := "KEGG"]
    enrich_all <- rbindlist(list(enrich_all, ekegg_df), fill=TRUE)
  }
}

# If empty, create template schema
if (nrow(enrich_all) == 0) {
  enrich_all <- data.table(
    ID = character(),
    Description = character(),
    GeneRatio = character(),
    BgRatio = character(),
    pvalue = numeric(),
    p.adjust = numeric(),
    qvalue = numeric(),
    geneID = character(),
    Count = integer(),
    source = character()
  )
}
fwrite(enrich_all, out_enrich, sep="\t")
message("Saved enrichment results to: ", out_enrich)

# 2) STRING PPI Network and Hub detection
hubs_dt <- tryCatch({
  suppressMessages({
    library(STRINGdb)
    library(igraph)
  })
  
  # Initialize STRINGdb
  # version 11.5, human (9606), confidence threshold 400
  string_db <- STRINGdb$new(version="11.5", species=9606, score_threshold=400, input_directory="")
  
  # Map genes
  mapped <- string_db$map(data.frame(gene = sig_genes), "gene", removeUnmappedRows=TRUE)
  
  if (!is.null(mapped) && nrow(mapped) > 0) {
    # Fetch protein interactions
    ppi <- string_db$get_interactions(mapped$STRING_id)
    
    if (!is.null(ppi) && nrow(ppi) > 0) {
      # Create undirected igraph
      g <- graph_from_data_frame(ppi[, c("from", "to")], directed=FALSE)
      g <- simplify(g)
      
      # Degree and Betweenness centralities
      deg <- degree(g)
      bet <- betweenness(g)
      
      hubs <- data.table(
        STRING_id = names(deg),
        degree = as.numeric(deg),
        betweenness = as.numeric(bet)
      )
      
      # Map STRING IDs back to HGNC gene symbols
      mapped_dt <- as.data.table(mapped)
      hubs <- merge(hubs, mapped_dt[, .(gene, STRING_id)], by="STRING_id")
      hubs <- hubs[order(-degree)]
      
      # Standard columns
      hubs <- hubs[, .(gene, degree, betweenness)]
      hubs
    } else {
      warning("No STRING PPI interactions found for these genes.")
      data.table(gene = sig_genes, degree = 0, betweenness = 0)
    }
  } else {
    warning("No genes could be mapped in STRINGdb.")
    data.table(gene = sig_genes, degree = 0, betweenness = 0)
  }
}, error = function(e) {
  warning("STRINGdb analysis failed or skipped (e.g. offline/network issue): ", e$message)
  data.table(gene = sig_genes, degree = 0, betweenness = 0)
})

fwrite(hubs_dt, out_hubs, sep="\t")
message("Saved PPI hubs to: ", out_hubs)

gc()
 