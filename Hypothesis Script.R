################################################################################
# SCRIPT 4 — Hypothesis Testing + Full GO Search + Novel Gene Discovery
# FIXES:
#   FIX 4: Universe GO map built ONCE outside run_full_go (was rebuilt 4x = ~8 min)
#          Now built once before function call = ~2 min total
#   Correct phyper K per GO term from universe (not circular)
#   K>=10 low count filter for GO terms
#   All poster numbers dynamic — reads kegg_summary.rds and deg_summary.rds
#   No hardcoded numbers anywhere
################################################################################

for (p in c("dplyr","tibble")) {
  if (!requireNamespace(p, quietly=TRUE)) install.packages(p)
  library(p, character.only=TRUE)
}
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
for (p in c("AnnotationDbi","org.At.tair.db","GO.db","DESeq2")) {
  if (!requireNamespace(p, quietly=TRUE)) BiocManager::install(p)
}
library(AnnotationDbi); library(org.At.tair.db); library(GO.db); library(DESeq2)
select <- dplyr::select

# ==============================================================================
# PATHS
# ==============================================================================
base_dir    <- "~/Downloads/deseq2_arabidopsis"
drought_dir <- file.path(base_dir, "arabidopsis_drought")
salt_dir    <- file.path(base_dir, "arabidopsis_salt")
out_dir     <- file.path(base_dir, "Hypothesis_Testing")
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

# ==============================================================================
# LOAD DESeq2 RESULTS
# ==============================================================================
safe_convert <- function(obj) {
  if (inherits(obj,"DESeqResults") || inherits(obj,"DataFrame")) {
    df <- as.data.frame(obj); df$gene <- rownames(df); rownames(df) <- NULL
  } else {
    df <- as.data.frame(obj)
    if (!"gene" %in% names(df)) { df$gene <- rownames(df); rownames(df) <- NULL }
  }
  names(df) <- gsub("log2foldchange","log2FoldChange",names(df),ignore.case=TRUE)
  df <- df[, c("gene",setdiff(names(df),"gene")), drop=FALSE]
  df <- df[!is.na(df$padj) & !is.na(df$log2FoldChange), ]
  df
}

drought_df <- safe_convert(readRDS(file.path(drought_dir, "deseq2_results.rds")))
salt_df    <- safe_convert(readRDS(file.path(salt_dir,    "deseq2_results.rds")))

padj_cut <- 0.05
lfc_cut  <- 1

drought_sig  <- drought_df[drought_df$padj < padj_cut &
                             abs(drought_df$log2FoldChange) > lfc_cut, ]
salt_sig     <- salt_df[   salt_df$padj    < padj_cut &
                             abs(salt_df$log2FoldChange)    > lfc_cut, ]
drought_up   <- drought_sig[drought_sig$log2FoldChange > 0, ]
drought_down <- drought_sig[drought_sig$log2FoldChange < 0, ]
salt_up      <- salt_sig[   salt_sig$log2FoldChange    > 0, ]
salt_down    <- salt_sig[   salt_sig$log2FoldChange    < 0, ]

cat(sprintf("Drought DEGs: %d (%d up / %d down)\n",
            nrow(drought_sig), nrow(drought_up), nrow(drought_down)))
cat(sprintf("Salt DEGs:    %d (%d up / %d down)\n\n",
            nrow(salt_sig), nrow(salt_up), nrow(salt_down)))

universe <- unique(c(drought_df$gene, salt_df$gene))
N        <- length(universe)
cat(sprintf("Universe: %d genes\n\n", N))

##############################################################################
# FIX 4: BUILD UNIVERSE GO MAP ONCE — outside run_full_go
#
# Previously the universe AnnotationDbi::select query ran INSIDE run_full_go
# which was called 4 times (drought_up, drought_down, salt_up, salt_down).
# Each call queried the full universe (~15k genes) against org.At.tair.db —
# approximately 2 minutes per call = ~8 minutes wasted on identical work.
#
# Now built once here, shared across all 4 calls via global objects.
##############################################################################
cat("Building universe GO map (runs once — ~2 minutes)...\n")

universe_go_raw <- tryCatch(
  AnnotationDbi::select(org.At.tair.db,
                        keys    = universe,
                        columns = c("GO","ONTOLOGY"),
                        keytype = "TAIR"),
  error=function(e) { cat(sprintf("ERROR: %s\n", e$message)); NULL }
)

if (is.null(universe_go_raw)) {
  stop("Cannot proceed: universe GO annotation lookup failed.")
}

universe_go_bp      <- universe_go_raw[
  !is.na(universe_go_raw$GO) & universe_go_raw$ONTOLOGY == "BP", ]
k_total_map_global  <- table(universe_go_bp$GO)

universe_go_terms <- tryCatch(
  AnnotationDbi::select(GO.db,
                        keys    = unique(universe_go_bp$GO),
                        columns = "TERM",
                        keytype = "GOID"),
  error=function(e) { cat(sprintf("ERROR: %s\n", e$message)); NULL }
)

if (is.null(universe_go_terms)) {
  stop("Cannot proceed: GO.db term lookup failed.")
}
names(universe_go_terms) <- c("GO","TERM")

universe_go_ann <- merge(universe_go_bp, universe_go_terms, by="GO", all.x=TRUE)

# Precompute per-GO-term gene lists and term name lists once
# These are reused across all 4 run_full_go calls
go_list_global   <- split(universe_go_ann$TAIR, universe_go_ann$GO)
term_list_global <- split(universe_go_ann$TERM,  universe_go_ann$GO)

cat(sprintf("Universe: %d unique BP terms across %d genes\n\n",
            length(k_total_map_global), length(unique(universe))))

##############################################################################
# FULL GO ENRICHMENT — CORRECTED HYPERGEOMETRIC TEST
# run_full_go no longer rebuilds the universe map — uses global objects above
##############################################################################
run_full_go <- function(gene_list, label) {
  cat(sprintf("  Running full GO for %s (%d genes)...\n",
              label, length(gene_list)))
  
  N_univ <- length(unique(universe))
  n_deg  <- length(gene_list)
  
  result_rows <- list()
  for (goid in names(go_list_global)) {
    
    # K_total specific to this GO term from the universe
    k_total <- as.integer(k_total_map_global[goid])
    
    # Low count filter: skip terms with fewer than 10 universe genes
    if (is.na(k_total) || k_total < 10) next
    
    go_genes <- unique(go_list_global[[goid]][!is.na(go_list_global[[goid]])])
    overlap  <- sum(gene_list %in% go_genes)
    if (overlap == 0) next
    
    term <- term_list_global[[goid]][1]
    # Skip uninformative root term
    if (!is.na(term) && term == "biological_process") next
    
    # Corrected hypergeometric test
    # phyper(q, m, n, k):
    #   q = overlap - 1
    #   m = k_total  (universe genes WITH this specific term)
    #   n = N_univ - k_total  (universe genes WITHOUT this term)
    #   k = n_deg  (DEG list size)
    pval <- phyper(
      overlap - 1,
      k_total,
      N_univ - k_total,
      n_deg,
      lower.tail = FALSE
    )
    
    result_rows[[length(result_rows)+1]] <- data.frame(
      GO_ID      = goid,
      TERM       = ifelse(is.na(term), "unknown", term),
      Gene_Count = overlap,
      K_total    = k_total,
      GeneRatio  = round(overlap / n_deg, 4),
      pvalue     = pval,
      stringsAsFactors = FALSE
    )
  }
  
  if (length(result_rows)==0) return(NULL)
  
  res         <- do.call(rbind, result_rows)
  res$padj_BH <- p.adjust(res$pvalue, method="BH")
  res         <- res[order(res$padj_BH, res$pvalue), ]
  res$rank    <- seq_len(nrow(res))
  res$label   <- label
  
  cat(sprintf("    Found %d GO terms | %d significant (padj<0.05)\n",
              nrow(res), sum(res$padj_BH < 0.05)))
  res
}

cat("Running full GO enrichment (no top-20 cutoff)...\n\n")

go_full_drought_up   <- run_full_go(drought_up$gene,   "drought_UP")
go_full_drought_down <- run_full_go(drought_down$gene, "drought_DOWN")
go_full_salt_up      <- run_full_go(salt_up$gene,      "salt_UP")
go_full_salt_down    <- run_full_go(salt_down$gene,    "salt_DOWN")

go_full_all <- dplyr::bind_rows(go_full_drought_up,  go_full_drought_down,
                                go_full_salt_up,     go_full_salt_down)

##############################################################################
# SEARCH FOR ABA / AQUAPORIN / ANTIPORTER
##############################################################################
aba_kw  <- c("abscisic acid","ABA","osmotic stress","water deprivation",
             "response to water","stomata","guard cell","desiccation",
             "abscisic acid-activated","snrk","pyr","pyl","pp2c",
             "response to abscisic","abscisic acid signaling")
aqua_kw <- c("aquaporin","water channel","water transport",
             "major intrinsic protein","MIP","PIP","TIP",
             "water homeostasis","transmembrane water")
anti_kw <- c("sodium","antiporter","NHX","SOS","salt tolerance",
             "cation transport","ion homeostasis","proton antiport",
             "sodium ion transport","CHX","monovalent cation","sodium:proton")

search_go_full <- function(df, keywords, family_name) {
  if (is.null(df) || nrow(df)==0) return(data.frame())
  pattern <- paste(keywords, collapse="|")
  hits    <- df[grepl(pattern, df$TERM, ignore.case=TRUE), ]
  if (nrow(hits)>0) {
    hits$family <- family_name
    cat(sprintf("  %-15s — %d terms | best rank: %d | best padj: %.2e\n",
                family_name, nrow(hits), min(hits$rank), min(hits$padj_BH)))
  } else {
    cat(sprintf("  %-15s — NOT FOUND\n", family_name))
  }
  hits
}

cat("\n══════════════════════════════════════════════\n")
cat("  SEARCHING FULL GO RESULTS (all ranks)\n")
cat("══════════════════════════════════════════════\n\n")
cat("  ABA terms:\n")
aba_go_full  <- search_go_full(go_full_all, aba_kw,  "ABA")
cat("\n  Aquaporin terms:\n")
aqua_go_full <- search_go_full(go_full_all, aqua_kw, "Aquaporin")
cat("\n  Antiporter terms:\n")
anti_go_full <- search_go_full(go_full_all, anti_kw, "Antiporter")

print_family_hits <- function(hits, family_name) {
  if (nrow(hits)==0) {
    cat(sprintf("  %s: NO TERMS FOUND\n", family_name)); return()
  }
  cat(sprintf("\n  %s:\n", family_name))
  for (i in seq_len(nrow(hits))) {
    sig <- if (hits$padj_BH[i]<0.05)  "*** SIGNIFICANT" else
      if (hits$padj_BH[i]<0.20)  "  ~ marginal"   else
        "    (not sig)"
    cat(sprintf("    Rank %4d  padj=%.2e  K=%4d  [%-12s]  %s  %s\n",
                hits$rank[i], hits$padj_BH[i], hits$K_total[i],
                hits$label[i], hits$TERM[i], sig))
  }
}

cat("\n══════════════════════════════════════════════\n")
cat("  DETAILED RESULTS\n")
cat("══════════════════════════════════════════════\n")
print_family_hits(aba_go_full,  "ABA SIGNALING")
print_family_hits(aqua_go_full, "AQUAPORINS")
print_family_hits(anti_go_full, "SODIUM ANTIPORTERS")

##############################################################################
# SAVE FULL GO RESULTS
##############################################################################
save_csv <- function(df, name) {
  if (is.data.frame(df) && nrow(df)>0) {
    write.csv(df, file.path(out_dir,name), row.names=FALSE)
    cat(sprintf("  SAVED (%d rows): %s\n", nrow(df), name))
  } else {
    cat(sprintf("  SKIPPED (empty): %s\n", name))
  }
}

cat("\nSaving full GO results:\n")
save_csv(go_full_drought_up,   "FULL_GO_drought_upregulated.csv")
save_csv(go_full_drought_down, "FULL_GO_drought_downregulated.csv")
save_csv(go_full_salt_up,      "FULL_GO_salt_upregulated.csv")
save_csv(go_full_salt_down,    "FULL_GO_salt_downregulated.csv")
save_csv(aba_go_full,          "SEARCH_ABA_all_ranks.csv")
save_csv(aqua_go_full,         "SEARCH_Aquaporin_all_ranks.csv")
save_csv(anti_go_full,         "SEARCH_Antiporter_all_ranks.csv")

##############################################################################
# HYPOTHESIS VERDICTS
##############################################################################
cat("\n")
cat("████████████████████████████████████████████████\n")
cat("█  WHAT YOUR DATA ACTUALLY SAYS               █\n")
cat("████████████████████████████████████████████████\n\n")

aba_sig  <- if (nrow(aba_go_full)>0)  sum(aba_go_full$padj_BH  < 0.05) else 0L
aqua_sig <- if (nrow(aqua_go_full)>0) sum(aqua_go_full$padj_BH < 0.05) else 0L
anti_sig <- if (nrow(anti_go_full)>0) sum(anti_go_full$padj_BH < 0.05) else 0L

cat("  HYPOTHESIS VERDICTS:\n\n")

cat("  H1 — ABA signaling UP in both drought AND salt:\n")
if (aba_sig > 0) {
  best_aba      <- aba_go_full[aba_go_full$padj_BH < 0.05, ]
  n_drought_aba <- sum(grepl("drought", best_aba$label))
  n_salt_aba    <- sum(grepl("salt",    best_aba$label))
  if (n_drought_aba > 0 && n_salt_aba > 0) {
    cat(sprintf("  SUPPORTED — %d sig ABA GO terms (%d drought, %d salt)\n",
                aba_sig, n_drought_aba, n_salt_aba))
  } else if (n_drought_aba > 0) {
    cat(sprintf("  PARTIAL — ABA significant in DROUGHT only (%d terms)\n",
                n_drought_aba))
  } else {
    cat(sprintf("  PARTIAL — ABA significant in SALT only (%d terms)\n",
                n_salt_aba))
  }
} else if (nrow(aba_go_full) > 0) {
  cat(sprintf("  NOT SUPPORTED — ABA terms exist but none padj<0.05\n"))
  cat(sprintf("     Best: rank %d, padj=%.3f\n",
              min(aba_go_full$rank), min(aba_go_full$padj_BH)))
} else {
  cat("  NOT SUPPORTED — no ABA terms found at any rank\n")
}

cat("\n  H3 — Aquaporins DOWN in drought:\n")
if (aqua_sig > 0) {
  best_aqua <- aqua_go_full[aqua_go_full$padj_BH < 0.05, ]
  n_down    <- sum(grepl("DOWN",    best_aqua$label))
  n_drought <- sum(grepl("drought", best_aqua$label))
  if (n_down > 0 && n_drought > 0) {
    cat(sprintf("  SUPPORTED — %d sig aquaporin terms (drought DOWN)\n", aqua_sig))
  } else {
    cat(sprintf("  PARTIAL — %d sig aquaporin terms but check direction/condition:\n",
                aqua_sig))
    for (i in seq_len(nrow(best_aqua)))
      cat(sprintf("     [%s] %s\n", best_aqua$label[i], best_aqua$TERM[i]))
  }
} else if (nrow(aqua_go_full) > 0) {
  cat(sprintf("  NOT SUPPORTED — terms exist but none padj<0.05\n"))
  cat(sprintf("     Best: rank %d, padj=%.3f\n",
              min(aqua_go_full$rank), min(aqua_go_full$padj_BH)))
} else {
  cat("  NOT SUPPORTED — no aquaporin terms found\n")
}

cat("\n  H2 — Sodium antiporters UP in salt:\n")
if (anti_sig > 0) {
  best_anti <- anti_go_full[anti_go_full$padj_BH < 0.05, ]
  n_salt    <- sum(grepl("salt", best_anti$label))
  if (n_salt > 0) {
    cat(sprintf("  SUPPORTED — %d sig antiporter terms in salt\n", n_salt))
  } else {
    cat(sprintf("  PARTIAL — %d sig antiporter terms but not in salt condition\n",
                anti_sig))
  }
} else if (nrow(anti_go_full) > 0) {
  cat(sprintf("  NOT SUPPORTED — terms exist but none padj<0.05\n"))
  cat(sprintf("     Best: rank %d, padj=%.3f\n",
              min(anti_go_full$rank), min(anti_go_full$padj_BH)))
} else {
  cat("  NOT SUPPORTED — no antiporter terms found\n")
}

cat("\n  TOP SIGNIFICANT GO TERMS (K>=10, corrected):\n\n")

for (res_obj in list(go_full_drought_up,  go_full_drought_down,
                     go_full_salt_up,     go_full_salt_down)) {
  if (is.null(res_obj) || nrow(res_obj)==0) next
  lbl  <- res_obj$label[1]
  sig5 <- head(res_obj[res_obj$padj_BH < 0.05, ], 5)
  if (nrow(sig5)==0) { cat(sprintf("  %s: no significant terms\n", lbl)); next }
  cat(sprintf("  %s top 5:\n", lbl))
  for (i in seq_len(nrow(sig5)))
    cat(sprintf("    rank%3d  padj=%.2e  K=%4d  %s\n",
                sig5$rank[i], sig5$padj_BH[i], sig5$K_total[i], sig5$TERM[i]))
  cat("\n")
}

##############################################################################
# NOVEL GENE DISCOVERY
##############################################################################
cat("████████████████████████████████████████████████\n")
cat("█  NOVEL GENE DEEP DIVE                        █\n")
cat("████████████████████████████████████████████████\n\n")

get_sym <- function(ids) {
  if (length(ids)==0) return(setNames(character(0),character(0)))
  tryCatch({
    s <- AnnotationDbi::mapIds(org.At.tair.db, keys=ids,
                               column="SYMBOL", keytype="TAIR", multiVals="first")
    s[is.na(s)] <- ids[is.na(s)]; s
  }, error=function(e) setNames(ids, ids))
}

get_desc <- function(ids) {
  if (length(ids)==0) return(setNames(character(0),character(0)))
  tryCatch({
    d <- AnnotationDbi::mapIds(org.At.tair.db, keys=ids,
                               column="GENENAME", keytype="TAIR", multiVals="first")
    d[is.na(d)] <- "no description"; d
  }, error=function(e) setNames(rep("lookup failed",length(ids)), ids))
}

get_go_for_gene <- function(gene_id) {
  tryCatch({
    res <- AnnotationDbi::select(org.At.tair.db,
                                 keys    = gene_id,
                                 columns = c("GO","ONTOLOGY"),
                                 keytype = "TAIR")
    res <- res[!is.na(res$GO) & res$ONTOLOGY=="BP", ]
    if (nrow(res)==0) return("NO_BP_GO_ANNOTATION")
    go_names <- AnnotationDbi::select(GO.db,
                                      keys    = unique(res$GO),
                                      columns = "TERM",
                                      keytype = "GOID")
    terms <- go_names$TERM[!is.na(go_names$TERM) &
                             go_names$TERM != "biological_process"]
    if (length(terms)==0) return("NO_SPECIFIC_BP_ANNOTATION")
    paste(unique(terms), collapse=" | ")
  }, error=function(e) "DATABASE_LOOKUP_FAILED")
}

# Novel = gene still identified only by TAIR locus ID (no common name assigned)
is_novel <- function(sym, gene) is.na(sym) | sym==gene | grepl("^AT[0-9]G[0-9]",sym)

drought_sig$symbol <- get_sym(drought_sig$gene)
drought_sig$desc   <- get_desc(drought_sig$gene)
salt_sig$symbol    <- get_sym(salt_sig$gene)
salt_sig$desc      <- get_desc(salt_sig$gene)

novel_drought <- drought_sig[is_novel(drought_sig$symbol, drought_sig$gene), ]
novel_salt    <- salt_sig[   is_novel(salt_sig$symbol,    salt_sig$gene),    ]
both_ids      <- intersect(novel_drought$gene, novel_salt$gene)

cat(sprintf("Uncharacterized genes significant in BOTH stresses: %d\n\n",
            length(both_ids)))

novel_full <- data.frame(gene=both_ids, symbol=get_sym(both_ids),
                         desc=get_desc(both_ids), stringsAsFactors=FALSE)

d_sub <- drought_sig[drought_sig$gene %in% both_ids,
                     c("gene","log2FoldChange","padj")]
names(d_sub)[2:3] <- c("lfc_drought","padj_drought")
s_sub <- salt_sig[salt_sig$gene %in% both_ids,
                  c("gene","log2FoldChange","padj")]
names(s_sub)[2:3] <- c("lfc_salt","padj_salt")

novel_full <- merge(novel_full, d_sub, by="gene", all.x=TRUE)
novel_full <- merge(novel_full, s_sub, by="gene", all.x=TRUE)

novel_full$drought_dir <- ifelse(novel_full$lfc_drought > 0, "UP", "DOWN")
novel_full$salt_dir    <- ifelse(novel_full$lfc_salt    > 0, "UP", "DOWN")
novel_full$pattern     <- paste0("Drought:", novel_full$drought_dir,
                                 " | Salt:", novel_full$salt_dir)

novel_full$high_conf <- abs(novel_full$lfc_drought) > lfc_cut &
  abs(novel_full$lfc_salt)    > lfc_cut &
  novel_full$padj_drought     < 0.01    &
  novel_full$padj_salt        < 0.01

novel_full$combined_lfc <- abs(novel_full$lfc_drought) + abs(novel_full$lfc_salt)
novel_full <- novel_full[order(!novel_full$high_conf, -novel_full$combined_lfc), ]

hc <- novel_full[novel_full$high_conf, ]
cat(sprintf("High-confidence novel candidates (padj<0.01, |lfc|>1 both): %d\n\n",
            nrow(hc)))

top_novel <- head(hc[order(-hc$combined_lfc), ], 10)

cat("Looking up GO annotations for top novel candidates...\n\n")
cat("══════════════════════════════════════════════\n")
cat("  TOP NOVEL GENES — FULL PROFILE\n")
cat("══════════════════════════════════════════════\n\n")

for (i in seq_len(nrow(top_novel))) {
  r      <- top_novel[i, ]
  go_ann <- get_go_for_gene(r$gene)
  cat(sprintf("  %d. %s (%s)\n",    i, r$symbol, r$gene))
  cat(sprintf("     Pattern:  %s\n", r$pattern))
  cat(sprintf("     Drought:  log2FC=%+.2f  padj=%.2e\n", r$lfc_drought, r$padj_drought))
  cat(sprintf("     Salt:     log2FC=%+.2f  padj=%.2e\n", r$lfc_salt,    r$padj_salt))
  cat(sprintf("     Desc:     %s\n", r$desc))
  cat(sprintf("     GO:       %s\n\n", go_ann))
}

up_both   <- hc[hc$drought_dir=="UP"   & hc$salt_dir=="UP",   ]
down_both <- hc[hc$drought_dir=="DOWN" & hc$salt_dir=="DOWN",  ]
up_both   <- up_both[order(-up_both$combined_lfc), ]
down_both <- down_both[order(-down_both$combined_lfc), ]

cat(sprintf("UP in both stresses:   %d\n", nrow(up_both)))
cat(sprintf("DOWN in both stresses: %d\n\n", nrow(down_both)))

##############################################################################
# NEW DISCOVERIES CSV EXPORT
##############################################################################
cat("Building new discoveries CSVs...\n\n")

novel_full$category <- "unknown"
novel_full$category[novel_full$drought_dir=="UP"   & novel_full$salt_dir=="UP"]   <- "UP_in_both"
novel_full$category[novel_full$drought_dir=="DOWN" & novel_full$salt_dir=="DOWN"] <- "DOWN_in_both"
novel_full$category[novel_full$drought_dir=="UP"   & novel_full$salt_dir=="DOWN"] <- "UP_drought_DOWN_salt"
novel_full$category[novel_full$drought_dir=="DOWN" & novel_full$salt_dir=="UP"]   <- "DOWN_drought_UP_salt"

cat("Looking up GO annotations for all high-confidence genes (~2 min)...\n\n")
hc_full <- novel_full[novel_full$high_conf, ]
hc_full$go_annotation <- sapply(hc_full$gene, function(g) {
  tryCatch({
    res <- AnnotationDbi::select(org.At.tair.db,
                                 keys    = g,
                                 columns = c("GO","ONTOLOGY"),
                                 keytype = "TAIR")
    res <- res[!is.na(res$GO) & res$ONTOLOGY=="BP", ]
    if (nrow(res)==0) return("NO_BP_GO_ANNOTATION")
    go_names <- AnnotationDbi::select(GO.db,
                                      keys    = unique(res$GO),
                                      columns = "TERM",
                                      keytype = "GOID")
    terms <- go_names$TERM[!is.na(go_names$TERM) &
                             go_names$TERM != "biological_process"]
    if (length(terms)==0) return("NO_SPECIFIC_BP_ANNOTATION")
    paste(unique(terms), collapse=" | ")
  }, error=function(e) "DATABASE_LOOKUP_FAILED")
})

transport_pattern <- paste(c(
  "transport","transporter","channel","carrier","permease",
  "efflux","influx","pump","exchanger","antiporter","symporter",
  "facilitator","exporter","importer"
), collapse="|")

hc_full$is_transport_related <- grepl(transport_pattern,
                                      hc_full$go_annotation,
                                      ignore.case=TRUE)
hc_full$combined_lfc_score   <- abs(hc_full$lfc_drought) + abs(hc_full$lfc_salt)
hc_full$concordant           <- hc_full$drought_dir == hc_full$salt_dir
hc_full$discordant           <- hc_full$drought_dir != hc_full$salt_dir
hc_full$has_specific_go      <- !hc_full$go_annotation %in%
  c("NO_BP_GO_ANNOTATION","NO_SPECIFIC_BP_ANNOTATION","DATABASE_LOOKUP_FAILED")

hc_full <- hc_full[order(
  -as.integer(hc_full$is_transport_related),
  -as.integer(hc_full$has_specific_go),
  -hc_full$combined_lfc_score
), ]
hc_full$discovery_rank <- seq_len(nrow(hc_full))

cols_ordered <- c(
  "discovery_rank","gene","category","pattern",
  "lfc_drought","padj_drought","lfc_salt","padj_salt",
  "combined_lfc_score","concordant","discordant",
  "is_transport_related","has_specific_go","go_annotation","desc"
)
cols_ordered <- intersect(cols_ordered, names(hc_full))
hc_export    <- hc_full[, cols_ordered]

cat("Saving discovery CSVs:\n")
save_csv(hc_export,
         "NEW_DISCOVERIES_all_high_confidence.csv")
save_csv(hc_export[hc_export$is_transport_related==TRUE, ],
         "NEW_DISCOVERIES_transport_related.csv")
save_csv(hc_export[hc_export$category=="UP_in_both", ],
         "NEW_DISCOVERIES_up_in_both.csv")
save_csv(hc_export[hc_export$category=="DOWN_in_both", ],
         "NEW_DISCOVERIES_down_in_both.csv")
save_csv(hc_export[hc_export$discordant==TRUE, ],
         "NEW_DISCOVERIES_discordant_stress_specific.csv")
save_csv(novel_full[novel_full$high_conf, ],
         "NOVEL_high_confidence_full_profile.csv")
save_csv(up_both,
         "NOVEL_up_in_both_stresses.csv")
save_csv(down_both,
         "NOVEL_down_in_both_stresses.csv")

##############################################################################
# DYNAMIC POSTER NUMBERS
##############################################################################
drought_kegg_summary <- tryCatch(
  readRDS(file.path(drought_dir, "kegg_summary.rds")),
  error=function(e) {
    warning("Could not load drought kegg_summary.rds — run Script 1 first")
    NULL
  }
)
salt_kegg_summary <- tryCatch(
  readRDS(file.path(salt_dir, "kegg_summary.rds")),
  error=function(e) {
    warning("Could not load salt kegg_summary.rds — run Script 2 first")
    NULL
  }
)

n_kegg_drought <- if (!is.null(drought_kegg_summary))
  drought_kegg_summary$n_pathways_total else NA
n_kegg_salt    <- if (!is.null(salt_kegg_summary))
  salt_kegg_summary$n_pathways_total else NA

cat("\n══════════════════════════════════════════════\n")
cat("  POSTER NUMBERS — all from data, no hardcoding\n")
cat("══════════════════════════════════════════════\n\n")

cat(sprintf("  Total DEGs drought:            %d  (UP:%d / DOWN:%d)\n",
            nrow(drought_sig), nrow(drought_up), nrow(drought_down)))
cat(sprintf("  Total DEGs salt:               %d  (UP:%d / DOWN:%d)\n",
            nrow(salt_sig), nrow(salt_up), nrow(salt_down)))
cat(sprintf("  ABA sig GO terms:              %d\n", aba_sig))
cat(sprintf("  Aquaporin sig GO terms:        %d\n", aqua_sig))
cat(sprintf("  Antiporter sig GO terms:       %d\n", anti_sig))
cat(sprintf("  Novel in both stresses:        %d\n", length(both_ids)))
cat(sprintf("  High-confidence novel:         %d\n", nrow(hc)))
cat(sprintf("  Novel UP in both:              %d\n", nrow(up_both)))
cat(sprintf("  Novel DOWN in both:            %d\n", nrow(down_both)))
cat(sprintf("  Novel discordant:              %d\n",
            sum(hc_export$discordant, na.rm=TRUE)))
cat(sprintf("  Novel transport-annotated:     %d\n",
            sum(hc_export$is_transport_related, na.rm=TRUE)))
cat(sprintf("  KEGG pathways drought:         %s\n",
            ifelse(is.na(n_kegg_drought), "RUN SCRIPT 1 FIRST", n_kegg_drought)))
cat(sprintf("  KEGG pathways salt:            %s\n",
            ifelse(is.na(n_kegg_salt),    "RUN SCRIPT 2 FIRST", n_kegg_salt)))

cat("\n  BEST NOVEL TRANSPORT GENE FOR POSTER:\n\n")
transport_novel <- hc_export[hc_export$is_transport_related==TRUE, ]
if (nrow(transport_novel) > 0) {
  best <- transport_novel[1, ]
  cat(sprintf("  Gene:    %s\n", best$gene))
  cat(sprintf("  Pattern: %s\n", best$pattern))
  cat(sprintf("  Drought: log2FC=%+.2f  padj=%.2e\n", best$lfc_drought, best$padj_drought))
  cat(sprintf("  Salt:    log2FC=%+.2f  padj=%.2e\n", best$lfc_salt,    best$padj_salt))
  cat(sprintf("  GO:      %s\n\n", best$go_annotation))
  cat("  POSTER SENTENCE TEMPLATE:\n")
  cat(sprintf("  '%s is %s in both drought (log2FC=%+.2f) and salt (log2FC=%+.2f)\n",
              best$gene,
              ifelse(best$lfc_drought > 0, "upregulated", "downregulated"),
              best$lfc_drought, best$lfc_salt))
  cat(sprintf("   stressed roots and is annotated with %s,\n", best$go_annotation))
  cat("   making it a high-priority candidate for future functional study.'\n")
} else {
  cat("  No transport-related novel genes found in high-confidence set.\n")
  if (sum(hc_export$has_specific_go, na.rm=TRUE) > 0) {
    best_go <- hc_export[hc_export$has_specific_go==TRUE, ][1, ]
    cat(sprintf("  Best candidate with GO annotation: %s\n", best_go$gene))
    cat(sprintf("  Pattern: %s | GO: %s\n", best_go$pattern, best_go$go_annotation))
  }
}

cat(sprintf("\n  All CSVs saved to: %s\n", out_dir))
cat("Script 4 done.\n")


