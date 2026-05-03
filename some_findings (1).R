###########################################################
# Things to Check
###########################################################

###########################################################
# 1. Hardcoded Poster Numbers (IMPORTANT)
###########################################################
cat("\n════════════════════════════════════════════════════════════════\n")
cat("  POSTER NUMBERS — paste directly\n")
cat("════════════════════════════════════════════════════════════════\n\n")
cat(sprintf("  Total DEGs drought:         %d  (↑%d / ↓%d)\n",
            nrow(drought_sig), nrow(drought_up), nrow(drought_down)))  # ✅ dynamic
cat(sprintf("  Total DEGs salt:            %d  (↑%d / ↓%d)\n",
            nrow(salt_sig), nrow(salt_up), nrow(salt_down)))           # ✅ dynamic
cat(sprintf("  ABA sig GO terms:           %d\n",  aba_sig))           # ✅ dynamic
cat(sprintf("  Aquaporin sig GO terms:     %d\n",  aqua_sig))          # ✅ dynamic
cat(sprintf("  Antiporter sig GO terms:    %d\n",  anti_sig))          # ✅ dynamic
cat(sprintf("  Novel in both stresses:     191\n"))  # ❌ HARDCODED
cat(sprintf("  High-confidence novel:      111\n"))  # ❌ HARDCODED
cat(sprintf("  KEGG pathways drought:      30\n"))   # ❌ HARDCODED
cat(sprintf("  KEGG pathways salt:         30\n"))   # ❌ HARDCODED

# Novel in both stresses — computed in the novel gene deep dive section
cat(sprintf("  Novel in both stresses:     %d\n", length(both_ids)))

# High-confidence novel — computed just below that
cat(sprintf("  High-confidence novel:      %d\n", nrow(hc)))

# KEGG pathways — these would need to be captured from the individual 
# pipeline outputs, e.g. by counting rows in the saved KEGG CSVs,
# or by passing the counts through as variables. Right now there is
# no variable in this script that holds KEGG pathway counts at all
# Hardcoding should always be avoided !!

# The Fix
# Add this after the existing write.csv call for KEGG at the end of the individual pipelines
kegg_summary <- data.frame(
  condition        = "drought",
  n_pathways_up   = ifelse(is.null(up_kegg),   0, nrow(up_kegg)),
  n_pathways_down = ifelse(is.null(down_kegg), 0, nrow(down_kegg)),
  n_pathways_total = ifelse(is.null(up_kegg) && is.null(down_kegg), 0, 
                            nrow(combined_kegg))
)
saveRDS(kegg_summary, "kegg_summary.rds")

# Add to Part 13 of both pipeline scripts
deg_summary <- data.frame(
  condition      = "drought",  # or "salt"
  n_up           = length(up_genes),
  n_down         = length(down_genes),
  n_sig_total    = length(up_genes) + length(down_genes)
)
saveRDS(deg_summary, "deg_summary.rds")

# Load saved summaries from each pipeline
drought_kegg_summary <- tryCatch(
  readRDS(file.path(drought_dir, "kegg_summary.rds")),
  error = function(e) { 
    warning("Could not load drought KEGG summary — run drought pipeline first")
    NULL 
  }
)
salt_kegg_summary <- tryCatch(
  readRDS(file.path(salt_dir, "kegg_summary.rds")),
  error = function(e) { 
    warning("Could not load salt KEGG summary — run salt pipeline first")
    NULL 
  }
)

# Extract KEGG counts safely
n_kegg_drought <- if (!is.null(drought_kegg_summary)) 
  drought_kegg_summary$n_pathways_total else NA
n_kegg_salt    <- if (!is.null(salt_kegg_summary))    
  salt_kegg_summary$n_pathways_total else NA

# Now the poster block — everything dynamic
cat("\n════════════════════════════════════════════════════════════════\n")
cat("  POSTER NUMBERS — paste directly\n")
cat("════════════════════════════════════════════════════════════════\n\n")
cat(sprintf("  Total DEGs drought:         %d  (↑%d / ↓%d)\n",
            nrow(drought_sig), nrow(drought_up), nrow(drought_down)))
cat(sprintf("  Total DEGs salt:            %d  (↑%d / ↓%d)\n",
            nrow(salt_sig), nrow(salt_up), nrow(salt_down)))
cat(sprintf("  ABA sig GO terms:           %d\n",  aba_sig))
cat(sprintf("  Aquaporin sig GO terms:     %d\n",  aqua_sig))
cat(sprintf("  Antiporter sig GO terms:    %d\n",  anti_sig))
cat(sprintf("  Novel in both stresses:     %d\n",  length(both_ids)))  # computed above
cat(sprintf("  High-confidence novel:      %d\n",  nrow(hc)))          # computed above
cat(sprintf("  KEGG pathways drought:      %s\n",  
            ifelse(is.na(n_kegg_drought), "RUN DROUGHT PIPELINE FIRST", 
                   n_kegg_drought)))
cat(sprintf("  KEGG pathways salt:         %s\n",  
            ifelse(is.na(n_kegg_salt),    "RUN SALT PIPELINE FIRST",    
                   n_kegg_salt)))

# Everything that appears on the posted MUST
# 1. Be computed directly from data loaded in the top of the script or
# 2. Loaded from the saved output of another script (readRDS, read.csv, etc.)

###########################################################
# 2. Hypergeometric Test Universe Step (IMPORTANT)
###########################################################
# Current and incorrect
go_counts$pvalue <- sapply(seq_len(nrow(go_counts)), function(i) {
  phyper(go_counts$Gene_Count[i] - 1,
         total_with_go,           # ← K: genes in universe WITH any GO annotation WRONG
         total_genes - total_with_go,  # ← N-K: genes WITHOUT any GO annotation WRONG
         length(gene_list),       # ← n: number drawn (your DEG list)
         lower.tail = FALSE)
})

# PROBLEM: total_with_go is used as K for every GO term, when it should be the number of genes
# in the unvierse annotated with that specific term. total_with_go are genes with ANY BP annotation
# This unvervalues the most biologically interesting terms, biasing against specificity

# For each GO term i, K should be how many universe genes have THAT term
# not how many universe genes have ANY GO annotation
phyper(k - 1,
       K_for_this_term,           # genes in universe with this specific GO term
       N - K_for_this_term,       # genes in universe WITHOUT this specific GO term  
       length(gene_list),         # your DEG list size
       lower.tail = FALSE)

# run_full_go in your hypothesis testing actually gets this right?
for (goid in names(go_list)) {
  go_genes <- unique(go_list[[goid]][!is.na(go_list[[goid]])])
  K <- length(go_genes)   # ✅ specific to this GO term
  k <- sum(gene_list %in% go_genes)
  if (k == 0) next
  pval <- phyper(k-1, K, N-K, n, lower.tail=FALSE)
}

###########################################################
# 3. No Low Count Filtering? (IMPORTANT)
###########################################################
# Filter low counts — PRESENT in mine
dds <- dds[rowSums(counts(dds)) > 10, ]
cat("Genes retained after filtering:", nrow(dds), "\n\n")

###########################################################
# 4. Manual Pipeline vs. cluterProfiler (IMPORTANT)
###########################################################
# From Original 
# Apply log-fold change shrinkage — PRESENT in original
resLFC <- lfcShrink(dds, coef="condition_drought_vs_control", type="apeglm")

# And then used throughout:
resOrdered <- as.data.frame(resLFC[order(resLFC$padj), ])

# New pipelines — raw results used directly, no shrinkage
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "drought", "control"))
# lfcShrink is never called — res is used directly for everything in Foster's

# The consequence of this is here:
hc$combined_lfc <- abs(hc$lfc_drought) + abs(hc$lfc_salt)
hc <- hc[order(-hc$combined_lfc),]
top_novel <- head(hc, 10)

# Without shrinkage, this ranking will be disproportionately influenced 
# by low-count genes with artificially extreme fold changes!!

###########################################################
# 5. Vestigial Code?
###########################################################
# BLOCK 1 — computed but result is never used after this
results <- lapply(split(all_go_ann$TAIR, all_go_ann$GO), function(go_genes) {
  go_genes <- unique(go_genes[!is.na(go_genes)])
  K  <- length(go_genes)
  k  <- sum(gene_list %in% go_genes)
  if (k == 0) return(NULL)
  pval <- phyper(k-1, K, N-K, n, lower.tail=FALSE)
  data.frame(GO_ID      = names(split(all_go_ann$TAIR, all_go_ann$GO))[1], # ← also a bug: always returns index [1]
             Gene_Count = k,
             ...)
})
# `results` is assigned here and then never referenced again ↓

# BLOCK 2 — immediately redoes the same computation correctly
go_list   <- split(all_go_ann$TAIR, all_go_ann$GO)
term_list <- split(all_go_ann$TERM, all_go_ann$GO)

result_rows <- list()
for (goid in names(go_list)) {
  # ... correct version of the same logic
}

# Would remove BLOCK 1

###########################################################
# 6. Consistency Across Parallel Scripts (Minor)
###########################################################
# Drought pipeline — no explicit namespace, vulnerable to conflicts
go_data <- select(org.At.tair.db, ...)

# Salt pipeline — explicit namespace, correct
go_data <- AnnotationDbi::select(org.At.tair.db, ...)

# I would make a script with your shared functions and source into both pipelines
# shared_functions.R
get_go_BP    <- function(...) { ... }
get_kegg_BH  <- function(...) { ... }
save_summary <- function(...) { ... }

# Then at the top of both drought and salt pipeline scripts:
source("shared_functions.R")

###########################################################
# 7. Circular Phyper Call (IMPORTANT)
###########################################################
# In Part 8 of the comparative analysis script:
# get_go_BP_comparison <- function(gene_list, label) {
#   
#   go_counts$pvalue <- sapply(seq_len(nrow(go_counts)), function(i) {
#     phyper(go_counts$Gene_Count[i] - 1,
#            length(gene_list),   # K set to DEG list size
#            N - length(gene_list), # complement of DEG list size
#            length(gene_list),   # n also set to DEG list size
#            lower.tail = FALSE)
#   })

# K and n are both length(gene_list)
# "given I drew n genes, what's the probability of getting k genes from a pool of size n"
# circular logic and will produce meaningless p-values
  
###########################################################
# 8. Fisher vs. Venn Diagrams Intersection Inconsistency (IMPORTANT)
###########################################################
# Part 4 of comparison script
# Fisher's test defines universe as genes tested in BOTH experiments
universe <- intersect(rownames(drought_df[!is.na(drought_df$padj), ]),
                      rownames(salt_df[!is.na(salt_df$padj), ]))
N <- length(universe)

# No restriction to shared universe for Venn Diagrams
drought_up   <- rownames(subset(drought_df, padj < 0.05 & log2FoldChange >  1))
salt_up      <- rownames(subset(salt_df,    padj < 0.05 & log2FoldChange >  1))

# Solution
# Define universe once at the top and use consistently everywhere
universe <- intersect(
  rownames(drought_df[!is.na(drought_df$padj), ]),
  rownames(salt_df[!is.na(salt_df$padj), ])
)

# Then filter all gene sets to this universe
drought_up   <- rownames(subset(drought_df[universe,], 
                                padj < 0.05 & log2FoldChange >  1))
drought_down <- rownames(subset(drought_df[universe,], 
                                padj < 0.05 & log2FoldChange < -1))
salt_up      <- rownames(subset(salt_df[universe,],    
                                padj < 0.05 & log2FoldChange >  1))
salt_down    <- rownames(subset(salt_df[universe,],    
                                padj < 0.05 & log2FoldChange < -1))

