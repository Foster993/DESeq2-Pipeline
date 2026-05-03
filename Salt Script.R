################################################################################
# SCRIPT 2 — Arabidopsis Salt RNA-seq Pipeline
# FIXES: dplyr in cran_pkgs | ashr in bioc_pkgs |
#        make.unique() on heatmap rownames | low-count filter |
#        lfcShrink | correct phyper K per term | K>=10 filter |
#        no hardcoding | saves kegg_summary.rds + deg_summary.rds
################################################################################

# ==============================================================================
# PARAMETERS — only section you ever change
# ==============================================================================
STRESS_IDS   <- c("GSM8375077","GSM8375078","GSM8375079")
CONTROL_IDS  <- c("GSM8375074","GSM8375075","GSM8375076")
CONDITION    <- "salt"
CONTROL      <- "control"
FILE_PATTERN <- "gene_count.*\\.txt\\.gz$"
FILE_TYPE    <- "counts"   # "kallisto" or "counts"
N_CONTROLS   <- length(CONTROL_IDS)
N_STRESS     <- length(STRESS_IDS)

# ==============================================================================
# PART 1: PACKAGES
# ==============================================================================
cat("========================================\n")
cat("PART 1: Environment Setup\n")
cat("========================================\n\n")

# FIX 1: dplyr added to cran_pkgs (was missing — caused crash on select <- dplyr::select)
cran_pkgs <- c("ggplot2","pheatmap","ggrepel","reshape2","dplyr")
for (pkg in cran_pkgs) {
  if (!pkg %in% rownames(installed.packages())) install.packages(pkg)
  else message("OK ", pkg)
}
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")

bioc_pkgs <- c("GEOquery","DESeq2","AnnotationDbi","org.At.tair.db",
               "EnhancedVolcano","GO.db","ashr","KEGGREST")
for (pkg in bioc_pkgs) {
  if (!pkg %in% rownames(installed.packages())) BiocManager::install(pkg)
  else message("OK ", pkg)
}

library(GEOquery);  library(DESeq2);   library(ggplot2)
library(pheatmap);  library(ggrepel);  library(reshape2)
library(AnnotationDbi); library(org.At.tair.db)
library(EnhancedVolcano); library(GO.db); library(KEGGREST)

# Must restore after AnnotationDbi masks dplyr::select
select <- dplyr::select

cat("\nAll libraries loaded\n\n")
setwd("~/Downloads/deseq2_arabidopsis")

# ==============================================================================
# COLORS — *** ONLY SECTION CHANGED ***
# grey = control / not significant / below cutoff
# teal = upregulated salt
# pink = downregulated salt
# ==============================================================================
COL_CONTROL <- "#999999"   # grey  (control bars, not-sig points)
COL_STRESS  <- "#008B8B"   # teal  (upregulated salt)
COL_DARK    <- "#004444"   # dark teal
COL_MID     <- "#40BFBF"   # medium teal
COL_LIGHT   <- "#B3E6E6"   # light teal
COL_PALE    <- "#eeeeee"
COL_GRAY1   <- "#666666"
COL_GRAY2   <- "#999999"   # grey  (not-significant volcano points)
COL_DOWN    <- "#FF69B4"   # pink  (downregulated salt)
COL_UP      <- "#008B8B"   # teal  (upregulated salt)

ANNOT_COLORS <- list(
  condition = setNames(c(COL_CONTROL, COL_STRESS), c(CONTROL, CONDITION)),
  Direction = c("Upregulated"=COL_UP, "Downregulated"=COL_DOWN)
)

theme_scheme <- function(base_size=11) {
  theme_bw(base_size=base_size) +
    theme(
      plot.title    = element_text(hjust=0.5, size=base_size+3,
                                   face="bold", color=COL_DARK),
      plot.subtitle = element_text(hjust=0.5, size=base_size, color=COL_GRAY1),
      axis.title    = element_text(size=base_size+1, color=COL_DARK),
      axis.text     = element_text(size=base_size-1, color=COL_GRAY1),
      panel.border  = element_rect(color=COL_STRESS, linewidth=0.8),
      legend.text   = element_text(size=base_size-1, color=COL_GRAY1),
      legend.title  = element_text(size=base_size, color=COL_DARK),
      strip.background = element_rect(fill=COL_STRESS),
      strip.text    = element_text(color="white", face="bold", size=base_size)
    )
}

# ==============================================================================
# PART 2: DOWNLOAD
# ==============================================================================
cat("========================================\n")
cat("PART 2: Data Download\n")
cat("========================================\n\n")

options(timeout=600)
all_ids  <- c(STRESS_IDS, CONTROL_IDS)
out_name <- paste0("arabidopsis_", CONDITION)
dir.create(out_name, showWarnings=FALSE)

for (id in all_ids) {
  message("Downloading: ", id)
  success <- FALSE; attempts <- 0
  while (!success && attempts < 3) {
    attempts <- attempts + 1
    tryCatch({
      getGEOSuppFiles(id, baseDir=out_name,
                      makeDirectory=TRUE, fetch_files=TRUE)
      success <- TRUE; message("  OK ", id)
    }, error=function(e) {
      message("  Attempt ", attempts, " failed: ", e$message)
      if (attempts < 3) Sys.sleep(5)
    })
  }
}
cat("\nDownload complete\n\n")

# ==============================================================================
# PART 3: COUNT MATRIX
# ==============================================================================
cat("========================================\n")
cat("PART 3: Count Matrix Construction\n")
cat("========================================\n\n")

setwd(out_name)

files <- list.files(".", pattern=FILE_PATTERN, recursive=TRUE, full.names=TRUE)
cat("Found", length(files), "files\n")

expected_n <- N_CONTROLS + N_STRESS
if (length(files) != expected_n)
  stop("Expected ", expected_n, " files, found ", length(files))

expected_names <- c(paste0(CONTROL,   "_", seq_len(N_CONTROLS)),
                    paste0(CONDITION, "_", seq_len(N_STRESS)))

count_list <- lapply(seq_along(files), function(i) {
  f <- files[i]
  cat("  Processing", i, "of", length(files), ":", basename(f), "\n")
  
  if (FILE_TYPE == "kallisto") {
    df         <- read.table(f, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    df$gene_id <- sub("\\.\\d+$", "", df$target_id)
    collapsed  <- aggregate(est_counts ~ gene_id, data=df, FUN=sum)
    counts_col <- round(collapsed$est_counts)
    gene_names <- collapsed$gene_id
  } else {
    df         <- read.table(f, header=FALSE, sep="\t", stringsAsFactors=FALSE)
    counts_col <- as.integer(df[, 2])
    gene_names <- df[, 1]
  }
  
  out_df           <- data.frame(count=counts_col, row.names=gene_names)
  colnames(out_df) <- expected_names[i]
  out_df
})

count_matrix <- count_list[[1]]
for (i in 2:length(count_list)) {
  count_matrix <- merge(count_matrix, count_list[[i]],
                        by="row.names", all=TRUE)
  rownames(count_matrix) <- count_matrix$Row.names
  count_matrix            <- count_matrix[, -1, drop=FALSE]
}
count_matrix[is.na(count_matrix)] <- 0
count_matrix[] <- lapply(count_matrix, as.numeric)
count_matrix   <- round(as.data.frame(count_matrix))
count_matrix   <- count_matrix[, expected_names, drop=FALSE]

cat("\nCount matrix:", nrow(count_matrix), "genes x", ncol(count_matrix), "samples\n\n")

# ==============================================================================
# PART 4: METADATA
# ==============================================================================
col_data <- data.frame(
  sample    = expected_names,
  condition = factor(c(rep(CONTROL,   N_CONTROLS),
                       rep(CONDITION, N_STRESS)),
                     levels=c(CONTROL, CONDITION)),
  row.names = expected_names
)
print(col_data)
stopifnot(identical(colnames(count_matrix), rownames(col_data)))
cat("\nMetadata validated\n\n")

# ==============================================================================
# PART 5: DESeq2 OBJECT + LOW COUNT FILTER
# ==============================================================================
cat("========================================\n")
cat("PART 5: DESeq2 Object + Low Count Filter\n")
cat("========================================\n\n")

dds <- DESeqDataSetFromMatrix(countData=count_matrix,
                              colData=col_data,
                              design=~condition)

n_before <- nrow(dds)
keep     <- rowSums(counts(dds)) >= 10
dds      <- dds[keep, ]
n_after  <- nrow(dds)
cat(sprintf("Low-count filter (rowSums >= 10): %d -> %d genes (%d removed)\n\n",
            n_before, n_after, n_before - n_after))

# ==============================================================================
# PART 6: SAVE INTERMEDIATE FILES
# ==============================================================================
write.csv(count_matrix[keep, ], "clean_gene_counts.csv", quote=FALSE)
write.csv(col_data, "sample_metadata.csv", quote=FALSE)
cat("Saved clean_gene_counts.csv + sample_metadata.csv\n\n")

# ==============================================================================
# PART 7: PRE-DESeq2 QC PLOTS (1-4)
# ==============================================================================
cat("========================================\n")
cat("PART 7: QC Plots\n")
cat("========================================\n\n")

vsd_blind  <- vst(dds, blind=TRUE)
pca_data   <- plotPCA(vsd_blind, intgroup="condition", returnData=TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))
cond_cols  <- setNames(c(COL_CONTROL, COL_STRESS), c(CONTROL, CONDITION))
gene_counts <- rowSums(counts(dds))

png("QC_gene_count_distribution.png", width=8, height=6, units="in", res=300)
ggplot(data.frame(x=log10(gene_counts+1)), aes(x=x)) +
  geom_histogram(bins=50, fill=COL_STRESS, color=COL_DARK, alpha=0.85) +
  geom_vline(xintercept=log10(10), linetype="dashed",
             color=COL_CONTROL, linewidth=1) +
  annotate("text", x=log10(10)+0.05, y=Inf, hjust=0, vjust=1.5,
           label="filter cutoff (10)", size=3.5, color=COL_CONTROL) +
  labs(title=sprintf("Gene Count Distribution — %s (post-filter)", CONDITION),
       x="log10(Total Counts + 1)", y="Number of Genes") +
  theme_scheme()
dev.off()
cat("1/7 QC_gene_count_distribution.png\n")

png("QC_PCA_v1_standard.png", width=10, height=8, units="in", res=300)
print(plotPCA(vsd_blind, intgroup="condition") +
        scale_color_manual(values=cond_cols) +
        ggtitle(sprintf("PCA: %s vs %s", CONDITION, CONTROL)) +
        theme_scheme(13))
dev.off()
cat("2/7 QC_PCA_v1_standard.png\n")

png("QC_PCA_v2_labeled.png", width=10, height=8, units="in", res=300)
ggplot(pca_data, aes(x=PC1, y=PC2, color=condition, label=name)) +
  geom_point(size=5, alpha=0.9) +
  geom_text_repel(size=4, color=COL_DARK) +
  scale_color_manual(values=cond_cols) +
  labs(title=sprintf("PCA Labeled: %s vs %s", CONDITION, CONTROL),
       x=paste0("PC1: ",percentVar[1],"% variance"),
       y=paste0("PC2: ",percentVar[2],"% variance")) +
  theme_scheme(13)
dev.off()
cat("3/7 QC_PCA_v2_labeled.png\n")

samp_dists <- dist(t(assay(vsd_blind)))
samp_mat   <- as.matrix(samp_dists)
tryCatch({
  png("QC_sample_distance_heatmap.png", width=8, height=7, units="in", res=150)
  pheatmap(samp_mat,
           clustering_distance_rows=samp_dists,
           clustering_distance_cols=samp_dists,
           annotation_col=col_data["condition"],
           annotation_row=col_data["condition"],
           annotation_colors=ANNOT_COLORS,
           color=colorRampPalette(c(COL_STRESS, COL_PALE, COL_CONTROL))(100),
           main="Sample Distance Heatmap (QC)")
  dev.off()
  cat("4/7 QC_sample_distance_heatmap.png\n\n")
}, error=function(e) {
  try(dev.off(), silent=TRUE)
  cat(sprintf("  SKIPPED QC_sample_distance_heatmap: %s\n\n", e$message))
})

# ==============================================================================
# PART 8: RUN DESeq2 + LFC SHRINKAGE + POST-QC (5-7)
# ==============================================================================
cat("========================================\n")
cat("PART 8: DESeq2 + lfcShrink\n")
cat("========================================\n\n")

dds     <- DESeq(dds)
res_raw <- results(dds, contrast=c("condition", CONDITION, CONTROL))

# lfcShrink stabilises fold changes for low/moderate count genes.
# Using ashr — reliable and does not require apeglm.
res <- tryCatch({
  cat("  Using ashr shrinkage\n")
  lfcShrink(dds, contrast=c("condition",CONDITION,CONTROL), type="ashr")
}, error=function(e) {
  cat(sprintf("  lfcShrink failed (%s) — using unshrunken\n", e$message))
  res_raw
})

cat(sprintf("\nDESeq2 complete. Genes tested: %d\n", nrow(res)))
summary(res)

png("QC_dispersion_plot.png", width=10, height=6, units="in", res=300)
plotDispEsts(dds, main="Dispersion Estimates")
dev.off()
cat("5/7 QC_dispersion_plot.png\n")

# p-value distribution uses res_raw — lfcShrink does not preserve raw p-values
png("QC_pvalue_distribution.png", width=8, height=6, units="in", res=300)
ggplot(data.frame(p=res_raw$pvalue[!is.na(res_raw$pvalue)]), aes(x=p)) +
  geom_histogram(bins=50, fill=COL_STRESS, color=COL_DARK,
                 alpha=0.85, boundary=0) +
  labs(title="P-value Distribution (unshrunken raw p-values)",
       x="Raw p-value", y="Frequency") +
  theme_scheme()
dev.off()
cat("6/7 QC_pvalue_distribution.png\n")

# Independent filtering uses res_raw — filterNumRej metadata only on unshrunken results
filt_df <- as.data.frame(metadata(res_raw)$filterNumRej)
opt_t   <- filt_df$theta[which.max(filt_df$numRej)]
png("QC_independent_filtering.png", width=8, height=6, units="in", res=300)
ggplot(filt_df, aes(x=theta, y=numRej)) +
  geom_line(color=COL_STRESS, linewidth=0.8) +
  geom_point(shape=21, color=COL_STRESS, fill=COL_LIGHT, size=3) +
  geom_vline(xintercept=opt_t, linetype="dashed",
             color=COL_CONTROL, linewidth=1) +
  labs(title="Independent Filtering",
       subtitle=paste("Optimal threshold:", round(opt_t, 2)),
       x="Filter Quantile", y="Rejections") +
  theme_scheme()
dev.off()
cat("7/7 QC_independent_filtering.png\n\n")

res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df      <- res_df[!is.na(res_df$padj) & !is.na(res_df$log2FoldChange), ]

res_df_filtered <- res_df[res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, ]
sig_genes  <- rownames(res_df_filtered)
up_genes   <- rownames(res_df_filtered[res_df_filtered$log2FoldChange >  1, ])
down_genes <- rownames(res_df_filtered[res_df_filtered$log2FoldChange < -1, ])

universe_genes <- res_df$gene
N_universe     <- length(universe_genes)

cat(sprintf("DEGs: %d total (%d UP / %d DOWN)\n\n",
            length(sig_genes), length(up_genes), length(down_genes)))

vsd <- varianceStabilizingTransformation(dds, blind=FALSE)

# ==============================================================================
# PART 9: RESULTS PLOTS (1-4)
# *** CHANGE: significant column now has 3 levels (Up / Down / Not Significant)
# *** CHANGE: scale_color_manual updated to use 3 colors
# ==============================================================================
cat("========================================\n")
cat("PART 9: Results Visualisation\n")
cat("========================================\n\n")

# *** CHANGED: 3-level significant column for direction-aware coloring
res_df$significant <- ifelse(
  res_df$padj < 0.05 & res_df$log2FoldChange >  1, "Upregulated",
  ifelse(
    res_df$padj < 0.05 & res_df$log2FoldChange < -1, "Downregulated",
    "Not Significant"))

png("RESULTS_volcano_plot.png", width=10, height=8, units="in", res=300)
ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), color=significant)) +
  geom_point(alpha=0.6, size=1.5) +
  # *** CHANGED: 3-color scale
  scale_color_manual(values=c("Not Significant"=COL_GRAY2,
                              "Upregulated"    =COL_UP,
                              "Downregulated"  =COL_DOWN)) +
  geom_vline(xintercept=c(-1,1), linetype="dashed", color=COL_DARK) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color=COL_DARK) +
  labs(title=sprintf("Volcano: %s vs %s (shrunken LFC)", CONDITION, CONTROL),
       x="log2 Fold Change (shrunken)", y="-log10(padj)") +
  theme_scheme()
dev.off()
cat("1/4 RESULTS_volcano_plot.png\n")

top15 <- rownames(res_df)[order(res_df$padj)][1:min(15, nrow(res_df))]
top15_syms <- tryCatch(
  AnnotationDbi::mapIds(org.At.tair.db, keys=top15,
                        column="SYMBOL", keytype="TAIR", multiVals="first"),
  error=function(e) setNames(top15, top15))
top15_syms[is.na(top15_syms)] <- top15[is.na(top15_syms)]
res_ev       <- res_df
res_ev$label <- ifelse(rownames(res_ev) %in% top15,
                       top15_syms[rownames(res_ev)], "")

# Build directional color vector: teal=UP, pink=DOWN, grey=not sig
ev_colors <- ifelse(
  res_ev$padj < 0.05 & res_ev$log2FoldChange >  1, COL_UP,
  ifelse(
    res_ev$padj < 0.05 & res_ev$log2FoldChange < -1, COL_DOWN,
    COL_GRAY2))
names(ev_colors) <- rownames(res_ev)

tryCatch({
  png("RESULTS_volcano_enhanced_labeled.png",
      width=14, height=10, units="in", res=150)
  print(
    ggplot(res_ev, aes(x=log2FoldChange, y=-log10(padj),
                       color=significant, label=label)) +
      geom_point(alpha=0.6, size=1.5) +
      scale_color_manual(
        name   = "Expression",
        values = c("Not Significant"=COL_GRAY2,
                   "Upregulated"    =COL_UP,
                   "Downregulated"  =COL_DOWN)) +
      geom_text_repel(
        data          = res_ev[res_ev$label != "", ],
        size          = 3.5,
        fontface      = "bold",
        color         = "black",
        box.padding   = 0.4,
        max.overlaps  = 30,
        segment.color = COL_GRAY1) +
      geom_vline(xintercept=c(-1, 1),
                 linetype="dashed", color=COL_DARK, linewidth=0.6) +
      geom_hline(yintercept=-log10(0.05),
                 linetype="dashed", color=COL_DARK, linewidth=0.6) +
      labs(title    = sprintf("%s vs %s — Labeled Volcano", CONDITION, CONTROL),
           subtitle = "Shrunken LFC | Top 15 genes labeled | padj<0.05, |FC|>1",
           x        = "log2 Fold Change (shrunken)",
           y        = "-log10(padj)") +
      theme_scheme(13)
  )
  dev.off()
  cat("2/4 RESULTS_volcano_enhanced_labeled.png\n")
}, error=function(e) {
  try(dev.off(), silent=TRUE)
  cat(sprintf("  SKIPPED labeled volcano: %s\n", e$message))
})

png("RESULTS_MA_plot.png", width=10, height=8, units="in", res=300)
plotMA(res,
       main      = sprintf("MA Plot: %s vs %s (shrunken LFC)", CONDITION, CONTROL),
       colNonSig = COL_GRAY2,
       colSig    = COL_STRESS,
       colLine   = COL_CONTROL,
       alpha     = 0.05)
dev.off()
cat("3/4 RESULTS_MA_plot.png\n")

top50_genes <- head(sig_genes[order(res_df[sig_genes, "padj"])], 50)
top50_genes <- top50_genes[top50_genes %in% rownames(assay(vsd))]
if (length(top50_genes) >= 3) {
  mat50    <- assay(vsd)[top50_genes, ]
  mat50_sc <- t(scale(t(mat50)))
  row_s    <- tryCatch(
    AnnotationDbi::mapIds(org.At.tair.db, keys=rownames(mat50_sc),
                          column="SYMBOL", keytype="TAIR", multiVals="first"),
    error=function(e) setNames(rownames(mat50_sc), rownames(mat50_sc)))
  row_s[is.na(row_s)] <- rownames(mat50_sc)[is.na(row_s)]
  rownames(mat50_sc)  <- make.unique(as.character(row_s))
  tryCatch({
    png("RESULTS_heatmap_top50.png", width=10, height=14, units="in", res=150)
    pheatmap(mat50_sc,
             annotation_col   = col_data["condition"],
             annotation_colors= ANNOT_COLORS,
             color=colorRampPalette(c(COL_CONTROL,COL_PALE,COL_STRESS))(100),
             border_color=NA,
             main=sprintf("Top 50 DEGs — %s", CONDITION))
    dev.off()
    cat("4/4 RESULTS_heatmap_top50.png\n\n")
  }, error=function(e) {
    try(dev.off(), silent=TRUE)
    cat(sprintf("  SKIPPED RESULTS_heatmap_top50: %s\n\n", e$message))
  })
}

# ==============================================================================
# PART 10: MASSIVE HEATMAP
# ==============================================================================
cat("========================================\n")
cat("PART 10: Massive Expression Heatmap\n")
cat("========================================\n\n")

genes_hm <- sig_genes[sig_genes %in% rownames(assay(vsd))]
if (length(genes_hm) >= 5) {
  MAX_HM_GENES <- 500
  if (length(genes_hm) > MAX_HM_GENES) {
    genes_hm <- head(genes_hm[order(res_df[genes_hm, "padj"])], MAX_HM_GENES)
    cat(sprintf("  Heatmap capped at top %d genes by padj\n", MAX_HM_GENES))
  }
  mat_hm <- assay(vsd)[genes_hm, ]
  mat_sc <- t(scale(t(mat_hm)))
  mat_sc <- pmin(pmax(mat_sc, -3), 3)
  row_ann    <- data.frame(
    Direction = ifelse(genes_hm %in% up_genes, "Upregulated", "Downregulated"),
    row.names = genes_hm)
  ann_col_hm <- data.frame(condition=col_data$condition,
                           row.names=colnames(mat_sc))
  dyn_h <- min(max(14, length(genes_hm) * 0.08), 120)
  tryCatch({
    png("EXPRESSION_massive_heatmap_clustered.png",
        width=12, height=dyn_h, units="in", res=150)
    pheatmap(mat_sc,
             annotation_col    = ann_col_hm,
             annotation_row    = row_ann,
             annotation_colors = ANNOT_COLORS,
             color = colorRampPalette(c(COL_CONTROL, COL_PALE, COL_STRESS))(100),
             border_color=NA,
             show_rownames = length(genes_hm) <= 200,
             fontsize_row=6, cutree_rows=4, cutree_cols=2,
             treeheight_row=60, treeheight_col=30,
             main=sprintf("Top %d DEGs — %s", length(genes_hm), CONDITION))
    dev.off()
    cat(sprintf("EXPRESSION_massive_heatmap_clustered.png (%d genes)\n\n",
                length(genes_hm)))
  }, error=function(e) {
    try(dev.off(), silent=TRUE)
    cat(sprintf("  SKIPPED EXPRESSION_massive_heatmap: %s\n\n", e$message))
  })
}

# ==============================================================================
# PART 11: GO ENRICHMENT — CORRECTED
# ==============================================================================
cat("========================================\n")
cat("PART 11: GO Enrichment (corrected)\n")
cat("========================================\n\n")

cat("  Building universe GO map...\n")
universe_go_ann <- tryCatch(
  AnnotationDbi::select(org.At.tair.db,
                        keys    = universe_genes,
                        columns = c("GO","ONTOLOGY"),
                        keytype = "TAIR"),
  error=function(e) { cat("  WARNING: universe GO lookup failed\n"); NULL }
)

k_total_go <- NULL
if (!is.null(universe_go_ann)) {
  univ_go_bp <- universe_go_ann[
    !is.na(universe_go_ann$GO) & universe_go_ann$ONTOLOGY=="BP", ]
  k_total_go <- table(univ_go_bp$GO)
  cat(sprintf("  Universe BP terms: %d\n\n", length(k_total_go)))
}

cat("  Building universe KEGG map...\n")
universe_kegg_ann <- tryCatch(
  AnnotationDbi::select(org.At.tair.db,
                        keys    = universe_genes,
                        columns = "PATH",
                        keytype = "TAIR"),
  error=function(e) { cat("  WARNING: universe KEGG lookup failed\n"); NULL }
)

k_total_kegg <- NULL
if (!is.null(universe_kegg_ann)) {
  univ_kegg_clean <- universe_kegg_ann[!is.na(universe_kegg_ann$PATH), ]
  k_total_kegg    <- table(univ_kegg_clean$PATH)
  cat(sprintf("  Universe KEGG pathways: %d\n\n", length(k_total_kegg)))
}

build_kegg_names_for_ids <- function(pathway_ids) {
  ids      <- unique(trimws(as.character(pathway_ids)))
  name_map <- setNames(ids, ids)
  
  if (!requireNamespace("KEGGREST", quietly=TRUE)) {
    cat("  KEGGREST not available — pathway IDs will be shown\n")
    return(name_map)
  }
  
  batches    <- split(ids, ceiling(seq_along(ids) / 10))
  n_resolved <- 0L
  
  for (batch in batches) {
    keys <- paste0("ath", batch)
    tryCatch({
      res_list <- KEGGREST::keggGet(keys)
      for (j in seq_along(res_list)) {
        r     <- res_list[[j]]
        input_id <- batch[j]
        pname <- tryCatch(trimws(as.character(r[["NAME"]])[1]),
                          error = function(e) NA_character_)
        pname <- sub(" - .*", "", pname)
        pname <- trimws(pname)
        if (!is.na(pname) && nchar(pname) > 0 && pname != input_id) {
          name_map[input_id] <- pname
          n_resolved <- n_resolved + 1L
        }
      }
    }, error = function(e) NULL)
  }
  
  cat(sprintf("  KEGG names resolved: %d / %d pathways\n",
              n_resolved, length(ids)))
  name_map
}

add_kegg_names <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)
  lkp              <- build_kegg_names_for_ids(df$KEGG_Pathway)
  clean_ids        <- trimws(as.character(df$KEGG_Pathway))
  df$Pathway_Name  <- ifelse(
    is.na(lkp[clean_ids]) | lkp[clean_ids] == clean_ids,
    clean_ids,
    lkp[clean_ids]
  )
  df
}

get_go_BP <- function(gene_list, category_name) {
  if (length(gene_list)==0 || is.null(k_total_go)) return(NULL)
  tryCatch({
    go_data <- AnnotationDbi::select(
      org.At.tair.db,
      keys    = as.character(gene_list),
      columns = c("GO","ONTOLOGY"),
      keytype = "TAIR")
    go_data <- go_data[!is.na(go_data$GO) & go_data$ONTOLOGY=="BP", ]
    if (nrow(go_data)==0) return(NULL)
    
    go_counts           <- as.data.frame(table(go_data$GO))
    colnames(go_counts) <- c("GO_ID","Gene_Count")
    
    go_terms <- AnnotationDbi::select(
      GO.db, keys=as.character(go_counts$GO_ID),
      columns="TERM", keytype="GOID")
    go_counts <- merge(go_counts, go_terms,
                       by.x="GO_ID", by.y="GOID", all.x=TRUE)
    go_counts$TERM[is.na(go_counts$TERM)] <- "Unknown BP term"
    
    go_counts$K_total <- as.integer(k_total_go[as.character(go_counts$GO_ID)])
    go_counts <- go_counts[!is.na(go_counts$K_total) & go_counts$K_total >= 10, ]
    if (nrow(go_counts)==0) return(NULL)
    
    n_deg <- length(gene_list)
    
    go_counts$pvalue <- mapply(
      function(overlap, k_total) {
        phyper(overlap - 1,
               k_total,
               N_universe - k_total,
               n_deg,
               lower.tail=FALSE)
      },
      go_counts$Gene_Count,
      go_counts$K_total
    )
    
    go_counts$padj_BH   <- p.adjust(go_counts$pvalue, method="BH")
    go_counts$log10padj <- -log10(go_counts$padj_BH + 1e-300)
    go_counts$GeneRatio <- go_counts$Gene_Count / n_deg
    go_counts$Category  <- category_name
    go_counts           <- go_counts[order(go_counts$padj_BH), ]
    
    cat(sprintf("  %s: %d sig GO terms (padj<0.05)\n",
                category_name, sum(go_counts$padj_BH < 0.05, na.rm=TRUE)))
    head(go_counts, 20)
  }, error=function(e) { warning("GO error: ", e$message); NULL })
}

up_go   <- get_go_BP(up_genes,   "Upregulated")
down_go <- get_go_BP(down_genes, "Downregulated")

go_dotplot_fn <- function(go_df, title_str, hi_color) {
  ggplot(go_df, aes(x=GeneRatio, y=reorder(TERM,GeneRatio),
                    size=Gene_Count, color=log10padj)) +
    geom_point(alpha=0.85) +
    scale_color_gradient(low=COL_LIGHT, high=hi_color, name="-log10(BH padj)") +
    scale_size_continuous(range=c(3,12), name="Gene Count") +
    labs(title=title_str,
         subtitle="BP | BH | K per term from universe | K>=10",
         x="Gene Ratio", y="GO Term") +
    theme_scheme(11) + theme(axis.text.y=element_text(size=9))
}

if (!is.null(up_go) && nrow(up_go)>0) {
  png("GO_BP_upregulated_dotplot.png",  width=14, height=10, units="in", res=300)
  print(go_dotplot_fn(up_go, sprintf("GO BP Upregulated — %s", CONDITION), COL_UP))
  dev.off()
  
  png("GO_BP_upregulated_gradient.png", width=14, height=10, units="in", res=300)
  print(ggplot(up_go, aes(x=Gene_Count, y=reorder(TERM,Gene_Count),
                          size=Gene_Count, color=log10padj)) +
          geom_point(alpha=0.85) +
          scale_color_gradient(low=COL_LIGHT, high=COL_UP, name="-log10(BH padj)") +
          scale_size_continuous(range=c(3,12)) +
          labs(title=sprintf("GO BP Gradient UP — %s", CONDITION),
               x="Gene Count", y="GO Term") +
          theme_scheme(11) + theme(axis.text.y=element_text(size=9)))
  dev.off()
  
  png("GO_BP_upregulated_lollipop.png", width=14, height=10, units="in", res=300)
  print(ggplot(up_go, aes(x=Gene_Count, y=reorder(TERM,Gene_Count))) +
          geom_segment(aes(x=0, xend=Gene_Count,
                           y=reorder(TERM,Gene_Count),
                           yend=reorder(TERM,Gene_Count)),
                       color=COL_LIGHT, linewidth=0.6) +
          geom_point(aes(size=Gene_Count, color=log10padj), alpha=0.9) +
          scale_color_gradient2(low=COL_MID, mid=COL_PALE, high=COL_UP,
                                midpoint=median(up_go$log10padj)) +
          scale_size_continuous(range=c(4,14)) +
          labs(title=sprintf("GO BP Lollipop UP — %s", CONDITION),
               x="Gene Count", y="") +
          theme_scheme(11) + theme(axis.text.y=element_text(size=9),
                                   panel.grid.major.y=element_blank()))
  dev.off()
  
  write.csv(up_go, "GO_BP_upregulated_top20.csv", row.names=FALSE)
  cat("GO UP plots saved (1-3/6)\n")
} else cat("No UP GO results\n")

if (!is.null(down_go) && nrow(down_go)>0) {
  png("GO_BP_downregulated_dotplot.png",  width=14, height=10, units="in", res=300)
  print(go_dotplot_fn(down_go, sprintf("GO BP Downregulated — %s", CONDITION), COL_DOWN))
  dev.off()
  
  png("GO_BP_downregulated_gradient.png", width=14, height=10, units="in", res=300)
  print(ggplot(down_go, aes(x=Gene_Count, y=reorder(TERM,Gene_Count),
                            size=Gene_Count, color=log10padj)) +
          geom_point(alpha=0.85) +
          scale_color_gradient(low=COL_LIGHT, high=COL_DOWN, name="-log10(BH padj)") +
          scale_size_continuous(range=c(3,12)) +
          labs(title=sprintf("GO BP Gradient DOWN — %s", CONDITION),
               x="Gene Count", y="GO Term") +
          theme_scheme(11) + theme(axis.text.y=element_text(size=9)))
  dev.off()
  
  png("GO_BP_downregulated_lollipop.png", width=14, height=10, units="in", res=300)
  print(ggplot(down_go, aes(x=Gene_Count, y=reorder(TERM,Gene_Count))) +
          geom_segment(aes(x=0, xend=Gene_Count,
                           y=reorder(TERM,Gene_Count),
                           yend=reorder(TERM,Gene_Count)),
                       color=COL_LIGHT, linewidth=0.6) +
          geom_point(aes(size=Gene_Count, color=log10padj), alpha=0.9) +
          scale_color_gradient2(low=COL_MID, mid=COL_PALE, high=COL_DOWN,
                                midpoint=median(down_go$log10padj)) +
          scale_size_continuous(range=c(4,14)) +
          labs(title=sprintf("GO BP Lollipop DOWN — %s", CONDITION),
               x="Gene Count", y="") +
          theme_scheme(11) + theme(axis.text.y=element_text(size=9),
                                   panel.grid.major.y=element_blank()))
  dev.off()
  
  write.csv(down_go, "GO_BP_downregulated_top20.csv", row.names=FALSE)
  cat("GO DOWN plots saved (4-6/6)\n\n")
} else cat("No DOWN GO results\n\n")

# ==============================================================================
# PART 12: KEGG ENRICHMENT — CORRECTED
# ==============================================================================
cat("========================================\n")
cat("PART 12: KEGG Enrichment (corrected)\n")
cat("========================================\n\n")

get_kegg_BH <- function(gene_list, direction) {
  if (length(gene_list)==0 || is.null(k_total_kegg)) return(NULL)
  tryCatch({
    kegg_data <- AnnotationDbi::select(
      org.At.tair.db,
      keys    = as.character(gene_list),
      columns = "PATH",
      keytype = "TAIR")
    kegg_data <- kegg_data[!is.na(kegg_data$PATH), ]
    if (nrow(kegg_data)==0) return(NULL)
    
    k_counts           <- as.data.frame(table(kegg_data$PATH))
    colnames(k_counts) <- c("KEGG_Pathway","Gene_Count")
    k_counts$K_total   <- as.integer(k_total_kegg[as.character(k_counts$KEGG_Pathway)])
    k_counts <- k_counts[!is.na(k_counts$K_total) & k_counts$K_total >= 10, ]
    if (nrow(k_counts)==0) return(NULL)
    
    n_deg <- length(gene_list)
    
    k_counts$pvalue <- mapply(
      function(overlap, k_total) {
        phyper(overlap - 1,
               k_total,
               N_universe - k_total,
               n_deg,
               lower.tail=FALSE)
      },
      k_counts$Gene_Count,
      k_counts$K_total
    )
    
    k_counts$padj_BH   <- p.adjust(k_counts$pvalue, method="BH")
    k_counts$log10padj <- -log10(k_counts$padj_BH + 1e-300)
    k_counts$GeneRatio <- k_counts$Gene_Count / n_deg
    k_counts$Direction <- direction
    k_counts           <- k_counts[order(k_counts$padj_BH), ]
    
    cat(sprintf("  KEGG %s: %d sig (padj<0.05)\n",
                direction, sum(k_counts$padj_BH < 0.05, na.rm=TRUE)))
    head(k_counts, 15)
  }, error=function(e) { warning("KEGG error: ", e$message); NULL })
}

up_kegg   <- get_kegg_BH(up_genes,   "Upregulated")
down_kegg <- get_kegg_BH(down_genes, "Downregulated")

up_kegg   <- add_kegg_names(up_kegg)
down_kegg <- add_kegg_names(down_kegg)

kegg_parts <- list()
if (!is.null(up_kegg)   && nrow(up_kegg)>0)   kegg_parts[["Upregulated"]]   <- up_kegg
if (!is.null(down_kegg) && nrow(down_kegg)>0)  kegg_parts[["Downregulated"]] <- down_kegg

if (length(kegg_parts)>0) {
  combined_kegg <- do.call(rbind, kegg_parts)
  
  png("KEGG_dotplot_combined.png", width=14, height=10, units="in", res=300)
  print(ggplot(combined_kegg,
               aes(x=GeneRatio, y=reorder(Pathway_Name,GeneRatio),
                   size=Gene_Count, color=log10padj)) +
          geom_point(alpha=0.85) +
          scale_color_gradient(low=COL_LIGHT, high=COL_DARK,
                               name="-log10(BH padj)") +
          scale_size_continuous(range=c(3,12)) +
          facet_wrap(~Direction, scales="free_y", ncol=2) +
          labs(title=sprintf("KEGG Pathways: %s", CONDITION),
               subtitle="BH correction | Top 15 per direction",
               x="Gene Ratio", y="") +
          theme_scheme(11) + theme(axis.text.y=element_text(size=9)))
  dev.off()
  
  png("KEGG_bubble_combined.png", width=14, height=10, units="in", res=300)
  print(ggplot(combined_kegg,
               aes(x=Gene_Count, y=reorder(Pathway_Name,Gene_Count),
                   color=Direction, size=Gene_Count)) +
          geom_point(alpha=0.75) +
          scale_color_manual(
            values=c("Upregulated"=COL_UP,"Downregulated"=COL_DOWN)) +
          scale_size_continuous(range=c(3,10)) +
          labs(title=sprintf("KEGG Bubble — %s", CONDITION),
               x="Gene Count", y="") +
          theme_scheme(11) + theme(axis.text.y=element_text(size=9)))
  dev.off()
  
  png("KEGG_bubble_faceted.png", width=16, height=10, units="in", res=300)
  print(ggplot(combined_kegg,
               aes(x=Gene_Count, y=reorder(Pathway_Name,Gene_Count),
                   color=Direction, size=Gene_Count)) +
          geom_point(alpha=0.85) +
          scale_color_manual(
            values=c("Upregulated"=COL_UP,"Downregulated"=COL_DOWN)) +
          scale_size_continuous(range=c(4,12)) +
          facet_wrap(~Direction, scales="free_y", ncol=2) +
          labs(title=sprintf("KEGG Faceted — %s", CONDITION),
               x="Gene Count", y="") +
          theme_scheme(11) + theme(axis.text.y=element_text(size=9),
                                   legend.position="bottom"))
  dev.off()
  
  write.csv(combined_kegg, "KEGG_pathways_BH_combined.csv", row.names=FALSE)
  cat("KEGG plots saved (1-3/3)\n\n")
} else cat("No KEGG results\n\n")

# ==============================================================================
# PART 13: SAVE RESULTS + SUMMARIES FOR SCRIPT 4
# ==============================================================================
cat("========================================\n")
cat("PART 13: Save Results\n")
cat("========================================\n\n")

saveRDS(res, "deseq2_results.rds")
saveRDS(dds, "dds_object.rds")
saveRDS(vsd, "vsd_object.rds")
write.csv(res_df, "deseq2_results_full.csv")

deg_summary <- data.frame(
  condition       = CONDITION,
  n_up            = length(up_genes),
  n_down          = length(down_genes),
  n_sig_total     = length(sig_genes),
  n_universe      = N_universe,
  stringsAsFactors= FALSE
)
saveRDS(deg_summary, "deg_summary.rds")

kegg_summary <- data.frame(
  condition        = CONDITION,
  n_pathways_up    = ifelse(is.null(up_kegg),   0L, nrow(up_kegg)),
  n_pathways_down  = ifelse(is.null(down_kegg), 0L, nrow(down_kegg)),
  n_pathways_total = ifelse(length(kegg_parts)==0, 0L,
                            nrow(do.call(rbind, kegg_parts))),
  stringsAsFactors = FALSE
)
saveRDS(kegg_summary, "kegg_summary.rds")

cat("Saved: deseq2_results.rds | dds_object.rds | vsd_object.rds\n")
cat("Saved: deg_summary.rds | kegg_summary.rds\n\n")

# ==============================================================================
# FINAL VERIFICATION
# ==============================================================================
cat("========================================\n")
cat("FINAL VERIFICATION\n")
cat("========================================\n\n")

expected_plots <- c(
  "QC_gene_count_distribution.png","QC_PCA_v1_standard.png",
  "QC_PCA_v2_labeled.png","QC_sample_distance_heatmap.png",
  "QC_dispersion_plot.png","QC_pvalue_distribution.png",
  "QC_independent_filtering.png","RESULTS_volcano_plot.png",
  "RESULTS_volcano_enhanced_labeled.png","RESULTS_MA_plot.png",
  "RESULTS_heatmap_top50.png","EXPRESSION_massive_heatmap_clustered.png",
  "GO_BP_upregulated_dotplot.png","GO_BP_upregulated_gradient.png",
  "GO_BP_upregulated_lollipop.png","GO_BP_downregulated_dotplot.png",
  "GO_BP_downregulated_gradient.png","GO_BP_downregulated_lollipop.png",
  "KEGG_dotplot_combined.png","KEGG_bubble_combined.png",
  "KEGG_bubble_faceted.png"
)

all_ok <- TRUE
for (f in expected_plots) {
  if (file.exists(f) && file.info(f)$size > 1000) {
    cat(sprintf("OK  %-45s %.1f KB\n", f, file.info(f)$size/1024))
  } else if (file.exists(f)) {
    cat(sprintf("SMALL  %s\n", f)); all_ok <- FALSE
  } else {
    cat(sprintf("MISSING  %s\n", f)); all_ok <- FALSE
  }
}

cat("\n")
if (all_ok) cat("ALL 21 PLOTS VERIFIED\n") else cat("Some files missing\n")
cat(sprintf("\nLow-count filter: %d -> %d genes\n", n_before, n_after))
cat(sprintf("DEGs: %d (%d UP / %d DOWN)\n",
            length(sig_genes), length(up_genes), length(down_genes)))