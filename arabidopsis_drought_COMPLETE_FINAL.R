################################################################################
# Arabidopsis Drought RNA-seq Complete Analysis Pipeline
# ALL 21 PLOTS | Color scheme applied | BP+BH GO & KEGG | FIXED
################################################################################

# ==============================================================================
# PART 1: ENVIRONMENT SETUP
# ==============================================================================

cat("========================================\n")
cat("PART 1: Environment Setup\n")
cat("========================================\n\n")

if (!"DESeq2" %in% rownames(installed.packages())) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("DESeq2")
}

library(GEOquery)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(org.At.tair.db)
library(EnhancedVolcano)
library(GO.db)

cat("✅ All libraries loaded!\n\n")

setwd("~/Downloads/deseq2_arabidopsis")

# ==============================================================================
# GLOBAL COLOR SCHEME
# Control = #0088cc (blue)  |  Drought = #cc4400 (orange)
# ==============================================================================

COL_CONTROL  <- "#0088cc"   # blue  — all control samples & bars
COL_DROUGHT  <- "#cc4400"   # orange — all drought samples & bars
COL_DARK     <- "#004D74"   # dark navy — titles, borders, outlines
COL_MID      <- "#48AADA"   # mid blue — lines, secondary points
COL_LIGHT    <- "#91CCE9"   # pale blue — fills, lollipop segments
COL_PALE     <- "#eeeeee"   # near-white — heatmap midpoint
COL_GRAY1    <- "#666666"   # dark gray — axis text
COL_GRAY2    <- "#979797"   # mid gray — non-significant points
COL_DOWN     <- "#87cc00"   # triad green — downregulated GO/KEGG
COL_UP       <- "#cc4400"   # orange — upregulated GO/KEGG (matches drought)

# Annotation colors for all pheatmap calls
ANNOT_COLORS <- list(
  condition = c("control" = COL_CONTROL, "drought" = COL_DROUGHT),
  Direction = c("Upregulated" = COL_UP,  "Downregulated" = COL_DOWN)
)

# Shared ggplot2 theme — applied to every plot
theme_scheme <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title       = element_text(hjust = 0.5, size = base_size + 3,
                                      face = "bold", color = COL_DARK),
      plot.subtitle    = element_text(hjust = 0.5, size = base_size,
                                      color = COL_GRAY1),
      axis.title       = element_text(size = base_size + 1, color = COL_DARK),
      axis.text        = element_text(size = base_size - 1, color = COL_GRAY1),
      panel.border     = element_rect(color = COL_CONTROL, linewidth = 0.8),
      legend.text      = element_text(size = base_size - 1, color = COL_GRAY1),
      legend.title     = element_text(size = base_size, color = COL_DARK),
      strip.background = element_rect(fill = COL_CONTROL),
      strip.text       = element_text(color = "white", face = "bold",
                                      size = base_size)
    )
}

# ==============================================================================
# PART 2: DATA DOWNLOAD
# ==============================================================================

cat("========================================\n")
cat("PART 2: Data Download\n")
cat("========================================\n\n")

options(timeout = 600)
options(download.file.method = "libcurl")

drought_ids <- c("GSM8346417","GSM8346418","GSM8346419","GSM8346420")
control_ids <- c("GSM8346402","GSM8346403","GSM8346404","GSM8346405")
all_ids     <- c(drought_ids, control_ids)

dir.create("arabidopsis_drought", showWarnings = FALSE)

for (id in all_ids) {
  message("Downloading: ", id)
  success <- FALSE; attempts <- 0
  while (!success && attempts < 3) {
    attempts <- attempts + 1
    tryCatch({
      getGEOSuppFiles(id, baseDir = "arabidopsis_drought",
                      makeDirectory = TRUE, fetch_files = TRUE)
      success <- TRUE
      message("  ✅ ", id, " downloaded")
    }, error = function(e) {
      message("  ⚠️  Attempt ", attempts, " failed: ", e$message)
      if (attempts < 3) { message("  Retrying in 5s..."); Sys.sleep(5) }
      else message("  ❌ All attempts failed for ", id)
    })
  }
}
cat("\n✅ Download complete!\n\n")

# ==============================================================================
# PART 3: COUNT MATRIX CONSTRUCTION
# ==============================================================================

cat("========================================\n")
cat("PART 3: Count Matrix Construction\n")
cat("========================================\n\n")

setwd("arabidopsis_drought")

abundance_files <- list.files(".", pattern = "\\.abundance\\.tsv\\.gz$",
                               recursive = TRUE, full.names = TRUE)
cat("Found", length(abundance_files), "files:\n")
print(basename(abundance_files)); cat("\n")

if (length(abundance_files) != 8)
  stop("❌ Expected 8 files, found ", length(abundance_files))

count_list <- lapply(seq_along(abundance_files), function(i) {
  f  <- abundance_files[i]
  cat("  Processing", i, "of", length(abundance_files), ":", basename(f), "\n")
  df <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  df$gene   <- sub("\\.\\d+$", "", df$target_id)
  collapsed <- aggregate(df$est_counts, by = list(gene = df$gene), FUN = sum)
  rownames(collapsed) <- collapsed$gene
  collapsed <- collapsed[, "x", drop = FALSE]
  colnames(collapsed) <- basename(dirname(f))
  collapsed
})

count_matrix <- count_list[[1]]
for (i in 2:length(count_list)) {
  count_matrix <- merge(count_matrix, count_list[[i]], by = "row.names", all = TRUE)
  rownames(count_matrix) <- count_matrix$Row.names
  count_matrix <- count_matrix[, -1, drop = FALSE]
}
count_matrix[is.na(count_matrix)] <- 0
count_matrix[] <- lapply(count_matrix, as.numeric)
count_matrix   <- round(as.data.frame(count_matrix))

expected_names         <- c(paste0("control_", 1:4), paste0("drought_", 1:4))
colnames(count_matrix) <- expected_names
cat("\n✅ Count matrix built!\n\n")

# ==============================================================================
# PART 4: METADATA CREATION
# ==============================================================================

cat("========================================\n")
cat("PART 4: Metadata Creation\n")
cat("========================================\n\n")

col_data <- data.frame(
  sample    = expected_names,
  condition = factor(c(rep("control", 4), rep("drought", 4)),
                     levels = c("control", "drought")),
  row.names = expected_names
)
print(col_data)

cat("\n--- Sanity Checks ---\n")
cat("Column names match metadata rownames:",
    identical(colnames(count_matrix), rownames(col_data)), "\n")
cat("Number of genes:", nrow(count_matrix), "\n")
cat("Number of samples:", ncol(count_matrix), "\n\n")

if (!identical(colnames(count_matrix), rownames(col_data)))
  stop("❌ Metadata doesn't match count matrix!")

cat("✅ All sanity checks passed!\n\n")

# ==============================================================================
# PART 5: CREATE DESEQ2 OBJECT
# ==============================================================================

cat("========================================\n")
cat("PART 5: DESeq2 Object Creation\n")
cat("========================================\n\n")

dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                               colData   = col_data,
                               design    = ~ condition)
cat("✅ DESeq2 dataset created\n\n")
print(dds); cat("\n")

# ==============================================================================
# PART 6: SAVE INTERMEDIATE FILES
# ==============================================================================

cat("========================================\n")
cat("PART 6: Saving Intermediate Files\n")
cat("========================================\n\n")

write.csv(count_matrix, "clean_gene_counts.csv", quote = FALSE)
write.csv(col_data,     "sample_metadata.csv",   quote = FALSE)
cat("✅ Saved clean_gene_counts.csv + sample_metadata.csv\n\n")

# ==============================================================================
# PART 7: PRE-DESEQ2 QC PLOTS (plots 1-4)
# NOTE: Dispersion plot (plot 5) is in Part 8 — requires DESeq() to run first
# ==============================================================================

cat("========================================\n")
cat("PART 7: Quality Control Plots (Pre-DESeq2)\n")
cat("========================================\n\n")

vsd <- vst(dds, blind = TRUE)

# --- QC Plot 1: Gene Count Distribution Histogram ---
gene_counts <- rowSums(counts(dds))
count_df    <- data.frame(total_counts = gene_counts)

pdf("QC_gene_count_distribution.pdf", width = 8, height = 6)
ggplot(count_df, aes(x = log10(total_counts + 1))) +
  geom_histogram(bins = 50, fill = COL_CONTROL, color = COL_DARK, alpha = 0.85) +
  geom_vline(xintercept = log10(11), linetype = "dashed",
             color = COL_DROUGHT, linewidth = 1) +
  labs(title = "Gene Count Distribution Before Filtering",
       x = "Log10(Total Counts + 1)", y = "Number of Genes") +
  theme_scheme()
dev.off()
cat("✅ 1/7 QC_gene_count_distribution.pdf saved\n")

# --- QC Plot 2: PCA V1 Standard ---
pdf("QC_PCA_v1_standard.pdf", width = 10, height = 8)
plotPCA(vsd, intgroup = "condition") +
  scale_color_manual(values = c("control" = COL_CONTROL,
                                "drought" = COL_DROUGHT)) +
  ggtitle("PCA: Drought vs Control") +
  theme_scheme(base_size = 13) +
  theme(legend.position = "right")
dev.off()
cat("✅ 2/7 QC_PCA_v1_standard.pdf saved\n")

# --- QC Plot 3: PCA V2 Labeled ---
pca_data   <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

pdf("QC_PCA_v2_labeled.pdf", width = 10, height = 8)
ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = name)) +
  geom_point(size = 5, alpha = 0.9) +
  geom_text(vjust = -1, hjust = 0.5, size = 4, color = COL_DARK) +
  scale_color_manual(values = c("control" = COL_CONTROL,
                                "drought" = COL_DROUGHT)) +
  labs(title = "PCA: Drought vs Control (Labeled)",
       x     = paste0("PC1: ", percentVar[1], "% variance"),
       y     = paste0("PC2: ", percentVar[2], "% variance"),
       color = "Condition") +
  theme_scheme(base_size = 13) +
  theme(legend.position = "right")
dev.off()
cat("✅ 3/7 QC_PCA_v2_labeled.pdf saved\n")

# --- QC Plot 4: Sample Distance Heatmap ---
sample_dists       <- dist(t(assay(vsd)))
sample_dist_matrix <- as.matrix(sample_dists)

pdf("QC_sample_distance_heatmap.pdf", width = 8, height = 7)
pheatmap(sample_dist_matrix,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         annotation_col    = col_data["condition"],
         annotation_row    = col_data["condition"],
         annotation_colors = ANNOT_COLORS,
         color    = colorRampPalette(c(COL_CONTROL, COL_PALE, COL_DROUGHT))(100),
         main     = "Sample-to-Sample Distance Heatmap (QC)",
         fontsize = 10)
dev.off()
cat("✅ 4/7 QC_sample_distance_heatmap.pdf saved\n\n")

# ==============================================================================
# PART 8: RUN DESEQ2 + POST-DESeq2 QC PLOTS (plots 5, 6, 7)
# ==============================================================================

cat("========================================\n")
cat("PART 8: DESeq2 Analysis\n")
cat("========================================\n\n")

dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "drought", "control"))

cat("✅ DESeq2 analysis complete!\n\n")
summary(res); cat("\n")

# Build res_df now — needed by GO/KEGG functions downstream
res_df <- as.data.frame(res)
res_df$significant <- ifelse(
  !is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1,
  "Significant", "Not Significant")
res_df_filtered <- res_df[!is.na(res_df$padj), ]

# --- QC Plot 5: Dispersion Plot ---
# Must come after DESeq() — dispersions are NULL on dds until then
pdf("QC_dispersion_plot.pdf", width = 10, height = 6)
plotDispEsts(dds, main = "Dispersion Estimates", legend = TRUE)
dev.off()
cat("✅ 5/7 QC_dispersion_plot.pdf saved\n")

# --- QC Plot 6: P-value Distribution ---
pval_df <- data.frame(pvalue = res$pvalue[!is.na(res$pvalue)])

pdf("QC_pvalue_distribution.pdf", width = 8, height = 6)
ggplot(pval_df, aes(x = pvalue)) +
  geom_histogram(bins = 50, fill = COL_CONTROL, color = COL_DARK,
                 alpha = 0.85, boundary = 0) +
  labs(title = "P-value Distribution", x = "P-value", y = "Frequency") +
  theme_scheme()
dev.off()
cat("✅ 6/7 QC_pvalue_distribution.pdf saved\n")

# --- QC Plot 7: Independent Filtering ---
metadata_res  <- metadata(res)
filter_df     <- data.frame(theta  = metadata_res$filterNumRej$theta,
                            numRej = metadata_res$filterNumRej$numRej)
optimal_theta <- filter_df$theta[which.max(filter_df$numRej)]

pdf("QC_independent_filtering.pdf", width = 8, height = 6)
ggplot(filter_df, aes(x = theta, y = numRej)) +
  geom_line(color = COL_CONTROL, linewidth = 0.8) +
  geom_point(shape = 21, color = COL_CONTROL, fill = COL_LIGHT,
             size = 3, stroke = 1.5) +
  geom_vline(xintercept = optimal_theta, linetype = "dashed",
             color = COL_DROUGHT, linewidth = 1) +
  labs(title    = "Independent Filtering",
       subtitle = paste("Optimal threshold:", round(optimal_theta, 2)),
       x = "Quantiles of Mean Normalized Counts",
       y = "Number of Rejections") +
  theme_scheme()
dev.off()
cat("✅ 7/7 QC_independent_filtering.pdf saved\n\n")

# ==============================================================================
# PART 9: RESULTS VISUALIZATION (4 plots)
# ==============================================================================

cat("========================================\n")
cat("PART 9: Results Visualization\n")
cat("========================================\n\n")

# --- Results Plot 1: Simple Volcano ---
pdf("RESULTS_volcano_plot.pdf", width = 10, height = 8)
ggplot(res_df_filtered,
       aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Not Significant" = COL_GRAY2,
                                "Significant"     = COL_DROUGHT)) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = COL_DARK) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = COL_DARK) +
  labs(title = "Volcano Plot: Drought vs Control",
       x     = "Log2 Fold Change",
       y     = "-Log10(Adjusted P-value)",
       color = "Status") +
  theme_scheme() +
  theme(legend.position = "top")
dev.off()
cat("✅ 1/4 RESULTS_volcano_plot.pdf saved\n")

# --- Results Plot 2: Enhanced Labeled Volcano ---
res_volcano        <- as.data.frame(res)
res_volcano$gene   <- rownames(res_volcano)
res_volcano$symbol <- mapIds(org.At.tair.db,
                             keys      = rownames(res_volcano),
                             column    = "SYMBOL",
                             keytype   = "TAIR",
                             multiVals = "first")
res_volcano$label  <- ifelse(is.na(res_volcano$symbol),
                             res_volcano$gene,
                             paste0(res_volcano$gene,
                                    "\n(", res_volcano$symbol, ")"))
top_genes_to_label <- head(rownames(res[order(res$padj), ]), 15)

pdf("RESULTS_volcano_enhanced_labeled.pdf", width = 14, height = 10)
EnhancedVolcano(res_volcano,
                lab             = res_volcano$label,
                x               = "log2FoldChange",
                y               = "padj",
                selectLab       = res_volcano$label[
                  res_volcano$gene %in% top_genes_to_label],
                title           = "Drought vs Control: Enhanced Volcano",
                subtitle        = "Top 15 genes labeled",
                pCutoff         = 0.05,
                FCcutoff        = 1.0,
                pointSize       = 2.5,
                labSize         = 4.0,
                labCol          = COL_DARK,
                labFace         = "bold",
                boxedLabels     = TRUE,
                colAlpha        = 0.5,
                col             = c(COL_GRAY2, COL_GRAY1, COL_MID, COL_DROUGHT),
                legendPosition  = "right",
                legendLabSize   = 12,
                legendIconSize  = 5.0,
                drawConnectors  = TRUE,
                widthConnectors = 0.5,
                colConnectors   = COL_DARK,
                arrowheads      = FALSE,
                gridlines.major = TRUE,
                gridlines.minor = FALSE,
                border          = "full",
                borderWidth     = 1.5,
                borderColour    = COL_CONTROL,
                xlim            = c(-12, 12))
dev.off()
cat("✅ 2/4 RESULTS_volcano_enhanced_labeled.pdf saved\n")

# --- Results Plot 3: MA Plot ---
pdf("RESULTS_MA_plot.pdf", width = 10, height = 8)
plotMA(res,
       ylim      = c(-10, 10),
       alpha     = 0.05,
       main      = "MA Plot: Drought vs Control",
       colNonSig = COL_GRAY2,
       colSig    = COL_DROUGHT,
       colLine   = COL_CONTROL,
       cex       = 0.8)
legend("topright",
       legend = c("Significant (padj < 0.05)", "Not Significant"),
       col = c(COL_DROUGHT, COL_GRAY2), pch = 16, cex = 0.9)
dev.off()
cat("✅ 3/4 RESULTS_MA_plot.pdf saved\n")

# --- Results Plot 4: Top 50 Heatmap ---
top50        <- head(order(res$padj), 50)
top50_counts <- assay(vsd)[top50, ]

pdf("RESULTS_heatmap_top50.pdf", width = 10, height = 12)
pheatmap(top50_counts,
         cluster_rows      = TRUE,
         cluster_cols      = TRUE,
         show_rownames     = TRUE,
         annotation_col    = col_data["condition"],
         annotation_colors = ANNOT_COLORS,
         scale             = "row",
         color  = colorRampPalette(c(COL_CONTROL, COL_PALE, COL_DROUGHT))(100),
         main   = "Top 50 Differentially Expressed Genes (Drought vs Control)",
         fontsize_row = 8,
         fontsize_col = 10)
dev.off()
cat("✅ 4/4 RESULTS_heatmap_top50.pdf saved\n\n")

# ==============================================================================
# PART 10: MASSIVE EXPRESSION PATTERN HEATMAP (1 plot)
# ==============================================================================

cat("========================================\n")
cat("PART 10: Massive Expression Heatmap\n")
cat("========================================\n\n")

sig_genes  <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
up_genes   <- rownames(sig_genes[sig_genes$log2FoldChange >  1, ])
down_genes <- rownames(sig_genes[sig_genes$log2FoldChange < -1, ])

cat("Upregulated:", length(up_genes), "\n")
cat("Downregulated:", length(down_genes), "\n\n")

mat        <- assay(vsd)[rownames(sig_genes), ]
mat_scaled <- t(scale(t(mat)))

row_annot  <- data.frame(
  Direction = ifelse(rownames(mat_scaled) %in% up_genes,
                     "Upregulated", "Downregulated"),
  row.names = rownames(mat_scaled)
)

plot_height <- min(max(14, nrow(mat_scaled) * 0.05), 60)
show_rows   <- nrow(mat_scaled) <= 200

pdf("EXPRESSION_massive_heatmap_clustered.pdf",
    width = 12, height = plot_height)
pheatmap(mat_scaled,
         cluster_rows      = TRUE,
         cluster_cols      = TRUE,
         show_rownames     = show_rows,
         show_colnames     = TRUE,
         annotation_col    = col_data["condition"],
         annotation_row    = row_annot,
         annotation_colors = ANNOT_COLORS,
         color  = colorRampPalette(c(COL_CONTROL, COL_PALE, COL_DROUGHT))(100),
         main   = paste0("All Significant DEGs — Clustered Expression Patterns\n(",
                         nrow(mat_scaled), " genes | padj<0.05 | |log2FC|>1)"),
         fontsize       = 9,
         fontsize_row   = 6,
         border_color   = NA,
         cutree_rows    = 4,
         cutree_cols    = 2,
         treeheight_row = 60,
         treeheight_col = 30)
dev.off()
cat("✅ EXPRESSION_massive_heatmap_clustered.pdf saved (",
    nrow(mat_scaled), "genes)\n\n")

# ==============================================================================
# PART 11: GO ENRICHMENT — BIOLOGICAL PROCESS ONLY + BH (6 plots)
# ==============================================================================

cat("========================================\n")
cat("PART 11: GO Enrichment (BP + BH)\n")
cat("========================================\n\n")

get_go_BP <- function(gene_list, category_name) {
  if (length(gene_list) == 0) { warning("No genes for ", category_name); return(NULL) }
  tryCatch({
    go_data <- select(org.At.tair.db,
                      keys    = as.character(gene_list),
                      columns = c("GO", "ONTOLOGY"),
                      keytype = "TAIR")
    # Biological Process only
    go_data <- go_data[!is.na(go_data$GO) & go_data$ONTOLOGY == "BP", ]
    if (nrow(go_data) == 0) { warning("No BP terms for ", category_name); return(NULL) }

    go_counts           <- as.data.frame(table(go_data$GO))
    colnames(go_counts) <- c("GO_ID", "Gene_Count")

    go_terms  <- select(GO.db,
                        keys    = as.character(go_counts$GO_ID),
                        columns = "TERM",
                        keytype = "GOID")
    go_counts <- merge(go_counts, go_terms,
                       by.x = "GO_ID", by.y = "GOID", all.x = TRUE)
    go_counts$TERM[is.na(go_counts$TERM)] <- "Unknown BP term"

    total_genes   <- nrow(res_df_filtered)
    total_with_go <- length(unique(go_data$TAIR))

    go_counts$pvalue <- sapply(seq_len(nrow(go_counts)), function(i) {
      phyper(go_counts$Gene_Count[i] - 1,
             total_with_go,
             total_genes - total_with_go,
             length(gene_list),
             lower.tail = FALSE)
    })

    # Benjamini-Hochberg correction
    go_counts$padj_BH   <- p.adjust(go_counts$pvalue, method = "BH")
    go_counts$log10padj <- -log10(go_counts$padj_BH + 1e-300)
    go_counts$GeneRatio <- go_counts$Gene_Count / length(gene_list)
    go_counts$Category  <- category_name

    go_counts <- go_counts[order(go_counts$padj_BH), ]
    head(go_counts, 20)

  }, error = function(e) {
    warning("GO BP error for ", category_name, ": ", e$message); NULL
  })
}

up_go   <- get_go_BP(up_genes,   "Upregulated")
down_go <- get_go_BP(down_genes, "Downregulated")

# Reusable GO dot plot function
go_dotplot <- function(go_df, title_str, hi_color) {
  ggplot(go_df,
         aes(x = GeneRatio, y = reorder(TERM, GeneRatio),
             size = Gene_Count, color = log10padj)) +
    geom_point(alpha = 0.85) +
    scale_color_gradient(low = COL_LIGHT, high = hi_color,
                         name = "-log10(BH padj)") +
    scale_size_continuous(range = c(3, 12), name = "Gene Count") +
    labs(title    = title_str,
         subtitle = "Biological Process | BH correction | Top 20",
         x = "Gene Ratio", y = "GO Term") +
    theme_scheme(base_size = 11) +
    theme(axis.text.y = element_text(size = 9))
}

if (!is.null(up_go) && nrow(up_go) > 0) {

  # GO Plot 1: Upregulated Dot Plot
  pdf("GO_BP_upregulated_dotplot.pdf", width = 14, height = 10)
  print(go_dotplot(up_go,
                   "GO Enrichment (BP): Upregulated Genes — Drought Stress",
                   COL_UP))
  dev.off()
  cat("✅ 1/6 GO_BP_upregulated_dotplot.pdf saved\n")

  # GO Plot 2: Upregulated Gradient
  pdf("GO_BP_upregulated_gradient.pdf", width = 14, height = 10)
  print(ggplot(up_go, aes(x = Gene_Count, y = reorder(TERM, Gene_Count),
                          size = Gene_Count, color = log10padj)) +
    geom_point(alpha = 0.85) +
    scale_color_gradient(low = COL_LIGHT, high = COL_UP,
                         name = "-log10(BH padj)") +
    scale_size_continuous(range = c(3, 12), name = "Gene Count") +
    labs(title    = "GO Enrichment (BP): Upregulated Genes — Drought Stress",
         subtitle = "Gradient Bubble | BH correction | Top 20",
         x = "Number of Genes", y = "GO Term") +
    theme_scheme(base_size = 11) +
    theme(axis.text.y = element_text(size = 9)))
  dev.off()
  cat("✅ 2/6 GO_BP_upregulated_gradient.pdf saved\n")

  # GO Plot 3: Upregulated Lollipop
  pdf("GO_BP_upregulated_lollipop.pdf", width = 14, height = 10)
  print(ggplot(up_go, aes(x = Gene_Count, y = reorder(TERM, Gene_Count))) +
    geom_segment(aes(x = 0, xend = Gene_Count,
                     y = reorder(TERM, Gene_Count),
                     yend = reorder(TERM, Gene_Count)),
                 color = COL_LIGHT, linewidth = 0.6) +
    geom_point(aes(size = Gene_Count, color = log10padj), alpha = 0.9) +
    scale_color_gradient2(low = COL_MID, mid = COL_PALE, high = COL_UP,
                          midpoint = median(up_go$log10padj),
                          name = "-log10(BH padj)") +
    scale_size_continuous(range = c(4, 14), name = "Gene Count") +
    labs(title    = "GO Enrichment (BP): Upregulated Genes — Drought Stress",
         subtitle = "Lollipop | BH correction | Top 20",
         x = "Number of Genes", y = "") +
    theme_scheme(base_size = 11) +
    theme(axis.text.y        = element_text(size = 9),
          panel.grid.major.y = element_blank()))
  dev.off()
  cat("✅ 3/6 GO_BP_upregulated_lollipop.pdf saved\n")

  write.csv(up_go, "GO_BP_upregulated_top20.csv", row.names = FALSE)
}

if (!is.null(down_go) && nrow(down_go) > 0) {

  # GO Plot 4: Downregulated Dot Plot
  pdf("GO_BP_downregulated_dotplot.pdf", width = 14, height = 10)
  print(go_dotplot(down_go,
                   "GO Enrichment (BP): Downregulated Genes — Drought Stress",
                   COL_DOWN))
  dev.off()
  cat("✅ 4/6 GO_BP_downregulated_dotplot.pdf saved\n")

  # GO Plot 5: Downregulated Gradient
  pdf("GO_BP_downregulated_gradient.pdf", width = 14, height = 10)
  print(ggplot(down_go, aes(x = Gene_Count, y = reorder(TERM, Gene_Count),
                            size = Gene_Count, color = log10padj)) +
    geom_point(alpha = 0.85) +
    scale_color_gradient(low = COL_LIGHT, high = COL_DOWN,
                         name = "-log10(BH padj)") +
    scale_size_continuous(range = c(3, 12), name = "Gene Count") +
    labs(title    = "GO Enrichment (BP): Downregulated Genes — Drought Stress",
         subtitle = "Gradient Bubble | BH correction | Top 20",
         x = "Number of Genes", y = "GO Term") +
    theme_scheme(base_size = 11) +
    theme(axis.text.y = element_text(size = 9)))
  dev.off()
  cat("✅ 5/6 GO_BP_downregulated_gradient.pdf saved\n")

  # GO Plot 6: Downregulated Lollipop
  pdf("GO_BP_downregulated_lollipop.pdf", width = 14, height = 10)
  print(ggplot(down_go, aes(x = Gene_Count, y = reorder(TERM, Gene_Count))) +
    geom_segment(aes(x = 0, xend = Gene_Count,
                     y = reorder(TERM, Gene_Count),
                     yend = reorder(TERM, Gene_Count)),
                 color = COL_LIGHT, linewidth = 0.6) +
    geom_point(aes(size = Gene_Count, color = log10padj), alpha = 0.9) +
    scale_color_gradient2(low = COL_MID, mid = COL_PALE, high = COL_DOWN,
                          midpoint = median(down_go$log10padj),
                          name = "-log10(BH padj)") +
    scale_size_continuous(range = c(4, 14), name = "Gene Count") +
    labs(title    = "GO Enrichment (BP): Downregulated Genes — Drought Stress",
         subtitle = "Lollipop | BH correction | Top 20",
         x = "Number of Genes", y = "") +
    theme_scheme(base_size = 11) +
    theme(axis.text.y        = element_text(size = 9),
          panel.grid.major.y = element_blank()))
  dev.off()
  cat("✅ 6/6 GO_BP_downregulated_lollipop.pdf saved\n\n")

  write.csv(down_go, "GO_BP_downregulated_top20.csv", row.names = FALSE)
}

# ==============================================================================
# PART 12: KEGG PATHWAY ANALYSIS + BH (3 plots)
# ==============================================================================

cat("========================================\n")
cat("PART 12: KEGG Pathway Analysis (BH)\n")
cat("========================================\n\n")

get_kegg_BH <- function(gene_list, direction) {
  if (length(gene_list) == 0) { warning("No genes for KEGG ", direction); return(NULL) }
  tryCatch({
    kegg_data <- select(org.At.tair.db,
                        keys    = as.character(gene_list),
                        columns = "PATH",
                        keytype = "TAIR")
    kegg_data <- kegg_data[!is.na(kegg_data$PATH), ]
    if (nrow(kegg_data) == 0) { warning("No KEGG for ", direction); return(NULL) }

    kegg_counts           <- as.data.frame(table(kegg_data$PATH))
    colnames(kegg_counts) <- c("KEGG_Pathway", "Gene_Count")

    total_genes     <- nrow(res_df_filtered)
    total_with_path <- length(unique(kegg_data$TAIR))

    kegg_counts$pvalue <- sapply(seq_len(nrow(kegg_counts)), function(i) {
      phyper(kegg_counts$Gene_Count[i] - 1,
             total_with_path,
             total_genes - total_with_path,
             length(gene_list),
             lower.tail = FALSE)
    })

    # Benjamini-Hochberg correction
    kegg_counts$padj_BH   <- p.adjust(kegg_counts$pvalue, method = "BH")
    kegg_counts$log10padj <- -log10(kegg_counts$padj_BH + 1e-300)
    kegg_counts$GeneRatio <- kegg_counts$Gene_Count / length(gene_list)
    kegg_counts$Direction <- direction

    kegg_counts <- kegg_counts[order(kegg_counts$padj_BH), ]
    head(kegg_counts, 15)

  }, error = function(e) {
    warning("KEGG error for ", direction, ": ", e$message); NULL
  })
}

up_kegg   <- get_kegg_BH(up_genes,   "Upregulated")
down_kegg <- get_kegg_BH(down_genes, "Downregulated")

if (!is.null(up_kegg) && !is.null(down_kegg)) {

  combined_kegg <- rbind(up_kegg, down_kegg)

  # KEGG Plot 1: Dot Plot Faceted
  pdf("KEGG_dotplot_combined.pdf", width = 14, height = 10)
  print(ggplot(combined_kegg,
               aes(x = GeneRatio, y = reorder(KEGG_Pathway, GeneRatio),
                   size = Gene_Count, color = log10padj)) +
    geom_point(alpha = 0.85) +
    scale_color_gradient(low = COL_LIGHT, high = COL_DARK,
                         name = "-log10(BH padj)") +
    scale_size_continuous(range = c(3, 12), name = "Gene Count") +
    facet_wrap(~ Direction, scales = "free_y", ncol = 2) +
    labs(title    = "KEGG Pathways: Drought Stress",
         subtitle = "BH correction | Top 15 per direction",
         x = "Gene Ratio", y = "KEGG Pathway") +
    theme_scheme(base_size = 11) +
    theme(axis.text.y = element_text(size = 9)))
  dev.off()
  cat("✅ 1/3 KEGG_dotplot_combined.pdf saved\n")

  # KEGG Plot 2: Bubble Combined
  pdf("KEGG_bubble_combined.pdf", width = 14, height = 10)
  print(ggplot(combined_kegg,
               aes(x = Gene_Count, y = reorder(KEGG_Pathway, Gene_Count),
                   color = Direction, size = Gene_Count)) +
    geom_point(alpha = 0.75) +
    scale_color_manual(values = c("Upregulated"   = COL_DROUGHT,
                                  "Downregulated" = COL_DOWN)) +
    scale_size_continuous(range = c(3, 10)) +
    labs(title    = "KEGG Pathways: Upregulated vs Downregulated",
         subtitle = "BH correction | Top 15 per direction",
         x = "Gene Count", y = "KEGG Pathway",
         color = "Direction", size = "Gene Count") +
    theme_scheme(base_size = 11) +
    theme(axis.text.y = element_text(size = 9)))
  dev.off()
  cat("✅ 2/3 KEGG_bubble_combined.pdf saved\n")

  # KEGG Plot 3: Bubble Faceted
  pdf("KEGG_bubble_faceted.pdf", width = 16, height = 10)
  print(ggplot(combined_kegg,
               aes(x = Gene_Count, y = reorder(KEGG_Pathway, Gene_Count),
                   color = Direction, size = Gene_Count)) +
    geom_point(alpha = 0.85) +
    scale_color_manual(values = c("Upregulated"   = COL_DROUGHT,
                                  "Downregulated" = COL_DOWN)) +
    scale_size_continuous(range = c(4, 12)) +
    facet_wrap(~ Direction, scales = "free_y", ncol = 2) +
    labs(title    = "KEGG Pathways: Side by Side",
         subtitle = "BH correction | Top 15 per direction",
         x = "Gene Count", y = "KEGG Pathway") +
    theme_scheme(base_size = 11) +
    theme(axis.text.y     = element_text(size = 9),
          legend.position = "bottom"))
  dev.off()
  cat("✅ 3/3 KEGG_bubble_faceted.pdf saved\n\n")

  write.csv(combined_kegg, "KEGG_pathways_BH_combined.csv", row.names = FALSE)
}

# ==============================================================================
# PART 13: SAVE ALL RESULTS
# ==============================================================================

cat("========================================\n")
cat("PART 13: Saving Results\n")
cat("========================================\n\n")

saveRDS(res, "deseq2_results.rds")
saveRDS(dds, "dds_object.rds")
saveRDS(vsd, "vsd_object.rds")
write.csv(as.data.frame(res), "deseq2_results_full.csv", row.names = TRUE)

cat("✅ Saved:\n")
cat("  - deseq2_results.rds\n")
cat("  - dds_object.rds\n")
cat("  - vsd_object.rds\n")
cat("  - deseq2_results_full.csv\n\n")

# ==============================================================================
# ANALYSIS COMPLETE
# ==============================================================================

cat("========================================\n")
cat("🎉 DROUGHT ANALYSIS COMPLETE!\n")
cat("========================================\n\n")

cat("Significant genes (padj<0.05, |log2FC|>1):",
    sum(res_df$significant == "Significant"), "\n")
cat("  Upregulated:", length(up_genes), "\n")
cat("  Downregulated:", length(down_genes), "\n\n")

cat("ALL PLOTS SAVED:\n")
cat("  QC (7 plots):\n")
cat("    ✅ QC_gene_count_distribution.pdf\n")
cat("    ✅ QC_PCA_v1_standard.pdf\n")
cat("    ✅ QC_PCA_v2_labeled.pdf\n")
cat("    ✅ QC_sample_distance_heatmap.pdf\n")
cat("    ✅ QC_dispersion_plot.pdf\n")
cat("    ✅ QC_pvalue_distribution.pdf\n")
cat("    ✅ QC_independent_filtering.pdf\n")
cat("  Results (4 plots):\n")
cat("    ✅ RESULTS_volcano_plot.pdf\n")
cat("    ✅ RESULTS_volcano_enhanced_labeled.pdf\n")
cat("    ✅ RESULTS_MA_plot.pdf\n")
cat("    ✅ RESULTS_heatmap_top50.pdf\n")
cat("  Expression (1 plot):\n")
cat("    ✅ EXPRESSION_massive_heatmap_clustered.pdf\n")
cat("  GO BP + BH (6 plots):\n")
cat("    ✅ GO_BP_upregulated_dotplot.pdf\n")
cat("    ✅ GO_BP_upregulated_gradient.pdf\n")
cat("    ✅ GO_BP_upregulated_lollipop.pdf\n")
cat("    ✅ GO_BP_downregulated_dotplot.pdf\n")
cat("    ✅ GO_BP_downregulated_gradient.pdf\n")
cat("    ✅ GO_BP_downregulated_lollipop.pdf\n")
cat("  KEGG + BH (3 plots):\n")
cat("    ✅ KEGG_dotplot_combined.pdf\n")
cat("    ✅ KEGG_bubble_combined.pdf\n")
cat("    ✅ KEGG_bubble_faceted.pdf\n")
cat("\n  TOTAL: 21 plots\n")
cat("\nReady for comparison with salt data!\n")
