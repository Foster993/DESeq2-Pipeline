################################################################################
# SCRIPT 3 — Arabidopsis Drought vs Salt Comparative Analysis
# All original plots kept | 4 directional Venns restored
# KEGG y-axis shows pathway names not IDs (via KEGGREST)
# GO y-axis shows term names (always was — no change needed)
# All fixes retained: correct phyper | universe consistency | PNG | no dead code
################################################################################
# ==============================================================================
# PART 1: SETUP
# ==============================================================================
cat("========================================\n")
cat("PART 1: Setup\n")
cat("========================================\n\n")
cran_pkgs <- c("ggplot2","ggvenn","ggrepel","pheatmap",
               "gridExtra","scales","RColorBrewer","UpSetR","dplyr")
for (pkg in cran_pkgs) {
  if (!pkg %in% rownames(installed.packages())) {
    message("Installing: ", pkg); install.packages(pkg)
  } else message("OK ", pkg)
}
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
bioc_pkgs <- c("DESeq2","AnnotationDbi","org.At.tair.db","GO.db","KEGGREST")
for (pkg in bioc_pkgs) {
  if (!pkg %in% rownames(installed.packages())) {
    message("Installing: ", pkg); BiocManager::install(pkg)
  } else message("OK ", pkg)
}
library(ggplot2); library(ggvenn);    library(ggrepel);  library(UpSetR)
library(pheatmap); library(gridExtra); library(scales);   library(dplyr)
library(RColorBrewer); library(DESeq2); library(KEGGREST)
library(AnnotationDbi); library(org.At.tair.db); library(GO.db)
select <- dplyr::select
cat("\nAll libraries loaded\n\n")

# ==============================================================================
# COLORS — *** ONLY SECTION CHANGED ***
# grey        = control / not significant
# red         = drought upregulated
# green       = drought downregulated
# teal        = salt upregulated
# pink        = salt downregulated
# ==============================================================================
COL_DROUGHT_UP   <- "#CC2200"   # red   — drought upregulated
COL_DROUGHT_DOWN <- "#33AA44"   # green — drought downregulated
COL_SALT_UP      <- "#008B8B"   # teal  — salt upregulated
COL_SALT_DOWN    <- "#FF69B4"   # pink  — salt downregulated
COL_DROUGHT      <- "#CC2200"   # red   — drought identity color (Venns, bars)
COL_SALT         <- "#008B8B"   # teal  — salt identity color (Venns, bars)
COL_CONTROL      <- "#999999"   # grey  — control / not significant
COL_DARK         <- "#1a1a2e"
COL_UP           <- "#CC2200"   # red   (drought UP default)
COL_DOWN         <- "#33AA44"   # green (drought DOWN default)
COL_OVERLAP      <- "#884488"   # purple — shared/overlap
COL_GRAY1        <- "#666666"
COL_GRAY2        <- "#999999"
COL_PALE         <- "#f5f5f5"

theme_compare <- function(base_size=12) {
  theme_bw(base_size=base_size) +
    theme(
      plot.title    = element_text(hjust=0.5, face="bold",
                                   size=base_size+3, color=COL_DARK),
      plot.subtitle = element_text(hjust=0.5, size=base_size, color=COL_GRAY1),
      axis.title    = element_text(size=base_size+1, color=COL_DARK),
      axis.text     = element_text(size=base_size-1, color=COL_GRAY1),
      panel.border  = element_rect(color=COL_DARK, linewidth=0.8),
      legend.title  = element_text(size=base_size, color=COL_DARK),
      legend.text   = element_text(size=base_size-1, color=COL_GRAY1),
      strip.background = element_rect(fill=COL_DARK),
      strip.text    = element_text(color="white", face="bold", size=base_size)
    )
}

# ==============================================================================
# PART 2: LOAD RESULTS
# ==============================================================================
cat("========================================\n")
cat("PART 2: Load Results\n")
cat("========================================\n\n")
base_dir    <- "~/Downloads/deseq2_arabidopsis"
drought_dir <- file.path(base_dir, "arabidopsis_drought")
salt_dir    <- file.path(base_dir, "arabidopsis_salt")
out_dir     <- file.path(base_dir, "Comparative_Analysis")
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)
drought_res <- readRDS(file.path(drought_dir, "deseq2_results.rds"))
salt_res    <- readRDS(file.path(salt_dir,    "deseq2_results.rds"))
drought_vsd <- readRDS(file.path(drought_dir, "vsd_object.rds"))
salt_vsd    <- readRDS(file.path(salt_dir,    "vsd_object.rds"))
cat(sprintf("Drought: %d genes | Salt: %d genes\n\n",
            nrow(drought_res), nrow(salt_res)))
drought_df           <- as.data.frame(drought_res)
drought_df$gene      <- rownames(drought_df)
drought_df$condition <- "Drought"
salt_df           <- as.data.frame(salt_res)
salt_df$gene      <- rownames(salt_df)
salt_df$condition <- "Salt"
shared_universe <- intersect(
  rownames(drought_df[!is.na(drought_df$padj), ]),
  rownames(salt_df[!is.na(salt_df$padj), ])
)
N <- length(shared_universe)
cat(sprintf("Shared universe (tested in both): %d genes\n\n", N))
drought_up   <- rownames(subset(drought_df[shared_universe, ],
                                padj < 0.05 & log2FoldChange >  1))
drought_down <- rownames(subset(drought_df[shared_universe, ],
                                padj < 0.05 & log2FoldChange < -1))
salt_up      <- rownames(subset(salt_df[shared_universe, ],
                                padj < 0.05 & log2FoldChange >  1))
salt_down    <- rownames(subset(salt_df[shared_universe, ],
                                padj < 0.05 & log2FoldChange < -1))
cat(sprintf("Drought UP: %d | DOWN: %d\n",   length(drought_up),   length(drought_down)))
cat(sprintf("Salt    UP: %d | DOWN: %d\n\n", length(salt_up),      length(salt_down)))
get_symbols <- function(gene_ids) {
  tryCatch({
    s <- AnnotationDbi::mapIds(org.At.tair.db, keys=gene_ids,
                               column="SYMBOL", keytype="TAIR", multiVals="first")
    s[is.na(s)] <- gene_ids[is.na(s)]; s
  }, error=function(e) setNames(gene_ids, gene_ids))
}

# ==============================================================================
# BUILD KEGG PATHWAY NAME LOOKUP
# ==============================================================================
cat("Building KEGG pathway name lookup...\n")
kegg_name_lookup <- tryCatch({
  pw   <- KEGGREST::keggList("pathway", "ath")
  ids  <- sub("path:ath", "", names(pw))
  nms  <- sub(" - Arabidopsis thaliana.*", "", as.character(pw))
  lkp  <- setNames(nms, ids)
  cat(sprintf("  KEGG names loaded: %d pathways\n\n", length(lkp)))
  lkp
}, error=function(e) {
  cat(sprintf("  KEGGREST unavailable (%s) — will use pathway IDs\n\n", e$message))
  NULL
})
kegg_label <- function(ids) {
  if (is.null(kegg_name_lookup)) return(as.character(ids))
  nms <- kegg_name_lookup[as.character(ids)]
  nms[is.na(nms)] <- as.character(ids)[is.na(nms)]
  nms
}

# ==============================================================================
# PART 3: SIMPLE VENN — All drought DEGs vs all salt DEGs
# Drought circle = red (COL_DROUGHT), Salt circle = teal (COL_SALT)
# ==============================================================================
cat("========================================\n")
cat("PART 3: Simple Venn Diagram\n")
cat("========================================\n\n")
all_drought_degs <- union(drought_up, drought_down)
all_salt_degs    <- union(salt_up,    salt_down)
shared_any       <- intersect(all_drought_degs, all_salt_degs)
drought_only     <- setdiff(all_drought_degs, all_salt_degs)
salt_only        <- setdiff(all_salt_degs,    all_drought_degs)
cat(sprintf("  Drought only: %d | Shared: %d | Salt only: %d\n\n",
            length(drought_only), length(shared_any), length(salt_only)))
p_simple <- ggvenn(
  list(Drought=all_drought_degs, Salt=all_salt_degs),
  fill_color    = c(COL_DROUGHT, COL_SALT),
  fill_alpha    = 0.45,
  stroke_size   = 1.2,
  set_name_size = 5,
  text_size     = 5
) +
  labs(title    = "All Drought DEGs vs All Salt DEGs",
       subtitle = sprintf("Drought only: %d  |  Shared: %d  |  Salt only: %d",
                          length(drought_only), length(shared_any),
                          length(salt_only))) +
  theme_void() +
  theme(
    plot.title      = element_text(hjust=0.5, face="bold", size=16,
                                   color=COL_DARK, margin=margin(b=6)),
    plot.subtitle   = element_text(hjust=0.5, size=12,
                                   color=COL_GRAY1, margin=margin(b=10)),
    plot.background = element_rect(fill="white", color=NA),
    plot.margin     = margin(20,20,20,20)
  )
ggsave(file.path(out_dir, "Venn_simple_drought_vs_salt.png"),
       p_simple, width=8, height=7, dpi=300)
cat("Venn_simple_drought_vs_salt.png saved\n\n")
write.csv(data.frame(gene=shared_any, symbol=get_symbols(shared_any)),
          file.path(out_dir, "genes_shared_any_direction.csv"), row.names=FALSE)

# ==============================================================================
# PART 4: FOUR DIRECTIONAL VENN DIAGRAMS
# Each Venn uses the DIRECTIONAL colors for its two circles:
#   Venn1 Down/Down  : green  (drought DOWN) vs pink (salt DOWN)
#   Venn2 Down/Up    : green  (drought DOWN) vs teal (salt UP)
#   Venn3 Up/Down    : red    (drought UP)   vs pink (salt DOWN)
#   Venn4 Up/Up      : red    (drought UP)   vs teal (salt UP)
# ==============================================================================
cat("========================================\n")
cat("PART 4: Four Directional Venn Diagrams\n")
cat("========================================\n\n")
make_venn <- function(list_a, list_b,
                      label_a, label_b,
                      title_str, subtitle_str,
                      filename,
                      col_a, col_b) {          # *** directional colors passed in
  venn_list        <- list(list_a, list_b)
  names(venn_list) <- c(label_a, label_b)
  
  overlap <- intersect(list_a, list_b)
  only_a  <- setdiff(list_a, list_b)
  only_b  <- setdiff(list_b, list_a)
  
  cat(sprintf("%s\n  Only %s: %d | Shared: %d | Only %s: %d\n\n",
              title_str, label_a, length(only_a),
              length(overlap), label_b, length(only_b)))
  
  p <- ggvenn(venn_list,
              fill_color    = c(col_a, col_b),   # *** directional colors
              fill_alpha    = 0.45,
              stroke_size   = 1.2,
              set_name_size = 5,
              text_size     = 5) +
    labs(title=title_str, subtitle=subtitle_str) +
    theme_void() +
    theme(
      plot.title      = element_text(hjust=0.5, face="bold", size=16,
                                     color=COL_DARK, margin=margin(b=6)),
      plot.subtitle   = element_text(hjust=0.5, size=12,
                                     color=COL_GRAY1, margin=margin(b=10)),
      plot.background = element_rect(fill="white", color=NA),
      plot.margin     = margin(20,20,20,20)
    )
  
  ggsave(file.path(out_dir, filename), p, width=8, height=7, dpi=300)
  cat(sprintf("  %s saved\n\n", filename))
  
  list(overlap=overlap, only_a=only_a, only_b=only_b)
}

# Venn 1: DOWN in both — green vs pink
venn1 <- make_venn(
  drought_down, salt_down,
  "Drought DOWN", "Salt DOWN",
  "Venn 1: Downregulated in Both Stresses",
  "Genes repressed by drought AND salt",
  "Venn1_Down_Down.png",
  COL_DROUGHT_DOWN, COL_SALT_DOWN        # green, pink
)
# Venn 2: Drought DOWN, Salt UP — green vs teal
venn2 <- make_venn(
  drought_down, salt_up,
  "Drought DOWN", "Salt UP",
  "Venn 2: Drought DOWN vs Salt UP",
  "Genes repressed by drought but induced by salt",
  "Venn2_Down_Up.png",
  COL_DROUGHT_DOWN, COL_SALT_UP          # green, teal
)
# Venn 3: Drought UP, Salt DOWN — red vs pink
venn3 <- make_venn(
  drought_up, salt_down,
  "Drought UP", "Salt DOWN",
  "Venn 3: Drought UP vs Salt DOWN",
  "Genes induced by drought but repressed by salt",
  "Venn3_Up_Down.png",
  COL_DROUGHT_UP, COL_SALT_DOWN          # red, pink
)
# Venn 4: UP in both — red vs teal
venn4 <- make_venn(
  drought_up, salt_up,
  "Drought UP", "Salt UP",
  "Venn 4: Upregulated in Both Stresses",
  "Genes induced by drought AND salt",
  "Venn4_Up_Up.png",
  COL_DROUGHT_UP, COL_SALT_UP            # red, teal
)
write.csv(data.frame(gene=venn1$overlap),
          file.path(out_dir, "genes_down_down.csv"), row.names=FALSE)
write.csv(data.frame(gene=venn2$overlap),
          file.path(out_dir, "genes_down_up.csv"),   row.names=FALSE)
write.csv(data.frame(gene=venn3$overlap),
          file.path(out_dir, "genes_up_down.csv"),   row.names=FALSE)
write.csv(data.frame(gene=venn4$overlap),
          file.path(out_dir, "genes_up_up.csv"),     row.names=FALSE)
cat("All 4 directional Venn diagrams + gene list CSVs saved\n\n")

# ==============================================================================
# PART 5: UPSET PLOT
# ==============================================================================
cat("========================================\n")
cat("PART 5: UpSet Plot\n")
cat("========================================\n\n")
all_deg_ids <- unique(c(drought_up, drought_down, salt_up, salt_down))
upset_mat <- data.frame(
  gene         = all_deg_ids,
  Drought_UP   = as.integer(all_deg_ids %in% drought_up),
  Drought_DOWN = as.integer(all_deg_ids %in% drought_down),
  Salt_UP      = as.integer(all_deg_ids %in% salt_up),
  Salt_DOWN    = as.integer(all_deg_ids %in% salt_down),
  stringsAsFactors = FALSE
)
cat(sprintf("UpSet matrix: %d total DEGs\n\n", nrow(upset_mat)))
write.csv(upset_mat,
          file.path(out_dir, "upset_membership_matrix.csv"),
          row.names=FALSE)
png(file.path(out_dir, "UpSet_DEG_intersections.png"),
    width=14, height=8, units="in", res=300)
upset(
  upset_mat,
  sets           = c("Drought_UP","Drought_DOWN","Salt_UP","Salt_DOWN"),
  sets.bar.color = c(COL_DROUGHT_UP, COL_DROUGHT_DOWN,
                     COL_SALT_UP,    COL_SALT_DOWN),  # *** directional colors
  order.by       = "freq",
  decreasing     = TRUE,
  keep.order     = FALSE,
  mb.ratio       = c(0.6, 0.4),
  number.angles  = 0,
  point.size     = 3.5,
  line.size      = 1.2,
  text.scale     = c(1.5, 1.3, 1.2, 1.1, 1.4, 1.1),
  mainbar.y.label = "Intersection Size",
  sets.x.label    = "Set Size",
  main.bar.color  = COL_DARK,
  matrix.color    = COL_DARK,
  shade.color     = COL_PALE,
  shade.alpha     = 0.25
)
grid::grid.text(
  "DEG Set Intersections: Drought vs Salt",
  x=0.65, y=0.98,
  gp=grid::gpar(fontsize=15, fontface="bold", col=COL_DARK),
  just=c("centre","top")
)
grid::grid.text(
  "Ordered by intersection size | Shared universe only",
  x=0.65, y=0.94,
  gp=grid::gpar(fontsize=10, col=COL_GRAY1),
  just=c("centre","top")
)
dev.off()
cat("UpSet_DEG_intersections.png saved\n\n")

# ==============================================================================
# PART 6: FISHER'S EXACT TESTS
# ==============================================================================
cat("========================================\n")
cat("PART 6: Fisher's Exact Tests\n")
cat("========================================\n\n")
cat(sprintf("Universe (same as all Venns): %d\n\n", N))
run_fisher <- function(set_a, set_b, label) {
  in_a    <- intersect(set_a, shared_universe)
  in_b    <- intersect(set_b, shared_universe)
  both    <- length(intersect(in_a, in_b))
  only_a  <- length(in_a) - both
  only_b  <- length(in_b) - both
  neither <- N - both - only_a - only_b
  mat <- matrix(c(both, only_a, only_b, neither), nrow=2,
                dimnames=list(c("In_A","Not_A"), c("In_B","Not_B")))
  ft  <- fisher.test(mat, alternative="greater")
  cat(sprintf("%s: OR=%.3f  p=%.3e  %s\n",
              label, ft$estimate, ft$p.value,
              ifelse(ft$p.value<0.05,"*","")))
  data.frame(
    Comparison  = label,
    N_universe  = N,
    N_set_A     = length(in_a),
    N_set_B     = length(in_b),
    N_overlap   = both,
    Odds_Ratio  = round(ft$estimate, 4),
    pvalue      = ft$p.value,
    CI_low      = ft$conf.int[1],
    CI_high     = ft$conf.int[2],
    Significant = ifelse(ft$p.value < 0.05, "YES", "NO"),
    stringsAsFactors = FALSE
  )
}
fisher_results <- rbind(
  run_fisher(drought_down, salt_down, "Down-Down (both repressed)"),
  run_fisher(drought_down, salt_up,   "Down-Up (drought_down_salt_up)"),
  run_fisher(drought_up,   salt_down, "Up-Down (drought_up_salt_down)"),
  run_fisher(drought_up,   salt_up,   "Up-Up (both induced)")
)
fisher_results$padj_BH     <- p.adjust(fisher_results$pvalue, method="BH")
fisher_results$neg_log10_p <- -log10(fisher_results$pvalue)
write.csv(fisher_results,
          file.path(out_dir, "Fishers_exact_test_results.csv"),
          row.names=FALSE)
cat("\nFishers_exact_test_results.csv saved\n\n")
fisher_plot_df            <- fisher_results
fisher_plot_df$Comparison <- factor(fisher_plot_df$Comparison,
                                    levels=rev(fisher_plot_df$Comparison))
fisher_plot_df$label_sig  <- ifelse(fisher_plot_df$padj_BH<0.001,"***",
                                    ifelse(fisher_plot_df$padj_BH<0.01, "**",
                                           ifelse(fisher_plot_df$padj_BH<0.05,"*","ns")))
png(file.path(out_dir, "Fisher_forest_plot.png"),
    width=10, height=6, units="in", res=300)
print(ggplot(fisher_plot_df,
             aes(x=log2(Odds_Ratio), y=Comparison, color=Significant)) +
        geom_vline(xintercept=0, linetype="dashed",
                   color=COL_GRAY1, linewidth=0.8) +
        geom_errorbarh(aes(xmin=log2(CI_low+0.001), xmax=log2(CI_high)),
                       height=0.2, linewidth=0.8) +
        geom_point(size=5) +
        geom_text(aes(label=label_sig), hjust=-0.5, vjust=0.3,
                  size=6, color=COL_DARK) +
        scale_color_manual(values=c("YES"=COL_OVERLAP, "NO"=COL_GRAY2),
                           name="Significant (p<0.05)") +
        labs(title="Fisher's Exact: Overlap Significance",
             subtitle="All four directional comparisons | Shared universe",
             x="log2(Odds Ratio)", y="") +
        theme_compare())
dev.off()
cat("Fisher_forest_plot.png saved\n\n")

# ==============================================================================
# PART 7: FISHER 2x2 INTERSECTION HEATMAP
# ==============================================================================
cat("========================================\n")
cat("PART 7: Fisher 2x2 Intersection Heatmap\n")
cat("========================================\n\n")
overlap_mat <- matrix(
  c(length(intersect(drought_up,   salt_up)),
    length(intersect(drought_down, salt_up)),
    length(intersect(drought_up,   salt_down)),
    length(intersect(drought_down, salt_down))),
  nrow=2, ncol=2,
  dimnames=list(c("Drought UP","Drought DOWN"),
                c("Salt UP",   "Salt DOWN"))
)
get_cell_text <- function(drought_set, salt_set) {
  key_map <- c(
    "Drought UP.Salt UP"     = "Up-Up (both induced)",
    "Drought UP.Salt DOWN"   = "Up-Down (drought_up_salt_down)",
    "Drought DOWN.Salt UP"   = "Down-Up (drought_down_salt_up)",
    "Drought DOWN.Salt DOWN" = "Down-Down (both repressed)"
  )
  key     <- paste0(drought_set, ".", salt_set)
  comp    <- key_map[key]
  row_idx <- which(fisher_results$Comparison == comp)
  if (length(row_idx)==0)
    return(sprintf("%d\n?", overlap_mat[drought_set, salt_set]))
  p   <- fisher_results$padj_BH[row_idx]
  or  <- fisher_results$Odds_Ratio[row_idx]
  sig <- if (p<0.001) "***" else if (p<0.01) "**" else if (p<0.05) "*" else "ns"
  sprintf("%d\nOR=%.1f\n%s", overlap_mat[drought_set, salt_set], or, sig)
}
cell_texts <- matrix(
  c(get_cell_text("Drought UP",   "Salt UP"),
    get_cell_text("Drought DOWN", "Salt UP"),
    get_cell_text("Drought UP",   "Salt DOWN"),
    get_cell_text("Drought DOWN", "Salt DOWN")),
  nrow=2, ncol=2, dimnames=dimnames(overlap_mat)
)
log_overlap <- log10(overlap_mat + 1)
n_cols      <- 100
pal         <- colorRampPalette(c("#f5f5f5", COL_OVERLAP, "#330033"))(n_cols)
breaks      <- seq(0, max(log_overlap)*1.05, length.out=n_cols+1)
col_idx     <- cut(as.vector(log_overlap), breaks=breaks, labels=FALSE)
col_idx[is.na(col_idx)] <- 1L
png(file.path(out_dir, "Fisher_intersection_heatmap.png"),
    width=9, height=7, units="in", res=300)
par(mar=c(3,9,5,9), bg="white")
plot(NULL, xlim=c(0,2), ylim=c(0,2), axes=FALSE, xlab="", ylab="", asp=1)
for (dr in 1:2) {
  for (sa in 1:2) {
    flat_idx <- (sa-1)*2 + dr
    rect(sa-1, 2-dr, sa, 3-dr,
         col=pal[col_idx[flat_idx]], border="white", lwd=2)
    text(sa-0.5, 2.5-dr, labels=cell_texts[dr, sa],
         cex=1.3, col=COL_DARK, font=2)
  }
}
# *** directional colors for axis labels
mtext(c("Salt UP","Salt DOWN"), side=3, at=c(0.5,1.5),
      line=1.2, cex=1.2,
      col=c(COL_SALT_UP, COL_SALT_DOWN), font=2)
mtext(c("Drought UP","Drought DOWN"), side=2, at=c(1.5,0.5),
      line=1.2, cex=1.2,
      col=c(COL_DROUGHT_UP, COL_DROUGHT_DOWN), font=2, las=1)
title(main="Directional Overlap: Drought vs Salt DEGs",
      sub ="n = overlap | OR = Odds Ratio | * padj<0.05 | ** <0.01 | *** <0.001",
      cex.main=1.3, col.main=COL_DARK, cex.sub=0.85, col.sub=COL_GRAY1)
par(xpd=TRUE)
lx  <- 2.15; lh <- 2/n_cols
for (i in seq_len(n_cols))
  rect(lx, 2-i*lh, lx+0.12, 2-(i-1)*lh, col=pal[i], border=NA)
text(lx+0.12, 2.0, labels=max(overlap_mat), adj=c(0,0.5), cex=0.85, col=COL_DARK)
text(lx+0.12, 0.0, labels="0",              adj=c(0,0.5), cex=0.85, col=COL_DARK)
text(lx+0.06, 2.18, labels="n",             adj=c(0.5,0), cex=0.85, col=COL_DARK, font=3)
par(xpd=FALSE)
dev.off()
cat("Fisher_intersection_heatmap.png saved\n\n")

# ==============================================================================
# PART 8: SCATTER PLOT — log2FC comparison
# Not Significant = grey
# Drought Only    = red
# Salt Only       = teal
# Shared DEG      = purple (overlap)
# ==============================================================================
cat("========================================\n")
cat("PART 8: log2FC Scatter Plot\n")
cat("========================================\n\n")
scatter_df <- data.frame(
  gene         = shared_universe,
  drought_lfc  = drought_df[shared_universe, "log2FoldChange"],
  salt_lfc     = salt_df[shared_universe,    "log2FoldChange"],
  drought_padj = drought_df[shared_universe, "padj"],
  salt_padj    = salt_df[shared_universe,    "padj"],
  stringsAsFactors = FALSE
)
scatter_df <- scatter_df[
  !is.na(scatter_df$drought_lfc) & !is.na(scatter_df$salt_lfc), ]
scatter_df$category <- "Not Significant"
scatter_df$category[
  scatter_df$drought_padj<0.05 & abs(scatter_df$drought_lfc)>1 &
    !(scatter_df$salt_padj<0.05 & abs(scatter_df$salt_lfc)>1)
] <- "Drought Only"
scatter_df$category[
  scatter_df$salt_padj<0.05   & abs(scatter_df$salt_lfc)>1 &
    !(scatter_df$drought_padj<0.05 & abs(scatter_df$drought_lfc)>1)
] <- "Salt Only"
scatter_df$category[
  scatter_df$drought_padj<0.05 & abs(scatter_df$drought_lfc)>1 &
    scatter_df$salt_padj<0.05   & abs(scatter_df$salt_lfc)>1
] <- "Shared DEG"
cat("Scatter categories:\n"); print(table(scatter_df$category)); cat("\n")
shared_sig      <- scatter_df[scatter_df$category=="Shared DEG", ]
shared_sig$dist <- abs(shared_sig$drought_lfc) + abs(shared_sig$salt_lfc)
top_ids         <- head(shared_sig[order(-shared_sig$dist), "gene"], 20)
top_syms <- tryCatch(
  AnnotationDbi::mapIds(org.At.tair.db, keys=top_ids,
                        column="SYMBOL", keytype="TAIR", multiVals="first"),
  error=function(e) setNames(top_ids, top_ids))
top_syms[is.na(top_syms)] <- top_ids[is.na(top_syms)]
scatter_df$label <- ifelse(scatter_df$gene %in% top_ids,
                           top_syms[scatter_df$gene], "")
r_val <- round(cor(scatter_df$drought_lfc, scatter_df$salt_lfc,
                   use="complete.obs"), 3)
# *** directional color scheme for scatter categories
cat_colors <- c("Not Significant"=COL_GRAY2,
                "Drought Only"   =COL_DROUGHT,
                "Salt Only"      =COL_SALT,
                "Shared DEG"     =COL_OVERLAP)
ggsave(
  file.path(out_dir, "Scatter_log2FC_drought_vs_salt.png"),
  ggplot(scatter_df[order(scatter_df$category=="Not Significant",
                          decreasing=TRUE), ],
         aes(x=drought_lfc, y=salt_lfc, color=category)) +
    annotate("rect", xmin=1,    xmax=Inf,  ymin=1,    ymax=Inf,
             fill=COL_OVERLAP,      alpha=0.05) +
    annotate("rect", xmin=-Inf, xmax=-1,   ymin=-Inf, ymax=-1,
             fill=COL_DROUGHT_DOWN, alpha=0.05) +
    annotate("rect", xmin=1,    xmax=Inf,  ymin=-Inf, ymax=-1,
             fill=COL_SALT_DOWN,    alpha=0.07) +
    annotate("rect", xmin=-Inf, xmax=-1,   ymin=1,    ymax=Inf,
             fill=COL_DROUGHT_DOWN, alpha=0.04) +
    annotate("text", x=9,  y=9,  label="Both UP",        size=3.5,
             color=COL_OVERLAP,      fontface="italic", alpha=0.8) +
    annotate("text", x=-9, y=-9, label="Both DOWN",       size=3.5,
             color=COL_DROUGHT_DOWN, fontface="italic", alpha=0.8) +
    annotate("text", x=9,  y=-9, label="D\u2191 S\u2193", size=3.5,
             color=COL_SALT_DOWN,    fontface="italic", alpha=0.8) +
    annotate("text", x=-9, y=9,  label="D\u2193 S\u2191", size=3.5,
             color=COL_SALT_UP,      fontface="italic", alpha=0.8) +
    geom_point(alpha=0.45, size=1.1) +
    geom_text_repel(
      data=scatter_df[scatter_df$label!="", ],
      aes(label=label), size=2.8, color=COL_DARK,
      max.overlaps=30, box.padding=0.5,
      segment.color=COL_GRAY1, segment.size=0.3) +
    scale_color_manual(values=cat_colors, name="Category") +
    geom_hline(yintercept=0, linetype="dashed", color=COL_GRAY1, linewidth=0.5) +
    geom_vline(xintercept=0, linetype="dashed", color=COL_GRAY1, linewidth=0.5) +
    geom_hline(yintercept=c(-1,1), linetype="dotted",
               color=COL_SALT,    linewidth=0.4, alpha=0.6) +
    geom_vline(xintercept=c(-1,1), linetype="dotted",
               color=COL_DROUGHT, linewidth=0.4, alpha=0.6) +
    geom_abline(slope=1, intercept=0, linetype="dotted",
                color=COL_DARK, linewidth=0.4, alpha=0.5) +
    annotate("text", x=Inf, y=-Inf, hjust=1.05, vjust=-0.5,
             label=paste0("r = ", r_val),
             size=4.5, color=COL_DARK, fontface="bold") +
    coord_cartesian(xlim=c(-12,12), ylim=c(-12,12)) +
    labs(title    = "log2FC: Drought vs Salt (shrunken LFC)",
         subtitle = sprintf("Shared universe (n=%d) | Top %d shared DEGs labeled",
                            nrow(scatter_df), length(top_ids)),
         x="log2FC (Drought vs Control)",
         y="log2FC (Salt vs Control)") +
    theme_compare() +
    guides(color=guide_legend(override.aes=list(size=4, alpha=1))),
  width=11, height=10, dpi=300
)
cat("Scatter_log2FC_drought_vs_salt.png saved\n\n")

# ==============================================================================
# PART 9: SIDE-BY-SIDE VOLCANO PLOTS
# Drought UP=red, Drought DOWN=green, Salt UP=teal, Salt DOWN=pink, grey=not sig
# ==============================================================================
cat("========================================\n")
cat("PART 9: Side-by-Side Volcano Plots\n")
cat("========================================\n\n")
make_volcano_df <- function(df, condition_name) {
  df2 <- df[!is.na(df$padj) & df$gene %in% shared_universe, ]
  df2$status <- "Not Significant"
  df2$status[df2$padj<0.05 & df2$log2FoldChange>1]    <- "Upregulated"
  df2$status[df2$padj<0.05 & df2$log2FoldChange<(-1)] <- "Downregulated"
  df2$condition <- condition_name
  df2
}
drought_vol <- make_volcano_df(drought_df, "Drought")
salt_vol    <- make_volcano_df(salt_df,    "Salt")
both_vol    <- rbind(drought_vol, salt_vol)
both_vol$condition <- factor(both_vol$condition, levels=c("Drought","Salt"))

# *** assign per-point directional colors
both_vol$point_color <- COL_GRAY2
both_vol$point_color[
  both_vol$condition=="Drought" & both_vol$status=="Upregulated"
] <- COL_DROUGHT_UP
both_vol$point_color[
  both_vol$condition=="Drought" & both_vol$status=="Downregulated"
] <- COL_DROUGHT_DOWN
both_vol$point_color[
  both_vol$condition=="Salt" & both_vol$status=="Upregulated"
] <- COL_SALT_UP
both_vol$point_color[
  both_vol$condition=="Salt" & both_vol$status=="Downregulated"
] <- COL_SALT_DOWN

count_annot <- data.frame(
  condition = factor(c("Drought","Drought","Salt","Salt"),
                     levels=c("Drought","Salt")),
  n         = c(length(drought_up), length(drought_down),
                length(salt_up),    length(salt_down)),
  x         = c(8, -8, 8, -8),
  y         = c(18, 18, 18, 18),
  col       = c(COL_DROUGHT_UP, COL_DROUGHT_DOWN,
                COL_SALT_UP,    COL_SALT_DOWN)   # *** directional colors
)
png(file.path(out_dir, "Volcano_side_by_side.png"),
    width=16, height=8, units="in", res=300)
print(ggplot(both_vol, aes(x=log2FoldChange, y=-log10(padj))) +
        geom_point(aes(color=point_color), alpha=0.5, size=1.2) +
        scale_color_identity() +
        geom_vline(xintercept=c(-1,1), linetype="dashed",
                   color=COL_DARK, linewidth=0.5) +
        geom_hline(yintercept=-log10(0.05), linetype="dashed",
                   color=COL_DARK, linewidth=0.5) +
        geom_text(data=count_annot,
                  aes(x=x, y=y, label=paste0("n=",n), color=col),
                  size=4.5, fontface="bold", inherit.aes=FALSE) +
        scale_color_identity() +
        facet_wrap(~condition, ncol=2) +
        coord_cartesian(xlim=c(-12,12), ylim=c(0,20)) +
        labs(title    = "Volcano Plots: Drought vs Salt (shared universe)",
             subtitle = "padj<0.05 | |log2FC|>1 | shrunken LFC",
             x="log2 Fold Change", y="-log10(padj)") +
        theme_compare())
dev.off()
cat("Volcano_side_by_side.png saved\n\n")

# ==============================================================================
# PART 10: COMBINED EXPRESSION HEATMAP
# ==============================================================================
cat("========================================\n")
cat("PART 10: Combined Heatmap\n")
cat("========================================\n\n")
shared_all <- unique(c(venn1$overlap, venn2$overlap,
                       venn3$overlap, venn4$overlap))
cat(sprintf("Total directionally-shared DEGs: %d\n", length(shared_all)))
if (length(shared_all) >= 5) {
  d_mat <- assay(drought_vsd)[
    shared_all[shared_all %in% rownames(assay(drought_vsd))], ]
  s_mat <- assay(salt_vsd)[
    shared_all[shared_all %in% rownames(assay(salt_vsd))], ]
  common_rows <- intersect(rownames(d_mat), rownames(s_mat))
  d_mat <- d_mat[common_rows, ]
  s_mat <- s_mat[common_rows, ]
  colnames(d_mat) <- paste0("D_", colnames(d_mat))
  colnames(s_mat) <- paste0("S_", colnames(s_mat))
  comb_mat    <- cbind(d_mat, s_mat)
  comb_scaled <- t(scale(t(comb_mat)))
  
  d_col_data <- as.data.frame(colData(drought_vsd))
  s_col_data <- as.data.frame(colData(salt_vsd))
  
  col_ann_hm <- data.frame(
    Experiment = c(rep("Drought", ncol(d_mat)), rep("Salt", ncol(s_mat))),
    Condition  = c(as.character(d_col_data$condition),
                   as.character(s_col_data$condition)),
    row.names  = colnames(comb_mat)
  )
  row_ann_hm <- data.frame(
    Direction = ifelse(common_rows %in% venn4$overlap, "Both UP",
                       ifelse(common_rows %in% venn1$overlap, "Both DOWN",
                              ifelse(common_rows %in% venn3$overlap,
                                     "D_up S_down","D_down S_up"))),
    row.names = common_rows)
  
  # *** directional colors in heatmap annotations
  ann_col_hm <- list(
    Experiment = c("Drought"=COL_DROUGHT,      "Salt"=COL_SALT),
    Direction  = c("Both UP"    =COL_OVERLAP,
                   "Both DOWN"  =COL_DROUGHT_DOWN,
                   "D_up S_down"=COL_DROUGHT_UP,
                   "D_down S_up"=COL_SALT_UP))
  
  ph <- min(max(10, length(common_rows)*0.12), 60)
  tryCatch({
    png(file.path(out_dir, "Combined_heatmap_shared_DEGs.png"),
        width=14, height=ph, units="in", res=150)
    pheatmap(comb_scaled,
             cluster_rows=TRUE, cluster_cols=FALSE,
             show_rownames=length(common_rows)<=80,
             annotation_col=col_ann_hm, annotation_row=row_ann_hm,
             annotation_colors=ann_col_hm,
             color=colorRampPalette(c(COL_SALT_UP,"#eeeeee",COL_DROUGHT_UP))(100),
             main=sprintf("Shared DEGs Combined Heatmap (%d genes)",
                          length(common_rows)),
             gaps_col=ncol(d_mat), fontsize=9, fontsize_row=7,
             border_color=NA, cutree_rows=4, treeheight_row=40)
    dev.off()
    cat(sprintf("Combined_heatmap_shared_DEGs.png (%d genes)\n\n",
                length(common_rows)))
  }, error=function(e) {
    try(dev.off(), silent=TRUE)
    cat(sprintf("  SKIPPED Combined_heatmap: %s\n\n", e$message))
  })
} else {
  cat("Too few shared genes for heatmap\n\n")
}

# ==============================================================================
# PART 11: GO ENRICHMENT COMPARISON GRID
# ==============================================================================
cat("========================================\n")
cat("PART 11: GO Enrichment Comparison Grid\n")
cat("========================================\n\n")
cat("  Building shared universe GO map...\n")
comp_go_ann <- tryCatch(
  AnnotationDbi::select(org.At.tair.db, keys=shared_universe,
                        columns=c("GO","ONTOLOGY"), keytype="TAIR"),
  error=function(e) NULL)
k_total_comp <- NULL
if (!is.null(comp_go_ann)) {
  comp_go_bp   <- comp_go_ann[!is.na(comp_go_ann$GO) &
                                comp_go_ann$ONTOLOGY=="BP", ]
  k_total_comp <- table(comp_go_bp$GO)
  cat(sprintf("  Shared universe BP terms: %d\n\n", length(k_total_comp)))
}
get_go_BP_comparison <- function(gene_list, label) {
  if (length(gene_list) < 3) {
    cat(sprintf("  Too few genes for %s (%d)\n", label, length(gene_list)))
    return(NULL)
  }
  if (is.null(k_total_comp)) {
    cat("  Skipping: universe map not available\n"); return(NULL)
  }
  tryCatch({
    go_data <- AnnotationDbi::select(
      org.At.tair.db, keys=as.character(gene_list),
      columns=c("GO","ONTOLOGY"), keytype="TAIR")
    go_data <- go_data[!is.na(go_data$GO) & go_data$ONTOLOGY=="BP", ]
    if (nrow(go_data)==0) return(NULL)
    go_counts           <- as.data.frame(table(go_data$GO))
    colnames(go_counts) <- c("GO_ID","Gene_Count")
    go_term_names <- AnnotationDbi::select(
      GO.db, keys=as.character(go_counts$GO_ID),
      columns="TERM", keytype="GOID")
    go_counts <- merge(go_counts, go_term_names,
                       by.x="GO_ID", by.y="GOID", all.x=TRUE)
    go_counts$TERM[is.na(go_counts$TERM)] <- "Unknown BP term"
    go_counts$K_total <- as.integer(k_total_comp[as.character(go_counts$GO_ID)])
    go_counts <- go_counts[!is.na(go_counts$K_total) & go_counts$K_total >= 10, ]
    if (nrow(go_counts)==0) return(NULL)
    n_deg <- length(gene_list)
    go_counts$pvalue <- mapply(
      function(overlap, k_total)
        phyper(overlap-1, k_total, N-k_total, n_deg, lower.tail=FALSE),
      go_counts$Gene_Count, go_counts$K_total)
    go_counts$padj_BH   <- p.adjust(go_counts$pvalue, method="BH")
    go_counts$log10padj <- -log10(go_counts$padj_BH + 1e-300)
    go_counts$GeneRatio <- go_counts$Gene_Count / n_deg
    go_counts$Set       <- label
    go_counts           <- go_counts[order(go_counts$padj_BH), ]
    cat(sprintf("  %-22s: %d terms | %d sig\n",
                label, nrow(go_counts), sum(go_counts$padj_BH<0.05, na.rm=TRUE)))
    head(go_counts, 15)
  }, error=function(e) { cat("GO error:", e$message, "\n"); NULL })
}
go_up_up     <- get_go_BP_comparison(venn4$overlap, "Both Upregulated")
go_down_down <- get_go_BP_comparison(venn1$overlap, "Both Downregulated")
go_d_up_only <- get_go_BP_comparison(setdiff(drought_up, salt_up), "Drought Only UP")
go_s_up_only <- get_go_BP_comparison(setdiff(salt_up, drought_up), "Salt Only UP")
go_sets <- Filter(Negate(is.null),
                  list(go_up_up, go_down_down, go_d_up_only, go_s_up_only))
if (length(go_sets) >= 2) {
  go_combined <- do.call(rbind, go_sets)
  go_combined$TERM_short <- ifelse(
    nchar(go_combined$TERM)>45,
    paste0(substr(go_combined$TERM,1,42),"..."),
    go_combined$TERM)
  go_combined$Set <- factor(go_combined$Set,
                            levels=c("Both Upregulated","Both Downregulated",
                                     "Drought Only UP","Salt Only UP"))
  # *** directional strip colors via manual scale on Set
  set_colors <- c("Both Upregulated"  =COL_OVERLAP,
                  "Both Downregulated"=COL_DROUGHT_DOWN,
                  "Drought Only UP"   =COL_DROUGHT_UP,
                  "Salt Only UP"      =COL_SALT_UP)
  
  png(file.path(out_dir, "GO_comparison_grid_dotplot.png"),
      width=18, height=12, units="in", res=300)
  print(ggplot(go_combined,
               aes(x=GeneRatio, y=reorder(TERM_short, GeneRatio),
                   size=Gene_Count, color=log10padj)) +
          geom_point(alpha=0.85) +
          scale_color_gradient(low="#cceecc", high=COL_DARK,
                               name="-log10(BH padj)") +
          scale_size_continuous(range=c(2,10), name="Gene Count") +
          facet_wrap(~Set, scales="free_y", ncol=2) +
          labs(title    = "GO Enrichment Comparison: Drought vs Salt",
               subtitle = "Biological Process | BH | K per term from shared universe | K>=10",
               x="Gene Ratio", y="GO Term") +
          theme_compare(10) +
          theme(axis.text.y=element_text(size=8),
                panel.spacing=unit(1.2,"lines"),
                strip.background=element_rect(
                  fill=set_colors[levels(go_combined$Set)[1]]),
                strip.text=element_text(color="white",face="bold",size=10)))
  dev.off()
  cat("GO_comparison_grid_dotplot.png saved\n\n")
  write.csv(go_combined,
            file.path(out_dir, "GO_comparison_all_sets.csv"),
            row.names=FALSE)
} else {
  cat("Not enough GO sets — skipping grid\n\n")
}

# ==============================================================================
# PART 12: BAR CHARTS
# ==============================================================================
cat("========================================\n")
cat("PART 12: Bar Charts\n")
cat("========================================\n\n")
deg_summary <- data.frame(
  Stress    = factor(rep(c("Drought","Salt"), each=2), levels=c("Drought","Salt")),
  Direction = rep(c("Upregulated","Downregulated"), 2),
  Count     = c(length(drought_up), length(drought_down),
                length(salt_up),    length(salt_down))
)
deg_summary$Count_signed <- ifelse(deg_summary$Direction=="Downregulated",
                                   -deg_summary$Count, deg_summary$Count)
# *** per-bar directional fill
deg_summary$bar_fill <- c(COL_DROUGHT_UP, COL_DROUGHT_DOWN,
                          COL_SALT_UP,    COL_SALT_DOWN)

ggsave(file.path(out_dir, "Barchart_DEG_counts.png"),
       ggplot(deg_summary, aes(x=Stress, y=Count_signed, fill=bar_fill)) +
         geom_bar(stat="identity", position="identity", width=0.5,
                  color=COL_DARK, linewidth=0.4) +
         geom_hline(yintercept=0, color=COL_DARK, linewidth=0.8) +
         geom_text(aes(label=Count,
                       y=ifelse(Count_signed>0,
                                Count_signed+30,
                                Count_signed-30)),
                   size=5, fontface="bold", color=COL_DARK) +
         scale_fill_identity(guide=guide_legend(
           title="Direction",
           override.aes=list(
             fill=c(COL_DROUGHT_UP, COL_DROUGHT_DOWN,
                    COL_SALT_UP,    COL_SALT_DOWN)),
           labels=c("Drought UP","Drought DOWN","Salt UP","Salt DOWN"))) +
         scale_y_continuous(labels=function(x) abs(x), name="Number of DEGs") +
         annotate("text", x=0.55, y=max(deg_summary$Count)*0.7,
                  label="UP",   size=4, color=COL_DARK) +
         annotate("text", x=0.55, y=-max(deg_summary$Count)*0.7,
                  label="DOWN", size=4, color=COL_DARK) +
         labs(title="DEG Counts: Drought vs Salt (shared universe)",
              subtitle="padj<0.05 | |log2FC|>1 | shrunken LFC",
              x="") +
         theme_compare(),
       width=9, height=7, dpi=300)
cat("Barchart_DEG_counts.png saved\n")

overlap_summary <- data.frame(
  Category = c("Both UP","Both DOWN",
               "Drought UP\nSalt DOWN","Drought DOWN\nSalt UP",
               "Drought Only UP","Salt Only UP",
               "Drought Only DOWN","Salt Only DOWN"),
  Count    = c(length(venn4$overlap), length(venn1$overlap),
               length(venn3$overlap), length(venn2$overlap),
               length(venn3$only_a),  length(venn3$only_b),
               length(venn1$only_a),  length(venn1$only_b)),
  Type     = c("Shared","Shared","Discordant","Discordant",
               "Unique","Unique","Unique","Unique"),
  bar_fill = c(COL_OVERLAP,      COL_DROUGHT_DOWN,
               COL_DROUGHT_UP,   COL_SALT_UP,
               COL_DROUGHT_UP,   COL_SALT_UP,
               COL_DROUGHT_DOWN, COL_SALT_DOWN)
)
overlap_summary$Category <- factor(overlap_summary$Category,
                                   levels=overlap_summary$Category)
ggsave(file.path(out_dir, "Barchart_overlap_breakdown.png"),
       ggplot(overlap_summary, aes(x=Category, y=Count, fill=bar_fill)) +
         geom_bar(stat="identity", width=0.65,
                  color=COL_DARK, linewidth=0.4) +
         geom_text(aes(label=Count), vjust=-0.4, size=4.5,
                   fontface="bold", color=COL_DARK) +
         scale_fill_identity() +
         labs(title="DEG Overlap Breakdown (shared universe)",
              x="", y="Number of Genes") +
         theme_compare() +
         theme(axis.text.x=element_text(size=9, angle=15, hjust=1)),
       width=12, height=7, dpi=300)
cat("Barchart_overlap_breakdown.png saved\n")

stack_df <- data.frame(
  Stress   = factor(c(rep("Drought",3), rep("Salt",3)),
                    levels=c("Drought","Salt")),
  Category = rep(c("Shared (both)","Stress-specific","Discordant"), 2),
  Count    = c(
    length(venn4$overlap) + length(venn1$overlap),
    length(setdiff(drought_up,salt_up)) + length(setdiff(drought_down,salt_down)),
    length(venn2$overlap) + length(venn3$overlap),
    length(venn4$overlap) + length(venn1$overlap),
    length(setdiff(salt_up,drought_up)) + length(setdiff(salt_down,drought_down)),
    length(venn2$overlap) + length(venn3$overlap)
  )
)
stack_df$Category <- factor(stack_df$Category,
                            levels=c("Discordant","Stress-specific","Shared (both)"))
ggsave(file.path(out_dir, "Barchart_proportions_stacked.png"),
       ggplot(stack_df, aes(x=Stress, y=Count, fill=Category)) +
         geom_bar(stat="identity", width=0.5,
                  color=COL_DARK, linewidth=0.4) +
         geom_text(aes(label=Count),
                   position=position_stack(vjust=0.5),
                   size=4.5, fontface="bold", color="white") +
         scale_fill_manual(
           values=c("Shared (both)"   =COL_OVERLAP,
                    "Stress-specific" =COL_GRAY2,
                    "Discordant"      =COL_SALT_DOWN),
           name="") +
         labs(title="DEG Composition: Drought vs Salt",
              x="", y="Number of DEGs") +
         theme_compare() +
         theme(legend.position="bottom"),
       width=8, height=7, dpi=300)
cat("Barchart_proportions_stacked.png saved\n\n")

# ==============================================================================
# FINAL VERIFICATION
# ==============================================================================
cat("========================================\n")
cat("FINAL VERIFICATION\n")
cat("========================================\n\n")
expected_files <- c(
  "Venn_simple_drought_vs_salt.png",
  "genes_shared_any_direction.csv",
  "Venn1_Down_Down.png",
  "Venn2_Down_Up.png",
  "Venn3_Up_Down.png",
  "Venn4_Up_Up.png",
  "genes_down_down.csv","genes_down_up.csv",
  "genes_up_down.csv","genes_up_up.csv",
  "upset_membership_matrix.csv",
  "UpSet_DEG_intersections.png",
  "Fishers_exact_test_results.csv",
  "Fisher_forest_plot.png",
  "Fisher_intersection_heatmap.png",
  "Scatter_log2FC_drought_vs_salt.png",
  "Volcano_side_by_side.png",
  "Combined_heatmap_shared_DEGs.png",
  "GO_comparison_grid_dotplot.png",
  "GO_comparison_all_sets.csv",
  "Barchart_DEG_counts.png",
  "Barchart_overlap_breakdown.png",
  "Barchart_proportions_stacked.png"
)
all_ok <- TRUE
for (f in expected_files) {
  fp <- file.path(out_dir, f)
  if (file.exists(fp) && file.info(fp)$size > 500) {
    cat(sprintf("OK  %-45s %.1f KB\n", f, file.info(fp)$size/1024))
  } else if (file.exists(fp)) {
    cat(sprintf("SMALL  %s\n", f)); all_ok <- FALSE
  } else {
    cat(sprintf("MISSING  %s\n", f)); all_ok <- FALSE
  }
}
cat("\n")
if (all_ok) cat("ALL FILES VERIFIED\n\n") else cat("Some files missing\n\n")
cat("SUMMARY (shared universe only):\n")
cat(sprintf("  Drought: %d UP | %d DOWN\n", length(drought_up),   length(drought_down)))
cat(sprintf("  Salt:    %d UP | %d DOWN\n", length(salt_up),      length(salt_down)))
cat(sprintf("  Simple Venn — Drought only: %d | Shared: %d | Salt only: %d\n",
            length(drought_only), length(shared_any), length(salt_only)))
cat(sprintf("  Both UP: %d | Both DOWN: %d | D-up/S-dn: %d | D-dn/S-up: %d\n",
            length(venn4$overlap), length(venn1$overlap),
            length(venn3$overlap), length(venn2$overlap)))
cat("Fisher results:\n")
for (i in seq_len(nrow(fisher_results)))
  cat(sprintf("  %s -> %s (OR=%.2f, p=%.2e)\n",
              fisher_results$Comparison[i], fisher_results$Significant[i],
              fisher_results$Odds_Ratio[i], fisher_results$pvalue[i]))