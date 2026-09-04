#######################################

#__Date:__ 7-24-2026
#__Author:__ John Lemas
# Reworked for the ~ line + trt (dds_full) / ~ line + trt + line:trt (dds_int) design
# List of packages required
packages <- c("DESeq2",
              "tidyverse",
              "pheatmap",
              "RColorBrewer")

# if not installed, install and load it
for (package in packages) {
  if (!(package %in% rownames(installed.packages()))) {
    BiocManager::install(package)
  }
  library(package, character.only = TRUE)
}

#################################################

pdf(NULL)

# Set the working directory:
setwd("/Users/jmlemas/Library/CloudStorage/OneDrive-MichiganStateUniversity/Dissertation/Chapter3/RNA_seq/03_outputs/04_DESeq2/")

# ===========>------- Load saved objects --------<=========== #

dds_full <- readRDS("dds_full.rds")
dds_int  <- readRDS("dds_int.rds")  # not directly plotted here, but loaded for reference/consistency

sig_trt_genes    <- read.table("sig_trt_genes.txt",    header = TRUE, sep = "\t", row.names = 1)
sig_drt_genes    <- read.table("sig_drt_genes.txt",    header = TRUE, sep = "\t", row.names = 1)
sig_herb_genes   <- read.table("sig_herb_genes.txt",   header = TRUE, sep = "\t", row.names = 1)
sig_baseline_drt <- read.table("sig_baseline_drt.txt", header = TRUE, sep = "\t", row.names = 1)
sig_baseline_hrb <- read.table("sig_baseline_hrb.txt", header = TRUE, sep = "\t", row.names = 1)

# ===========>------- QC: dispersion estimates --------<=========== #

pdf(file = "disp_estimate.pdf", width = 8, height = 8, family = "Times")
plotDispEsts(dds_full)
dev.off()

plotDispEsts(dds_int)

# ===========>------- VST transform --------<=========== #
# NOTE: transformed into a separate object (vsd_full) instead of overwriting
# dds_full, so the original counts/dds object is still available later if needed.

vsd_full <- varianceStabilizingTransformation(dds_full, blind = TRUE)

# ===========>------- Sample-distance heatmap: all samples --------<=========== #

sampleDists <- dist(t(assay(vsd_full)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(vsd_full)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

p <- pheatmap(sampleDistMatrix,
              clustering_distance_rows = sampleDists,
              clustering_distance_cols = sampleDists,
              col = colors)

pdf("heatmap.pdf", height = 6, width = 7)
p
dev.off()

# ===========>------- PCA: all samples, all genes --------<=========== #

var_all <- prcomp(t(assay(vsd_full)))
perc_var <- (var_all$sdev)^2 / sum(var_all$sdev^2) * 100
x_label <- paste0("PC1 (", round(perc_var[1], 1), "%)")
y_label <- paste0("PC2 (", round(perc_var[2], 1), "%)")

pcaData <- plotPCA(vsd_full, intgroup = c("line", "trt"), returnData = TRUE)

plot_pca_full <- ggplot(pcaData, aes(x = PC1, y = PC2, color = line, shape = trt)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "B.hybridum Recurrent Selection PCA",
       x = x_label,
       y = y_label,
       color = "Line",
       shape = "Treatment"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold", family = "Times"),
    axis.title = element_text(size = 12, family = "Times"),
    legend.title = element_text(size = 10, family = "Times"),
    legend.text = element_text(size = 9, family = "Times")
  ) +
  scale_color_brewer(palette = "Set2")

pdf(file = "PCA_plot.pdf", width = 8, height = 8, family = "Times")
plot_pca_full
dev.off()

loadings <- var_all$rotation
write.csv(loadings, "pca_loadings.csv", row.names = TRUE)

# ===========>------- PCA: top 500 variable genes --------<=========== #

gene_vars <- rowVars(assay(vsd_full))
top_500_indices <- order(gene_vars, decreasing = TRUE)[1:500]
vst_matrix_top500 <- assay(vsd_full)[top_500_indices, ]

var_top500 <- prcomp(t(vst_matrix_top500))
perc_var_top500 <- (var_top500$sdev)^2 / sum(var_top500$sdev^2) * 100
x_label <- paste0("PC1 (", round(perc_var_top500[1], 1), "%)")
y_label <- paste0("PC2 (", round(perc_var_top500[2], 1), "%)")

pcaData_top500 <- as.data.frame(var_top500$x[, 1:2])
colnames(pcaData_top500) <- c("PC1", "PC2")
pcaData_top500$line <- vsd_full$line
pcaData_top500$trt  <- vsd_full$trt

plot_pca_top500 <- ggplot(pcaData_top500, aes(x = PC1, y = PC2, color = line, shape = trt)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "B.hybridum Recurrent Selection PCA (Top 500 Variable Genes)",
       x = x_label,
       y = y_label,
       color = "Line",
       shape = "Treatment"
  ) +
  scale_color_brewer(palette = "Set2")

pdf(file = "PCA_plot_top500.pdf", width = 8, height = 8, family = "Times")
plot_pca_top500
dev.off()

# ===========>------- PCA / heatmap: untreated samples only --------<=========== #
# Replaces the old dds_split object. Since we no longer fit a separate dds
# for the untreated subset, we subset the already-transformed full-dataset
# VST matrix to Control-treatment columns instead.

untrt_samples <- colnames(vsd_full)[vsd_full$trt == "Control"]
vst_untrt <- assay(vsd_full)[, untrt_samples]

sampleDists_untrt <- dist(t(vst_untrt))
sampleDistMatrix_untrt <- as.matrix(sampleDists_untrt)

p_untrt <- pheatmap(sampleDistMatrix_untrt,
                     clustering_distance_rows = sampleDists_untrt,
                     clustering_distance_cols = sampleDists_untrt,
                     col = colors)

pdf("heatmap_untrt.pdf", height = 6, width = 7)
p_untrt
dev.off()

var_untrt <- prcomp(t(vst_untrt))
perc_var_untrt <- (var_untrt$sdev)^2 / sum(var_untrt$sdev^2) * 100
x_label <- paste0("PC1 (", round(perc_var_untrt[1], 1), "%)")
y_label <- paste0("PC2 (", round(perc_var_untrt[2], 1), "%)")

pcaData_untrt <- as.data.frame(var_untrt$x[, 1:2])
colnames(pcaData_untrt) <- c("PC1", "PC2")
pcaData_untrt$line <- vsd_full$line[colnames(vsd_full) %in% untrt_samples]

plot_pca_untrt <- ggplot(pcaData_untrt, aes(x = PC1, y = PC2, color = line)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "B.hybridum Recurrent Selection PCA (Untreated Samples Only)",
       x = x_label,
       y = y_label,
       color = "Line"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold", family = "Times"),
    axis.title = element_text(size = 12, family = "Times"),
    legend.title = element_text(size = 10, family = "Times"),
    legend.text = element_text(size = 9, family = "Times")
  ) +
  scale_color_brewer(palette = "Set2")

pdf(file = "PCA_plot_untrt.pdf", width = 8, height = 8, family = "Times")
plot_pca_untrt
dev.off()

write.csv(var_untrt$rotation, "pca_loadings_untrt.csv", row.names = TRUE)

# ===========>------- Shared setup for DEG heatmaps --------<=========== #

vst_mat <- assay(vsd_full)

annotation_col <- as.data.frame(colData(dds_full)[, c("line", "trt")])
rownames(annotation_col) <- colnames(dds_full)

annotation_col$line <- factor(as.character(annotation_col$line),
                              levels = c("Control", "Drought", "Herbicide"))
annotation_col$trt <- factor(as.character(annotation_col$trt),
                             levels = c("Control", "Pinoxaden"))

# Sanity check
sum(is.na(annotation_col$line))  # should be 0
sum(is.na(annotation_col$trt))   # should be 0

# Sort once, up front -- every heatmap below draws its sample order from this
annotation_col <- annotation_col[order(annotation_col$line, annotation_col$trt), , drop = FALSE]

ann_colors <- list(
  line = setNames(c("grey70", "goldenrod", "firebrick"), levels(annotation_col$line)),
  trt  = setNames(c("steelblue", "orange"), levels(annotation_col$trt))
)

# ------------------------------------------------------------
# Manual gene ID -> functional annotation lookup
# Add to this as you curate new Diamond hits for the updated gene lists.
# Any gene NOT listed here just gets labeled with its gene ID -- the script
# won't break, it'll print a reminder of which ones still need curating.
# ------------------------------------------------------------
gene_annotations <- c(
  "Bradi3g31750.v3.2" = "Tau-Class GST",
  "Bradi3g48770.v3.2" = "Unc. LRR-like",
  "Bradi2g43330.v3.2" = "APTH",
  "Bradi2g13010.v3.2" = "Phi-Class GST",
  "Bradi3g48800.v3.2" = "Unc. LRR-like"
)

# ------------------------------------------------------------
# Helper: build row labels for a set of gene IDs, falling back
# to the bare gene ID when no curated annotation exists yet
# ------------------------------------------------------------
build_row_labels <- function(gene_ids, annotations = gene_annotations) {
  labels <- annotations[gene_ids]
  missing <- gene_ids[is.na(labels)]
  if (length(missing) > 0) {
    message("No curated annotation yet for: ", paste(missing, collapse = ", "))
  }
  labels[is.na(labels)] <- gene_ids[is.na(labels)]
  paste0(labels, " (", gene_ids, ")")
}

# ------------------------------------------------------------
# Helper: z-score heatmap for a gene set, over a chosen set of samples
# (pass all samples for the full picture, or a subset -- e.g. untreated
# samples only -- when the underlying contrast only tested that subset)
# ------------------------------------------------------------
plot_zscore_heatmap <- function(gene_ids, samples, filename, title) {
  ann <- annotation_col[rownames(annotation_col) %in% samples, , drop = FALSE]
  mat <- vst_mat[gene_ids, rownames(ann), drop = FALSE]

  pdf(file = filename, width = 12, height = 8, family = "Times")
  pheatmap(mat,
           cluster_rows      = FALSE,
           cluster_cols      = FALSE,
           show_rownames     = TRUE,
           labels_row        = build_row_labels(gene_ids),
           show_colnames     = TRUE,
           annotation_col    = ann,
           annotation_colors = ann_colors,
           scale             = "row",
           color             = colorRampPalette(c("blue", "white", "red"))(100),
           main              = title)
  dev.off()
}

# ------------------------------------------------------------
# Helper: within-line LFC heatmap -- each line's Pinoxaden samples relative
# to that SAME line's own Control-treatment samples. This directly
# visualizes the quantity an interaction term (line:trt) is testing, so
# it's the natural plot for resistance-specific candidate genes.
# ------------------------------------------------------------
plot_within_line_lfc_heatmap <- function(gene_ids, filename, title) {
  mat <- vst_mat[gene_ids, rownames(annotation_col), drop = FALSE]
  lines <- levels(annotation_col$line)

  within_line_lfc <- mat
  within_line_lfc[] <- NA

  for (ln in lines) {
    ref_samples <- rownames(annotation_col)[
      annotation_col$line == ln & annotation_col$trt == "Control"
    ]
    line_samples <- rownames(annotation_col)[annotation_col$line == ln]

    if (length(ref_samples) == 0) {
      warning(paste("No Control-treatment reference samples found for line:", ln))
      next
    }

    ref_means <- rowMeans(mat[, ref_samples, drop = FALSE])
    within_line_lfc[, line_samples] <- sweep(mat[, line_samples, drop = FALSE],
                                             1, ref_means, FUN = "-")
  }

  max_val <- max(abs(within_line_lfc), na.rm = TRUE)

  pdf(file = filename, width = 12, height = 8, family = "Times")
  pheatmap(within_line_lfc,
           cluster_rows      = FALSE,
           cluster_cols      = FALSE,
           show_rownames     = TRUE,
           labels_row        = build_row_labels(gene_ids),
           show_colnames     = TRUE,
           annotation_col    = annotation_col,
           annotation_colors = ann_colors,
           scale             = "none",
           color             = colorRampPalette(c("blue", "white", "red"))(100),
           breaks            = seq(-max_val, max_val, length.out = 101),
           main              = title)
  dev.off()
}

# ===========>------- Heatmap set 1: shared pinoxaden response --------<=========== #
# sig_trt_genes -- from dds_full, the treatment effect shared/averaged
# across all three lines. Answers: "what responds to pinoxaden broadly?"

deg_ids_trt <- rownames(sig_trt_genes)

if (length(deg_ids_trt) > 0) {
  plot_zscore_heatmap(
    gene_ids = deg_ids_trt,
    samples  = rownames(annotation_col),
    filename = "vst_heatmap_shared_response.pdf",
    title    = "Z-score of Shared Pinoxaden-Responsive DEGs\nAcross Lines and Treatment"
  )

  plot_within_line_lfc_heatmap(
    gene_ids = deg_ids_trt,
    filename = "lfc_heatmap_shared_response.pdf",
    title    = "Within-Line Log2 Fold Change of Shared Pinoxaden-Responsive DEGs"
  )
} else {
  message("No shared-response DEGs to plot -- skipping this heatmap pair.")
}

# ===========>------- Heatmap set 2: resistance-specific response --------<=========== #
# sig_drt_genes / sig_herb_genes -- interaction terms (line:trt) from
# dds_int. Genes whose response to pinoxaden differs from Control
# specifically in the Drought/Herbicide lines -- your TF-search candidates.
# Within-line LFC is the primary visual, since it's literally the quantity
# the interaction term tests.

deg_ids_resist <- union(rownames(sig_drt_genes), rownames(sig_herb_genes))

if (length(deg_ids_resist) > 0) {
  plot_within_line_lfc_heatmap(
    gene_ids = deg_ids_resist,
    filename = "lfc_heatmap_resistance_specific.pdf",
    title    = "Within-Line Log2 Fold Change of Resistance-Specific DEGs\n(Drought/Herbicide Response vs. Control Line Response)"
  )

  plot_zscore_heatmap(
    gene_ids = deg_ids_resist,
    samples  = rownames(annotation_col),
    filename = "vst_heatmap_resistance_specific.pdf",
    title    = "Z-score of Resistance-Specific DEGs\nAcross Lines and Treatment"
  )
} else {
  message("No resistance-specific (interaction-term) DEGs to plot -- skipping this heatmap pair.")
}

# ===========>------- Heatmap set 3: baseline (untreated) line differences --------<=========== #
# sig_baseline_drt / sig_baseline_hrb -- the line effect at trt = Control.
# Restricted to untreated columns only: this contrast never tested the
# treated samples, so including them would visualize variance the test
# wasn't about.

deg_ids_baseline <- union(rownames(sig_baseline_drt), rownames(sig_baseline_hrb))

if (length(deg_ids_baseline) > 0) {
  plot_zscore_heatmap(
    gene_ids = deg_ids_baseline,
    samples  = untrt_samples,
    filename = "vst_heatmap_baseline.pdf",
    title    = "Z-score of Baseline (Untreated) Line-Differentiated DEGs"
  )
} else {
  message("No baseline DEGs to plot -- skipping this heatmap.")
}

message("All figures written to working directory.")
