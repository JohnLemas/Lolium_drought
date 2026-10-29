#######################################
# Figures for Lolium R/S x drought RNA-seq
# Run after DESeq2_resistance_x_drought.R (uses its saved outputs)
# Makes: PCA, heatmaps, volcano plots, UpSet and Venn plots
#######################################

packages <- c("DESeq2", "tidyverse", "ggrepel", "pheatmap", "EnhancedVolcano",
              "UpSetR", "ggVennDiagram", "limma", "RColorBrewer")
for (package in packages) {
  if (!(package %in% rownames(installed.packages()))) BiocManager::install(package)
  library(package, character.only = TRUE)
}

setwd("/Users/jmlemas/Library/CloudStorage/OneDrive-MichiganStateUniversity/Dissertation/Chapter3/Lolium_samples/RNA_seq")
indir  <- "03_outputs/04_DESeq2/"
figdir <- "03_outputs/05_figures/"
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

alpha   <- 0.05
lfc_cut <- 1

# ---------------- Load ----------------
dds <- readRDS(paste0(indir, "dds_int.rds"))
vsd <- vst(dds, blind = FALSE)          # blind = FALSE: design-aware, better for heatmaps
cd  <- as.data.frame(colData(dds))
cd$group <- factor(paste(cd$status, cd$trt, sep = "_"),
                   levels = c("S_untrt", "S_trt", "R_untrt", "R_trt"))

contrast_names <- c("drought_in_S", "drought_in_R", "interaction_RxS",
                    "baseline_R_vs_S", "drought_R_vs_S")
res_list <- setNames(lapply(contrast_names, function(x)
  as.data.frame(readRDS(paste0(indir, "res_", x, ".rds")))), contrast_names)

# Colours used throughout
col_group <- c(S_untrt = "#9ecae1", S_trt = "#08519c", R_untrt = "#fdae6b", R_trt = "#a63603")
col_species <- c(LOLMU = "#1b9e77", LOLRI = "#7570b3")

# ---------------- 1. PCA ----------------
pca_plot <- function(mat, title, file) {
  mat <- mat[order(rowVars(mat), decreasing = TRUE)[1:min(500, nrow(mat))], ]
  pc  <- prcomp(t(mat))
  pv  <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
  df  <- data.frame(pc$x[, 1:3], cd[rownames(pc$x), c("species", "population", "achl")],
                    group = cd[rownames(pc$x), "group"], name = rownames(pc$x))
  p <- ggplot(df, aes(PC1, PC2, colour = group, shape = species)) +
    geom_point(size = 3.5, alpha = 0.9) +
    geom_text_repel(aes(label = name), size = 2.3, max.overlaps = 8, show.legend = FALSE) +
    scale_colour_manual(values = col_group) +
    labs(x = paste0("PC1 (", pv[1], "%)"), y = paste0("PC2 (", pv[2], "%)"), title = title) +
    theme_bw(base_size = 12)
  ggsave(file, p, width = 8, height = 6)
  # PC2 vs PC3 and PC1 coloured by chlorophyll fluorescence
  p2 <- ggplot(df, aes(PC1, PC2, colour = achl, shape = group)) +
    geom_point(size = 3.5) + scale_colour_viridis_c(name = "a-chl (Fv/Fm)") +
    labs(x = paste0("PC1 (", pv[1], "%)"), y = paste0("PC2 (", pv[2], "%)"),
         title = paste(title, "- coloured by fluorescence")) + theme_bw(base_size = 12)
  ggsave(sub("\\.pdf$", "_achl.pdf", file), p2, width = 8, height = 6)
}

pca_plot(assay(vsd), "PCA, top 500 variable genes (vst)", paste0(figdir, "pca_vst.pdf"))

# Same PCA with the species effect removed, to see status/drought structure
vsd_nosp <- removeBatchEffect(assay(vsd), batch = cd$species,
                              design = model.matrix(~ status + trt, cd))
pca_plot(vsd_nosp, "PCA, species effect removed (visualization only)",
         paste0(figdir, "pca_vst_species_removed.pdf"))

# ---------------- 2. Heatmaps ----------------
ann_col <- data.frame(group = cd$group, species = cd$species, population = cd$population,
                      row.names = rownames(cd))
ann_colors <- list(group = col_group, species = col_species)
col_order <- order(cd$group, cd$species, cd$population)

make_heatmap <- function(genes, title, file, max_genes = 60) {
  genes <- head(genes, max_genes)
  genes <- genes[genes %in% rownames(vsd)]
  if (length(genes) < 2) { message("Skipping heatmap: ", title); return(invisible()) }
  mat <- assay(vsd)[genes, col_order, drop = FALSE]
  pheatmap(mat, scale = "row", cluster_cols = FALSE,
           annotation_col = ann_col, annotation_colors = ann_colors,
           color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
           show_colnames = TRUE, fontsize_row = 6, fontsize_col = 7,
           gaps_col = cumsum(table(cd$group[col_order])),
           main = title, filename = file, width = 11, height = max(5, 0.14 * length(genes) + 3))
}

top_genes <- function(df, n = 60) {
  df <- df[which(df$padj < alpha & abs(df$log2FoldChange) > lfc_cut), ]
  rownames(df)[order(df$padj)][seq_len(min(n, nrow(df)))]
}

make_heatmap(top_genes(res_list$interaction_RxS),
             "Top genes: R x drought interaction", paste0(figdir, "heatmap_interaction.pdf"))
make_heatmap(top_genes(res_list$baseline_R_vs_S),
             "Top genes: constitutive R vs S (untreated)", paste0(figdir, "heatmap_baseline.pdf"))
make_heatmap(top_genes(res_list$drought_in_R),
             "Top genes: drought response in R", paste0(figdir, "heatmap_drought_R.pdf"))

# Candidates that are both constitutive and drought-responsive in R
cand <- readRDS(paste0(indir, "candidate_gene_sets.rds"))
make_heatmap(cand$constitutive_and_responsive,
             "Constitutive R difference AND drought-responsive in R",
             paste0(figdir, "heatmap_constitutive_and_responsive.pdf"))

# ---------------- 3. Volcano plots ----------------
for (nm in contrast_names) {
  df <- res_list[[nm]]
  df <- df[!is.na(df$padj), ]
  p <- EnhancedVolcano(df, lab = rownames(df), x = "log2FoldChange", y = "padj",
                       pCutoff = alpha, FCcutoff = lfc_cut,
                       title = nm, subtitle = paste0("padj < ", alpha, ", |log2FC| > ", lfc_cut,
                                                     " (ashr-shrunk LFC)"),
                       ylab = bquote(-log[10]~"adj. p"), pointSize = 1.6, labSize = 3,
                       drawConnectors = TRUE, max.overlaps = 15, caption = NULL)
  ggsave(paste0(figdir, "volcano_", nm, ".pdf"), p, width = 7.5, height = 7)
}

# ---------------- 4. UpSet and Venn ----------------
sig_sets <- lapply(res_list, function(df)
  rownames(df)[which(df$padj < alpha & abs(df$log2FoldChange) > lfc_cut)])
sig_sets <- sig_sets[lengths(sig_sets) > 0]
print(lengths(sig_sets))

if (length(sig_sets) >= 2) {
  pdf(paste0(figdir, "upset_all_contrasts.pdf"), width = 9, height = 6)
  print(upset(fromList(sig_sets), nsets = length(sig_sets), nintersects = 30,
              order.by = "freq", mb.ratio = c(0.6, 0.4),
              text.scale = 1.2, sets.bar.color = "grey30"))
  dev.off()
}

# Up vs down split for the interaction and baseline contrasts
dir_sets <- list()
for (nm in c("interaction_RxS", "baseline_R_vs_S", "drought_in_R", "drought_in_S")) {
  df <- res_list[[nm]]
  s <- which(df$padj < alpha & abs(df$log2FoldChange) > lfc_cut)
  dir_sets[[paste0(nm, "_up")]]   <- rownames(df)[s[df$log2FoldChange[s] > 0]]
  dir_sets[[paste0(nm, "_down")]] <- rownames(df)[s[df$log2FoldChange[s] < 0]]
}
dir_sets <- dir_sets[lengths(dir_sets) > 0]
if (length(dir_sets) >= 2) {
  pdf(paste0(figdir, "upset_direction.pdf"), width = 11, height = 6.5)
  print(upset(fromList(dir_sets), nsets = length(dir_sets), nintersects = 30,
              order.by = "freq", text.scale = 1.1))
  dev.off()
}

# Venn: drought response in S vs R, plus the baseline and interaction sets (max 4 sets)
venn_sets <- sig_sets[intersect(c("drought_in_S", "drought_in_R", "baseline_R_vs_S",
                                  "interaction_RxS"), names(sig_sets))]
if (length(venn_sets) >= 2) {
  pv <- ggVennDiagram(venn_sets, label_alpha = 0) +
    scale_fill_gradient(low = "white", high = "#a63603") +
    theme(legend.position = "none")
  ggsave(paste0(figdir, "venn_contrasts.pdf"), pv, width = 7, height = 6)
}

message("Figures written to ", figdir)