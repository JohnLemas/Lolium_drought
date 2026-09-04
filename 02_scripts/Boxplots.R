################################################################################
# Date: 5-14-2026
# Author: John Lemas (Optimized via Custom Functions via Gemini)
# Description: Modularized DEG Boxplot Engine for L. corniculatus RNA-Seq
################################################################################

# 1. SETUP & PACKAGE LOADING
required_packages <- c("tidyverse", "DESeq2", "RColorBrewer", "ggpubr")

for (package in required_packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    if (package == "DESeq2") {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(package)
    } else {
      install.packages(package)
    }
  }
  library(package, character.only = TRUE)
}

# Set working directory
setwd("/Users/jmlemas/Library/CloudStorage/OneDrive-MichiganStateUniversity/Dissertation/Chapter3/RNA_seq/03_outputs/04_DESeq2/")

# Whole dataset contrast
dds           <- readRDS("dds.rds")
resLFC_trt    <- readRDS("resLFC_trt.rds")
sig_trt_genes <- readRDS("sig_trt_genes.rds")

sig_trt_genes

plot_data <- plotCounts(dds, gene = "Bradi3g31827.v3.2", 
                        intgroup = c("trt", "line"), returnData = TRUE)

p <- ggplot(plot_data, aes(x = line, y = count)) + 
  geom_boxplot(alpha = 0.7, outlier.shape = NA, position = position_dodge(0.8), show.legend = FALSE) + 
  geom_jitter(position = position_dodge(0.8), size = 2, color = "black", alpha = 0.7, show.legend = FALSE) +
  scale_y_log10() + 
  scale_fill_brewer(palette = "Set2") +
  labs(title = bquote(italic("B. hybridum") ~ "gene1"),
       subtitle = "ID = Bradi2g43330.v3.2",
       y = "Normalized Counts (log10 scale)",
       x = "Line") + 
  theme_bw() + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", family = "Times"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, family = "Times"),
    axis.title = element_text(size = 12, family = "Times"),
    axis.text = element_text(family = "Times"),
    legend.title = element_text(size = 10, family = "Times"),
    legend.text = element_text(size = 9, family = "Times"),
    strip.text = element_text(size = 11, face = "bold", family = "Times") 
  ) +
  stat_compare_means(
    aes(group = line),
    method = "t.test",        
    label = "p.format"        
  )

p

pdf("UNK_plot.pdf", width = 8, height = 8, family = "Times")
p
dev.off()
