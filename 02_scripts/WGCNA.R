################################################################################
# Date: 5-14-2026
# Author: John Lemas (Optimized & Debugged using Gemini)
# Description: BFT project WGCNA with Integrated DESeq2 Output
################################################################################

# 1. SETUP & PACKAGE LOADING
options(stringsAsFactors = FALSE)

required_packages <- c("DESeq2", "tidyverse", "WGCNA", "ggdendro")
for (package in required_packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    if (package %in% c("DESeq2", "WGCNA")) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(package)
    } else {
      install.packages(package)
    }
  }
  library(package, character.only = TRUE)
}

# Set environment
setwd("/Users/jmlemas/Library/CloudStorage/OneDrive-MichiganStateUniversity/Dissertation/Chapter3/RNA_seq/03_outputs/04_DESeq2/")

# 2. LOAD PRE-COMPUTED DATA
dds           <- readRDS("dds.rds")
resLFC_trt    <- readRDS("resLFC_trt.rds")
resLFC_line   <- readRDS("resLFC_line.rds")
sig_trt_genes <- readRDS("sig_trt_genes.rds")
coldata       <- as.data.frame(colData(dds))

# 3. TRANSFORMATION & QUALITY CHECK
dds_vst   <- varianceStabilizingTransformation(dds, blind = FALSE)
geneCount <- t(assay(dds_vst))

gsg <- goodSamplesGenes(geneCount, verbose = 3)
if (!gsg$allOK) {
  geneCount <- geneCount[gsg$goodSamples, gsg$goodGenes]
}

# 4. SAMPLE CLUSTERING DENDROGRAM
sampleTree    <- hclust(dist(geneCount), method = "average")
sampleTree_gg <- dendro_data(sampleTree)

pdf(file = "Tree_DESeq2_means_combined.pdf", width = 18, height = 9, family = "Times")
print(
  ggplot() +
    geom_segment(data = segment(sampleTree_gg), aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_text(data = label(sampleTree_gg), aes(x = x, y = y, label = label), hjust = -0.1, size = 3) +
    coord_flip() + scale_y_reverse() + theme_minimal() +
    labs(title = "Hierarchical Clustering Dendrogram")
)
dev.off()

# 5. CHOOSE SOFT-THRESHOLDING POWER
powers <- c(1:10, seq(from = 12, to = 24, by = 2))
disableWGCNAThreads()
sft <- pickSoftThreshold(geneCount, powerVector = powers, verbose = 3)
write.csv(sft$fitIndices, "power_table.csv", row.names = FALSE)

pdf(file = "soft_power.pdf", width = 12, height = 5, family = "Times")
par(mfrow = c(1,2), cex = 0.9)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit,signed R^2", type="n", main="Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, col="red")
abline(h=0.90, col="red"); abline(h=0.80, col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n", main="Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, col="red")
dev.off()

# 6. BLOCKWISE NETWORK CONSTRUCTION
softPower <- 5
enableWGCNAThreads(nThreads = 12)

# Protect and override namespace safely
tmp_cor <- cor; cor <- WGCNA::cor

net <- blockwiseModules(geneCount, 
                        power = softPower, 
                        corType = "pearson",
                        networkType = "signed",
                        TOMType = "signed", 
                        minModuleSize = 30,
                        maxBlockSize = 20000,
                        reassignThreshold = 0, 
                        mergeCutHeight = 0.25,
                        numericLabels = TRUE, 
                        pamRespectsDendro = FALSE,
                        saveTOMs = TRUE,
                        saveTOMFileBase = "B.hyb_TOM",
                        verbose = 3,
                        nThreads = 12)

cor <- tmp_cor # Restore original cor function
mergedColors <- labels2colors(net$colors)

# Plot Final Network Dendrogram
pdf("dendrogram_final_means.pdf", width = 12, height = 9, family = "Times")
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

# Save simple tracking module key
tibble(gene_id = names(net$colors), colors = mergedColors) %>% 
  write_delim("gene_modules.tsv", delim = "\t")

# 7. MODULE-TRAIT CORRELATION HEATMAP
MEs <- orderMEs(moduleEigengenes(geneCount, colors = mergedColors)$eigengenes)
module_order <- gsub("ME", "", names(MEs))

correlation_metadata <- read.table("../../01_inputs/correlation_metadata.txt", header = TRUE)

moduleTraitcor    <- cor(MEs, correlation_metadata, use = "p")
nSamples          <- nrow(geneCount)
moduleTraitPvalue <- corPvalueStudent(moduleTraitcor, nSamples)

testMatrix <- matrix(paste0(signif(moduleTraitcor, 2), "\n(", signif(moduleTraitPvalue, 1), ")"),
                     nrow = nrow(moduleTraitcor), ncol = ncol(moduleTraitcor), dimnames = dimnames(moduleTraitcor))

pdf(file = "MEheatmap.pdf", width = 12, height = 12, family = "Times", pointsize = 10)
labeledHeatmap(Matrix = moduleTraitcor, xLabels = colnames(correlation_metadata), yLabels = names(MEs), 
               ySymbols = module_order, colorLabels = FALSE, colors = blueWhiteRed(50), 
               textMatrix = testMatrix, setStdMargins = TRUE, yLabelsPosition = "left", cex.text = 1.0, zlim = c(-1, 1))
dev.off()

# 8. EIGENGENE EXPRESSION PLOTS
coldata_augmented <- coldata %>%
  rownames_to_column("Sample") %>%
  mutate(treatment_group = factor(paste(line, trt, sep = "_")))

# Order treatment groups (UT then T)
ordered_groups <- coldata_augmented %>%
  distinct(treatment_group, trt) %>%
  arrange(factor(trt, levels = c("Control", "Pinoxaden"))) %>%
  pull(treatment_group)

MEs_meta <- MEs %>%
  rownames_to_column("Sample") %>%
  left_join(coldata_augmented, by = "Sample")

# Boxplot of Target Modules
modules_of_interest <- c("MEyellow", "MEgreenyellow", "MEblue", "MEtan") %>% intersect(colnames(MEs_meta))
if(length(modules_of_interest) > 0) {
  pdf(file = "ME_dotplots.pdf", width = 8, height = 6, family = "Times")
  print(
    MEs_meta %>%
      dplyr::select(treatment_group, all_of(modules_of_interest)) %>%
      pivot_longer(cols = -treatment_group, names_to = "Modules", values_to = "Eigenvalues") %>%
      mutate(treatment_group = factor(treatment_group, levels = ordered_groups)) %>%
      ggplot(aes(x = treatment_group, y = Eigenvalues, color = treatment_group)) +
      geom_boxplot(outliers = FALSE, show.legend = FALSE) +
      geom_point(size = 1.0, alpha = 0.5, position = position_jitter(width = 0.2), show.legend = FALSE) +
      facet_wrap(~Modules, scales = "free_y") + theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 14)) +
      labs(title = "Eigengene Expression by Treatment Group", x = "Treatment", y = "Eigengene Value")
  )
  dev.off()
}

# 9. ADVANCED METRIC COMPUTATIONS (MM & GS)
gene_ids <- colnames(geneCount)
modNames <- substring(names(MEs), 3)

# Module Membership (MM) across all modules
geneModuleMembership <- as.data.frame(cor(geneCount, MEs, use = "p"))
MMPvalue             <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
colnames(geneModuleMembership) <- paste0("MM.", modNames)
colnames(MMPvalue)             <- paste0("p.MM.", modNames)

# FIX: Robust dynamic assignment loop handling "grey" module cases safely
MM_own   <- numeric(length(gene_ids))
p.MM_own <- numeric(length(gene_ids))

for (i in seq_along(gene_ids)) {
  target_mod <- mergedColors[i]
  col_idx    <- match(target_mod, modNames)
  
  if (!is.na(col_idx)) {
    MM_own[i]   <- abs(geneModuleMembership[i, col_idx])
    p.MM_own[i] <- MMPvalue[i, col_idx]
  } else {
    MM_own[i]   <- 0  # Grey module exceptions default gracefully
    p.MM_own[i] <- 1
  }
}

# FIX: Correcting the undefined `names(Temp)` bug for Gene Significance
# Using your metadata metric column directly (assumed "R" column represents Genotype)
trait_vector <- as.data.frame(correlation_metadata$Treated)
colnames(trait_vector) <- "Genotype"

geneTraitSignificance <- as.data.frame(cor(geneCount, trait_vector, use = "p"))
GSPvalue             <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
colnames(geneTraitSignificance) <- "GS.Genotype"
colnames(GSPvalue)             <- "p.GS.Genotype"

# 10. GENERATING THE COMPREHENSIVE COMBINED TABLE
wgcna_base <- tibble(
  gene_id = gene_ids,
  WGCNA_module = mergedColors,
  MM_OwnModule = MM_own,
  p.MM_OwnModule = p.MM_own,
  GS_Genotype = abs(geneTraitSignificance$GS.Genotype),
  p.GS_Genotype = GSPvalue$p.GS.Genotype
)

# Coerce and prefix DESeq2 tables to keep clean column traces
deseq_line_df <- as.data.frame(resLFC_line) %>% 
  rownames_to_column("gene_id") %>% 
  dplyr::select(gene_id, baseMean, log2FoldChange, pvalue, padj) %>% 
  rename_with(~paste0("DESeq2_gen_", .), -gene_id)

deseq_trt_df <- as.data.frame(resLFC_trt) %>% 
  rownames_to_column("gene_id") %>% 
  dplyr::select(gene_id, log2FoldChange, pvalue, padj) %>% 
  rename_with(~paste0("DESeq2_trt_", .), -gene_id)

# Perform clean sequential joins
combined_master_table <- wgcna_base %>%
  left_join(deseq_gen_df, by = "gene_id") %>%
  left_join(deseq_trt_df, by = "gene_id") %>%
  mutate(is_hub_gene = if_else(MM_OwnModule >= 0.8 & GS_Genotype >= 0.8, TRUE, FALSE)) %>% 
  arrange(DESeq2_gen_padj)

# Save Master Table Output
write.csv(combined_master_table, "DESeq2_WGCNA_master_metrics.csv", row.names = FALSE)
cat("\nPipeline successful! Master metrics saved to 'DESeq2_WGCNA_master_metrics.csv'\n")