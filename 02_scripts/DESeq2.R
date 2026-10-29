#######################################
# DESeq2: herbicide resistance (R/S) x drought, Lolium
# Revised template based on DESeq2.R
#######################################

packages <- c("DESeq2", "tidyverse", "ashr")
for (package in packages) {
  if (!(package %in% rownames(installed.packages()))) BiocManager::install(package)
  library(package, character.only = TRUE)
}

setwd("/Users/jmlemas/Library/CloudStorage/OneDrive-MichiganStateUniversity/Dissertation/Chapter3/Lolium_samples/RNA_seq")
out <- "03_outputs/04_DESeq2/"

alpha   <- 0.05   # FDR
lfc_cut <- 1      # post-hoc |log2FC| cut (2 is very strict for interaction terms)

# ---------------- Metadata ----------------
meta <- read.table("01_inputs/Lol_metadata.txt", header = FALSE, sep = "\t")
colnames(meta) <- c("fasta1", "name1", "rep", "trt", "status", "species", "achl")  # no hyphen in achl
rownames(meta) <- meta$name1

# population = line ID without replicate number (KB, OR, PR, EMR, GMR, ARG)
meta$population <- sub("[0-9]+(\\.[0-9]+)?-.*$", "", meta$name1)

coldata <- meta[, c("status", "trt", "species", "population", "rep", "achl")]
coldata$status  <- factor(coldata$status,  levels = c("S", "R"))          # S = reference
coldata$trt     <- factor(coldata$trt,     levels = c("untrt", "trt"))    # untrt = reference
coldata$species <- factor(coldata$species, levels = c("LOLMU", "LOLRI"))
coldata$population <- factor(coldata$population)
coldata$rep <- factor(coldata$rep)

# Sanity check of the design
print(with(coldata, table(population, trt)))
print(with(coldata, table(status, species, trt)))

# ---------------- Counts ----------------
countsData <- read.table("03_outputs/03_featureCounts/counts_multi.txt",
                         header = FALSE, row.names = 1, skip = 2)
n_samples <- ncol(countsData) - 5
if (n_samples != nrow(meta)) {
  stop("counts has ", n_samples, " sample columns but metadata has ", nrow(meta),
       " rows. Match columns to metadata by BAM name before continuing.")
}
colnames(countsData) <- c("chr", "start", "stop", "strand", "length", meta$name1)
cts <- round(as.matrix(countsData[, 6:ncol(countsData)]))
stopifnot(all(rownames(coldata) == colnames(cts)))

# ---------------- Model ----------------
# species is crossed with status, so it is adjusted for.
# Coefficients: species_LOLRI_vs_LOLMU, status_R_vs_S, trt_trt_vs_untrt, statusR.trttrt
dds <- DESeqDataSetFromMatrix(cts, coldata,
                              design = ~ species + status + trt + status:trt)

keep <- rowSums(counts(dds) >= 10) >= 5   # >= smallest group size (S untreated n = 5)
dds <- dds[keep, ]
dds <- DESeq(dds)
saveRDS(dds, file.path(out, "dds_int.rds"))
print(resultsNames(dds))

# ---------------- QC: look at PR1-0D and population structure ----------------
vsd <- vst(dds, blind = TRUE)
pca <- plotPCA(vsd, intgroup = c("trt", "status", "species"), returnData = TRUE)
pca$population <- coldata$population
pca$name <- rownames(pca)
p <- ggplot(pca, aes(PC1, PC2, colour = trt, shape = status)) +
  geom_point(size = 3) + geom_text(aes(label = name), size = 2, vjust = -1)
ggsave(file.path(out, "pca_vst.pdf"), p, width = 8, height = 6)

# ---------------- Contrasts ----------------
run_contrast <- function(spec, label) {
  if (is.character(spec)) {
    res <- results(dds, name = spec, alpha = alpha)
    shr <- lfcShrink(dds, coef = spec, res = res, type = "ashr")
  } else {
    res <- results(dds, contrast = spec, alpha = alpha)
    shr <- lfcShrink(dds, contrast = spec, res = res, type = "ashr")
  }
  message("---- ", label)
  print(summary(res))
  shr$padj <- res$padj   # keep the test's padj, shrunken LFC for effect size
  df  <- as.data.frame(shr)
  sig <- df[which(df$padj < alpha & abs(df$log2FoldChange) > lfc_cut), ]
  sig <- sig[order(sig$padj), ]
  write.table(df,  paste0(out, "all_",  label, ".txt"), quote = FALSE, sep = "\t")
  write.table(sig, paste0(out, "sig_",  label, ".txt"), quote = FALSE, sep = "\t")
  saveRDS(shr, paste0(out, "res_", label, ".rds"))
  invisible(df)
}

contrasts <- list(
  # 1. Drought response in susceptible plants
  drought_in_S      = "trt_trt_vs_untrt",
  # 2. Drought response in resistant plants
  drought_in_R      = list(c("trt_trt_vs_untrt", "statusR.trttrt")),
  # 3. KEY: does the drought response differ between R and S? (interaction)
  interaction_RxS   = "statusR.trttrt",
  # 4. Constitutive R vs S differences (untreated)
  baseline_R_vs_S   = "status_R_vs_S",
  # 5. R vs S under drought
  drought_R_vs_S    = list(c("status_R_vs_S", "statusR.trttrt"))
)

res_list <- imap(contrasts, ~ run_contrast(.x, .y))

# ---------------- Candidate classes ----------------
get_sig <- function(df) rownames(df)[which(df$padj < alpha & abs(df$log2FoldChange) > lfc_cut)]
s_int  <- get_sig(res_list$interaction_RxS)
s_base <- get_sig(res_list$baseline_R_vs_S)
s_dR   <- get_sig(res_list$drought_in_R)
s_dS   <- get_sig(res_list$drought_in_S)

candidates <- list(
  differential_drought_response = s_int,                   # R responds differently than S
  constitutive_and_responsive   = intersect(s_base, s_dR), # R-baseline difference AND drought-responsive in R
  R_specific_drought            = setdiff(s_dR, s_dS),     # responds in R only
  shared_drought                = intersect(s_dR, s_dS)
)
print(lengths(candidates))
saveRDS(candidates, file.path(out, "candidate_gene_sets.rds"))

# ---------------- Robustness across populations ----------------
# Status is tested with only 2 S vs 4 R populations; replicates within a population
# are not independent for the status effect. For top hits, check that the direction
# is consistent in each population (mean vst by population x trt):
top <- head(s_int, 50)
if (length(top) > 0) {
  m <- assay(vsd)[top, , drop = FALSE]
  popmeans <- sapply(split(seq_len(ncol(m)), paste(coldata$population, coldata$trt)),
                     function(i) rowMeans(m[, i, drop = FALSE]))
  write.table(popmeans, paste0(out, "top_interaction_pop_means.txt"), quote = FALSE, sep = "\t")
}

# ---------------- Optional: chlorophyll fluorescence ----------------
# achl is almost collinear with trt, so only use it WITHIN treated samples:
# genes whose expression tracks the severity of photosystem damage.
trt_idx <- coldata$trt == "trt"
dds_t <- DESeqDataSetFromMatrix(cts[, trt_idx], coldata[trt_idx, ],
                                design = ~ species + achl)
dds_t$achl <- scale(dds_t$achl)[, 1]
dds_t <- dds_t[rowSums(counts(dds_t) >= 10) >= 5, ]
dds_t <- DESeq(dds_t)
res_achl <- results(dds_t, name = "achl", alpha = alpha)
summary(res_achl)
write.table(as.data.frame(res_achl), paste0(out, "all_achl_within_treated.txt"),
            quote = FALSE, sep = "\t")
