#######################################

#__Date:__ 8-12-2026
#__Author:__ John Lemas
# List of packages required
packages <- c("ggVennDiagram",
              "ggplot2")

# if not installed, install and load it
for (package in packages) {
  if (!(package %in% rownames(installed.packages()))) {
    BiocManager::install(package)
  }
  library(package, character.only = TRUE)
}

#################################################

# Set the working directory:
setwd("/Users/jmlemas/Library/CloudStorage/OneDrive-MichiganStateUniversity/Dissertation/Chapter3/RNA_seq/03_outputs/04_DESeq2/")

# ===========>------- Load saved objects --------<=========== #
sig_trt_genes    <- read.table("sig_trt_genes.txt",    header = TRUE, sep = "\t", row.names = 1)
sig_drt_genes    <- read.table("sig_drt_genes.txt",    header = TRUE, sep = "\t", row.names = 1)
sig_herb_genes   <- read.table("sig_herb_genes.txt",   header = TRUE, sep = "\t", row.names = 1)
sig_baseline_drt <- read.table("sig_baseline_drt.txt", header = TRUE, sep = "\t", row.names = 1)
sig_baseline_hrb <- read.table("sig_baseline_hrb.txt", header = TRUE, sep = "\t", row.names = 1)

# each element is a vector of gene IDs (e.g., rownames of your sig DEG subset)
gene_lists <- list(
  drt_trt = rownames(sig_drt_genes),
  herb_trt = rownames(sig_herb_genes),
  drt_base = rownames(sig_baseline_drt),
  hrb_base = rownames(sig_baseline_hrb)
)

p <- ggVennDiagram(gene_lists, label = "count") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  ggtitle("DEG overlap")

p + scale_fill_gradient(
  low = "white", high = "steelblue",
  trans = "log1p"
)
