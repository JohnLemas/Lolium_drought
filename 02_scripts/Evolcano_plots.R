#######################################

#__Date:__ 5-14-2026
#__Author:__ John Lemas
# List of packages required
packages <- c("tidyverse",
              "EnhancedVolcano")

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

# Whole dataset contrast
res_trt    <- readRDS("res_trt.rds")

# Interaction datasets
resLFC_interaction_hrb <- readRDS("resLFC_interaction_hrb.rds")
resLFC_interaction_drt <- readRDS("resLFC_interaction_drt.rds")
resLFC_baseline_dr <- readRDS("resLFC_baseline_drt.rds")
resLFC_baseline_hrb <- readRDS("resLFC_baseline_hrb.rds")

# =========>-------- Whole dataset plot -------<========== #
# treated and untreated samples
# deseq2: design ~ line + trt
# restuls: contrast = c("trt", "Pinoxaden", "Control")

trt_plot <- EnhancedVolcano(toptable = res_trt, lab = rownames(res_trt), 
                                x = 'log2FoldChange',
                                y = "pvalue",
                                title = 'Treatment Results',
                                subtitle = 'contrast = c("trt", "Pinoxaden", "Control")',
                                caption = bquote(~Log[2]~ "fold change cutoff = ±2; p-value cutoff = 0.01"),
                                xlab = bquote(~Log[2]~ 'fold change'),
                                selectLab = c('Bradi3g31750.v3.2',
                                              'Bradi3g48770.v3.2',
                                              'Bradi2g43330.v3.2',
                                              'Bradi2g13010.v3.2',
                                              'Bradi3g48800.v3.2'),
                                pCutoff = 0.01,
                                FCcutoff = 2,
                                cutoffLineType = 'twodash',
                                cutoffLineWidth = 0.8,
                                pointSize = 4.0,
                                labSize = 4.5,
                                labFace = 'bold',
                                boxedLabels = TRUE,
                                colAlpha = 0.5,
                                legendPosition = 'none',
                                drawConnectors = TRUE,
                                widthConnectors = 0.5,
                                colConnectors = 'black')

trt_plot

dev.off()

pdf("trt_Evolcano.pdf", 
    width = 8, height = 8)

par(mfrow = c(1,1))

trt_plot

dev.off()

# =========>-------- Herbicide selected lines plot -------<========== #
# Herbicide lineage response to treatment
# deseq2: design = ~ line + trt + line:trt
# results: name = "lineHerbicide.trtPinoxaden"

hrb_plot <- EnhancedVolcano(toptable = resLFC_interaction_hrb, lab = rownames(resLFC_interaction_hrb), 
                            x = 'log2FoldChange',
                            y = "pvalue",
                            title = 'Within-line analysis:\nRecurrent Pinoxaden Selection',
                            subtitle = 'name = "lineHerbicide.trtPinoxaden"',
                            caption = bquote(~Log[2]~ "fold change cutoff = ±2; p-value cutoff = 0.01"),
                            xlab = bquote(~Log[2]~ 'fold change'),
                            pCutoff = 0.01,
                            FCcutoff = 2,
                            selectLab = c('Bradi3g46720.v3.2',
                                          'Bradi4g39920.v3.2',
                                          'Bradi4g35342.v3.2',
                                          'Bradi2g45080.v3.2',
                                          'Bradi2g31950.v3.2',
                                          'Bradi1g12690.v3.2'),
                            cutoffLineType = 'twodash',
                            cutoffLineWidth = 0.8,
                            pointSize = 4.0,
                            labSize = 4.5,
                            labFace = 'bold',
                            boxedLabels = TRUE,
                            colAlpha = 0.5,
                            legendPosition = 'none',
                            drawConnectors = TRUE,
                            widthConnectors = 0.5,
                            colConnectors = 'black')

hrb_plot

dev.off()

pdf("hrb_Evolcano.pdf", 
    width = 8, height = 8)

par(mfrow = c(1,1))

hrb_plot

dev.off()



# =========>-------- Drought selected lines plot -------<========== #
# Drought lineage response to treatment
# deseq2: design = ~ line + trt + line:trt
# results: name = "lineDrought.trtPinoxaden"

drt_plot <- EnhancedVolcano(toptable = resLFC_interaction_drt, lab = rownames(resLFC_interaction_drt), 
                            x = 'log2FoldChange',
                            y = "pvalue",
                            title = 'Within-line analysis:\nRecurrent Drought Selection',
                            subtitle = 'name = "lineDrought.trtPinoxaden"',
                            caption = bquote(~Log[2]~ "fold change cutoff = ±2; p-value cutoff = 0.01"),
                            xlab = bquote(~Log[2]~ 'fold change'),
                            pCutoff = 0.01,
                            FCcutoff = 2,
                            selectLab = c('Bradi2g45080.v3.2',
                                          'Bradi1g12690.v3.2',
                                          'Bradi3g52060.v3.2',
                                          'Bradi4g29246.v3.2',
                                          'Bradi3g56340.v3.2',
                                          'Bradi2g58989.v3.2',
                                          'Bradi2g17530.v3.2'),
                            cutoffLineType = 'twodash',
                            cutoffLineWidth = 0.8,
                            pointSize = 4.0,
                            labSize = 4.5,
                            labFace = 'bold',
                            boxedLabels = TRUE,
                            colAlpha = 0.5,
                            legendPosition = 'none',
                            drawConnectors = TRUE,
                            widthConnectors = 0.5,
                            colConnectors = 'black')

drt_plot

dev.off()

pdf("drt_Evolcano.pdf", 
    width = 8, height = 8)

par(mfrow = c(1,1))

drt_plot

dev.off()
