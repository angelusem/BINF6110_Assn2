#ASSN2 BINF6110: DE in yeast biofilm

#loading required packages
library(readr)
library(dplyr)
library(tximport)
library(DESeq2)
library(GenomicFeatures)
library(AnnotationDbi)
library(ggplot2)
library(pheatmap)


# inputting metadata
meta <- read_tsv("metadata/metadata.tsv", show_col_types = FALSE)
meta <- meta %>%
  mutate(stage = factor(stage, levels = c("early","thin","mature")))

# pointing to Salmon quant.sf files
files <- file.path("results/quants", meta$srr, "quant.sf")
names(files) <- meta$sample_id

stopifnot(all(file.exists(files)))

# building tx2gene from the GTF. (specifiying namespaces to make sure appropriate packages are used)

txdb <- GenomicFeatures::makeTxDbFromGFF("ref/genes.gtf")

k <- AnnotationDbi::keys(txdb, keytype = "TXNAME")

tx2gene <- AnnotationDbi::select(txdb, keys = k, keytype = "TXNAME", columns = "GENEID") |>
  dplyr::distinct(TXNAME, GENEID) |>
  dplyr::rename(tx = TXNAME, gene = GENEID)

head(tx2gene)

# importing transcript counts and summarize to genes
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

dds <- DESeqDataSetFromTximport(txi, colData = as.data.frame(meta), design = ~ stage)
dds <- DESeq(dds)

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# making pairwise contrasting comparisons
res_early_thin   <- results(dds, contrast = c("stage","early","thin"))
res_thin_mature  <- results(dds, contrast = c("stage","thin","mature"))
res_early_mature <- results(dds, contrast = c("stage","early","mature"))

write.csv(as.data.frame(res_early_thin),   "results/tables/DE_early_vs_thin.csv")
write.csv(as.data.frame(res_thin_mature),  "results/tables/DE_thin_vs_mature.csv")
write.csv(as.data.frame(res_early_mature), "results/tables/DE_early_vs_mature.csv")

# Figure 1: PCA (for overall structure)
vsd <- vst(dds, blind = FALSE)
pca <- plotPCA(vsd, intgroup = "stage", returnData = TRUE)
percentVar <- round(100 * attr(pca, "percentVar"))

p_pca <- ggplot(pca, aes(PC1, PC2, color = stage, label = name)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw()

ggsave("results/figures/PCA_stage.png", p_pca, width = 6, height = 4, dpi = 300)

# Figure 2: Volcano plot with the top DE genes for one contrast
volc <- as.data.frame(res_early_mature) %>%
  mutate(gene = rownames(.),
         neglog10padj = -log10(padj))

p_volc <- ggplot(volc, aes(x = log2FoldChange, y = neglog10padj)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  xlab("log2 fold change (early vs mature)") +
  ylab("-log10 adjusted p-value")

ggsave("results/figures/Volcano_early_vs_mature.png", p_volc, width = 6, height = 4, dpi = 300)

# Figure 3: heatmap of top genes
top <- head(order(res_early_mature$padj), 50)
mat <- assay(vsd)[top, ]
mat <- mat - rowMeans(mat)
pheatmap(mat, annotation_col = data.frame(stage = meta$stage, row.names = meta$sample_id),
         filename = "results/figures/Heatmap_top50_early_vs_mature.png",
         width = 7, height = 8)
