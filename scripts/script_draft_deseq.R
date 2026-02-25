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
  library(clusterProfiler)
  library(org.Sc.sgd.db)
  library(enrichplot)

## PART 1: DE Analysis using DeSeq2 ----

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
# generating DE strength summary between contrasts made

summ<-function(res){
  df<-as.data.frame(res)
  df<-df[!is.na(df$padj), ]
  c(
    n_tested=nrow(df),
    n_FDR05=sum(df$padj<0.05),
    n_FDR05_LFC1=sum(df$padj<0.05&abs(df$log2FoldChange)>=1),
    med_absLFC_sig=median(abs(df$log2FoldChange[df$padj<0.05]),na.rm=TRUE)
  )
}

rbind(
  early_thin=summ(res_early_thin),
  thin_mature=summ(res_thin_mature),
  early_mature=summ(res_early_mature)
)

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
# Is early vs mature the largest shift? 
# checking centroid distances between groups to verify

pca$stage<-meta$stage[match(pca$name,meta$sample_id)]

centroids<-aggregate(cbind(PC1,PC2)~stage,pca,mean)

dist_early_thin<-with(centroids,sqrt((PC1[stage=="early"]-PC1[stage=="thin"])^2+
                                       (PC2[stage=="early"]-PC2[stage=="thin"])^2))

dist_thin_mature<-with(centroids,sqrt((PC1[stage=="thin"]-PC1[stage=="mature"])^2+
                                        (PC2[stage=="thin"]-PC2[stage=="mature"])^2))

dist_early_mature<-with(centroids,sqrt((PC1[stage=="early"]-PC1[stage=="mature"])^2+
                                         (PC2[stage=="early"]-PC2[stage=="mature"])^2))

c(early_thin=dist_early_thin,thin_mature=dist_thin_mature,early_mature=dist_early_mature)

# Figure 2: Volcano plot with the top DE genes for early vs mature contrast 
volc <- as.data.frame(res_early_mature) %>%
  mutate(gene = rownames(.),
         neglog10padj = -log10(padj))

p_volc <- ggplot(volc, aes(x = log2FoldChange, y = neglog10padj)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  xlab("log2 fold change (early vs mature)") +
  ylab("-log10 adjusted p-value")

ggsave("results/figures/Volcano_early_vs_mature.png", p_volc, width = 6, height = 4, dpi = 300)

# Figure 3: heatmap of top genes in early vs mature contrast
top <- head(order(res_early_mature$padj), 50)
mat <- assay(vsd)[top, ]
mat <- mat - rowMeans(mat)
pheatmap(mat, annotation_col = data.frame(stage = meta$stage, row.names = meta$sample_id),
         filename = "results/figures/Heatmap_top50_early_vs_mature.png",
         width = 7, height = 8)
## PART 2: Functional Enrichment Analysis using GSEA ----
# building my ranked list for GO GSEA enrichment analysis, from my DeSeq2 output
df <- as.data.frame(res_early_mature) %>%
  mutate(gene = rownames(.)) %>%
  filter(!is.na(stat))

geneList <- df$stat
names(geneList) <- df$gene
geneList <- sort(geneList, decreasing = TRUE)

# running GSEA
gsea_bp <- gseGO(
  geneList = geneList,
  OrgDb    = org.Sc.sgd.db,
  keyType  = "ORF",
  ont      = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  seed = TRUE
)

dir.create("results/tables", recursive=TRUE, showWarnings=FALSE)
write.csv(as.data.frame(gsea_bp), "results/tables/GSEA_GO_BP_early_vs_mature.csv", row.names=FALSE)
# plotting a dotplot
p1 <- dotplot(gsea_bp, showCategory = 20) + ggtitle("GSEA GO BP: early vs mature")
ggsave("results/figures/GSEA_GO_BP_dotplot_early_vs_mature.png", p1, width=8, height=6, dpi=300)
# plotting an enrichment curve for the top term

p2 <- gseaplot2(gsea_bp, geneSetID = 1, title = gsea_bp@result$Description[1])
ggsave("results/figures/GSEA_GO_BP_topTerm_curve_early_vs_mature.png", p2, width=8, height=5, dpi=300)

# making a table for results
gsea_tbl <- as.data.frame(gsea_bp@result)[, c("Description","NES","p.adjust","setSize","core_enrichment")]
head(gsea_tbl, 10)

# connecting pathways back to actual genes
top1 <- gsea_bp@result[1, ]
core <- unlist(strsplit(top1$core_enrichment, "/"))

core_df <- as.data.frame(res_early_mature)[core, c("log2FoldChange","padj","baseMean")]
core_df$gene <- rownames(core_df)
core_df <- core_df[order(core_df$padj), ]
head(core_df, 10)
# to get gene names and descriptions for discussion
AnnotationDbi::select(org.Sc.sgd.db,
                      keys = core_df$gene[1:10],
                      keytype = "ORF",
                      columns = c("GENENAME","DESCRIPTION"))
