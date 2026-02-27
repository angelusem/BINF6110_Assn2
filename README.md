# Assignment 2 — Bulk RNA-seq differential expression and functional analysis on Yeast biofilm

## Introduction
## Methods

## **Results**

RNA-seq libraries from three velum developmental stages: early, thin, mature ( *n* \= 3 per stage) were quantified and summarized to gene-level counts for differential expression (DE) analysis.

Variance-stabilized expression values showed strong stage-associated structure. In performing PCA, it was observed that samples clustered by stage, with **PC1** explaining **72%** and **PC2** explaining **24%** of the variance. Early and mature samples were clearly separated along PC1, while thin samples formed a distinct cluster intermediate along PC1 but separated along PC2, consistent with a progressive transcriptome shift across development of biofilm (Figure 1).

**Figure 1\. PCA of VST-transformed gene expression by stage: global expression patterns observed to be separate by developmental stage**   
PCA was computed on VST-transformed DESeq2 expression values and colored by stage. Points here represent individual samples. Replicates from the same stage are seen to be clustering more closely to one another than to other stages.

(results/figures/PCA\_stage.png)

### **Differential expression across stages**

DESeq2 identified widespread differential expression across all contrasts (FDR \< 0.05), with the largest number and strongest effect sizes in early vs mature. Using a more stringent subset (**FDR \< 0.05 and |log2FC| ≥ 1**), early vs mature still retained the most DE genes, consistent with the greatest transcriptional remodeling occurring between the earliest and latest stages of velum biofilm development. The results for each contrast are summarized in Table 1\.

**Table 1\. Differential expression summary by contrast**

|  contrast | Number of genes tested | Number of genes with  FDR \< 0.05 | Number of genes with  FDR \< 0.05 and |log2FC| ≥ 1 | med\_absLFC\_sig  (rounded to 3 decimal places) |
| :---- | ----- | ----- | ----- | ----- |
| early vs thin   | 5957 | 2194 | 1119 | 1.009 |
| thin vs mature | 5957 | 2352 | 1320 | 1.059 |
| early vs mature | 5957 | 2968 | 1866 | 1.193 |

To visualize the distribution of effect sizes and statistical support in the early vs mature contrast, a volcano plot was generated (Figure 2). A heatmap of the 50 most significant genes by adjusted p-value further highlights coordinated, stage-associated expression patterns across samples (Figure 3).

**Figure 2\. Volcano plot for early vs mature differential expression.**  
 Each point represents a gene; x-axis is log2 fold change (early vs mature) and y-axis is −log10(FDR). Observed is broad differential expression spanning both positive and negative log2 fold changes, with many genes strongly supported after multiple-testing correction. (results/figures/Volcano\_early\_vs\_mature.png).

**Figure 3\. Heatmap of the top 50 DE genes in early vs mature contrast.**  
 Heatmap shows VST expression centered per gene for the 50 lowest-FDR, most significant,  genes. Columns are annotated by stage. Depicted is stage-specific expression patterns: early samples clustered together and were distinct from mature samples, while thin samples formed their own group with intermediate profiles across many genes (results/figures/Heatmap\_top50\_early\_vs\_mature.png)

### **Functional shifts captured by GSEA (GO Biological Process; early vs mature)**

Multiple metabolic and nucleotide-related categories showed strong **positive enrichment (NES \> 0\)**, indicating their member genes were concentrated toward the “early-high” end of the ranked list. In contrast, translation/mitochondria-associated categories showed **negative enrichment (NES \< 0\)**, indicating their genes were concentrated toward the “mature-high” end (Figure 4).

Two representative, highly significant terms were:

* **GO:0009185 – ribonucleoside diphosphate metabolic process**: NES \= **2.756**, FDR (**p.adjust**) \= **6.02×10⁻⁹**, leading-edge **21/35 genes** (GeneRatio \= **0.60**)

* **GO:0032543 – mitochondrial translation**: NES \= **−2.531**, FDR (**p.adjust**) \= **6.02×10⁻⁹**, leading-edge **84/135 genes** (GeneRatio \= **0.62**)

The enrichment curve for the top term (GO:0009185) shows a strong positive running enrichment score early in the ranked list, consistent with many leading-edge genes appearing among the most early-upregulated signals (Figure 5).

**Figure 4\. GSEA dotplot (GO BP) for early vs mature.** (GSEA\_GO\_BP\_dotplot\_early\_vs\_mature.png)  
 **Figure 5\. Enrichment plot for GO:0009185.** (GSEA\_GO\_BP\_topTerm\_curve\_early\_vs\_mature.png)

### **Genes linking GSEA terms to expression**

Representative “leading-edge” genes from the positively enriched nucleotide/metabolic term showed decreasing expression from early → thin → mature, consistent with positive early vs mature log2 fold changes. For example, **YJL052W** (log2FC \= **5.35**, padj \= **7.62×10⁻¹⁶²**) and **YCR012W** (log2FC \= **4.23**, padj \= **1.78×10⁻⁹⁰**) both decline across stages (Figure 6). Conversely, representative leading-edge genes from the negatively enriched mitochondrial translation term increased across development, consistent with negative log2 fold changes (higher in mature), including **YJL102W** (log2FC \= **−1.79**, padj \= **4.80×10⁻²⁶**) and **YML009C** (log2FC \= **−2.69**, padj \= **2.35×10⁻¹⁹**) (Figure 6).

**Figure 6\. VST expression trajectories for selected leading-edge genes across stages.** (chosenGenesFA\_stage.png)

## Discussion
## References 


