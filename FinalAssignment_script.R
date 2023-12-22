# Yuehui He Assignment Script

# Change working directory.
path  = "C:/Users/heyue/Desktop/R/bio" 

file_name = "brca_tcga_pan_can_atlas_2018.tar.gz"

# extract the files
untar(file_name)

# change directory to the extracted folders
setwd(paste(getwd() , "/brca_tcga_pan_can_atlas_2018", sep = ""))

# read the files
clinical = read.delim("data_clinical_patient.txt")

rnaseq = read.delim("data_mrna_seq_v2_rsem.txt")

# delete the genes for which there's more than one Hugo Symbol
# These are typically genes with no Hugo Symbol ("" as an entry) or pseudogenes.

keep = !duplicated(rnaseq[,1])

rnaseq = rnaseq[keep,]

# set rownames of rnaseq to hugo symbols

rownames(rnaseq)  = rnaseq[,1]

# read CNA Data

cna = read.delim('data_cna.txt')

# find ERBB2 in cna

erbb2_indx = which(cna[,1] == 'ERBB2')

# match patients in rnaseq to patients in cna.

rna_cna_id = which(is.element(colnames(rnaseq[,-c(1,2)]), colnames(cna[,-c(1,2)])))

# select only the rna cases which have cna data.

rna_cna_sub = rnaseq[,2+rna_cna_id]

# pre-allocate memory for ERBB2

meta_erbb2 = matrix(0,length(rna_cna_id),1)
for (i in 1:length(rna_cna_id)){
  # access the colnames of i
  col_i = colnames(rna_cna_sub)[i]
  # get the index in cna for the same patient
  col_cna = which(colnames(cna)==col_i)
  # store if they're amplified.
  meta_erbb2[i,] = 1*(cna[erbb2_indx,col_cna]>0)
}


# add a title to the metadata.

colnames(meta_erbb2) = 'ERBB2Amp'

# transform into integers

rna_cna_sub = round(rna_cna_sub)



library(DESeq2)


# convert meta_erbb2 to a dataframe
meta_erbb2_df <- data.frame(ERBB2Amp = factor(meta_erbb2))

# perform normalization
dds <- DESeqDataSetFromMatrix(countData = rna_cna_sub, colData = meta_erbb2_df, design = ~ ERBB2Amp)
dds <- DESeq(dds)
res <- results(dds)

summary(res)


# get top 10 genes based on fold change
topGenes <- head(res[order(-abs(res$log2FoldChange)), ], 10)
print(topGenes)

# get the variance stabilised transformed expression values.
vsd <- vst(dds, blind=FALSE)
vst_values <- assay(vsd)

# perform PCA
pca_data <- prcomp(t(vst_values))

plot(pca_data$x[,1], pca_data$x[,2], 
     col = as.factor(meta_erbb2_df$ERBB2Amp),
     pch = 16, 
     main = "PCA Plot with VST Values",
     xlab = "PC1",
     ylab = "PC2",
     cex = 0.7
)
legend("topright", legend = levels(as.factor(meta_erbb2_df$ERBB2Amp)), col = 1:2, pch = 16, title = "ERBB2Amp")



# Pathway Enrichment Analysis
library(clusterProfiler)

# get significantly expressed gene
signif = which(res$padj<0.05)
deg = res[signif,]
# separate them 
dup = deg[deg[,2]>0.,]
ddown = deg[deg[,2]<0.,]

entrez_ids = rnaseq[keep,2]

entrez_all = entrez_ids[signif]
entrez_up = entrez_all[signif[deg[,2]>0.]]
entrez_down = entrez_all[signif[deg[,2]<0.]]

all_paths =   enrichKEGG(gene = entrez_all, organism = 'hsa', pvalueCutoff = 0.05)
head(all_paths)

up_paths = enrichKEGG(gene = entrez_up, organism = 'hsa', pvalueCutoff = 0.05)
head(up_paths)



# Heatmap

install.packages("pheatmap")
library(pheatmap)

signif2 = which(res$padj<0.05 & abs(res$log2FoldChange)>1)
deg2 = res[signif2,]
deg_names <- rownames(deg2)
heatmap_data <- res[deg_names,]

pheatmap(heatmap_data,
        col = colorRampPalette(c("blue", "white", "red"))(50),
        scale = "row",
        main = "Heatmap",
        cluster_cols = FALSE,
        show_colnames = TRUE,
        show_rownames = FALSE) 


# volcano plot for all differentially expressed genes
plot(res$log2FoldChange, -log10(res$pvalue), col = ifelse(res$padj < 0.05, "red", "black"), main = "Volcano Plot", xlab="Log2 Fold Change",ylab="-Log10 (p-value)")

# highlight up-regulated genes in red
points(dup$log2FoldChange, -log10(dup$pvalue), col = "red", pch = 16)

# highlight down-regulated genes in blue
points(ddown$log2FoldChange, -log10(ddown$pvalue), col = "blue", pch = 16)

# add labels and legend
legend("topleft", legend = c("All Genes", "Up-regulated", "Down-regulated"), col = c("black", "red", "blue"), pch = 16)
threshold_x <- 1
threshold_y <- 200
text(dup$log2FoldChange[dup$log2FoldChange > threshold_x & -log10(dup$pvalue) > threshold_y],
     -log10(dup$pvalue[dup$log2FoldChange > threshold_x & -log10(dup$pvalue) > threshold_y]),
     labels = rownames(dup)[dup$log2FoldChange > threshold_x & -log10(dup$pvalue) > threshold_y],
     pos = 1, col = "blue", cex = 0.8)



# perform Clustering

# calculate distance matrix
dist_matrix <- dist(rna_cna_sub)

# perform hierarchical clustering
hc <- hclust(dist_matrix)

# cut the dendrogram to form clusters
clusters <- cutree(hc, k = 4)

# PCA
pca_result <- prcomp(rna_cna_sub)

# extract PC scores
pc_scores <- pca_result$x

# create a data frame with PC scores and cluster assignments
pca_df <- data.frame(PC1 = pc_scores[, 1], PC2 = pc_scores[, 2], Cluster = as.factor(clusters))

# Plot PCA
plot(pca_df$PC1, pca_df$PC2, col = as.factor(pca_df$Cluster), pch = 16, 
     main = "PCA Plot with Clustering", xlab = "PC1", ylab = "PC2")

# add legend
legend("topright", legend = levels(pca_df$Cluster), col = 1:4, pch = 16, title = "Cluster")












