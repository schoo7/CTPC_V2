### CTPC V3 RNAseq analysis
library(Rtsne)
library(plotly)
library(ggplot2)
library(preprocessCore)
library(sva)
library(limma)
library(viridisLite)
library(ggrepel)
library(dplyr)
library(htmlwidgets)
library(tidyr)
library(Seurat)
library(singscore)
library(GSEABase)
library(ExperimentHub)
library(msigdb)
library(pheatmap)
library(reshape2)
library(biomaRt)
library(ggforce)
library(decoupleR)
library(SCpubr)
library(vcfR)
library(SNPRelate)
#
setwd("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/CTPC_V3")
# Define the main directory containing the subfolders
main_dir <- "/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/CTPC_V3/rnaseq"
# List of subfolders
subfolders <- paste0(main_dir, "/CTPC_result_", 1:32)
# Initialize a list to store data frames
tpm_list <- list()
# Function to read and store each data frame
read_tpm_files <- function(subfolder) {
  file_path <- file.path(subfolder, "star_salmon/salmon.merged.gene_tpm.tsv")
  df <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  return(df)
}
# Read all the files and store them in tpm_list with appropriate names
for (i in 1:length(subfolders)) {
  tpm_list[[i]] <- read_tpm_files(subfolders[i])
  names(tpm_list)[i] <- paste0("TPM_", i)
}
# Check if the first two columns are identical
identical_columns <- function(df_list) {
  reference <- df_list[[1]][, 1:2]
  for (i in 2:length(df_list)) {
    if (!all(reference == df_list[[i]][, 1:2])) {
      return(FALSE)
    }
  }
  return(TRUE)
}
# Combine data frames if the first two columns are identical
if (identical_columns(tpm_list)) {
  combined_tpm <- tpm_list[[1]][, 1:2]  # Start with gene_id and gene_name
  for (i in 1:length(tpm_list)) {
    combined_tpm <- cbind(combined_tpm, tpm_list[[i]][, -c(1, 2)])  # Append sample columns
  }
  sample_names <- unlist(lapply(tpm_list, function(x) names(x)[-c(1, 2)]))
  names(combined_tpm) <- c("gene_id", "gene_name", sample_names)
  TPM <- combined_tpm
} else {
  stop("The first two columns (gene_id and gene_name) are not identical in all files.")
}

index=duplicated(TPM$gene_name)
TPM=TPM[!index,]
rownames(TPM)=TPM$gene_name
TPM=TPM[,-c(1:2)]
saveRDS(TPM,"CTPC_TPM_original.rds")
#TPM=readRDS("CTPC_TPM_original.rds")
# combine with meta
meta=read.csv("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/CTPC_V3/CTPC_meta.csv")
rownames(meta)=meta$ID
meta=meta[colnames(TPM),]
table(meta$group)
# Clean up TPM and log2-transform
TPM_nonzero = TPM[rowSums(TPM != 0) > 0, ]
TPM_log2 = log2(TPM_nonzero + 1)  # Add 1 to avoid log(0)
# Select the most variable genes
gene_variance = apply(TPM_log2, 1, var)
top_genes = names(sort(gene_variance, decreasing = TRUE))[1:5000]  # Select top 5000 most variable genes
TPM_top_var_genes = TPM_log2[top_genes, ]
# Initial PCA analysis with the most variable genes
matrix = as.matrix(t(TPM_top_var_genes))
pca_result = prcomp(matrix, center = TRUE, scale. = TRUE)
# Extract PCA results
pca_df = data.frame(Sample = rownames(pca_result$x), 
                    PC1 = pca_result$x[,1], 
                    PC2 = pca_result$x[,2])
# Merge with metadata
pca_df = merge(pca_df, meta, by.x = "Sample", by.y = "row.names")
# Plot the PCA results
ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 1) +
  labs(title = "PCA of Log2-Transformed RNAseq Data (Top Variable Genes)", x = "PC1", y = "PC2") +
  theme_minimal()

# Run t-SNE on the PCA results (using the top principal components)
set.seed(123)  # For reproducibility
tsne_result = Rtsne(pca_result$x[, 1:50], pca=F,dims = 2, perplexity = 30, verbose = TRUE, max_iter = 1000)
# Ensure the rownames are retained
tsne_df = data.frame(Sample = rownames(matrix), 
                     tSNE1 = tsne_result$Y[,1], 
                     tSNE2 = tsne_result$Y[,2])
# Merge with metadata
tsne_df = merge(tsne_df, meta, by.x = "Sample", by.y = "row.names")
# Plot the 2D t-SNE results
ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = group)) +
  geom_point(size = 1) +
  labs(title = "2D t-SNE of PCA Results", x = "t-SNE 1", y = "t-SNE 2") +
  theme_minimal()
# Run 3D t-SNE on the PCA results
tsne_result_3d = Rtsne(pca_result$x[, 1:50], pca=FALSE, dims = 3, perplexity = 30, verbose = TRUE, max_iter = 1000)
# Ensure the rownames are retained
tsne_3d_df = data.frame(Sample = rownames(matrix), 
                        tSNE1 = tsne_result_3d$Y[,1], 
                        tSNE2 = tsne_result_3d$Y[,2], 
                        tSNE3 = tsne_result_3d$Y[,3])
# Merge with metadata
tsne_3d_df = merge(tsne_3d_df, meta, by.x = "Sample", by.y = "row.names")
# Plot the 3D t-SNE results with adjusted point size
plot_ly(tsne_3d_df, x = ~tSNE1, y = ~tSNE2, z = ~tSNE3, color = ~group, colors = c('#1f77b4', '#ff7f0e', '#2ca02c')) %>%
  add_markers(marker = list(size = 3)) %>%  # Change the size to your desired value
  layout(title = "3D t-SNE of PCA Results",
         scene = list(xaxis = list(title = 't-SNE 1'),
                      yaxis = list(title = 't-SNE 2'),
                      zaxis = list(title = 't-SNE 3')))

### samples are still scattered, now remove low quality samples
# cut1 duplication level < 70%
dup=read.delim("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/CTPC_V3/rnaseq/multiqc_data/picard_deduplication.txt",header = T)
dup$dup_rate_unpaired=dup$Duplicate.Unpaired/(dup$Duplicate.Unpaired+dup$Unique.Unpaired)
dup$dup_rate_paired=dup$Duplicate.Pairs.Nonoptical/(dup$Duplicate.Pairs.Nonoptical+dup$Unique.Pairs)
dup_keep=c(subset(dup,(dup$Unique.Pairs==0)&(dup$dup_rate_unpaired<0.7))$Sample,
           subset(dup,(dup$Unique.Pairs!=0)&(dup$dup_rate_paired<0.7))$Sample)
# cut2 star alignment number of unique reads > 10M and >70%
star=read.delim("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/CTPC_V3/rnaseq/multiqc_data/star_alignment_plot.txt",header = T)
star$rate=star$Uniquely.mapped/(star$Uniquely.mapped+star$Mapped.to.multiple.loci+star$Mapped.to.too.many.loci+
                                  star$Unmapped..too.short+star$Unmapped..other)
star_keep=subset(star,(star$Uniquely.mapped>10000000)&(star$rate>0.7))$Sample

### start to filter samples
keep=intersect(dup_keep,star_keep)
TPM_clean=TPM[,keep]
meta_clean=meta[keep,]
### refine the meta information
meta_clean$group=gsub("LNCaP_42D","LNCaP-42D",meta_clean$group)
meta_clean$group=factor(meta_clean$group,levels = c("LNCaP","LNCaP-abl","LNCaP-95","C4-2","C4-2B","VCaP","22RV1","PC3","DU145","H660","LASCPC1","LNCaP-42D"))
saveRDS(meta_clean,"meta_clean.rds")
saveRDS(TPM_clean,"TPM_clean.rds")
### now re-do the PCA analysis
# Clean up TPM and log2-transform
TPM_nonzero = TPM_clean[rowSums(TPM_clean != 0) > 0, ]
# Select the most variable genes
gene_variance = apply(TPM_nonzero, 1, var)
top_genes = names(sort(gene_variance, decreasing = TRUE))[1:5000]  # Select top 5000 most variable genes
TPM_top_var_genes = TPM_nonzero[top_genes, ]
# Initial PCA analysis with the most variable genes
matrix = as.matrix(t(TPM_top_var_genes))
pca_result = prcomp(matrix, center = TRUE, scale. = TRUE)
# Extract PCA results
pca_df = data.frame(Sample = rownames(pca_result$x), 
                    PC1 = pca_result$x[,1], 
                    PC2 = pca_result$x[,2])
# Merge with metadata
pca_df = merge(pca_df, meta_clean,by.x = "Sample", by.y = "ID")
# Plot the PCA results
ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 1) +
  labs(title = "PCA of Log2-Transformed RNAseq Data (Top Variable Genes)", x = "PC1", y = "PC2") +
  theme_minimal()
# Run t-SNE on the PCA results (using the top principal components)
set.seed(123)  # For reproducibility
tsne_result = Rtsne(pca_result$x[, 1:50], pca=F,dims = 2, perplexity = 30, verbose = TRUE, max_iter = 1000)
# Ensure the rownames are retained
tsne_df = data.frame(Sample = rownames(matrix), 
                     tSNE1 = tsne_result$Y[,1], 
                     tSNE2 = tsne_result$Y[,2])
# Merge with metadata
tsne_df = merge(tsne_df, meta, by.x = "Sample", by.y = "row.names")
# Plot the 2D t-SNE results
ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = group)) +
  geom_point(size = 1) +
  labs(title = "2D t-SNE of PCA Results", x = "t-SNE 1", y = "t-SNE 2") +
  theme_minimal()
# Run 3D t-SNE on the PCA results
tsne_result_3d = Rtsne(pca_result$x[, 1:50], pca=FALSE, dims = 3, perplexity = 30, verbose = TRUE, max_iter = 1000)
# Ensure the rownames are retained
tsne_3d_df = data.frame(Sample = rownames(matrix), 
                        tSNE1 = tsne_result_3d$Y[,1], 
                        tSNE2 = tsne_result_3d$Y[,2], 
                        tSNE3 = tsne_result_3d$Y[,3])
# Merge with metadata
tsne_3d_df = merge(tsne_3d_df, meta_clean, by.x = "Sample", by.y = "ID")
# Plot the 3D t-SNE results with adjusted point size
plot_ly(tsne_3d_df, x = ~tSNE1, y = ~tSNE2, z = ~tSNE3, color = ~group, colors = c('#1f77b4', '#ff7f0e', '#2ca02c')) %>%
  add_markers(marker = list(size = 3)) %>%  # Change the size to your desired value
  layout(title = "3D t-SNE of PCA Results",
         scene = list(xaxis = list(title = 't-SNE 1'),
                      yaxis = list(title = 't-SNE 2'),
                      zaxis = list(title = 't-SNE 3')))

### try to remove batch effect
TPM=readRDS("TPM_clean.rds")
meta_clean=readRDS("meta_clean.rds")
identical(colnames(TPM),rownames(meta_clean))
# Clean up TPM and log2-transform
TPM_nonzero = TPM[rowSums(TPM != 0) > 0, ]
TPM_log2 = log2(TPM_nonzero + 1)
TPM_log2 = normalize.quantiles(as.matrix(TPM_log2),keep.names = T)
meta_clean$Batch=gsub("\t", "", meta_clean$Batch)


### the PCA plot showed the NEPCa cells H660, LASCPC1 and LNCaP-42D are far from AdPCa cells
### I used the batch effect removal on all samples and it generated artifical effect
### To overcome this, I splite the samples into NEPCa and AdPCa then performed the batch effect removal and the data is meaningful
# Split the data into two groups
group_1 = c("H660", "LASCPC1", "LNCaP-42D")
idx_group_1 = meta_clean$group %in% group_1
TPM_log2_1 = TPM_log2[, idx_group_1]
TPM_log2_2 = TPM_log2[, !idx_group_1]
# Split the meta_clean data correspondingly
meta_clean_1 = meta_clean[idx_group_1, ]
meta_clean_2 = meta_clean[!idx_group_1, ]
# Perform batch effect removal on each group separately
TPM_cor_1 = removeBatchEffect(TPM_log2_1, batch = meta_clean_1$Batch, group = meta_clean_1$group)
TPM_cor_2 = removeBatchEffect(TPM_log2_2, batch = meta_clean_2$Batch, group = meta_clean_2$group)
# Combine the batch-corrected data back into a single dataset
TPM_cor = cbind(TPM_cor_1, TPM_cor_2)
# Reorder the columns to match the original order
TPM_cor=as.data.frame(TPM_cor)
TPM_cor = TPM_cor[, colnames(TPM)]
#TPM_cor=removeBatchEffect(TPM_log2, batch=meta_clean$Batch, group=meta_clean$group)
#TPM_cor=as.data.frame(TPM_cor)
#rownames(TPM_cor)=rownames(TPM_nonzero)
#colnames(TPM_cor)=colnames(TPM_nonzero)
### save batch corrected data, the TPM data has been log2 transformed
saveRDS(TPM_cor,"CTPC_TPM_batchCorrected.rds")
saveRDS(meta_clean,"CTPC_meta_cleaned.rds")




### now re-do the PCA analysis
# Select the most variable genes
TPM_cor=readRDS("CTPC_TPM_batchCorrected.rds")
meta_clean=readRDS("CTPC_meta_cleaned.rds")
gene_variance = apply(TPM_cor, 1, var)
top_genes = names(sort(gene_variance, decreasing = TRUE))[1:5000]  # Select top 5000 most variable genes
TPM_top_var_genes = TPM_cor[top_genes, ]
# Initial PCA analysis with the most variable genes
matrix = as.matrix(t(TPM_top_var_genes))
pca_result = prcomp(matrix, center = TRUE, scale. = TRUE)
# Extract PCA results
pca_df = data.frame(Sample = rownames(pca_result$x), 
                    PC1 = pca_result$x[,1], 
                    PC2 = pca_result$x[,2])
# Merge with metadata
pca_df = merge(pca_df, meta_clean, by.x = "Sample", by.y = "ID")
# Plot the PCA results
ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 0.1) +
  labs(title = "PCA of Log2-Transformed RNAseq Data (Top Variable Genes)", x = "PC1", y = "PC2") +
  theme_classic() +
  scale_color_hue(l = 50)
ggsave("CTPC_corrected_PCA.pdf",width = 24,height = 16,units = "cm")
# Run t-SNE on the PCA results (using the top principal components)
set.seed(123)  # For reproducibility
tsne_result = Rtsne(pca_result$x[, 1:50], pca=F,dims = 2, perplexity = 30, verbose = TRUE, max_iter = 1000)
# Ensure the rownames are retained
tsne_df = data.frame(Sample = rownames(matrix), 
                     tSNE1 = tsne_result$Y[,1], 
                     tSNE2 = tsne_result$Y[,2])
# Merge with metadata
tsne_df = merge(tsne_df, meta_clean, by.x = "Sample", by.y = "ID")
# Plot the 2D t-SNE results
ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = group)) +
  geom_point(size = 0.5) +
  labs(title = "PCA of Log2-Transformed RNAseq Data (Top Variable Genes)", x = "PC1", y = "PC2") +
  theme_classic() +
  scale_color_hue(l = 50)
ggsave("CTPC_corrected_2Dtsne.pdf",width = 18,height = 18,units = "cm")
# Run 3D t-SNE on the PCA results
tsne_result_3d = Rtsne(pca_result$x[, 1:50], pca=FALSE, dims = 3, perplexity = 30, verbose = TRUE, max_iter = 1000)
# Ensure the rownames are retained
tsne_3d_df = data.frame(Sample = rownames(matrix), 
                        tSNE1 = tsne_result_3d$Y[,1], 
                        tSNE2 = tsne_result_3d$Y[,2], 
                        tSNE3 = tsne_result_3d$Y[,3])
# Merge with metadata
tsne_3d_df = merge(tsne_3d_df,meta_clean, by.x = "Sample", by.y = "ID")
# Plot the 3D t-SNE results with adjusted point size and custom color map
color_palette <- scale_color_hue(l = 50)$palette(12)
group_colors <- setNames(color_palette, unique(tsne_3d_df$group))
plot_ly(tsne_3d_df, x = ~tSNE1, y = ~tSNE2, z = ~tSNE3, color = ~group, colors = group_colors) %>%
  add_markers(marker = list(size = 2)) %>%  # Change the size to your desired value
  layout(title = "3D t-SNE of PCA Results",
         scene = list(xaxis = list(title = 't-SNE 1', showgrid = FALSE),  # Remove x-axis grid
                      yaxis = list(title = 't-SNE 2', showgrid = FALSE),  # Remove y-axis grid
                      zaxis = list(title = 't-SNE 3', showgrid = FALSE)))  # Remove z-axis grid
# Calculate centroids for each group to position the labels
centroids <- tsne_3d_df %>%
  group_by(group) %>%
  summarize(tSNE1 = mean(tSNE1), tSNE2 = mean(tSNE2), tSNE3 = mean(tSNE3))

# Adjust label positions if needed (e.g., slightly offset from the centroid)
centroids <- centroids %>%
  mutate(tSNE1_label = tSNE1 + 0.5,
         tSNE2_label = tSNE2 + 0.5,
         tSNE3_label = tSNE3 + 0.5)

# Plot the 3D t-SNE results with group labels
p <- plot_ly(tsne_3d_df, x = ~tSNE1, y = ~tSNE2, z = ~tSNE3, color = ~group, colors = group_colors) %>%
  add_markers(marker = list(size = 2)) %>%  # Change the size to your desired value
  add_text(data = centroids, x = ~tSNE1_label, y = ~tSNE2_label, z = ~tSNE3_label, text = ~group, 
           textfont = list(size = 20), showlegend = FALSE) %>%
  layout(title = "3D t-SNE of PCA Results",
         scene = list(xaxis = list(title = 't-SNE 1', showgrid = FALSE),  # Remove x-axis grid
                      yaxis = list(title = 't-SNE 2', showgrid = FALSE),  # Remove y-axis grid
                      zaxis = list(title = 't-SNE 3', showgrid = FALSE)))  # Remove z-axis grid
p




# now combine the CTPC with ProAtlas to determine the relative clinical status of prostate cancer cell lines
ctpc=readRDS("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/CTPC_V3/CTPC_TPM_batchCorrected.rds")
ctpc_meta=readRDS("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/CTPC_V3/CTPC_meta_cleaned.rds")
proatlas=readRDS("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/HuPSA/HuPSA_integration/Proatlas.rds")
# the proatlas TPM was not log2 transformed
# convert proatlas TPM to log2+1
proatlas <- proatlas %>%
  pivot_wider(
    names_from = Gene,  # Column to become new column names
    values_from = TPM   # Column with values to fill the new columns
  )
proatlas=as.data.frame(proatlas)
rownames(proatlas)=proatlas$ID
proatlas=proatlas[,-1]
proatlas_meta=proatlas[,1,drop=F]
proatlas=proatlas[,-1]
proatlas = log2(proatlas + 1)
proatlas=t(proatlas)
proatlas=as.data.frame(proatlas)
# the CTPC was batch corrected, it generated values lower than 0
ctpc_meta_new=ctpc_meta[,2,drop=F]
ctpc_meta_new$source="cellline"
proatlas_meta_new=proatlas_meta
colnames(proatlas_meta_new)="group"
proatlas_meta_new$source="patient"
# now combine ctpc and proatlas and perform batch effect removal again
# it seems the batch effect removal didn't help
# keep the original and do the PCA analysis
gene=intersect(rownames(ctpc),rownames(proatlas))
saveRDS(gene,"CTPC_ProAtlas_commen_genes.rds")
# now only select the common genes
data=cbind(ctpc[gene,],proatlas[gene,])
meta=rbind(ctpc_meta_new,proatlas_meta_new)
# now start to batch removal to remove difference between in vivo and in vitro samples
meta=meta[colnames(data),]
meta$ID=rownames(meta)


#####
#skip this part because the markers has been generated
# the PCA analysis will not work well for the combination of cell line and patient samples
# I will use ssGSEA to run single-cell derived sigantures
# generate single-cell cancer population marker gene lists
#hupsa=readRDS("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/HuPSA/HuPSA_integration/data_annotiated_refined_pathwayscored_renamed.rds")
#hupsa$refined_histo=hupsa$histo
#hupsa$refined_histo[which(hupsa$sample=="S2")]="KRT7"
#hupsa$refined_histo[which(hupsa$sample=="S5")]="NEPC"
#hupsa$refined_histo[which(hupsa$sample=="S6")]="NEPC"
#hupsa$refined_histo[which(hupsa$sample=="S14")]="NEPC"
#hupsa$refined_histo[which(hupsa$sample=="S16")]="NEPC"
#hupsa$refined_histo[which(hupsa$sample=="S12")]="Progenitor"
#hupsa$refined_histo[which(hupsa$sample=="S19")]="Progenitor"
#hupsa$refined_histo=gsub("mCRPC","CRPC",hupsa$refined_histo)
# combine all non-cancer samples
#hupsa$refined_histo=gsub("Normal_adj","Normal",hupsa$refined_histo)
#hupsa$refined_histo=gsub("Benign","Normal",hupsa$refined_histo)
#table(hupsa$refined_histo)
# do not include normal samples
#hupsa$cell_type4=Idents(hupsa)
#Idents(hupsa)=hupsa$refined_histo
#hupsa=subset(x=hupsa,idents=c("Normal","PCa_Cribriform"),invert=T)
#Idents(hupsa)=hupsa$cell_type4
#table(Idents(hupsa))
#hupsa=subset(x=hupsa,idents=c("AdPCa_AR+_1","AdPCa_AR+_2","AdPCa_ARhi","AdPCa_ARlo","KRT7","Progenitor_like","NEPCa"))
#markers=FindAllMarkers(hupsa, only.pos = TRUE,logfc.threshold = 1,min.pct = 0.7)
#markers %>%
#  group_by(cluster) %>%
#  dplyr::filter(avg_log2FC > 1.2) %>%
#  slice_head(n = 50) %>%
#  ungroup() -> top50
#top50=as.data.frame(top50)
#gene_list_by_cluster <- top50 %>%
#  split(.$cluster) %>%                      # Split the data frame by 'cluster' into a named list
#  lapply(function(df) df$gene) 
#gene_list_by_cluster <- gene_list_by_cluster[!names(gene_list_by_cluster) %in% "AdPCa_ARlo"]
#gene_list_by_cluster <- gene_list_by_cluster[!names(gene_list_by_cluster) %in% "AdPCa_AR+_2"]
#gene_list_by_cluster <- gene_list_by_cluster[!names(gene_list_by_cluster) %in% "AdPCa_ARhi"]
#saveRDS(gene_list_by_cluster,"Cell_type_marker_from_hupsa.rds")
# now start Signature score
gene_list_by_cluster=readRDS("Cell_type_marker_from_hupsa.rds")
rankData=rankGenes(as.matrix(data))
score=multiScore(rankData = rankData,upSetColc = gene_list_by_cluster)
score=as.data.frame(t(score$Scores))
score$ID=rownames(score)
score=merge(score,meta,by="ID")
score=melt(score,id=c("ID","group","source"))
ggplot(score, aes(x = group, y = value, color = group, fill = group)) +
  geom_boxplot(outliers = FALSE,width=0.3) +
  geom_jitter(alpha =0.2,size=0.3) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_hue(l = 30) +
  facet_wrap(~ variable, ncol = 1)
ggsave("CTPC_celltype_score.pdf",width = 6,height = 9)
# now perform more signatures using gsea base
eh = ExperimentHub()
query(eh , 'msigdb')
msigdb.hs = getMsigdb(org = 'hs', id = 'SYM', version = '7.4')
listSubCollections(msigdb.hs)
# now perform the hallmark first
msigdb.hs_h=subsetCollection(msigdb.hs, 'h')
gsea=multiScore(rankData = rankData,upSetColc = msigdb.hs_h)
gsea=as.data.frame(t(as.data.frame(gsea$Scores)))
gsea$ID=rownames(gsea)
gsea=merge(gsea,meta,by="ID")
gsea=melt(gsea,id=c("ID","group","source"))
ggplot(gsea, aes(x = group, y = value, color = group, fill = group)) +
  geom_boxplot(outliers = FALSE,width=0.3) +
  geom_jitter(alpha =0.2,size=0.3) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_hue(l = 30) +
  facet_wrap(~ variable, ncol = 3)
ggsave("CTPC_gsea_score.pdf",width = 18,height = 40)
### also generate heatmap
gsea_wide <- gsea %>%
  pivot_wider(names_from = variable, values_from = value)
gsea_sorted <- gsea_wide %>%
  arrange(group)
gsea_sorted=as.data.frame(gsea_sorted)
gsea_sorted=subset(gsea_sorted,gsea_sorted$source=="cellline")
rownames(gsea_sorted)=gsea_sorted$ID
hm=as.data.frame(t(gsea_sorted[,c(4:53)]))
anno=as.data.frame(gsea_sorted[,c(2:3)])
# Set the color scale from white to red
color_palette <- colorRampPalette(c("navy", "white", "firebrick3"))(100)
# Define the breaks for the color scale
breaks <- seq(-2, 2, length.out = 100)
pdf("PCTA_hallmark.pdf", width = 36, height = 18)
pheatmap(hm, 
         cluster_cols = FALSE,
         annotation_col = anno,
         scale = "row",
         show_colnames = FALSE,
         gaps_col = c(470, 486, 614, 671, 722, 885, 981, 1159, 1212, 1225, 1243, 1253),
         breaks = breaks,
         color = color_palette)
dev.off()
#gaps_col = c(470, 486, 614, 671, 722, 885, 981, 1159, 1212, 1225, 1243, 1253, 1301, 1339, 1552, 1601, 1700, 1707, 1711, 1734),



### run TF activity analysis
# Retrieve prior knowledge network.
ctpc=readRDS("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/CTPC_V3/CTPC_TPM_batchCorrected.rds")
ctpc_meta=readRDS("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/CTPC_V3/CTPC_meta_cleaned.rds")
network <- decoupleR::get_dorothea(organism = "human",
                                   levels = c("A"))
# Run weighted means algorithm.
activities <- decoupleR::run_wmean(mat = as.matrix(ctpc),
                                   network = network,
                                   .source = "source",
                                   .targe = "target",
                                   .mor = "mor",
                                   times = 1000,
                                   minsize = 20)




wide_tf <- activities %>%
  filter(statistic == 'norm_wmean') %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source')


sorted_ctpc_meta <- ctpc_meta %>%
  arrange(group)

wide_tf=wide_tf[,rownames(sorted_ctpc_meta)]
breaks <- seq(-2, 2, length.out = 100)
color_palette <- colorRampPalette(c("navy", "white", "firebrick3"))(100)
pdf("PCTA_TFactivity.pdf", width = 36, height = 18)
pheatmap(wide_tf,cluster_cols = F,show_colnames = F,
         scale = "row",annotation_col = sorted_ctpc_meta[,2,drop=F],
         gaps_col = c(470, 486, 614, 671, 722, 885, 981, 1159, 1212, 1225, 1243, 1253),
         breaks = breaks,
         color = color_palette)
dev.off()


### generate data for CTPC-v2 app
ctpc=readRDS("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/CTPC_V3/CTPC_TPM_batchCorrected.rds")
ctpc_meta=readRDS("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/CTPC_V3/CTPC_meta_cleaned.rds")
### now only keep protein-coding genes
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
protein_coding_genes <- getBM(attributes = c("hgnc_symbol"), 
                              filters = "biotype", 
                              values = "protein_coding", 
                              mart = mart)
ctpc <- ctpc[rownames(ctpc) %in% protein_coding_genes$hgnc_symbol, ]
ctpc_app=as.data.frame(t(ctpc))
ctpc_app$ID=rownames(ctpc_app)
ctpc_app=melt(ctpc_app,id="ID")
colnames(ctpc_app)=c("ID","Gene","TPM")
saveRDS(ctpc_app,"CTPC_log2corrected_forApp.rds")



### now generate individual gene expression data for paper
ctpc=readRDS("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/CTPC_V3/CTPC_TPM_batchCorrected.rds")
ctpc_meta=readRDS("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/CTPC_V3/CTPC_meta_cleaned.rds")
#### only keep the protein coding genes, same as previous ProAtlas dataset
gene=readRDS("CTPC_ProAtlas_commen_genes.rds")
ctpc=ctpc[gene,]
### now only keep protein-coding genes
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
protein_coding_genes <- getBM(attributes = c("hgnc_symbol"), 
                              filters = "biotype", 
                              values = "protein_coding", 
                              mart = mart)
ctpc <- ctpc[rownames(ctpc) %in% protein_coding_genes$hgnc_symbol, ]

ctpc=as.data.frame(t(ctpc))
ctpc$ID=rownames(ctpc)
ctpc=merge(ctpc,ctpc_meta[,-4],by="ID")
ctpc=melt(ctpc,id=c("ID","group","Treatment","Batch"))
genes=c("KRT7","AR","INSM1","NKX2-1","CHGA","FOXA2","SOX2")
for (i in genes) {
  tem=subset(ctpc,ctpc$variable==i)
  tem$group=factor(tem$group,levels = rev(levels(tem$group)))
  ggplot(tem, aes(x = group, y = value, color = group, fill = group)) +
    geom_boxplot(outliers = FALSE,width=0.5) +
    geom_jitter(alpha =0.2,size=0.3,width = 0.3) +
    theme_classic() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_color_hue(l = 30)+coord_flip()
  ggsave(paste0(i,"_expression_PCTA.pdf"),width = 4,height = 6)
}


### now calculate outliers to identify treatment to gene expression relationship
# Define a function to identify outliers and calculate required statistics
find_outliers <- function(data) {
  data %>%
    group_by(group, variable) %>%
    mutate(
      Q1 = quantile(value, 0.25),
      Q3 = quantile(value, 0.75),
      IQR = Q3 - Q1,
      LowerBound = Q1 - 1.5 * IQR,
      UpperBound = Q3 + 1.5 * IQR,
      Outlier = ifelse(value < LowerBound | value > UpperBound, TRUE, FALSE),
      MedianValue = median(value),
      DifferenceToMedian = value - MedianValue
    ) %>%
    filter(Outlier == TRUE) %>%
    ungroup()
}


# Apply the function to your dataset
outliers <- find_outliers(ctpc)
# Calculate p-values for outliers
outliers <- outliers %>%
  group_by(group, variable) %>%
  mutate(
    p_value = 2 * pnorm(-abs((value - MedianValue) / sd(value)))
  ) %>%
  ungroup()
outliers_clean=subset(outliers,abs(outliers$DifferenceToMedian)>1.5)

colnames(outliers_clean)[5:6]=c("Gene","log2TPM")

saveRDS(outliers_clean,"outliers_clean.rds")
write.csv(outliers_clean,"outliers_clean.csv")

# visulize the outlier data
outliers_clean=readRDS("outliers_clean.rds")

table(outliers_clean$group)




ggplot(outliers_clean, aes(x = log2TPM, y = DifferenceToMedian, color = group)) +
  geom_point(size = 0.3) +
  labs(title = "log2TPM vs Difference to Median by Group", x = "log2TPM", y = "Difference to Median") +
  facet_wrap(~ group, ncol = 3) +
  coord_polar() +
  xlim(-1, 15) +  # Set x-axis range
  theme_minimal() +
  theme(
    # X-axis label black and positioned outside
    axis.title.x = element_text(size = 10, color = "black", vjust = 2), 
    
    # Y-axis label black
    axis.title.y = element_text(size = 10, color = "black"),
    
    # X and Y axis text size and color
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    
    # Gradient effect on y-axis grid lines
    panel.grid.major.y = element_line(color = scales::seq_gradient_pal("grey", "black")(seq(0, 1, length.out = 5))),
    panel.grid.major.x = element_line(color = "black"),  # X-axis grid lines in black
    panel.grid.minor = element_line(color = NA),  # No minor grid lines
    
    # Facet strip text settings
    strip.text = element_text(size = 8, angle = 0, hjust = 0.5),
    
    # Background and axis ticks settings
    panel.background = element_blank(),
    axis.ticks = element_line(color = "black")
  ) +
  scale_color_hue(l = 30) +
  # Adding an arrow for the x-axis
  annotate(
    "segment",
    x = 0, xend = 15,
    y = 0, yend = 0,
    color = "black",
    linewidth = 0.5,
    arrow = arrow(length = unit(0.2, "cm"), type = "closed")
  )
ggsave("outlier_distribution.pdf",width = 9,height = 11)

### generate sample IDs for mutation study
star=read.delim("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/CTPC_V3/rnaseq/multiqc_data/star_alignment_plot.txt",header = T)

star=merge(star,ctpc_meta,by.x="Sample",by.y="ID",all.x=F,all=F)

top_star <- star %>%
  group_by(group) %>%
  arrange(desc(Uniquely.mapped)) %>%  # Sort within each group by Uniquely.mapped in descending order
  slice_head(n = 20) %>%  # Keep the top 20 rows per group
  ungroup()  # Ungroup the data
write.csv(top_star,"CTPC_mutation_samples.csv")



#### VCF file PCA analysis
### genomics study
vcf_files <- list.files(path = ".", pattern = "^GSM.*\\.haplotypecaller\\.filtered\\.vcf\\.gz$", 
                        recursive = TRUE, full.names = TRUE)
# Convert each VCF file to GDS format
gds_files <- c()
for (i in 1:length(vcf_files)) {
  vcf.fn <- vcf_files[i]  # Use the file path directly
  
  # Extract the GSM identifier from the filename
  gsm_id <- sub("\\.haplotypecaller\\.filtered\\.vcf\\.gz$", "", basename(vcf.fn))
  
  # Create the GDS file name using the GSM identifier
  gds.fn <- paste0(gsm_id, ".gds")
  
  # Convert to GDS
  if (!file.exists(gds.fn)) {
    snpgdsVCF2GDS(vcf.fn, gds.fn, method = "biallelic.only")
  }
  
  # Store the GDS file name
  gds_files <- c(gds_files, gds.fn)
}

# Ensure gds_files is unique
gds_files <- unique(gds_files)

# Merge all GDS files into a single GDS file
merge_gds.fn <- "merged_data.gds"
snpgdsCombineGeno(gds_files, merge_gds.fn)

# Perform PCA on the merged GDS file
genofile <- snpgdsOpen(merge_gds.fn)
pca <- snpgdsPCA(genofile, autosome.only = F,bayesian = T)

# Extract and plot the PCA results
pca_df <- data.frame(PC1 = pca$eigenvect[,1], PC2 = pca$eigenvect[,2])
rownames(pca_df)=pca$sample.id
pca_df$ID=rownames(pca_df)
ctpc_meta=readRDS("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/CTPC_V3/CTPC_meta_cleaned.rds")
pca_df=merge(pca_df,ctpc_meta,by="ID")
snpgdsClose(genofile)
ggplot(pca_df, aes(x = PC1, y = PC2,color=group)) + geom_point(size=1) + theme_classic()+scale_color_hue(l = 50)
ggsave("LNCaP_celllines_genomic_PCA.pdf",width = 18,height = 16,units = "cm")


### now perform transcriptome PCA for these samples
### generate data for CTPC-v2 app
ctpc=readRDS("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/CTPC_V3/CTPC_TPM_batchCorrected.rds")
ctpc_meta=readRDS("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/CTPC_V3/CTPC_meta_cleaned.rds")
# Define the groups to keep
groups_to_keep <- c("LNCaP", "LNCaP-abl", "LNCaP-95", "C4-2", "C4-2B")
# Subset the table
ctpc_meta <- ctpc_meta[ctpc_meta$group %in% groups_to_keep, ]
ctpc=ctpc[,rownames(ctpc_meta)]
gene_variance = apply(ctpc, 1, var)
top_genes = names(sort(gene_variance, decreasing = TRUE))[1:5000]  # Select top 5000 most variable genes
TPM_top_var_genes = ctpc[top_genes, ]
# Initial PCA analysis with the most variable genes
matrix = as.matrix(t(TPM_top_var_genes))
pca_result = prcomp(matrix, center = TRUE, scale. = TRUE)
# Extract PCA results
pca_df = data.frame(Sample = rownames(pca_result$x), 
                    PC1 = pca_result$x[,1], 
                    PC2 = pca_result$x[,2])
# Merge with metadata
pca_df = merge(pca_df, ctpc_meta, by.x = "Sample", by.y = "ID")
# Plot the PCA results
ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 0.5) +
  labs(title = "PCA of Log2-Transformed RNAseq Data (Top Variable Genes)", x = "PC1", y = "PC2") +
  theme_classic() +
  scale_color_hue(l = 50)

### also try to use the non-batch-corrected data
ctpc=readRDS("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/CTPC_V3/TPM_clean.rds")
ctpc_meta=readRDS("/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/CTPC_V3/CTPC_meta_cleaned.rds")
# Define the groups to keep
groups_to_keep <- c("LNCaP", "LNCaP-abl", "LNCaP-95", "C4-2", "C4-2B")
# Subset the table
ctpc_meta <- ctpc_meta[ctpc_meta$group %in% groups_to_keep, ]
ctpc=ctpc[,rownames(ctpc_meta)]
gene_variance = apply(ctpc, 1, var)
top_genes = names(sort(gene_variance, decreasing = TRUE))[1:5000]  # Select top 5000 most variable genes
TPM_top_var_genes = ctpc[top_genes, ]
# Initial PCA analysis with the most variable genes
matrix = as.matrix(t(TPM_top_var_genes))
pca_result = prcomp(matrix, center = TRUE, scale. = TRUE)
# Extract PCA results
pca_df = data.frame(Sample = rownames(pca_result$x), 
                    PC1 = pca_result$x[,1], 
                    PC2 = pca_result$x[,2])
# Merge with metadata
pca_df = merge(pca_df, ctpc_meta, by.x = "Sample", by.y = "ID")
# Plot the PCA results
ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 0.5) +
  labs(title = "PCA of Log2-Transformed RNAseq Data (Top Variable Genes)", x = "PC1", y = "PC2") +
  theme_classic() +
  scale_color_hue(l = 50)
ggsave("LNCaP_celllines_uncorrectedTranscriptome_PCA.pdf",width = 18,height = 16,units = "cm")
