library(SpatialDecon)
library(Rtsne)
library(ggplot2)

# My reverse deconvoluted object
load(url("https://github.com/preventionhub/MyTestingPlace/raw/main/Reverse%20Deconvolution.RData"))

# My matrix
rev.mat <- as.matrix(rdecon$coefs)

# ================================================================ #
#                             tSNE
# ================================================================ #
# Use the Rtsne() function to compute the t-SNE embedding of the data
set.seed(1234) # for reproducibility
tsne.obj <- Rtsne(rev.mat, perplexity = 35, dims = 2) # you can adjust the perplexity and dims as needed

# Convert the t-SNE object to a data frame and add the sample labels:
tsne.mat <- as.data.frame(tsne.obj$Y)
tsne.mat$Genes <- rownames(rev.mat)

# Plot the t-SNE embedding using ggplot2 :
ggplot(tsne.mat, aes(x = V1, y = V2)) + 
  geom_point(size = 0.5, color = "lightblue1") + 
  theme_bw()


# ================================================================ #
#                             UMAP
# In addition to regular stuff related to UMAP, here I would like to 
# try to "convert" reverse-deconvolution object to Seurat object
# ================================================================ #
library(umap)
library(Seurat)

# Create a Seurat object from reverseDecon object:
seurat.obj <- CreateSeuratObject(counts = rdecon$coefs)
str(seurat.obj)

# variable features
seurat.obj =SCTransform(seurat.obj)
seurat.obj = FindVariableFeatures(seurat.obj,nfeatures=20)

# Run normalization and scaling steps:
seurat.obj <- NormalizeData(object = seurat.obj)
seurat.obj <- ScaleData(object = seurat.obj)

# Perform Principal Component Analysis (PCA):
seurat.obj <- RunPCA(object = seurat.obj, ndims = 30, npcs = 11, approx=FALSE, verbose = FALSE)

# Generate UMAP projection:
seurat.obj <- RunUMAP(object = seurat.obj, dims = 1:11, n.neighbors = 11, seed.use = 12345)

# Store my UMAP data
umap_data <- as.data.frame(seurat.obj@reductions$umap@cell.embeddings)
umap_data$CellsType <- row.names(umap_data)

# Plot UMAP using ggplot2 :
ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = CellsType)) +
  geom_point(size = 1.5)+
  theme_bw()

