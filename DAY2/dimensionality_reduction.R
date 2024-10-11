"*******************************************************************************
DAY 2 Lecture 1/4
-----------------
After having completed this chapter you will be able to:

  - Describe and perform standard procedures for normalization 
    and scaling with the package Seurat
  - Select the most variable genes from a Seurat object for downstream analyses

*******************************************************************************"
seu <- readRDS("seu_day1-3.rds")

seu <- Seurat::RunPCA(seu)
Seurat::DimPlot(seu, reduction = 'pca')
# We can colour the PCA plot according to any factor that is present in @meta.data,
# or for any gene. 
# For example we can take the column percent.globin:
Seurat::FeaturePlot(seu, reduction = "pca", features = "percent.globin") 
# Generate a PCA plot where color is according to counts of a gene (i.e. gene expression). 
Seurat::FeaturePlot(seu, reduction = 'pca', features = "HBA1")

# We can generate heatmaps according to their principal component scores calculated in the rotation matrix:
Seurat::DimHeatmap(seu, dims = 1:12, cells = 500, balanced = TRUE)

# The elbowplot can help you in determining how many PCs to use for downstream analysis such as UMAP:
Seurat::ElbowPlot(seu, ndims = 40)
# 25 seems to capture most of the variance
"The elbow plot ranks principle components based on the percentage of variance 
explained by each one. 
Where we observe an “elbow” or flattening curve, 
the majority of true signal is captured by this number of PCs, 
eg around 25 PCs for the seu dataset."

# UMAP
# UMAP: The goal of these algorithms is to learn the underlying manifold of the data 
# in order to place similar cells together in low-dimensional space.

seu <- Seurat::RunUMAP(seu, dims= 1:25)
Seurat::DimPlot(seu, reduction = "umap")

# Exercise
# A: Color the dots in the UMAP according to a variable (e.g. percent.globin or HBA1). 
# Any idea where the erythrocytes probably are in the UMAP?
Seurat::FeaturePlot(seu, reduction = 'umap', features = c('HBA1','percent.globin'))

#B: change number of neigbors look at the effects
#The default number of neighbours is 30. 
seu <- Seurat::RunUMAP(seu, dims = 1:25, n.neighbors = 5)
Seurat::DimPlot(seu, reduction = "umap")
# C:
# The number of dims to extremes dims = 1:5 or dims = 1:50 
# how did it affect the output ? 
# In your opinion better few PCAs too much or too few ?
seu <- Seurat::RunUMAP(seu, dims = 1:5)
Seurat::DimPlot(seu, reduction = 'umap')
     
seu <- Seurat::RunUMAP(seu, dims = 1:50)
Seurat::DimPlot(seu, reduction = 'umap')               
# Taking dims = 1:100 does not work as in the step RunPCA by default only 50pcs are calculated, 

seu <- Seurat::RunUMAP(seu, dims = 1:25)
saveRDS(seu, "seu_after_dim_red.rds")
# END of Dimensionality Reduction
