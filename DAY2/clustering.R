library(Seurat)
seu <- readRDS("seu_day2-2.rds")



"The method implemented in Seurat first constructs a SNN graph based on the 
euclidean distance in PCA space, and refine the edge weights between any two 
cells based on the shared overlap in their local neighborhoods (Jaccard similarity). 
This step is performed using the FindNeighbors() function, 
and takes as input the previously defined dimensionality of the dataset."
seu <- Seurat::FindNeighbors(seu, dims = 1:25, reduction = "integrated.cca")
"To cluster the cells, Seurat next implements modularity optimization techniques 
such as the Louvain algorithm (default) 
to iteratively group cells together, 
with the goal of optimizing the standard modularity function."

"The FindClusters() function implements this procedure, 
and contains a resolution parameter that sets the ‘granularity’ of the downstream
clustering, with increased values leading to a greater number of clusters."
seu <- Seurat::FindClusters(seu, resolution = seq(0.1,0.8, by=0.1))
# Cluster id of each cell is added to the metadata object, 
# as a new column for each resolution tested:
head(seu@meta.data)
library(clustree)
clustree::clustree(seu@meta.data[,grep("RNA_snn_res", colnames(seu@meta.data))],
                   prefix = "RNA_snn_res.")
Seurat::DimPlot(seu, group.by = "RNA_snn_res.0.1")

Seurat::DimPlot(seu, group.by = "RNA_snn_res.0.3")
