library(Seurat)
"*******************************************************************************
DAY 2 Lecture 2/4
-----------------

*******************************************************************************"
seu<-readRDS("seu_after_dim_red.rds")

# Let’s have a look at the UMAP again. 
Seurat::DimPlot(seu, reduction = "umap")
# Although cells of different samples are shared amongst ‘clusters’, 
# you can still see separation within the clusters

"To perform the integration, we split our object by sample, 
resulting into a set of layers within the RNA assay. 
The layers are integrated and stored in the reduction slot - 
in our case we call it integrated.cca. Then, we re-join the layers"
seu[["RNA"]] <- split(seu[["RNA"]], f = seu$orig.ident)
# integrate layers
seu <- Seurat::IntegrateLayers(object = seu,
                               method = Seurat::CCAIntegration,
                               orig.reduction = "pca",
                               new.reduction = "integrated.cca",
                               )
# re-join layers after integration
seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
# We can then use this new integrated matrix for clustering and visualization.

# Now, we can re-run and visualize the results with UMAP.
seu <- RunUMAP(seu, dims = 1:30, reduction = "integrated.cca")
Seurat::DimPlot(seu, reduction = "umap")

saveRDS(seu, "seu_day2-2.rds")
# END of integration

