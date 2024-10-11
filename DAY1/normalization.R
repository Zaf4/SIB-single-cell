"*******************************************************************************
DAY 1 Lecture 4/4
-----------------
After having completed this chapter you will be able to:

  - Describe and perform standard procedures for normalization 
    and scaling with the package Seurat
  - Select the most variable genes from a Seurat object for downstream analyses

*******************************************************************************"
seu <- readRDS("DAY1/seu_after_QC.rds")

# NORMALIZATION
## After removing unwanted cells from the dataset, the next step is to normalize the data
Seurat::GetAssayData(seu, layer = "counts") #slot is depr. use `layer`
seu <- Seurat::NormalizeData(seu,
                             scale.factor = 10000,
                             normalization.method = "LogNormalize")

# VARIABLE FEATURES
## We next calculate a subset of features that exhibit high cell-to-cell variation in the dataset.
## Focusing on these genes in downstream analysis helps to highlight biological signal in single-cell datasets.
seu <- Seurat::FindVariableFeatures(seu,
                                    selection.method = "vst", # default is vst
                                    nfeatures = 2000) # default is 2000 

# Identify the 10 most highly variable genes
top10 <- head(Seurat::VariableFeatures(seu), 10)
top10

# variable feature plot to visualize gene expression levels and variances
vf_plot <- Seurat::VariableFeaturePlot(seu)
Seurat::LabelPoints(plot = vf_plot, # labels points the data with given point names
                    points = top10,
                    repel = TRUE)
vf_plot

# SCALING
## Next, we apply scaling, a linear transformation that is a standard
## pre-processing step prior to dimensional reduction techniques like PCA.
"*******************************************************************************
The ScaleData() function
  1. shifts the expression of each gene, so that the mean expression across cells is 0
  2. scales the expression of each gene, so that the variance across cells is 1

This step gives equal weight in downstream analyses, 
so that highly-expressed genes do not dominate. 
The results of this are stored in seu$RNA@scale.data
*******************************************************************************"
seu <- Seurat::ScaleData(seu)

seu <- Seurat::SCTransform(seu) 
# Covers all prior functions (NormalizeData, VariableFeatures and ScaleData)
names(seu@assays)
Seurat::DefaultAssay(seu) <- "RNA"

# save day1 output
saveRDS(seu, "seu_day1-3.rds")

# END of DAY1
