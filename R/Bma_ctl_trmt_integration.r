## Integration of tBM and utBM based on Seurat integration tutorial 

# load required libraries
library(ZamanianLabThemes)
library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(pals)
library(ggtext)



## Perform the integration starting with the utBM_subset and tBM_subset processed Seurat objects
utBM_subset <- readRDS("/Users/chenthorn/Library/CloudStorage/Box-Box/ZamanianLab/SeqLibraries/Mapping/singlecell/scRNAesq_Bma_multi/RDS_objects/scmulti_utBM.rds")
tBM_subset <- readRDS("/Users/chenthorn/Library/CloudStorage/Box-Box/ZamanianLab/SeqLibraries/Mapping/singlecell/scRNAesq_Bma_multi/RDS_objects/scmulti_tBM.rds")

# Identify anchors in both datasets with the FindIntegrationAnchors function 
# Completes canonical correlation analysis (CCA) by using PCA to find the greaest sources of shared variation across the groups and then identifies the anchors across the datasets
# The anchors or mutual nearest neighbors (MNNs) are cellular 'best buddies' --> looking for a cells closest neighbor in the other condition based on gene expression and if the reciprocal cell also identifies the original cell as a neighbor, the cells are marked as anchors
# incorrect anchors are removed by comparing the similarity between anchor pairs by the overlap in the local neighborhoods - do the adjacent cell shave "best buddies" that are adjacent to each other?

BM_anchors <- FindIntegrationAnchors(object.list = list(utBM_subset, tBM_subset), dims = 1:50) #69995 anchors identified, filtered down to 14207

# Integrate the two datasets based on the anchors identified (If cell types are present in one dataset but not the other, then the cells will appear as a separate sample-specific cluster)
BM_combined <- IntegrateData(anchorset = BM_anchors, dims = 1:50, features.to.integrate = )


# Now we can run dimensionality reduction to visualize
BM_combined <- ScaleData(BM_combined, verbose = TRUE)
BM_combined <- RunPCA(BM_combined, verbose = TRUE)
BM_combined <- RunUMAP(BM_combined, reduction = "pca", dims = 1:50)
BM_combined <- FindNeighbors(BM_combined, reduction = "pca", dims = 1:50)
BM_combined <- FindClusters(BM_combined, resolution = 0.5)

p1 <- DimPlot(BM_combined, reduction= "umap", group.by = "orig.ident") + ggtitle("Integrated Bma Datasets")
p2 <- DimPlot(BM_combined, reduction = "umap", label = TRUE, cols = "polychrome")
plot_grid(p1, p2)



# how many cells are in each cluster and what proportion of those cells is represented by tBM/utBM
table(BM_combined@active.ident, group_by = BM_combined@meta.data$orig.ident)


saveRDS(BM_combined, "/Users/chenthorn/Library/CloudStorage/Box-Box/ZamanianLab/LabMembers/Clair/project_singlecell/other/BM_combined.RDS")
