#Using filtered, processed seurat objects of B. malayi mf tBM (1 uM IVM treated) and utBM (untreated) as input for cluster analysis. 

#Load libraries
library(tidyverse)
library(ggplot2)
library(cowplot)
library(DropletUtils)
library(scCB2)
library(dplyr)
library(methods)
library(Matrix)
library(Seurat)
library(SoupX)
library(ZamanianLabThemes)

# see package vignettes
vignette("scCB2")
browseVignettes("scDblFinder")

# read in subsetted tBM and utBM Seurat objects
tBM_subset <- readRDS("/Users/clairhenthorn/Library/CloudStorage/Box-Box/ZamanianLab/SeqLibraries/Mapping/singlecell/scRNAesq_Bma_multi/RDS_objects/scmulti_tBM_subset.rds")
utBM_subset <- readRDS("/Users/clairhenthorn/Library/CloudStorage/Box-Box/ZamanianLab/SeqLibraries/Mapping/singlecell/scRNAesq_Bma_multi/RDS_objects/scmulti_utBM_subset.rds")


# Normalize the data (requires subsetted Seurat object from pre-processing)
tBM_subset <-NormalizeData(tBM_subset, normalization.method = "CLR")
utBM_subset <-NormalizeData(utBM_subset, normalization.method = "CLR")


# Identify genes in cells that vary in expression (high in some and low in others)
tBM_hvfeatures <- FindVariableFeatures(tBM_subset, selection.method = "vst", nfeatures = 2000)
utBM_hvfeatures <- FindVariableFeatures(utBM_subset, selection.method = "vst", nfeatures = 2000)


# Scaling the data (linear transformation) prior to PCA analysis.

tBM_all.genes <- rownames(tBM_subset)
tBM_subset<- ScaleData(tBM_subset, features = tBM_all.genes)

utBM_all.genes <- rownames(utBM_subset)
utBM_subset <- ScaleData(utBM_subset, features = utBM_all.genes)


# Run PCA analysis (linear dimensional reduction) on the scaled data (defult dims = 50, change to 100)
tBM_subset <- RunPCA(tBM_subset, npcs = 100, features = VariableFeatures(object = tBM_hvfeatures))
utBM_subset <- RunPCA(utBM_subset, npcs = 100, features = VariableFeatures(object = utBM_hvfeatures))




# Look at the variance explained by each PC to identify what PCs should be included in downstream annalysis
ElbowPlot(tBM_subset, reduction = "pca", ndims = 100)
ElbowPlot(utBM_subset, reduction = "pca", ndims = 100)



# How many PCs should be included in order to capture the robustness of the dataset? Significant dropoff in significance at PC44

tBM_subset <- ScoreJackStraw(tBM_subset, dims = 1:50)
JackStrawPlot(tBM_subset, dims = 1:50)

utBM_subset <- ScoreJackStraw(utBM_subset, dims = 1:50)
JackStrawPlot(utBM_subset, dims = 1:50)


# clustering the cells
tBM_subset <- FindNeighbors(tBM_subset, dims = 1:50)
tBM_subset <- FindClusters(tBM_subset, resolution = 0.5)

utBM_subset<- FindNeighbors(utBM_subset, dims = 1:50)
utBM_subset<- FindClusters(utBM_subset, resolution = 0.5)



# UMAP

tBM_subset <- RunUMAP(tBM_subset, seed.use=123, dims = 1:50)
utBM_subset <- RunUMAP(utBM_subset,seed.use=123, dims = 1:50)


DimPlot(tBM_subset, reduction = "umap", label = TRUE, cols = "polychrome") + ggtitle("tBM_UMAP, res = 0.5")
DimPlot(utBM_subset, reduction = "umap", label = TRUE, cols = "polychrome") + ggtitle("utBM_UMAP, res = 0.5")


# tSNE
tBM_subset <- RunTSNE(
  tBM_subset,
  dims=1:44,
  seed.use = 10403,
  perplexity=50
)

DimPlot(tBM_subset,reduction = "tsne", pt.size = 0.5, split.by = "orig.ident", label = TRUE) + ggtitle("tBM_tSNE with Perplexity 50")

utBM_subset <- RunTSNE(
  utBM_subset,
  dims=1:44,
  seed.use = 10403, 
  perplexity=50
) 

DimPlot(utBM_subset,reduction = "tsne", pt.size = 0.5, split.by = "orig.ident", label = TRUE) + ggtitle("utBM_tSNE with Perplexity 50")





#save RDS objects for import used in downstream analysis

saveRDS(tBM_subset, "~/Box/ZamanianLab/SeqLibraries/Mapping/singlecell/scRNAesq_Bma_multi/RDS_objects/scmulti_tBM_subset.rds")

saveRDS(utBM_subset, "~/Box/ZamanianLab/SeqLibraries/Mapping/singlecell/scRNAesq_Bma_multi/RDS_objects/scmulti_utBM_subset.rds")
