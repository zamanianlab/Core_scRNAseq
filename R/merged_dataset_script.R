# Pre-processing and analysis of tBM and utBM samples based on the CellRanger filtered dataset


# load libraries
library(tidyverse)
library(ggplot2)
library(cowplot)
library(DropletUtils)
library(scCB2)
library(dplyr)
library(expss)
library(multtest)
library(tidyquant)
library(SummarizedExperiment)
library(scDblFinder)
library(monocle)
library(methods)
library(Matrix)
library(scater)
library(patchwork)
library(Seurat)
library(SoupX)
library(ggplot2) 
library(ZamanianLabThemes)

# set directories for tBM and utBM data sets
data.dir.tBM <- ("~/Box/ZamanianLab/SeqLibraries/Mapping/singlecell/210518_BH7FNCDRXY/outs_tBM_only/filtered_feature_bc_matrix/")
data.dir.utBM <- ("~/Box/ZamanianLab/SeqLibraries/Mapping/singlecell/210518_BH7FNCDRXY/outs_utBM_only/filtered_feature_bc_matrix/")


# read in the data
tBM <- Read10X(data.dir.tBM)
utBM <- Read10X(data.dir.utBM)


# create Seurat objects
tBM <- CreateSeuratObject(counts = tBM, project = "20210413_bma_tBM")
utBM <- CreateSeuratObject(utBM, project = "20210413_bma_utBM")


# create merged Seurat object to combine tBM and utBM prior to filtering the dataset
merged <- merge(tBM, utBM, add.cell.ids = c("tBM", "utBM"), project = "20210413_bma")
view(merged@meta.data)


# What percentage of all counts for each cell belong to mitochondrial genes WBGene00225418, WBGene00225415, WBGene00220387?

merged[["percent.mt"]] <- PercentageFeatureSet(merged, features = c("WBGene00225418" , "WBGene00225415", "WBGene00220387"))
view(merged@meta.data)

metadata <- merged@meta.data

# view the percent.mt distribution for each treatment group via density plot 
(mtrna_dis_plot <- metadata %>% 
    ggplot(aes(color=orig.ident, x=percent.mt, fill= orig.ident)) + 
    geom_density(alpha = 0.2) + 
    theme_classic() +
    scale_x_continuous(breaks = seq(0,50, 10), expand = c(0,0)) +
    ylab("Cell density") +
    NULL)


# look at gene representation for each cell and identiy cells that are largly represented by a single or few genes
merged$Percent.Largest.Gene <- apply(
  merged@assays$RNA@counts,
  2,
  function(x)(100*max(x))/sum(x)
)

view(merged@meta.data)


# Visualize QC metrics as a violin plots
VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3) + scale_y_log10()


#plotting different features with correlation as the plot title 
FeatureScatter(merged, feature1 = "nCount_RNA", feature2 = "percent.mt") # no relationship

FeatureScatter(merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") # strong positive (0.91)



# Conservative subsetting of the data
  merged_subset <- subset(merged,
      nCount_RNA < 60000 & 
      nFeature_RNA < 6000 & 
      percent.mt < 20 & 
      Percent.Largest.Gene < 40)
  
view(merged_subset@meta.data)


# Does the data look any better? Plotting different features with correlation as the plot title again
FeatureScatter(merged_subset, feature1 = "nCount_RNA", feature2 = "percent.mt") # no change

FeatureScatter(merged_subset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") # strong positive (0.94)


# Normalize the data
norm <- NormalizeData(merged_subset, normalization.method = "CLR")

as.tibble(
  norm@assays$RNA@data[,1:100]
) %>%
  pivot_longer(
    cols=everything(),
    names_to="cell",
    values_to="expression"
  ) %>%
  ggplot(aes(x=expression, group=cell)) +
  geom_density() +
  coord_cartesian(ylim=c(0,0.6), xlim=c(0,3))


# Identify genes in cells that vary in expression (high in some and low in others)
hvfeatures <- FindVariableFeatures(norm, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(hvfeatures), 10)
plot1 <- VariableFeaturePlot(hvfeatures)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Scaling the data (linear transformation) prior to PCA analysis. 
all.genes <- rownames(norm)
scaled <- ScaleData(norm, features = all.genes)


# Run PCA analysis (linear dimensional reduction) on the scaled data
PCA <- RunPCA(scaled, features = VariableFeatures(object = hvfeatures))

PCA_plots <- plot_grid(ncol = 3, DimPlot(PCA, reduction = "pca", group.by = "orig.ident", dims = 1:2), DimPlot(PCA, reduction = "pca", group.by = "orig.ident", dims = 3:4), DimPlot(PCA, reduction = "pca", group.by = "orig.ident", dims = 5:6))


# Look at the variance explained by each PC to identify if what PCs should be included in downstream annalysis
ElbowPlot(PCA, reduction = "pca", ndims = 50)



# clustering the cells
knn <- FindNeighbors(PCA, dims = 1:10)
clusters <- FindClusters(knn, resolution = 0.5)



# UMAP
UMAP <- RunUMAP(clusters, dims = 1:10)
UMAP_plot <-DimPlot(UMAP, reduction = "umap", label = TRUE, split.by = "orig.ident")

# tSNE
tSNE <- RunTSNE(
  clusters,
  dims=1:15,
  seed.use = 10403, 
  perplexity=30
) 

DimPlot(tSNE,reduction = "tsne", pt.size = 1, split.by = "orig.ident") + ggtitle("tSNE with Perplexity 30")


# identify cluster biomarkers 
markers <- FindAllMarkers(tSNE, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)



