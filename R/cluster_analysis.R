#Using filtered, processed seurat objects of B. malayi mf tBM (1 uM IVM treated) and utBM (untreated) as input for cluster analysis. 

#Install packages
BiocManager::install("scCB2")
install.packages("Seurat")
BiocManager::install("DropletUtils")
BiocManager::install("scDblFinder")
install.packages("SoupX")

#Load libraries
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

# see package vignettes
vignette("scCB2")
browseVignettes("scDblFinder")



# Normalize the data (requires subsetted Seurat object from pre-processing)
norm <- NormalizeData(merged, normalization.method = "CLR")
tBM_norm <-NormalizeData(tBM_subset, normalization.method = "CLR")
utBM_norm <-NormalizeData(utBM_subset, normalization.method = "CLR")


# Identify genes in cells that vary in expression (high in some and low in others)
hvfeatures <- FindVariableFeatures(norm, selection.method = "vst", nfeatures = 2000)
tBM_hvfeatures <- FindVariableFeatures(norm, selection.method = "vst", nfeatures = 2000)
utBM_hvfeatures <- FindVariableFeatures(norm, selection.method = "vst", nfeatures = 2000)


# Scaling the data (linear transformation) prior to PCA analysis. 
all.genes <- rownames(norm)
scaled <- ScaleData(norm, features = all.genes)

tBM_all.genes <- rownames(tBM_norm)
tBM_scaled <- ScaleData(tBM_norm, features = tBM_all.genes)

utBM_all.genes <- rownames(utBM_norm)
utBM_scaled <- ScaleData(utBM_norm, features = utBM_all.genes)


# Run PCA analysis (linear dimensional reduction) on the scaled data (defult dims = 50, change to 150)
PCA <- RunPCA(scaled, npcs = 100, features = VariableFeatures(object = hvfeatures))
tBM_PCA <- RunPCA(tBM_scaled, npcs = 100, features = VariableFeatures(object = tBM_hvfeatures))
utBM_PCA <- RunPCA(utBM_scaled, npcs = 100, features = VariableFeatures(object = utBM_hvfeatures))

PCA_plots <- plot_grid(ncol = 3, DimPlot(PCA, reduction = "pca", group.by = "orig.ident", dims = 1:2), DimPlot(PCA, reduction = "pca", group.by = "orig.ident", dims = 3:4), DimPlot(PCA, reduction = "pca", group.by = "orig.ident", dims = 5:6))



# Look at the variance explained by each PC to identify if what PCs should be included in downstream annalysis
ElbowPlot(PCA, reduction = "pca", ndims = 100)



# How many PCs should be included in order to capture the robustness of the dataset? Significant dropoff in significance at PC44
jackstraw <- JackStraw(PCA, num.replicate = 100, dims = 50)
jackstraw <- ScoreJackStraw(jackstraw, dims = 1:50)
JackStrawPlot(jackstraw, dims = 1:50)



# clustering the cells
knn <- FindNeighbors(PCA, dims = 1:50)
clusters <- FindClusters(knn, resolution = 0.5)

tBM_knn <- FindNeighbors(tBM_PCA, dims = 1:50)
tBM_clusters <- FindClusters(tBM_knn, resolution = 0.5)

utBM_knn <- FindNeighbors(utBM_PCA, dims = 1:50)
utBM_clusters <- FindClusters(utBM_knn, resolution = 0.5)



# UMAP
UMAP <- RunUMAP(clusters, dims = 1:50) # merged
tBM_UMAP <- RunUMAP(tBM_clusters, dims = 1:50)
utBM_UMAP <- RunUMAP(utBM_clusters, dims = 1:50)

DimPlot(UMAP, reduction = "umap", label = TRUE, split.by = "orig.ident", cols = "polychrome") + ggtitle("UMAP, res = 0.5")

DimPlot(tBM_UMAP, reduction = "umap", label = TRUE, cols = "polychrome") + ggtitle("tBM_UMAP, res = 0.5")
DimPlot(utBM_UMAP, reduction = "umap", label = TRUE, cols = "polychrome") + ggtitle("utBM_UMAP, res = 0.5")


# tSNE
tSNE <- RunTSNE(
  clusters,
  dims=1:44,
  seed.use = 10403, 
  perplexity=50
) 

DimPlot(tSNE,reduction = "tsne", pt.size = 0.5, split.by = "orig.ident", label = TRUE) + ggtitle("tSNE with Perplexity 50")



# identify cluster biomarkers based on UMAP. only.pos = TRUE indicating it will only report positive markers
markers <- FindAllMarkers(UMAP, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
tBM_markers <- FindAllMarkers(tBM_UMAP, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
utBM_markers <- FindAllMarkers(utBM_UMAP, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)



# import dataframe of B. malayi and C. elegans orthologs
orthologs <- read.csv("~/Desktop/ortholog.csv")
colnames(orthologs) <- c("Bma_gene_id", "Bma_gene_name", "cel_ortho")



# add the C. elegans ortholog to the markers object for each cluster
utBM_markers <- rownames_to_column(utBM_markers, var = "Bma_gene_id")

(utBM_markers_ortho <- utBM_markers %>% 
    left_join(orthologs, utBM_markers, by = "Bma_gene_id"))


(tBM_markers_ortho <- tBM_markers %>% 
    left_join(orthologs, tBM_markers, by = "Bma_gene_id"))




# load cengen enriched genes by cell type
exc_cell <- read.csv("~/Desktop/cengen_markers/CellMarkers_Excretory_cell.csv", header = TRUE, sep = ",")
intestine <- read.csv("~/Desktop/cengen_markers/CellMarkers_Intestine.csv", header = TRUE, sep = ",")
seam <- read.csv("~/Desktop/cengen_markers/CellMarkers_Seam_cell.csv", header = TRUE, sep = ",")
coel <- read.csv("~/Desktop/cengen_markers/CellMarkers_Coelomocyte.csv", header = TRUE, sep = ",")
bwm <- read.csv("~/Desktop/cengen_markers/CellMarkers_Body_wall_muscle.csv", header = TRUE, sep = ",")
epidermis <- read.csv("~/Desktop/cengen_markers/CellMarkers_Epidermis.csv", header = TRUE, sep = ",")
pharyngeal_mus <- read.csv("~/Desktop/cengen_markers/CellMarkers_pharyngeal_muscle.csv", header = TRUE, sep = ",")
head_mc <- read.csv("~/Desktop/cengen_markers/CellMarkers_hmc.csv", header = TRUE, sep = ",")
microarray <- read.csv("~/Desktop/cengen_markers/microarray.csv", header = TRUE, sep = ",")
germ <- read.csv("~/Desktop/cengen_markers/CellMarkers_Germline.csv", header = TRUE, sep = ",")



# add the C. elegans ortholog to the markers object for each cluster
(utBM_cellident <- utBM_markers_ortho %>% 
    left_join(exc_cell, utBM_markers_ortho, by = "cel_ortho") %>% 
    left_join(intestine, utBM_markers_ortho, by = "cel_ortho") %>% 
    left_join(seam, utBM_markers_ortho, by = "cel_ortho") %>% 
    left_join(coel, utBM_markers_ortho, by = "cel_ortho") %>% 
    left_join(bwm, utBM_markers_ortho, by = "cel_ortho") %>% 
    left_join(pharyngeal_mus, utBM_markers_ortho, by = "cel_ortho") %>%
    left_join(head_mc,utBM_markers_ortho, by = "cel_ortho") %>% 
    left_join(epidermis, utBM_markers_ortho, by = "cel_ortho") %>% 
    left_join(microarray, utBM_markers_ortho, by = "cel_ortho") %>% 
    left_join(germ, utBM_markers_ortho, by = "cel_ortho"))

(tBM_cellident <- tBM_markers_ortho %>% 
    left_join(exc_cell, tBM_markers_ortho, by = "cel_ortho") %>% 
    left_join(intestine, tBM_markers_ortho, by = "cel_ortho") %>% 
    left_join(seam, tBM_markers_ortho, by = "cel_ortho") %>% 
    left_join(coel, tBM_markers_ortho, by = "cel_ortho") %>% 
    left_join(bwm, tBM_markers_ortho, by = "cel_ortho") %>% 
    left_join(pharyngeal_mus, tBM_markers_ortho, by = "cel_ortho") %>%
    left_join(head_mc,tBM_markers_ortho, by = "cel_ortho") %>% 
    left_join(epidermis, tBM_markers_ortho, by = "cel_ortho") %>% 
    left_join(microarray, utBM_markers_ortho, by = "cel_ortho") %>% 
    left_join(germ, utBM_markers_ortho, by = "cel_ortho"))



# Identificatin of coelomocytes using markers cup-4 ("WBGene00268467") and lgc-26 ("WBGene00227182")
FeaturePlot(utBM_UMAP, features = c("WBGene00268467","WBGene00227182")) # cluster 9
FeaturePlot(tBM_UMAP, features = c("WBGene00268467","WBGene00227182")) # cluster 7

# Identification of neurons
# DVA tail interneuron (marker nlp-12)
FeaturePlot(tBM_UMAP, features = "WBGene00225297") 
FeaturePlot(utBM_UMAP, features = "WBGene00225297") 

#cholinergic neurons (cha-1, unc-17, cho-1)
FeaturePlot(utBM_UMAP, features = "WBGene00227648") #cha-1
FeaturePlot(tBM_UMAP, features = "WBGene00227648") #cha-1

FeaturePlot(tBM_UMAP, features = "WBGene00223095") #cho-1
FeaturePlot(utBM_UMAP, features = "WBGene00223095") #cho-1

FeaturePlot(utBM_UMAP, features = "WBGene00231904") #acr-15
FeaturePlot(tBM_UMAP, features = "WBGene00231904") #acr-15

#GABAergic neurons (unc-25)
FeaturePlot(tBM_UMAP, features = "WBGene00226652") 
FeaturePlot(utBM_UMAP, features = "WBGene00226652") 

#Dopaminergic neurons (dat-1, cat-2)

#Pharyngeal neurons (flr-2, ser-7, eya-1, pha-4)
FeaturePlot(tBM_UMAP, features = "WBGene00229074") #eya-1
FeaturePlot(utBM_UMAP, features = "WBGene00229074")

#AWA (odr-10), ASG (gcy-15), ASE (ASER= gcy-3, gcy-5; ASEL= gcy-6, gcy-7), AFD (gcy-8), and ASK (snet-1 and zig-4) neurons. 


# ASI/ASJ (ins-6, daf-28), AWB/AWC (daf-11, not ins-6 or snet-1), BAG, URX, SDQ, and other ciliated 3 sensory neurons.

#ARX, AQR, PQR, SDQ, ALN, PLN O2-sensory neurons, as well as the AVM 13 and  BDU  neurons  (gcy-35)
FeaturePlot(tBM_UMAP, features = "WBGene00223194") #gcy-35
FeaturePlot(utBM_UMAP, features = "WBGene00223194")#??? undefined


#ADQ/ALN/PLN O2-sensory neurons.(mec-1, lad-2, gcy-35)






#Add new cell type identities to the clusters
utBM_new_idents <- c("0", "1", "2", "Muscle", "Cholinergic neurons", "5", "Neuron", "7", "Muscle", "Coelomocytes", "Canal-associated neurons", "11", "12", "13", "14", "15", "DVA tail interneuron", "17", "18", "19", "GABAergic neurons")

names(utBM_new_idents) <- levels(utBM_UMAP)
utBM_UMAP.ID <- RenameIdents(utBM_UMAP, utBM_new_idents)
DimPlot(utBM_UMAP.ID, reduction = "umap", pt.size = 0.5, cols = "polychrome")



tBM_new_idents <- c("0", "1", "2", "Muscle", "Cholinergic neurons", "5", "Ciliated sensory neurons", "Coelomocyotes", "8", "Muscle", "10", "Canal-associated neurons", "12", "13", "14", "15", "16", "DVA tail interneuron", "18", "19", "20", "GABAergic neurons", "Oxygen sensory, AVM, BDU neurons", "23", "24")

names(tBM_new_idents) <- levels(tBM_UMAP)
tBM_UMAP.ID <- RenameIdents(tBM_UMAP, tBM_new_idents)
DimPlot(tBM_UMAP.ID, reduction = "umap", pt.size = 0.5, cols = "polychrome")






# feature plot looking at all the GluCls in the dataset 
# avr-14
FeaturePlot(UMAP, features ="WBGene00221971", split.by = "orig.ident")

# glc-2
FeaturePlot(UMAP, features = "WBGene00223839", split.by = "orig.ident")

#glc-3
FeaturePlot(UMAP, features ="WBGene00228311", split.by = "orig.ident")

#glc-4
FeaturePlot(UMAP, features = "WBGene00222703", split.by = "orig.ident")

#ben-1
FeaturePlot(UMAP, features = "WBGene00224994", split.by = "orig.ident")




# Cel body wall muscle markers: myo-3, pat-10, egl-20, php-3 (respectively)
FeaturePlot(UMAP, features = c("WBGene00231447","WBGene00224604", "WBGene00226493", "WBGene00229993"))


FeaturePlot(tBM_UMAP, features = c("WBGene00223381", "WBGene00223147", "WBGene00221982", "WBGene00225764"), cols = "blue") #egl-21, egl-3, ida-1, sbt-1 for distinguishing neurons, glia, and excretory cells




## THIS WOULD BE COOL!!!
plot <- FeaturePlot(tBM_UMAP, features = "WBGene00223381")
HoverLocator(plot = tBM_UMAP, information = FetchData(tBM_UMAP, vars = c("")))



