#10X Genomics scRNAseq analysis of B. malayi mf +/- 1 uM Ivermectin treatment
#Reading and filtering raw feature barcode matrix results from 10X Genomics Cellranger output. Filtration completed with scCB2 R package followed by SoupX to remove cell free mRNA contamination. Script completes with filtred, processed seurat objects ready for cluster analysis.

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


# Set PSOCK cluster setup_strategt to "sequential". Defult is "parallel". Required to use CB2CellFinder. 
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")


# set directories for tBM and utBM data sets
data.dir.tBM <- ("~/Box/ZamanianLab/SeqLibraries/Mapping/singlecell/scRNAesq_Bma_multi/scmulti_tBM/outs/raw_feature_bc_matrix")
data.dir.utBM <- ("~/Box/ZamanianLab/SeqLibraries/Mapping/singlecell/scRNAesq_Bma_multi/scmulti_utBM/outs/raw_feature_bc_matrix")


# Read in raw bc matrix from 10X Genomics Cellranger output
# read in raw data as a dgCMatrix
raw_tBM <- Read10xRaw(data.dir.tBM)
raw_utBM <- Read10xRaw(data.dir.utBM)



### Begin filtration of the raw data to remove empty droplets using scCB2

# run CB2 to distinguish real cells from empty droplets. Default lower=100
CBout_tBM <- CB2FindCell(raw_tBM, FDR_threshold = 0.01, lower = 100, Ncores = 2)
str(assay(CBout_tBM)) 
str(metadata(CBout_tBM))

CBout_utBM <- CB2FindCell(raw_utBM, FDR_threshold = 0.01, lower = 100, Ncores = 2)
str(assay(CBout_utBM)) 
str(metadata(CBout_utBM))



# extract the real cell matrix
RealCell_tBM <- GetCellMat(CBout_tBM)
str(RealCell_tBM)

RealCell_utBM <- GetCellMat(CBout_utBM)
str(RealCell_utBM)



# Create a Seurat object from the scCB2 filtered cells. Remove genes detected in less than 3 cells
tBM <- Seurat::CreateSeuratObject(counts = RealCell_tBM,  min.cells = 3, project = "tBM")
view(tBM@meta.data)

utBM <- Seurat::CreateSeuratObject(counts = RealCell_utBM, min.cells = 3,  project = "utBM")
view(utBM@meta.data)




# Estimating cell free mRNA contamination using SoupX package. SoupX requires a raw and a filtered bc matrix to estimate the % contamination. Because the cellranger filtered bc matrix removes many cells, the raw bc matrix output is filtered using scCB2 and then a new filtered bc matrix is made from this output for input into soupx for each sample. Additionally, SoupX can better parse out contamination when cluster information is supplied. Each filtered treatment group will be taken all the way through UMAP to retrieve cluster information for SoupX input. Because each sample is taken from the same sc prep, the background should be relatively the same...

# Normalize the data (non-subsetted)
norm_tBM <- NormalizeData(tBM, normalization.method = "CLR")
norm_utBM <- NormalizeData(utBM, normalization.method = "CLR")


# Identify genes in cells that vary in expression (high in some and low in others)
hvfeatures_tBM <- FindVariableFeatures(norm_tBM, selection.method = "vst", nfeatures = 2000)
hvfeatures_utBM <- FindVariableFeatures(norm_utBM, selection.method = "vst", nfeatures = 2000)


# Scaling the data (linear transformation) prior to PCA analysis. 
all.genes_tBM <- rownames(norm_tBM)
scaled_tBM <- ScaleData(norm_tBM, features = all.genes_tBM)

all.genes_utBM <- rownames(norm_utBM)
scaled_utBM <- ScaleData(norm_utBM, features = all.genes_utBM)


# Run PCA analysis (linear dimensional reduction) on the scaled data
PCA_tBM <- RunPCA(scaled_tBM, features = VariableFeatures(object = hvfeatures_tBM))
PCA_utBM <- RunPCA(scaled_utBM, features = VariableFeatures(object = hvfeatures_utBM))


# Look at the variance explained by each PC to identify if what PCs should be included in downstream annalysis
ElbowPlot(PCA_tBM, reduction = "pca", ndims = 50)
ElbowPlot(PCA_utBM, reduction = "pca", ndims = 50)


# clustering the cells
knn_tBM <- FindNeighbors(PCA_tBM, dims = 1:44)
clusters_tBM <- FindClusters(knn_tBM, resolution = 0.5)

knn_utBM <- FindNeighbors(PCA_utBM, dims = 1:44)
clusters_utBM <- FindClusters(knn_utBM, resolution = 1.0)


# UMAP
UMAP_tBM <- RunUMAP(clusters_tBM, dims = 1:44)
DimPlot(UMAP_tBM, reduction = "umap", label = TRUE, cols = "polychrome") + ggtitle("UMAP, res = 0.5")

UMAP_utBM <- RunUMAP(clusters_utBM, dims = 1:44)
DimPlot(UMAP_utBM, reduction = "umap", label = TRUE, cols = "polychrome") + ggtitle("UMAP, res = 1.0")



# Load the raw 10X output and the scCB2 filtered bc matrices. The estimateSoup function to estimate the expression profile of the cell free mRNA (soup) is computed automatically when using the SoupChannel function.
souped_tBM = SoupChannel(raw_tBM, RealCell_tBM, calcSoupProfile = FALSE)
souped_tBM = estimateSoup(souped_tBM, soupRange=c(0,25))
souped_tBM

souped_utBM = SoupChannel(raw_utBM, RealCell_utBM, calcSoupProfile = FALSE)
souped_utBM = estimateSoup(souped_utBM, soupRange=c(0,25))
souped_utBM

#SoupX can better understand contamination is dimension reduction data is supplied.

Embeddings(UMAP_tBM)[,1:2]
UMAP_DR_tBM <- as.data.frame(Embeddings(UMAP_tBM)[,1:2]) # extract the DR from the UMAP
head(UMAP_DR_tBM)

souped_tBM = setClusters(souped_tBM, Idents(UMAP_tBM)) #add cluster ID to the metadata
souped_tBM = setDR(souped_tBM, UMAP_DR_tBM, c("RD1", "RD2")) #add DR to the metadata



Embeddings(UMAP_utBM)[,1:2]
UMAP_DR_utBM <- as.data.frame(Embeddings(UMAP_utBM)[,1:2]) # extract the DR from the UMAP
head(UMAP_DR_utBM)

souped_utBM = setClusters(souped_utBM, Idents(UMAP_utBM)) #add cluster ID to the metadata
souped_utBM = setDR(souped_utBM, UMAP_DR_utBM, c("RD1", "RD2")) #add DR to the metadata

# Estimating the contamination fraction for removal (automated method) (CenGen manually estimated using cell-specific marers and estimated a global contamination fraction of 6.45%)
souped_tBM = autoEstCont(souped_tBM) #estimated ~10% contamination (rho 0.10)

souped_utBM = autoEstCont(souped_utBM) #estimated ~13% contamination (rho 0.13) but when cluster resolution was set to 0.5 like the tBM the soup was estimated to be 0.8...soo......



# We have now calculated the contamination fraction for each cell and can use this to remove the contamination from the original count matrix. 
out_tBM = adjustCounts(souped_tBM)
out_utBM = adjustCounts(souped_utBM)


# Create Seurat objects from the SoupX output
pro_tBM = CreateSeuratObject(out_tBM, min.cells = 3, project = "tBM")
pro_utBM = CreateSeuratObject(out_utBM, min.cells = 3, project = "utBM")





### Begin QC analysis to remove additional low quality cells

# Add mitochondrial representation to the metadata of each data set and view distribution
pro_tBM[["percent.mt"]] <- PercentageFeatureSet(pro_tBM, features = c("WBGene00225418" , "WBGene00225415", "WBGene00220387"))
view(pro_tBM@meta.data)
pro_tBM_metadata <- pro_tBM@meta.data
(tBM_mtrna_dis_plot <- pro_tBM_metadata %>% 
    ggplot(aes(color=orig.ident, x=percent.mt, fill= orig.ident)) + 
    geom_density(alpha = 0.2) + 
    theme_classic() +
    scale_x_continuous(breaks = seq(0,100,10)) +
    ylab("Cell density") +
    NULL)


pro_utBM[["percent.mt"]] <- PercentageFeatureSet(pro_utBM, features = c("WBGene00225418" , "WBGene00225415", "WBGene00220387"))
view(pro_utBM@meta.data)
pro_utBM_metadata <- pro_utBM@meta.data
(utBM_mtrna_dis_plot <- pro_utBM_metadata %>% 
    ggplot(aes(color=orig.ident, x=percent.mt, fill= orig.ident)) + 
    geom_density(alpha = 0.2) + 
    theme_classic() +
    scale_x_continuous(breaks = seq(0,100,10)) +
    ylab("Cell density") +
    NULL)

# Remove cells that express <= 10% mitochondrial representation in downstream subset function

# Quick overview glance via VlnPlot
VlnPlot(pro_tBM, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3) 
VlnPlot(pro_utBM, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3) 

# look at UMIs per cell across the dataset
(tBM_UMI_plot <- pro_tBM_metadata %>% 
    ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) + 
    geom_histogram(binwidth = 10) + 
    theme_classic() +
    xlab("Number of UMIs") +
    NULL)


(utBM_UMI_plot <- pro_utBM_metadata %>% 
    ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) + 
    geom_histogram(binwidth = 10) + 
    theme_classic() +
    xlab("Number of UMIs") +
    NULL)

# Remove cells that express x> 10,000 UMIs in downstream subset function


# look at gene representation for each cell and identify cells that are largely represented by a single or few genes
pro_tBM$Percent.Largest.Gene <- apply(
  tBM@assays$RNA@counts,
  2,
  function(x)(100*max(x))/sum(x))

tBM_metadata <- pro_tBM@meta.data

(lgst_gene_plot <- tBM_metadata %>% 
    ggplot(aes(color=orig.ident, x=Percent.Largest.Gene, fill= orig.ident)) + 
    geom_bar(stat = "count") + 
    theme_classic() +
    scale_x_continuous(breaks = seq(0,80,5)) +
    ylab("Cell Count") +
    NULL)


pro_utBM$Percent.Largest.Gene <- apply(
  utBM@assays$RNA@counts,
  2,
  function(x)(100*max(x))/sum(x))

utBM_metadata <- pro_utBM@meta.data

(lgst_gene_plot <- utBM_metadata %>% 
    ggplot(aes(color=orig.ident, x=Percent.Largest.Gene, fill= orig.ident)) + 
    geom_bar(stat = "count") + 
    theme_classic() +
    scale_x_continuous(breaks = seq(0,80,5)) +
    ylab("Cell Count") +
    NULL)

# Remove cells where a single gene represents >= 15% of the total cell genotype

# Visualize the distribution of genes detected per cell via histogram
(tBM_features <- tBM_metadata %>% 
    ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) + 
    geom_bar(stat = "count", show.legend = TRUE) + 
    theme_classic() +
    scale_x_continuous(limits = c(0,10000)) +
    ylab("Cell Count") +
    xlab("Total genes") +
    NULL)

(utBM_features <- utBM_metadata %>% 
    ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) + 
    geom_bar(stat = "count", show.legend = TRUE) + 
    theme_classic() +
    scale_x_continuous(limits = c(0,10000)) +
    ylab("Cell Count") +
    xlab("Total genes") +
    NULL)


# Combination of reads, genes, and percent.mt for each  cell
ggplot(tBM_metadata, aes( x = nCount_RNA, y = nFeature_RNA, color = percent.mt))+
  geom_point(size = 0.5, position = position_jitter(w = 4, h = 0.5))+
  scale_color_gradient(low = "blue", high = "red")+
  theme_classic()+
  facet_grid(~orig.ident)+
  xlab("Count depth")+
  ylab("Number of genes")+
  NULL


ggplot(utBM_metadata, aes( x = nCount_RNA, y = nFeature_RNA, color = percent.mt))+
  geom_point(size = 0.5, position = position_jitter(w = 4, h = 0.5))+
  scale_color_gradient(low = "blue", high = "red")+
  theme_classic()+
  facet_grid(~orig.ident)+
  xlab("Count depth")+
  ylab("Number of genes")+
  NULL


# Subset each seurat objected on the following parameters:
tBM_subset <- subset(tBM,
                     nCount_RNA < 10000 & 
                       nFeature_RNA < 2500 & 
                       percent.mt < 10 & 
                       Percent.Largest.Gene < 15)
view(tBM_subset@meta.data)

utBM_subset <- subset(utBM,
                      nCount_RNA < 10000 & 
                        nFeature_RNA < 2500 & 
                        percent.mt < 10 & 
                        Percent.Largest.Gene < 15)
view(utBM_subset@meta.data)


# use the tBM_subset and utBM_subsest seurat objects for further normalization and selection of informative genes 


