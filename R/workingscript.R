# install packages
BiocManager::install("scCB2")
install.packages("Seurat")
BiocManager::install("DropletUtils")
BiocManager::install("scDblFinder")
install.packages("SoupX")

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

# see package vignettes
vignette("scCB2")
browseVignettes("scDblFinder")

# fix PSOCK cluster error to make CB2CellFinder work. default is "parallel"
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")


# set directories for each raw ext data set
data.dir <- ("~/Box/ZamanianLab/SeqLibraries/Mapping/singlecell/210518_BH7FNCDRXY/outs_optext/outs/raw_feature_bc_matrix")


#######---------Analysis for 0 bp ext. dataset --------------------

### -------scCB2 analysis

# read in raw data as a dgCMatrix
raw <- Read10xRaw(data.dir)

# run CB2 to distinguish real cells from empty droplets.Change lower=50 over the default lower=100
CBout <- CB2FindCell(raw, FDR_threshold = 0.01, lower = 100, Ncores = 2)
str(assay(CBout)) 
str(metadata(CBout))

# extract the real cell matrix
RealCell <- GetCellMat(CBout)
str(RealCell)

# option to create a Seurat object from the scCB2 filtered cells 
bma<- Seurat::CreateSeuratObject(counts = RealCell, project = "20210413_bma")

view(bma_0bp@meta.data)

# write a csv file containing the count matrix from the Seurat object
write.table(as.matrix(GetAssayData(object = bma_0bp, slot = "counts")), 
            '~/Desktop/raw_0bp_counts.csv', 
            sep = ',', row.names = T, col.names = T, quote = F)

CB2_cells <- colnames(filt_0bp_mxt)

# ------emptyDrops analysis for comparison

# read in data as SingleCellExperiment object
ed_raw <- read10xCounts(data.dir)
colnames(ed_raw) <- ed_raw$Barcode
rowData(ed_raw) <- rowData(ed_raw)[,1:2]
gene_list_0bp <- read.table("~/Box/ZamanianLab/SeqLibraries/Mapping/singlecell/210518_BH7FNCDRXY/outs_0bp_ext/raw_feature_bc_matrix/features.tsv.gz", sep = "\t", stringsAsFactors = F, header = F)

rowData(raw_0bp)$Symbol <- as.character(VLOOKUP(rowData(raw_0bp)$ID, gene_list_0bp))

bcrank <- barcodeRanks(counts(raw_0bp))
uniq <- !duplicated(bcrank$rank)
plot(bcrank$rank[uniq], bcrank$total[uniq], log = "xy", xlab = "Rank", ylab = "Total UMI count", cex_lab = 1.2)

# # Using droplets with fewer than 50 UMIs as the background. Call cells at FDR of 0.1%

set.seed(100)

ed_filtered <- emptyDrops(counts(ed_raw), lower = 100, test.ambient = TRUE)
summary(ed_filtered$FDR <= 0.001, na.rm = TRUE)

# subset raw_xxxbp object to retain only the detected cells using which() to remove NAs
ed_filt <- ed_raw[, which(ed_filtered$FDR <= 0.001)]

#emptydrops() assumes cells with very low UMI counts are empty so we can look at the distribution of p-values for low-total barcodes. Should be uniform distribution


ed_filtered_check <- data.frame(PValue=ed_filtered@listData$PValue, Total=ed_filtered@listData$Total, FDR=ed_filtered@listData$FDR) %>%  # to create data frame for ggplot input
  drop_na() %>% #remove NAs from dataset
 
subset <- filter(ed_filtered_check, Total >0 & Total <= 100) # filter out any rows with greater than 100 UMIs (ambient)

# plot p-values of "ambient" RNA cells
(subset_plot <- subset %>% 
  ggplot(aes(x =PValue)) + 
  geom_histogram() +
  theme_nw()+
  NULL)

# processing to this point distinguishes real cells from empty cells containing ambient RNA but doesnt tell us anything about the quality of the real cells. 

















# add number of genes per UMI for each cell to metadata
bma$GenesPerUMI <- (bma$nFeature_RNA / bma$nCount_RNA)
view(bma@meta.data)




# --------------Begin analysis using Seurat--------------------------------

# What percentage of all counts for each cell belong to mitochondrial genes WBGene00225418, WBGene00225415, WBGene00220387?

bma[["percent.mt"]] <- PercentageFeatureSet(bma, features = c("WBGene00225418" , "WBGene00225415", "WBGene00220387"))
view(bma@meta.data)

metadata <- bma@meta.data

# plot the distribution of percent.mt per cell
(mtrna_dis_plot <- metadata %>% 
    ggplot(aes(color=orig.ident, x=percent.mt, fill= orig.ident)) + 
    geom_density(alpha = 0.2, show.legend = FALSE) + 
    theme_classic() +
    scale_x_continuous(breaks = seq(0,50, 10), expand = c(0,0)) +
    ylab("Cell density") +
    NULL)





# Remove cells that have mitochondrial genes >= 20%
bma_nomt <- bma[,which(bma@meta.data$percent.mt <= 20)]
view(bma_nomt@meta.data)




## looking at QC metrics-------------

metadata_nomt <- bma_nomt@meta.data

metadata_nomt$cells <- rownames(metadata_nomt)

metadata_nomt <- metadata_nomt %>%
  dplyr::rename(nUMI = nCount_RNA,
                nGene = nFeature_RNA)

metadata_nomt$log10GenesPerUMI <- log10(bma_nomt$nFeature_RNA) / log10(bma_nomt$nCount_RNA)

# visualize percent_mt per cell in filtered dataset
(mt_filtered_plot <- metadata_nomt %>% 
    ggplot(aes(color=orig.ident, x=percent.mt, fill= orig.ident)) + 
    geom_density(alpha = 0.2, show.legend = FALSE) + 
    theme_classic() +
    scale_x_continuous(breaks = seq(0,22)) +
    ylab("Cell density") +
    NULL)


# Visualizing # transcripts per cell
(metadata_plot <- metadata_nomt %>% 
    ggplot(aes(color=orig.ident, x=nUMI, fill= orig.ident)) + 
    geom_density(alpha = 0.2, show.legend = FALSE) + 
    theme_classic() +
    scale_x_log10() +
    ylab("Cell density") +
    NULL)

# Visualize the distribution of genes detected per cell via histogram
(metadata_plot2 <- metadata_nomt %>% 
    ggplot(aes(color=orig.ident, x=nGene, fill= orig.ident)) + 
    geom_density(alpha = 0.2, show.legend = FALSE) + 
    theme_classic() +
    scale_x_log10() + 
    NULL)

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
(metadata_plot3 <- metadata_nomt %>% 
    ggplot(aes(x=nUMI, y=nGene)) + 
    geom_point(show.legend = FALSE) + 
    scale_colour_gradient(low = "gray90", high = "black") +
    stat_smooth(method=lm) +
    scale_x_log10() + 
    scale_y_log10() + 
    theme_classic() +
    geom_vline(xintercept = 500) +
    geom_hline(yintercept = 250) +
    facet_wrap(~orig.ident) +
    NULL)


# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI

(metadata_plot4 <- metadata_nomt %>%
    ggplot(aes(x=log10GenesPerUMI, color = orig.ident, fill=orig.ident)) +
    geom_density(alpha = 0.2, show.legend = FALSE) +
    theme_classic() +
    geom_vline(xintercept = 0.8)+
    NULL)








## Visualize QC metrics as a violin plot
VlnPlot(bma, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)

#plotting different features with correlation as the plot title 
FeatureScatter(bma, feature1 = "nCount_RNA", feature2 = "percent.mt") 
FeatureScatter(bma, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 



## Remove cells that have over 5,000 unique feature counts. Filters 39500 to 37066. Super conservative! Seurat recommends 200-2,500 range. 
bma<- subset(bma, subset = nFeature_RNA > 100 & nFeature_RNA < 5000)



## Normalizing the datausing a global-scaling normalization method "LogNormalize"  
bma <- NormalizeData(bma, normalization.method = "LogNormalize", scale.factor = 10000)



## Identification of highly variable features (feature selection)
bma <- FindVariableFeatures(bma, selection.method = "vst", nfeatures = 2000)
## Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(bma), 10)

## plot variable features with and without labels
plot1 <- VariableFeaturePlot(bma)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

## scaling the data
all.genes <- rownames(bma)
bma <- ScaleData(bma, features = all.genes)

##PCA on the scaled data
bma <- RunPCA(bma, features = VariableFeatures(object = bma))

## Visualization
PCA <- DimPlot(bma, reduction  = "pca")



write.table(as.matrix(GetAssayData(object = bma_0bp, slot = "counts")), 
            '~/Desktop/raw_0bp_counts.csv', 
            sep = ',', row.names = T, col.names = T, quote = F)
