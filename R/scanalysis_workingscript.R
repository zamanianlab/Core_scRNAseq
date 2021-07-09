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
data.dir.0 <- ("~/Box/ZamanianLab/SeqLibraries/Mapping/singlecell/210518_BH7FNCDRXY/outs_0bp_ext/raw_feature_bc_matrix/")
data.dir.50 <- ("~/Box/ZamanianLab/SeqLibraries/Mapping/singlecell/210518_BH7FNCDRXY/outs_50bp_ext/raw_feature_bc_matrix/")
data.dir.100 <- ("~/Box/ZamanianLab/SeqLibraries/Mapping/singlecell/210518_BH7FNCDRXY/outs_100bp_ext/raw_feature_bc_matrix/")
data.dir.150 <- ("~/Box/ZamanianLab/SeqLibraries/Mapping/singlecell/210518_BH7FNCDRXY/outs_150bp_ext/raw_feature_bc_matrix/")
data.dir.200 <- ("~/Box/ZamanianLab/SeqLibraries/Mapping/singlecell/210518_BH7FNCDRXY/outs_200bp_ext/raw_feature_bc_matrix/")
data.dir.250 <- ("~/Box/ZamanianLab/SeqLibraries/Mapping/singlecell/210518_BH7FNCDRXY/outs_250bp_ext/raw_feature_bc_matrix/")
data.dir.300 <- ("~/Box/ZamanianLab/SeqLibraries/Mapping/singlecell/210518_BH7FNCDRXY/outs_300bp_ext/raw_feature_bc_matrix/")
data.dir.400 <- ("~/Box/ZamanianLab/SeqLibraries/Mapping/singlecell/210518_BH7FNCDRXY/outs_400bp_ext/raw_feature_bc_matrix/")
data.dir.500 <- ("~/Box/ZamanianLab/SeqLibraries/Mapping/singlecell/210518_BH7FNCDRXY/outs_500bp_ext/raw_feature_bc_matrix/")


#######---------Analysis for 0 bp ext. dataset --------------------

### -------scCB2 analysis

# read in raw data as a dgCMatrix
raw_0 <- Read10xRaw(data.dir.0)

# run CB2 to distinguish real cells from empty droplets.Change lower=50 over the default lower=100
CBout_0bp <- CB2FindCell(raw_0, FDR_threshold = 0.01, lower = 50, Ncores = 2)
str(assay(CBout_0bp)) 
str(metadata(CBout_0bp))

# extract the real cell matrix
RealCell_0bp <- GetCellMat(CBout_0bp)
str(RealCell_0bp)

# option to create a Seurat object from the scCB2 filtered cells 
bma_0bp<- Seurat::CreateSeuratObject(counts = RealCell_0bp, project = "20210413_bma")

view(bma_0bp@meta.data)

# write a csv file containing the count matrix from the Seurat object
write.table(as.matrix(GetAssayData(object = bma_0bp, slot = "counts")), 
            '~/Desktop/raw_0bp_counts.csv', 
            sep = ',', row.names = T, col.names = T, quote = F)

CB2_cells <- colnames(filt_0bp_mxt)

# ------emptyDrops analysis for comparison

# read in data as SingleCellExperiment object
raw_0bp <- read10xCounts(data.dir.0)
colnames(raw_0bp) <- raw_0bp$Barcode
rowData(raw_0bp) <- rowData(raw_0bp)[,1:2]
gene_list_0bp <- read.table("~/Box/ZamanianLab/SeqLibraries/Mapping/singlecell/210518_BH7FNCDRXY/outs_0bp_ext/raw_feature_bc_matrix/features.tsv.gz", sep = "\t", stringsAsFactors = F, header = F)

rowData(raw_0bp)$Symbol <- as.character(VLOOKUP(rowData(raw_0bp)$ID, gene_list_0bp))

bcrank <- barcodeRanks(counts(raw_0bp))
uniq <- !duplicated(bcrank$rank)
plot(bcrank$rank[uniq], bcrank$total[uniq], log = "xy", xlab = "Rank", ylab = "Total UMI count", cex_lab = 1.2)

# # Using droplets with fewer than 50 UMIs as the background. Call cells at FDR of 0.1%

set.seed(100)

ed.filtered_0bp <- emptyDrops(counts(raw_0bp), lower = 50)
summary(ed.filtered_0bp$FDR <= 0.001, na.rm = TRUE)

# subset raw_xxxbp object to retain only the detected cells using which() to remove NAs
raw_0bp <- raw_0bp[, which(filtered_0bp$FDR <= 0.001)]

#emptydrops() assumes cells with very low UMI counts are empty so we can look at the distribution of p-values for low-total barcodes. Should be uniform distribution


filtered_500bp_check = data.frame(PValue=filtered_500bp@listData$PValue, Total=filtered_500bp@listData$Total, FDR=filtered_500bp@listData$FDR) %>%  # to create data frame for ggplot input
  drop_na() %>% #remove NAs from dataset
  subset<- filter(filtered_500bp_check, FDR > 0 & FDR <= 0.001)

??????
  
  ggplot(subset, aes(x = PValue)) + geom_bar()




#######---------Analysis for 50 bp ext. dataset --------------------

### -------scCB2 analysis

# read in raw data as a dgCMatrix
raw_50 <- Read10xRaw(data.dir.50)

# run CB2 to distinguish real cells from empty droplets.Change lower=50 over the default lower=100
CBout_50bp <- CB2FindCell(raw_50, FDR_threshold = 0.01, lower = 50, Ncores = 2)
str(assay(CBout_50bp)) 
str(metadata(CBout_50bp))

# extract the real cell matrix
RealCell_50bp <- GetCellMat(CBout_50bp)
str(RealCell_50bp)

# option to create a Seurat object from the scCB2 filtered cells 
bma_50bp<- Seurat::CreateSeuratObject(counts = RealCell_50bp, project = "20210413_bma")

view(bma_50bp@meta.data)


# ------emptyDrops analysis for comparison

# read in data as SingleCellExperiment object
raw_50bp <- read10xCounts(data.dir.50)
colnames(raw_50bp) <- raw_50bp$Barcode
rowData(raw_50bp) <- rowData(raw_50bp)[,1:2]
gene_list_50bp <- read.table("~/Box/ZamanianLab/SeqLibraries/Mapping/singlecell/210518_BH7FNCDRXY/outs_50bp_ext/raw_feature_bc_matrix/features.tsv.gz", sep = "\t", stringsAsFactors = F, header = F)

rowData(raw_50bp)$Symbol <- as.character(VLOOKUP(rowData(raw_50bp)$ID, gene_list_50bp))

bcrank <- barcodeRanks(counts(raw_50bp))
uniq <- !duplicated(bcrank$rank)
plot(bcrank$rank[uniq], bcrank$total[uniq], log = "xy", xlab = "Rank", ylab = "Total UMI count", cex_lab = 1.2)

# # Using droplets with fewer than 50 UMIs as the background. Call cells at FDR of 0.1%

set.seed(100)

ed.filtered_50bp <- emptyDrops(counts(raw_50bp), lower = 50)
summary(ed.filtered_50bp$FDR <= 0.001, na.rm = TRUE)

# subset raw_xxxbp object to retain only the detected cells using which() to remove NAs
raw_50bp <- raw_50bp[, which(filtered_50bp$FDR <= 0.001)]

#emptydrops() assumes cells with very low UMI counts are empty so we can look at the distribution of p-values for low-total barcodes. Should be uniform distribution

#######---------Analysis for 100 bp ext. dataset --------------------

### -------scCB2 analysis

# read in raw data as a dgCMatrix
raw_100 <- Read10xRaw(data.dir.100)

# run CB2 to distinguish real cells from empty droplets.Change lower=50 over the default lower=100
CBout_100bp <- CB2FindCell(raw_100, FDR_threshold = 0.01, lower = 50, Ncores = 2)
str(assay(CBout_100bp)) 
str(metadata(CBout_100bp))

# extract the real cell matrix
RealCell_100bp <- GetCellMat(CBout_100bp)
str(RealCell_100bp)

# option to create a Seurat object from the scCB2 filtered cells 
bma_100bp<- Seurat::CreateSeuratObject(counts = RealCell_100bp, project = "20210413_bma")

view(bma_100bp@meta.data)


# ------emptyDrops analysis for comparison

# read in data as SingleCellExperiment object
raw_100bp <- read10xCounts(data.dir.100)
colnames(raw_100bp) <- raw_100bp$Barcode
rowData(raw_100bp) <- rowData(raw_100bp)[,1:2]
gene_list_100bp <- read.table("~/Box/ZamanianLab/SeqLibraries/Mapping/singlecell/210518_BH7FNCDRXY/outs_100bp_ext/raw_feature_bc_matrix/features.tsv.gz", sep = "\t", stringsAsFactors = F, header = F)

rowData(raw_100bp)$Symbol <- as.character(VLOOKUP(rowData(raw_100bp)$ID, gene_list_100bp))

bcrank <- barcodeRanks(counts(raw_100bp))
uniq <- !duplicated(bcrank$rank)
plot(bcrank$rank[uniq], bcrank$total[uniq], log = "xy", xlab = "Rank", ylab = "Total UMI count", cex_lab = 1.2)

# # Using droplets with fewer than 50 UMIs as the background. Call cells at FDR of 0.1%

set.seed(100)

ed.filtered_100bp <- emptyDrops(counts(raw_100bp), lower = 50)
summary(ed.filtered_100bp$FDR <= 0.001, na.rm = TRUE)

# subset raw_xxxbp object to retain only the detected cells using which() to remove NAs
raw_100bp <- raw_100bp[, which(filtered_100bp$FDR <= 0.001)]

#emptydrops() assumes cells with very low UMI counts are empty so we can look at the distribution of p-values for low-total barcodes. Should be uniform distribution















# add number of genes per UMI for each cell to metadata
bma$GenesPerUMI <- (bma$nFeature_RNA / bma$nCount_RNA)
view(bma@meta.data)


#calculate barcode rank and plot for raw and filtered
## knee plot of raw data

empty_thresh <- 100
bc_ranks <- barcodeRanks(counts(raw_50bp), lower = empty_thresh)

colData(raw_50bp)$BarcodeRank   <- bc_ranks$rank
colData(raw_50bp)$BarcodeTotal  <- bc_ranks$total
colData(raw_50bp)$BarcodeFitted <- bc_ranks$fitted

bc_data <- colData(raw_50bp) %>%
  as.data.frame() %>%
  select(Cell, Kept = CellRangerFilt, Rank = BarcodeRank,
         Total = BarcodeTotal, Fitted = BarcodeFitted) %>%
  arrange(Rank)

ggplot(bc_data, aes(x = Rank, y = Total)) +
  geom_point(shape = 1, aes(colour = Kept)) +
  geom_line(aes(y = Fitted), colour = "red") +
  geom_hline(yintercept = bc_ranks$knee,
             colour = "dodgerblue", linetype = "dashed") +
  annotate("text", x = 0, y = bc_ranks$knee, label = "Knee",
           colour = "dodgerblue", hjust = 0, vjust = -1) +
  geom_hline(yintercept = bc_ranks$inflection,
             colour = "forestgreen", linetype = "dashed") +
  annotate("text", x = 0, y = bc_ranks$inflection, label = "Inflection",
           colour = "forestgreen", hjust = 0, vjust = -1) +
  geom_hline(yintercept = empty_thresh,
             colour = "darkorchid", linetype = "dashed") +
  annotate("text", x = 0, y = empty_thresh, label = "Empty threshold",
           colour = "darkorchid", hjust = 0, vjust = -1) +
  scale_x_log10(labels = scales::number) +
  scale_y_log10(labels = scales::number) +
  scale_colour_manual(values = c("black", "violet")) +
  ylab("Total counts") +
  theme_minimal()

## knee plot of filtered RealCell data from scCB2
br.out.filt<- barcodeRanks(RealCell)

(br.out.filt_plot <-plot(br.out.filt$rank, br.out.filt$total, log="xy", xlab="Rank", ylab="Total") %>% 
    abline(h=metadata(br.out.filt)$knee, col="dodgerblue", lty=2)+
    abline(h=metadata(br.out.filt)$inflection, col="forestgreen", lty=2)+
    legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"),
           legend=c("knee", "inflection"))+
    NULL)



# --------------Begin analysis using Seurat--------------------------------

# What percentage of all counts for each cell belong to mitochondrial genes WBGene00225418, WBGene00225415, WBGene00220387?

bma[["percent.mt"]] <- PercentageFeatureSet(bma, features = c("WBGene00225418" , "WBGene00225415", "WBGene00220387"))
view(bma@meta.data)

# Remove cells that have mitochondrial genes 


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
