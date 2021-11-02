# Using C. elegans orthologs to define cluster cell types and identify 

# load libraries
library(tidyverse)
library(ggplot2)
library(cowplot)
library(ZamanianLabThemes)



# import DEG matrix for each treatment with the C. elegans orthologs appended to the table
tBM_markers_ortho <- readRDS("~/Box/ZamanianLab/SeqLibraries/Mapping/singlecell/scRNAesq_Bma_multi/RDS_objects/scmulti_tBM_markers_ortho.rds")  # 14,863 marker genes in total
utBM_markers_ortho <- readRDS("~/Box/ZamanianLab/SeqLibraries/Mapping/singlecell/scRNAesq_Bma_multi/RDS_objects/scmulti_utBM_markers_ortho.rds") # 13,804 marker genes in total


# create new column in dataframe that corresponds to the difference of pct.1 and pct.2 for each marker
tBM_markers_ortho$pct.df <- tBM_markers_ortho$pct.1 - tBM_markers_ortho$pct.2
utBM_markers_ortho$pct.df <- utBM_markers_ortho$pct.1 - utBM_markers_ortho$pct.2


# remove rows that do not have a Bma_gene_name (NA)
tBM_clean <- tBM_markers_ortho[!is.na(tBM_markers_ortho$Bma_gene_name),]  # 4,381 marker genes remaining
utBM_clean <- utBM_markers_ortho[!is.na(utBM_markers_ortho$Bma_gene_name),] # 4,010 marker genes remaining


# create a new matrix with cluster ID as rows and bma gene names as columns. The number in the corresponding box will be the difference between pct.1 and pct.2 indicating the significane of cells where the gene is detected in the cluster indicated. Where Bma_gene_names is duplicated due to multiple C. elegans orthologs


## Generate matrix using C. elegans orthologs
# wide matrix
tBM_cel_wide <- tBM_clean %>% 
  select(-Bma_gene_id, -p_val, -avg_log2FC, -pct.1, -pct.2, -p_val_adj, -gene, -Bma_gene_name) %>%
  group_by(tBM_clean$cel_ortho) %>% 
  data.frame(cel_ortho = make.unique(as.character(tBM_clean$cel_ortho))) %>% 
  ungroup() %>% 
  pivot_wider(id_cols = "cluster", names_from = "cel_ortho", values_from = "pct.df", values_fn = list, names_sort = TRUE) %>% 
  subset("pct.df" >= 0.2)

#long matrix
tBM_cel_long <- tBM_clean %>% 
  group_by(tBM_clean$cel_ortho) %>% 
  data.frame(cel_ortho = make.unique(as.character(tBM_clean$cel_ortho))) %>%
  select(-Bma_gene_id, -p_val, -avg_log2FC, -pct.1, -pct.2, -p_val_adj, -gene, -Bma_gene_name, -tBM_clean.cel_ortho) 
  subset(pct.df >= 0.2)


  
## utBM matrices  
  # wide matrix
  utBM_new_cel <- utBM_clean %>% 
    select(-Bma_gene_id, -p_val, -avg_log2FC, -pct.1, -pct.2, -p_val_adj, -gene, -Bma_gene_name) %>%
    group_by(utBM_clean$cel_ortho) %>% 
    data.frame(cel_ortho = make.unique(as.character(utBM_clean$cel_ortho))) %>% 
    ungroup() %>% 
    pivot_wider(id_cols = "cluster", names_from = "cel_ortho", values_from = "pct.df", values_fn = list, names_sort = TRUE) %>% 
    subset(pct.df >= 0.5)
  
  #long for heatmap
  utBM_hm_cel <- utBM_clean %>% 
    group_by(utBM_clean$cel_ortho) %>% 
    data.frame(cel_ortho = make.unique(as.character(utBM_clean$cel_ortho))) %>%
    select(-Bma_gene_id, -p_val, -avg_log2FC, -pct.1, -pct.2, -p_val_adj, -gene, -Bma_gene_name, -utBM_clean.cel_ortho)
  
  subset(pct.df >= 0.15)
  
  
 
## Using cengen marker lists to make inferences about our cell type clusters  
 
  # import cengen marker list for major cell types and neurons independently
  cengen_marker <- read.csv("~/Desktop/cengen_markers.csv")
  cengen_neurons <- read.csv("~/Desktop/cengen_neurons.csv")
  

  
  
  
  
  
  
  
  

# plots 
  
# tBM heatmap plot
pdf("~/Desktop/tBM_hm_cel.pdf", width = 24)
tBM_cel_long %>% 
  ggplot(aes(cel_ortho.1, cluster)) +
  geom_tile(aes(fill = pct.df))+
  scale_fill_gradient(low = "blue", high = "red") +
  theme(axis.text.x = element_text(angle = 90))
dev.off()


#utBM hm plot
pdf("~/Desktop/utBM_hm_cel.pdf")
utBM_hm_cel %>% 
  ggplot(aes(cel_ortho.1, cluster)) +
  geom_tile(aes(fill = pct.df) +
  scale_fill_gradient(low = "blue", high = "red") +
  theme(axis.text.x = element_text(angle = 90, size=6))
dev.off()










# save RDS objects
saveRDS(tBM_clean, "~/Box/ZamanianLab/SeqLibraries/Mapping/singlecell/scRNAesq_Bma_multi/RDS_objects/scmulti_tBM_clean.rds")
saveRDS(utBM_clean, "~/Box/ZamanianLab/SeqLibraries/Mapping/singlecell/scRNAesq_Bma_multi/RDS_objects/scmulti_utBM_clean.rds")

# objects for import
tBM_clean <- readRDS("~/Box/ZamanianLab/SeqLibraries/Mapping/singlecell/scRNAesq_Bma_multi/RDS_objects/scmulti_tBM_clean.rds")
utBM_clean <- readRDS("~/Box/ZamanianLab/SeqLibraries/Mapping/singlecell/scRNAesq_Bma_multi/RDS_objects/scmulti_utBM_clean.rds")

