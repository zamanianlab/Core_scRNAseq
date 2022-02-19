## Plotting targets of interest within the ut/tBM integrated dataset

# load required libraries
library(ZamanianLabThemes)
library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(pals)
library(ggtext)


# read in combined ut/t BM seurat object
BM_combined <- readRDS("/Users/clairhenthorn/Library/CloudStorage/Box-Box/ZamanianLab/LabMembers/Clair/project_singlecell/other/BM_combined.RDS")


# pull out the normalized counts, metadata, and UMAP coordinates into dataframes for plotting in ggplot2
data <- as_tibble(BM_combined@reductions$umap@cell.embeddings, rownames = 'index') # UMAP coordinates for each cell

md <- as_tibble(BM_combined@meta.data, rownames = 'index') # metadata detailing ut/t identity and cluster information

counts <- as_tibble(BM_combined@assays[["RNA"]]@data, rownames = "gene_id") %>%  # gene expression matrix of normalized counts
  pivot_longer(!gene_id, names_to = "index", values_to = "counts") %>%
  filter(counts > 0)


# Glutamate-gated chloride channels (GluCls)
#subset data to only glucls
glucls <- counts %>%
  filter(gene_id == c("WBGene00221971", "WBGene00222703", "WBGene00223839", "WBGene00228311")) %>%
  left_join(data)

#plot
glucls %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2))+
  geom_point(data = data, size = 0.5, alpha = 0.1, color = "grey")+
  geom_point(aes(color = gene_id), size = 0.7)+
  scale_color_manual(values = c("#f8d364", "purple", "#1c493d", "#c22e43"), name = "Markers", labels = c("*avr-14*", "*glc-4*", "*glc-2*", "*glc-3*"))+
  labs(title = "GluCls")+
  theme(legend.text = element_markdown(),
        panel.background = element_blank(),
        axis.line = element_line (colour = "black"),
        legend.background=element_blank(),
        legend.key = element_blank())+
  guides(color = guide_legend(override.aes = list(size=1.5)))+
  NULL



  ## coexpression of GluCl subunits

  #identify what barcodes have more than one gene assigned
  tmp<- glucls[duplicated(glucls[4]),]
  tmp <- tmp$index

  # use the output to subset glucl to include only the duplicates
  glucl_co <- glucls[glucls$index %in% tmp,] %>%
    select(-UMAP_1, -UMAP_2)
           
# plot
  glucl_co %>%
    ggplot(aes(x = gene_id, y = index, fill = counts))+
    geom_tile()+
  scale_x_discrete(labels = c("*avr-14*", "*glc-4*", "*glc-3*"))+
    labs( x = "Genes", y = "Cells", title = "GluCl Subunit Co-expression")+
    theme(axis.text.x = element_markdown())
  


## Betatubulins 
#create dataframe subset to only the btubs
btubs <- counts %>%
  filter(gene_id == c("WBGene00229959", "WBGene00224994", "WBGene00228922", "WBGene00233027")) %>%
  left_join(data)

#plot
btubs %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2))+
  geom_point(data = data, size = 0.5, alpha = 0.1, color = "grey")+
  geom_point(data = subset(btubs, gene_id == list("WBGene00224994", "WBGene00228922", "WBGene00233027")), aes(color = gene_id), size = 0.7, position = "jitter")+
  geom_point(data = subset(btubs, gene_id == "WBGene00229959" & counts >= 1), aes(color = gene_id), size = 0.6, position = "jitter")+
  scale_color_manual(values = c("#f8d364", "blue", "#1c493d", "red"), name = "Markers", labels = c("*ben-1*", "*mec-7*", "*tbb-2*", "*tbb-4*"))+
  labs(title = "Betatubulins")+
  theme(legend.text = element_markdown(),
        panel.background = element_blank(),
        axis.line = element_line (colour = "black"),
        legend.background=element_blank(),
        legend.key = element_blank())+
  guides(color = guide_legend(override.aes = list(size=1.5)))+
  NULL



#identify what barcodes have more than one gene assigned
tmp<- btubs[duplicated(btubs[4]),]
tmp <- tmp$index

# use the output to subset glucl to include only the duplicates
btubs_co <- btubs[btubs$index %in% tmp,] %>%
  select(-UMAP_1, -UMAP_2)

#plot
btubs_co %>%
  ggplot(aes(x = gene_id, y = index, fill = counts))+
  geom_tile()+
  scale_x_discrete(labels = c("*ben-1*", "*mec-7*", "*tbb-2*", "*tbb-4*"))+
  labs( x = "Genes", y = "Cells", title = "Betatubulin Co-expression")+
  theme(axis.text.x = element_markdown())


# slo-1 (Emodepside target)
slo1 <- counts %>%
  filter(gene_id == "WBGene00226980") %>%
  left_join(data)


slo1 %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2))+
  geom_point(data = data, size = 0.5, alpha = 0.1, color = "grey")+
  geom_point(data = subset(slo1, gene_id == "WBGene00226980" & counts >= 1), aes(color = gene_id), size = 0.7, position = "jitter")+
  scale_color_manual(values = "blue", name = "Marker", labels = "*slo-1*")+
  labs(title = "Emodepside Target")+
  theme(legend.text = element_markdown(),
        panel.background = element_blank(),
        axis.line = element_line (colour = "black"),
        legend.background=element_blank(),
        legend.key = element_blank())+
  guides(color = guide_legend(override.aes = list(size=1.5)))+
  NULL



### Combined plots of glucls, btubs, and slo-1
#create column in each dataframe with the target label
glucls$target <- "glucls"
btubs$target <- "btubs"
slo1$target <- "slo1"

# bind all 3 target dataframes and append the metadata for cluster analysis
targets <- rbind(glucls, slo1, btubs) %>%
  left_join(md)
  

# plot UMAP with targets as separate colors
targets %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2))+
  geom_point(data = data, size = 0.5, alpha = 0.1, color = "grey")+
  geom_point(data = subset(targets, target == "btubs"), aes(color = target))+
  geom_point(data = subset(targets, target == "slo1"), aes(color = target))+
  geom_point(data = subset(targets, target == "glucls"), aes(color = target))+
  scale_size(breaks = c(1, 4, 8, 12, 16))+
  theme(legend.text = element_markdown(),
        panel.background = element_blank(),
        axis.line = element_line (colour = "black"),
        legend.background=element_blank(),
        legend.key = element_blank())+
  guides(color = guide_legend(override.aes = list(size=1.5)))+
  NULL




## Expression of btubs, glucls, and slo1 by dotplot instead of UMAP

tmp <- c("WBGene00226980","WBGene00229959", "WBGene00224994", "WBGene00228922", "WBGene00233027", "WBGene00221971", "WBGene00222703", "WBGene00223839", "WBGene00228311")

# extract percent expressing of each gene in the assigned cluster from dotplot and join to the targets dataframe
a <- DotPlot(BM_combined, features = tmp, assay = "RNA") #seurat dotplot for extraction of the percent cells expressing in each cluster


dotplot<- a$data %>%
  rownames_to_column(var = "gene_id") %>%
  select(c("gene_id", "pct.exp")) %>%
  subset(gene_id == tmp) %>%
  left_join(targets)

# create the dotplot
dotplot %>%
  ggplot(aes(y = integrated_snn_res.0.5, x = gene_id, color = pct.exp, size= counts))+
  geom_point()+
  scale_x_discrete(labels = c("*avr-14*", "*glc-4*", "*glc-2", "*ben-1*", "*slo-1*", "*glc-3*",  "*mec-7*", "*tbb-2*", "*tbb-4*"))+
  theme_classic()+
  coord_flip()
        
  

########################################
## Complete Drug Targets
########################################
drugs <- read.csv("/Users/clairhenthorn/Desktop/drug_targets.csv")


# plotting GPCR targets
gpcr <- subset(drugs, type == "GPCR") %>% 
  left_join(counts) %>% 
  left_join(data)
  
  
gpcr %>% 
  ggplot(aes(x = UMAP_1, y = UMAP_2))+
  geom_point(data= data, size = 0.5, alpha = 0.1, color = "grey")+
  geom_point(aes(color = gene_id))+
  theme_classic()







