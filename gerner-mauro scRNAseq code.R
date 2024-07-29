library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
library(patchwork)
library(cowplot)
library(dplyr)
library(VGAM)
library(grid)
library(EnhancedVolcano)
library(SeuratWrappers)
library(scCustomize)

#####FIGURE 4#####
#object to make fig 4A
E13_2<-readRDS(file="/Users/knger/Downloads/E13_Chicken_Analysis/pbmc.R")
DefaultAssay(E13_2) <- "RNA"




pbmc <- E13_2
pbmc <- FindVariableFeatures(pbmc)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA (pbmc, npcs=30) %>% RunUMAP (dims=1:30) %>% FindNeighbors() %>% FindClusters()
VlnPlot(pbmc, features = c("nFeature_RNA"))
VlnPlot(pbmc, features = c("nCount_RNA"))
DimPlot(pbmc, label=F)
DefaultAssay(pbmc) <- "RNA"
VlnPlot(pbmc, features = c("CDH1","CDH5","PTPRC","COL3A1","PCNA"))

VlnPlot(pbmc, features = c("nFeature_RNA"))

#remove cluster 9 (low gen expression)
pbmc <- subset(pbmc, idents=c(0:8,10:15))
pbmc <- FindVariableFeatures(pbmc)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA (pbmc, npcs=30) %>% RunUMAP (dims=1:30) %>% FindNeighbors() %>% FindClusters()
DimPlot(pbmc, label=F)
DefaultAssay(pbmc) <- "RNA"
VlnPlot(pbmc, features = c("CDH1","CDH5","PTPRC","COL3A1","PCNA"))

#endo
Cluster_Highlight_Plot(seurat_object = pbmc, cluster_name = c("11"), highlight_color = "#6e016b",
                       background_color = "lightgray")

#immune
Cluster_Highlight_Plot(seurat_object = pbmc, cluster_name = c("8","12"), highlight_color = "#cb181d",
                       background_color = "lightgray")

#epi 
Cluster_Highlight_Plot(seurat_object = pbmc, cluster_name = c("0","1","2","5","6","7"), highlight_color = "#41ab5d",
                       background_color = "lightgray")
#mes
Cluster_Highlight_Plot(seurat_object = pbmc, cluster_name = c("3","4","9","10","13","14"), highlight_color = "#225ea8",
                       background_color = "lightgray")


pbmc@meta.data$cell_type <- pbmc@meta.data$seurat_clusters

pbmc@meta.data$cell_type <- plyr::mapvalues(x = pbmc@meta.data$cell_type,
                                            from = c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14"),
                                            to = c("epi","epi","epi","mes","mes","epi","epi","epi","imm","epi","mes","mes","endo","imm","mes"))

Idents(pbmc) <- "cell_type"
table(pbmc@active.ident)

DimPlot(pbmc, label=F, group.by = "cell_type")
DimPlot_scCustom(pbmc, reduction="umap", colors_use = c("#41ab5d","#225ea8","#cb181d","#6e016b"),label = F)

Stacked_VlnPlot(seurat_object = pbmc, features = c("CDH1","COL3A1","PTPRC","CDH5"), x_lab_rotate = TRUE,
                colors_use = c("#41ab5d","#225ea8","#cb181d","#6e016b"), split.by = "cell_type")


#read epithelial only data for the rest of figure 4
E13<-readRDS(file="/Users/knger/Downloads/E13_Chicken_Analysis/epi.R")
DefaultAssay(E13) <- "RNA"

epi <- E13
epi <- FindVariableFeatures(epi)
epi <- ScaleData(epi)
epi <- RunPCA (epi, npcs=30) %>% RunUMAP (dims=1:30) %>% FindNeighbors() %>% FindClusters()
VlnPlot(epi, features = c("nFeature_RNA"))
VlnPlot(epi, features = c("nCount_RNA"))
DimPlot(epi, label=T)
DefaultAssay(epi) <- "RNA"
VlnPlot(epi, features = c("CDH1","CDH5","PTPRC","COL3A1","PCNA"))
VlnPlot(epi, features = c("SOX9","SOX2","NKX2-1","ASCL1", "CCNO", "KRT14"))

#label experiment 
epi$geno <- substr(colnames(epi@assays$RNA@counts), 18, 18) #cell barcode

epi$geno <- plyr::mapvalues(x = epi$geno,
                            from = "1",
                            to = "E13 Chicken")

#remove cluster 7 (low gen expression of CDH1 and NKX)
epi <- subset(epi, idents=c(0:6,8:10))
epi <- FindVariableFeatures(epi)
epi <- ScaleData(epi)
epi <- RunPCA (epi, npcs=30) %>% RunUMAP (dims=1:30) %>% FindNeighbors() %>% FindClusters()
DimPlot(epi, label=T)
DefaultAssay(epi) <- "RNA"
VlnPlot(epi, features = c("CDH1","CDH5","PTPRC","COL3A1","PCNA"))
VlnPlot(epi, features = c("SOX9","SOX2","NKX2-1","FOXJ1", "CCNO", "KRT14"))

#remove 9 because low gene expression 
epi <- subset(epi, idents=c(0:8,10))
epi <- FindVariableFeatures(epi)
epi <- ScaleData(epi)
epi <- RunPCA (epi, npcs=30) %>% RunUMAP (dims=1:30) %>% FindNeighbors() %>% FindClusters(resolution=1.5)
#DimPlot(E15, label=F, group.by = "orig.ident")
DimPlot(epi, label=T)
DefaultAssay(epi) <- "RNA"
VlnPlot(epi, features = c("CDH1","CDH5","GAPDH","COL3A1","PCNA"))
VlnPlot(epi, features = c("SOX9","SOX2","NKX2-1","FOXJ1", "CCNO", "KRT14"))

#label cell types
epi@meta.data$cell_type <- epi@meta.data$seurat_clusters

epi@meta.data$cell_type <- plyr::mapvalues(x = epi@meta.data$cell_type,
                                           from = c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16"),
                                           to = c("SOX9","SOX9","SOX9","SOX9","SOX9","SOX9","Polif","SOX9","SOX9","SOX9","KRT14","Polif","SOX9","SOX9","KRT14","SOX9","Ciliated"))

E13 <- epi
Idents(E13) <- "cell_type"
table(E13@active.ident)
DimPlot(epi, label=F, group.by = "cell_type")



#make umap for cell types for fig 4B
DimPlot_scCustom(E13, reduction="umap", colors_use = c("#0cb702","#0571b0","#ca0020","#ed68ed"),label = F)

#make heatmap for fig 4C
Idents(E13) <- "cell_type"
top25_13 <- FindAllMarkers(E13, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.05)
tmp_13 <- group_by(top25_13, cluster) %>% top_n(n=20, wt=avg_log2FC)
DoHeatmap(E13, features = tmp_13$gene, group.by = "cell_type", group.colors = c("#0cb702","#0571b0","#ca0020","#ed68ed")) + NoLegend()

#make featureplots for fig 4D
FeaturePlot_scCustom(seurat_object = epi, colors_use = viridis_inferno_dark_high, features = "TOP2A")

FeaturePlot(epi, features="NKX2-1", split.by="orig.ident", cols = c("lightgrey", "red"))

#make scatterplot for fig 4E
# Load libraries
library(tidyverse)
library(cowplot)
library(edgeR)
library(Matrix)
library(reshape2)
library(S4Vectors)
library(SingleCellExperiment)
library(pheatmap)
library(png)
library(DESeq2)
library(RColorBrewer)
library(data.table)
library(Seurat)
library(SeuratObject)

###SOX9 chicken vs mouse ######
chicken <- readRDS("/Users/rnayak/Downloads/E13 2.rds")
mouse <- readRDS("/Users/rnayak/Downloads/E15.rds")
chickmouse.list <- list(chicken, mouse)
chickmouse.features <- SelectIntegrationFeatures(object.list = chickmouse.list)
chickmouse.anchors <- FindIntegrationAnchors(object.list = chickmouse.list, anchor.features = chickmouse.features)
#integrate
chickmouse <- IntegrateData(anchorset = chickmouse.anchors)
DefaultAssay(chickmouse) <- "integrated"
chickmouse <- ScaleData(chickmouse) %>% RunPCA(npcs = 30) %>% RunUMAP(dims=1:20) %>% FindNeighbors(dims=1:20) %>% FindClusters(resolution = 0.1)

chickmouse <- readRDS("/Users/rnayak/OneDrive - UTHealth Houston/kamryn/chickmouse.rds")
DefaultAssay(chickmouse) <- "RNA"

# Aggregate expression data by cell type and genotype ## ln CPM from here  https://github.com/satijalab/seurat/issues/2496
avg_chickmouse_expr <- AggregateExpression(chickmouse, group.by = c("cell_type", "geno"), 
                                           assays = "RNA", normalization.method = "LogNormalize", 
                                           scale.factor = 1e6) 

# Convert aggregated expression to a Seurat object
avg_chickmouse_expr_sobj <- AggregateExpression(chickmouse, group.by = c("cell_type", "geno"), 
                                                assays = "RNA", normalization.method = "LogNormalize", 
                                                scale.factor = 1e6, return.seurat = TRUE)

log_data <- GetAssayData(avg_chickmouse_expr_sobj)
#log2_data <- log(expm1(log_data) + 1, 2)

# Extract average expression for chicken and mouse SOX9 cells 
###calculated for all cell types so you can just start from here for other cells 
chicken_sox9_expr <- log_data[,"SOX9_E13 Chicken"]
mouse_sox9_expr <- log_data[,"SOX9_E15 Mouse"]

# Combine into a data frame 
chickmouse_sox9_exp_df <- data.frame(
  gene = rownames(avg_chickmouse_expr_sobj$RNA),
  chicken_exp = as.numeric(chicken_sox9_expr),
  mouse_exp = as.numeric(mouse_sox9_expr)
)

#remove rows with 0 in both
chickmouse_sox9_exp_df <- chickmouse_sox9_exp_df[!apply(chickmouse_sox9_exp_df[, 2:3]==0, 1, all),]


divergent_genes <- chickmouse_sox9_exp_df %>%
  filter((chicken_exp / (mouse_exp + 1) > 2.5) | (mouse_exp / (chicken_exp + 1) > 2.5))
# Add the above information to the dataframe
chickmouse_sox9_exp_df$divergent <- ifelse(chickmouse_sox9_exp_df$gene %in% divergent_genes$gene, "divergent", "non-divergent")

# Adjust the axis limits
x_limit <- max(chickmouse_sox9_exp_df$chicken_exp, na.rm = TRUE)
y_limit <- max(chickmouse_sox9_exp_df$mouse_exp, na.rm = TRUE)
##remove rownames with ENSG
chickmouse_sox9_exp_df <- chickmouse_sox9_exp_df[!grepl("ENSG", chickmouse_sox9_exp_df$gene),]
# Plot the scatter plot
ggplot(chickmouse_sox9_exp_df, aes(x = mouse_exp, y = chicken_exp)) +
  geom_point(aes(color = divergent), alpha = 1) +
  geom_text_repel(data = subset(chickmouse_sox9_exp_df, divergent == "divergent"),
                  aes(label = gene), size = 3, box.padding = 0.3) +
  scale_color_manual(values = c("non-divergent" = "grey", "divergent" = "darkred")) +
  labs(x = "Expression in Mouse (ln(CPM + 1))",
       y = "Expression in Chicken (ln(CPM + 1))",
       subtitle = paste("R =", round(cor(chickmouse_sox9_exp_df$mouse_exp, chickmouse_sox9_exp_df$chicken_exp, method = "pearson"), 2))) +
  xlim(0, x_limit) +
  ylim(0, y_limit) +
  theme_minimal()

goi <-  c("CPED1", "ROBO2", "TFF2", "BEX1", "BEX3", "BEX4","CLDN6","NKX2-1","SOX9")
#SOX9_genes <- c("CPED1", "ROBO2", "TFF2", "BEX1", "BEX3", "BEX4","CLDN6","NKX2-1","SOX9")
chickmouse_sox9_exp_df$GOI <- ifelse(chickmouse_sox9_exp_df$gene %in% goi, "goi", "other")
# Create a new column for highlight- add both goi and divergent info
chickmouse_sox9_exp_df$Highlight <- "other" 
chickmouse_sox9_exp_df$Highlight[chickmouse_sox9_exp_df$gene %in% goi] <- "goi"
chickmouse_sox9_exp_df$Highlight[chickmouse_sox9_exp_df$gene %in% divergent_genes$gene] <- "divergent"
# if you need prioritize labelling of goi over divergent in case of overlap
chickmouse_sox9_exp_df$Highlight[chickmouse_sox9_exp_df$gene %in% goi & chickmouse_sox9_exp_df$gene %in% divergent_genes$gene] <- "goi"
highlight_colors <- c("other" = "grey", "goi" = "black", "divergent" = "darkgreen")

#highlight diveregent and goi
ggplot(chickmouse_sox9_exp_df, aes(x = mouse_exp, y = chicken_exp, color = Highlight)) +
  geom_point(alpha = 1) +
  scale_color_manual(values = highlight_colors) +
  geom_text_repel(data = subset(chickmouse_sox9_exp_df, Highlight == "divergent", color = "darkred"),
                  aes(label = gene), size = 3, box.padding = 0.3) + 
  geom_text_repel(data = subset(chickmouse_sox9_exp_df, Highlight == "goi"),
                  aes(label = gene), size = 3, box.padding = 0.3, color = "darkgreen",  max.overlaps = 1000) +  
  labs(x = "Expression in Mouse (ln(CPM + 1))",
       y = "Expression in Chicken (ln(CPM + 1))",
       subtitle = paste("R =", round(cor(chickmouse_sox9_exp_df$mouse_exp, chickmouse_sox9_exp_df$chicken_exp, method = "pearson"), 2))) +
  xlim(0, x_limit) +
  ylim(0, y_limit) +
  theme_minimal() 

# Identify divergent genes, I set the cutoff to 2.5 here
divergent_genes <- chickmouse_sox9_exp_df %>%
  filter((chicken_exp / (mouse_exp + 1) > 1.5) | (mouse_exp / (chicken_exp + 1) > 2.5))
convergent_genes <- chickmouse_sox9_exp_df %>%
  filter(chicken_exp > 2.5 & mouse_exp > 2.5)
# Add the above information to the dataframe
chickmouse_sox9_exp_df$divergent <- ifelse(chickmouse_sox9_exp_df$gene %in% divergent_genes$gene, "divergent", "non-divergent")
chickmouse_sox9_exp_df$convergent <- ifelse(chickmouse_sox9_exp_df$gene %in% convergent_genes$gene, "convergent", "non-convergent")
# Adjust the axis limits
x_limit <- max(chickmouse_sox9_exp_df$chicken_exp, na.rm = TRUE)
y_limit <- max(chickmouse_sox9_exp_df$mouse_exp, na.rm = TRUE)
##remove rownames with ENSG
chickmouse_sox9_exp_df <- chickmouse_sox9_exp_df[!grepl("ENSG", chickmouse_sox9_exp_df$gene),]

goi <-  c("CPED1", "ROBO2", "TFF2", "BEX1", "BEX3", "BEX4","CLDN6","NKX2-1","SOX9")
#SOX9_genes <- c("CPED1", "ROBO2", "TFF2", "BEX1", "BEX3", "BEX4","CLDN6","NKX2-1","SOX9")
chickmouse_sox9_exp_df$GOI <- ifelse(chickmouse_sox9_exp_df$gene %in% goi, "goi", "other")
# Create a new column for highlight- add both goi and divergent info
chickmouse_sox9_exp_df$Highlight <- "other" 
chickmouse_sox9_exp_df$Highlight[chickmouse_sox9_exp_df$gene %in% goi] <- "goi"
chickmouse_sox9_exp_df$Highlight[chickmouse_sox9_exp_df$gene %in% divergent_genes$gene] <- "divergent"
chickmouse_sox9_exp_df$Highlight[chickmouse_sox9_exp_df$gene %in% convergent_genes$gene] <- "convergent"
# if you need prioritize labelling of goi over divergent in case of overlap
chickmouse_sox9_exp_df$Highlight[chickmouse_sox9_exp_df$gene %in% goi & chickmouse_sox9_exp_df$gene %in% divergent_genes$gene] <- "goi"
highlight_colors <- c("other" = "grey", "goi" = "black", "divergent" = "darkgreen", "convergent" = "grey")

#highlight diveregent and goi
ggplot(chickmouse_sox9_exp_df, aes(x = mouse_exp, y = chicken_exp, color = Highlight)) +
  geom_point(alpha = 1) +
  scale_color_manual(values = highlight_colors) +
  geom_text_repel(data = subset(chickmouse_sox9_exp_df, Highlight == c("divergent", "convergent")),
                  aes(label = gene), size = 3, box.padding = 0.3) + 
  geom_text_repel(data = subset(chickmouse_sox9_exp_df, Highlight == "goi"),
                  aes(label = gene), size = 3, box.padding = 0.3, color = "black", max.overlaps = 1000) +  
  labs(x = "Expression in Mouse (ln(CPM + 1))",
       y = "Expression in Chicken (ln(CPM + 1))",
       subtitle = paste("R =", round(cor(chickmouse_sox9_exp_df$mouse_exp, chickmouse_sox9_exp_df$chicken_exp, method = "pearson"), 2))) +
  xlim(0, x_limit) +
  ylim(0, y_limit) +
  theme_minimal() 

#####FIGURE 5#####
E20chicken.data <- Read10X(data.dir = "/Users/knger/Box/scRNA-seq/E20-Chicken-RNA/outs/galGal6/")

E20 = CreateSeuratObject(counts = E20chicken.data, project = "E20 Chicken", min.cells = 3, min.features = 200) %>% PercentageFeatureSet(pattern = "^mt-", col.name = "percent.mt") %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=F)
DefaultAssay(E20) <- "RNA"
E20 <- RunPCA (E20, npcs=30) %>% RunUMAP (dims=1:30) %>% FindNeighbors() %>% FindClusters()
DimPlot(E20, label=T)

FeaturePlot(E20, features="HBBA", split.by="orig.ident", cols = c("lightgrey", "red"))
VlnPlot(E20, features = c("CDH1","CDH5","PTPRC","COL3A1"))

#subset epi clusters 2,3,5,9,10,14 from reso 1.2
E20 <- subset(E20, idents=c(1,2,5,8,12,14))
E20 <- FindVariableFeatures(E20)
E20 <- ScaleData(E20)
E20 <- RunPCA (E20, npcs=30) %>% RunUMAP (dims=1:30) %>% FindNeighbors() %>% FindClusters()
DimPlot(E20, label=T)
VlnPlot(E20, features = c("CDH1","CDH5","PTPRC","COL3A1"))

#subset out doublet cluster 6
E20 <- subset(E20, idents=c(0:5,7,8,9))
E20 <- FindVariableFeatures(E20)
E20 <- ScaleData(E20)
E20 <- RunPCA (E20, npcs=30) %>% RunUMAP (dims=1:30) %>% FindNeighbors() %>% FindClusters()
#DimPlot(E15, label=F, group.by = "orig.ident")
DimPlot(E20, label=T)
VlnPlot(E20, features = c("CDH1","CDH5","PTPRC","COL3A1"))
VlnPlot(E20, features = c("SOX9","VEGFA","NKX2-1", "LAMP3", "ASCL1", "KRT14"))

#subset out doublet cluster 7
E20 <- subset(E20, idents=c(0:6,8:11))
E20 <- FindVariableFeatures(E20)
E20 <- ScaleData(E20)
E20 <- RunPCA (E20, npcs=30) %>% RunUMAP (dims=1:30) %>% FindNeighbors() %>% FindClusters()
#DimPlot(E15, label=F, group.by = "orig.ident")
DimPlot(E20, label=T)
VlnPlot(E20, features = c("CDH1","CDH5","PTPRC","COL3A1"))
VlnPlot(E20, features = c("SOX9","VEGFA","NKX2-1", "LAMP3", "CCNO", "KRT14"))

#subset out doublet cluster 9+10
E20 <- subset(E20, idents=c(0:8))
E20 <- FindVariableFeatures(E20)
E20 <- ScaleData(E20)
E20 <- RunPCA (E20, npcs=30) %>% RunUMAP (dims=1:30) %>% FindNeighbors() %>% FindClusters()
#DimPlot(E15, label=F, group.by = "orig.ident")
DimPlot(E20, label=T)
VlnPlot(E20, features = c("CDH1","CDH5","PTPRC","COL3A1"))
VlnPlot(E20, features = c("SOX9","VEGFA","PCNA", "LAMP3", "CCNO", "KRT14"))

#name cell types
E20@meta.data$cell_type <- E20@meta.data$seurat_clusters

E20@meta.data$cell_type <- plyr::mapvalues(x = E20@meta.data$cell_type,
                                           from = c("0","1","2","3","4","5","6","7","8"),
                                           to = c("AT1","AT1","KRT14","AT2","AT2","AT2","AT2","AT2","KRT14"))

E20$geno <- substr(colnames(E20@assays$RNA@counts), 18, 18) #cell barcode

E20$geno <- plyr::mapvalues(x = E20$geno,
                            from = "1",
                            to = "E20 Chicken")

Idents(E20) <- "cell_type"
table(E20@active.ident)

#make umap of celltypes for fig 5B
DimPlot_scCustom(E20, reduction="umap", colors_use = c("#dd3497","#ca0020","#66bd63"),label = F)

#make heatmap for fig 5C
top25 <- FindAllMarkers(E20, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.05)
tmp <- group_by(top25, cluster) %>% top_n(n=20, wt=avg_log2FC)
DoHeatmap(E20, features = tmp$gene, group.by = "cell_type", group.colors = c("#dd3497","#ca0020","#66bd63"))+ NoLegend()  #+ scale_fill_gradientn(colors = c("#7b3294", "white", "#008837"))

write.csv(top25, file = "/Users/knger/Box/Chicken paper/Single cell plots/e20-celltype-markers.csv")

#make featureplots for fig 5D
FeaturePlot_scCustom(seurat_object = E20, features = "PCNA", colors_use = viridis_inferno_dark_high)

#make volcano plot for fig 5E
## Krt14 only 
VlnPlot(E13, features = c("SOX9","VEGFA","NKX2-1", "LAMP3", "ASCL1", "KRT14"))
VlnPlot(E20, features = c("SOX9","VEGFA","NKX2-1", "LAMP3", "ASCL1", "KRT14"))

Idents(E20) <- "cell_type"
Idents(E13) <- "cell_type"
table(E13@active.ident)
E13_krt14 <- subset(E13, idents=c("KRT14"))
E13_krt14 <- FindVariableFeatures(E13_krt14)
E13_krt14 <- ScaleData(E13_krt14)
E13_krt14 <- RunPCA (E13_krt14, npcs=30) %>% RunUMAP (dims=1:30) %>% FindNeighbors() %>% FindClusters()
#DimPlot(E15, label=F, group.by = "orig.ident")
DimPlot(E13_krt14, label=T)
#VlnPlot(E20_krt14, features = c("CDH1","CDH5","PTPRC","COL3A1"))
VlnPlot(E13_krt14, features = c("SOX9","VEGFA","NKX2-1", "LAMP3", "ASCL1", "KRT14"))

Idents(E20) <- "cell_type"
table(E20@active.ident)
E20_krt14 <- subset(E20, idents=c("KRT14"))
E20_krt14 <- FindVariableFeatures(E20_krt14)
E20_krt14 <- ScaleData(E20_krt14)
E20_krt14 <- RunPCA (E20_krt14, npcs=30) %>% RunUMAP (dims=1:30) %>% FindNeighbors() %>% FindClusters()
#DimPlot(E15, label=F, group.by = "orig.ident")
DimPlot(E20_krt14, label=T)
#VlnPlot(E20_krt14, features = c("CDH1","CDH5","PTPRC","COL3A1"))
VlnPlot(E20_krt14, features = c("SOX9","VEGFA","NKX2-1", "LAMP3", "ASCL1", "KRT14"))

Idents(E20_krt14) <- "seurat_clusters"
table(E20_krt14@active.ident)
E20_krt14 <- subset(E20_krt14, idents=c(2))
E20_krt14 <- FindVariableFeatures(E20_krt14)
E20_krt14 <- ScaleData(E20_krt14)
E20_krt14 <- RunPCA (E20_krt14, npcs=30) %>% RunUMAP (dims=1:30) %>% FindNeighbors() %>% FindClusters()
#DimPlot(E15, label=F, group.by = "orig.ident")
DimPlot(E20_krt14, label=T)
#VlnPlot(E20_krt14, features = c("CDH1","CDH5","PTPRC","COL3A1"))
VlnPlot(E20_krt14, features = c("SOX9","VEGFA","NKX2-1", "LAMP3", "ASCL1", "KRT14"))


##krt14 vol
krt14.anchors <- FindIntegrationAnchors(object.list = list(E20,E13))
#7878 anchors, 671 retained
integrated_krt14 <- IntegrateData(anchorset = krt14.anchors)

DefaultAssay(integrated_krt14) <- "integrated"
integrated_krt14 <- ScaleData(integrated_krt14, verbose = FALSE)
integrated_krt14 <- RunPCA(integrated_krt14, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
integrated_krt14 <- RunUMAP(integrated_krt14, reduction = "pca", dims = 1:30)
integrated_krt14 <- FindNeighbors(integrated_krt14, reduction = "pca", dims = 1:30)
integrated_krt14 <- FindClusters(integrated_krt14, resolution = 1.5)
# Visualization
p1 <- DimPlot(integrated_krt14, reduction = "umap", group.by = "cell_type")
p2 <- DimPlot(integrated_krt14, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

integrated_krt14$geno.cell_type <- paste(integrated_krt14$geno,integrated_krt14$cell_type,sep = "_")

Idents(integrated_krt14) <- "geno.cell_type"
table(Idents(integrated_krt14))
diffvol <- FindMarkers(integrated_krt14, ident.1 = "E20 Chicken_KRT14", ident.2 = "E13 Chicken_KRT14", min.pct = 0.05, logfc.threshold = 0.05, pseudocount.use = 0.001)
#diffvol <- diffvol[!(rownames(diffvol) %in% boring.genes),]
EnhancedVolcano(diffvol, lab = rownames(diffvol), x='avg_log2FC', y='p_val_adj', FCcutoff = 1, pCutoff = 10e-6, 
                xlim=c(-25,25), 
                gridlines.major = F, gridlines.minor = F, xlab = "ln(fold change)", colAlpha = 1,
                title=NULL, subtitle=NULL, col=c('black','black','black','red'))


write.csv(diffvol, file = "/Users/knger/Box/Chicken Paper/Single cell plots/krt-volcano.csv")

#for fig 5A
#Run 297-299 and then immediately name celltypes
#subset out blood
E20 <- subset(E20, idents=c(1:14,16:18))
E20 <- FindVariableFeatures(E20)
E20 <- ScaleData(E20)
E20 <- RunPCA (E20, npcs=30) %>% RunUMAP (dims=1:30) %>% FindNeighbors() %>% FindClusters()
#DimPlot(E15, label=F, group.by = "orig.ident")
DimPlot(E20, label=T)
VlnPlot(E20, features = c("CDH1","CDH5","PTPRC","COL3A1"))
#subset out low count
E20 <- subset(E20, idents=c(0:8,10:17))
E20 <- FindVariableFeatures(E20)
E20 <- ScaleData(E20)
E20 <- RunPCA (E20, npcs=30) %>% RunUMAP (dims=1:30) %>% FindNeighbors() %>% FindClusters()
#DimPlot(E15, label=F, group.by = "orig.ident")
DimPlot(E20, label=T)
VlnPlot(E20, features = c("CDH1","CDH5","PTPRC","COL3A1"))
#subset out low count
E20 <- subset(E20, idents=c(0:3,5:8,10:14,16))
E20 <- FindVariableFeatures(E20)
E20 <- ScaleData(E20)
E20 <- RunPCA (E20, npcs=30) %>% RunUMAP (dims=1:30) %>% FindNeighbors() %>% FindClusters()
#DimPlot(E15, label=F, group.by = "orig.ident")
DimPlot(E20, label=T)
VlnPlot(E20, features = c("CDH1","CDH5","PTPRC","COL3A1"))
#subset out low count
E20 <- subset(E20, idents=c(0:13))
E20 <- FindVariableFeatures(E20)
E20 <- ScaleData(E20)
E20 <- RunPCA (E20, npcs=30) %>% RunUMAP (dims=1:30) %>% FindNeighbors() %>% FindClusters()
#DimPlot(E15, label=F, group.by = "orig.ident")
DimPlot(E20, label=F)
VlnPlot(E20, features = c("CDH1","CDH5","PTPRC","COL3A1"))

top5 <- FindAllMarkers(E20, min.pct = 0.05, logfc.threshold = 0.05)
tmp <- group_by(top5, cluster) %>% top_n(n=25, wt=avg_log2FC)

E20@meta.data$cell_type <- E20@meta.data$seurat_clusters

E20@meta.data$cell_type <- plyr::mapvalues(x = E20@meta.data$cell_type,
                                           from = c("0","1","2","3","4","5","6","7","8","9","10","11","12","13"),
                                           to = c("Mes","Epi","Epi","Epi","Endo","Mes","Epi","Endo","Mes","Imm","Mes","Endo","Mes","Endo"))

E20$geno <- substr(colnames(E20@assays$RNA@counts), 18, 18) #cell barcode

E20$geno <- plyr::mapvalues(x = E20$geno,
                            from = "1",
                            to = "E20 Chicken")

Idents(E20) <- "cell_type"
table(E20@active.ident)

#endo
Cluster_Highlight_Plot(seurat_object = E20, cluster_name = c("Endo"), highlight_color = "#6e016b",
                       background_color = "lightgray")

#immune
Cluster_Highlight_Plot(seurat_object = E20, cluster_name = c("Imm"), highlight_color = "#cb181d",
                       background_color = "lightgray")

#epi 
Cluster_Highlight_Plot(seurat_object = E20, cluster_name = c("Epi"), highlight_color = "#41ab5d",
                       background_color = "lightgray")

#mes
Cluster_Highlight_Plot(seurat_object = E20, cluster_name = c("Mes"), highlight_color = "#225ea8",
                       background_color = "lightgray")


Stacked_VlnPlot(seurat_object = E20, features = c("CDH1","CDH5","COL3A1","PTPRC"), x_lab_rotate = TRUE,
                colors_use = c("#41ab5d","#225ea8","#cb181d","#6e016b"), split.by = "cell_type")

DimPlot_scCustom(seurat_object = E20, colors_use = c("#225ea8","#41ab5d","#6e016b","#cb181d"), label = F)


saveRDS(E20, file="/Users/knger/Box/scRNA-seq/E20-Chicken-RNA/E20_allcells.rds")

#for fig 5G
#for monocle
library(monocle3)

E13_m <- subset(E13, idents= c("SOX9", "KRT14"))

m.anchors <- FindIntegrationAnchors(object.list = list(E20,E13_m))
#7938 anchors, 1932 retained
integrated_m <- IntegrateData(anchorset = m.anchors)

DefaultAssay(integrated_m) <- "integrated"
integrated_m <- ScaleData(integrated_m, verbose = FALSE)
integrated_m <- RunPCA(integrated_m, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
integrated_m <- RunUMAP(integrated_m, reduction = "pca", dims = 1:30)
integrated_m <- FindNeighbors(integrated_m, reduction = "pca", dims = 1:30)
integrated_m <- FindClusters(integrated_m, resolution = 1.5)
# Visualization
p1 <- DimPlot(integrated_m, reduction = "umap", group.by = "cell_type")
p2 <- DimPlot(integrated_m, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

integrated_m$geno.cell_type <- paste(integrated_m$geno,integrated_m$cell_type,sep = "_")


integrated_m@active.assay = 'RNA'
cds <- as.cell_data_set(integrated_m)
cds <- cluster_cells(cds, resolution=1e-3)
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(integrated_m[["RNA"]])



p1 <- plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = T)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = T)
wrap_plots(p1, p2)

cds <- learn_graph(cds, use_partition = F, verbose = FALSE)
plot_cells(cds,
           color_cells_by = "cell_type",
           label_groups_by_cluster=FALSE,
           label_leaves=F,
           label_branch_points=FALSE,
           label_roots=F)

cds_graph_test_results <- graph_test(cds,
                                     neighbor_graph = "principal_graph",
                                     cores = 8)

rowData(cds)$gene_short_name <- row.names(rowData(cds))

head(cds_graph_test_results, error=FALSE, message=FALSE, warning=FALSE)

deg_ids <- rownames(subset(cds_graph_test_results[order(cds_graph_test_results$morans_I, decreasing = TRUE),], q_value < 0.05))

rownames(cds@preprocess_aux$gene_loadings) 
rownames(cds) 




gene_module_df <- find_gene_modules(cds[deg_ids,], resolution=1e-2)

plot_cells(cds,
           color_cells_by = "geno.cell_type",
           label_cell_groups=FALSE,
           label_leaves=F,
           label_branch_points=F,
           graph_label_size=1.5)

cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5, label_roots=F)

# Add module score chicken E13E20
DefaultAssay(integrated_m) <- "RNA"

AT1 <- subset(tmp, cluster == "AT1")
AT2 <- subset(tmp, cluster == "AT2")
KRT14 <- subset(tmp, cluster == "KRT14")

AT1.list <- list(AT1$gene)
integrated_m <- AddModuleScore(integrated_m, features = AT1.list, name = "AT1")


AT2.list <- list(AT2$gene)
integrated_m <- AddModuleScore(integrated_m, features = AT2.list, name = "AT2")


KRT14.list <- list(KRT14$gene)
integrated_m <- AddModuleScore(integrated_m, features = KRT14.list, name = "KRT14")

FeaturePlot_scCustom(seurat_object = integrated_m, features = "AT11")

FeaturePlot_scCustom(seurat_object = integrated_m, features = "AT21")

FeaturePlot_scCustom(seurat_object = integrated_m, features = "KRT141")

#for Fig S4C
E20_m <- subset(E20, idents= c("AT2", "KRT14"))
DefaultAssay(E20_m) <- "RNA"
VlnPlot(E20_m, features = c("KRT14", "LAMP3", "SFTPC"),group.by = "cell_type", cols=c("#e6f5d0","#fcbba1"), pt.size = 0)
