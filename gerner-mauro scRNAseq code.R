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

#####FIGURE 5#####
#object to make fig 5A
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


#read epithelial only data for the rest of figure 5
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



#make umap for cell types for fig 5B
DimPlot_scCustom(E13, reduction="umap", colors_use = c("#0cb702","#0571b0","#ca0020","#ed68ed"),label = F)

#make heatmap for fig 5C
Idents(E13) <- "cell_type"
top25_13 <- FindAllMarkers(E13, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.05)
tmp_13 <- group_by(top25_13, cluster) %>% top_n(n=20, wt=avg_log2FC)
DoHeatmap(E13, features = tmp_13$gene, group.by = "cell_type", group.colors = c("#0cb702","#0571b0","#ca0020","#ed68ed")) + NoLegend()

#make featureplots for fig 5D
FeaturePlot_scCustom(seurat_object = epi, colors_use = viridis_inferno_dark_high, features = "CLDN10")

FeaturePlot(epi, features="CLDN10", split.by="orig.ident", cols = c("lightgrey", "red"))

#make scatterplot for fig 5E
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
chickmouse.features <- SelectIntegratichionFeatures(object.list = chickmouse.list)
chickmouse.anchors <- FindIntegrationAnchors(object.list = chickmouse.list, anchor.features = chickmouse.features)
#integrate
chickmouse <- IntegrateData(anchorset = chickmouse.anchors)
DefaultAssay(chickmouse) <- "integrated"
chickmouse <- ScaleData(chickmouse) %>% RunPCA(npcs = 30) %>% RunUMAP(dims=1:20) %>% FindNeighbors(dims=1:20) %>% FindClusters(resolution = 0.1)

chickmouse <- readRDS("/Users/rnayak/OneDrive - UTHealth Houston/kamryn/chickmouse.rds")
DefaultAssay(chickmouse) <- "RNA"

#####Aggregate expression data by cell type and genotype #####
#figured out the ln CPM from here  https://github.com/satijalab/seurat/issues/2496
# Convert aggregated expression to a Seurat object #see more here about normalization help("NormalizeData")
avg_chickmouse_expr_sobj <- AggregateExpression(chickmouse, group.by = c("cell_type", "geno"), 
                                                assays = "RNA", normalization.method = "LogNormalize", 
                                                scale.factor = 1e6, return.seurat = TRUE) ##Log Normalize - Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p
##scale factor - scale.factor = 1e6 For counts per million (CPM) 
sox9_log_data <- GetAssayData(avg_chickmouse_expr_sobj)
#### Extract average expression for chicken and mouse SOX9 cells ####
###calculated for all cell types so you can just start from here for other cells 
chicken_sox9_expr <- sox9_log_data[,"SOX9_E13 Chicken"]
mouse_sox9_expr <- sox9_log_data[,"SOX9_E15 Mouse"]

#### Combine into a data frame ####
chickmouse_sox9_exp_df <- data.frame(
  gene = rownames(avg_chickmouse_expr_sobj$RNA),
  chicken_exp = as.numeric(chicken_sox9_expr),
  mouse_exp = as.numeric(mouse_sox9_expr)
)

####remove rows with 0 in both ###
chickmouse_sox9_exp_df <- chickmouse_sox9_exp_df[!apply(chickmouse_sox9_exp_df[, 2:3]==0, 1, all),]
## Remove rows with ENSG gene names
chickmouse_sox9_exp_df <- chickmouse_sox9_exp_df[!grepl("ENSG", chickmouse_sox9_exp_df$gene),]

####Filtering out genes other than shared genes between the chick and mouse and plotting again #####
e15_mouse <- read.csv('E15_mouse.csv', header = F)
e15_chick <- read.csv('E15_chick.csv', header =F)
e15_mouse <- as.character(e15_mouse$V1)
e15_chick <- as.character(e15_chick$V1)
e15_chickmouse_shared <- intersect(e15_mouse, e15_chick) ##9754 genes left

#### Filter the data frame to keep only shared genes ####
chickmouse_sox9_exp_filt_df <- chickmouse_sox9_exp_df %>%
  filter(gene %in% e15_chickmouse_shared)

#### Highlight genes of interest (GOI) ####
goi <- c("CPED1", "ROBO2", "TFF2", "BEX1", "BEX3", "BEX4", "CLDN6", "NKX2-1", "SOX9", "SFTPC", "RPS17", "SLIT3",
         "SELENOW", "MIF", "NME1", "DYNLL1")
chickmouse_sox9_exp_filt_df$GOI <- ifelse(chickmouse_sox9_exp_filt_df$gene %in% goi, "goi", "other")

#### Update the axis limits based on filtered data ####
x_limit <- max(chickmouse_sox9_exp_filt_df$chicken_exp, na.rm = TRUE)
y_limit <- max(chickmouse_sox9_exp_filt_df$mouse_exp, na.rm = TRUE)
highlight_colors <- c("other" = "grey", "goi" = "black")

#### Plot the data ####
ggplot(chickmouse_sox9_exp_filt_df, aes(x = mouse_exp, y = chicken_exp, color = GOI)) +
  geom_point(alpha = 1) +
  scale_color_manual(values = highlight_colors) +
  geom_text_repel(data = subset(chickmouse_sox9_exp_filt_df, GOI == "goi"),
                  aes(label = gene), size = 3, box.padding = 0.3, color = "black", max.overlaps = 10000) +  
  labs(x = "Expression in Mouse (ln(CPM))",
       y = "Expression in Chicken (ln(CPM))",
       subtitle = paste("R =", round(cor(chickmouse_sox9_exp_filt_df$mouse_exp, chickmouse_sox9_exp_filt_df$chicken_exp, method = "pearson"), 2))) +
  xlim(0, x_limit) +
  ylim(0, 3.5)

#####FIGURE 6#####
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

#make umap of celltypes for fig 6B
DimPlot_scCustom(E20, reduction="umap", colors_use = c("#dd3497","#ca0020","#66bd63"),label = F)

#make heatmap for fig 6C
top25_20<- FindAllMarkers(E20, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.05)
tmp_20 <- group_by(top25, cluster) %>% top_n(n=20, wt=avg_log2FC)
DoHeatmap(E20, features = tmp$gene, group.by = "cell_type", group.colors = c("#dd3497","#ca0020","#66bd63"))+ NoLegend()  #+ scale_fill_gradientn(colors = c("#7b3294", "white", "#008837"))

write.csv(top25, file = "/Users/knger/Box/Chicken paper/Single cell plots/e20-celltype-markers.csv")

#make featureplots for fig 6D
FeaturePlot_scCustom(seurat_object = E20, features = "PCNA", colors_use = viridis_inferno_dark_high)

#make volcano plot for fig 6E
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
DimPlot(E20_krt14, label=T)
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
EnhancedVolcano(diffvol, lab = rownames(diffvol), x='avg_log2FC', y='p_val_adj', FCcutoff = 1, pCutoff = 10e-6, 
                xlim=c(-25,25), 
                gridlines.major = F, gridlines.minor = F, xlab = "ln(fold change)", colAlpha = 1,
                title=NULL, subtitle=NULL, col=c('black','black','black','red'))


write.csv(diffvol, file = "/Users/knger/Box/Chicken Paper/Single cell plots/krt-volcano.csv")

#for fig 6A
#Run 297-299 and then immediately name celltypes
#subset out blood
E20 <- subset(E20, idents=c(1:14,16:18))
E20 <- FindVariableFeatures(E20)
E20 <- ScaleData(E20)
E20 <- RunPCA (E20, npcs=30) %>% RunUMAP (dims=1:30) %>% FindNeighbors() %>% FindClusters()
#DimPlot(E15, label=F, group.by = "orig.ident")
DimPlot(E20, label=T)
VlnPlot(E20, features = c("CDH1","CDH5","PTPRC","COL3A1","HPSE"))
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


Idents(E20) <- "seurat_clusters"
E20@meta.data$cell_type <- E20@meta.data$seurat_clusters

E20@meta.data$cell_type <- plyr::mapvalues(x = E20@meta.data$cell_type,
                                           from = c("0","1","2","3","4","5","6","7","8","9","10","11","12","13"),
                                           to = c("Blood","Epi","Epi","Epi","Endo","Mes","Epi","Endo","Blood","Imm","Mes","Endo","Mes","Endo"))

E20$geno <- substr(colnames(E20@assays$RNA@counts), 18, 18) #cell barcode

E20$geno <- plyr::mapvalues(x = E20$geno,
                            from = "1",
                            to = "E20 Chicken")

Idents(E20) <- "cell_type"
table(E20@active.ident)


Stacked_VlnPlot(seurat_object = E20, features = c("HBBA","CDH1","CDH5","COL3A1","PTPRC"), x_lab_rotate = TRUE,
                colors_use = c("#67001f","#41ab5d","#6e016b","#225ea8","#cb181d"), split.by = "cell_type")

DimPlot_scCustom(seurat_object = E20, colors_use = c("#67001f","#41ab5d","#6e016b","#225ea8","#cb181d"), label = F)



saveRDS(E20, file="/Users/knger/Box/scRNA-seq/E20-Chicken-RNA/E20_allcells.rds")

#for fig 6H
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

# Add module score chicken E13E20 for Fig 6H
DefaultAssay(integrated_m) <- "RNA"

AT1 <- subset(tmp_20, cluster == "AT1")
AT2 <- subset(tmp_20, cluster == "AT2")
KRT14 <- subset(tmp_20, cluster == "KRT14")
AT0 <- subset(tmp, cluster == "AT0")

AT1.list <- list(AT1$gene)
integrated_m <- AddModuleScore(integrated_m, features = AT1.list, name = "AT1")


AT2.list <- list(AT2$gene)
integrated_m <- AddModuleScore(integrated_m, features = AT2.list, name = "AT2")


KRT14.list <- list(KRT14$gene)
integrated_m <- AddModuleScore(integrated_m, features = KRT14.list, name = "KRT14")

AT0.list <- list(AT0$gene)
integrated_m <- AddModuleScore(integrated_m, features = AT0.list, name = "AT0")

FeaturePlot_scCustom(seurat_object = integrated_m, features = "AT11")

FeaturePlot_scCustom(seurat_object = integrated_m, features = "AT21")

FeaturePlot_scCustom(seurat_object = integrated_m, features = "KRT141")

FeaturePlot_scCustom(seurat_object = integrated_m, features = "AT01")

#for Fig 6E
###AT1 AT2 chicken vs mouse####
e19e20 <- readRDS("/Users/rnayak/OneDrive - UTHealth Houston/kamryn/e19e20.rds")
DefaultAssay(e19e20) <- "RNA"

# Convert aggregated expression to a Seurat object
avg_e19e20_expr_sobj <- AggregateExpression(e19e20, group.by = c("cell_type", "geno"), 
                                            assays = "RNA", normalization.method = "LogNormalize", 
                                            scale.factor = 1e6, return.seurat = TRUE)

e19e20_log_data <- GetAssayData(avg_e19e20_expr_sobj)

# Extract average expression for chicken and mouse AT1 and AT2 cells 
chicken_AT1_expr <- e19e20_log_data[,"AT1_E20 Chicken"]
mouse_AT1_expr <- e19e20_log_data[,"AT1_E19 Mouse"]

chicken_AT2_expr <- e19e20_log_data[,"AT2_E20 Chicken"]
mouse_AT2_expr <- e19e20_log_data[,"AT2_E19 Mouse"]

# Combine into a data frame 
e19e20_AT1_exp_df <- data.frame(
  gene = rownames(avg_e19e20_expr_sobj$RNA),
  chicken_exp = as.numeric(chicken_AT1_expr),
  mouse_exp = as.numeric(mouse_AT1_expr)
)

e19e20_AT2_exp_df <- data.frame(
  gene = rownames(avg_e19e20_expr_sobj$RNA),
  chicken_exp = as.numeric(chicken_AT2_expr),
  mouse_exp = as.numeric(mouse_AT2_expr)
)
#remove rows with 0 in both
e19e20_AT1_exp_df <- e19e20_AT1_exp_df[!apply(e19e20_AT1_exp_df[, 2:3]==0, 1, all),]
e19e20_AT2_exp_df <- e19e20_AT2_exp_df[!apply(e19e20_AT2_exp_df[, 2:3]==0, 1, all),]
##remove rownames with ENSG
e19e20_AT1_exp_df <- e19e20_AT1_exp_df[!grepl("ENSG", e19e20_AT1_exp_df$gene),]
e19e20_AT2_exp_df <- e19e20_AT2_exp_df[!grepl("ENSG", e19e20_AT2_exp_df$gene),]


####Filtering out genes other than shared genes between the chick and mouse and plotting again #####
e19_mouse <- read.csv('E19_mouse.csv', header = F)
e19_chick <- read.csv('E19_chick.csv', header =F)
e19_chickmouse_shared <- intersect(e19_mouse, e19_chick) ##9548 genes left
colnames(e19_chickmouse_shared) <- 'gene'
e19_chickmouse_shared <- as.character(e19_chickmouse_shared$gene)

#### Filter the data frame to keep only shared genes ####
e19e20_AT1_exp_filt_df <- e19e20_AT1_exp_df %>%
  filter(gene %in% e19_chickmouse_shared)

e19e20_AT2_exp_filt_df <- e19e20_AT2_exp_df %>%
  filter(gene %in% e19_chickmouse_shared)


###AT1 GOI
goi <- c("KRT19", "KRT8", "AGER", "AQP5", "CAV1", "LY6E",
         "B3GAT2","CYTB", "NKX2-1", "VEGFA", "CDH1", "ACTB", "RPS28")
e19e20_AT1_exp_filt_df$GOI <- ifelse(e19e20_AT1_exp_filt_df$gene %in% goi, "goi", "other")

goi <- c("LYG2", "SFTPA2", "LYZ1", "LYZ2", "BEX2", "BEX4", "NKX2-1", "CDH1", "SFTPC", "SFTPA1", "RPL29", 
         "SLC34A2", "FABP5", "DLK1", "RPL39L", "FABP3", "HSPA2")
e19e20_AT2_exp_filt_df$GOI <- ifelse(e19e20_AT2_exp_filt_df$gene %in% goi, "goi", "other")

x_limit <- max(e19e20_AT2_exp_filt_df$chicken_exp, na.rm = TRUE)
y_limit <- max(e19e20_AT2_exp_filt_df$mouse_exp, na.rm = TRUE)

#### Plot the data ####
ggplot(e19e20_AT1_exp_filt_df, aes(x = mouse_exp, y = chicken_exp, color = GOI)) +
  geom_point(alpha = 1) +
  scale_color_manual(values = highlight_colors) +
  geom_text_repel(data = subset(e19e20_AT1_exp_filt_df, GOI == "goi"),
                  aes(label = gene), size = 3, box.padding = 0.3, color = "black", max.overlaps = 10000) +  
  labs(x = "Expression in Mouse (ln(CPM))",
       y = "Expression in Chicken (ln(CPM))",
       subtitle = paste("R =", round(cor(e19e20_AT1_exp_filt_df$mouse_exp, e19e20_AT1_exp_filt_df$chicken_exp, method = "pearson"), 2))) +
  xlim(0, x_limit) +
  ylim(0, y_limit)

ggplot(e19e20_AT2_exp_filt_df, aes(x = mouse_exp, y = chicken_exp, color = GOI)) +
  geom_point(alpha = 1) +
  scale_color_manual(values = highlight_colors) +
  geom_text_repel(data = subset(e19e20_AT2_exp_filt_df, GOI == "goi"),
                  aes(label = gene), size = 3, box.padding = 0.3, color = "black", max.overlaps = 10000) +  
  labs(x = "Expression in Mouse (ln(CPM))",
       y = "Expression in Chicken (ln(CPM))",
       subtitle = paste("R =", round(cor(e19e20_AT2_exp_filt_df$mouse_exp, e19e20_AT2_exp_filt_df$chicken_exp, method = "pearson"), 2))) +
  xlim(0, x_limit) +
  ylim(0, y_limit)

####divergent genes####
SOX9_divergent_genes <- chickmouse_sox9_exp_filt_df %>%
  filter((chicken_exp / (mouse_exp + 1) > 2) | (mouse_exp / (chicken_exp + 1) > 2))
AT2_divergent_genes <- e19e20_AT2_exp_filt_df %>%
  filter((chicken_exp / (mouse_exp + 1) > 2) | (mouse_exp / (chicken_exp + 1) > 2))


#for Fig S5B
E20_cd45.data <- Read10X(data.dir = "/Users/knger/Downloads/E20-chicken-CD5neg_analysis/outs/galGal6/filtered_feature_bc_matrix")
E20_cd45 = CreateSeuratObject(counts = E20_cd45.data, project = "E20 Chicken", min.cells = 3, min.features = 200) %>% PercentageFeatureSet(pattern = "^mt-", col.name = "percent.mt") %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=F)
DefaultAssay(E20_cd45) <- "RNA"
E20_cd45 <- RunPCA (E20_cd45, npcs=30) %>% RunUMAP (dims=1:30) %>% FindNeighbors() %>% FindClusters()
DimPlot(E20_cd45, label=T)

FeaturePlot(E20_cd45, features=c("CDH1","CDH5","PTPRC","COL3A1","HBBA"), cols = c("lightgrey", "red"))
VlnPlot(E20_cd45, features = c("CDH1","CDH5","PTPRC","COL3A1","HBBA","HPSE"))

#subset dbl clusters 4,11,13
E20_cd45 <- subset(E20_cd45, idents=c(0:3,5:10,12))
E20_cd45 <- FindVariableFeatures(E20_cd45)
E20_cd45 <- ScaleData(E20_cd45)
E20_cd45 <- RunPCA (E20_cd45, npcs=30) %>% RunUMAP (dims=1:30) %>% FindNeighbors() %>% FindClusters()
DimPlot(E20_cd45, label=T)
FeaturePlot(E20_cd45, features=c("CDH1","CDH5","PTPRC","COL3A1","HBBA"), cols = c("lightgrey", "red"))


#name the cell types
E20_cd45@meta.data$cell_type <- E20_cd45@meta.data$seurat_clusters

E20_cd45@meta.data$cell_type <- plyr::mapvalues(x = E20_cd45@meta.data$cell_type,
                                                from = c("0","1","2","3","4","5","6","7","8","9","10","11","12"),
                                                to = c("blood","blood","endo","blood","mes","mes","mes","mes","epi","endo","endo","blood","endo"))

E20_cd45$geno <- substr(colnames(E20_cd45@assays$RNA@counts), 18, 18) #cell barcode

E20_cd45$geno <- plyr::mapvalues(x = E20_cd45$geno,
                                 from = "1",
                                 to = "E20 Chicken")

Idents(E20_cd45) <- "cell_type"
table(E20_cd45@active.ident)

DimPlot(E20_cd45, label=T)
DimPlot_scCustom(E20_cd45, reduction="umap", colors_use = c("#67001f","#6e016b","#225ea8","#41ab5d"),label = F)

Stacked_VlnPlot(seurat_object = E20_cd45, features = c("HBBA","CDH1","CDH5","COL3A1","PTPRC"), x_lab_rotate = TRUE,
                colors_use = c("#67001f","#6e016b","#225ea8","#41ab5d"), split.by = "cell_type")
#for Fig S5D
E20_m <- subset(E20, idents= c("AT2", "KRT14"))
DefaultAssay(E20_m) <- "RNA"
VlnPlot(E20_m, features = c("KRT14", "LAMP3", "SFTPC"),group.by = "cell_type", cols=c("#e6f5d0","#fcbba1"), pt.size = 0)

#for Fig S5E
mouse.anchors <- FindIntegrationAnchors(object.list = list(E15,E19), dims = 1:20)
#9560 anchors, 1822 retained
mouse_integrated <- IntegrateData(anchorset = mouse.anchors)

DefaultAssay(mouse_integrated) <- "integrated"
mouse_integrated <- ScaleData(mouse_integrated, verbose = FALSE)
mouse_integrated <- RunPCA(mouse_integrated, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
mouse_integrated <- RunUMAP(mouse_integrated, reduction = "pca", dims = 1:30)
mouse_integrated <- FindNeighbors(mouse_integrated, reduction = "pca", dims = 1:30)
mouse_integrated <- FindClusters(mouse_integrated)
# Visualization
p1 <- DimPlot(mouse_integrated, reduction = "umap", group.by = "cell_type")
p2 <- DimPlot(mouse_integrated, reduction = "umap", split.by = "geno")
plot_grid(p1, p2)


DimPlot(mouse_integrated, reduction = "umap", group.by = "cell_type", split.by="geno")
DimPlot(mouse_integrated, reduction = "umap",group.by = "cell_type", cols=c("#dd3497","#4dac26","#3288bd","#ed68ed","#d73027","#addd8e"))

mouse_integrated$geno.cell_type <- paste(mouse_integrated$geno,mouse_integrated$cell_type,sep = "_")
Idents(mouse_integrated) <- "geno.cell_type"
table(mouse_integrated@active.ident)



# Add module score for revision
tmp_13<-read.csv(file="/Users/knger/Box/Chicken Paper/CellRef_signature_Epi_2024-11-20.csv")
tmp_AT0 <-read.csv(file="/Users/knger/Box/Chicken Paper/HLCA_AT0_signature_top50_2024-11-20.csv")

DefaultAssay(E20) <- "RNA"
head(tmp_13)
head(tmp_AT0)

table(tmp_13$ï..celltype)

AT1 <- subset(tmp_13, ï..celltype == "AT1")
AT2 <- subset(tmp_13, ï..celltype == "AT2")
Basal <- subset(tmp_13, ï..celltype == "Basal")
Ciliated <- subset(tmp_13, ï..celltype == "Ciliated")
Deuterosomal <- subset(tmp_13, ï..celltype == "Deuterosomal")
Ionocyte <- subset(tmp_13, ï..celltype == "Ionocyte")
MEC <- subset(tmp_13, ï..celltype == "MEC")
Mucous <- subset(tmp_13, ï..celltype == "Mucous")
PNEC <- subset(tmp_13, ï..celltype == "PNEC")
RAS <- subset(tmp_13, ï..celltype == "RAS")
Secretory <- subset(tmp_13, ï..celltype == "Secretory")
SMG_Basal <- subset(tmp_13, ï..celltype == "SMG_Basal")
Tuft <- subset(tmp_13, ï..celltype == "Tuft")
Goblet <- subset(tmp_13, ï..celltype == "Goblet")
AT0 <- subset(tmp_AT0, ï..celltype == "AT0")
Serous <- subset(tmp_13, ï..celltype == "Serous")


AT1.list <- list(AT1$gene)
E20 <- AddModuleScore(E20, features = AT1.list, name = "AT1")

AT2.list <- list(AT2$gene)
E20 <- AddModuleScore(E20, features = AT2.list, name = "AT2")

Basal.list <- list(Basal$gene)
E20 <- AddModuleScore(E20, features = Basal.list, name = "Basal")

Ciliated.list <- list(Ciliated$gene)
E20 <- AddModuleScore(E20, features = Ciliated.list, name = "Ciliated")

Deuterosomal.list <- list(Deuterosomal$gene)
E20 <- AddModuleScore(E20, features = Deuterosomal.list, name = "Deuterosomal")

Ionocyte.list <- list(Ionocyte$gene)
E20 <- AddModuleScore(E20, features = Ionocyte.list, name = "Ionocyte")

MEC.list <- list(MEC$gene)
E20 <- AddModuleScore(E20, features = MEC.list, name = "MEC")

PNEC.list <- list(PNEC$gene)
E20 <- AddModuleScore(E20, features = PNEC.list, name = "PNEC")

RAS.list <- list(RAS$gene)
E20 <- AddModuleScore(E20, features = RAS.list, name = "RAS")

Secretory.list <- list(Secretory$gene)
E20 <- AddModuleScore(E20, features = Secretory.list, name = "Secretory")

SMG_Basal.list <- list(SMG_Basal$gene)
E20 <- AddModuleScore(E20, features = SMG_Basal.list, name = "SMG_Basal")

Tuft <- list(Tuft$gene)
E20 <- AddModuleScore(E20, features = Tuft.list, name = "Tuft")

Serous.list <- list(Serous$gene)
E20 <- AddModuleScore(E20, features = Serous.list, name = "Serous")

AT0.list <- list(AT0$gene)
E20 <- AddModuleScore(E20, features = AT0.list, name = "AT0")

FeaturePlot_scCustom(seurat_object = E20, features = "AT11")

FeaturePlot_scCustom(seurat_object = E20, features = "AT21")

FeaturePlot_scCustom(seurat_object = E20, features = "Basal1")

FeaturePlot_scCustom(seurat_object = E20, features = "Ciliated1")

FeaturePlot_scCustom(seurat_object = E20, features = "Deuterosomal1")

FeaturePlot_scCustom(seurat_object = E20, features = "Goblet1")

FeaturePlot_scCustom(seurat_object = E20, features = "Ionocyte1")

FeaturePlot_scCustom(seurat_object = E20, features = "Suprabasal1")

FeaturePlot_scCustom(seurat_object = E20, features = "MEC1")

FeaturePlot_scCustom(seurat_object = E20, features = "Mucous1")

FeaturePlot_scCustom(seurat_object = E20, features = "PNEC1")

FeaturePlot_scCustom(seurat_object = E20, features = "RAS1")

FeaturePlot_scCustom(seurat_object = E20, features = "Secretory1")

FeaturePlot_scCustom(seurat_object = E20, features = "Serous1")

FeaturePlot_scCustom(seurat_object = E20, features = "SMG_Basal1")

FeaturePlot_scCustom(seurat_object = E20, features = "Tuft1")

FeaturePlot_scCustom(seurat_object = E20, features = "AT01")

