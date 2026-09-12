library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
library(data.table)
library(ggrepel)
library(ggpubr)
library(clustree)
library(celda)
library(SingleCellExperiment)
library(DoubletFinder)
theme_set(theme_cowplot())
setwd("C:/Users/rohin/Desktop/Pre_doc_projects/Project_RB/Huang_Lab")

dge_merged <- readRDS("dge_mergedall_090726.rds")
ncol(dge_merged)
colnames(dge_merged@meta.data)
dge_merged <- JoinLayers(dge_merged)


##Ambient RNA decontamination


sce <- SingleCellExperiment(list(counts = LayerData(dge_merged, assay = "RNA", layer = "counts" )))
sce <- decontX(sce)
sce_count <- decontXcounts(sce)

dge <- CreateSeuratObject(counts = sce_count, project = "Bladder_Final", 
                          min.cells = 0, min.features = 0,
                          meta.data = dge_merged@meta.data)
dim(dge)
colnames(dge@meta.data)
saveRDS( dge, "dge_postDecontX10thSept.rds")


##Assign patient IDs (explicit lookup)

sample_to_patient <- c(
  "GSM8878307_HY_20847_N2_S16" = "20847_N",
  "GSM8878308_HY_20847_T1-1_S19" = "20847_T",
  "GSM8878309_HY_20847_T1-2_S24" = "20847_T",
  "GSM8878310_HY_20847_T2-1_S20" = "20847_T",
  "GSM8878311_HY_20847_T2-2_S25" = "20847_T",
  "GSM8878312_HY_20847_T3-1_S21" = "20847_T",
  "GSM8878313_HY_20847_T3-2_S26" = "20847_T",
  "GSM8878314_HY_20847x_N1-2_S27" = "20847_N",
  "GSM8878315_HY_21032_N1_S10" = "21032_N",
  "GSM8878316_HY_21032_N2_S11" = "21032_N",
  "GSM8878317_HY_21032_N3_S12" = "21032_N",
  "GSM8878318_HY_21032_T1_S9" = "21032_T",
  "GSM8878319_HY_21217_N_48" = "21217_N",
  "GSM8878320_HY_21217_T1_37" = "21217_T",
  "GSM8878321_HY_21217_T1_48" = "21217_T",
  "GSM8878322_HY_21217_T2_48" = "21217_T",
  "GSM8878323_HY_21217_T3_48" = "21217_T",
  "GSM8878324_HY_21222_N_37" = "21222_N",
  "GSM8878325_HY_21222_N_48" = "21222_N",
  "GSM8878326_HY_21222_T1_37" = "21222_T",
  "GSM8878327_HY_21222_T1_48" = "21222_T",
  "GSM8878328_HY_21222_T2_37" = "21222_T",
  "GSM8878329_HY_21222_T2_48" = "21222_T",
  "GSM8878330_HY_21222_T3_48" = "21222_T",
  "GSM8878331_HY_21226_N_37" = "21226_N",
  "GSM8878332_HY_21226_N_48" = "21226_N",
  "GSM8878333_HY_21226_T1_48" = "21226_T",
  "GSM8878334_HY_21226_T2_37" = "21226_T",
  "GSM8878335_HY_21226_T2_48" = "21226_T",
  "GSM8878336_HY_21226_T3_48" = "21226_T",
  "GSM8878337_HY_21262_1_S39" = "21262_T",
  "GSM8878338_HY_21262_2_S40" = "21262_T",
  "GSM8878339_HY_FG_CIS_37" = "FG_CIS",
  "GSM8878340_HY_FG_CIS_48-1" = "FG_CIS",
  "GSM8878341_HY_FG_CIS_48-2" = "FG_CIS",
  "GSM8878342_HY_FG_N_48" = "FG_N",
  "GSM8878343_HY_FG_T1_37" = "FG_T",
  "GSM8878344_HY_FG_T1_48-1" = "FG_T",
  "GSM8878345_HY_FG_T1_48-2" = "FG_T",
  "GSM8878346_HY_FG_T2_37" = "FG_T",
  "GSM8878347_HY_FG_T2_48" = "FG_T",
  "GSM8878348_HY_HG_T1_37" = "HG_T",
  "GSM8878349_HY_HG_T1_48" = "HG_T",
  "GSM8878350_HY_HG_T2_37" = "HG_T",
  "GSM8878351_HY_JM_T1_S13" = "JM_T",
  "GSM8878352_HY_JM_T2_S14" = "JM_T",
  "GSM8878353_HY_JMx_N1_S15" = "JMx_N",
  "GSM8878354_HY_PG_MN1_S7" = "PG_MN",
  "GSM8878355_HY_PG_MN2_S8" = "PG_MN",
  "GSM8878356_HY_PG_MT1_S3" = "PG_MT",
  "GSM8878357_HY_PG_MT2_S4" = "PG_MT",
  "GSM8878358_HY_PG_SN1_S5" = "PG_SN",
  "GSM8878359_HY_PG_SN2_S6" = "PG_SN",
  "GSM8878360_HY_PG_ST1_S1" = "PG_ST",
  "GSM8878361_HY_PG_ST2_S2" = "PG_ST",
  "GSM8878362_HY_PS_T1-1_S17" = "PS_T",
  "GSM8878363_HY_PS_T1-2_S22" = "PS_T",
  "GSM8878364_HY_PS_T2-1_S18" = "PS_T",
  "GSM8878365_HY_PS_T2-2_S23" = "PS_T",
  "GSM8878366_PA_U67-11734-2_Pool_1_2_3_S21_L001" = "11734",
  "GSM8878367_PA_U67-12041-2_Pool_1_2_3_S22_L001" = "12041",
  "GSM8878368_PA_U67-12049_Pool_1_2_3_S20_L001" = "12049",
  "GSM8878369_PA_U67_12050N_Pool_1_2_3_S45_L002" = "12050_N",
  "GSM8878370_PA_U67-12050T1_Pool_1_2_3_S72_L003" = "12050_T",
  "GSM8878371_PA_U67_12050T2_Pool_1_2_3_S103_L004" = "12050_T",
  "GSM8878372_SG-1_S29" = "SG",
  "GSM8878373_SG-2_S30" = "SG"
)

#identical(colnames(sce_count), colnames(dge_merged))
ncol(sce_count) == ncol(dge_merged)

dge_merged$patient <- unname(sample_to_patient[dge_merged$sample])
table(dge_merged$patient, useNA = "always")


##Assign patient-level Name groups(paper-wise mapping)


dge_merged$Name <- "Unknown"
dge_merged$Name[dge_merged$patient %in% c("21217_N", "21217_T")] <- "VAR11"
dge_merged$Name[dge_merged$patient %in% c("21222_N", "21222_T")] <- "VAR1"
dge_merged$Name[dge_merged$patient %in% c("21226_N", "21226_T")] <- "VAR7"
dge_merged$Name[dge_merged$patient %in% c("21262_T")] <- "VAR9"
dge_merged$Name[dge_merged$patient %in% c("FG_CIS", "FG_N", "FG_T")] <- "VAR5"
dge_merged$Name[dge_merged$patient %in% c("HG_T")] <- "UC3"
dge_merged$Name[dge_merged$patient %in% c("20847_T")] <- "VAR6"
dge_merged$Name[dge_merged$patient %in% c("JM_T", "JMx_N")] <- "VAR3"
dge_merged$Name[dge_merged$patient %in% c("PG_SN", "PG_MT", "PG_ST")] <- "PG"
dge_merged$Name[dge_merged$patient %in% c("PS_T")] <- "VAR2"
dge_merged$Name[dge_merged$patient %in% c("11734")] <- "VAR8"
dge_merged$Name[dge_merged$patient %in% c("12041")] <- "VAR4"
dge_merged$Name[dge_merged$patient %in% c("12050_N", "12050_T")] <- "VAR10"
dge_merged$Name[dge_merged$patient %in% c("12049")] <- "UC1"
dge_merged$Name[dge_merged$patient %in% c("SG")] <- "UC2"

table(dge_merged$Name, useNA = "always")
ncol(dge_merged)


##Normalization -LogNormalize, scale factor 10000 (paper's value)


dge <- NormalizeData(object = dge, normalization.method = "LogNormalize" , scale.factor = 10000 )


##Variable features- top 2000(paper's value)


dge <- FindVariableFeatures( object = dge, selection.method = "vst", nfeatures = 2000)
top10 <- head(x= VariableFeatures(object = dge), 10)
top10

plot1 <- VariableFeaturePlot(object = dge)
plot2 <- LabelPoints( plot = plot1 , points = top10 ,repel =TRUE, xnudge = 0, ynudge = 0 )
combined_plot <- plot1 + plot2
ggsave(file = "VariableFeature_2000_10thSept.pdf", plot = combined_plot, width = 20 , height = 20 , units = "cm")


## Scale + PCA - 100 PCs stored(paper's value)

dge <- ScaleData(object = dge, features = VariableFeatures(object = dge))
dge <- RunPCA(object = dge, features = VariableFeatures(object = dge), npcs = 100)

slot(dge[["pca"]],"misc")
print(x= dge[["pca"]], dims = 1:5, nfeatures= 5)

mat <- Seurat::GetAssayData(dge, assay = "RNA", layer = "scale.data")
pca <- dge[["pca"]]
total_variance <- sum(matrixStats::rowVars(mat))
eigValues <- (pca@stdev)^2
varExplained <- eigValues / total_variance
sum(varExplained)


## Attach patient/Name metadata onto post-decontX dge


identical(colnames(dge),colnames(dge_merged))
dge$patient <- dge_merged$patient
dge$Name <- dge_merged$Name
colnames(dge@meta.data)
table(dge$Name, useNA = "always")
p <- DimPlot( object = dge, reduction = "pca", group.by = "Name")
ggsave(file = "PCA_byName10thSept.pdf", plot = p , width = 20 , height = 20 , units = "cm")
png("PCA_heatmap10thSept.png", width = 1200, height = 1000)
DimHeatmap(object = dge, dims = 1:10 , cells = 500, balanced = TRUE)
dev.off()


png("Elbowplot_dge10thSept.png")
ElbowPlot(dge, ndims = 100)
dev.off()


##Clustering - 75 PCs used for graph/clustering/UMAP(paper's value)

n_pc <- 75
dge <- FindNeighbors(dge, dims = 1:n_pc, k.param = 30)
dge <- FindClusters(dge, resolution=c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0))

pdf("Stability10thSept.pdf", width = 15, height = 25)
clustree(dge, prefix = "RNA_snn_res.")
dev.off()

## Paper's final resolution = 0.5 
dge <- FindClusters(dge, resolution = 0.5)
length(unique(Idents(dge)))


## t-SNE(using same 75 PCs)


dge <- RunTSNE(dge, dims = 1:n_pc, perplexity = 40, seed.use = 10, check_duplicates = FALSE)
p <- DimPlot(dge, reduction = "tsne", label = TRUE, pt.size = 1) + NoLegend()
ggsave( file = "TSNE_raw10thSept.pdf" , plot = p , width = 20 , height = 20, units = "cm")


## UMAP(using same 75 PCs)


dge <- RunUMAP (dge, dims = 1:n_pc)

p <- DimPlot (dge, reduction= "umap", label = TRUE) + NoLegend()
ggsave(file= "Umap_raw10thSept.pdf", plot = p, width = 20 , height = 20 , units = "cm")
p <- DimPlot (dge, reduction = "umap", group.by = "Name", label = FALSE)
ggsave(file = "Umap_byPatient10thSept.pdf", plot = p , width = 25 , height = 20 , units = "cm")


## Marker genes for cell-type annotation


marker_genes <- list(
  Epithelial    = c("EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"),
  Immune        = c("PTPRC", "CD3E", "CD68", "CD79A", "LYZ"),
  Fibroblast    = c("COL1A1", "COL1A2", "COL3A1", "DCN", "PDGFRB"),
  Endothelial   = c("PECAM1", "VWF", "CDH5", "CLDN5"),
  Smooth_Muscle = c("ACTA2", "MYH11", "TAGLN", "DES")
)

all_markers <- unlist(marker_genes)
present_markers <- all_markers[all_markers %in% rownames(dge)]
missing_markers <- setdiff(all_markers, present_markers)
missing_markers

DotPlot(dge, features = present_markers, group.by = "seurat_clusters") + 
  RotatedAxis() +
  theme(axis.text.x = element_text(size = 8))
ggsave("Marker_DotPlot.pdf", width = 16 , height = 10)

FeaturePlot( dge, features =c("EPCAM", "PTPRC", "COL1A1", "PECAM1", "ACTA2"),
             reduction = "umap", ncol = 3)
ggsave("Marker_Feature10thSept.pdf", width = 15 , height = 10)
avg_exp <- AverageExpression(dge, features = present_markers, group.by = "seurat_clusters", assays = "RNA")
avg_exp_mat <- avg_exp$RNA
View(avg_exp_mat)


## Assign cell-type labels per cluster(filled on the basis of DotPlot inspection)

cluster_to_celltype <- c(
  "0"  = "Immune",   "1"  = "Fibroblast",    "2"  = "Epithelial",
  "3"  = "Smooth_Muscle",    "4"  = "Fibroblast",     "5"  = "Immune",
  "6"  = "Immune",       "7"  = "Epithelial",     "8"  = "Fibroblast",
  "9"  = "Immune",       "10" = "Epithelial",     "11" = "Epithelial",
  "12" = "Epithelial",   "13" = "Fibroblast",  "14" = "Fibroblast",
  "15" = "Immune",       "16" = "Epithelial",     "17" = "Epithelial",
  "18" = "Epithelial",   "19" = "Endothelial",    "20" = "Immune",
  "21" = "Epithelial",   "22" = "Endothelial",    "23" = "Epithelial",
  "24" = "Endothelial",  "25" = "Immune",         "26" = "Epithelial",
  "27" = "Epithelial",   "28" = "Immune",         "29" = "Immune",
  "30" = "Endothelial",    "31" = "Immune",         "32" = "Epithelial",
  "33" = "Fibroblast",   "34" = "Endothelial",  "35" = "Smooth_Muscle"
)

dge$cell_type <- unname(cluster_to_celltype[as.character(dge$seurat_clusters)])
table(dge$cell_type, useNA = "always")

celltype_colors <- c(
  "Endothelial" = "#F8766D",
  "Epithelial"    = "#B79F00",
  "Fibroblast"    = "#00BA38",
  "Immune"        = "#00BFC4",
  "Smooth_Muscle" = "#F564E3"
  
)

p <- DimPlot(dge, reduction = "umap", group.by = "cell_type", cols = celltype_colors ) +
  ggtitle(paste0("Full dataset (N=", ncol(dge),")"))
ggsave("UMAP_by_celltype_10thSeptv2.pdf", plot = p, width = 8, height = 6 )

table(dge$cell_type)


