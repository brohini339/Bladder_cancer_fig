These are two R scripts I have used to analyse the scRNA seq data publicly available under accession number GSE293189 from the paper entitled "Bladder cancer variants share aggressive
features including a CA125+ cell state and targetable TM4SF1 expression"(https://doi.org/10.1038/s41467-025-59888-8) and the corresponding figures generated from this dataset. 
I have generated Fig 2B of the supplementary data provided with this paper , i.e. the UMAP of whole dataset , colour-coded by broad cell type.
I have also generated the raw UMAP of the whole dataset and colour-coded by tumor specimen type.

The following folders are present in this  repository:

QC- Contains QC_DGE_070926.R , the script for QC filtering according to the paper's thresholds , doublet removal and saving merged samples in a single .rds  file.

UMAP_wholedata - Contains UMAP_generation_11thSept.R , the script for  ambient RNA decontamination , normalization, and clustering pipeline for the QC'd, doublet-filtered bladder cancer scRNA-seq dataset , reproducing the published methods and annotating clusters into five major cell types via canonical marker expression.

mainpaper_and_supplementarydata -Contains the following files:

s41467-025-59888-8.pdf- the main paper

41467_2025_59888_MOESM1_ESM.pdf - the supplementary data provided with the paper

Results- Contains the following generated results:
VariableFeature_2000_11thSept.pdf - Scatter plot of the top 2000 variable genes (mean expression vs. dispersion), with the 10 most variable genes labelled.
PCA_byName11thSept.pdf- Cells in PCA space (PC1 vs PC2), coloured by patient group (Name). 
PCA_heatmap11thSept.png -Heatmap of the top genes driving each of the first 10 principal components for 500 representative cells. 
Elbowplot_dge11thSept.png-   Standard deviation explained by each of the 100 stored PCs. 
Stability11thSept.pdf - Clustree plot showing how cluster assignments shift across resolutions (0.2–1.0).
TSNE_raw11thSept.pdf - t-SNE embedding of all cells (75 PCs), labelled by cluster number, no cell-type colouring. 
Umap_raw11thSept.pdf- UMAP embedding of all cells (75 PCs), labelled by cluster number. 
Umap_byPatient11thSept.pdf- Same UMAP, coloured by patient group (Name) instead of cluster.
Marker_DotPlot11thSept.pdf- Dot plot of canonical marker gene expression (percent expressed + average expression) across all 36 clusters, used to assign cell-type identity. 
Marker_Feature11thSept.pdf -UMAP feature plots showing expression of one representative marker per cell type (EPCAM, PTPRC, COL1A1, PECAM1, ACTA2). 
UMAP_by_celltype_11thSep.pdf - Final annotated UMAP coloured by assigned cell type (Epithelial, Immune, Fibroblast, Endothelial, Smooth Muscle) , the main result figure, comparable to Fig 2B of the supplementary data of the paper.






