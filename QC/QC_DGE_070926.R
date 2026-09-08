library(Seurat)
library(data.table)
library(DoubletFinder)
library(dplyr)

dge_files <- list.files("C:/Users/rohin/Desktop/Pre_doc_projects/Project_RB/Huang_Lab/GSE293189_RAW_extracted",
                        pattern = "_dge.txt$", full.names = TRUE)

seurat_list <- list()

for (f in dge_files) {
  sample_id <- gsub("_dge.txt$", "", basename(f))
  mat <- fread( f, sep = "\t", header =TRUE, data.table= FALSE)
  rownames(mat) <- mat$GENE
  mat$GENE <- NULL
  mat <- as.matrix(mat)
  obj <- CreateSeuratObject(counts = mat, project = sample_id)
  obj$sample <- sample_id
  seurat_list[[sample_id]] <- obj
}
length(seurat_list)
sapply(seurat_list , ncol)

for (id in names(seurat_list)) {
  obj <- seurat_list[[id]]
  obj$percent.mt <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj <- subset(obj, subset = nFeature_RNA > 300 &
                  nCount_RNA > 500 &
                  percent.mt < 20)
  seurat_list[[id]] <- obj
}

sapply(seurat_list, ncol)

for (id in names(seurat_list)) {
  obj <- seurat_list[[id]]
  
  if (ncol(obj) < 50) {
    message(id, " has only ", ncol(obj), " cells after QC — skipping DoubletFinder, keeping as-is")
    next
}
    obj <- NormalizeData(obj)
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures= 2000 )
    obj <- ScaleData(obj, features= VariableFeatures(obj))
    obj <- RunPCA(obj, features = VariableFeatures(obj),npcs = min(30, ncol(obj) -1) )

    sweep.res <- paramSweep(obj, PCs= 1:10, sct = FALSE)
    sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    pK_opt <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

    nExp <- round(0.075 *ncol(obj))
    obj <- doubletFinder(obj, PCs = 1:10, pN = 0.25, pK = pK_opt, nExp = nExp)

    df_col <- grep("^DF.classifications", colnames(obj@meta.data), value = TRUE )
    obj <- subset(obj, cells = colnames(obj)[obj[[df_col]] == "Singlet"])
    seurat_list[[id]] <- obj
}

sapply(seurat_list, ncol)

dge_merged <- merge(seurat_list[[1]], y= seurat_list[-1],
                     add.cell.ids = names(seurat_list))
ncol(dge_merged)
table(dge_merged$sample)

cols_to_drop <-grep("^pANN_|^DF.classifications_", colnames(dge_merged@meta.data), value = TRUE)
length(cols_to_drop)

dge_merged@meta.data <- dge_merged@meta.data[,!colnames(dge_merged@meta.data) %in% cols_to_drop]
colnames(dge_merged@meta.data)

saveRDS(dge_merged, "C:/Users/rohin/Desktop/Pre_doc_projects/Project_RB/Huang_Lab/dge_mergedall_090726.rds" )
sapply(seurat_list, ncol)


