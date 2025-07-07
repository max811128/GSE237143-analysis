library(Seurat)
library(tidyverse)
NCD<- Read10X(data.dir = "F:/GSE237143_RAW/NCD/")
HFD<- Read10X(data.dir = "F:/GSE237143_RAW/HFD/")

NCD<- CreateSeuratObject(counts = NCD, project = "NCD", min.cells = 3, min.features = 200)
NCD<- PercentageFeatureSet(NCD, pattern = "^mt-", col.name = "percent.mt")
HFD<- CreateSeuratObject(counts = HFD, project = "HFD", min.cells = 3, min.features = 200)
HFD<- PercentageFeatureSet(HFD, pattern = "^mt-", col.name = "percent.mt")

NCD <- subset(NCD, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 10)
HFD<- subset(HFD, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 10)
#VlnPlot(NCD, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#VlnPlot(HFD, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

All.list <- list(NCD=NCD, HFD=HFD) 

All.list <- lapply(X = All.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)})

features <- SelectIntegrationFeatures(object.list = All.list)
All.anchors <- FindIntegrationAnchors(object.list = All.list, anchor.features = features)

All.combined <- IntegrateData(anchorset = All.anchors)
DefaultAssay(All.combined) <- "integrated"

All.combined <- ScaleData(All.combined, verbose = FALSE)
All.combined <- RunPCA(All.combined, npcs = 50, verbose = FALSE)
All.combined <- RunUMAP(All.combined, reduction = "pca", dims = 1:30)
All.combined <- FindNeighbors(All.combined, reduction = "pca", dims = 1:30)
All.combined <- FindClusters(All.combined, resolution = 0.1)

#DimPlot(All.combined, reduction = "umap", label = TRUE)
#DimPlot(All.combined, reduction = "umap", group.by = "orig.ident")
#DimPlot(All.combined, reduction = "umap", split.by = "orig.ident")

## SingleR annotation
library(SingleR)
library(SingleCellExperiment)
library(celldex)
ref<-celldex::MouseRNAseqData()
data.input <- GetAssayData(All.combined, assay = "integrated", layer = "data")
meta <- data.frame(row.names = colnames(All.combined))
sce <- SingleCellExperiment(assays = list(logcounts = data.input), colData = meta)
results <- SingleR(test = sce, ref = ref, labels = ref$label.main)
All.combined$SingleR.label <- results$labels

#DimPlot(All.combined, group.by = "SingleR.label", label = TRUE, reduction = "umap")
#DimPlot(All.combined, group.by = "SingleR.label", label = TRUE, reduction = "umap", split.by = "orig.ident")



