library(Seurat)
library(tidyverse)

HFD<- Read10X(data.dir = "F:/GSE237143_RAW/HFD/")
HFD<- CreateSeuratObject(counts = HFD, project = "HFD", min.cells = 3, min.features = 200)
HFD<- PercentageFeatureSet(HFD, pattern = "^mt-", col.name = "percent.mt")
HFD<- subset(HFD, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 10)

HFD <- NormalizeData(HFD)
HFD <- FindVariableFeatures(HFD, selection.method = "vst", nfeatures = 2000)

HFD <- ScaleData(HFD, features = VariableFeatures(HFD))

HFD <- RunPCA(HFD, features = VariableFeatures(HFD))
HFD <- RunUMAP(HFD, reduction = "pca", dims = 1:30)
HFD <- FindNeighbors(HFD, reduction = "pca", dims = 1:30)
HFD <- FindClusters(HFD, resolution = 0.1)

#Dimplot
DimPlot(HFD, reduction = "umap", label = TRUE)

## SingleR annotation
library(SingleR)
library(SingleCellExperiment)
library(celldex)
ref<-celldex::MouseRNAseqData()
data.input <- GetAssayData(HFD, assay = "RNA", layer = "data")
meta <- data.frame(row.names = colnames(HFD))
sce <- SingleCellExperiment(assays = list(logcounts = data.input), colData = meta)
results <- SingleR(test = sce, ref = ref, labels = ref$label.main)
HFD$SingleR.label <- results$labels

# Dimplot
DimPlot(HFD, group.by = "SingleR.label", label = TRUE, reduction = "umap")

# CellChat
library(CellChat)
library(reticulate)
library(patchwork)
library(circlize)
meta <- data.frame(group = HFD$SingleR.label, row.names = names(HFD$SingleR.label))
data.input <- GetAssayData(HFD, assay = "RNA", layer = "data")
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")

# Load the ligand-receptor interaction database
CellChatDB.mouse <- CellChatDB.mouse
# ShowDatabaseCategory(CellChatDB.mouse)
CellChatDB.use <- CellChatDB.mouse
cellchat@DB <- CellChatDB.use

# Subset and pre-processing the expression data 
# Subset the expression data to use less RAM
cellchat <- subsetData(cellchat)

# Pre-processing the expression data
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- smoothData(cellchat, adj = PPI.mouse)
cellchat <- computeCommunProb(cellchat, raw.use = FALSE)

cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)
cellchat@net$count
cellchat@net$weight

#Circle plot
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1, 2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# Signaling sent from each cell group
mat <- cellchat@net$weight
par(mfrow = c(2, 6), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = F, 
                   edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

# Interaction between specific cell group
mat <- cellchat@net$weight
groupSize <- as.numeric(table(cellchat@idents))
target.groups <- c("Fibroblasts", "Adipocytes", "Macrophages")
target.indices <- which(rownames(mat) %in% target.groups)
par(mfrow = c(1, length(target.indices)), xpd = TRUE)
for (i in target.indices) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  
  netVisual_circle(mat2, 
                   vertex.weight = groupSize, 
                   weight.scale = FALSE, 
                   edge.weight.max = max(mat),)}
