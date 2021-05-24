#remind that you should read the sentences explaining what you are doing

library(dplyr)
library(Seurat)
library(patchwork)

#Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "/home/rstudio/r-test/filtered_gene_bc_matrices/hg19/")
#Initialize the Seurat object with the raw (non-normalized data)
pbmc <- CreateSeuratObject(counts=pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

#The[[ operator can add columns to object matadata. This is a great place to stash QC stats]]
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

#Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#FeatureScatter is typically used to visualize feature-feature relationships, but can be used
#for anything calculated by the object , i.e.columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt <5)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#but the same result will be gained with the line below
#pbmc <- NormalizeData(pbmc)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

#Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

#plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#a standard pre-processing prior to dimensional reduction techniques like PCA
all.genes <- rownames(pbmc)
pbmc<- ScaleData(pbmc, features = all.genes)

#PCA - dimentional reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

#Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

#NOTE: This process can take a long time for big datasets, comment out for expendiency. More
#approximate techniques such as those implemented in ElbowPlot() can be used to reduce
#computation time
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)
ElbowPlot(pbmc)

#Cluster the cells
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

#Look at cluster IDs of the first 5 cells
head(Idents(pbmc),5)

#If you haven't installed UMAP, you can do via reticulate::py_install(packages =
#'umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:10)

#note that you can set 'label = TRUE' or use the LabelClusters function to help label
#indivisual clusters
DimPlot(pbmc, reduction = "umap")
saveRDS(pbmc, file = "/home/rstudio/r-test/output/pbmc_tutorial.rds")

#done here, but you might do it again every time------------
