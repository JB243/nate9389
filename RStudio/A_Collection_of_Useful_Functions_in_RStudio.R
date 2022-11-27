library(Seurat)
library(dplyr)
library(ggplot2)

load('cancer_lung.rdata')
library(Matrix)
cancer_lung_sp= as(as.matrix(cancer_lung), "sparseMatrix")

writeMM(cancer_lung_sp,file='cancer_lung_sp.mtx')
#readMM(file='cancer_lung_sp.mtx')
load('cancer_lung_id.rdata')


#Seurat
cancer_seurat <- CreateSeuratObject(counts = cancer_lung_sp,
                             min.cells = 3,
                             project = "cancer_lung")
cancer_seurat[["percent.mt"]] <- PercentageFeatureSet(cancer_seurat, pattern = "^MT")
FeatureScatter(cancer_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(cancer_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

cancer_seurat <- NormalizeData(object = cancer_seurat,
                        normalization.method = "LogNormalize",
                        scale.factor = 1e4)

cancer_seurat <- FindVariableFeatures(object = cancer_seurat,
                               selection.method = "vst", nfeatures = 2000)
VariableFeaturePlot(object = cancer_seurat)


cancer_seurat <- ScaleData(object=cancer_seurat)
cancer_seurat <- RunPCA(object = cancer_seurat, features = VariableFeatures(object = cancer_seurat))
DimPlot(cancer_seurat, reduction = "pca")

cancer_seurat <- FindNeighbors(object = cancer_seurat, dims = 1:5)
cancer_seurat <- FindClusters(object = cancer_seurat, 
                       resolution = 0.1)

cancer_seurat <- RunTSNE(object = cancer_seurat, dims = 1:10, do.fast= TRUE, check_duplicates = FALSE)
cancer_seurat <- RunUMAP(cancer_seurat, dims = 1:10)
DimPlot(cancer_seurat, reduction = "umap")
DimPlot(cancer_seurat, reduction = "tsne")

save(cancer_seurat, file = 'cancer_seurat.rdata')

load('cancer_seurat.rdata')
cancer.markers <- FindAllMarkers(object = cancer_seurat,
                                 only.pos = TRUE,
                                 min.pct = 0.25,
                                 thresh.use = 0.25)
save(cancer.markers, file = 'cancer.markers.rdata')

load('cancer_lung_cell.rdata')

cancer_seurat <- CreateSeuratObject(counts = cancer_lung_cell,
                                    min.cells = 3,
                                    project = "cancer_lung")
cancer_seurat[["percent.mt"]] <- PercentageFeatureSet(cancer_seurat, pattern = "^MT")
FeatureScatter(cancer_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(cancer_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

cancer_seurat <- NormalizeData(object = cancer_seurat,
                               normalization.method = "LogNormalize",
                               scale.factor = 1e4)

cancer_seurat <- FindVariableFeatures(object = cancer_seurat,
                                      selection.method = "vst", nfeatures = 2000)
VariableFeaturePlot(object = cancer_seurat)


cancer_seurat <- ScaleData(object=cancer_seurat)
cancer_seurat <- RunPCA(object = cancer_seurat, features = VariableFeatures(object = cancer_seurat))
DimPlot(cancer_seurat, reduction = "pca")

cancer_seurat <- FindNeighbors(object = cancer_seurat, dims = 1:5)
cancer_seurat <- FindClusters(object = cancer_seurat, 
                              resolution = 0.2)

cancer_seurat <- RunTSNE(object = cancer_seurat, dims = 1:10, do.fast= TRUE, check_duplicates = FALSE)
cancer_seurat <- RunUMAP(cancer_seurat, dims = 1:10)
DimPlot(cancer_seurat, reduction = "umap")
DimPlot(cancer_seurat, reduction = "tsne")

cancer_epi_seurat = cancer_seurat
save(cancer_epi_seurat, file = 'cancer_epi_seurat.rdata')

load('cancer_epi_seurat.rdata')
cancer.epi.markers <- FindAllMarkers(object = cancer_epi_seurat,
                                 only.pos = TRUE,
                                 min.pct = 0.25,
                                 thresh.use = 0.25)
save(cancer.epi.markers, file = 'cancer_epi_markers.rdata')


###
load('normal_lung.rdata')

normal_seurat <- CreateSeuratObject(counts = normal_lung,
                                    min.cells = 3,
                                    project = "normal_lung")
normal_seurat[["percent.mt"]] <- PercentageFeatureSet(normal_seurat, pattern = "^MT")
FeatureScatter(normal_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(normal_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

normal_seurat <- NormalizeData(object = normal_seurat,
                               normalization.method = "LogNormalize",
                               scale.factor = 1e4)

normal_seurat <- FindVariableFeatures(object = normal_seurat,
                                      selection.method = "vst", nfeatures = 2000)
VariableFeaturePlot(object = normal_seurat)


normal_seurat <- ScaleData(object=normal_seurat)
normal_seurat <- RunPCA(object = normal_seurat, features = VariableFeatures(object = normal_seurat))
DimPlot(normal_seurat, reduction = "pca")

normal_seurat <- FindNeighbors(object = normal_seurat, dims = 1:5)
normal_seurat <- FindClusters(object = normal_seurat, 
                              resolution = 0.1)

normal_seurat <- RunTSNE(object = normal_seurat, dims = 1:10, do.fast= TRUE, check_duplicates = FALSE)
normal_seurat <- RunUMAP(normal_seurat, dims = 1:10)
DimPlot(normal_seurat, reduction = "umap")
DimPlot(normal_seurat, reduction = "tsne")

save(normal_seurat, file = 'normal_seurat.rdata')

load('normal_seurat.rdata')
normal.markers <- FindAllMarkers(object = normal_seurat,
                                 only.pos = TRUE,
                                 min.pct = 0.25,
                                 thresh.use = 0.25)
save(normal.markers, file = 'normal.markers.rdata')

load('normal_lung_cell.rdata')

normal_seurat <- CreateSeuratObject(counts = normal_lung_cell,
                                    min.cells = 3,
                                    project = "normal_lung")
normal_seurat[["percent.mt"]] <- PercentageFeatureSet(normal_seurat, pattern = "^MT")
FeatureScatter(normal_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(normal_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

normal_seurat <- NormalizeData(object = normal_seurat,
                               normalization.method = "LogNormalize",
                               scale.factor = 1e4)

normal_seurat <- FindVariableFeatures(object = normal_seurat,
                                      selection.method = "vst", nfeatures = 2000)
VariableFeaturePlot(object = normal_seurat)


normal_seurat <- ScaleData(object=normal_seurat)
normal_seurat <- RunPCA(object = normal_seurat, features = VariableFeatures(object = normal_seurat))
DimPlot(normal_seurat, reduction = "pca")

normal_seurat <- FindNeighbors(object = normal_seurat, dims = 1:5)
normal_seurat <- FindClusters(object = normal_seurat, 
                              resolution = 0.5)

normal_seurat <- RunTSNE(object = normal_seurat, dims = 1:10, do.fast= TRUE, check_duplicates = FALSE)
normal_seurat <- RunUMAP(normal_seurat, dims = 1:10)
DimPlot(normal_seurat, reduction = "umap")
DimPlot(normal_seurat, reduction = "tsne")

normal_epi_seurat = normal_seurat
save(normal_epi_seurat, file = 'normal_epi_seurat.rdata')

normal.epi.markers <- FindAllMarkers(object = normal_seurat,
                                 only.pos = TRUE,
                                 min.pct = 0.25,
                                 thresh.use = 0.25)
save(normal.epi.markers, file = 'normal_epi_markers.rdata')
