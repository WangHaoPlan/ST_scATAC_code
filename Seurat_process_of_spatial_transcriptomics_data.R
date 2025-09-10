####
# Clear workspace and garbage collect
rm(list = ls())
gc()

# Load required libraries
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(harmony)

# Set random seed for reproducibility
set.seed(123)

# Function to load and process each spatial dataset
load_spatial_data <- function(expr_path, img_path, project_name, slice_num) {
  # Read expression data
  expr_data <- Read10X_h5(filename = expr_path)
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(
    counts = expr_data, 
    project = project_name, 
    assay = 'Spatial'
  )
  
  # Add metadata
  seurat_obj$slice <- slice_num
  seurat_obj$region <- project_name
  
  # Read and process spatial image
  img <- Read10X_Image(image.dir = img_path)
  DefaultAssay(object = img) <- 'Spatial'
  img <- img[colnames(x = seurat_obj)]
  seurat_obj[['image']] <- img
  
  return(seurat_obj)
}

# Load four spatial datasets
mydata1 <- load_spatial_data(
  "f:/DJ_spital/plan5d/DJ_plan1_A_20230503/outs/filtered_feature_bc_matrix.h5",
  "f:/DJ_spital/plan5d/DJ_plan1_A_20230503/outs/spatial/",
  'DJ_1_A',
  1
)

mydata2 <- load_spatial_data(
  "f:/DJ_spital/plan5d/DJ_plan1_B_20230503/outs/filtered_feature_bc_matrix.h5",
  "f:/DJ_spital/plan5d/DJ_plan1_B_20230503/outs/spatial/",
  'DJ_1_B',
  2
)

mydata3 <- load_spatial_data(
  "f:/DJ_spital/plan5d/DJ_plan1_C_20230503/outs/filtered_feature_bc_matrix.h5",
  "f:/DJ_spital/plan5d/DJ_plan1_C_20230503/outs/spatial/",
  'DJ_1_C',
  3
)

mydata4 <- load_spatial_data(
  "f:/DJ_spital/plan5d/DJ_plan1_D_20230503/outs/filtered_feature_bc_matrix.h5",
  "f:/DJ_spital/plan5d/DJ_plan1_D_20230503/outs/spatial/",
  'DJ_1_D',
  4
)

# Merge all datasets
mydata <- merge(mydata1, c(mydata2, mydata3, mydata4))

# Quality control plots
plot1 <- VlnPlot(mydata, features = "nCount_Spatial", pt.size = 0.1)
plot2 <- SpatialFeaturePlot(mydata, features = "nCount_Spatial") + theme(legend.position = "right")
plot_grid(plot1, plot2)

plot1 <- VlnPlot(mydata, features = "nFeature_Spatial", pt.size = 0.1) 
plot2 <- SpatialFeaturePlot(mydata, features = "nFeature_Spatial") + theme(legend.position = "right")
plot_grid(plot1, plot2)

FeatureScatter(mydata, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")

# Filter cells based on QC metrics
mydata_sub <- subset(mydata, 
                    subset = nFeature_Spatial > 300 & nFeature_Spatial < 7500 &
                             nCount_Spatial > 500 & nCount_Spatial < 60000)

# Normalize data using SCTransform
mydata_sub <- SCTransform(mydata_sub, assay = "Spatial", verbose = FALSE)

# Dimensionality reduction and clustering
mydata_sub <- RunPCA(mydata_sub, assay = "SCT")
mydata_sub <- FindNeighbors(mydata_sub, reduction = "pca", dims = 1:30)
mydata_sub <- FindClusters(mydata_sub, resolution = 0.5)
mydata_sub <- RunUMAP(mydata_sub, reduction = "pca", dims = 1:30)
mydata_sub <- RunTSNE(mydata_sub, reduction = "pca", dims = 1:30)

# Batch correction using Harmony
seuratObj <- RunHarmony(mydata_sub, "orig.ident")
seuratObj <- RunUMAP(seuratObj, dims = 1:30, reduction = "harmony")
seuratObj <- RunTSNE(seuratObj, dims = 1:30, reduction = "harmony")

# Final clustering
sce <- FindNeighbors(seuratObj, reduction = "harmony", dims = 1:25)
sce <- FindClusters(sce, resolution = 0.4)

# Find marker genes
data.markers <- FindAllMarkers(sce, 
                              only.pos = TRUE, 
                              min.pct = 0.3, 
                              logfc.threshold = 0.3,
                              test.use = "wilcox")

# Rename clusters with cell type annotations
new.cluster.ids <- c("Intestine", "Neoblast", "Epidermal", "Muscle", "Intestine", "Neural")
names(new.cluster.ids) <- levels(sce)
sce <- RenameIdents(sce, new.cluster.ids)
sce$celltype <- Idents(sce)

# Save final data
save(sce, data.markers, file = "DJ_st_5d_whgtf_final.rdata")
saveRDS(sce, file = "dj_STwhgtf_5D_final.rds")

# Export marker genes
write.csv(data.markers, 
          file = "DJ_st_5d_markers_celltype.csv", 
          row.names = FALSE, 
          quote = FALSE)
