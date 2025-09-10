# Clear workspace and free up memory
rm(list = ls())
gc()

# Load required libraries
library(Seurat)
library(dplyr)
library(stringr)

# Set random seed for reproducibility
set.seed(123)

# Load and process single-cell RNA-seq data
sce_dj_singlecell <- Read10X('F:/DJ_singlecell/result/outs/count/filtered_feature_bc_matrix/')
sce_dj_singlecell <- CreateSeuratObject(
  counts = sce_dj_singlecell,
  project = "dj_reg",
  min.cells = 10,
  min.features = 200
)

# Quality control: Check relationship between UMI counts and gene counts
FeatureScatter(sce_dj_singlecell, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Extract time information from cell barcodes and assign to metadata
sce_dj_singlecell$time <- str_split(colnames(sce_dj_singlecell), "-", simplify = TRUE)[, 2]
table(sce_dj_singlecell$time)

# Rename time points for clarity
Idents(sce_dj_singlecell) <- sce_dj_singlecell$time
Idents(sce_dj_singlecell, cells = WhichCells(sce_dj_singlecell, idents = "1")) <- 'intact'
Idents(sce_dj_singlecell, cells = WhichCells(sce_dj_singlecell, idents = "2")) <- '5dpa'
sce_dj_singlecell$time <- Idents(sce_dj_singlecell)
table(sce_dj_singlecell$time)

# Visualize QC metrics
VlnPlot(sce_dj_singlecell, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

# Calculate percentage of hemoglobin genes (potential red blood cell markers)
HB.genes <- c("DJKW.19663", "DJKW.1695")
HB_m <- match(HB.genes, rownames(sce_dj_singlecell@assays$RNA))
HB.genes <- rownames(sce_dj_singlecell@assays$RNA)[HB_m]
HB.genes <- HB.genes[!is.na(HB.genes)]
sce_dj_singlecell[["percent.HB"]] <- PercentageFeatureSet(sce_dj_singlecell, features = HB.genes)

# Visualize QC metrics including hemoglobin percentage
VlnPlot(sce_dj_singlecell, features = c("nFeature_RNA", "nCount_RNA", "percent.HB"), ncol = 3)

# Filter cells based on QC metrics
sce_sub <- subset(sce_dj_singlecell, 
                  subset = (nFeature_RNA > 300 & nFeature_RNA < 4000 & 
                             nCount_RNA < 20000 & nCount_RNA > 2000))
table(sce_sub$time)

# Normalize data
sce_sub <- NormalizeData(sce_sub, 
                         normalization.method = "LogNormalize",
                         scale.factor = 10000)

# Identify highly variable features
sce_sub <- FindVariableFeatures(sce_sub, 
                                selection.method = "vst", 
                                nfeatures = 2000)

# Scale data
sce_sub <- ScaleData(sce_sub)

# Perform PCA
sce_sub <- RunPCA(sce_sub, features = VariableFeatures(object = sce_sub))

# Find neighbors and clusters
sce_sub <- FindNeighbors(sce_sub, dims = 1:30)

# Remove doublets using DoubletFinder
library(DoubletFinder)
keloid <- sce_sub

# Parameter sweep for doublet detection
sweep.res.list_keloid <- paramSweep_v3(keloid, PCs = 1:30, sct = FALSE)
sweep.stats_keloid <- summarizeSweep(sweep.res.list_keloid, GT = FALSE)
bcmvn_keloid <- find.pK(sweep.stats_keloid)
mpK <- as.numeric(as.vector(bcmvn_keloid$pK[which.max(bcmvn_keloid$BCmetric)]))

# Estimate doublet rate
annotations <- keloid@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
DoubletRate = ncol(keloid) * 8 * 1e-6  # 8 doublets per 1000 cells
nExp_poi <- round(DoubletRate * ncol(keloid))
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

# Run DoubletFinder
keloid <- doubletFinder_v3(keloid, 
                           PCs = 1:30, 
                           pN = 0.25, 
                           pK = mpK, 
                           nExp = nExp_poi, 
                           reuse.pANN = FALSE, 
                           sct = FALSE)

# Process doublet classification results
keloid@meta.data[, "DF_hi.lo"] <- keloid@meta.data$DF.classifications_0.25_0.28_1016
keloid@meta.data$DF_hi.lo[which(keloid@meta.data$DF_hi.lo == "Doublet")] <- "Doublet-High Confidence"
table(keloid@meta.data$DF_hi.lo)

# Visualize doublets
DimPlot(keloid, reduction = "umap", group.by = "DF_hi.lo", cols = c("black", "red"))

# Remove doublets
pool_sub_remove_double <- subset(keloid, subset = DF_hi.lo == "Singlet")
table(pool_sub_remove_double$time)

# Re-cluster after doublet removal
pool_sub_remove_double <- FindNeighbors(pool_sub_remove_double, dims = 1:30)
pool_sub_remove_double <- RunTSNE(pool_sub_remove_double, dims.use = 1:30)
pool_sub_remove_double <- RunUMAP(pool_sub_remove_double, dims = 1:30)
pool_sub_remove_double <- FindClusters(pool_sub_remove_double, resolution = 0.3)

# Visualize clusters by time and identity
p1 <- DimPlot(pool_sub_remove_double, reduction = "umap", group.by = "time")
p2 <- DimPlot(pool_sub_remove_double, reduction = "umap", label = TRUE)
CombinePlots(plots = list(p1, p2))

# Find marker genes for each cluster
Idents(pool_sub_remove_double) <- pool_sub_remove_double$seurat_clusters
sce.markers <- FindAllMarkers(pool_sub_remove_double, 
                              only.pos = TRUE, 
                              min.pct = 0.25, 
                              logfc.threshold = 0.25)

# Get top markers
top10 <- sce.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

# Annotate cell types based on marker expression
Idents(pool_sub_remove_double) <- pool_sub_remove_double$seurat_clusters
new.cluster.ids <- c("Neural", "Muscle", "Epidermal", "Neural", "Neural",
                     "Cathepsin+ cell", "Cathepsin+ cell", "Protonephridia", "Epidermal", "Neural",
                     "Neoblast", "Intestine", "Neural", "Epidermal", "Intestine",
                     "Cathepsin+ cell", "Neural", "Cathepsin+ cell", "Cathepsin+ cell", "Parenchymal")

names(new.cluster.ids) <- levels(pool_sub_remove_double)
pool_sub_remove_double <- RenameIdents(pool_sub_remove_double, new.cluster.ids)
pool_sub_remove_double$celltype <- Idents(pool_sub_remove_double)

# Visualize cell types
DimPlot(pool_sub_remove_double, reduction = "umap", label = TRUE, pt.size = 0.5)
DimPlot(pool_sub_remove_double, reduction = "umap", group.by = "time", label = TRUE, pt.size = 1)

# Find markers for annotated cell types
pool_sub_remove_double.markers_celltype <- FindAllMarkers(pool_sub_remove_double, 
                                                         only.pos = FALSE, 
                                                         min.pct = 0.2, 
                                                         logfc.threshold = 0.2)

# Rename samples for clarity
Idents(pool_sub_remove_double) <- pool_sub_remove_double$orig.ident
new.cluster.ids <- c("5dpa", "intact")
names(new.cluster.ids) <- levels(pool_sub_remove_double)
pool_sub_remove_double <- RenameIdents(pool_sub_remove_double, new.cluster.ids)
levels(pool_sub_remove_double) <- c("intact", "5dpa")
pool_sub_remove_double$sample <- Idents(pool_sub_remove_double)

# Final visualization
DimPlot(pool_sub_remove_double, reduction = "umap", label = TRUE, pt.size = 0.5)
DimPlot(pool_sub_remove_double, split.by = "sample", reduction = "umap", label = TRUE, pt.size = 0.5)

# Save final object
dj_singcell <- pool_sub_remove_double
dj_singcell_markers_celltype <- pool_sub_remove_double.markers_celltype
save(dj_singcell, dj_singcell_markers_celltype, file = "DJ_singlecell_processed.rdata")
