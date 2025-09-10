# Clear workspace and free memory
rm(list = ls())
gc()

# Load required libraries
library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(TxDb.Oaries.Dugesiajaponica.DJwhgtf)

# Set seed for reproducibility
set.seed(123)

# Load genome annotation
TxDb <- TxDb.Oaries.Dugesiajaponica.DJwhgtf
genome <- genes(TxDb, columns = c("tx_id", "tx_name", "gene_id", "TXTYPE"))
genome$gene_name <- as.character(genome$gene_id)
genome$gene_biotype <- "protein_coding"

# Define TSS regions
tss <- promoters(TxDb, upstream = 1, downstream = 0, columns = c("tx_id", "tx_name", "TXTYPE"))
tss$gene_name <- tss$tx_name
tss$gene_id <- tss$gene_name
tss$gene_biotype <- tss$TXTYPE

# Read ATAC-seq data
counts <- Read10X_h5(filename = "../../DJ_ATACaline_zhengda/two/aggre/outs/filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "../../DJ_ATACaline_zhengda/two/aggre/outs/singlecell.csv",
  header = TRUE,
  row.names = 1
)

# Create chromatin assay
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  annotation = genome,
  fragments = '../../DJ_ATACaline_zhengda/two/aggre/outs/fragments.tsv.gz',
  min.cells = 1,
  min.features = 200
)

# Create Seurat object
dj_atac <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

# Extract sample information from cell barcodes
dj_atac$sample <- str_split(colnames(dj_atac), "-", simplify = TRUE)[, 2]
table(dj_atac$sample)

# Calculate QC metrics
dj_atac <- NucleosomeSignal(object = dj_atac)
dj_atac$pct_reads_in_peaks <- dj_atac$peak_region_fragments / dj_atac$passed_filters * 100
dj_atac <- TSSEnrichment(object = dj_atac, tss.positions = tss, fast = FALSE)

# Visualize QC metrics
p1 <- VlnPlot(dj_atac, c('pct_reads_in_peaks', 'peak_region_fragments'), pt.size = 0.1)
p2 <- VlnPlot(dj_atac, c('nucleosome_signal'), pt.size = 0.1) & scale_y_log10()
p1 | p2

# Filter cells based on QC metrics
dj_atac <- subset(
  x = dj_atac,
  subset = peak_region_fragments > 1000 &
    peak_region_fragments < 10000 &
    pct_reads_in_peaks > 20 &
    nucleosome_signal < 5 &
    TSS.enrichment > 1.5
)

# Normalize data and perform dimensionality reduction
dj_atac <- RunTFIDF(dj_atac)
dj_atac <- FindTopFeatures(dj_atac, min.cutoff = 'q5')
dj_atac <- RunSVD(dj_atac)
dj_atac <- RunUMAP(object = dj_atac, reduction = 'lsi', dims = 2:20)
dj_atac <- FindNeighbors(object = dj_atac, reduction = 'lsi', dims = 2:20)
dj_atac <- FindClusters(object = dj_atac, verbose = TRUE, algorithm = 3)

# Classify samples as normal or regeneration
dj_atac$sample_nr <- ifelse(dj_atac$sample == "1" | dj_atac$sample == "2", "normal", "reg5d")
table(dj_atac$sample_nr)

# Calculate gene activity scores
gene.activities <- GeneActivity(dj_atac, assay = "peaks")
dj_atac[['RNA']] <- CreateAssayObject(counts = gene.activities)
dj_atac <- NormalizeData(
  object = dj_atac,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(dj_atac$nCount_RNA)
)

# Find marker genes for each cluster
DefaultAssay(dj_atac) <- 'RNA'
ATAC_singlecell_maker <- FindAllMarkers(dj_atac, only.pos = TRUE, logfc.threshold = 0.3, test.use = "wilcox")

# Annotate cell types based on marker expression
Idents(dj_atac) <- dj_atac$seurat_clusters
new.cluster.ids <- c("Neural", "Neural", "Neural", "Muscle", "Muscle",
                     "Neural", "Epidermal", "Intestine", "Muscle", "Intestine",
                     "Epidermal", "Neural", "Epidermal", "Epidermal", "Muscle",
                     "Parenchymal", "Neoblast", "Intestine", "Intestine", "Neural",
                     "Intestine", "Neoblast")

names(new.cluster.ids) <- levels(dj_atac)
dj_atac <- RenameIdents(dj_atac, new.cluster.ids)
dj_atac$celltype <- Idents(dj_atac)

# Save processed data
save(dj_atac, ATAC_singlecell_maker, file = "DJ_ATAC_processed_data.rdata")
