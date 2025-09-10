# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)
library(org.DDugesiajaponica.eg.db)

# Load spatial transcriptomics data
load("DJ_st_5d_whgtf.rdata")

# Visualize spatial distribution of cells
SpatialDimPlot(sce, alpha = 0.8)

# Highlight intestinal cells in spatial plot
SpatialDimPlot(sce,
               alpha = 0.8,
               cells.highlight = CellsByIdentities(object = sce, idents = c("Intestine")),
               facet.highlight = TRUE,
               ncol = 1) +
  scale_fill_manual(values = c('red', "white"))

# Identify cells in specific spatial region (bottom-left quadrant)
idx <- WhichCells(sce, expression = image_imagerow < 600 | image_imagecol < 0)

# Create new identity classification based on spatial location
Idents(sce) <- 'post_blastema'
Idents(sce, cells = idx) <- 'blastema'
table(Idents(sce))

# Find marker genes between blastema and post_blastema regions
data.markers_blast <- FindAllMarkers(sce, 
                                    only.pos = TRUE, 
                                    logfc.threshold = 0.3,
                                    test.use = "wilcox")
table(data.markers_blast$cluster)

# Extract top marker genes
topgene <- data.markers_blast %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = -p_val_adj)

# Create heatmap of top marker genes
png('DJ_ST_blastema_heatmap.png', 
    width = 7, 
    height = 5, 
    units = "in", 
    res = 300)
DoHeatmap(sce, 
          features = topgene$gene, 
          size = 5) + 
  NoLegend() +
  scale_fill_gradientn(colors = c("white", "grey", "firebrick3"))
dev.off()

# Visualize expression of top markers
top1 <- data.markers_blast %>% 
  group_by(cluster) %>% 
  top_n(n = 5, wt = -p_val_adj)

VlnPlot(sce, 
        features = top1$gene, 
        pt.size = 0, 
        ncol = 5)
ggsave("DJ_ST_blastema_violinplot.pdf", scale = 2)

# Functional enrichment analysis for blastema region
blastema_genes <- data.markers_blast %>%
  filter(cluster == "blastema", 
         p_val_adj < 0.001, 
         avg_log2FC > 0.5) %>%
  pull(gene)

# Load pathway and gene annotation data
load("F:/DJ_tran_rawdata/Go和kegg注释/DJ_KEGG_pathway.rdata")
load("../../../DJ_tran_rawdata/Go和kegg注释/gene_annotation.rdata")

# KEGG pathway enrichment
de_ekp <- enricher(blastema_genes,
                   TERM2GENE = gene2pathway,
                   TERM2NAME = pathway2name,
                   qvalueCutoff = 1,
                   pvalueCutoff = 1)
barplot(de_ekp, showCategory = 20)

# GO enrichment analysis
ego_blast <- enrichGO(gene = blastema_genes,
                      OrgDb = org.DDugesiajaponica.eg.db,
                      keyType = "GID",
                      ont = "all",
                      qvalueCutoff = 1,
                      pvalueCutoff = 1)

# Functional enrichment analysis for post_blastema region
post_blastema_genes <- data.markers_blast %>%
  filter(cluster == "post_blastema",
         p_val_adj < 0.001,
         avg_log2FC > 0.5) %>%
  pull(gene)

ego_nomal <- enrichGO(gene = post_blastema_genes,
                      OrgDb = org.DDugesiajaponica.eg.db,
                      keyType = "GID",
                      ont = "all",
                      qvalueCutoff = 1,
                      pvalueCutoff = 1)

# Prepare data for comparative visualization
ego_blast_dat <- ego_blast@result %>%
  head(5) %>%
  select(Description, pvalue) %>%
  mutate(dat = "blast")

ego_nomal_dat <- ego_nomal@result %>%
  head(5) %>%
  select(Description, pvalue) %>%
  mutate(dat = "post_blastema")

de_ekp_dat <- rbind(ego_blast_dat, ego_nomal_dat)

# Create comparative bar plot
p <- ggplot(de_ekp_dat, 
            aes(x = reorder(Description, -log10(pvalue)),
                y = -log10(pvalue), 
                fill = dat)) +
  geom_bar(stat = 'identity') +
  theme_classic() +
  coord_flip() +
  labs(x = "GO Term", y = "-log10(p-value)") +
  theme(axis.text = element_text(size = 10))

ggsave("DJ_ST_blastema_post_blastema_GO_enrichment.pdf", p, scale = 2)

# Save all results to RData file
save(sce, data.markers_blast, ego_blast, ego_nomal, de_ekp_dat, 
     file = "DJ_ST_spatial_analysis_results.rdata")
