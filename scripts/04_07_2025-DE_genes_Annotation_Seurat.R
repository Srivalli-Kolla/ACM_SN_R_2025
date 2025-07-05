library(Seurat)
library(dplyr)
library(ggplot2)
library(SeuratDisk)

#### Load processed integrated object ####
integrated <- readRDS("/home/gruengroup/srivalli/Github/ACM_SN_R_2025/data/integrated_processed.rds")

#### Perform DE analysis ####
markers <- FindAllMarkers(
  integrated,
  only.pos = TRUE,
  test.use = "wilcox",
  min.pct = 0.25,
  logfc.threshold = 0.25
)

top5 <- markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)

# Save DE marker results
write.csv(markers, "/home/gruengroup/srivalli/Github/ACM_SN_R_2025/data/de_all_markers.csv", row.names = FALSE)
write.csv(top10, "/home/gruengroup/srivalli/Github/ACM_SN_R_2025/data/de_top_10_clusters.csv", row.names = FALSE)
saveRDS(markers, file = "/home/gruengroup/srivalli/Github/ACM_SN_R_2025/data/de_all_markers.rds")

# Dot plot for top 5 marker genes
p_dot <- DotPlot(integrated, features = unique(top5$gene)) + RotatedAxis()
ggsave("/home/gruengroup/srivalli/Github/ACM_SN_R_2025/data/dotplot_top5_markers.png", p_dot, width = 10, height = 6)

#### Manually set human-readable cluster names ####
cluster_names <- c(
  "Fibroblasts",                        # 0
  "Endothelial/Progenitor Cells",      # 1
  "Cardiomyocytes (ventricular)",      # 2
  "Mast Cells",                         # 3
  "Epicardial/Endothelial Cells",      # 4
  "Cardiomyocytes (possibly atrial)",  # 5
  "Pericytes/Vascular SMC",            # 6
  "Lymphatic Endothelial Cells",       # 7
  "Endocardial/Conduction Cells",      # 8
  "Immune Cells (Macrophages)",        # 9
  "Pericytes/Fibroblasts",             # 10
  "B Lymphocytes",                     # 11
  "T Lymphocytes",                     # 12
  "Vascular Smooth Muscle Cells",     # 13
  "Proliferating/Cycling Cells",      # 14
  "Neural Crest/Schwann Cells",       # 15
  "Epicardial/Mesothelial Cells"      # 16
)

# Map readable cluster names
integrated$cluster_named <- factor(
  cluster_names[as.numeric(as.character(integrated$seurat_clusters)) + 1],
  levels = cluster_names
)

#### UMAP Plots ####
dir.create("/home/gruengroup/srivalli/Github/ACM_SN_R_2025/plots", showWarnings = FALSE)

# 1. Seurat Cluster ID
p0 <- DimPlot(integrated, group.by = "seurat_clusters", label = TRUE, repel = TRUE) +
  ggtitle("UMAP by Seurat Cluster ID")
ggsave("/home/gruengroup/srivalli/Github/ACM_SN_R_2025/plots/umap_seurat_clusters.png", p0, width = 10, height = 7)

# 2. Named clusters
p1 <- DimPlot(integrated, group.by = "cluster_named", label = TRUE, repel = TRUE) +
  ggtitle("UMAP by Annotated Cluster")
ggsave("/home/gruengroup/srivalli/Github/ACM_SN_R_2025/plots/umap_named_clusters.png", p1, width = 10, height = 7)

# 3. Sample-wise UMAP
if ("sample" %in% colnames(integrated@meta.data)) {
  p2 <- DimPlot(integrated, group.by = "sample", label = FALSE) +
    ggtitle("UMAP by Sample")
  ggsave("/home/gruengroup/srivalli/Github/ACM_SN_R_2025/plots/umap_by_sample.png", p2, width = 10, height = 7)
}

# 4. Genotype-wise UMAP
if ("genotype" %in% colnames(integrated@meta.data)) {
  p3 <- DimPlot(integrated, group.by = "genotype", label = FALSE) +
    ggtitle("UMAP by Genotype")
  ggsave("/home/gruengroup/srivalli/Github/ACM_SN_R_2025/plots/umap_by_genotype.png", p3, width = 10, height = 7)
}

#### Save annotated Seurat object ####
saveRDS(integrated, "/home/gruengroup/srivalli/Github/ACM_SN_R_2025/data/integrated_named_clusters.rds")
SaveH5Seurat(integrated, filename = "/home/gruengroup/srivalli/Github/ACM_SN_R_2025/data/integrated_named_clusters.h5Seurat", overwrite = TRUE)
Convert("/home/gruengroup/srivalli/Github/ACM_SN_R_2025/data/integrated_named_clusters.h5Seurat", dest = "h5ad", overwrite = TRUE)