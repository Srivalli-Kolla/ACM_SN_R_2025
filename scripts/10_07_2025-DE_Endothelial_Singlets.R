library(Seurat)
library(dplyr)
library(ggplot2)
library(SeuratDisk)

#### Load processed integrated object and Select Endothelial only ####
integrated_singlets <- readRDS("/Volumes/srivalli/dzhi/ACM_SN_R_2025/data/integrated_named_clusters_singlets.rds")
DefaultAssay(integrated_singlets) <- "integrated"
en_labels <- c(
  "Capillary Endothelial Cells",
  "Lymphatic Endothelial Cells",
  "Endocardial Cells (conduction system)")
en_cells <- colnames(integrated_singlets)[integrated_singlets$cluster_named %in% en_labels ]
en_singlets <- subset(integrated_singlets, cells = en_cells)

#### Reclustering Endothelial ####
DefaultAssay(en_singlets) <- "RNA"
en_singlets <- NormalizeData(en_singlets)
en_singlets <- FindVariableFeatures(en_singlets)
en_singlets <- ScaleData(en_singlets)
en_singlets <- RunPCA(en_singlets)
en_singlets <- RunUMAP(en_singlets, dims = 1:10)
en_singlets <- FindNeighbors(en_singlets, dims = 1:10)
en_singlets <- FindClusters(en_singlets, resolution = 0.3) 

p0 <- DimPlot(en_singlets, label = TRUE, group.by = "seurat_clusters") + ggtitle("Endothelial Subclusters")
p0
ggsave("/Volumes/srivalli/dzhi/ACM_SN_R_2025/plots/umap_en_clusters.png", p0, width = 10, height = 7)

#### DE analysis ####
en_markers <- FindAllMarkers(
  en_singlets,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

top5 <- en_markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
top10 <- en_markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)

# Save DE marker results
write.csv(en_markers, "/Volumes/srivalli/dzhi/ACM_SN_R_2025/data/de_en_markers_singlets.csv", row.names = FALSE)
write.csv(top10, "/Volumes/srivalli/dzhi/ACM_SN_R_2025/data/de_top_10_en_singlets.csv", row.names = FALSE)

#### Manually set human-readable cluster names ####
cluster_names <- c(
  "Ventricular Cardiomyocytes",
  "Capillary Endothelial Cells",
  "Lymphatic Endothelial Cells",
  "Endocardial Cells (conduction system)",
  "Fibroblasts (matrix/adventitial)",
  "Erythroid-like Cells (contaminant)",
  "Fibroblasts (matricellular/fibrotic/epicardial-like)",
  "Cardiac Macrophages/Monocytes"
)

# Map readable cluster names
en_singlets$cluster_named <- factor(
  cluster_names[as.numeric(as.character(en_singlets$seurat_clusters)) + 1],
  levels = cluster_names
)
# UMAPS
p_dot <- DotPlot(en_singlets, features = unique(top5$gene)) + RotatedAxis()
p_dot
ggsave("/Volumes/srivalli/dzhi/ACM_SN_R_2025/plots/dotplot_top5_en_markers_singlets.png", p_dot, width = 20, height = 10)

# 2. Named clusters
p1 <- DimPlot(en_singlets, group.by = "cluster_named", repel = TRUE,raster = FALSE) +
  ggtitle("UMAP by Annotated Cluster - Endothelial")
p1
ggsave("/Volumes/srivalli/dzhi/ACM_SN_R_2025/plots/umap_en_named_clusters_singlets.png", p1, width = 10, height = 7)

# 3. Sample-wise UMAP
if ("Sample_ID" %in% colnames(en_singlets@meta.data)) {
  p2 <- DimPlot(en_singlets, group.by = "Sample_ID", label = FALSE,raster = FALSE) +
    ggtitle("UMAP by Sample")
  p2
  ggsave("/Volumes/srivalli/dzhi/ACM_SN_R_2025/plots/umap_en_by_sample_singlets.png", p2, width = 10, height = 7)
}

# 4. Genotype-wise UMAP
if ("Genotype" %in% colnames(en_singlets@meta.data)) {
  p3 <- DimPlot(en_singlets, group.by = "Genotype", label = FALSE,raster = FALSE) +
    ggtitle("UMAP by Genotype")
  p3
  ggsave("/Volumes/srivalli/dzhi/ACM_SN_R_2025/plots/umap_en_by_genotype_singlets.png", p3, width = 10, height = 7)
}

# 5. Condition-wise UMAP
if ("Condition" %in% colnames(en_singlets@meta.data)) {
  p4 <- DimPlot(en_singlets, group.by = "Condition", label = FALSE,raster = FALSE) +
    ggtitle("UMAP by Condition")
  p4
  ggsave("/Volumes/srivalli/dzhi/ACM_SN_R_2025/plots/umap_by_en_Condition_singlets.png", p4, width = 10, height = 7)
}

# Split UMAP by genotype
p_umap_split <- DimPlot(
  en_singlets,
  group.by = "cluster_named",
  split.by = "Condition",,raster = FALSE
) +
  ggtitle("UMAP Split by Condition")
p_umap_split
ggsave(
  "/Volumes/srivalli/dzhi/ACM_SN_R_2025/plots/umap_en_split_by_condition_singlets.png",
  p_umap_split, width = 14, height = 6
)

# Create summary table
bar_data <- en_singlets@meta.data %>%
  group_by(Condition, cluster_named) %>%
  summarise(Cell_Count = n(), .groups = "drop")
p_bar <- ggplot(bar_data, aes(x = cluster_named, y = Cell_Count, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(
    aes(label = Cell_Count),
    position = position_dodge(width = 0.9),
    vjust = -0.3,
    size = 2
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", size = 14)
  ) +
  ylab("Number of Cells") +
  xlab("Annotated Cluster") +
  ggtitle("Cluster Cell Counts by Condition")
p_bar
ggsave(
  "/Volumes/srivalli/dzhi/ACM_SN_R_2025/plots/barplot_en_cluster_by_genotype_singlets.png",
  p_bar, width = 12, height = 6
)
write.csv(bar_data, "/Volumes/srivalli/dzhi/ACM_SN_R_2025/data/counts_en_clusters_condition_singlets.csv", row.names = FALSE)

#### Integrate Singlets in Overall ####
# Extract the subcluster labels from the en_singlets object
en_subclusters <- en_singlets$cluster_named
names(en_subclusters) <- colnames(en_singlets)
integrated_singlets$en_subcluster <- NA

# Assign subcluster labels to matching barcodes
matching_barcodes <- intersect(names(en_subclusters), colnames(integrated_singlets))
integrated_singlets$en_subcluster[matching_barcodes] <- as.character(en_subclusters[matching_barcodes])
p5 <- DimPlot(integrated_singlets, group.by = "en_subcluster", raster = FALSE) +
  ggtitle("Endothelial Subclusters in Integrated Object")
p5
ggsave("/Volumes/srivalli/dzhi/ACM_SN_R_2025/plots/umap_en_subcluster_integrated.png", p5, width = 10, height = 7)

#### Save annotated Seurat object ####
saveRDS(en_singlets, "/Volumes/srivalli/dzhi/ACM_SN_R_2025/data/integrated_en_clusters_singlets.rds")
