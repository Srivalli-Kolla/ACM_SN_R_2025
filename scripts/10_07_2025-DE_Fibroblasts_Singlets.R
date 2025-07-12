library(Seurat)
library(dplyr)
library(ggplot2)
library(SeuratDisk)

#### Load processed integrated object and Select Fibroblasts only ####
integrated_singlets <- readRDS("/Volumes/srivalli/dzhi/ACM_SN_R_2025/data/integrated_named_clusters_singlets.rds")
DefaultAssay(integrated_singlets) <- "integrated"
fb_labels <- c(
    "Cardiac Fibroblasts (ECM remodeling)",
    "Adventitial Fibroblasts/Pericytes",
    "Pericytes (microvascular mural cells)")
fb_cells <- colnames(integrated_singlets)[integrated_singlets$cluster_named %in% fb_labels ]
fb_singlets <- subset(integrated_singlets, cells = fb_cells)

#### Reclustering Fibroblasts ####
DefaultAssay(fb_singlets) <- "RNA"
fb_singlets <- NormalizeData(fb_singlets)
fb_singlets <- FindVariableFeatures(fb_singlets)
fb_singlets <- ScaleData(fb_singlets)
fb_singlets <- RunPCA(fb_singlets)
fb_singlets <- RunUMAP(fb_singlets, dims = 1:10)
fb_singlets <- FindNeighbors(fb_singlets, dims = 1:10)
fb_singlets <- FindClusters(fb_singlets, resolution = 0.3) 

p0 <- DimPlot(fb_singlets, label = TRUE, group.by = "seurat_clusters") + ggtitle("Fibroblasts Subclusters")
p0
ggsave("/Volumes/srivalli/dzhi/ACM_SN_R_2025/plots/umap_fb_clusters.png", p0, width = 10, height = 7)

#### DE analysis ####
fb_markers <- FindAllMarkers(
  fb_singlets,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

top5 <- fb_markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
top10 <- fb_markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)

# Save DE marker results
write.csv(fb_markers, "/Volumes/srivalli/dzhi/ACM_SN_R_2025/data/de_fb_markers_singlets.csv", row.names = FALSE)
write.csv(top10, "/Volumes/srivalli/dzhi/ACM_SN_R_2025/data/de_top_10_fb_singlets.csv", row.names = FALSE)

#### Manually set human-readable cluster names ####
cluster_names <- c(
"Ventricular Cardiomyocytes (contractile)",
"Arterial/Vascular Endothelial Cells",
"Cardiac Fibroblasts (ECM remodeling)",
"Cardiac Fibroblasts (matrix remodeling)",
"Adventitial Fibroblasts/Pericytes",
"Microvascular Endothelial Cells",
"Erythroid-like Cells (contaminant)",
"Pericytes (microvascular mural cells)",
"Cardiac Macrophages/Monocytes"
)

# Map readable cluster names
fb_singlets$cluster_named <- factor(
  cluster_names[as.numeric(as.character(fb_singlets$seurat_clusters)) + 1],
  levels = cluster_names
)
# UMAPS
p_dot <- DotPlot(fb_singlets, features = unique(top5$gene)) + RotatedAxis()
p_dot
ggsave("/Volumes/srivalli/dzhi/ACM_SN_R_2025/plots/dotplot_top5_fb_markers_singlets.png", p_dot, width = 20, height = 10)

# 2. Named clusters
p1 <- DimPlot(fb_singlets, group.by = "cluster_named", repel = TRUE,raster = FALSE) +
  ggtitle("UMAP by Annotated Cluster - Fibroblasts")
p1
ggsave("/Volumes/srivalli/dzhi/ACM_SN_R_2025/plots/umap_fb_named_clusters_singlets.png", p1, width = 10, height = 7)

# 3. Sample-wise UMAP
if ("Sample_ID" %in% colnames(fb_singlets@meta.data)) {
  p2 <- DimPlot(fb_singlets, group.by = "Sample_ID", label = FALSE,raster = FALSE) +
    ggtitle("UMAP by Sample")
  p2
  ggsave("/Volumes/srivalli/dzhi/ACM_SN_R_2025/plots/umap_fb_by_sample_singlets.png", p2, width = 10, height = 7)
}

# 4. Genotype-wise UMAP
if ("Genotype" %in% colnames(fb_singlets@meta.data)) {
  p3 <- DimPlot(fb_singlets, group.by = "Genotype", label = FALSE,raster = FALSE) +
    ggtitle("UMAP by Genotype")
  p3
  ggsave("/Volumes/srivalli/dzhi/ACM_SN_R_2025/plots/umap_fb_by_genotype_singlets.png", p3, width = 10, height = 7)
}

# 5. Condition-wise UMAP
if ("Condition" %in% colnames(fb_singlets@meta.data)) {
  p4 <- DimPlot(fb_singlets, group.by = "Condition", label = FALSE,raster = FALSE) +
    ggtitle("UMAP by Condition")
  p4
  ggsave("/Volumes/srivalli/dzhi/ACM_SN_R_2025/plots/umap_by_fb_Condition_singlets.png", p4, width = 10, height = 7)
}

# Split UMAP by genotype
p_umap_split <- DimPlot(
  fb_singlets,
  group.by = "cluster_named",
  split.by = "Condition",,raster = FALSE
) +
  ggtitle("UMAP Split by Condition")
p_umap_split
ggsave(
  "/Volumes/srivalli/dzhi/ACM_SN_R_2025/plots/umap_fb_split_by_condition_singlets.png",
  p_umap_split, width = 14, height = 6
)

# Create summary table
bar_data <- fb_singlets@meta.data %>%
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
  "/Volumes/srivalli/dzhi/ACM_SN_R_2025/plots/barplot_fb_cluster_by_genotype_singlets.png",
  p_bar, width = 12, height = 6
)
write.csv(bar_data, "/Volumes/srivalli/dzhi/ACM_SN_R_2025/data/counts_fb_clusters_condition_singlets.csv", row.names = FALSE)

#### Integrate Singlets in Overall ####
# Extract the subcluster labels from the fb_singlets object
fb_subclusters <- fb_singlets$cluster_named
names(fb_subclusters) <- colnames(fb_singlets)
integrated_singlets$fb_subcluster <- NA

# Assign subcluster labels to matching barcodes
matching_barcodes <- intersect(names(fb_subclusters), colnames(integrated_singlets))
integrated_singlets$fb_subcluster[matching_barcodes] <- as.character(fb_subclusters[matching_barcodes])
p5 <- DimPlot(integrated_singlets, group.by = "fb_subcluster", raster = FALSE) +
  ggtitle("Fibroblasts Subclusters in Integrated Object")
p5
ggsave("/Volumes/srivalli/dzhi/ACM_SN_R_2025/plots/umap_fb_subcluster_integrated.png", p5, width = 10, height = 7)

#### Save annotated Seurat object ####
saveRDS(fb_singlets, "/Volumes/srivalli/dzhi/ACM_SN_R_2025/data/integrated_fb_clusters_singlets.rds")
