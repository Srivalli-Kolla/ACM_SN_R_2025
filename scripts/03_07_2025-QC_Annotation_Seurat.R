#### Load required libraries ####
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(tibble)
library(reshape2)

#### 1. Define sample names and Add Metadata ####
sample_names <- c(
  "B3_Lib2", "B4_Lib2", "B5_Lib2", "B6_Lib2", "B7_Lib2",
  "B8_Lib2", "B9_Lib2", "B10_Lib2", "B11_Lib2"
)
data_paths <- paste0("/Volumes/srivalli/dzhi/ACM_sn_2025/data/Library2/", sample_names, "/filtered_feature_bc_matrix/")

#### 2. Load and QC-filter Seurat objects ####
seurat_list <- list()
for (i in seq_along(sample_names)) {
  counts <- Read10X(data.dir = data_paths[i])
  obj <- CreateSeuratObject(counts = counts, project = sample_names[i])
  obj$sample <- sample_names[i]
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
  
  # Filter low-quality cells
  obj <- subset(obj, subset = nFeature_RNA > 200 & percent.mt < 10)
  seurat_list[[sample_names[i]]] <- obj
}

#### 3. Merge for visualization & Add Metadata ####
combined_filtered <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list), project = "after_filter")
combined_filtered[["percent.mt"]] <- PercentageFeatureSet(combined_filtered, pattern = "^mt-")
combined_filtered$sample_base <- gsub("_Lib2$", "", combined_filtered$sample)

meta_data <- read.csv("/Volumes/srivalli/dzhi/ACM_sn_2025/data/acm_sn_metadata.csv", stringsAsFactors = FALSE)

seurat_meta <- combined_filtered@meta.data %>%
  rownames_to_column("cell") %>%
  left_join(meta_data, by = c("sample_base" = "sample")) %>%
  column_to_rownames("cell")
combined_filtered@meta.data <- seurat_meta

#### 4. Violin QC Plot ####
p_violin <- VlnPlot(
  combined_filtered,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "Sample_ID",
  pt.size = 0.1,
  ncol = 3
)
p_violin
ggsave("./Documents/Github/ACM_SN_R_2025/plots/qc_violin_after_filtered.png", p_violin, width = 12, height = 6)

#### 5. Barplot: Cell counts after filtering ####
cell_counts <- as.data.frame(table(combined_filtered$Sample_ID))
colnames(cell_counts) <- c("Sample_ID", "CellCount")

p_bar <- ggplot(cell_counts, aes(x = Sample_ID, y = CellCount)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = CellCount), vjust = -0.3, size = 3) +
  theme_minimal() +
  labs(title = "Cell Counts After Filtering") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_bar
ggsave("./Documents/Github/ACM_SN_R_2025/plots/cell_counts_barplot_filtered.png", p_bar, width = 10, height = 6)

#### 6. Normalize and Identify Variable Features ####
seurat_list <- lapply(seurat_list, function(obj) {
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  return(obj)
})

#### 5. Integration ####
features <- SelectIntegrationFeatures(object.list = seurat_list)
anchors <- FindIntegrationAnchors(object.list = seurat_list, anchor.features = features,dims = 1:30)
integrated <- IntegrateData(anchorset = anchors, dims = 1:30)

#### 6. Scaling, PCA, UMAP, Clustering ####
DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated)
integrated <- RunPCA(integrated, npcs = 30)

# Save Elbow plot
elbow <- ElbowPlot(integrated, ndims = 30)
ggsave("./Documents/Github/ACM_SN_R_2025/plots/elbow_plot.png", elbow)

integrated <- FindNeighbors(integrated, dims = 1:30)
integrated <- FindClusters(integrated, resolution = 0.3)
integrated <- RunUMAP(integrated, dims = 1:30)

#### 7. UMAP Visualizations ####
p1 <- DimPlot(integrated, reduction = "umap", group.by = "sample") + ggtitle("By Sample")
p1
ggsave("./Documents/Github/ACM_SN_R_2025/plots/umap_by_sample.png", p1)

p2 <- DimPlot(integrated, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + ggtitle("By Cluster")
p2
ggsave("./Documents/Github/ACM_SN_R_2025/plots/umap_by_cluster.png", p2)

#### 8. DE genes ####
cluster_markers <- FindAllMarkers(
  integrated,
  only.pos = TRUE,
  test.use = "wilcox",
  min.pct = 0.25,
  logfc.threshold = 0.25
)
top5 <- cluster_markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
top10 <- cluster_markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)

write.csv(top10, "../Documents/Github/ACM_SN_R_2025/de_top_10_clusters.csv", row.names = FALSE)

DotPlot(integrated, features = unique(top5$gene)) + RotatedAxis()
ggsave("../Documents/Github/ACM_SN_R_2025/dotplot_top5_markers.png")

#### 9. Annotation ####
##### Marker Gene Annotation #####
marker_ref <- read.csv("../Documents/Github/ACM_SN_R_2025/marker_genes_heart- marker_genes.csv", stringsAsFactors = FALSE)
marker_ref$Marker.Genes <- toupper(marker_ref$Marker.Genes)
marker_ref$marker_list <- strsplit(marker_ref$Marker.Genes, "[ ,;]+")
top_markers <- cluster_markers %>%
  group_by(cluster) %>%
  top_n(30, avg_log2FC)
top_markers$gene <- toupper(top_markers$gene)

annotated_clusters <- data.frame(
  cluster = unique(top_markers$cluster),
  CellType_Subtype = NA,
  BroadType = NA
)

for (i in seq_along(annotated_clusters$cluster)) {
  clust <- annotated_clusters$cluster[i]
  clust_genes <- unique(top_markers$gene[top_markers$cluster == clust])
  
  best_match <- NULL
  best_score <- 0
  
  for (j in 1:nrow(marker_ref)) {
    ref_genes <- marker_ref$marker_list[[j]]
    score <- length(intersect(clust_genes, ref_genes))
    
    if (score > best_score) {
      best_score <- score
      best_match <- j
    }
  }
  
  if (!is.null(best_match)) {
    annotated_clusters$CellType_Subtype[i] <- paste0(
      marker_ref$Cell.Type[best_match], "_", marker_ref$Subtype...State[best_match]
    )
    annotated_clusters$BroadType[i] <- marker_ref$Cell.Type[best_match]
  }
}

ct_map <- setNames(annotated_clusters$CellType_Subtype, annotated_clusters$cluster)
broad_map <- setNames(annotated_clusters$BroadType, annotated_clusters$cluster)

integrated$celltype_marker_subtype <- plyr::mapvalues(
  x = as.character(integrated$seurat_clusters),
  from = names(ct_map),
  to = ct_map
)

integrated$celltype_marker_broad <- plyr::mapvalues(
  x = as.character(integrated$seurat_clusters),
  from = names(broad_map),
  to = broad_map
)

DimPlot(integrated, group.by = "celltype_marker_subtype", label = TRUE, repel = TRUE) +
  ggtitle("Clusters Annotated by CellType_Subtype")
ggsave("../Documents/Github/ACM_SN_R_2025/umap_celltype_subtype.png", width = 10, height = 7)

DimPlot(integrated, group.by = "celltype_marker_broad", label = TRUE, repel = TRUE) +
  ggtitle("Clusters Annotated by Broad Cell Type")
ggsave("../Documents/Github/ACM_SN_R_2025/umap_celltype_broad.png", width = 10, height = 7)

