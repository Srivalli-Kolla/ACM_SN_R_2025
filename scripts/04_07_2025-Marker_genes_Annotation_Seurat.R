library(Seurat)
library(dplyr)
library(plyr)
library(ggplot2)
library(SeuratDisk)

#### Load integrated object ####
integrated <- readRDS("/home/gruengroup/srivalli/Github/ACM_SN_R_2025/data/integrated_processed.rds")

#### Load marker gene reference ####
marker_ref <- read.csv("/home/gruengroup/srivalli/Github/ACM_SN_R_2025/data/marker_genes_heart_copy - marker_genes.csv", stringsAsFactors = FALSE)
marker_ref$Marker.Genes <- toupper(marker_ref$Marker.Genes)
marker_ref$marker_list <- strsplit(marker_ref$Marker.Genes, "[ ,;]+")

#### Load DEGs or use top N from scratch ####
markers <- readRDS("/home/gruengroup/srivalli/Github/ACM_SN_R_2025/data/de_all_markers.rds")
top_markers <- markers %>% group_by(cluster) %>% top_n(30, avg_log2FC)
top_markers$gene <- toupper(top_markers$gene)

#### Matching markers to annotation ####
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

# Map annotations
ct_map <- setNames(annotated_clusters$CellType_Subtype, annotated_clusters$cluster)
broad_map <- setNames(annotated_clusters$BroadType, annotated_clusters$cluster)

integrated$marker_subtype <- plyr::mapvalues(
  as.character(integrated$seurat_clusters),
  names(ct_map),
  ct_map
)
integrated$marker_broad <- plyr::mapvalues(
  as.character(integrated$seurat_clusters),
  names(broad_map),
  broad_map
)

# UMAPs with annotations
p1 <- DimPlot(integrated, group.by = "marker_subtype", label = TRUE, repel = TRUE) +
  ggtitle("Clusters Annotated by CellType - Subtype")
ggsave("/home/gruengroup/srivalli/Github/ACM_SN_R_2025/plots/umap_marker_subtype.png", p1, width = 10, height = 7)

p2 <- DimPlot(integrated, group.by = "marker_broad", label = TRUE, repel = TRUE) +
  ggtitle("Clusters Annotated by Broad Cell Type")
ggsave("/home/gruengroup/srivalli/Github/ACM_SN_R_2025/plots/umap_marker_broad.png", p2, width = 10, height = 7)

# Save annotated version
SaveH5Seurat(integrated, filename = "/home/gruengroup/srivalli/Github/ACM_SN_R_2025/plots/integrated_annotated.h5Seurat", overwrite = TRUE)
Convert("/home/gruengroup/srivalli/Github/ACM_SN_R_2025/data/integrated_annotated.h5Seurat", dest = "h5ad", overwrite = TRUE)
saveRDS(integrated, file = "/home/gruengroup/srivalli/Github/ACM_SN_R_2025/data/integrated_annotated.rds")