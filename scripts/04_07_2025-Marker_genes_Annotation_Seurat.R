library(Seurat)
library(dplyr)
library(ggplot2)
library(SeuratDisk)

#-------------------------#
# Load Data
#-------------------------#
integrated <- readRDS("/home/gruengroup/srivalli/Github/ACM_SN_R_2025/data/integrated_processed.rds")

marker_ref <- read.csv("/home/gruengroup/srivalli/Github/ACM_SN_R_2025/data/marker_genes_heart.csv", stringsAsFactors = FALSE)
colnames(marker_ref)[1:3] <- c("Cell.Type", "Subtype", "Marker.Genes")
marker_ref$Marker.Genes <- toupper(marker_ref$Marker.Genes)
marker_ref$marker_list <- strsplit(marker_ref$Marker.Genes, "[ ,;]+")

#-------------------------#
# Compute Avg Expression per Cluster
#-------------------------#
DefaultAssay(integrated) <- "RNA"
integrated <- NormalizeData(integrated)
avg_exp <- AverageExpression(integrated, return.seurat = TRUE, assays = "RNA", slot = "data")
avg_matrix <- GetAssayData(avg_exp, assay = "RNA", slot = "data")
clusters <- levels(integrated$seurat_clusters)

#-------------------------#
# Scoring: Top 2 Matching Cell Types per Cluster
#-------------------------#
annotated_clusters <- data.frame(
  cluster = clusters,
  CellType_Subtype_1 = NA,
  BroadType_1 = NA,
  Score_1 = NA,
  CellType_Subtype_2 = NA,
  BroadType_2 = NA,
  Score_2 = NA
)

for (clust in clusters) {
  gene_expr <- avg_matrix[, clust]

  # Mean expression over marker genes
  scores <- sapply(seq_len(nrow(marker_ref)), function(j) {
    ref_genes <- intersect(marker_ref$marker_list[[j]], names(gene_expr))
    if (length(ref_genes) == 0) return(0)
    mean(gene_expr[ref_genes], na.rm = TRUE)
  })

  ranked <- order(scores, decreasing = TRUE)
  top1 <- ranked[1]
  top2 <- ranked[2]

  # Assign top hits
  annotated_clusters[annotated_clusters$cluster == clust, ] <- list(
    clust,
    paste0(marker_ref$Cell.Type[top1], "_", marker_ref$Subtype[top1]),
    marker_ref$Cell.Type[top1],
    scores[top1],
    paste0(marker_ref$Cell.Type[top2], "_", marker_ref$Subtype[top2]),
    marker_ref$Cell.Type[top2],
    scores[top2]
  )
}

#-------------------------#
# Add Top Matches to Seurat Object
#-------------------------#
map1 <- setNames(annotated_clusters$CellType_Subtype_1, annotated_clusters$cluster)
map2 <- setNames(annotated_clusters$CellType_Subtype_2, annotated_clusters$cluster)

integrated$marker_subtype_top1 <- dplyr::recode(as.character(integrated$seurat_clusters), !!!map1)
integrated$marker_subtype_top2 <- dplyr::recode(as.character(integrated$seurat_clusters), !!!map2)

#-------------------------#
# Create Plot Output Directory
#-------------------------#
#dir.create("/home/gruengroup/srivalli/Github/ACM_SN_R_2025/plots", showWarnings = FALSE)

#-------------------------#
# UMAP Plots
#-------------------------#
# 1. Seurat Clusters
p0 <- DimPlot(integrated, group.by = "seurat_clusters", label = TRUE, repel = TRUE) +
  ggtitle("UMAP: Seurat Clusters")
ggsave("/home/gruengroup/srivalli/Github/ACM_SN_R_2025/plots/umap_seurat_clusters.png", p0, width = 10, height = 7)

# 2. Marker Top 1
p1 <- DimPlot(integrated, group.by = "marker_subtype_top1", label = TRUE, repel = TRUE) +
  ggtitle("UMAP: Marker Gene Subtype (Top 1 Match)")
ggsave("/home/gruengroup/srivalli/Github/ACM_SN_R_2025/plots/umap_marker_subtype_top1.png", p1, width = 10, height = 7)

# 3. Marker Top 2
p2 <- DimPlot(integrated, group.by = "marker_subtype_top2", label = TRUE, repel = TRUE) +
  ggtitle("UMAP: Marker Gene Subtype (Top 2 Match)")
ggsave("/home/gruengroup/srivalli/Github/ACM_SN_R_2025/plots/umap_marker_subtype_top2.png", p2, width = 10, height = 7)

# 4. Sample ID
if ("sample_id" %in% colnames(integrated@meta.data)) {
  p3 <- DimPlot(integrated, group.by = "sample_id", label = TRUE, repel = TRUE) +
    ggtitle("UMAP: Sample ID")
  ggsave("/home/gruengroup/srivalli/Github/ACM_SN_R_2025/plots/umap_sample_id.png", p3, width = 10, height = 7)
}

# 5. Genotype
if ("genotype" %in% colnames(integrated@meta.data)) {
  p4 <- DimPlot(integrated, group.by = "genotype", label = TRUE, repel = TRUE) +
    ggtitle("UMAP: Genotype")
  ggsave("/home/gruengroup/srivalli/Github/ACM_SN_R_2025/plots/umap_genotype.png", p4, width = 10, height = 7)
}

#-------------------------#
# Save Annotated Seurat Object
#-------------------------#
saveRDS(integrated, "/home/gruengroup/srivalli/Github/ACM_SN_R_2025/data/integrated_marker_annotated.rds")

SaveH5Seurat(integrated,
  filename = "/home/gruengroup/srivalli/Github/ACM_SN_R_2025/data/integrated_marker_annotated.h5Seurat",
  overwrite = TRUE
)

Convert("/home/gruengroup/srivalli/Github/ACM_SN_R_2025/data/integrated_marker_annotated.h5Seurat",
        dest = "h5ad", overwrite = TRUE)

#-------------------------#
# Save Annotation Table
#-------------------------#
write.csv(annotated_clusters,
          "/home/gruengroup/srivalli/Github/ACM_SN_R_2025/data/marker_gene_cluster_annotations.csv",
          row.names = FALSE)