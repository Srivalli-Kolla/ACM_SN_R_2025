library(Seurat)
library(dplyr)

# Load data
integrated <- readRDS("/Volumes/srivalli/dzhi/ACM_SN_R_2025/data/integrated_processed.rds")

marker_ref <- read.csv("/Volumes/srivalli/dzhi/ACM_SN_R_2025/data/marker_genes_heart.csv", stringsAsFactors = FALSE)
colnames(marker_ref)[1:3] <- c("Cell.Type", "Subtype", "Marker.Genes")

# Split marker genes into lists exactly as they appear (no uppercasing)
marker_ref$marker_list <- strsplit(marker_ref$Marker.Genes, "[ ,;]+")

# Find marker genes present in the Seurat object (exact match)
all_marker_genes <- unique(unlist(marker_ref$marker_list))
genes_in_data <- intersect(all_marker_genes, rownames(integrated))
cat("Marker genes matched in data:", length(genes_in_data), "\n")

# Subset Seurat object to marker genes found
subset_integrated <- subset(integrated, features = genes_in_data)

# Calculate average expression per cluster
DefaultAssay(subset_integrated) <- "RNA"
avg_exp <- AverageExpression(subset_integrated, return.seurat = FALSE, slot = "data")
avg_matrix <- avg_exp$RNA  
clusters <- colnames(avg_matrix)

# Initialize score matrix: rows=clusters, cols=marker sets
score_matrix <- matrix(0, nrow = length(clusters), ncol = nrow(marker_ref))
rownames(score_matrix) <- clusters
colnames(score_matrix) <- paste0(marker_ref$Cell.Type, "_", marker_ref$Subtype)

# Score clusters against each marker set by averaging expression of marker genes
for (clust in clusters) {
  gene_expr <- avg_matrix[, clust]
  
  scores <- sapply(seq_len(nrow(marker_ref)), function(j) {
    ref_genes <- intersect(marker_ref$marker_list[[j]], rownames(avg_matrix))
    if (length(ref_genes) == 0) return(0)
    mean(gene_expr[ref_genes], na.rm = TRUE)
  })
  
  score_matrix[clust, ] <- scores
}

# Annotate clusters with top 3 marker sets by score
annotated_clusters <- data.frame(cluster = clusters)

top_matches <- apply(score_matrix, 1, function(x) {
  ord <- order(x, decreasing = TRUE)
  top1 <- ord[1]
  top2 <- ord[2]
  top3 <- ord[3]
  c(
    paste0(marker_ref$Cell.Type[top1], "_", marker_ref$Subtype[top1]),
    x[top1],
    paste0(marker_ref$Cell.Type[top2], "_", marker_ref$Subtype[top2]),
    x[top2],
    paste0(marker_ref$Cell.Type[top3], "_", marker_ref$Subtype[top3]),
    x[top3]
  )
})

annotated_clusters$CellType_Subtype_1 <- top_matches[1, ]
annotated_clusters$Score_1 <- as.numeric(top_matches[2, ])
annotated_clusters$CellType_Subtype_2 <- top_matches[3, ]
annotated_clusters$Score_2 <- as.numeric(top_matches[4, ])
annotated_clusters$CellType_Subtype_3 <- top_matches[5, ]
annotated_clusters$Score_3 <- as.numeric(top_matches[6, ])

# Map annotations back to Seurat object metadata
map1 <- setNames(annotated_clusters$CellType_Subtype_1, annotated_clusters$cluster)
map2 <- setNames(annotated_clusters$CellType_Subtype_2, annotated_clusters$cluster)
map3 <- setNames(annotated_clusters$CellType_Subtype_3, annotated_clusters$cluster)
# Ensure cluster IDs in Seurat metadata are character
integrated$seurat_clusters <- as.character(integrated$seurat_clusters)

# Make sure the names in the map are character as well (usually they are)
names(map1) <- sub("^g", "", names(map1))
names(map2) <- sub("^g", "", names(map2))
names(map3) <- sub("^g", "", names(map3))

# Now recode the cluster IDs to your annotation names
integrated$marker_subtype_top1 <- dplyr::recode(integrated$seurat_clusters, !!!map1, .default = NA_character_)
integrated$marker_subtype_top2 <- dplyr::recode(integrated$seurat_clusters, !!!map2, .default = NA_character_)
integrated$marker_subtype_top3 <- dplyr::recode(integrated$seurat_clusters, !!!map3, .default = NA_character_)

# Check that the mapping worked:
table(integrated$marker_subtype_top1, useNA = "ifany")
p1 <- DimPlot(integrated, group.by = "marker_subtype_top1", repel = TRUE) +
  ggtitle("UMAP: Marker Gene Subtype (Top 1 Match)")
p1
ggsave("/Volumes/srivalli/dzhi/ACM_SN_R_2025/plots/umap_marker_subtype_top1.png", p1, width = 10, height = 7)


# Save results 
write.csv(score_matrix, "/Volumes/srivalli/dzhi/ACM_SN_R_2025/data/marker_score_matrix.csv", row.names = TRUE)
write.csv(annotated_clusters, "/Volumes/srivalli/dzhi/ACM_SN_R_2025/data/marker_gene_cluster_annotations.csv", row.names = FALSE)
saveRDS(integrated, "/Volumes/srivalli/dzhi/ACM_SN_R_2025/data/integrated_marker_annotated.rds")
