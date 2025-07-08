library(Seurat)
library(DoubletFinder)
library(tidyverse)


# Load your Seurat object
integrated <- readRDS("/Volumes/srivalli/dzhi/ACM_SN_R_2025/data/integrated_named_clusters.rds")
DefaultAssay(integrated) <- "RNA"

# Split by sample
sample_list <- SplitObject(integrated, split.by = "Sample_ID")

# Run DoubletFinder per sample
for (i in seq_along(sample_list)) {
  sample_obj <- sample_list[[i]]
  cat("Processing:", unique(sample_obj$Sample_ID), "with", ncol(sample_obj), "cells\n")
  if (ncol(sample_obj) < 100) next
  
  # Preprocessing
  sample_obj <- NormalizeData(sample_obj)
  sample_obj <- FindVariableFeatures(sample_obj)
  sample_obj <- ScaleData(sample_obj)
  sample_obj <- RunPCA(sample_obj, npcs = 20)
  
  # Remove any old DoubletFinder columns
  df_cols <- grep("^(pANN|DF.classifications)", colnames(sample_obj@meta.data), value = TRUE)
  if (length(df_cols) > 0) {
    sample_obj@meta.data[, df_cols] <- NULL
  }
  
  # Parameter sweep
  sweep.res.list <- paramSweep(sample_obj, PCs = 1:20, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  optimal.pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  nExp <- round(0.05 * ncol(sample_obj))
  if (nExp < 1) next
  
  tryCatch({
    sample_obj <- doubletFinder(
      sample_obj,
      PCs = 1:20,
      pN = 0.25,
      pK = optimal.pK,
      nExp = nExp,
      reuse.pANN = FALSE, 
      sct = FALSE
    )
    # Check for data frames after DoubletFinder
    meta_types <- sapply(sample_obj@meta.data, class)
    if (any(meta_types == "data.frame")) {
      print(meta_types[meta_types == "data.frame"])
      stop("Data frame found in meta.data after DoubletFinder!")
    }
    sample_list[[i]] <- sample_obj
    cat("DoubletFinder classifications added for sample:", unique(sample_obj$Sample_ID), "\n")
  }, error = function(e) {
    cat("DoubletFinder failed for sample:", unique(sample_obj$Sample_ID), "\nError:", e$message, "\n")
  })
}

# Combine DoubletFinder results back to integrated object
doublet_calls <- unlist(lapply(sample_list, function(obj) {
  df_col <- grep("^DF.classifications", colnames(obj@meta.data), value = TRUE)
  if (length(df_col) == 0) return(NULL)
  calls <- obj@meta.data[[df_col[length(df_col)]]]
  names(calls) <- rownames(obj@meta.data)
  return(calls)
}))
names(doublet_calls) <- sub("^.*\\.", "", names(doublet_calls))

integrated$DoubletFinder <- NA
matching_cells <- intersect(names(doublet_calls), colnames(integrated))
integrated$DoubletFinder[matching_cells] <- doublet_calls[matching_cells]

## Visualization
DefaultAssay(integrated) <- "RNA"
if (!"umap" %in% names(integrated@reductions)) {
  integrated <- RunPCA(integrated, npcs = 20)
  integrated <- RunUMAP(integrated, dims = 1:20)
}

# UMAP plot colored by DoubletFinder classification
umap_plot <- DimPlot(
  integrated,
  group.by = "DoubletFinder",
  reduction = "umap",
  pt.size = 0.5
)
umap_plot
ggsave("/Volumes/srivalli/dzhi/ACM_SN_R_2025/plots/umap_doubletfinder.png", plot = umap_plot, width = 7, height = 5, dpi = 300)

# Doublets per Sample 
df_sample <- integrated@meta.data %>%
  group_by(Sample_ID, DoubletFinder) %>%
  tally() %>%
  filter(DoubletFinder == "Doublet")

# Identify the sample with the highest doublet count
max_sample_row <- df_sample[which.max(df_sample$n), ]
sample_to_annotate <- max_sample_row$Sample_ID
y_to_annotate <- max_sample_row$n + 2  # Adjust as needed for spacing

# Bar plot: Doublets per Sample, with annotation
p_sample <- ggplot(df_sample, aes(x = Sample_ID, y = n, fill = Sample_ID)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = n), vjust = -0.5) +  # Add count labels above bars
  labs(title = "Number of Doublets per Sample",
       x = "Sample",
       y = "Doublet Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none") +
  annotate("text", x = sample_to_annotate, y = y_to_annotate,
           , label= NA,size = 5, fontface = "bold")
p_sample
ggsave("/Volumes/srivalli/dzhi/ACM_SN_R_2025/plots/doublets_per_sample.png", plot = p_sample, width = 7, height = 5, dpi = 300)

# Doublets per Cluster
df_cluster <- integrated@meta.data %>%
  group_by(seurat_clusters, DoubletFinder) %>%
  tally() %>%
  filter(DoubletFinder == "Doublet")

# Identify the cluster with the highest doublet count
max_cluster_row <- df_cluster[which.max(df_cluster$n), ]
cluster_to_annotate <- max_cluster_row$seurat_clusters
y_cluster_annotate <- max_cluster_row$n + 2

# Bar plot: Doublets per Cluster, with annotation
p_cluster <- ggplot(df_cluster, aes(x = seurat_clusters, y = n, fill = seurat_clusters)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = n), vjust = -0.5) +
  labs(title = "Number of Doublets per Cluster",
       x = "Cluster",
       y = "Doublet Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none") +
  annotate("text", x = cluster_to_annotate, y = y_cluster_annotate,
           label = NA, size = 5, fontface = "bold")
p_cluster
ggsave("/Volumes/srivalli/dzhi/ACM_SN_R_2025/plots/doublets_per_cluster.png", plot = p_cluster, width = 7, height = 5, dpi = 300)

table(integrated$cluster_named,integrated$DoubletFinder,integrated$Condition)

# Save the object
saveRDS(integrated, "/Volumes/srivalli/dzhi/ACM_SN_R_2025/data/integrated_with_doubletfinder.rds")
