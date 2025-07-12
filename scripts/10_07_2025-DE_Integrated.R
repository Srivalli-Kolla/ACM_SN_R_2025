library(Seurat)
library(ggplot2)
library(dplyr)

# 1. Load integrated object
integrated <- readRDS("/Volumes/srivalli/dzhi/ACM_SN_R_2025/data/integrated_named_clusters_singlets.rds")

# 2. Initialize new column with current cluster names
integrated$combined_subclusters <- integrated$cluster_named
integrated$combined_subclusters <- as.character(integrated$combined_subclusters)

# 3. Define your subcluster files and the corresponding column name in each
subcluster_files <- list(
  cm = "/Volumes/srivalli/dzhi/ACM_SN_R_2025/data/integrated_cm_clusters_singlets.rds",
  fb = "/Volumes/srivalli/dzhi/ACM_SN_R_2025/data/integrated_fb_clusters_singlets.rds",
  en = "/Volumes/srivalli/dzhi/ACM_SN_R_2025/data/integrated_en_clusters_singlets.rds",
  others = "/Volumes/srivalli/dzhi/ACM_SN_R_2025/data/integrated_others_clusters_singlets.rds"
)

# 4. For each file, update the combined_subclusters column by barcode
for (file in subcluster_files) {
  obj <- readRDS(file)
  meta <- as.character(obj$cluster_named)
  names(meta) <- colnames(obj)
  common_barcodes <- intersect(names(meta), colnames(integrated))
  integrated$combined_subclusters[common_barcodes] <- meta[common_barcodes]
  rm(obj, meta)
  gc()
}

# 5. Visualize
integrated$combined_subclusters <- as.character(integrated$combined_subclusters)
p0 <- DimPlot(integrated, group.by = "combined_subclusters", raster = FALSE)
p0
ggsave("/Volumes/srivalli/dzhi/ACM_SN_R_2025/plots/umap_combined_subcluster_integrated.png", p0, width = 15, height = 7)

bar_data <- integrated@meta.data %>%
  group_by(combined_subclusters, Condition) %>%
  tally() %>%
  group_by(Condition) %>%
  mutate(Percent = 100 * n / sum(n))
bar_data
write.csv(bar_data,"/Volumes/srivalli/dzhi/ACM_SN_R_2025/data/subcluster_counts.csv", row.names = FALSE)

# Specify the unwanted cluster labels
clusters_to_remove <- c("Atrial Cardiomyocytes","Erythroid-like Cells (contaminant)","Unknown/Doublets (likely non-CM)","Erythroid-like/Contaminant or Progenitor")

# Subset the Seurat object, removing those clusters
integrated_filtered <- subset(
  integrated,
  subset = !(combined_subclusters %in% clusters_to_remove)
)
integrated_filtered$combined_subclusters <- as.character(integrated_filtered$combined_subclusters)
p1 <- DimPlot(integrated_filtered, group.by = "combined_subclusters", raster = FALSE) + 
  ggtitle('Subclusters - Filtered')
p1
ggsave("/Volumes/srivalli/dzhi/ACM_SN_R_2025/plots/umap_combined_subcluster_filtered_integrated.png", p1, width = 15, height = 7)

bar_data <- integrated_filtered@meta.data %>%
  group_by(combined_subclusters, Condition) %>%
  tally() %>%
  group_by(Condition) %>%
  mutate(Percent = 100 * n / sum(n))
bar_data
ggplot(bar_data, aes(x = combined_subclusters, y = Percent, fill = Condition)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_text(aes(label = sprintf("%.1f%%", Percent)),
            position = position_dodge(width = 0.9),
            vjust = -0.2, size = 2)+
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(
    x = "Subcluster",
    y = "Percent of Cells",
    title = "Cell Type Composition by Condition"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
write.csv(bar_data,"/Volumes/srivalli/dzhi/ACM_SN_R_2025/data/subcluster_filtered_counts.csv", row.names = FALSE)

#### Major grouping #####
integrated_filtered$MajorGroup <- NA

integrated_filtered$MajorGroup[grepl("Cardiomyocyte", integrated_filtered$combined_subclusters, ignore.case = TRUE)] <- "Cardiomyocytes"
integrated_filtered$MajorGroup[grepl("Fibroblast", integrated_filtered$combined_subclusters, ignore.case = TRUE)] <- "Fibroblasts"
integrated_filtered$MajorGroup[grepl("Endothelial", integrated_filtered$combined_subclusters, ignore.case = TRUE)] <- "Endothelial"
integrated_filtered$MajorGroup[grepl("Endocardial", integrated_filtered$combined_subclusters, ignore.case = TRUE)] <- "Endocardial"
integrated_filtered$MajorGroup[grepl("Macrophage|Monocyte", integrated_filtered$combined_subclusters, ignore.case = TRUE)] <- "Macrophages"
integrated_filtered$MajorGroup[grepl("Mast Cell", integrated_filtered$combined_subclusters, ignore.case = TRUE)] <- "Mast Cells"
integrated_filtered$MajorGroup[grepl("T Lymphocyte|CD4|CD8", integrated_filtered$combined_subclusters, ignore.case = TRUE)] <- "T Lymphocytes"
integrated_filtered$MajorGroup[grepl("B Lymphocyte", integrated_filtered$combined_subclusters, ignore.case = TRUE)] <- "B Lymphocytes"
integrated_filtered$MajorGroup[grepl("Pericyte", integrated_filtered$combined_subclusters, ignore.case = TRUE)] <- "Pericytes"
integrated_filtered$MajorGroup[grepl("Neuronal|Glial", integrated_filtered$combined_subclusters, ignore.case = TRUE)] <- "Neuronal"
integrated_filtered$MajorGroup[grepl("Smooth Muscle", integrated_filtered$combined_subclusters, ignore.case = TRUE)] <- "Smooth Muscle Cells"
integrated_filtered$MajorGroup[grepl("Proliferating|Proliferative", integrated_filtered$combined_subclusters, ignore.case = TRUE)] <- "Proliferating/Dividing Cells"
integrated_filtered$MajorGroup[integrated_filtered$combined_subclusters == "Purkinje/Conduction System-like Cells"] <- "Cardiomyocytes"
integrated_filtered$MajorGroup[integrated_filtered$combined_subclusters == "Purkinje/Conduction System-like CMs"] <- "Cardiomyocytes"

meta_integrated_filtered <- integrated_filtered@meta.data

summary_counts <- meta_integrated_filtered %>%
  group_by(MajorGroup, Condition) %>%
  summarise(total_cells = n(), .groups = 'drop')

p <- ggplot(summary_counts, aes(x = MajorGroup, y = total_cells, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(aes(label = total_cells), 
            position = position_dodge(width = 0.9), 
            vjust = -0.2, size = 2) +
  theme_minimal() +
  labs(
    x = "Major Cell Group",
    y = "Cell Count",
    title = "Cell Counts per Major Cell Group and Condition _ After Filtering"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/Volumes/srivalli/dzhi/ACM_SN_R_2025/plots/barplot_after_filtering_major.png", plot = p, width = 8, height = 5)

bar_data2 <- meta_df %>%
  group_by(MajorGroup, combined_subclusters, Condition) %>%
  tally(name = "n")
bar_data2 <- bar_data2 %>%
  filter(!is.na(MajorGroup))
major_groups <- unique(bar_data2$MajorGroup)

# Loop through each major group and plot/save
for (grp in major_groups) {
  sub_data <- bar_data2 %>% filter(MajorGroup == grp)
  
  p <- ggplot(sub_data, aes(x = combined_subclusters, y = n, fill = Condition)) +
    geom_col(position = position_dodge(width = 0.9)) +
    geom_text(
      aes(label = n),
      position = position_dodge(width = 0.9),
      vjust = -0.2, size = 2
    )+
    labs(
      x = "Subcluster",
      y = "Cell Count",
      title = paste(grp, "- Cell Counts per Subcluster and Condition")
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold")
    )
  
  # Safe filename
  filename <- paste0("plot_", gsub("[^A-Za-z0-9]", "_", grp), ".png")
  fullpath <- file.path("/Volumes/srivalli/dzhi/ACM_SN_R_2025/plots/", filename)
  
  ggsave(filename = fullpath, plot = p, width = 8, height = 5)
}
#### Save annotated Seurat object ####
saveRDS(integrated_filtered, "/Volumes/srivalli/dzhi/ACM_SN_R_2025/data/integrated_filtered.rds")
