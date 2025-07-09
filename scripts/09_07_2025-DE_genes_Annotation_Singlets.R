library(Seurat)
library(dplyr)
library(ggplot2)
library(SeuratDisk)

#### Load processed integrated object ####
integrated <- readRDS("/Volumes/srivalli/dzhi/ACM_SN_R_2025/data/integrated_with_doubletfinder.rds")
DefaultAssay(integrated) <- "integrated"
singlet_col <- grep("DoubletFinder", colnames(integrated@meta.data), value = TRUE)[1]
singlets_logical <- integrated@meta.data[[singlet_col]] == "Singlet"
names(singlets_logical) <- rownames(integrated@meta.data)
singlet_cells <- names(singlets_logical)[singlets_logical]
singlet_cells <- intersect(singlet_cells, colnames(integrated)); integrated@graphs <- list()
integrated@neighbors <- list()
DefaultAssay(integrated) <- "integrated"  

integrated_singlets <- subset(integrated, cells = singlet_cells)
remove(integrated)

#### Perform DE analysis ####
Idents(integrated_singlets) <- "seurat_clusters"
markers <- FindAllMarkers(
  integrated_singlets,
  only.pos = TRUE,
  test.use = "wilcox",
  min.pct = 0.25,
  logfc.threshold = 0.25
)

top5 <- markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)

# Save DE marker results
write.csv(markers, "/Volumes/srivalli/dzhi/ACM_SN_R_2025/data/de_all_markers_singlets.csv", row.names = FALSE)
write.csv(top10, "/Volumes/srivalli/dzhi/ACM_SN_R_2025/data/de_top_10_clusters_singlets.csv", row.names = FALSE)

#### Manually set human-readable cluster names ####
cluster_names <- c(
  "Cardiac Fibroblasts (ECM remodeling)",        # 0
  "Capillary Endothelial Cells",                 # 1
  "Ventricular Cardiomyocytes (contractile)",    # 2
  "Cardiac Mast Cells",                          # 3
  "Epicardial Cells",                            # 4
  "Ventricular Cardiomyocytes (Metabolic/Structural)",      # 5
  "Pericytes (microvascular mural cells)",       # 6
  "Lymphatic Endothelial Cells",                 # 7
  "Endocardial Cells (conduction system)",       # 8
  "Cardiac Macrophages/Monocytes",               # 9
  "Adventitial Fibroblasts/Pericytes",           # 10
  "B Lymphocytes (mature B cells)",              # 11
  "T Lymphocytes (CD4+/CD8+ mix)",               # 12
  "Vascular Smooth Muscle Cells (contractile)",  # 13
  "Cycling/Proliferating Cells (mitotic)",       # 14
  "Schwann Cells (neural crest-derived glia)",   # 15
  "Epicardial Mesothelial Cells (Wt1+)"          # 16
)

# Map readable cluster names
integrated_singlets$cluster_named <- factor(
  cluster_names[as.numeric(as.character(integrated_singlets$seurat_clusters)) + 1],
  levels = cluster_names
)

# UMAPS
p_dot <- DotPlot(integrated_singlets, features = unique(top5$gene)) + RotatedAxis()
p_dot
ggsave("/Volumes/srivalli/dzhi/ACM_SN_R_2025/plots/dotplot_top5_markers_singlets.png", p_dot, width = 20, height = 10)

# Create summary table
bar_data <- integrated_singlets@meta.data %>%
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
  "/Volumes/srivalli/dzhi/ACM_SN_R_2025/plots/barplot_cluster_by_genotype.png",
  p_bar, width = 12, height = 6
)
write.csv(bar_data, "/Volumes/srivalli/dzhi/ACM_SN_R_2025/data/counts_clusters_condition_singlets.csv", row.names = FALSE)

bar_data_prop <- bar_data %>%
  group_by(Condition) %>%
  mutate(Proportion = Cell_Count / sum(Cell_Count))
p_bar_prop <- ggplot(bar_data_prop, aes(x = cluster_named, y = Proportion, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(
    aes(label = round(Proportion, 1)),  # Round to 2 decimal places
    position = position_dodge(width = 0.9),
    vjust = -0.3,
    size = 2
  ) +
  theme_minimal() +
  ylab("Proportion of Cells") +
  xlab("Annotated Cluster") +
  ggtitle("Cluster Proportions by Condition") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_bar_prop
ggsave(
  "/Volumes/srivalli/dzhi/ACM_SN_R_2025/plots/barplot_cluster_proportion_by_condition_singlets.png",
  p_bar_prop, width = 12, height = 6
)
