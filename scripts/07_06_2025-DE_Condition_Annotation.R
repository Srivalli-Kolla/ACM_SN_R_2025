library(Seurat)
library(dplyr)
library(ggplot2)
library(SeuratDisk)

#### Load processed integrated object ####
integrated <- readRDS("/Volumes/srivalli/dzhi/ACM_SN_R_2025/data/integrated_processed.rds")
# Metadata table
meta_info <- data.frame(
  sample = paste0("B", 3:11),
  Sample_Name = c(
    "20241021_AC_Pkp2-MCMV_Mouse_Single_Heart_1",
    "20241021_AC_Pkp2-MCMV_Mouse_Single_Heart_2",
    "20241021_AC_Pkp2-MCMV_Mouse_Single_Heart_3",
    "20241021_AC_Pkp2-MCMV_Mouse_Single_Heart_4",
    "20241021_AC_Pkp2-MCMV_Mouse_Single_Heart_5",
    "20241021_AC_Pkp2-MCMV_Mouse_Single_Heart_6",
    "20241021_AC_Pkp2-MCMV_Mouse_Single_Heart_7",
    "20241021_AC_Pkp2-MCMV_Mouse_Single_Heart_8",
    "20241021_AC_Pkp2-MCMV_Mouse_Single_Heart_9"
  ),
  Sex = rep("male", 9),
  Genotype = c("Ctr", "het_KO", "het_KO", "Ctr", "het_KO", "het_KO", "Ctr", "het_KO", "het_KO"),
  Treatment = c("noninf", "noninf", "MCMV", "noninf", "noninf", "MCMV", "noninf", "noninf", "MCMV"),
  Condition = c("Ctr_noninf", "het_KO_noninf", "het_KO_MCMV", "Ctr_noninf", "het_KO_noninf", "het_KO_MCMV", "Ctr_noninf", "het_KO_noninf", "het_KO_MCMV"),
  Sample_ID = c(
    "B3_Ctr_noninf", "B4_het_KO_noninf", "B5_het_KO_MCMV",
    "B6_Ctr_noninf", "B7_het_KO_noninf", "B8_het_KO_MCMV",
    "B9_Ctr_noninf", "B10_het_KO_noninf", "B11_het_KO_MCMV"
  ),
  stringsAsFactors = FALSE
)
# Merge metadata
integrated$sample_short <- toupper(sub("_.*", "", integrated$sample))
meta_merged <- left_join(integrated@meta.data, meta_info, by = c("sample_short" = "sample"))

# Update the Seurat object's metadata
integrated@meta.data <- meta_merged
correct_cells <- rownames(Embeddings(integrated, "umap"))
rownames(integrated@meta.data) <- correct_cells

#### Perform DE analysis ####
Idents(integrated) <- "Condition"
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
write.csv(markers, "/Volumes/srivalli/dzhi/ACM_SN_R_2025/data/de_condition_all_markers.csv", row.names = FALSE)
write.csv(top10, "/Volumes/srivalli/dzhi/ACM_SN_R_2025/data/de_top_10_condition_clusters.csv", row.names = FALSE)

p_dot <- DotPlot(integrated, features = unique(top10$gene)) + RotatedAxis()
p_dot
ggsave("/Volumes/srivalli/dzhi/ACM_SN_R_2025/plots/dotplot_top5_condition_markers.png", p_dot, width = 20, height = 10)
