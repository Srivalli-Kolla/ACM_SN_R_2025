library(Seurat)
library(dplyr)
library(ggplot2)

#### Load processed object ####
integrated <- readRDS("/home/gruengroup/srivalli/Github/ACM_SN_R_2025/data/integrated_processed.rds")

#### DE analysis ####
markers <- FindAllMarkers(
  integrated,
  only.pos = TRUE,
  test.use = "wilcox",
  min.pct = 0.25,
  logfc.threshold = 0.25
)

top5 <- markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)

write.csv(markers, "/home/gruengroup/srivalli/Github/ACM_SN_R_2025/data/de_all_markers.csv", row.names = FALSE)
write.csv(top10, "/home/gruengroup/srivalli/Github/ACM_SN_R_2025/data/de_top_10_clusters.csv", row.names = FALSE)

# Dot plot for top 5
p_dot <- DotPlot(integrated, features = unique(top5$gene)) + RotatedAxis()
ggsave("/home/gruengroup/srivalli/Github/ACM_SN_R_2025/data/dotplot_top5_markers.png", p_dot, width = 10, height = 6)


saveRDS(markers, file = "/home/gruengroup/srivalli/Github/ACM_SN_R_2025/data/de_all_markers.rds")