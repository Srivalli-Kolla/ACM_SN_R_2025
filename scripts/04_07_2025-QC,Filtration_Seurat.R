#### Load libraries ####
library(Seurat)
library(dplyr)
library(ggplot2)
library(tibble)

#### Define samples and paths ####
sample_names <- c("B3_Lib2", "B4_Lib2", "B5_Lib2", "B6_Lib2", "B7_Lib2",
                  "B8_Lib2", "B9_Lib2", "B10_Lib2", "B11_Lib2")
data_paths <- paste0("/home/gruengroup/srivalli/Github/ACM_SN_R_2025/data/Library2/", sample_names, "/filtered_feature_bc_matrix/")

#### Load Seurat objects (unfiltered) ####
seurat_raw_list <- list()
for (i in seq_along(sample_names)) {
  counts <- Read10X(data.dir = data_paths[i])
  obj <- CreateSeuratObject(counts = counts, project = sample_names[i])
  obj$sample <- sample_names[i]
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
  seurat_raw_list[[sample_names[i]]] <- obj
}

#### Create merged object for before-filtering violin plot ####
combined_raw <- merge(seurat_raw_list[[1]], y = seurat_raw_list[-1], add.cell.ids = sample_names)
combined_raw$Sample_ID <- combined_raw$sample

#### Create output directory ####
dir.create("/home/gruengroup/srivalli/Github/ACM_SN_R_2025/plots", showWarnings = FALSE)

#### Violin Plot Before Filtering ####
p_violin_before <- VlnPlot(combined_raw,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "Sample_ID", pt.size = 0.1, ncol = 3
)
ggsave("/home/gruengroup/srivalli/Github/ACM_SN_R_2025/plots/qc_violin_before_filtered.png", p_violin_before, width = 12, height = 6)

#### QC filter the raw objects ####
seurat_list <- lapply(seurat_raw_list, function(obj) {
  subset(obj, subset = nFeature_RNA > 200 & percent.mt < 10)
})

#### Add metadata ####
combined <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = sample_names)
combined$sample_base <- gsub("_Lib2$", "", combined$sample)
meta_data <- read.csv("/home/gruengroup/srivalli/Github/ACM_SN_R_2025/data/acm_sn_metadata.csv", stringsAsFactors = FALSE)

combined@meta.data <- combined@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  left_join(meta_data, by = c("sample_base" = "sample")) %>%
  tibble::column_to_rownames("cell")

#### Violin Plot After Filtering ####
combined$Sample_ID <- combined$sample
p_violin_after <- VlnPlot(combined,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "Sample_ID", pt.size = 0.1, ncol = 3
)
ggsave("/home/gruengroup/srivalli/Github/ACM_SN_R_2025/plots/qc_violin_after_filtered.png", p_violin_after, width = 12, height = 6)

#### Barplot ####
cell_counts <- as.data.frame(table(combined$Sample_ID))
colnames(cell_counts) <- c("Sample_ID", "CellCount")
p_bar <- ggplot(cell_counts, aes(x = Sample_ID, y = CellCount)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = CellCount), vjust = -0.3, size = 3) +
  theme_minimal() +
  labs(title = "Cell Counts After Filtering") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/home/gruengroup/srivalli/Github/ACM_SN_R_2025/plots/cell_counts_barplot_filtered.png", p_bar, width = 10, height = 6)

#### SCTransform ####
seurat_list <- lapply(seurat_list, function(obj) {
  SCTransform(obj, verbose = FALSE)
})

#### Save processed list ####
saveRDS(seurat_list, "/home/gruengroup/srivalli/Github/ACM_SN_R_2025/data/seurat_list_sctransformed.rds")