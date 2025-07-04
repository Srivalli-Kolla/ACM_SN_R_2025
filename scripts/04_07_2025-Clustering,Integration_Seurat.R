library(Seurat)
library(ggplot2)
library(future)

# Parallel processing setup
plan("multicore", workers = 4)  
options(future.globals.maxSize = 100 * 1024^3)  

# Load SCTransformed Seurat list
seurat_list <- readRDS("/home/gruengroup/srivalli/Github/ACM_SN_R_2025/data/seurat_list_sctransformed.rds")

# Feature selection
features <- SelectIntegrationFeatures(seurat_list, nfeatures = 3000)

# Prep for SCT integration
seurat_list <- PrepSCTIntegration(seurat_list, anchor.features = features)

# Find anchors
anchors <- FindIntegrationAnchors(seurat_list, normalization.method = "SCT", anchor.features = features)

# Integrate data
integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

# Continue downstream
DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated)
integrated <- RunPCA(integrated, npcs = 30)

# Elbow plot
elbow <- ElbowPlot(integrated)
ggsave("/home/gruengroup/srivalli/Github/ACM_SN_R_2025/plots/elbow_plot.png", elbow, width = 7, height = 5)

# Clustering and UMAP
integrated <- FindNeighbors(integrated, dims = 1:30)
integrated <- FindClusters(integrated, resolution = 0.3)
integrated <- RunUMAP(integrated, dims = 1:30)

# Save final integrated object
saveRDS(integrated, "/home/gruengroup/srivalli/Github/ACM_SN_R_2025/data/integrated_processed.rds")