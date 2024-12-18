# Integrate data using RPCA.


library(Seurat)
library(ggplot2)
library(Matrix)
library(data.table)

source("./utils/io.R")

options(future.globals.maxSize = 1000 * 1024^3)


#### Main ####

data_dir <- "./outs/mtx/concatenated/"
seu_file <- "./outs/seu/rpca_integrated.seu.rds"
mtx_dir <- "./outs/mtx/integrated/"
fig_dir <- "./outs/figs/integrate/"
project_name <- "vs"


# load data
data <- read_mtx(data_dir, use_symbol = FALSE)
obj <- CreateSeuratObject(counts = data[[1]], meta.data = data[[2]],
                          project = project_name)

obj[["RNA"]] <- split(obj[["RNA"]], f = obj$sample)

# preprocess
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)

obj <- FindNeighbors(obj, dims = 1:30, reduction = "pca")
obj <- FindClusters(obj, resolution = 2, cluster.name = "unintegrated_cluster")

obj <- RunUMAP(obj, dims = 1:30, reduction = "pca",
               reduction.name = "umap.unintegrated")

# rpca
obj <- IntegrateLayers(object = obj, method = RPCAIntegration,
                       orig.reduction = "pca",
                       new.reduction = "rpca",
                       verbose = FALSE)

obj <- FindNeighbors(obj, reduction = "rpca", dims = 1:30)
obj <- FindClusters(obj, resolution = 2, cluster.name = "rpca_cluster")

obj <- RunUMAP(obj, reduction = "rpca", dims = 1:30,
               reduction.name = "umap.rpca")

# plot
p_unintegrated <- DimPlot(obj, group.by = c("sample", "unintegrated_cluster"),
                          reduction = "umap.unintegrated", shuffle = TRUE)
p_rpca <- DimPlot(obj, group.by = c("sample", "rpca_cluster"),
                  reduction = "umap.rpca", shuffle = TRUE)

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(paste0(fig_dir, "/", "unintegrated_umap.png"), p_unintegrated,
       width = 30, height = 10, units = "cm")
ggsave(paste0(fig_dir, "/", "rpca_umap.png"), p_rpca,
       width = 30, height = 10, units = "cm")

obj <- JoinLayers(obj)

# save
dir.create(dirname(seu_file), recursive = TRUE, showWarnings = FALSE)
saveRDS(obj, seu_file)
seu2mtx(obj, mtx_dir)
# save pca and rpca
fwrite(obj[["pca"]]@cell.embeddings,
       paste0(mtx_dir, "/", "pca_cell_embeddings.tsv"),
       sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
fwrite(obj[["pca"]]@feature.loadings,
       paste0(mtx_dir, "/", "pca_cell_loadings.tsv"),
       sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
fwrite(obj[["rpca"]]@cell.embeddings,
       paste0(mtx_dir, "/", "rpca_cell_embeddings.tsv"),
       sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
fwrite(obj[["rpca"]]@feature.loadings,
       paste0(mtx_dir, "/", "rpca_cell_loadings.tsv"),
       sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

print("Done!")
