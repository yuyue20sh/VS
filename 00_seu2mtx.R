# Convert Seurat object from AJP dataset to mtx format data


library(Seurat)
library(Matrix)


seu2mtx <- function(seu_obj, save_dir) {
  #
  # Convert seurat object to mtx format data
  #
  # Args:
  #   seu_obj: Seurat object
  #   save_dir: str, path to save mtx data
  #
  # Returns:
  #   None
  #
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

  matrix <- GetAssayData(seu_obj, assay = "RNA", layer = "counts")
  barcodes <- data.frame(colnames(seu_obj))
  features <- data.frame(rownames(seu_obj))
  meta <- seu_obj@meta.data

  write.table(barcodes, paste0(save_dir, "/barcodes.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(features, paste0(save_dir, "/features.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  writeMM(matrix, paste0(save_dir, "/matrix.mtx"))
  write.table(meta, paste0(save_dir, "/metadata.tsv"),
              sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
}


#### Main ####

data_path <- "./data/ajp.rds"
save_dir <- "./outs/mtx/ajp/"


# load data
ajp <- readRDS(data_path)

# create new ajp object
ajp_matrix <- GetAssayData(ajp, assay = "RNA", layer = "counts")

ajp_meta <- data.frame(row.names = colnames(ajp))
ajp_meta["gsm"] <- "-"
ajp_meta["sample"] <- ajp@meta.data[["orig.ident"]]
ajp_meta["gse"] <- "-"
ajp_meta["dataset"] <- "AJP"

ajp <- CreateSeuratObject(ajp_matrix, meta.data = ajp_meta, project = "AJP")
ajps <- SplitObject(ajp, split.by = "sample")

# convert to mtx
for (i in seq_along(ajps)) {
  out_dir <- paste0(save_dir, "/", names(ajps)[[i]])
  seu2mtx(ajps[[i]], out_dir)
}

print("Done!")
