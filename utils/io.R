library(Seurat)
library(Matrix)
library(data.table)


seu2mtx <- function(seu_obj, save_dir) {
  #
  # Convert seurat object to mtx format data.
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

  fwrite(barcodes, paste0(save_dir, "/barcodes.tsv"),
         sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  fwrite(features, paste0(save_dir, "/features.tsv"),
         sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  writeMM(matrix, paste0(save_dir, "/matrix.mtx"))
  fwrite(meta, paste0(save_dir, "/metadata.tsv"),
         sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
}


read_mtx <- function(data_dir, use_symbol = TRUE) {
  #
  # Read mtx format data.
  # input files:  matrix.mtx, barcodes.tsv, features.tsv, and metadata.tsv
  #
  # Args:
  #     data_dir: str, path to the directory containing the matrix.mtx,
  #         barcodes.tsv, features.tsv, and metadata.tsv files
  #     use_symbol: bool, whether to use gene symbol (the second column) as row
  #          names
  #
  # Returns:
  #     list, the matrix and meta data
  #
  matrix <- readMM(file.path(data_dir, "matrix.mtx"))
  barcodes <- read.table(file.path(data_dir, "barcodes.tsv"), header = FALSE)
  features <- read.table(file.path(data_dir, "features.tsv"), header = FALSE)
  colnames(matrix) <- barcodes[[1]]
  if (use_symbol) {
    rownames(matrix) <- features[[2]]
  } else {
    rownames(matrix) <- features[[1]]
  }
  meta <- read.table(file.path(data_dir, "metadata.tsv"), sep = "\t",
                     header = TRUE, row.names = 1)

  return(list(matrix, meta))
}
