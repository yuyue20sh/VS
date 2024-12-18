# Convert Seurat object from AJP dataset to mtx format data.


library(Seurat)

source("./utils/io.R")


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
