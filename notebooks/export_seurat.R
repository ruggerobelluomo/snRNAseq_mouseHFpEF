
suppressPackageStartupMessages({library(Seurat); library(Matrix)})
source("_common.R")
path_export  <- file.path(PATHS$objects, "for_seurat")
path_out_rds <- file.path(PATHS$objects, "integrated_LV_snRNAseq_seurat.rds")
path_out_meta <- file.path(PATHS$objects, "atlas_metadata.csv.gz")

cat("reading counts...\n")
counts   <- as(readMM(gzfile(file.path(path_export, "counts.mtx.gz"))), "CsparseMatrix")
genes    <- read.csv(gzfile(file.path(path_export, "genes.csv.gz")))$gene
barcodes <- read.csv(gzfile(file.path(path_export, "barcodes.csv.gz")))$barcode
dimnames(counts) <- list(genes, barcodes)
cat(sprintf("  %d genes x %d nuclei\n", nrow(counts), ncol(counts)))

nucleus_meta <- read.csv(gzfile(file.path(path_export, "metadata.csv.gz")), row.names = 1, check.names = FALSE)
stopifnot(all(barcodes %in% rownames(nucleus_meta)), ncol(counts) == nrow(nucleus_meta))
nucleus_meta <- nucleus_meta[barcodes, , drop = FALSE]
nucleus_meta$cellbender_predicted_doublets <- NULL
for (col in c("sample", "individual", "sex", "batch", "seq", "qc_flag")) nucleus_meta[[col]] <- factor(nucleus_meta[[col]])
nucleus_meta$leiden    <- factor(nucleus_meta$leiden, levels = as.character(sort(unique(as.integer(nucleus_meta$leiden)))))
nucleus_meta$disease   <- factor(nucleus_meta$disease, levels = GROUP_ORDER)
nucleus_meta$age       <- factor(nucleus_meta$age, levels = c("young", "old"))
nucleus_meta$condition <- factor(nucleus_meta$condition, levels = c("control", "HFpEF"))
nucleus_meta$cell_type <- factor(nucleus_meta$cell_type, levels = c(CELL_TYPE_ORDER, setdiff(unique(nucleus_meta$cell_type), CELL_TYPE_ORDER)))
stopifnot(!anyNA(nucleus_meta$disease), !anyNA(nucleus_meta$cell_type))

atlas_obj <- CreateSeuratObject(counts = counts, meta.data = nucleus_meta, project = "LV_snRNAseq", assay = "RNA")

atlas_obj <- NormalizeData(atlas_obj, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)

for (reduction in c("pca", "harmony", "umap")) {
  embedding <- as.matrix(read.csv(gzfile(file.path(path_export, paste0(reduction, ".csv.gz"))), header = FALSE, row.names = 1))
  embedding <- embedding[barcodes, , drop = FALSE]
  reduction_key <- switch(reduction, pca = "PC_", harmony = "harmony_", umap = "UMAP_")
  colnames(embedding) <- paste0(reduction_key, seq_len(ncol(embedding)))
  atlas_obj[[reduction]] <- CreateDimReducObject(embeddings = embedding, key = reduction_key, assay = "RNA")
}
Idents(atlas_obj) <- "cell_type"

cat("\n=== object ===\n")
cat("dim:", dim(atlas_obj), "\n")
cat("layers:", paste(Layers(atlas_obj[["RNA"]]), collapse = ", "), "\n")
cat("reductions:", paste(names(atlas_obj@reductions), collapse = ", "), "\n")
cat("metadata columns:", ncol(atlas_obj@meta.data), "\n")
print(table(atlas_obj$cell_type))
print(table(atlas_obj$sample, atlas_obj$qc_flag))
data_check <- GetAssayData(atlas_obj, layer = "data")[, 1:200]; counts_check <- counts[, 1:200]
cat("max |data - log1p(CP10k)| on 200 nuclei:", max(abs(data_check - log1p(t(t(counts_check) / colSums(counts_check)) * 1e4))), "\n")

saveRDS(atlas_obj, path_out_rds)
nucleus_meta_out <- atlas_obj@meta.data[, c("sample", "individual", "disease", "age", "condition", "sex", "batch", "seq", "qc_flag",
                              "cellbender_ncount", "cellbender_ngenes", "percent_mito", "percent_mito_cellbender", "exon_prop",
                              "cellranger_doublet_scores", "leiden", "cell_type")]
write.csv(nucleus_meta_out, gzfile(path_out_meta), row.names = TRUE)
cat("wrote", path_out_rds, "\nwrote", path_out_meta, "\n")
session_footer()
