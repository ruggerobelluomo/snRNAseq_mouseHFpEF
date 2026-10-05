
install.packages(c("BiocManager", "remotes"))
install.packages(c("rmarkdown", "knitr", "Seurat", "SeuratObject", "harmony", "msigdbr", "NMF", "data.table", "dplyr", "tidyr", "tibble", "readr", "stringr", "ggplot2", "ggrepel", "patchwork", "scales", "scattermore", "svglite", "DT", "Matrix", "digest", "future"))
BiocManager::install(c("MAST", "SingleCellExperiment", "limma", "edgeR", "speckle", "fgsea", "GenomeInfoDbData"))
remotes::install_github("jinworks/CellChat")
