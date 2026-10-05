
suppressPackageStartupMessages({ library(dplyr); library(readr); library(msigdbr) })
path_resources <- file.path("data", "resources")

uniprot_query <- paste0("https://rest.uniprot.org/uniprotkb/stream?format=tsv&fields=accession,gene_primary,cc_subcellular_location",
                        "&query=(taxonomy_id:10090)%20AND%20(reviewed:true)%20AND%20(keyword:KW-0964)")
raw <- readr::read_tsv(uniprot_query, show_col_types = FALSE) |> setNames(c("uniprot", "gene", "location"))

raw$primary <- vapply(strsplit(raw$location, "SUBCELLULAR LOCATION: "), function(location_parts) {
  location_parts <- location_parts[nzchar(location_parts) & !grepl("^\\[", location_parts)]; if (length(location_parts)) location_parts[1] else "" }, character(1))
secretome_tab <- raw |> filter(!is.na(gene), grepl("^Secreted", primary)) |>
  arrange(gene, nchar(uniprot), uniprot) |> distinct(gene, .keep_all = TRUE) |>
  transmute(gene, uniprot, primary_location = sub(" \\{.*$|[.;].*$", "", primary),
            retrieved_on = as.character(Sys.Date()), source_query = uniprot_query)
readr::write_csv(secretome_tab, file.path(path_resources, "mouse_secretome_swissprot.csv"))

fetch_msigdb <- function(...) { msig_df <- msigdbr(db_species = "HS", species = "Mus musculus", ...); lapply(split(msig_df$gene_symbol, msig_df$gs_name), unique) }
gs_snapshot <- list(sets = list(Hallmark = fetch_msigdb(collection = "H"),
                                Reactome = fetch_msigdb(collection = "C2", subcollection = "CP:REACTOME"),
                                `GO:BP`  = fetch_msigdb(collection = "C5", subcollection = "GO:BP")),
                    db_version = unique(msigdbr(db_species = "HS", species = "Mus musculus", collection = "H")$db_version),
                    msigdbr = as.character(packageVersion("msigdbr")), retrieved_on = Sys.Date(),
                    db_species = "HS", species = "Mus musculus")
saveRDS(gs_snapshot, file.path(path_resources, "ec_gsea_gene_sets.rds"))

collectri_url <- "https://zenodo.org/records/8192729/files/CollecTRI_regulons.csv?download=1"
title_case <- function(symbol) paste0(toupper(substr(symbol, 1, 1)), tolower(substring(symbol, 2)))
collectri_human <- read.csv(collectri_url)
collectri_mouse <- unique(data.frame(source = title_case(collectri_human$source), target = title_case(collectri_human$target),
                                     weight = collectri_human$weight))
write.csv(collectri_mouse, file.path(path_resources, "collectri_mouse_network.csv"), row.names = FALSE)

cat("wrote", file.path(path_resources, c("mouse_secretome_swissprot.csv", "ec_gsea_gene_sets.rds", "collectri_mouse_network.csv")), sep = "\n")
print(sessionInfo())
