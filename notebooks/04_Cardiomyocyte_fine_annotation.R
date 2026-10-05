# 04_Cardiomyocyte_fine_annotation.R
# Code-only export of the analysis pipeline. Data are not included; see README.md.

params <- list(rds_in = "results/02_subclustering/Cardiomyocyte_subcluster.rds", out_dir = "results/04_cm_annotation")


knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = TRUE, fig.align = "center", dpi = 150)

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyr)
})
source("_common.R")
knitr::opts_chunk$set(cache = TRUE, autodep = TRUE, cache.lazy = FALSE, cache.extra = input_md5(params))

path_input_rds    <- resolve_path(params$rds_in)
path_markers_csv  <- sub("\\.rds$", "_markers.csv", path_input_rds)
nb2_cleaning_tab <- read.csv(file.path(dirname(path_input_rds), "atlas_cleaning_clusters.csv"))
path_out      <- resolve_path(params$out_dir)
path_tables      <- file.path(path_out, "tables");  dir.create(path_tables, recursive = TRUE, showWarnings = FALSE)
path_figures      <- file.path(path_out, "figures"); dir.create(path_figures, recursive = TRUE, showWarnings = FALSE)
FILE_PREFIX  <- "cm"
path_table      <- function(name) file.path(path_tables, paste0(FILE_PREFIX, "_", name))
path_figure      <- function(name) file.path(path_figures, paste0(FILE_PREFIX, "_", name))
CLUSTER_COL  <- "seurat_clusters"
MERGED_COL   <- "merged_cluster"
DOUBLET_COL  <- "cellbender_doublet_scores"

cm_obj    <- readRDS(path_input_rds)
stopifnot(all(c("sample", "age", "condition", "batch", DESIGN_COVARIATES) %in% colnames(cm_obj@meta.data)))
subclustering_run <- cm_obj@misc$subclustering

LABEL_MAP     <- c("0" = "Ventricular CM",
                   "1" = "Low-content ambient CM",
                   "2" = "IFN-GTPase CM",
                   "3" = "Conduction-like CM",
                   "4" = "Nppa-high stress CM",
                   "5" = "Suspected artefact (Gm10800 family)")
LABEL_KIND    <- c("Ventricular CM"           = "subtype",
                   "Low-content ambient CM"   = "technical",
                   "IFN-GTPase CM"            = "subtype",
                   "Conduction-like CM"       = "subtype",
                   "Nppa-high stress CM"      = "subtype",
                   "Suspected artefact (Gm10800 family)" = "suspected artefact")
LABEL_MERGE   <- character()
MERGE_EXCLUDE <- character()
stopifnot(setequal(names(LABEL_KIND), LABEL_MAP))

GM10800_FAMILY <- c("Gm10800", "Gm10801", "Gm10718", "Gm21738")

GM10800_LOCI <- data.frame(gene = GM10800_FAMILY, chrom = c("chr2", "chr2", "chr9", "chr14"),
                           start = c(98666547, 98662237, 3023547, 19415857), end = c(98667301, 98664083, 3025218, 19418930),
                           strand = c("-", "+", "+", "-"), gene_type = "protein_coding", tsl = 5,
                           masked_bp = c(755, 1487, 1672, 3074), gsat_bp = c(755, 872, 1672, 3074),
                           gap_kb = c(NA, NA, 23.5, 0.8), gap_type = c(NA, NA, "centromere", "contig"))
PANEL_EXCLUDED <- c(Mural_contam = "Abcc9")
NB2_AMBIENT <- with(nb2_cleaning_tab, sort(table(best_other_lineage[removed & cell_type != subclustering_run$cell_type]), decreasing = TRUE))
IFN_ANIMAL <- names(ANIMAL_NOTES)[grepl("interferon", ANIMAL_NOTES)]
LOW_YIELD  <- names(QC_FLAGGED_ANIMALS)[grepl("^low_yield", QC_FLAGGED_ANIMALS)]

merged_clusters_of   <- function(label) names(LABEL_MAP)[LABEL_MAP %in% label]
subclusters_of   <- function(label) names(merge_of)[merge_of %in% merged_clusters_of(label)]
qc_value   <- function(label, col) qc_tab[[col]][qc_tab$cluster %in% subclusters_of(label)]
module_score   <- function(label, module) unlist(score_tab[score_tab$module == module, paste0("c", subclusters_of(label))])
module_score_others  <- function(label, module) unlist(score_tab[score_tab$module == module, paste0("c", setdiff(names(merge_of), subclusters_of(label)))])
top_cluster_of  <- function(module) sub("^c", "", names(which.max(unlist(score_tab[score_tab$module == module, -1]))))
detection_all <- function(gene) setNames(unlist(gene_tab[gene_tab$gene == gene, paste0("pct_c", names(merge_of))]), unname(LABEL_MAP[merge_of]))
detection_of  <- function(label, gene) unname(detection_all(gene)[names(detection_all(gene)) == label])
detection_others <- function(label, gene) detection_all(gene)[names(detection_all(gene)) != label]
detection_by_subcluster <- function(gene) setNames(unlist(gene_tab[gene_tab$gene == gene, paste0("pct_c", levels(cluster_ids))]), levels(cluster_ids))
pct_range     <- function(x, number_format = "%.1f") paste0(paste(unique(sprintf(number_format, range(x))), collapse = "-"), "%")
top_animals <- function(label, k = 2) head(sort(setNames(rowSums(support[, subclusters_of(label), drop = FALSE]), rownames(support)), decreasing = TRUE), k)
format_count   <- function(x) formatC(x, big.mark = ",", format = "d")
qc_others  <- function(label, col) qc_tab[[col]][!qc_tab$cluster %in% subclusters_of(label)]
own_share  <- function(label, animal) sum(support[animal, subclusters_of(label)]) / sum(support[animal, ])
best_auc   <- function(label) max(nb2_all$myAUC[nb2_all$cluster %in% subclusters_of(label)])
top_markers <- function(label, k = 3) paste(head(nb2_markers$gene[nb2_markers$cluster %in% subclusters_of(label)], k), collapse = ", ")

SUBTYPE_ORDER <- unname(LABEL_MAP[LABEL_KIND[LABEL_MAP] == "subtype"])
SUBTYPE_SHORT <- c("Ventricular CM" = "Ventricular", "IFN-GTPase CM" = "IFN-GTPase",
                   "Conduction-like CM" = "Conduction-like", "Nppa-high stress CM" = "Nppa-high stress")
LABEL_ORDER   <- c(SUBTYPE_ORDER, setdiff(unname(LABEL_MAP), SUBTYPE_ORDER))
QC_STRATA     <- setdiff(LABEL_ORDER, SUBTYPE_ORDER)
MERGED_GREY   <- "grey50"
SUBTYPE_COLS  <- subtype_palette(SUBTYPE_ORDER)
LABEL_COLS    <- c(SUBTYPE_COLS, setNames(ifelse(QC_STRATA %in% names(LABEL_MERGE), MERGED_GREY, TECH_GREY), QC_STRATA))
TEXT_COLS     <- setNames(colorspace::darken(LABEL_COLS, 0.3), names(LABEL_COLS))

print_versions(c("ggplot2", "dplyr", "Seurat", "tidyr"))

DefaultAssay(cm_obj) <- "RNA"
Idents(cm_obj) <- CLUSTER_COL
cm_obj$sample      <- as.character(cm_obj$sample)
cm_obj$qc_flag     <- unname(coalesce(sub(":.*", "", QC_FLAGGED_ANIMALS[cm_obj$sample]), "pass"))
cm_obj$animal_note <- unname(coalesce(ANIMAL_NOTES[cm_obj$sample], "none"))

module_panels <- list(
  PanCM          = c("Tnnt2","Ttn","Actn2","Myh6","Ryr2","Tnni3","Pln","Actc1"),
  Ventricular    = c("Myl2","Myl3","Hey2","Irx4","Fhl2"),
  Atrial         = c("Myl7","Myl4","Sln","Kcnj3","Nr2f2","Bmp10"),
  Conduction_VCS = c("Irx3","Gja5","Cntn2","Cpne5","Etv1"),
  IEG            = c("Fos","Fosb","Jun","Junb","Egr1","Atf3","Ier2","Nr4a1"),
  Remodelling    = c("Xirp1","Xirp2","Ankrd1","Flnc","Hspb1","Ccn2","Rcan1","Serpine1"),
  Fetal_stress   = c("Nppa","Nppb","Myh7","Acta1"),
  IFN_GTPase     = c("Iigp1","Igtp","Irgm1","Gbp2","Gbp6","Gbp7","Gbp10"),
  IFN_ISG        = c("Isg15","Rsad2","Ifit1","Ifit3","Oasl2","Stat1","Irf7"),
  FAO            = c("Cpt1b","Cpt2","Acadm","Acadvl","Hadha","Hadhb","Slc25a20","Decr1"),
  OXPHOS         = c("Ndufa4","Ndufb7","Cox7a1","Cox6c","Uqcrh","Atp5f1a","Chchd10"),
  Metabolic_regulators = c("Ppargc1a","Esrra","Esrrg","Ppara","Pdk4","Hk2"),
  Cell_cycle     = c("Mki67","Top2a","Ccnb1","Birc5","Cdk1"),

  Mito_tech      = c("mt-Co1","mt-Co2","mt-Co3","mt-Nd1","mt-Nd4","mt-Atp6","mt-Cytb"),
  Ribo_tech      = c("Rpl13a","Rplp0","Rpl28","Rplp1","Rpl17","Rps10","Rps3","Rps19","Rpl11","Rps11"),

  EC_contam      = c("Pecam1","Cdh5","Flt1","Ptprb","Egfl7","Fabp4","Cyyr1"),
  Fibro_contam   = c("Dcn","Pdgfra","Col1a1","Col1a2","Col3a1"),
  Mural_contam   = c("Rgs5","Pdgfrb","Kcnj8","Myh11","Notch3"),
  Myeloid_contam = c("Ptprc","F13a1","Mrc1","Cd163","Lyz2","Tyrobp","Lst1","Fcer1g")
)
stopifnot(!anyDuplicated(unlist(module_panels)))
module_panels <- lapply(module_panels, intersect, y = rownames(cm_obj))
module_panels <- module_panels[lengths(module_panels) > 0]
print(lengths(module_panels))

cm_obj <- AddModuleScore(cm_obj, features = module_panels, name = "MS_", ctrl = 50, seed = 1)
module_score_cols <- paste0("MS.", names(module_panels))
colnames(cm_obj@meta.data)[match(paste0("MS_", seq_along(module_panels)), colnames(cm_obj@meta.data))] <- module_score_cols

score_tab <- cm_obj@meta.data |>
  group_by(cluster = Idents(cm_obj)) |>
  summarise(across(all_of(module_score_cols), mean), .groups = "drop") |>
  pivot_longer(-cluster, names_to = "module", values_to = "score") |>
  mutate(module = sub("^MS\\.", "", module)) |>
  pivot_wider(names_from = cluster, values_from = score, names_prefix = "c")
write.csv(score_tab, path_table("module_scores_by_cluster.csv"), row.names = FALSE)
knitr::kable(score_tab, digits = 3,
             caption = sprintf("Mean module score per cluster (%s). Scores are read within panel (across clusters), not across panels.", CLUSTER_COL))

key <- intersect(unique(c(unlist(module_panels), GM10800_FAMILY, PANEL_EXCLUDED)), rownames(cm_obj))
expr_data   <- GetAssayData(cm_obj, layer = "data")[key, , drop = FALSE]
cluster_ids  <- Idents(cm_obj)
gene_tab <- data.frame(gene = key)
for (k in levels(cluster_ids)) {
  sub <- expr_data[, cluster_ids == k, drop = FALSE]
  gene_tab[[paste0("mean_c", k)]] <- round(Matrix::rowMeans(sub), 3)
  gene_tab[[paste0("pct_c",  k)]] <- round(100 * Matrix::rowMeans(sub > 0), 1)
}
write.csv(rename_with(gene_tab, ~ paste0(sub("^pct_", "", .x), "_pct"), starts_with("pct_")), path_table("canonical_markers_by_cluster.csv"), row.names = FALSE)
DT::datatable(rename_with(gene_tab, ~ paste(sub("^pct_", "", .x), "(%)"), starts_with("pct_")), rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE),
              caption = "Mean log-normalised expression and % positive nuclei per cluster for all panel genes, the Gm10800 family and the gene left out of the mural panel (Abcc9).")

nb2_all <- read.csv(path_markers_csv)
below_auc <- nb2_all |> group_by(cluster) |> summarise(best_AUC = max(myAUC)) |> filter(best_AUC < subclustering_run$auc_stop)
gm10800_top  <- nb2_all |> filter(gene %in% GM10800_FAMILY) |> slice_max(myAUC, n = 1, with_ties = FALSE)
nb2_markers <- nb2_all |>
  group_by(cluster) |>
  slice_max(myAUC, n = 5, with_ties = FALSE) |>
  ungroup() |>
  transmute(cluster, gene,
            AUC = round(myAUC, 3), avg_log2FC = round(avg_log2FC, 2), `pct.1 (proportion)` = pct.1, `pct.2 (proportion)` = pct.2)
knitr::kable(nb2_markers,
             caption = sprintf("Top-5 ROC markers per sub-cluster from the Notebook 2 run (all sub-clusters shown; without a marker at the Notebook 2 threshold AUC >= %s: %s).",
                               subclustering_run$auc_stop, if (nrow(below_auc)) paste("sub-cluster", below_auc$cluster, collapse = ", ") else "none"))

qc_tab <- cm_obj@meta.data |>
  group_by(cluster = Idents(cm_obj)) |>
  summarise(n = n(),
            median_nFeature  = median(nFeature_RNA),
            median_pct_mito  = median(percent_mito),
            median_exon_prop = median(exon_prop),
            median_doublet   = median(.data[[DOUBLET_COL]]),
            n_animals        = n_distinct(sample),
            max_animal       = names(which.max(table(sample))),
            max_animal_share = max(table(sample)) / n(),
            .groups = "drop") |>
  mutate(single_animal_dominated = max_animal_share > 0.5)

panel_detection <- sapply(module_panels[grepl("_contam$", names(module_panels))], function(panel_genes) Matrix::colSums(expr_data[panel_genes, , drop = FALSE] > 0) > 0)
qc_tab <- left_join(qc_tab, as.data.frame(panel_detection) |> mutate(cluster = cluster_ids) |> group_by(cluster) |>
                      summarise(across(everything(), mean, .names = "any_{.col}"), .groups = "drop"), by = "cluster")
knitr::kable(qc_tab |> rename(`median_pct_mito (%)` = median_pct_mito, `median_exon_prop (proportion)` = median_exon_prop,
                             `max_animal_share (proportion)` = max_animal_share) |> rename_with(~ paste(.x, "(proportion)"), starts_with("any_")), digits = 3,
             caption = "Per-cluster QC and animal support. exon_prop = exon-only UMIs / intron-inclusive UMIs (proxy for cytoplasmic carry-over); max_animal_share = largest single-animal share of the cluster; any_[panel]_contam = share of nuclei detecting at least one gene of that contamination panel.")

support <- as.data.frame.matrix(table(cm_obj$sample, Idents(cm_obj)))
write.csv(support, path_table("cluster_support_by_animal.csv"))
knitr::kable(support, caption = "Nuclei per animal (rows) and cluster (columns).")

nb2_removed <- nb2_cleaning_tab |> filter(cell_type == subclustering_run$cell_type, removed)
ward <- ward_merge(cm_obj, CLUSTER_COL, VariableFeatures(cm_obj), exclude = union(as.character(nb2_removed$cluster), MERGE_EXCLUDE))
merge_of <- ward$merge_of
cm_obj$merged_cluster <- factor(unname(merge_of[as.character(cm_obj[[CLUSTER_COL, drop = TRUE]])]), levels = unique(merge_of))
svglite::svglite(path_figure("ward_merge_dendrogram.svg"), width = 7, height = 4.5)
plot_ward_merge(ward, "Cardiomyocyte sub-clusters: Ward assessment")
invisible(dev.off())
plot_ward_merge(ward, "Cardiomyocyte sub-clusters: Ward assessment")
merge_tab <- data.frame(subcluster = names(merge_of), merged_cluster = unname(merge_of),
                        removed_before_merge = names(merge_of) %in% c(as.character(nb2_removed$cluster), MERGE_EXCLUDE),
                        n_nuclei = as.integer(table(cm_obj[[CLUSTER_COL, drop = TRUE]])[names(merge_of)]))
write.csv(merge_tab, path_table("ward_merge.csv"), row.names = FALSE)
knitr::kable(merge_tab, caption = sprintf("Ward assessment: Notebook 2 sub-cluster -> merged cluster (cut height %.2f = %s x %.2f).",
                                          ward$cut, WARD_MERGE_FRACTION, ward$max_dist))

label_map <- data.frame(merged_cluster = names(LABEL_MAP), label = unname(LABEL_MAP), kind = unname(LABEL_KIND[LABEL_MAP]),
                        merged_into = unname(coalesce(LABEL_MERGE[LABEL_MAP], "")),
                        n_nuclei = as.integer(table(cm_obj$merged_cluster)[names(LABEL_MAP)])) |>
  mutate(fate = case_when(kind == "subtype" ~ "retained", merged_into != "" ~ "merged", TRUE ~ "removed"))
knitr::kable(label_map, caption = "Label map: merged cluster, label, kind (subtype or QC stratum), parent subtype of a merged technical stratum, nuclei and fate before composition.")

merged_ids <- as.character(cm_obj@meta.data[[MERGED_COL]])
if (!all(merged_ids %in% names(LABEL_MAP))) stop("unmapped clusters: ", paste(setdiff(merged_ids, names(LABEL_MAP)), collapse = ", "))
cm_obj$CM_label <- factor(unname(LABEL_MAP[merged_ids]), levels = LABEL_ORDER)
cm_obj$CM_kind  <- unname(LABEL_KIND[as.character(cm_obj$CM_label)])
Idents(cm_obj) <- "CM_label"

qc_tab <- qc_tab |> mutate(merged_cluster = unname(merge_of[as.character(cluster)]), label = unname(LABEL_MAP[merged_cluster]),
                           kind = unname(LABEL_KIND[label]), .after = cluster)
write.csv(qc_tab |> rename(median_mito_pct = median_pct_mito, max_animal_share_prop = max_animal_share) |> rename_with(~ paste0(.x, "_prop"), starts_with("any_")),
          path_table("cluster_qc.csv"), row.names = FALSE)
animal_tab <- cm_obj@meta.data |>
  group_by(sample, age, condition, sex, seq, batch, qc_flag, animal_note) |>
  summarise(n_nuclei_input = n(),
            merged_share  = mean(CM_label %in% names(LABEL_MERGE)),
            removed_share = mean(!CM_label %in% c(names(LABEL_MERGE), names(LABEL_KIND)[LABEL_KIND == "subtype"])),
            median_exon_prop = median(exon_prop), .groups = "drop")

annotation_export <- data.frame(
  barcode   = colnames(cm_obj),
  cluster   = as.character(cm_obj@meta.data[[CLUSTER_COL]]),
  merged_cluster = merged_ids,
  CM_label  = as.character(cm_obj$CM_label),
  CM_kind   = cm_obj$CM_kind,
  sample    = cm_obj$sample,
  age       = cm_obj$age,
  condition = cm_obj$condition,
  qc_flag   = cm_obj$qc_flag,
  animal_note = cm_obj$animal_note,
  row.names = NULL
)
write.csv(annotation_export, file.path(path_tables, sub("\\.rds$", "_barcode_annotations.csv", basename(path_input_rds))),
          row.names = FALSE)

dot_panels <- list(
  "Pan-CM"       = c("Tnnt2","Ttn","Ryr2","Myh6"),
  "Ventr."       = c("Myl2","Myl3","Fhl2"),
  "Atrial"       = c("Myl7","Myl4","Sln"),
  "VCS"          = c("Irx3","Gja5","Cntn2"),
  "IEG"          = c("Fos","Jun","Atf3","Nr4a1"),
  "Remodelling"  = c("Xirp1","Ankrd1","Flnc","Hspb1","Ccn2","Rcan1"),
  "Fetal"        = c("Nppa","Nppb","Myh7","Acta1"),
  "IFN-\u03b3"   = c("Iigp1","Igtp","Gbp6"),
  "ISG"          = c("Isg15","Rsad2","Ifit1","Stat1"),
  "Mito"         = c("mt-Co1","mt-Co3","mt-Atp6"),
  "Ribo"         = c("Rplp0","Rpl13a","Rps3"),
  "Gm10800"      = GM10800_FAMILY,
  "Contam"       = c("Flt1","Ptprb","Dcn","Pdgfra","Rgs5","F13a1")
)
dot_panels <- lapply(dot_panels, intersect, y = rownames(cm_obj))
p_dot <- DotPlot(cm_obj, features = dot_panels, cols = c("lightgrey", "firebrick"), dot.scale = 6) +
  guides(colour = guide_colourbar(title = "Avg. expr.\n(scaled)"), size = guide_legend(title = "% expressing")) +
  labs(x = NULL, y = NULL, title = "Canonical marker panels by cluster label") +
  theme_house() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9),
        axis.text.y = element_text(size = 10, colour = TEXT_COLS[LABEL_ORDER]),
        plot.margin = margin(5, 5, 5, 15))
save_svg(p_dot, path_figure("annotation_dotplot.svg"), width = 17, height = 5)

umap_stratum_keys <- c(merged  = sprintf("cluster(s) %s: merged into %s", paste(merged_clusters_of(names(LABEL_MERGE)), collapse = ", "), paste(unique(LABEL_MERGE), collapse = ", ")),
             removed = sprintf("cluster(s) %s: removed", paste(merged_clusters_of(QC_STRATA[!QC_STRATA %in% names(LABEL_MERGE)]), collapse = ", ")))
umap_df <- as.data.frame(Embeddings(cm_obj, "umap")) |> setNames(c("UMAP_1", "UMAP_2")) |>
  mutate(label = cm_obj$CM_label,
         stratum = case_when(label %in% SUBTYPE_ORDER ~ as.character(label),
                             label %in% names(LABEL_MERGE) ~ umap_stratum_keys[["merged"]], TRUE ~ umap_stratum_keys[["removed"]]))
umap_label_pos <- umap_df |> filter(label %in% SUBTYPE_ORDER) |> group_by(label) |> summarise(across(c(UMAP_1, UMAP_2), median), .groups = "drop")
p_umap <- ggplot(umap_df, aes(UMAP_1, UMAP_2, colour = stratum)) +
  scattermore::geom_scattermore(pointsize = 2.2, alpha = .7, pixels = c(1400, 1400)) +
  ggrepel::geom_text_repel(data = umap_label_pos, aes(label = label), colour = "black", size = 3, seed = 1, min.segment.length = Inf) +
  scale_colour_manual(values = c(SUBTYPE_COLS, setNames(c(MERGED_GREY, TECH_GREY), umap_stratum_keys)), breaks = unname(umap_stratum_keys), name = "QC strata") +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  labs(title = "CM cluster labels (Harmony UMAP)", x = "UMAP 1", y = "UMAP 2") +
  theme_house() + theme(legend.position = "bottom", legend.direction = "vertical", axis.text = element_blank(), axis.ticks = element_blank())
save_svg(p_umap, path_figure("annotated_umap.svg"), width = 6.5, height = 6.5)

p_umap
p_dot

FILL_CAP <- 1.2
evidence <- cm_obj@meta.data |>
  group_by(CM_label) |>
  summarise(across(all_of(module_score_cols), mean), .groups = "drop") |>
  pivot_longer(-CM_label, names_to = "module", values_to = "score") |>
  mutate(module = sub("^MS\\.", "", module),
         kind = case_when(grepl("_contam$", module) ~ "contamination",
                          grepl("_tech$", module)   ~ "technical",
                          TRUE                      ~ "signature"))

p_heat <- ggplot(evidence, aes(CM_label, module, fill = score)) +
  geom_tile(colour = "white", linewidth = .5) +
  geom_text(aes(label = sprintf("%.2f", score), colour = abs(score) > FILL_CAP / 2), size = 3.5) +
  facet_grid(kind ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0,
                       limits = c(-FILL_CAP, FILL_CAP), oob = scales::squish,
                       name = sprintf("mean\nmodule score\n(fill capped\nat \u00b1%.1f)", FILL_CAP)) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey20"), guide = "none") +
  labs(x = NULL, y = NULL, title = "Mean module score per cluster label") +
  theme_house() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, colour = TEXT_COLS[LABEL_ORDER]),
        axis.text.y = element_text(size = 9), strip.text.y = element_text(angle = 0), panel.grid = element_blank())
save_svg(p_heat, path_figure("module_score_heatmap.svg"), width = 13, height = 6)
p_heat

keep_labels     <- c(names(LABEL_KIND)[LABEL_KIND == "subtype"], names(LABEL_MERGE))
remove_clusters <- names(LABEL_MAP)[!LABEL_MAP %in% keep_labels]
n_nuclei_before <- ncol(cm_obj)
cm_obj <- subset(cm_obj, cells = colnames(cm_obj)[!cm_obj[[MERGED_COL, drop = TRUE]] %in% remove_clusters])
cm_obj$CM_label   <- droplevels(cm_obj$CM_label)
cm_obj$CM_subtype <- factor(coalesce(unname(LABEL_MERGE[as.character(cm_obj$CM_label)]), as.character(cm_obj$CM_label)), levels = SUBTYPE_ORDER)
Idents(cm_obj) <- "CM_subtype"

clean_tab <- as.data.frame(table(CM_label = cm_obj$CM_label, CM_subtype = cm_obj$CM_subtype)) |> filter(Freq > 0)
knitr::kable(clean_tab, caption = sprintf("Nuclei per cluster label and inferential subtype after removing merged cluster(s) %s (%s of %s nuclei removed).",
                                          paste(remove_clusters, collapse = ", "), format_count(n_nuclei_before - ncol(cm_obj)), format_count(n_nuclei_before)))

animal_meta <- animal_tab |>
  left_join(count(cm_obj@meta.data, sample, name = "n_nuclei"), by = "sample") |>
  mutate(group = make_group(age, condition))
write.csv(rename(animal_meta, merged_share_prop = merged_share, removed_share_prop = removed_share), path_table("sample_metadata.csv"), row.names = FALSE)
knitr::kable(rename(animal_meta, `merged_share (proportion)` = merged_share, `removed_share (proportion)` = removed_share, `median_exon_prop (proportion)` = median_exon_prop), digits = 3,
             caption = "One row per animal: input and retained CM nuclei, share merged (technical stratum with parent) and removed (other QC strata), library-quality flag and biological note.")
sex_by_group    <- table(animal_meta$group, animal_meta$sex)
df_resid <- nrow(animal_meta) - qr(model.matrix(DESIGN_FORMULA, animal_meta))$rank

subtype_props <- speckle::getTransformedProps(clusters = cm_obj$CM_subtype, sample = cm_obj$sample, transform = "logit")
animal_props <- as.data.frame(subtype_props$Proportions) |>
  setNames(c("CM_subtype", "sample", "proportion")) |>
  left_join(animal_meta, by = "sample") |>
  mutate(CM_subtype = factor(CM_subtype, levels = SUBTYPE_ORDER),
         group = factor(as.character(group), levels = GROUP_ORDER))
write.csv(rename(animal_props, merged_share_prop = merged_share, removed_share_prop = removed_share), path_table("subtype_proportions_per_sample.csv"), row.names = FALSE)
low_min <- animal_props |> filter(sample == LOW_YIELD) |> slice_min(proportion, n = 1)

dominance <- cm_obj@meta.data |>
  count(CM_subtype, sample) |>
  group_by(CM_subtype) |>
  summarise(n_nuclei = sum(n), n_animals_present = n(), max_animal = sample[which.max(n)],
            max_animal_share = max(n) / sum(n), .groups = "drop") |>
  mutate(single_animal_dominated = max_animal_share > 0.5,
         max_animal_note = unname(coalesce(ANIMAL_NOTES[max_animal], "none")))
knitr::kable(rename(dominance, `max_animal_share (proportion)` = max_animal_share), digits = 3, caption = "Per-subtype animal support after cleaning.")

composition_res <- fit_composition(subtype_props, animal_meta, covariates = DESIGN_COVARIATES) |>
  mutate(subtype = factor(subtype, levels = SUBTYPE_ORDER)) |>
  left_join(dominance |> select(subtype = CM_subtype, n_animals_present, max_animal, max_animal_share), by = "subtype") |>
  arrange(contrast, subtype)
write.csv(rename(composition_res, max_animal_share_prop = max_animal_share), path_table("propeller_results.csv"), row.names = FALSE)
top_result <- composition_res[which.min(composition_res$FDR), ]
se_ratio <- composition_res |> group_by(subtype) |>
  summarise(se_ratio = SE[contrast == "Interaction"] / mean(SE[contrast != "Interaction"]), .groups = "drop") |> pull(se_ratio)

group_pct <- function(subtype_label, group_id, drop = character()) 100 * animal_props$proportion[animal_props$CM_subtype == subtype_label & animal_props$group == group_id & !animal_props$sample %in% drop]
animal_pct <- function(subtype_label, animal) 100 * animal_props$proportion[animal_props$CM_subtype == subtype_label & animal_props$sample == animal]
ifn_group  <- as.character(animal_meta$group[animal_meta$sample == IFN_ANIMAL])
nppa_animal <- as.character(dominance$max_animal[dominance$CM_subtype == "Nppa-high stress CM"])
result_of  <- function(subtype_label, contrast_name, column) composition_res[[column]][composition_res$subtype == subtype_label & composition_res$contrast == contrast_name]
design_matrix     <- model.matrix(DESIGN_FORMULA, animal_meta[match(colnames(subtype_props$TransformedProps), animal_meta$sample), ])
rss_ifn <- apply(residuals(limma::lmFit(subtype_props$TransformedProps, design_matrix), subtype_props$TransformedProps)^2, 1,
                 function(animal_rss) animal_rss[IFN_ANIMAL] / sum(animal_rss))

noted_subtypes <- dominance |> filter(max_animal_note != "none")
noted_text <- paste(mapply(function(subtype_label, animal_id, animal_share) {
  own_pct <- 100 * animal_props$proportion[animal_props$CM_subtype == subtype_label & animal_props$sample == animal_id]
  other_pct <- 100 * animal_props$proportion[animal_props$CM_subtype == subtype_label & animal_props$sample != animal_id]
  sprintf("%s: %.0f%% of its nuclei come from %s (%s); the subtype is %.1f%% of that animal's CM nuclei versus %.1f-%.1f%% in the other animals",
          subtype_label, 100 * animal_share, animal_id, ANIMAL_NOTES[animal_id], own_pct, min(other_pct), max(other_pct))
}, as.character(noted_subtypes$CM_subtype), noted_subtypes$max_animal, noted_subtypes$max_animal_share), collapse = "; ")

knitr::kable(
  composition_res |> select(-n_animals) |>
    mutate(across(c(logit_effect, SE, CI_low, CI_high, max_animal_share), ~ round(.x, 3)),
           across(c(p_value, FDR, FDR_all_contrasts), ~ signif(.x, 3))) |>
    rename(`max_animal_share (proportion)` = max_animal_share),
  caption = paste("Animal-level limma results on logit-transformed CM subtype proportions, adjusted for sex and sequencing run (~ 0 + group + sex + seq).",
                  "Positive values indicate higher relative abundance in the first-named group; a positive",
                  "interaction means a more positive HFpEF effect in old than in young animals. FDR is within",
                  "contrast; FDR_all_contrasts is across all tests. Intervals are pointwise 95% intervals.",
                  "logit_effect, SE and CI in logit units; max_animal / max_animal_share (proportion): largest single-animal contributor to the subtype.")
)

effect_df <- composition_res |>
  mutate(subtype = factor(as.character(subtype), levels = rev(SUBTYPE_ORDER)),
         contrast_label = factor(CONTRAST_LABEL[as.character(contrast)], levels = unname(CONTRAST_LABEL)),
         star = fdr_stars(FDR))

p_effects <- ggplot(effect_df, aes(logit_effect, subtype, colour = subtype)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = .4) +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high), width = .25, orientation = "y", linewidth = .45) +
  geom_point(size = 2.2) +
  geom_text(aes(label = star), colour = "black", nudge_y = .3, fontface = "bold", size = 4) +
  facet_wrap(~ contrast_label, nrow = 1) +
  scale_colour_manual(values = TEXT_COLS[SUBTYPE_ORDER], guide = "none") +
  scale_y_discrete(labels = SUBTYPE_SHORT) +
  labs(x = "logit-proportion effect (95% CI)", y = NULL,
       title = "Animal-level CM subtype composition contrasts",
       subtitle = "limma on logit proportions (~ 0 + group + sex + seq). Positive = higher in the first-named group; stars: within-contrast BH-FDR (* <.05, ** <.01, *** <.001)") +
  theme_house()

strip_lab <- setNames(paste0(SUBTYPE_SHORT[SUBTYPE_ORDER],
                             ifelse(dominance$single_animal_dominated[match(SUBTYPE_ORDER, dominance$CM_subtype)], " \u2020", "")),
                      SUBTYPE_ORDER)
per_animal_df <- animal_props |> mutate(x = as.numeric(group), x_j = x + withr::with_seed(1, runif(nrow(animal_props), -.12, .12)),
                   point_lab = ifelse(qc_flag != "pass", sprintf("%s (n = %s)", sample, format_count(n_nuclei)), sample))
p_per_animal <- ggplot(per_animal_df, aes(y = 100 * proportion)) +
  stat_summary(aes(x = x, fill = group, colour = CM_subtype), fun = mean, geom = "col", width = .7, linewidth = .8) +
  geom_point(aes(x = x_j, shape = qc_flag == "pass"), size = 1.5, alpha = .85) +
  ggrepel::geom_text_repel(data = filter(per_animal_df, qc_flag != "pass" | animal_note != "none"), aes(x = x_j, label = point_lab),
                           size = 2.4, min.segment.length = 0, seed = 1) +
  facet_wrap(~ CM_subtype, scales = "free_y", nrow = 1, labeller = as_labeller(strip_lab)) +
  scale_x_continuous(breaks = seq_along(GROUP_ORDER), labels = GROUP_ORDER) +
  scale_fill_manual(values = GROUP_FILL, labels = GROUP_LABEL, name = NULL) +
  scale_colour_manual(values = TEXT_COLS[SUBTYPE_ORDER], guide = "none") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "qc_flag = pass", `FALSE` = "library-quality flag"), name = NULL) +
  labs(x = NULL, y = "% of recovered CM nuclei",
       title = "Per-animal CM subtype proportions",
       subtitle = "Bars: group means; points: animals (labelled: qc_flag, with retained nuclei, or animal_note). \u2020 more than 50% of the subtype's nuclei come from one animal.") +
  theme_house() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom")

save_svg(p_effects, path_figure("propeller_effectsizes.svg"), width = 13, height = 4.2)
save_svg(p_per_animal, path_figure("propeller_composition.svg"), width = 12, height = 4.2)
p_effects
p_per_animal

cm_obj@misc$annotation <- list(cluster_col = CLUSTER_COL, merged_col = MERGED_COL, merge_of = merge_of,
                           ward_merge_fraction = WARD_MERGE_FRACTION, ward_cut = ward$cut, label_map = LABEL_MAP, lab_kind = LABEL_KIND, lab_merge = LABEL_MERGE,
                           removed_clusters = remove_clusters, model = deparse(DESIGN_FORMULA),
                           input_rds = basename(path_input_rds))
saveRDS(cm_obj, file.path(path_out, sub("\\.rds$", "_annotated_clean.rds", basename(path_input_rds))))

session_footer()
