# 09_HFpEF_model_validation.R
# Code-only export of the analysis pipeline. Data are not included; see README.md.

params <- list(atlas_rds = "data/objects/integrated_LV_snRNAseq_seurat_clean.rds", cm_barcodes_csv = "results/04_cm_annotation/tables/Cardiomyocyte_subcluster_barcode_annotations.csv", 
    ec_barcodes_csv = "results/03_ec_annotation/tables/Endothelial_subcluster_barcode_annotations.csv", 
    cm_props_csv = "results/04_cm_annotation/tables/cm_subtype_proportions_per_sample.csv", 
    cm_propeller_csv = "results/04_cm_annotation/tables/cm_propeller_results.csv", 
    mast_csv = "results/05_ec_de/ec_DE_all_MAST_annotated.csv", universe_csv = "results/05_ec_de/ec_DE_all_MAST_gsea_universe.csv", 
    out_dir = "results/09_hfpef_validation")


knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = TRUE,
                      fig.align = "center", dpi = 150)
options(width = 110)

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(limma); library(speckle); library(tidyr); library(patchwork)
})
source("_common.R")
knitr::opts_chunk$set(cache = TRUE, autodep = TRUE, cache.lazy = FALSE, cache.extra = input_md5(params))
path_atlas        <- resolve_path(params$atlas_rds)
path_cm_barcodes  <- resolve_path(params$cm_barcodes_csv)
path_ec_barcodes  <- resolve_path(params$ec_barcodes_csv)
path_cm_props     <- resolve_path(params$cm_props_csv)
path_cm_propeller <- resolve_path(params$cm_propeller_csv)
path_mast_csv     <- resolve_path(params$mast_csv)
path_universe_csv <- resolve_path(params$universe_csv)
path_out     <- results_dir("09_hfpef_validation")
path_figures <- file.path(path_out, "figures"); dir.create(path_figures, recursive = TRUE, showWarnings = FALSE)
format_count <- function(x) formatC(x, big.mark = ",", format = "d")

SEED          <- 1
MIN_PB_NUCLEI <- 25
MODULE_CTRL   <- 50
SCALE_FACTOR  <- 1e4
FDR_CUT       <- 0.05
PADJ_CUT      <- 0.05
LFC_CUT       <- 0.2
LINEAGES <- c("Cardiomyocyte", "Endothelial", "Fibroblast", "Macrophage")
NPPA_CM  <- "Nppa-high stress CM"

GENE_SETS <- list(
  hypertrophy   = c("Nppa", "Nppb", "Myh7", "Acta1", "Ankrd1", "Xirp2"),
  calcium       = c("Atp2a2", "Pln", "Ryr2", "Slc8a1", "Camk2d"),
  metabolic     = c("Pdk4", "Cd36", "Ucp3", "Hmgcs2", "Angptl4", "Acot1", "Cpt1a"),
  fb_activation = c("Postn", "Col1a1", "Col3a1", "Cthrc1", "Comp", "Fn1", "Lox", "Timp1", "Ltbp2"),
  senescence    = c("Cdkn2a", "Cdkn1a", "Trp53"),
  inflammatory  = c("Il1b", "Tnf", "Ccl2", "Cxcl2", "Nlrp3"),
  upr           = c("Hspa5", "Atf4", "Ddit3", "Xbp1"))
EC_ACT_GENES <- c("Icam1", "Vcam1", "Sele", "Pecam1", "Nos3")
MODULE_LABEL <- c(hypertrophy = "hypertrophy module", calcium = "calcium-handling module", metabolic = "metabolic module",
                  fb_activation = "FB activation module", senescence = "senescence module", inflammatory = "inflammatory module",
                  upr = "UPR module")

LINEAGE_PLAN <- list(
  Cardiomyocyte = list(modules = c("hypertrophy", "calcium", "metabolic", "senescence", "upr"),
                       genes = unlist(GENE_SETS[c("hypertrophy", "calcium", "metabolic", "senescence", "upr")], use.names = FALSE)),
  Endothelial   = list(modules = "senescence", genes = c(GENE_SETS$senescence, EC_ACT_GENES)),
  Fibroblast    = list(modules = c("fb_activation", "senescence"), genes = unlist(GENE_SETS[c("fb_activation", "senescence")], use.names = FALSE)),
  Macrophage    = list(modules = c("inflammatory", "senescence"), genes = c(GENE_SETS$inflammatory, "Nos2", GENE_SETS$senescence)))
SECTION_LABEL <- c(S1 = "1 Hypertrophy / fetal genes", S2 = "2 Calcium handling", S3 = "3 Metabolic (HFD)", S4 = "4 Fibrosis",
                   S5 = "5 Senescence", S7 = "7 Inflammation and UPR")

READOUTS <- bind_rows(
  tibble(section = "S1", lineage = "Cardiomyocyte", kind = c("module", rep("gene", 6), "composition"),
         feature = c("hypertrophy", GENE_SETS$hypertrophy, "Nppa-high CM share")),
  tibble(section = "S2", lineage = "Cardiomyocyte", kind = c("module", rep("gene", 5)), feature = c("calcium", GENE_SETS$calcium)),
  tibble(section = "S3", lineage = "Cardiomyocyte", kind = "module", feature = "metabolic"),
  tibble(section = "S4", lineage = c("Fibroblast", "atlas"), kind = c("module", "composition"), feature = c("fb_activation", "fibroblast proportion")),
  tibble(section = "S5", lineage = LINEAGES, kind = "module", feature = "senescence"),
  expand_grid(section = "S5", lineage = LINEAGES, kind = "gene", feature = GENE_SETS$senescence),
  tibble(section = "S7", lineage = c("Macrophage", "Macrophage", "Cardiomyocyte"), kind = c("module", "gene", "module"),
         feature = c("inflammatory", "Nos2", "upr"))) |>
  mutate(lineage_short = coalesce(unname(CELL_TYPE_SHORT[lineage]), "atlas"),
         label = paste0(ifelse(kind == "module", MODULE_LABEL[feature], feature), " (", lineage_short, ")"),
         unit = case_when(kind == "module" ~ "module score", kind == "gene" ~ "mean log-normalised expression",
                          TRUE ~ "proportion (0-1); model on logit scale"))
ANIMAL_CLASS_SHAPE <- c(pass = 21, flag = 1, note = 24)
ANIMAL_CLASS_LABEL <- c(pass = "qc pass",
                        flag = sprintf("library-quality flag (%s)", paste(names(QC_FLAGGED_ANIMALS), collapse = ", ")),
                        note = sprintf("interferon note (%s)", paste(names(ANIMAL_NOTES), collapse = ", ")))
theme_fig <- theme_house() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   axis.line = element_line(colour = "grey60", linewidth = .3))
save_fig <- function(fig, name, width, height) {
  save_svg(fig, paste0(name, ".svg"), width = width, height = height, dir = path_figures)
  ggsave(file.path(path_figures, paste0(name, ".png")), fig, width = width, height = height, dpi = 200, bg = "white")
}
fmt_fdr <- function(p) formatC(p, format = "g", digits = 2)
md_escape <- function(x) gsub("*", "\\*", x, fixed = TRUE)
contrast_key <- paste(CONTRAST_CODE, CONTRAST_LABEL, sep = " = ", collapse = "; ")

print_versions(c("ggplot2", "dplyr", "Seurat", "Matrix", "limma", "speckle", "tidyr", "patchwork"))

# [Removed from the public code export: hard-coded lab phenotyping results. Read from a private file instead.]
LAB_FINDINGS <- readr::read_csv(file.path(PROJECT_ROOT, "config", "lab_findings.csv"), show_col_types = FALSE)
readr::write_csv(LAB_FINDINGS, file.path(path_out, "hfpef_lab_findings_reported.csv"))
knitr::kable(LAB_FINDINGS |> mutate(reported = md_escape(reported)), row.names = FALSE,
             caption = "Reported findings of the phenotyping of the same model (manuscript Figures 1-3; larger cohorts, sequenced animals not identifiable). Stars as reported: */**/***/**** pairwise p < 0.05/0.01/0.001/0.0001.")

atlas <- readRDS(path_atlas)
atlas_meta <- atlas@meta.data |>
  transmute(barcode = colnames(atlas), sample = as.character(sample), cell_type = as.character(cell_type),
            age = as.character(age), condition = as.character(condition), sex = as.character(sex), seq = as.character(seq))
atlas_counts <- GetAssayData(atlas, assay = "RNA", layer = "counts")
rm(atlas); invisible(gc())
animal_meta <- atlas_meta |> distinct(sample, age, condition, sex, seq) |>
  mutate(group = make_group(age, condition),
         qc_flag = unname(coalesce(sub(":.*", "", QC_FLAGGED_ANIMALS[sample]), "pass")),
         note = unname(coalesce(ANIMAL_NOTES[sample], "")),
         animal_class = factor(case_when(qc_flag != "pass" ~ "flag", nzchar(note) ~ "note", TRUE ~ "pass"), levels = names(ANIMAL_CLASS_SHAPE))) |>
  arrange(match(as.character(group), GROUP_ORDER), sample)
stopifnot(!anyDuplicated(animal_meta$sample), nrow(animal_meta) == 16)
cm_clean <- readr::read_csv(path_cm_barcodes, show_col_types = FALSE) |> filter(CM_kind == "subtype")
ec_clean <- readr::read_csv(path_ec_barcodes, show_col_types = FALSE) |> filter(EC_kind %in% c("subtype", "technical"))
barcode_check <- tibble(
  lineage = c("Cardiomyocyte", "Endothelial"),
  atlas_nuclei = c(sum(atlas_meta$cell_type == "Cardiomyocyte"), sum(atlas_meta$cell_type == "Endothelial")),
  clean_object_nuclei = c(nrow(cm_clean), nrow(ec_clean)),
  identical_barcodes = c(setequal(atlas_meta$barcode[atlas_meta$cell_type == "Cardiomyocyte"], cm_clean$barcode),
                         setequal(atlas_meta$barcode[atlas_meta$cell_type == "Endothelial"], ec_clean$barcode)))
stopifnot(all(barcode_check$identical_barcodes))
knitr::kable(barcode_check, caption = "Atlas lineages against the clean objects of notebooks 04 (CM, kind = subtype) and 03 (EC, subtypes with merged technical strata).")
lineage_counts <- as.data.frame.matrix(table(factor(atlas_meta$sample, levels = animal_meta$sample), factor(atlas_meta$cell_type, levels = LINEAGES)))
animal_meta |> select(sample, group, sex, seq, qc_flag, note) |>
  bind_cols(lineage_counts[animal_meta$sample, ], all_nuclei = as.vector(table(atlas_meta$sample)[animal_meta$sample])) |>
  knitr::kable(row.names = FALSE, caption = sprintf("Animals and nuclei per lineage (atlas); a lineage profile with fewer than %d nuclei leaves that animal out of that lineage's models only.", MIN_PB_NUCLEI))

score_lineage <- function(lineage_name) {
  in_lineage   <- atlas_meta$cell_type == lineage_name
  counts_l     <- atlas_counts[, in_lineage]
  obj <- NormalizeData(CreateSeuratObject(counts = counts_l), normalization.method = "LogNormalize", scale.factor = SCALE_FACTOR, verbose = FALSE)
  plan <- LINEAGE_PLAN[[lineage_name]]
  candidate_genes <- intersect(unique(c(unlist(GENE_SETS[plan$modules]), plan$genes)), rownames(counts_l))
  nuclei_detecting <- Matrix::rowSums(counts_l[candidate_genes, , drop = FALSE] > 0)
  set_genes <- lapply(GENE_SETS[plan$modules], function(g) g[g %in% names(nuclei_detecting)[nuclei_detecting > 0]])
  obj <- AddModuleScore(obj, features = set_genes, name = "module_", ctrl = MODULE_CTRL, seed = SEED)
  module_mat <- t(as.matrix(obj@meta.data[, paste0("module_", seq_along(set_genes)), drop = FALSE]))
  rownames(module_mat) <- paste("module", names(set_genes), sep = "|")
  gene_ids <- intersect(plan$genes, rownames(obj))
  gene_mat <- GetAssayData(obj, layer = "data")[gene_ids, , drop = FALSE]
  indicator <- Matrix::t(Matrix::fac2sparse(factor(atlas_meta$sample[in_lineage], levels = animal_meta$sample)))
  n_nuclei  <- Matrix::colSums(indicator)
  means     <- sweep(rbind(as.matrix(module_mat %*% indicator), as.matrix(gene_mat %*% indicator) |> `rownames<-`(paste("gene", gene_ids, sep = "|"))),
                     2, n_nuclei, "/")
  detection <- sweep(as.matrix((counts_l[gene_ids, , drop = FALSE] > 0) %*% indicator), 2, n_nuclei, "/")
  list(values = as.data.frame(as.table(means), stringsAsFactors = FALSE) |> setNames(c("key", "sample", "value")) |>
         separate(key, c("kind", "feature"), sep = "\\|") |> mutate(lineage = lineage_name, n_nuclei = n_nuclei[sample]),
       set_genes = set_genes,
       detection = tibble(lineage = lineage_name, gene = gene_ids, nuclei_detecting = unname(nuclei_detecting[gene_ids]),
                          detected_prop = nuclei_detecting / sum(in_lineage),
                          animals_detecting = rowSums(detection > 0, na.rm = TRUE), min_animal_detected_prop = apply(detection, 1, min, na.rm = TRUE),
                          max_animal_detected_prop = apply(detection, 1, max, na.rm = TRUE)))
}
lineage_scores <- lapply(setNames(LINEAGES, LINEAGES), score_lineage)
broad_props <- speckle::getTransformedProps(clusters = atlas_meta$cell_type, sample = atlas_meta$sample, transform = "logit")
broad_composition <- fit_composition(broad_props, animal_meta, covariates = DESIGN_COVARIATES)
rm(atlas_counts); invisible(gc())
set_sizes <- bind_rows(lapply(LINEAGES, function(l) bind_rows(lapply(names(lineage_scores[[l]]$set_genes), function(m)
  tibble(module = m, lineage = l, genes_in_set = length(GENE_SETS[[m]]), genes_scored = length(lineage_scores[[l]]$set_genes[[m]]),
         genes = paste(lineage_scores[[l]]$set_genes[[m]], collapse = ", "),
         not_detected = paste(setdiff(GENE_SETS[[m]], lineage_scores[[l]]$set_genes[[m]]), collapse = ", "))))))
readr::write_csv(set_sizes, file.path(path_out, "hfpef_gene_set_sizes.csv"))
gene_detection <- bind_rows(lapply(lineage_scores, `[[`, "detection"))
readr::write_csv(gene_detection, file.path(path_out, "hfpef_gene_detection.csv"))
knitr::kable(set_sizes, caption = "Module gene sets: genes in the set, genes detected in at least one nucleus of the lineage (scored) and genes without any count.")

animal_values <- bind_rows(lapply(lineage_scores, `[[`, "values")) |>
  bind_rows(tibble(sample = colnames(broad_props$Proportions), value = as.vector(broad_props$Proportions["Fibroblast", ]),
                   kind = "composition", feature = "fibroblast proportion", lineage = "atlas",
                   n_nuclei = as.vector(colSums(broad_props$Counts)))) |>
  bind_rows(readr::read_csv(path_cm_props, show_col_types = FALSE) |> filter(CM_subtype == NPPA_CM) |>
              transmute(sample, value = proportion, kind = "composition", feature = "Nppa-high CM share", lineage = "Cardiomyocyte", n_nuclei)) |>
  mutate(in_model = n_nuclei >= MIN_PB_NUCLEI) |>
  left_join(animal_meta, by = "sample") |>
  left_join(READOUTS |> select(lineage, kind, feature, section, label), by = c("lineage", "kind", "feature"))
readr::write_csv(animal_values |> transmute(sample, group, sex, seq, qc_flag, note, lineage, kind, feature, readout = label, n_nuclei, in_model,
                                            value_unit = case_when(kind == "module" ~ "module score", kind == "gene" ~ "mean log-normalised expression",
                                                                   TRUE ~ "proportion (0-1)"), value),
                 file.path(path_out, "hfpef_per_animal_values.csv"))
below_min <- animal_values |> distinct(sample, lineage, n_nuclei, in_model) |> filter(!in_model)

fit_animal_model <- function(values, profiles) {
  design_matrix <- model.matrix(DESIGN_FORMULA, data = profiles)
  colnames(design_matrix) <- sub("^group", "", colnames(design_matrix))
  fit <- contrasts.fit(lmFit(values, design_matrix), make_contrast_matrix(design_matrix))
  std_error <- fit$stdev.unscaled * fit$sigma
  bind_rows(lapply(colnames(fit$coefficients), function(contrast_name) {
    t_stat <- fit$coefficients[, contrast_name] / std_error[, contrast_name]
    data.frame(key = rownames(values), contrast = contrast_name, effect = fit$coefficients[, contrast_name],
               SE = std_error[, contrast_name], t = t_stat, p_value = 2 * pt(-abs(t_stat), fit$df.residual),
               CI_low = fit$coefficients[, contrast_name] - qt(.975, fit$df.residual) * std_error[, contrast_name],
               CI_high = fit$coefficients[, contrast_name] + qt(.975, fit$df.residual) * std_error[, contrast_name],
               n_animals = nrow(profiles), row.names = NULL)
  }))
}
tested_values <- animal_values |> filter(!is.na(label), kind != "composition", in_model)
expression_res <- bind_rows(lapply(split(tested_values, tested_values$lineage), function(d) {
  profiles <- d |> distinct(sample, group, sex, seq) |> arrange(sample)
  value_mat <- d |> select(label, sample, value) |> pivot_wider(names_from = sample, values_from = value) |>
    tibble::column_to_rownames("label") |> as.matrix()
  fit_animal_model(value_mat[, profiles$sample, drop = FALSE], profiles)
})) |> dplyr::rename(label = key) |> mutate(FDR_source = NA_real_, source = "limma, this notebook")
nppa_composition <- readr::read_csv(path_cm_propeller, show_col_types = FALSE) |> filter(subtype == NPPA_CM)
composition_res <- bind_rows(
  broad_composition |> filter(subtype == "Fibroblast") |> mutate(label = "fibroblast proportion (atlas)", source = "fit_composition, this notebook"),
  nppa_composition |> mutate(label = "Nppa-high CM share (CM)", source = "notebook 04 cm_propeller_results.csv")) |>
  transmute(label, contrast = as.character(contrast), effect = logit_effect, SE, t = logit_effect / SE, p_value, CI_low, CI_high,
            n_animals, FDR_source = FDR, source)
animal_res <- bind_rows(expression_res, composition_res) |>
  group_by(contrast) |> mutate(FDR = p.adjust(p_value, "BH")) |> ungroup() |>
  left_join(READOUTS, by = "label") |>
  mutate(contrast = factor(contrast, levels = CONTRAST_ORDER), label = factor(label, levels = READOUTS$label)) |>
  arrange(label, contrast)
stopifnot(nrow(animal_res) == nrow(READOUTS) * length(CONTRAST_ORDER))
readr::write_csv(animal_res |> select(section, label, lineage, kind, feature, unit, contrast, effect, SE, CI_low, CI_high, t, p_value, FDR, FDR_source, source, n_animals),
                 file.path(path_out, "hfpef_animal_model.csv"))
group_means <- animal_values |> filter(!is.na(label), in_model) |> group_by(label, group) |> summarise(mean = mean(value), .groups = "drop") |>
  pivot_wider(names_from = group, values_from = mean)
res_of   <- function(lab, cn) filter(animal_res, label == lab, contrast == cn)
res_line <- function(lab, cn) with(res_of(lab, cn), sprintf("%+.3g (95%% CI %.3g to %.3g, FDR = %s)", effect, CI_low, CI_high, fmt_fdr(FDR)))
gm_text  <- function(lab, digits = 3) with(filter(group_means, label == lab), sprintf(paste0("YC %.", digits, "f, OC %.", digits, "f, YH %.", digits, "f, OH %.", digits, "f"), YC, OC, YH, OH))
hits_text <- function(sections) {
  hits <- filter(animal_res, section %in% sections, FDR < FDR_CUT)
  if (!nrow(hits)) return("none")
  paste(sprintf("%s, %s: %+.3g, FDR = %s", hits$label, CONTRAST_LABEL[as.character(hits$contrast)], hits$effect, fmt_fdr(hits$FDR)), collapse = "; ")
}
nominal_text <- function(sections) {
  nominal <- filter(animal_res, section %in% sections, p_value < 0.05)
  if (!nrow(nominal)) return("none")
  paste(sprintf("%s %s (%+.3g, p = %s)", nominal$label, CONTRAST_CODE[as.character(nominal$contrast)], nominal$effect, fmt_fdr(nominal$p_value)), collapse = "; ")
}
results_table <- function(sections, caption) {
  filter(animal_res, section %in% sections) |>
    transmute(readout = as.character(label), contrast = CONTRAST_LABEL[as.character(contrast)], effect = signif(effect, 3),
              `95% CI` = sprintf("%.3g to %.3g", CI_low, CI_high), p = signif(p_value, 2), FDR = signif(FDR, 2), n_animals) |>
    knitr::kable(row.names = FALSE, caption = caption)
}
contrast_summary <- animal_res |> group_by(contrast) |>
  summarise(tests = n(), FDR05 = sum(FDR < FDR_CUT), nominal05 = sum(p_value < 0.05), min_FDR = min(FDR), .groups = "drop")
knitr::kable(contrast_summary |> mutate(contrast = CONTRAST_LABEL[as.character(contrast)], min_FDR = signif(min_FDR, 2)),
             caption = "Animal-level tests per contrast: readouts at FDR < 0.05 (BH within contrast over all readouts) and at nominal p < 0.05.")

animal_point_layers <- function() list(
  stat_summary(aes(group = group), fun = mean, fun.min = mean, fun.max = mean, geom = "errorbar", width = .45, linewidth = .35, colour = "grey45"),
  geom_point(aes(fill = group, colour = group, shape = animal_class), size = 2.1, stroke = .8,
             position = position_jitter(width = .12, height = 0, seed = SEED)),
  scale_fill_manual(values = GROUP_FILL, labels = GROUP_LABEL, name = NULL, drop = FALSE),
  scale_colour_manual(values = GROUP_FILL, labels = GROUP_LABEL, name = NULL, drop = FALSE),
  scale_shape_manual(values = ANIMAL_CLASS_SHAPE, labels = ANIMAL_CLASS_LABEL, name = NULL, drop = FALSE),
  scale_x_discrete(limits = GROUP_ORDER),
  guides(fill = guide_legend(override.aes = list(shape = 21)), shape = guide_legend(override.aes = list(fill = "grey60", colour = "grey30"))),
  theme_fig, theme(legend.position = "bottom", legend.box = "vertical", strip.background = element_blank(), strip.text = element_text(size = 7.5, face = "plain")))
strip_text <- function(labels) vapply(labels, function(lab) {
  hits <- filter(animal_res, label == lab, FDR < FDR_CUT)
  paste0(lab, "\n", if (nrow(hits)) paste(sprintf("%s FDR %s", CONTRAST_CODE[as.character(hits$contrast)], fmt_fdr(hits$FDR)), collapse = "; ") else "no FDR < 0.05")
}, "")
per_animal_plot <- function(labels, title, ncol, y_lab = NULL) {
  plot_df <- animal_values |> filter(label %in% labels, in_model) |>
    mutate(panel = factor(strip_text(as.character(label)), levels = strip_text(labels)), group = factor(as.character(group), levels = GROUP_ORDER))
  ggplot(plot_df, aes(group, value)) + animal_point_layers() + facet_wrap(~ panel, scales = "free_y", ncol = ncol) +
    labs(x = NULL, y = y_lab, title = title,
         subtitle = paste0("Points = animals (profiles with >= ", MIN_PB_NUCLEI, " nuclei), bars = group means. Strip: contrasts at FDR < 0.05 (BH within contrast over all ",
                           nrow(READOUTS), " readouts).\n", contrast_key))
}
OVERVIEW <- READOUTS$label[READOUTS$kind != "gene" | READOUTS$feature == "Nos2"]
fig <- per_animal_plot(OVERVIEW, "Molecular readouts of the HFpEF phenotype per animal",
                       ncol = 4, y_lab = "module score | proportion (0-1) | mean log-normalised expression (Nos2)")
save_fig(fig, "hfpef_readouts_per_animal", width = 13, height = 10)
fig

results_table("S1", "Readout 1, animal-level model. Effect: module score / mean log-normalised expression (natural log) / logit proportion difference (pairwise) or difference of differences (Age x HFpEF); FDR: BH within contrast over all readouts.")
nppa_share <- filter(group_means, label == "Nppa-high CM share (CM)")

fig <- per_animal_plot(READOUTS$label[READOUTS$section %in% c("S1", "S2") & READOUTS$kind == "gene"],
                       "Fetal-gene and calcium-handling genes in cardiomyocytes per animal", ncol = 6,
                       y_lab = "mean log-normalised expression (natural log)")
save_fig(fig, "hfpef_cm_genes_per_animal", width = 13, height = 8)
fig
cm_detection <- gene_detection |> filter(lineage == "Cardiomyocyte", gene %in% c(GENE_SETS$hypertrophy, GENE_SETS$calcium))

results_table("S2", "Readout 2, animal-level model. Effect: module score / mean log-normalised expression difference (pairwise) or difference of differences; FDR: BH within contrast over all readouts.")

results_table("S3", "Readout 3, animal-level model. Effect: module score difference (pairwise) or difference of differences; FDR: BH within contrast over all readouts.")

results_table("S4", "Readout 4, animal-level model. Effect: module score difference or logit-proportion difference (pairwise) or difference of differences; FDR: BH within contrast over all readouts.")
fb_prop <- filter(group_means, label == "fibroblast proportion (atlas)")

results_table("S5", "Readout 5, animal-level model. Effect: module score / mean log-normalised expression difference (pairwise) or difference of differences; FDR: BH within contrast over all readouts.")
sen_detection <- gene_detection |> filter(gene %in% GENE_SETS$senescence) |> arrange(gene, match(lineage, LINEAGES))
knitr::kable(sen_detection |> transmute(gene, lineage, nuclei_detecting, `detected (proportion of nuclei)` = signif(detected_prop, 2), animals_detecting,
                                        `per-animal detected (proportion), min-max` = sprintf("%.3f-%.3f", min_animal_detected_prop, max_animal_detected_prop)),
             row.names = FALSE, caption = "Detection of the senescence genes per lineage.")

fig <- per_animal_plot(as.vector(outer(GENE_SETS$senescence, paste0(" (", CELL_TYPE_SHORT[LINEAGES], ")"), paste0)) |> intersect(READOUTS$label),
                       "Senescence genes per lineage and animal", ncol = 3, y_lab = "mean log-normalised expression (natural log)")
save_fig(fig, "hfpef_senescence_genes_per_animal", width = 12, height = 9)
fig
cdkn2a_det <- filter(sen_detection, gene == "Cdkn2a")

EC_SUBTYPE_ORDER <- names(sort(table(ec_clean$EC_label |> replace(ec_clean$EC_kind == "technical", "Capillary EC")), decreasing = TRUE))
mast_hits     <- readr::read_csv(path_mast_csv, show_col_types = FALSE) |> filter(gene %in% EC_ACT_GENES)
mast_universe <- readr::read_csv(path_universe_csv, show_col_types = FALSE) |> filter(gene %in% EC_ACT_GENES)
stopifnot(setequal(unique(mast_hits$celltype), EC_SUBTYPE_ORDER))
flag_mark <- function(sig, sad, low, bkg) paste0(ifelse(sig & sad %in% TRUE, "\u2020", ""), ifelse(sig & low %in% TRUE, "\u00b0", ""),
                                               ifelse(sig & bkg %in% TRUE, "\u2021", ""))
ec_mast <- expand_grid(celltype = EC_SUBTYPE_ORDER, contrast = CONTRAST_ORDER, gene = EC_ACT_GENES) |>
  left_join(mast_hits |> select(celltype, contrast, gene, avg_log2FC, p_val_adj, pct.1, pct.2, sig, reported, single_animal_driven,
                                low_detection, background_suspect, n_animals_detected_i1, n_animals_detected_i2, max_animal_share, max_share_animal),
            by = c("celltype", "contrast", "gene")) |>
  left_join(mast_universe |> select(celltype, contrast, gene, universe_effect = avg_log2FC, universe_p_adj = p_val_adj), by = c("celltype", "contrast", "gene")) |>
  mutate(status = case_when(!is.na(p_val_adj) ~ "hit table (>= 5% detection)", !is.na(universe_p_adj) ~ "ranking universe only (1-5% detection)",
                            TRUE ~ "not tested (< 1% detection)"),
         effect_scale = ifelse(contrast == "Interaction", "interaction_log2", "avg_log2FC"),
         mark = paste0(ifelse(sig %in% TRUE, as.character(fdr_stars(p_val_adj)), ""),
                       flag_mark(sig %in% TRUE, single_animal_driven, low_detection, background_suspect)))
readr::write_csv(ec_mast |> dplyr::rename(det_i1_prop = pct.1, det_i2_prop = pct.2, max_animal_share_prop = max_animal_share),
                 file.path(path_out, "hfpef_ec_activation_mast.csv"))
ec_mast_summary <- ec_mast |> filter(!is.na(avg_log2FC)) |> group_by(gene, contrast) |>
  summarise(n_tested = n(), median_effect = median(avg_log2FC), n_lower = sum(reported & avg_log2FC < 0), n_higher = sum(reported & avg_log2FC > 0),
            lower_subtypes = paste(short_subtype(celltype[reported & avg_log2FC < 0]), collapse = ", "),
            higher_subtypes = paste(short_subtype(celltype[reported & avg_log2FC > 0]), collapse = ", "),
            min_p_lower = suppressWarnings(min(p_val_adj[avg_log2FC < 0])), min_p_higher = suppressWarnings(min(p_val_adj[avg_log2FC > 0])),
            n_flagged = sum(reported & (single_animal_driven | low_detection | background_suspect)), .groups = "drop") |>
  complete(gene = EC_ACT_GENES, contrast = CONTRAST_ORDER, fill = list(n_tested = 0L, n_lower = 0L, n_higher = 0L, n_flagged = 0L)) |>
  mutate(gene = factor(gene, levels = EC_ACT_GENES), contrast = factor(contrast, levels = CONTRAST_ORDER)) |> arrange(gene, contrast)
readr::write_csv(ec_mast_summary, file.path(path_out, "hfpef_ec_activation_mast_summary.csv"))
ec_mast_summary |>
  transmute(gene, contrast = CONTRAST_LABEL[as.character(contrast)], `subtypes tested` = n_tested, `median effect` = round(median_effect, 3),
            `reported lower` = ifelse(n_lower > 0, sprintf("%d (%s)", n_lower, lower_subtypes), "0"),
            `reported higher` = ifelse(n_higher > 0, sprintf("%d (%s)", n_higher, higher_subtypes), "0"), `reported, flagged` = n_flagged) |>
  knitr::kable(row.names = FALSE, caption = sprintf("Endothelial activation and NO genes in the MAST screen of notebook 05, per gene x contrast over the %d subtypes. Median effect over the subtypes in the hit table (avg_log2FC; Age x HFpEF interaction_log2); reported = Bonferroni p < 0.05, |effect| > 0.2, non-technical; flagged = reported and single-animal-driven, low-detection or background-suspect.", length(EC_SUBTYPE_ORDER)))
mast_count <- function(g, cn) with(filter(ec_mast_summary, gene == g, contrast == cn),
  if (n_tested == 0) "not in the hit table" else sprintf("median %+.2f over %d subtypes; reported lower in %d%s, higher in %d%s", median_effect, n_tested,
                                                         n_lower, ifelse(n_lower > 0, paste0(" (", lower_subtypes, ")"), ""), n_higher, ifelse(n_higher > 0, paste0(" (", higher_subtypes, ")"), "")))

heatmap_df <- ec_mast |>
  mutate(gene = factor(gene, levels = EC_ACT_GENES), celltype = factor(celltype, levels = rev(EC_SUBTYPE_ORDER)),
         contrast = factor(contrast, levels = CONTRAST_ORDER), effect_capped = pmax(pmin(avg_log2FC, 1.2), -1.2),
         label = ifelse(is.na(avg_log2FC), "\u2013", mark))
fig <- ggplot(heatmap_df, aes(contrast, celltype, fill = effect_capped)) +
  geom_tile(colour = "white", linewidth = .6) +
  geom_text(aes(label = label, colour = is.na(avg_log2FC)), size = 3, show.legend = FALSE) +
  facet_wrap(~ gene, nrow = 1) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, limits = c(-1.2, 1.2), na.value = "grey90",
                       breaks = c(-1.2, -.6, 0, .6, 1.2), name = "effect\navg_log2FC |\ninteraction_log2") +
  scale_colour_manual(values = c(`FALSE` = "grey15", `TRUE` = "grey55")) +
  scale_x_discrete(labels = CONTRAST_CODE) + scale_y_discrete(labels = short_subtype) +
  labs(x = NULL, y = NULL, title = "Endothelial activation and NO genes in the MAST screen (notebook 05)",
       subtitle = paste0("*/**/*** Bonferroni p < 0.05/0.01/0.001 with |effect| > ", LFC_CUT, "; \u2020 single-animal-driven, \u00b0 low detection, \u2021 background suspect; ",
                         "grey dash = below the 5% detection floor, not in the hit table.\n", contrast_key,
                         ". Pairwise fill = avg_log2FC; Age x HFpEF fill = interaction_log2 (positive = larger HFpEF effect in old).")) +
  theme_fig + theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_text(size = 7),
                    strip.background = element_blank(), strip.text = element_text(face = "bold.italic"))
save_fig(fig, "hfpef_ec_activation_mast_heatmap", width = 13, height = 4.6)
fig

ec_codes <- function(g) with(filter(ec_mast_summary, gene == g, n_lower + n_higher > 0),
  if (length(gene)) paste(sprintf("%s %s%s", CONTRAST_CODE[as.character(contrast)], ifelse(n_lower > 0, paste0("\u2193", n_lower), ""),
                                  ifelse(n_higher > 0, paste0("\u2191", n_higher), "")), collapse = "; ") else "none reported")
ec_plot_df <- animal_values |> filter(lineage == "Endothelial", kind == "gene", feature %in% EC_ACT_GENES, in_model) |>
  mutate(panel = factor(sprintf("%s (all EC)\nMAST: %s", feature, vapply(feature, ec_codes, "")),
                        levels = sprintf("%s (all EC)\nMAST: %s", EC_ACT_GENES, vapply(EC_ACT_GENES, ec_codes, ""))),
         group = factor(as.character(group), levels = GROUP_ORDER))
fig <- ggplot(ec_plot_df, aes(group, value)) + animal_point_layers() + facet_wrap(~ panel, scales = "free_y", nrow = 1) +
  labs(x = NULL, y = "mean log-normalised expression (natural log)", title = "Endothelial activation and NO genes per animal (all endothelial nuclei)",
       subtitle = paste0("Points = animals, bars = group means; shown for visualisation, tested by MAST per subtype (notebook 05). Strip: number of subtypes reported lower (\u2193) / higher (\u2191) per contrast.\n", contrast_key))
save_fig(fig, "hfpef_ec_activation_per_animal", width = 13, height = 4.8)
fig
ec_gm <- ec_plot_df |> group_by(feature, group) |> summarise(mean = mean(value), .groups = "drop") |> pivot_wider(names_from = group, values_from = mean)

results_table("S7", "Readout 7, animal-level model. Effect: module score / mean log-normalised expression difference (pairwise) or difference of differences; FDR: BH within contrast over all readouts.")
infl_detection <- gene_detection |> filter((lineage == "Macrophage" & gene %in% c("Nos2", GENE_SETS$inflammatory)) | (lineage == "Cardiomyocyte" & gene %in% GENE_SETS$upr))

infl_genes <- animal_values |> filter(kind == "gene", in_model,
                                      (lineage == "Macrophage" & feature %in% GENE_SETS$inflammatory) | (lineage == "Cardiomyocyte" & feature %in% GENE_SETS$upr)) |>
  mutate(panel = factor(paste0(feature, " (", CELL_TYPE_SHORT[lineage], ")\nshown, not tested"),
                        levels = paste0(c(GENE_SETS$inflammatory, GENE_SETS$upr), " (", rep(c("Mac", "CM"), c(5, 4)), ")\nshown, not tested")),
         group = factor(as.character(group), levels = GROUP_ORDER))
fig <- ggplot(infl_genes, aes(group, value)) + animal_point_layers() + facet_wrap(~ panel, scales = "free_y", nrow = 2) +
  labs(x = NULL, y = "mean log-normalised expression (natural log)", title = "Inflammatory macrophage genes and UPR genes in cardiomyocytes per animal",
       subtitle = paste0("Points = animals, bars = group means; the genes enter the tested inflammatory and UPR modules (overview figure).\n", contrast_key))
save_fig(fig, "hfpef_inflammation_upr_genes_per_animal", width = 13, height = 6.5)
fig

# [Removed from the public code export: hard-coded lab phenotyping results. Read from a private file instead.]
CONCORDANCE_MAP <- readr::read_csv(file.path(PROJECT_ROOT, "config", "lab_concordance_map.csv"), show_col_types = FALSE)
animal_cells <- animal_res |> transmute(readout = as.character(label), contrast = as.character(contrast), effect, stat = FDR, sig = FDR < FDR_CUT,
                                        method = ifelse(kind == "composition", "logit-proportion model", "limma, per-animal means"),
                                        cell_text = sprintf("%+.2g\nFDR %s", effect, fmt_fdr(FDR)))

mast_cells <- ec_mast_summary |>
  transmute(readout = paste0(gene, " (EC, MAST)"), contrast = as.character(contrast), effect = median_effect,
            n_lower, n_higher, min_p_lower, min_p_higher, method = "MAST (notebook 05), per subtype",
            cell_text = sprintf("%+.2f\n%d\u2193 %d\u2191 /%d", median_effect, n_lower, n_higher, n_tested))
concordance_cells <- bind_rows(animal_cells, mast_cells)
concordance <- CONCORDANCE_MAP |> mutate(row_id = row_number()) |>
  left_join(concordance_cells, by = "readout", relationship = "many-to-many") |>
  mutate(lab_contrast = mapply(function(cs, cn) !is.na(cn) && cn %in% strsplit(cs, ";")[[1]], lab_contrasts, contrast),
         lab_sign = c(up = 1, down = -1)[lab_direction],
         is_mast = grepl("MAST", method),
         mast_side = case_when(!is_mast ~ NA_real_, lab_contrast ~ lab_sign, n_higher > 0 & n_lower == 0 ~ 1, n_lower > 0 & n_higher == 0 ~ -1,
                               TRUE ~ sign(effect)),
         stat = ifelse(is_mast, ifelse(mast_side > 0, min_p_higher, min_p_lower), stat),
         sig  = ifelse(is_mast, ifelse(mast_side > 0, n_higher > 0, n_lower > 0), sig),
         opposite_mast = is_mast & lab_contrast & ifelse(lab_sign > 0, n_lower > 0, n_higher > 0) %in% TRUE,
         effect_dir = ifelse(is_mast & sig %in% TRUE, mast_side, sign(effect)),
         same_direction = coalesce(lab_contrast & sign(effect) == lab_sign, FALSE),
         same_direction_sig = same_direction & sig %in% TRUE,
         opposite_sig = ifelse(is_mast, opposite_mast, coalesce(lab_contrast & sign(effect) == -lab_sign & sig, FALSE)))
concordance_table <- concordance |> group_by(row_id, phenotype, lab_finding, basis, lab_direction, lab_contrasts, readout) |>
  summarise(method = first(method),
            per_contrast = list(setNames(ifelse(is.na(effect), NA, sprintf("%+.3g (%s %s)%s", effect, ifelse(is_mast, paste0("min p_Bonf, ", ifelse(mast_side > 0, "higher", "lower"), " side"), "FDR"),
                                                                          ifelse(is.finite(stat), fmt_fdr(stat), "none"),
                                                                          ifelse(sig %in% TRUE, " sig", ""))), contrast)),
            n_lab_contrasts = sum(lab_contrast), concordant_direction = sum(same_direction, na.rm = TRUE), significant = sum(same_direction_sig, na.rm = TRUE),
            opposite_significant = sum(opposite_sig), .groups = "drop") |>
  mutate(concordant_direction_yes_no = case_when(is.na(readout) | n_lab_contrasts == 0 ~ "n/a", concordant_direction == n_lab_contrasts ~ "yes",
                                                 concordant_direction > 0 ~ "partial", TRUE ~ "no"),
         significant_yes_no = case_when(is.na(readout) | n_lab_contrasts == 0 ~ "n/a", significant == n_lab_contrasts ~ "yes",
                                        significant > 0 ~ "partial", TRUE ~ "no"),
         concordant_direction = ifelse(n_lab_contrasts > 0, sprintf("%d/%d", concordant_direction, n_lab_contrasts), "n/a"),
         significant = ifelse(n_lab_contrasts > 0, sprintf("%d/%d", significant, n_lab_contrasts), "n/a"),
         opposite_significant = ifelse(n_lab_contrasts > 0, sprintf("%d/%d", opposite_significant, n_lab_contrasts), "n/a"))
for (cn in CONTRAST_ORDER) concordance_table[[cn]] <- vapply(concordance_table$per_contrast, function(x) if (is.null(x) || !cn %in% names(x)) NA_character_ else unname(x[cn]), "")
concordance_table <- concordance_table |> select(-per_contrast) |> arrange(row_id)
readr::write_csv(concordance_table |> select(-row_id), file.path(path_out, "hfpef_concordance_table.csv"))
concordance_table |>
  mutate(across(all_of(CONTRAST_ORDER), \(x) coalesce(x, "-"))) |>
  transmute(phenotype, `lab finding` = md_escape(lab_finding), basis, readout = coalesce(readout, "none"),
            `lab contrasts` = vapply(strsplit(lab_contrasts, ";"), function(x) paste(CONTRAST_CODE[x], collapse = ", "), ""),
            across(all_of(CONTRAST_ORDER)),
            `concordant direction` = paste0(concordant_direction_yes_no, " (", concordant_direction, ")"), significant = paste0(significant_yes_no, " (", significant, ")"),
            `opposite, significant` = opposite_significant) |>
  rename_with(\(x) CONTRAST_CODE[x], all_of(CONTRAST_ORDER)) |>
  knitr::kable(row.names = FALSE, caption = paste("Concordance of the snRNA readouts with the phenotyping of the same model (manuscript Figures 1-3). Per contrast: effect (FDR, BH within contrast over all readouts; MAST: median avg_log2FC over subtypes and the smallest Bonferroni p on the side read - the lab side in a lab contrast, otherwise the side holding a reported subtype (or the median side); 'sig' = passes its threshold on that side. Concordant direction / significant / opposite, significant: number of lab-significant (or model-expected) contrasts in which the snRNA effect has the lab direction / has it and passes its threshold / has the opposite direction and passes its threshold.", contrast_key))

CAT_FILL <- c("higher, significant" = "#b2182b", "higher, n.s." = "#f4c6bd", "no change (effect 0)" = "grey85", "lower, n.s." = "#c3d8ea", "lower, significant" = "#2166ac", "not tested" = "grey96")
fig_df <- concordance |> filter(!is.na(readout)) |>
  mutate(effect_dir = ifelse(opposite_sig, -lab_sign, effect_dir), sig = sig %in% TRUE | opposite_sig,
         category = factor(case_when(is.na(effect) ~ "not tested", sig & effect_dir > 0 ~ "higher, significant", sig & effect_dir < 0 ~ "lower, significant",
                                     effect > 0 ~ "higher, n.s.", effect == 0 ~ "no change (effect 0)", TRUE ~ "lower, n.s."), levels = names(CAT_FILL)),
         contrast = factor(contrast, levels = CONTRAST_ORDER),
         row_label = paste0(readout, "  [", ifelse(nchar(lab_finding) > 46, paste0(substr(lab_finding, 1, 45), "\u2026"), lab_finding), "]"),
         row_label = factor(row_label, levels = rev(unique(row_label))),
         phenotype = factor(phenotype, levels = unique(CONCORDANCE_MAP$phenotype)),
         lab_mark = ifelse(lab_contrast, ifelse(lab_direction == "up", "\u25b2", "\u25bc"), ""))
agree_df <- concordance_table |> filter(!is.na(readout)) |>
  transmute(row_id, text = ifelse(concordant_direction == "n/a", "n/a", sprintf("dir %s | sig %s | opp %s", concordant_direction, significant, opposite_significant))) |>
  left_join(distinct(fig_df, row_id, row_label, phenotype), by = "row_id")
fig <- ggplot(fig_df, aes(contrast, row_label)) +
  geom_tile(aes(fill = category, colour = lab_contrast, linewidth = lab_contrast), width = .94, height = .9) +
  geom_text(aes(label = cell_text), size = 2.1, lineheight = .85, colour = "grey10") +
  geom_text(aes(label = lab_mark), size = 2.2, nudge_x = -.4, nudge_y = .25, colour = "grey10") +
  geom_text(data = agree_df, aes(x = 6, y = row_label, label = text), inherit.aes = FALSE, size = 2.5, hjust = 0, colour = "grey20") +
  facet_grid(phenotype ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_manual(values = CAT_FILL, name = NULL, drop = FALSE) +
  scale_colour_manual(values = c(`FALSE` = "white", `TRUE` = "grey15"), guide = "none") +
  scale_linewidth_manual(values = c(`FALSE` = .3, `TRUE` = .7), guide = "none") +
  scale_x_discrete(labels = CONTRAST_LABEL_WRAP, expand = expansion(add = c(.6, 1.9))) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = NULL, title = "Concordance of snRNA readouts with the phenotyping of the same model (manuscript Figures 1-3)",
       subtitle = paste0("Cell: effect and FDR (animal level, BH within contrast over all readouts) or MAST median effect over EC subtypes with subtypes reported lower/higher/tested.\n",
                         "Outlined cell with \u25b2/\u25bc = contrast in which the lab reports a significant increase/decrease (or the model-expected direction). ",
                         "Right column: lab contrasts with the lab direction | with it and passing the threshold | opposite and passing the threshold.")) +
  theme_fig + theme(axis.line = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), legend.position = "bottom",
                    strip.placement = "outside", strip.background = element_blank(), strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold", size = 8),
                    axis.text.y = element_text(size = 7.5), panel.spacing.y = unit(2, "pt"), plot.margin = margin(5, 110, 5, 5),
                    plot.title.position = "plot")
save_fig(fig, "hfpef_concordance_table", width = 14, height = 12.5)
fig
conc_rows <- filter(concordance_table, !is.na(readout), concordant_direction != "n/a")
row_list <- function(rows) if (nrow(rows)) paste(sprintf("%s (%s)", rows$readout, rows$significant), collapse = "; ") else "none"

hfpef_cells <- concordance |> filter(lab_contrast, contrast %in% c("YH_vs_YC", "OH_vs_OC"))
age_cells   <- concordance |> filter(lab_contrast, contrast == "OC_vs_YC")
cell_summary <- function(d) sprintf("%d of %d cells with the lab direction, %d passing the threshold", sum(d$same_direction, na.rm = TRUE), nrow(d), sum(d$same_direction_sig, na.rm = TRUE))

session_footer()
