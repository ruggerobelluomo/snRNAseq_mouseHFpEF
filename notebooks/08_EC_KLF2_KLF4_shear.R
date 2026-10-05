# 08_EC_KLF2_KLF4_shear.R
# Code-only export of the analysis pipeline. Data are not included; see README.md.

params <- list(rds_in = "results/03_ec_annotation/Endothelial_subcluster_annotated_clean.rds", 
    mast_csv = "results/05_ec_de/ec_DE_all_MAST_annotated.csv", universe_csv = "results/05_ec_de/ec_DE_all_MAST_gsea_universe.csv", 
    net_csv = "data/resources/collectri_mouse_network.csv", out_dir = "results/08_ec_klf")


knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = TRUE,
                      fig.align = "center", dpi = 150)
options(width = 110)

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(edgeR); library(limma); library(decoupleR); library(tidyr); library(patchwork); library(ggrepel)
})
source("_common.R")
knitr::opts_chunk$set(cache = TRUE, autodep = TRUE, cache.lazy = FALSE, cache.extra = input_md5(params))
path_rds_in       <- resolve_path(params$rds_in)
path_mast_csv     <- resolve_path(params$mast_csv)
path_universe_csv <- resolve_path(params$universe_csv)
path_net_csv      <- resolve_path(params$net_csv)
path_out     <- results_dir("08_ec_klf")
path_figures <- file.path(path_out, "figures"); dir.create(path_figures, recursive = TRUE, showWarnings = FALSE)
format_count <- function(x) formatC(x, big.mark = ",", format = "d")

SEED          <- 1
MIN_PB_NUCLEI <- 25
MIN_TARGETS   <- 5
FDR_CUT       <- 0.05
PADJ_CUT      <- 0.05
LFC_CUT       <- 0.2
MODULE_CTRL   <- 50
ULM_CHUNK     <- 5000

SHEAR_PROG <- c("Klf2", "Klf4", "Nos3", "Thbd", "Nqo1", "Mef2c", "Nppc", "Adamts1", "Ptgis", "S1pr1", "Timp3", "Ackr1")
KLF_GENES  <- c("Klf2", "Klf4")
FOCUS_CT   <- "Capillary EC"
ALL_EC     <- "All EC"
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

print_versions(c("ggplot2", "dplyr", "Seurat", "Matrix", "edgeR", "limma", "decoupleR", "tidyr", "patchwork", "ggrepel"))

ec_obj <- readRDS(path_rds_in)
DefaultAssay(ec_obj) <- "RNA"
SUBTYPE_ORDER <- levels(droplevels(factor(ec_obj$EC_subtype)))
LEVEL_ORDER   <- c(SUBTYPE_ORDER, ALL_EC)
module_detecting <- Matrix::rowSums(GetAssayData(ec_obj, layer = "counts")[intersect(SHEAR_PROG, rownames(ec_obj)), ] > 0)
MODULE_GENES <- names(module_detecting)[module_detecting > 0]
ec_obj <- AddModuleScore(ec_obj, features = list(MODULE_GENES), name = "shear_module", ctrl = MODULE_CTRL, seed = SEED)
nucleus_meta <- data.frame(sample = as.character(ec_obj$sample), EC_subtype = factor(ec_obj$EC_subtype, levels = SUBTYPE_ORDER),
                           module_score = ec_obj$shear_module1)
animal_meta <- ec_obj@meta.data |> distinct(sample, age, condition, sex, seq) |>
  mutate(across(where(is.factor), as.character), group = make_group(age, condition),
         qc_flag = unname(coalesce(sub(":.*", "", QC_FLAGGED_ANIMALS[sample]), "pass")),
         note = unname(coalesce(ANIMAL_NOTES[sample], "")),
         animal_class = factor(case_when(qc_flag != "pass" ~ "flag", nzchar(note) ~ "note", TRUE ~ "pass"), levels = names(ANIMAL_CLASS_SHAPE)))
counts  <- GetAssayData(ec_obj, layer = "counts")
lognorm <- GetAssayData(ec_obj, layer = "data")
rm(ec_obj); invisible(gc())
animal_meta |> arrange(match(as.character(group), GROUP_ORDER), sample) |>
  left_join(count(nucleus_meta, sample, name = "EC_nuclei"), by = "sample") |>
  select(sample, group, age, condition, sex, seq, EC_nuclei, qc_flag, note) |>
  knitr::kable(row.names = FALSE, caption = "Animals (biological replicates); every animal enters every model.")

fit_animal_model <- function(values, profiles, moderated = FALSE) {
  design_matrix <- model.matrix(DESIGN_FORMULA, data = profiles)
  colnames(design_matrix) <- sub("^group", "", colnames(design_matrix))
  fit <- contrasts.fit(lmFit(values, design_matrix), make_contrast_matrix(design_matrix))
  if (moderated) fit <- eBayes(fit, robust = TRUE)
  residual_sd <- if (moderated) sqrt(fit$s2.post) else fit$sigma
  test_df     <- if (moderated) fit$df.total else fit$df.residual
  std_error   <- fit$stdev.unscaled * residual_sd
  bind_rows(lapply(colnames(fit$coefficients), function(contrast_name) {
    t_stat <- fit$coefficients[, contrast_name] / std_error[, contrast_name]
    data.frame(feature = rownames(values), contrast = contrast_name, effect = fit$coefficients[, contrast_name],
               SE = std_error[, contrast_name], t = t_stat, p_value = 2 * pt(-abs(t_stat), test_df),
               CI_low = fit$coefficients[, contrast_name] - qt(.975, test_df) * std_error[, contrast_name],
               CI_high = fit$coefficients[, contrast_name] + qt(.975, test_df) * std_error[, contrast_name],
               n_animals = nrow(profiles), row.names = NULL)
  }))
}

profile_sums <- function(mat, level_of_nucleus)
  mat %*% Matrix::t(Matrix::fac2sparse(factor(paste(nucleus_meta$sample, level_of_nucleus, sep = "|"))))
all_profile_sums <- function(mat) cbind(profile_sums(mat, nucleus_meta$EC_subtype), profile_sums(mat, ALL_EC))
profile_meta <- bind_rows(count(nucleus_meta, sample, level = as.character(EC_subtype), name = "n_nuclei"),
                          count(nucleus_meta, sample, name = "n_nuclei") |> mutate(level = ALL_EC)) |>
  mutate(profile = paste(sample, level, sep = "|"), in_model = n_nuclei >= MIN_PB_NUCLEI) |>
  left_join(animal_meta, by = "sample") |> mutate(level = factor(level, levels = LEVEL_ORDER))
pb_counts <- all_profile_sums(counts)[, profile_meta$profile]
profile_means <- function(mat) sweep(as.matrix(all_profile_sums(mat)[, profile_meta$profile, drop = FALSE]), 2, profile_meta$n_nuclei, "/")
readout_means <- profile_means(rbind(as.matrix(lognorm[KLF_GENES, ]), module = nucleus_meta$module_score))
profile_meta |> transmute(sample, level, nuclei = paste0(format_count(n_nuclei), ifelse(in_model, "", " (below minimum)"))) |>
  pivot_wider(names_from = level, values_from = nuclei) |>
  knitr::kable(caption = sprintf("Nuclei per animal x subtype profile; profiles below %d nuclei do not enter that subtype's models.", MIN_PB_NUCLEI))
below_min <- filter(profile_meta, !in_model)

mast_hits     <- readr::read_csv(path_mast_csv, show_col_types = FALSE) |> filter(gene %in% SHEAR_PROG)
mast_universe <- readr::read_csv(path_universe_csv, show_col_types = FALSE) |> filter(gene %in% SHEAR_PROG)
flag_mark <- function(sig, sad, low, bkg) paste0(ifelse(sig & sad %in% TRUE, "\u2020", ""), ifelse(sig & low %in% TRUE, "\u00b0", ""),
                                               ifelse(sig & bkg %in% TRUE, "\u2021", ""))
module_mast <- expand_grid(celltype = SUBTYPE_ORDER, contrast = CONTRAST_ORDER, gene = SHEAR_PROG) |>
  left_join(mast_hits |> select(celltype, contrast, gene, avg_log2FC, p_val_adj, pct.1, pct.2, sig, reported, single_animal_driven,
                                low_detection, background_suspect, n_animals_detected_i1, n_animals_detected_i2, max_animal_share, max_share_animal),
            by = c("celltype", "contrast", "gene")) |>
  left_join(mast_universe |> select(celltype, contrast, gene, universe_effect = avg_log2FC, universe_p_adj = p_val_adj),
            by = c("celltype", "contrast", "gene")) |>
  mutate(status = case_when(!is.na(p_val_adj) ~ "hit table (>= 5% detection)",
                            !is.na(universe_p_adj) ~ "ranking universe only (1-5% detection)",
                            TRUE ~ "not tested (< 1% detection)"),
         effect_scale = ifelse(contrast == "Interaction", "interaction_log2", "avg_log2FC"),
         mark = paste0(ifelse(sig %in% TRUE, as.character(fdr_stars(p_val_adj)), ""),
                       flag_mark(sig %in% TRUE, single_animal_driven, low_detection, background_suspect)))
readr::write_csv(module_mast |> dplyr::rename(det_i1_prop = pct.1, det_i2_prop = pct.2, max_animal_share_prop = max_animal_share),
                 file.path(path_out, "ec_klf_mast_module_genes.csv"))
klf_mast <- module_mast |> filter(gene %in% KLF_GENES)
klf_mast |> arrange(gene, match(celltype, SUBTYPE_ORDER), match(contrast, CONTRAST_ORDER)) |>
  transmute(gene, subtype = short_subtype(celltype), contrast = CONTRAST_LABEL[contrast], effect_scale, effect = round(avg_log2FC, 3),
            p_bonferroni = format_pvalue(p_val_adj), mark = md_escape(mark), reported, single_animal_driven, low_detection,
            `detection i1 / i2 (proportion)` = paste(pct.1, pct.2, sep = " / ")) |>
  knitr::kable(row.names = FALSE, caption = "Klf2 and Klf4 in the MAST screen of notebook 05. Effect: avg_log2FC (pairwise) or interaction_log2 (Age x HFpEF). Mark: */**/*** Bonferroni p < 0.05/0.01/0.001 with |effect| > 0.2; \u2020 single-animal-driven, \u00b0 low detection, \u2021 background suspect.")
klf_mast_summary <- klf_mast |> group_by(gene, contrast) |>
  summarise(n_reported = sum(reported), n_lower = sum(reported & avg_log2FC < 0), n_higher = sum(reported & avg_log2FC > 0),
            subtypes = paste(short_subtype(celltype[reported]), collapse = ", "), n_flagged = sum(reported & (single_animal_driven | low_detection)),
            .groups = "drop")
mast_count <- function(g, cn) with(filter(klf_mast_summary, gene == g, contrast == cn),
                                   ifelse(n_reported == 0, "no subtype",
                                          sprintf("%d of %d subtypes (%d lower, %d higher: %s)", n_reported, length(SUBTYPE_ORDER), n_lower, n_higher, subtypes)))
mast_focus <- function(g, cn) with(filter(klf_mast, gene == g, contrast == cn, celltype == FOCUS_CT),
                                   sprintf("%+.2f, Bonferroni p = %s", avg_log2FC, format_pvalue(p_val_adj)))

profile_values <- as.data.frame(t(readout_means)) |> mutate(profile = colnames(readout_means)) |>
  left_join(profile_meta, by = "profile")
module_res <- bind_rows(lapply(LEVEL_ORDER, function(level_name) {
  profiles <- filter(profile_meta, level == level_name, in_model)
  fit_animal_model(readout_means["module", profiles$profile, drop = FALSE], profiles) |> mutate(level = level_name)
})) |>
  group_by(contrast) |> mutate(FDR = p.adjust(p_value, "BH")) |> ungroup() |>
  mutate(level = factor(level, levels = LEVEL_ORDER), contrast = factor(contrast, levels = CONTRAST_ORDER))
readr::write_csv(module_res, file.path(path_out, "ec_klf_module_animal_model.csv"))
module_group_means <- profile_values |> filter(in_model) |> group_by(level, group) |> summarise(mean = mean(module), .groups = "drop") |>
  pivot_wider(names_from = group, values_from = mean) |> select(level, all_of(GROUP_ORDER))
module_res |> arrange(contrast, level) |>
  transmute(contrast = CONTRAST_LABEL[as.character(contrast)], level, effect = round(effect, 4), CI = sprintf("%.4f to %.4f", CI_low, CI_high),
            p = signif(p_value, 2), FDR = signif(FDR, 2), n_animals) |>
  knitr::kable(row.names = FALSE, caption = "Module score, animal-level model. Effect: difference in mean module score (pairwise) or difference of differences (Age x HFpEF), 95% CI; BH within contrast over the six levels.")
knitr::kable(module_group_means |> mutate(across(all_of(GROUP_ORDER), \(x) round(x, 4))), caption = "Group means of the per-animal module score.")
module_best <- module_res |> slice_min(FDR, n = 1, with_ties = FALSE)
module_line <- function(cn, lv) with(filter(module_res, contrast == cn, level == lv), sprintf("%+.3f (95%% CI %.3f to %.3f, FDR = %s)", effect, CI_low, CI_high, fmt_fdr(FDR)))

heatmap_df <- module_mast |>
  mutate(gene = factor(gene, levels = SHEAR_PROG), celltype = factor(celltype, levels = rev(SUBTYPE_ORDER)),
         contrast = factor(contrast, levels = CONTRAST_ORDER), effect_capped = pmax(pmin(avg_log2FC, 1.2), -1.2),
         label = ifelse(is.na(avg_log2FC), "\u2013", mark))
fig <- ggplot(heatmap_df, aes(contrast, celltype, fill = effect_capped)) +
  geom_tile(colour = "white", linewidth = .6) +
  geom_text(aes(label = label, colour = is.na(avg_log2FC)), size = 3, show.legend = FALSE) +
  facet_wrap(~ gene, nrow = 2) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, limits = c(-1.2, 1.2), na.value = "grey90",
                       breaks = c(-1.2, -.6, 0, .6, 1.2), name = "effect\navg_log2FC |\ninteraction_log2") +
  scale_colour_manual(values = c(`FALSE` = "grey15", `TRUE` = "grey55")) +
  scale_x_discrete(labels = CONTRAST_CODE) + scale_y_discrete(labels = short_subtype) +
  labs(x = NULL, y = NULL, title = "Shear-protective module genes in the MAST screen (notebook 05)",
       subtitle = paste0("*/**/*** Bonferroni p < 0.05/0.01/0.001 with |effect| > ", LFC_CUT, "; \u2020 single-animal-driven, \u00b0 low detection, \u2021 background suspect; ",
                         "grey dash = below the 5% detection floor, not in the hit table.\n", paste(CONTRAST_CODE, CONTRAST_LABEL, sep = " = ", collapse = "; "),
                         ". Pairwise fill = avg_log2FC; Age x HFpEF fill = interaction_log2 (positive = larger HFpEF effect in old).")) +
  theme_fig + theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_text(size = 7),
                    strip.background = element_blank(), strip.text = element_text(face = "bold.italic"))
save_fig(fig, "ec_klf_module_genes_mast_heatmap", width = 13, height = 5.2)
fig
module_gene_counts <- module_mast |> group_by(contrast) |>
  summarise(lower = sum(reported %in% TRUE & avg_log2FC < 0), higher = sum(reported %in% TRUE & avg_log2FC > 0), .groups = "drop")
module_gene_text <- function(cn) with(filter(module_gene_counts, contrast == cn), sprintf("%d lower / %d higher", lower, higher))

animal_point_layers <- function() list(
  stat_summary(aes(group = group), fun = mean, fun.min = mean, fun.max = mean, geom = "errorbar", width = .45, linewidth = .35, colour = "grey45"),
  geom_point(aes(fill = group, colour = group, shape = animal_class), size = 2.3, stroke = .8,
             position = position_jitter(width = .12, height = 0, seed = SEED)),
  scale_fill_manual(values = GROUP_FILL, labels = GROUP_LABEL, name = NULL, drop = FALSE),
  scale_colour_manual(values = GROUP_FILL, labels = GROUP_LABEL, name = NULL, drop = FALSE),
  scale_shape_manual(values = ANIMAL_CLASS_SHAPE, labels = ANIMAL_CLASS_LABEL, name = NULL, drop = FALSE),
  scale_x_discrete(limits = GROUP_ORDER),
  guides(fill = guide_legend(override.aes = list(shape = 21)), shape = guide_legend(override.aes = list(fill = "grey60", colour = "grey30"))),
  theme_fig, theme(legend.position = "bottom", legend.box = "vertical", strip.background = element_blank()))
readout_label <- c(Klf2 = "Klf2 (mean log-normalised expression)", Klf4 = "Klf4 (mean log-normalised expression)", module = "shear-protective module score")
module_code <- function(lv) with(filter(module_res, level == lv, FDR < FDR_CUT), paste(sprintf("%s FDR %s", CONTRAST_CODE[as.character(contrast)], fmt_fdr(FDR)), collapse = "; "))
mast_code <- function(g) with(filter(klf_mast, gene == g, celltype == FOCUS_CT, reported), paste0(CONTRAST_CODE[contrast], fdr_stars(p_val_adj), collapse = " "))
expression_strips <- expand_grid(level = c(FOCUS_CT, ALL_EC), readout = names(readout_label)) |>
  mutate(note_text = case_when(readout == "module" ~ paste("limma:", vapply(level, module_code, "")),
                               level == FOCUS_CT ~ paste("MAST:", vapply(readout, mast_code, "")),
                               TRUE ~ "MAST: tested per subtype"),
         note_text = sub(": $", ": no FDR < 0.05", note_text),
         panel = paste0(level, "\n", readout_label[readout], "\n", note_text))
expression_plot_df <- profile_values |> filter(level %in% c(FOCUS_CT, ALL_EC), in_model) |>
  pivot_longer(c(Klf2, Klf4, module), names_to = "readout", values_to = "value") |>
  mutate(level = as.character(level), group = factor(as.character(group), levels = GROUP_ORDER)) |>
  left_join(expression_strips, by = c("level", "readout")) |>
  mutate(panel = factor(panel, levels = expression_strips$panel))
fig <- ggplot(expression_plot_df, aes(group, value)) + animal_point_layers() +
  facet_wrap(~ panel, scales = "free_y", nrow = 2) +
  labs(x = NULL, y = NULL, title = "Klf2, Klf4 and the shear-protective module per animal",
       subtitle = paste0("Points = animals, bars = group means. Strip: contrasts passing the MAST screen (Klf2/Klf4, ", FOCUS_CT,
                         ") or animal-level FDR < 0.05 (module).\n", paste(CONTRAST_CODE, CONTRAST_LABEL, sep = " = ", collapse = "; "))) +
  theme(strip.text = element_text(size = 7.5, face = "plain"))
save_fig(fig, "ec_klf_expression_per_animal", width = 11, height = 6.2)
fig
expression_group_means <- expression_plot_df |> group_by(level, readout, group) |> summarise(mean = mean(value), .groups = "drop") |>
  pivot_wider(names_from = group, values_from = mean)
flagged_rank <- expression_plot_df |> filter(level == FOCUS_CT) |> group_by(readout, group) |>
  mutate(rank_in_group = rank(value, ties.method = "first"), n_group = n()) |> ungroup() |> filter(animal_class != "pass")
flagged_text <- function(r) with(filter(flagged_rank, readout == r), paste(sprintf("%s %d/%d", sample, rank_in_group, n_group), collapse = ", "))
old_young_text <- function(lv, r) with(filter(expression_group_means, level == lv, readout == r),
                                       sprintf("OC %.3f vs YC %.3f, OH %.3f vs YH %.3f", OC, YC, OH, YH))

regulons <- readr::read_csv(path_net_csv, show_col_types = FALSE)
klf_targets <- lapply(setNames(KLF_GENES, KLF_GENES), function(tf) regulons$target[regulons$source == tf])
data.frame(TFs = n_distinct(regulons$source), interactions = nrow(regulons), repressive = sum(regulons$weight < 0),
           Klf2_targets = length(klf_targets$Klf2), Klf4_targets = length(klf_targets$Klf4),
           Klf2_targets_in_EC = sum(klf_targets$Klf2 %in% rownames(counts)), Klf4_targets_in_EC = sum(klf_targets$Klf4 %in% rownames(counts))) |>
  knitr::kable(caption = "CollecTRI regulon network (mouse symbols); *_in_EC = targets present in the endothelial assay.")

logcpm_level <- function(level_name) {
  profiles <- filter(profile_meta, level == level_name, in_model)
  dge <- DGEList(as.matrix(pb_counts[, profiles$profile]))
  dge <- calcNormFactors(dge[filterByExpr(dge, group = profiles$group), , keep.lib.sizes = FALSE])
  cpm(dge, log = TRUE)
}
logcpm <- lapply(setNames(LEVEL_ORDER, LEVEL_ORDER), logcpm_level)
tf_scores <- lapply(logcpm, function(expr_mat)
  run_ulm(mat = expr_mat, network = regulons, .source = source, .target = target, .mor = weight, minsize = MIN_TARGETS) |>
    select(source, condition, score) |> pivot_wider(names_from = condition, values_from = score) |>
    tibble::column_to_rownames("source") |> as.matrix())
tf_screen <- bind_rows(lapply(LEVEL_ORDER, function(level_name) {
  scores <- tf_scores[[level_name]][, filter(profile_meta, level == level_name, in_model)$profile]
  fit_animal_model(scores, filter(profile_meta, level == level_name, in_model), moderated = TRUE) |> mutate(level = level_name)
})) |> dplyr::rename(TF = feature) |>
  group_by(level, contrast) |>
  mutate(FDR = p.adjust(p_value, "BH"), n_TF = n(), rank_p = rank(p_value, ties.method = "min"),
         rank_in_direction = ifelse(t < 0, rank(t, ties.method = "min"), rank(-t, ties.method = "min"))) |> ungroup() |>
  mutate(level = factor(level, levels = LEVEL_ORDER), contrast = factor(contrast, levels = CONTRAST_ORDER))
readr::write_csv(tf_screen, file.path(path_out, "ec_klf_tf_activity_screen.csv"))
universe_tab <- data.frame(level = LEVEL_ORDER, profiles = vapply(logcpm, ncol, 1L), genes = vapply(logcpm, nrow, 1L),
                           TFs_scored = vapply(tf_scores, nrow, 1L),
                           Klf2_targets = vapply(logcpm, function(m) sum(klf_targets$Klf2 %in% rownames(m)), 1L),
                           Klf4_targets = vapply(logcpm, function(m) sum(klf_targets$Klf4 %in% rownames(m)), 1L))
knitr::kable(universe_tab, row.names = FALSE, caption = "Expression universe and scored TFs per level.")
screen_summary <- tf_screen |> filter(level %in% c(FOCUS_CT, ALL_EC)) |> group_by(level, contrast) |>
  summarise(TFs = n(), FDR05 = sum(FDR < FDR_CUT), lower_FDR05 = sum(FDR < FDR_CUT & effect < 0), higher_FDR05 = sum(FDR < FDR_CUT & effect > 0),
            min_FDR = min(FDR), top_TFs = paste(TF[order(p_value)][1:5], collapse = ", "), .groups = "drop")
screen_summary |> mutate(contrast = CONTRAST_LABEL[as.character(contrast)], min_FDR = signif(min_FDR, 2)) |>
  knitr::kable(caption = "TF screen per contrast: TFs at FDR < 0.05 (BH within contrast over all scored TFs) and the five TFs with the smallest p.")

klf_activity <- tf_screen |> filter(TF %in% KLF_GENES) |>
  group_by(contrast) |> mutate(FDR_klf = p.adjust(p_value, "BH")) |> ungroup()
readr::write_csv(klf_activity, file.path(path_out, "ec_klf_activity_pseudobulk.csv"))
klf_activity |> arrange(TF, level, contrast) |>
  transmute(TF = toupper(TF), level, contrast = CONTRAST_LABEL[as.character(contrast)], effect = round(effect, 3), CI = sprintf("%.2f to %.2f", CI_low, CI_high),
            p = signif(p_value, 2), FDR_targeted = signif(FDR_klf, 2), FDR_screen = signif(FDR, 2),
            rank_p = paste0(rank_p, "/", n_TF), rank_in_direction = paste0(rank_in_direction, "/", n_TF, ifelse(effect < 0, " (lower)", " (higher)")), n_animals) |>
  knitr::kable(row.names = FALSE, caption = "KLF2/KLF4 pseudobulk activity. Effect: difference in ULM score (pairwise) or difference of differences (Age x HFpEF). FDR_targeted: BH within contrast over KLF2/KLF4 x six levels; FDR_screen: BH within contrast over all TFs of the level; rank_p: rank by p among all TFs of the level (1 = smallest); rank_in_direction: rank among TFs moving the same way (1 = strongest).")
klf_targeted_hits <- filter(klf_activity, FDR_klf < FDR_CUT)
klf_rank_line <- function(tf, cn, lv) with(filter(klf_activity, TF == tf, contrast == cn, level == lv),
  sprintf("%+.2f (p = %s, FDR = %s; rank %d/%d by p)", effect, fmt_fdr(p_value), fmt_fdr(FDR_klf), rank_p, n_TF))

klf_regulons <- filter(regulons, source %in% KLF_GENES)
nucleus_universe <- rownames(logcpm[[ALL_EC]])
nucleus_chunks <- split(seq_len(ncol(lognorm)), ceiling(seq_len(ncol(lognorm)) / ULM_CHUNK))
klf_nucleus_scores <- do.call(cbind, lapply(nucleus_chunks, function(idx)
  run_ulm(mat = as.matrix(lognorm[nucleus_universe, idx]), network = klf_regulons, .source = source, .target = target, .mor = weight,
          minsize = MIN_TARGETS) |>
    select(source, condition, score) |> pivot_wider(names_from = condition, values_from = score) |>
    tibble::column_to_rownames("source") |> as.matrix() |> (\(m) m[KLF_GENES, colnames(lognorm)[idx]])()))
nucleus_activity_means <- profile_means(klf_nucleus_scores)
klf_nucleus_res <- bind_rows(lapply(LEVEL_ORDER, function(level_name) {
  profiles <- filter(profile_meta, level == level_name, in_model)
  fit_animal_model(nucleus_activity_means[, profiles$profile], profiles) |> mutate(level = level_name)
})) |> dplyr::rename(TF = feature) |>
  group_by(contrast) |> mutate(FDR = p.adjust(p_value, "BH")) |> ungroup() |>
  mutate(level = factor(level, levels = LEVEL_ORDER), contrast = factor(contrast, levels = CONTRAST_ORDER))
readr::write_csv(klf_nucleus_res, file.path(path_out, "ec_klf_activity_nucleus_mean.csv"))
profile_readouts <- profile_values |>
  left_join(data.frame(profile = colnames(nucleus_activity_means), Klf2_activity_nucleus = nucleus_activity_means["Klf2", ],
                       Klf4_activity_nucleus = nucleus_activity_means["Klf4", ]), by = "profile") |>
  left_join(bind_rows(lapply(tf_scores, function(m) data.frame(profile = colnames(m), Klf2_activity_pseudobulk = m["Klf2", ],
                                                                Klf4_activity_pseudobulk = m["Klf4", ]))), by = "profile")
readr::write_csv(profile_readouts |> select(profile, sample, level, group, sex, seq, qc_flag, note, n_nuclei, in_model, Klf2_mean_expr = Klf2,
                                          Klf4_mean_expr = Klf4, module_score = module, ends_with("pseudobulk"), ends_with("nucleus")),
                 file.path(path_out, "ec_klf_profile_readouts.csv"))
profile_readouts |> filter(in_model, level %in% c(FOCUS_CT, ALL_EC)) |>
  group_by(level, group) |>
  summarise(animals = n(), KLF2_nucleus_mean = mean(Klf2_activity_nucleus), KLF4_nucleus_mean = mean(Klf4_activity_nucleus),
            KLF2_pseudobulk = mean(Klf2_activity_pseudobulk), KLF4_pseudobulk = mean(Klf4_activity_pseudobulk), .groups = "drop") |>
  mutate(group = factor(as.character(group), levels = GROUP_ORDER), across(where(is.double), \(x) round(x, 3))) |> arrange(level, group) |>
  knitr::kable(caption = "Group means of the per-animal KLF2/KLF4 activity: mean of the per-nucleus ULM scores and pseudobulk ULM score.")
klf_nucleus_res |> transmute(TF = toupper(TF), level, contrast = CONTRAST_LABEL[as.character(contrast)], effect = round(effect, 3),
                             CI = sprintf("%.3f to %.3f", CI_low, CI_high), p = signif(p_value, 2), FDR = signif(FDR, 2), n_animals) |>
  arrange(TF, level) |>
  knitr::kable(row.names = FALSE, caption = "Per-nucleus KLF2/KLF4 activity averaged per animal, animal-level model; BH within contrast over KLF2/KLF4 x six levels.")
nucleus_line <- function(tf, cn, lv) with(filter(klf_nucleus_res, TF == tf, contrast == cn, level == lv),
  sprintf("%+.3f (95%% CI %.3f to %.3f, FDR = %s)", effect, CI_low, CI_high, fmt_fdr(FDR)))
nucleus_best <- klf_nucleus_res |> slice_min(FDR, n = 1, with_ties = FALSE)

volcano_df <- tf_screen |> filter(level %in% c(FOCUS_CT, ALL_EC)) |>
  mutate(class = case_when(TF %in% KLF_GENES ~ "KLF2 / KLF4", FDR < FDR_CUT ~ sprintf("FDR < %s", FDR_CUT), TRUE ~ "other TFs"),
         class = factor(class, levels = c("other TFs", sprintf("FDR < %s", FDR_CUT), "KLF2 / KLF4")),
         label = ifelse(TF %in% KLF_GENES, sprintf("%s (%d)", toupper(TF), rank_p), NA))
volcano_top <- volcano_df |> filter(!TF %in% KLF_GENES) |> group_by(level, contrast) |> slice_min(p_value, n = 3) |> ungroup()
fig <- ggplot(volcano_df, aes(effect, -log10(p_value))) +
  geom_vline(xintercept = 0, colour = "grey75", linewidth = .3) +
  geom_point(aes(colour = class, size = class)) +
  geom_text_repel(data = volcano_top, aes(label = toupper(TF)), size = 2.3, colour = "grey45", min.segment.length = 0, segment.size = .2, seed = SEED) +
  geom_text_repel(aes(label = label), size = 2.8, colour = "#D55E00", fontface = "bold", min.segment.length = 0, segment.size = .25,
                  na.rm = TRUE, seed = SEED, box.padding = .5) +
  facet_grid(level ~ contrast, scales = "free", labeller = labeller(contrast = CONTRAST_LABEL)) +
  scale_colour_manual(values = c("other TFs" = "grey80", setNames("#0072B2", sprintf("FDR < %s", FDR_CUT)), "KLF2 / KLF4" = "#D55E00"), name = NULL, drop = FALSE) +
  scale_size_manual(values = c(.7, 1.2, 2.4), guide = "none") +
  labs(x = "TF activity effect (difference in pseudobulk ULM score; Age x HFpEF: difference of differences)", y = expression(-log[10] ~ p ~ "(moderated t)"),
       title = "Unbiased TF-activity screen (CollecTRI regulons, ULM on per-animal pseudobulk)",
       subtitle = sprintf("Points = TFs; KLF2/KLF4 in orange with their rank by p in brackets; grey labels = three smallest p per panel; BH within contrast over %s TFs (%s) and %s (%s).",
                          universe_tab$TFs_scored[universe_tab$level == FOCUS_CT], FOCUS_CT, universe_tab$TFs_scored[universe_tab$level == ALL_EC], ALL_EC)) +
  theme_fig + theme(legend.position = "bottom", strip.background = element_blank())
save_fig(fig, "ec_klf_tf_activity_volcano", width = 14, height = 6.4)
fig

VIEW_LABEL <- c(pseudobulk = "pseudobulk ULM score", nucleus = "mean per-nucleus ULM score")
activity_strips <- bind_rows(klf_activity |> transmute(TF, level, contrast, FDR = FDR_klf, view = "pseudobulk"),
                             klf_nucleus_res |> transmute(TF, level, contrast, FDR, view = "nucleus")) |>
  filter(level %in% c(FOCUS_CT, ALL_EC)) |> mutate(level = as.character(level)) |>
  group_by(view, TF, level) |>
  summarise(fdr_text = paste(sprintf("%s FDR %s", CONTRAST_CODE[as.character(contrast[FDR < FDR_CUT])], fmt_fdr(FDR[FDR < FDR_CUT])), collapse = "; "),
            .groups = "drop") |>
  mutate(fdr_text = ifelse(nzchar(fdr_text), fdr_text, "no FDR < 0.05"),
         panel = paste0(toupper(TF), ", ", level, "\n", VIEW_LABEL[view], "\n", fdr_text)) |>
  arrange(match(view, names(VIEW_LABEL)), TF, match(level, c(FOCUS_CT, ALL_EC)))
activity_plot_df <- profile_readouts |> filter(level %in% c(FOCUS_CT, ALL_EC), in_model) |>
  pivot_longer(c(Klf2_activity_pseudobulk, Klf4_activity_pseudobulk, Klf2_activity_nucleus, Klf4_activity_nucleus), names_to = "readout", values_to = "value") |>
  mutate(TF = sub("_.*", "", readout), view = sub(".*_", "", readout), level = as.character(level),
         group = factor(as.character(group), levels = GROUP_ORDER)) |>
  left_join(activity_strips, by = c("view", "TF", "level")) |>
  mutate(panel = factor(panel, levels = activity_strips$panel))
fig <- ggplot(activity_plot_df, aes(group, value)) + animal_point_layers() +
  facet_wrap(~ panel, scales = "free_y", nrow = 2) +
  labs(x = NULL, y = "regulon activity (ULM score)", title = "KLF2 and KLF4 regulon activity per animal",
       subtitle = paste0("Points = animals, bars = group means. Top: pseudobulk ULM score; bottom: mean of the per-nucleus ULM scores. Strip: contrasts at FDR < 0.05 (BH over KLF2/KLF4 x six levels).\n",
                         paste(CONTRAST_CODE, CONTRAST_LABEL, sep = " = ", collapse = "; "))) +
  theme(strip.text = element_text(size = 7.5, face = "plain"))
save_fig(fig, "ec_klf_activity_per_animal", width = 12, height = 6.4)
fig

activity_group_means <- activity_plot_df |> group_by(panel, group) |> summarise(mean = mean(value), .groups = "drop") |>
  pivot_wider(names_from = group, values_from = mean)
activity_flagged <- activity_plot_df |> filter(level == FOCUS_CT, view == "pseudobulk") |> group_by(TF, group) |>
  mutate(rank_in_group = rank(value, ties.method = "first"), n_group = n()) |> ungroup() |> filter(animal_class != "pass")
activity_flagged_text <- function(tf) with(filter(activity_flagged, TF == tf), paste(sprintf("%s %d/%d", sample, rank_in_group, n_group), collapse = ", "))
activity_overlap <- activity_plot_df |> group_by(panel) |>
  summarise(overlap_control = max(value[group == "OC"]) > min(value[group == "YC"]), overlap_hfpef = max(value[group == "OH"]) > min(value[group == "YH"]), .groups = "drop")

fmt_limma <- function(effect, fdr, digits) sprintf(paste0("%+.", digits, "f (FDR %s)%s"), effect, fmt_fdr(fdr), ifelse(fdr < FDR_CUT, " *", ""))
READOUT_ORDER <- c("Klf2 mRNA (MAST)", "Klf4 mRNA (MAST)", "module score", "KLF2 activity, pseudobulk", "KLF4 activity, pseudobulk",
                   "KLF2 activity, per-nucleus mean", "KLF4 activity, per-nucleus mean")
mast_answer <- klf_mast |> filter(celltype == FOCUS_CT) |> left_join(klf_mast_summary, by = c("gene", "contrast")) |>
  transmute(contrast, readout = paste0(gene, " mRNA (MAST)"),
            capillary = sprintf("%+.2f (p_Bonf %s)%s", avg_log2FC, format_pvalue(p_val_adj), ifelse(nzchar(mark), paste0(" ", mark), "")),
            all_ec = sprintf("reported in %d/%d subtypes (%d lower, %d higher)", n_reported, length(SUBTYPE_ORDER), n_lower, n_higher),
            significant_capillary = reported, significant_all_ec = n_reported > 0)
limma_answer <- bind_rows(
  module_res |> transmute(contrast, level, readout = "module score", text = fmt_limma(effect, FDR, 4), sig = FDR < FDR_CUT),
  klf_activity |> transmute(contrast, level, readout = paste0(toupper(TF), " activity, pseudobulk"),
                            text = sprintf("%s; rank %d/%d", fmt_limma(effect, FDR_klf, 2), rank_p, n_TF), sig = FDR_klf < FDR_CUT),
  klf_nucleus_res |> transmute(contrast, level, readout = paste0(toupper(TF), " activity, per-nucleus mean"), text = fmt_limma(effect, FDR, 3), sig = FDR < FDR_CUT)) |>
  filter(level %in% c(FOCUS_CT, ALL_EC)) |>
  mutate(contrast = as.character(contrast), level = ifelse(level == FOCUS_CT, "capillary", "all_ec")) |>
  pivot_wider(names_from = level, values_from = c(text, sig)) |>
  transmute(contrast, readout, capillary = text_capillary, all_ec = text_all_ec, significant_capillary = sig_capillary, significant_all_ec = sig_all_ec)
answer <- bind_rows(mast_answer, limma_answer) |>
  mutate(contrast = factor(contrast, levels = CONTRAST_ORDER), readout = factor(readout, levels = READOUT_ORDER)) |>
  arrange(contrast, readout)
readr::write_csv(answer, file.path(path_out, "ec_klf_results_per_contrast.csv"))
answer |> transmute(contrast = CONTRAST_LABEL[as.character(contrast)], readout, `Capillary EC` = md_escape(capillary), `All EC` = md_escape(all_ec)) |>
  knitr::kable(row.names = FALSE, caption = "Per contrast: effect and multiplicity-adjusted p. MAST (notebook 05): avg_log2FC (pairwise) or interaction_log2 (Age x HFpEF), Bonferroni p, */**/*** passing the screen, flags \u2020/\u00b0/\u2021 as in Step 2.1. Module and activity (animal-level limma): difference in per-animal mean (pairwise) or difference of differences (Age x HFpEF), BH FDR, * FDR < 0.05; activity rank = rank by p among all TFs of the level.")
describe_contrast <- function(cn) {
  rows <- filter(answer, contrast == cn)
  cell <- function(r, col) md_escape(rows[[col]][rows$readout == r])
  changed <- as.character(rows$readout[rows$significant_capillary | rows$significant_all_ec])
  sprintf(paste("**%s.** *Klf2* mRNA %s; %s %s. *Klf4* mRNA %s; %s %s. Module score: %s %s, all EC %s.",
                "KLF2 activity (pseudobulk): %s %s, all EC %s; KLF4: %s %s, all EC %s.",
                "Per-nucleus mean activity, %s: KLF2 %s, KLF4 %s. Readouts passing their threshold (MAST: reported in any subtype; animal-level: FDR < 0.05 in %s or all EC): %s."),
          CONTRAST_LABEL[[cn]], cell("Klf2 mRNA (MAST)", "all_ec"), FOCUS_CT, cell("Klf2 mRNA (MAST)", "capillary"),
          cell("Klf4 mRNA (MAST)", "all_ec"), FOCUS_CT, cell("Klf4 mRNA (MAST)", "capillary"),
          FOCUS_CT, cell("module score", "capillary"), cell("module score", "all_ec"),
          FOCUS_CT, cell("KLF2 activity, pseudobulk", "capillary"), cell("KLF2 activity, pseudobulk", "all_ec"),
          FOCUS_CT, cell("KLF4 activity, pseudobulk", "capillary"), cell("KLF4 activity, pseudobulk", "all_ec"),
          FOCUS_CT, cell("KLF2 activity, per-nucleus mean", "capillary"), cell("KLF4 activity, per-nucleus mean", "capillary"), FOCUS_CT,
          ifelse(length(changed), paste(changed, collapse = ", "), "none"))
}

session_footer()
