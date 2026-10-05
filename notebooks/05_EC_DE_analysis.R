# 05_EC_DE_analysis.R
# Code-only export of the analysis pipeline. Data are not included; see README.md.

params <- list(rds_in = "results/03_ec_annotation/Endothelial_subcluster_annotated_clean.rds", 
    atlas_rds = "data/objects/integrated_LV_snRNAseq_seurat_clean.rds", out_dir = "results/05_ec_de")


knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = TRUE,
                      fig.align = "center", dpi = 150)
options(width = 110)

suppressPackageStartupMessages({
  library(Seurat); library(MAST); library(Matrix); library(tidyr); library(patchwork); library(ggrepel)
})
source("_common.R")
knitr::opts_chunk$set(cache = TRUE, autodep = TRUE, cache.lazy = FALSE, cache.extra = input_md5(params))
path_rds_in  <- resolve_path(params$rds_in)
path_atlas_rds <- resolve_path(params$atlas_rds)
path_out <- resolve_path(params$out_dir)
path_figures <- file.path(path_out, "figures"); dir.create(path_figures, recursive = TRUE, showWarnings = FALSE)

PADJ_CUT     <- 0.05
LFC_CUT      <- 0.2
MIN_CELLS    <- 30
MIN_ANIMALS  <- 3
MAST_COVARIATES <- c(DESIGN_COVARIATES, "cngeneson", "exon_prop")
GSEA_MIN_PCT <- 0.01
DETECT_FRAC  <- 0.10
MAX_ANIMAL_SHARE <- 0.5

VOLCANO_CT   <- "Capillary EC"
LEAD_N       <- 20
BKG_CUT      <- 0.4
FILL_SHAPES  <- c(21, 22, 24, 23, 25)
FLAG_COL     <- c(none = "grey25", `single-animal driven` = "#E69F00", `low detection` = "#CC79A7",
                  `background suspect` = "#009E73", `several flags` = "#D55E00")
LFC_SHOW     <- 3

print_versions(c("ggplot2", "dplyr", "Seurat", "MAST", "Matrix", "tidyr", "patchwork", "ggrepel"))

ec_obj <- readRDS(path_rds_in)
DefaultAssay(ec_obj) <- "RNA"
required <- c("EC_subtype", "sample", "individual", "age", "condition", "disease", DESIGN_COVARIATES, "batch", "nFeature_RNA", "exon_prop")
stopifnot(all(required %in% colnames(ec_obj@meta.data)))
SUBTYPE_ORDER <- levels(droplevels(factor(ec_obj$EC_subtype)))
PAIRWISE      <- CONTRASTS
CONTRASTS     <- c(PAIRWISE, list(list(name = "Interaction", i1 = NA, i2 = NA, meaning = "age x HFpEF interaction, (OH - OC) - (YH - YC)")))
ec_obj$EC_subtype <- factor(ec_obj$EC_subtype, levels = SUBTYPE_ORDER)
ec_obj$qc_flag    <- unname(coalesce(QC_FLAGGED_ANIMALS[as.character(ec_obj$sample)], "pass"))
ec_obj$note       <- unname(coalesce(ANIMAL_NOTES[as.character(ec_obj$sample)], ""))
animal_meta <- ec_obj@meta.data |> distinct(individual, sample, age, condition, disease, sex, batch, seq, qc_flag, note) |>
  mutate(across(where(is.factor), as.character))
if (anyDuplicated(animal_meta$individual)) stop("Expected one sample and one metadata row per animal")
knitr::kable(animal_meta |> arrange(match(disease, GROUP_ORDER), individual), row.names = FALSE,
             caption = "Animal-level metadata. qc_flag and note are annotations; every animal enters every fit.")

design <- ec_obj@meta.data |>
  group_by(age, condition) |>
  summarise(nuclei = n(), animals = n_distinct(individual),
            females = n_distinct(individual[sex == "F"]), males = n_distinct(individual[sex == "M"]),
            median_genes = median(nFeature_RNA), median_exon_prop = round(median(exon_prop), 3),
            flagged = paste(unique(individual[qc_flag != "pass"]), collapse = ", "), .groups = "drop") |>
  mutate(nuclei_per_animal = round(nuclei / animals))
knitr::kable(dplyr::rename(design, `median_exon_prop (proportion)` = median_exon_prop), caption = "Nuclei and animal representation by experimental group.")
knitr::kable(as.data.frame.matrix(table(ec_obj$individual, ec_obj$EC_subtype)),
             caption = "Nuclei per animal and subtype: single-animal dominance of a group is visible here.")
sex_tab   <- as.data.frame.matrix(table(factor(animal_meta$disease, levels = GROUP_ORDER), animal_meta$sex))
unbalanced_groups     <- rownames(sex_tab)[sex_tab$M != sex_tab$F]

knitr::kable(bind_rows(CONTRASTS) |> dplyr::rename(ident.1 = i1, ident.2 = i2),
             caption = "Contrasts. Pairwise avg_log2FC is positive when expression is higher in ident.1; Interaction interaction_log2 is positive when the HFpEF effect is larger in old than in young animals.")

fit_one <- function(sub, contrast_def) {
  contrast_obj  <- subset(sub, cells = colnames(sub)[sub$disease %in% c(contrast_def$i1, contrast_def$i2)])
  n_nuclei_i1 <- sum(contrast_obj$disease == contrast_def$i1); n_nuclei_i2 <- sum(contrast_obj$disease == contrast_def$i2)
  n_animals_i1 <- n_distinct(contrast_obj$individual[contrast_obj$disease == contrast_def$i1]); n_animals_i2 <- n_distinct(contrast_obj$individual[contrast_obj$disease == contrast_def$i2])
  if (n_nuclei_i1 < MIN_CELLS || n_nuclei_i2 < MIN_CELLS) return(list(skip = sprintf("too few nuclei (%d vs %d)", n_nuclei_i1, n_nuclei_i2)))
  if (n_animals_i1 < MIN_ANIMALS || n_animals_i2 < MIN_ANIMALS) return(list(skip = sprintf("too few represented animals (%d vs %d)", n_animals_i1, n_animals_i2)))
  Idents(contrast_obj) <- "disease"
  contrast_obj$disease <- droplevels(factor(contrast_obj$disease))
  contrast_obj$sex <- droplevels(factor(contrast_obj$sex)); contrast_obj$seq <- droplevels(factor(contrast_obj$seq))
  contrast_obj$cngeneson <- as.numeric(scale(contrast_obj$nFeature_RNA))
  latent_vars <- MAST_COVARIATES
  design_matrix <- model.matrix(reformulate(c("disease", latent_vars)), data = contrast_obj@meta.data)
  if (qr(design_matrix)$rank < ncol(design_matrix)) return(list(skip = "group and covariates are not jointly estimable"))

  run_mast <- function(latent, features = NULL, min_pct = 0.05) {
    mast_res <- FindMarkers(contrast_obj, ident.1 = contrast_def$i1, ident.2 = contrast_def$i2, assay = "RNA", test.use = "MAST",
                     features = features, latent.vars = latent, min.pct = min_pct, logfc.threshold = 0,
                     only.pos = FALSE, verbose = FALSE)
    data.frame(gene = rownames(mast_res), mast_res, row.names = NULL)
  }
  mast_adjusted <- run_mast(latent_vars); mast_unadjusted <- run_mast(NULL); mast_universe <- run_mast(latent_vars, min_pct = GSEA_MIN_PCT)
  count_hits <- function(mast_res) sum(mast_res$p_val_adj < PADJ_CUT & abs(mast_res$avg_log2FC) > LFC_CUT, na.rm = TRUE)
  list(adj = mast_adjusted, raw = mast_unadjusted, uni = mast_universe,
       meta = data.frame(n_ident1 = n_nuclei_i1, n_ident2 = n_nuclei_i2, n_animals_ident1 = n_animals_i1, n_animals_ident2 = n_animals_i2,
                         covariates = paste(latent_vars, collapse = "+"), n_hits_adj = count_hits(mast_adjusted), n_hits_raw = count_hits(mast_unadjusted)))
}

fit_interaction <- function(sub) {
  nucleus_meta <- sub@meta.data |>
    mutate(age = factor(age, levels = c("young", "old")), condition = factor(condition, levels = c("control", "HFpEF")),
           sex = droplevels(factor(sex)), seq = droplevels(factor(seq)), cngeneson = as.numeric(scale(nFeature_RNA)))
  n_nuclei_group <- table(nucleus_meta$disease)[GROUP_ORDER]; n_animals_group <- tapply(nucleus_meta$individual, nucleus_meta$disease, n_distinct)[GROUP_ORDER]
  if (any(is.na(n_nuclei_group) | n_nuclei_group < MIN_CELLS)) return(list(skip = sprintf("too few nuclei in a group (%s)", paste(coalesce(n_nuclei_group, 0L), collapse = "/"))))
  if (any(is.na(n_animals_group) | n_animals_group < MIN_ANIMALS)) return(list(skip = sprintf("too few represented animals in a group (%s)", paste(coalesce(n_animals_group, 0L), collapse = "/"))))
  interaction_formula <- reformulate(c("age * condition", MAST_COVARIATES))
  design_matrix <- model.matrix(interaction_formula, data = nucleus_meta)
  if (qr(design_matrix)$rank < ncol(design_matrix)) return(list(skip = "interaction and covariates are not jointly estimable"))
  expr_data   <- GetAssayData(sub, layer = "data"); detection_frac <- rowMeans(expr_data > 0)
  is_hfpef  <- nucleus_meta$condition == "HFpEF"
  screened <- names(detection_frac)[detection_frac >= 0.05]

  at_mean  <- colMeans(design_matrix)
  group_design_row  <- function(is_old, is_hfpef) replace(at_mean, c("ageold", "conditionHFpEF", "ageold:conditionHFpEF"), c(is_old, is_hfpef, is_old * is_hfpef))
  fit <- function(genes) {
    sca_obj <- FromMatrix(as.matrix(expr_data[genes, , drop = FALSE]), cData = nucleus_meta[, c("age", "condition", MAST_COVARIATES)], check_sanity = FALSE)
    zlm_fit   <- zlm(interaction_formula, sca_obj, method = "bayesglm", ebayes = TRUE, silent = TRUE)
    lrt_p   <- lrTest(zlm_fit, "age:condition")[, "hurdle", "Pr(>Chisq)"]
    interaction_lfc  <- logFC(zlm_fit, rbind(YC = group_design_row(0, 0)), cbind(YH = group_design_row(0, 1), OC = group_design_row(1, 0), OH = group_design_row(1, 1)))$logFC
    data.frame(gene = genes, p_val = lrt_p[genes], avg_log2FC = ((interaction_lfc[genes, "OH"] - interaction_lfc[genes, "OC"]) - interaction_lfc[genes, "YH"]) / log(2),
               pct.1 = round(rowMeans(expr_data[genes, is_hfpef, drop = FALSE] > 0), 3), pct.2 = round(rowMeans(expr_data[genes, !is_hfpef, drop = FALSE] > 0), 3),
               p_val_adj = pmin(lrt_p[genes] * nrow(sub), 1), row.names = NULL)
  }
  mast_adjusted <- fit(screened); mast_universe <- fit(names(detection_frac)[detection_frac >= GSEA_MIN_PCT])
  list(adj = mast_adjusted, uni = mast_universe,
       meta = data.frame(n_ident1 = sum(is_hfpef), n_ident2 = sum(!is_hfpef), n_animals_ident1 = n_distinct(nucleus_meta$individual[is_hfpef]),
                         n_animals_ident2 = n_distinct(nucleus_meta$individual[!is_hfpef]), covariates = paste(MAST_COVARIATES, collapse = "+"),
                         n_hits_adj = sum(mast_adjusted$p_val_adj < PADJ_CUT & abs(mast_adjusted$avg_log2FC) > LFC_CUT, na.rm = TRUE), n_hits_raw = NA_integer_))
}

path_de_files <- file.path(path_out, c("ec_DE_all_MAST.csv", "ec_DE_all_MAST_unadj.csv", "ec_DE_all_summary.csv",
                                 "ec_DE_all_skipped.csv", "ec_DE_all_MAST_gsea_universe.csv"))
adjusted_list <- unadjusted_list <- universe_list <- summary_list <- skip_list <- list()
for (subtype_name in SUBTYPE_ORDER) {
  sub <- subset(ec_obj, cells = colnames(ec_obj)[ec_obj$EC_subtype == subtype_name])
  for (contrast_def in CONTRASTS) {
    key <- paste(subtype_name, contrast_def$name, sep = " | "); cat(format(Sys.time(), "%H:%M"), "fitting", key, "\n", file = stderr())
    block_fit <- if (contrast_def$name == "Interaction") fit_interaction(sub) else fit_one(sub, contrast_def)
    if (!is.null(block_fit$skip)) { skip_list[[key]] <- data.frame(celltype = subtype_name, contrast = contrast_def$name, reason = block_fit$skip); next }
    add_block_keys <- function(block_df) cbind(celltype = subtype_name, contrast = contrast_def$name, block_df)
    adjusted_list[[key]] <- add_block_keys(block_fit$adj); universe_list[[key]] <- add_block_keys(block_fit$uni); summary_list[[key]] <- add_block_keys(block_fit$meta)
    if (!is.null(block_fit$raw)) unadjusted_list[[key]] <- add_block_keys(block_fit$raw)
  }
}
skipped <- bind_rows(data.frame(celltype = character(), contrast = character(), reason = character()), bind_rows(skip_list))
de_all <- bind_rows(adjusted_list); de_all_raw <- bind_rows(unadjusted_list); de_universe <- bind_rows(universe_list); summary_all <- bind_rows(summary_list)
readr::write_csv(de_all, path_de_files[1]); readr::write_csv(de_all_raw, path_de_files[2])
readr::write_csv(summary_all, path_de_files[3]); readr::write_csv(skipped, path_de_files[4])
readr::write_csv(de_universe, path_de_files[5])
classify <- function(block_df) block_df |>
  mutate(celltype = factor(celltype, levels = SUBTYPE_ORDER), contrast = factor(contrast, levels = CONTRAST_ORDER),
         sig = is.finite(p_val_adj) & p_val_adj < PADJ_CUT & abs(avg_log2FC) > LFC_CUT,
         gene_flag = case_when(is_technical_gene(gene) ~ "technical", is_sex_gene(gene) ~ "sex", TRUE ~ ""))
de_all <- classify(de_all)
if (nrow(skipped)) knitr::kable(skipped, caption = "Comparisons not fitted")
de_universe |> mutate(celltype = factor(celltype, levels = SUBTYPE_ORDER), contrast = factor(contrast, levels = CONTRAST_ORDER)) |>
  count(celltype, contrast, name = "universe_n") |>
  pivot_wider(names_from = contrast, values_from = universe_n) |> arrange(celltype) |>
  knitr::kable(caption = sprintf("Genes in the GSEA ranking universe (detected in >= %s%% of the subtype's nuclei) per block.", GSEA_MIN_PCT * 100))

atlas_counts <- GetAssayData(readRDS(path_atlas_rds), assay = "RNA", layer = "counts")
gene_share <- Matrix::rowSums(atlas_counts) / sum(atlas_counts)
bkg_prob   <- cumsum(sort(gene_share))[names(gene_share)]
bkg_genes  <- intersect(unique(de_all$gene), rownames(atlas_counts))
bkg_counts <- atlas_counts[bkg_genes, ]
rm(atlas_counts)
equal_prevalence_ppv <- function(detected, in_subtype) {
  tpr <- Matrix::rowMeans(detected[, in_subtype, drop = FALSE]); fpr <- Matrix::rowMeans(detected[, !in_subtype, drop = FALSE])
  coalesce(tpr / (tpr + fpr), 0)
}
background <- bind_rows(lapply(SUBTYPE_ORDER, function(subtype_name) {
  in_subtype <- colnames(bkg_counts) %in% colnames(ec_obj)[ec_obj$EC_subtype == subtype_name]
  data.frame(celltype = subtype_name, gene = bkg_genes, bkg_prob = unname(bkg_prob[bkg_genes]),
             ppv0 = equal_prevalence_ppv(bkg_counts > 0, in_subtype), ppv1 = equal_prevalence_ppv(bkg_counts > 1, in_subtype))
})) |>
  mutate(nontarget_prob = 1 - (ppv0 + ppv1) / 2, background_score = bkg_prob * nontarget_prob,
         background_suspect = background_score > BKG_CUT)
readr::write_csv(background, file.path(path_out, "ec_DE_background_contamination.csv"))
background |> mutate(celltype = factor(celltype, levels = SUBTYPE_ORDER)) |> group_by(celltype) |> summarise(genes = n(), background_suspect = sum(background_suspect), .groups = "drop") |>
  knitr::kable(caption = sprintf("Genes of the MAST tables flagged background_suspect (bkg_prob x nontarget_prob > %s) per subtype.", BKG_CUT))

animal_stats <- function(subtype_name, genes) {
  sub   <- subset(ec_obj, cells = colnames(ec_obj)[ec_obj$EC_subtype == subtype_name])
  genes <- intersect(genes, rownames(sub))
  expr_data <- GetAssayData(sub, layer = "data")[genes, , drop = FALSE]
  animal_indicator <- fac2sparse(droplevels(sub$individual))
  n  <- rowSums(animal_indicator)
  to_long <- function(m, name) as.data.frame(as.table(t(as.matrix(animal_indicator %*% t(m)) / n))) |> setNames(c("gene", "individual", name))
  to_long(expr_data, "mean_expr") |> left_join(to_long(expr_data > 0, "det_frac"), by = c("gene", "individual")) |>
    mutate(across(c(gene, individual), as.character), celltype = subtype_name, n_nuclei = n[individual])
}
animal_expr <- bind_rows(lapply(SUBTYPE_ORDER, function(subtype_name)
  animal_stats(subtype_name, c(de_universe$gene[de_universe$celltype == subtype_name], de_all$gene[de_all$celltype == subtype_name])))) |>
  left_join(animal_meta |> select(individual, disease, qc_flag), by = "individual")

support_of <- function(keys) {
  keys <- keys |> mutate(across(c(celltype, contrast), as.character))
  pairwise_keys <- keys |> inner_join(bind_rows(PAIRWISE) |> select(contrast = name, i1, i2), by = "contrast") |>
    mutate(hi_sign = sign(avg_log2FC))

  grp_means <- animal_expr |> group_by(celltype, gene, disease) |> summarise(m = mean(mean_expr), .groups = "drop") |>
    pivot_wider(names_from = disease, values_from = m)
  interaction_keys <- keys |> filter(contrast == "Interaction") |> inner_join(grp_means, by = c("celltype", "gene")) |>
    mutate(old_larger = abs(OH - OC) >= abs(YH - YC), i1 = ifelse(old_larger, "OH", "YH"), i2 = ifelse(old_larger, "OC", "YC"),
           hi_sign = sign(ifelse(old_larger, OH - OC, YH - YC))) |>
    select(all_of(names(pairwise_keys)))
  bind_rows(pairwise_keys, interaction_keys) |>
    inner_join(animal_expr, by = c("celltype", "gene"), relationship = "many-to-many") |>
    filter(disease == i1 | disease == i2) |>
    group_by(celltype, contrast, gene) |>
    summarise(n_animals_detected_i1 = sum(disease == i1[1] & det_frac >= DETECT_FRAC),
              n_animals_detected_i2 = sum(disease == i2[1] & det_frac >= DETECT_FRAC),
              direction_concordance = if (contrast[1] == "Interaction") NA_real_ else {
                mean_i1 <- mean(mean_expr[disease == i1[1]]); mean_i2 <- mean(mean_expr[disease == i2[1]])
                mean(sign(ifelse(disease == i1[1], mean_expr - mean_i2, mean_i1 - mean_expr)) == sign(avg_log2FC[1])) },
              max_animal_share = { is_high_side <- disease == (if (hi_sign[1] >= 0) i1[1] else i2[1]); max(mean_expr[is_high_side]) / sum(mean_expr[is_high_side]) },
              max_share_animal = { is_high_side <- disease == (if (hi_sign[1] >= 0) i1[1] else i2[1]); individual[is_high_side][which.max(mean_expr[is_high_side])] },
              low_detection = (if (hi_sign[1] >= 0) n_animals_detected_i1 else n_animals_detected_i2) < MIN_ANIMALS,
              .groups = "drop") |>
    mutate(single_animal_driven = is.finite(max_animal_share) & max_animal_share > MAX_ANIMAL_SHARE)
}
support <- support_of(de_all |> distinct(celltype, contrast, gene, avg_log2FC))
add_support <- function(animal_df) animal_df |>
  left_join(support |> mutate(celltype = factor(celltype, levels = SUBTYPE_ORDER), contrast = factor(contrast, levels = CONTRAST_ORDER)),
            by = c("celltype", "contrast", "gene")) |>
  left_join(background |> select(celltype, gene, bkg_prob, nontarget_prob, background_suspect) |>
              mutate(celltype = factor(celltype, levels = SUBTYPE_ORDER)), by = c("celltype", "gene")) |>
  mutate(reported = sig & gene_flag == "")
de_all <- add_support(de_all)
readr::write_csv(de_all, file.path(path_out, "ec_DE_all_MAST_annotated.csv"))
readr::write_csv(animal_expr, file.path(path_out, "ec_DE_per_animal_expression.csv"))

reported_hits    <- de_all |> filter(reported)
single_animal_by      <- sort(table(reported_hits$max_share_animal[reported_hits$single_animal_driven]), decreasing = TRUE)
low_detection_by      <- sort(table(as.character(reported_hits$celltype[reported_hits$low_detection])), decreasing = TRUE)
note_animal <- names(ANIMAL_NOTES)[1]
note_group  <- animal_meta$disease[animal_meta$individual == note_animal]
ifn_subtype      <- grep("Interferon", SUBTYPE_ORDER, value = TRUE)
ifn_group_nuclei     <- table(as.character(ec_obj$individual[ec_obj$EC_subtype == ifn_subtype & ec_obj$disease == note_group]))
note_contrasts   <- vapply(Filter(function(contrast_def) note_group %in% c(contrast_def$i1, contrast_def$i2), PAIRWISE), `[[`, "", "name")
other_contrasts  <- setdiff(vapply(PAIRWISE, `[[`, "", "name"), note_contrasts)
ifn_note_hits    <- reported_hits |> filter(celltype == ifn_subtype, contrast %in% note_contrasts)
ifn_note_driven  <- sum(ifn_note_hits$single_animal_driven & ifn_note_hits$max_share_animal == note_animal)
ifn_note_pct     <- round(100 * ifn_note_driven / nrow(ifn_note_hits))
ifn_note_top     <- ifn_note_hits |> group_by(contrast) |> summarise(top = sum(max_share_animal == note_animal), n = n(), .groups = "drop")
ifn_note_top_pct <- round(100 * sum(ifn_note_top$top) / nrow(ifn_note_hits))

summary_ann <- summary_all |>
  left_join(de_all |> group_by(celltype, contrast) |>
              summarise(n_hits_reported = sum(reported), n_hits_technical_or_sex = sum(sig & gene_flag != ""),
                        n_single_animal_driven = sum(reported & single_animal_driven),
                        n_low_detection = sum(reported & low_detection),
                        n_background_suspect = sum(reported & background_suspect %in% TRUE), .groups = "drop") |>
              mutate(across(c(celltype, contrast), as.character)),
            by = c("celltype", "contrast"))
summary_ann |>
  transmute(subtype = short_subtype(celltype), contrast = CONTRAST_LABEL[contrast], n_ident1, n_ident2,
            n_animals_ident1, n_animals_ident2, n_hits_raw, n_hits_adj, n_hits_reported,
            n_hits_technical_or_sex, n_single_animal_driven, n_low_detection, n_background_suspect) |>
  knitr::kable(caption = "MAST hits per block: unadjusted fit, covariate-adjusted fit, reported (adjusted minus technical/sex genes), and reported hits flagged single-animal-driven, low-detection or background-suspect. Interaction rows: nuclei/animals = HFpEF vs control of the subtype; no unadjusted fit.")
readr::write_csv(summary_ann, file.path(path_out, "ec_DE_all_summary_annotated.csv"))

pseudobulk_lfc <- readr::read_csv(file.path(PATHS$results, "05b_ec_pseudobulk", "ec_pseudobulk_DE.csv"), show_col_types = FALSE) |>
  select(celltype, contrast, gene, pseudobulk_log2FC = logFC)
direction_rows <- de_all |> filter(reported, contrast != "Interaction") |> mutate(across(c(celltype, contrast), as.character)) |>
  group_by(celltype, contrast) |> arrange(p_val, desc(abs(avg_log2FC)), .by_group = TRUE) |> mutate(leading = row_number() <= LEAD_N) |> ungroup() |>
  inner_join(pseudobulk_lfc, by = c("celltype", "contrast", "gene")) |>
  mutate(same_sign = sign(avg_log2FC) == sign(pseudobulk_log2FC))
direction_tab <- direction_rows |> group_by(contrast = factor(contrast, levels = CONTRAST_ORDER)) |>
  summarise(hits_compared = n(), agreement_all = mean(same_sign), leading_compared = sum(leading), agreement_leading = mean(same_sign[leading]), .groups = "drop")
direction_tab |> transmute(contrast = CONTRAST_LABEL[as.character(contrast)], hits_compared, `sign agreement, all hits (proportion)` = round(agreement_all, 3),
                           leading_compared, `sign agreement, leading hits (proportion)` = round(agreement_leading, 3)) |>
  knitr::kable(caption = sprintf("Sign of the Seurat avg_log2FC vs sign of the pseudobulk DESeq2 log2 fold change (notebook 05b) for the reported pairwise hits and the %d leading hits per block, summed over subtypes.", LEAD_N))

HIT_CLASS <- c("not significant" = "grey88", "technical gene" = "grey40", "sex-chromosome gene" = "#8C564B",
               "background suspect" = FLAG_COL[["background suspect"]], "low detection" = FLAG_COL[["low detection"]],
               "single-animal driven" = FLAG_COL[["single-animal driven"]],
               "reported hit" = "#0072B2")
volcano_df <- de_all |> filter(celltype == VOLCANO_CT) |>
  mutate(class = factor(case_when(!sig ~ "not significant", gene_flag == "technical" ~ "technical gene",
                                  gene_flag == "sex" ~ "sex-chromosome gene", single_animal_driven ~ "single-animal driven",
                                  low_detection ~ "low detection", background_suspect %in% TRUE ~ "background suspect",
                                  TRUE ~ "reported hit"), levels = names(HIT_CLASS)),
         nlp = -log10(pmax(p_val_adj, 1e-300))) |>
  arrange(class)
volcano_labels <- volcano_df |> filter(reported) |> group_by(contrast) |> slice_min(p_val_adj, n = 8)
fig <- ggplot(volcano_df, aes(avg_log2FC, nlp, colour = class)) +
  geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), colour = "grey70", linetype = "dotted") +
  geom_hline(yintercept = -log10(PADJ_CUT), colour = "grey70", linetype = "dotted") +
  geom_point(size = .7, alpha = .8) +
  geom_text_repel(data = volcano_labels, aes(label = gene), size = 2.5, fontface = "italic", max.overlaps = 30, show.legend = FALSE) +
  facet_wrap(~ contrast, nrow = 1, labeller = as_labeller(CONTRAST_LABEL)) +
  scale_colour_manual(values = HIT_CLASS, name = NULL) +
  labs(x = "effect: avg_log2FC (pairwise) | interaction_log2 (Age x HFpEF)", y = expression(-log[10]~"Bonferroni p"),
       title = sprintf("%s: covariate-adjusted MAST screen", VOLCANO_CT),
       subtitle = sprintf("Hit: Bonferroni p < %s & |effect| > %s; labels = 8 smallest p per contrast. Pairwise: avg_log2FC = Seurat log2 ratio of group means (pseudocount 1, unadjusted).\nAge x HFpEF: interaction_log2 = ((OH - OC) - (YH - YC)) of MAST expected log expression / ln 2 (covariate-adjusted; positive = larger HFpEF effect in old)", PADJ_CUT, LFC_CUT)) +
  theme(legend.position = "bottom")
save_svg(fig, "ec_volcano_capillary.svg", width = 13, height = 4.2, dir = path_figures)
fig

flag_mark <- function(sig, sad, low, bkg) paste0(ifelse(sig & sad %in% TRUE, "\u2020", ""), ifelse(sig & low %in% TRUE, "\u00b0", ""),
                                               ifelse(sig & bkg %in% TRUE, "\u2021", ""))
animal_dotplot <- function(genes, subtype_name, title, nrow = 1, expr_df = animal_expr) {
  dot_df <- expr_df |> filter(celltype == subtype_name, gene %in% genes) |>
    mutate(gene = factor(gene, levels = genes), disease = factor(disease, levels = GROUP_ORDER))

  strip_codes <- de_all |> filter(celltype == subtype_name, gene %in% genes, sig) |>
    distinct(gene, contrast, .keep_all = TRUE) |> arrange(contrast) |>
    group_by(gene) |>
    summarise(lab = paste(strwrap(paste0(CONTRAST_CODE[as.character(contrast)], fdr_stars(p_val_adj),
                                         flag_mark(sig, single_animal_driven, low_detection, background_suspect), collapse = " "), width = 14),
                          collapse = "\n"))
  strip_labels <- setNames(genes, genes); strip_labels[strip_codes$gene] <- paste0(strip_codes$gene, "\n", strip_codes$lab)
  ggplot(dot_df, aes(disease, mean_expr)) +
    geom_point(aes(fill = disease, shape = qc_flag == "pass", size = n_nuclei), colour = "grey25", stroke = .4,
               position = position_jitter(width = .12, height = 0, seed = 1)) +
    facet_wrap(~ gene, scales = "free_y", nrow = nrow, labeller = as_labeller(strip_labels)) +
    scale_fill_manual(values = GROUP_FILL, labels = GROUP_LABEL, name = NULL) +
    scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 24), labels = c(`TRUE` = "qc pass", `FALSE` = "qc flagged"), name = NULL) +
    scale_size_continuous(range = c(1.5, 4), name = "nuclei") +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 3))) +
    labs(x = NULL, y = "mean log-normalised expression per animal", title = title,
         subtitle = paste0("Points = animals. Strip: contrasts passing the MAST screen (*/**/*** Bonferroni p < 0.05/0.01/0.001; ",
                           "\u2020 single-animal-driven, \u00b0 low detection, \u2021 background suspect)\n",
                           paste(CONTRAST_CODE, CONTRAST_LABEL, sep = " = ", collapse = "; "))) +
    theme(legend.position = "right", strip.text = element_text(size = 7, face = "bold"),
          axis.text.x = element_text(size = 7))
}
top_hits <- de_all |> filter(celltype == VOLCANO_CT, reported) |> group_by(contrast) |> slice_min(p_val_adj, n = 3) |>
  ungroup() |> distinct(gene) |> pull(gene)
fig <- animal_dotplot(top_hits, VOLCANO_CT, sprintf("Top reported MAST hits (3 per contrast), %s", VOLCANO_CT),
                    nrow = ceiling(length(top_hits) / 8))
save_svg(fig, "ec_top_hits_per_animal.svg", width = 13, height = 6.8, dir = path_figures)
fig

path_secretome  <- file.path(PATHS$resources, "mouse_secretome_swissprot.csv")
secretome <- readr::read_csv(path_secretome, show_col_types = FALSE)
secretome_genes <- secretome$gene
secreted_hits <- de_all |> filter(reported, gene %in% secretome_genes) |>
  left_join(secretome |> select(gene, uniprot), by = "gene")
cat(sprintf("Secretome: %d Swiss-Prot genes (retrieved %s); %d present in the assay; %d screened (tested in >= 1 block)\n",
            length(secretome_genes), secretome$retrieved_on[1], sum(secretome_genes %in% rownames(ec_obj)), n_distinct(de_all$gene[de_all$gene %in% secretome_genes])))
cat(sprintf("%d reported secreted gene-block hits across %d unique genes\n", nrow(secreted_hits), n_distinct(secreted_hits$gene)))

secreted_counts <- expand_grid(celltype = SUBTYPE_ORDER, contrast = CONTRAST_ORDER) |>
  semi_join(summary_all, by = c("celltype", "contrast")) |>
  left_join(secreted_hits |> mutate(across(c(celltype, contrast), as.character)) |> count(celltype, contrast, name = "n_secreted"),
            by = c("celltype", "contrast")) |>
  mutate(n_secreted = coalesce(n_secreted, 0L), celltype = factor(celltype, levels = SUBTYPE_ORDER),
         contrast = factor(contrast, levels = CONTRAST_ORDER))
p_counts <- ggplot(secreted_counts, aes(contrast, n_secreted, fill = contrast)) +
  geom_col(colour = "grey30", linewidth = .2) +
  geom_text(aes(label = n_secreted), vjust = -.3, size = 2.8) +
  facet_wrap(~ celltype, nrow = 1, labeller = as_labeller(short_subtype)) +
  scale_fill_manual(values = CONTRAST_FILL, labels = CONTRAST_LABEL, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, .15))) +
  labs(x = NULL, y = "reported secreted-gene hits", title = "Secreted-protein gene hits per comparison") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom")

recurrent_genes <- secreted_hits |> count(gene, name = "n_blocks") |> filter(n_blocks >= 2)
recurrent_df <- secreted_hits |> semi_join(recurrent_genes, by = "gene") |> left_join(recurrent_genes, by = "gene") |>
  mutate(gene = reorder(gene, n_blocks), lfc_plot = pmax(pmin(avg_log2FC, LFC_SHOW), -LFC_SHOW),
         n_flags = single_animal_driven + low_detection + (background_suspect %in% TRUE),
         flag = factor(case_when(n_flags > 1 ~ "several flags", single_animal_driven ~ "single-animal driven",
                                 low_detection ~ "low detection", background_suspect %in% TRUE ~ "background suspect",
                                 TRUE ~ "none"), levels = names(FLAG_COL)))
p_recurrent <- ggplot(recurrent_df, aes(lfc_plot, gene, fill = contrast, shape = celltype)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = .4) +
  geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), colour = "grey85", linetype = "dotted", linewidth = .4) +
  geom_point(aes(colour = flag), size = 2.6, stroke = .6, show.legend = TRUE) +
  geom_text(data = recurrent_df |> filter(abs(avg_log2FC) > LFC_SHOW), aes(label = sprintf("%.1f", avg_log2FC)),
            nudge_y = .45, size = 2.2, colour = "grey30", show.legend = FALSE) +
  scale_x_continuous(limits = c(-LFC_SHOW, LFC_SHOW)) +
  scale_shape_manual(values = setNames(FILL_SHAPES[seq_along(SUBTYPE_ORDER)], SUBTYPE_ORDER), labels = short_subtype, name = "subtype") +
  scale_fill_manual(values = CONTRAST_FILL, labels = CONTRAST_LABEL, name = "comparison") +
  scale_colour_manual(values = FLAG_COL, name = "outline: flag", drop = FALSE) +
  guides(fill = guide_legend(override.aes = list(shape = 21)), shape = guide_legend(override.aes = list(fill = "grey80")),
         colour = guide_legend(override.aes = list(shape = 21, fill = "white", size = 3))) +
  labs(x = "effect: avg_log2FC (pairwise) | interaction_log2 (Age x HFpEF)", y = NULL, title = "Secreted-protein genes reported in >= 2 blocks",
       subtitle = sprintf("|effect| > %s drawn at the axis limit and labelled with its value", LFC_SHOW)) +
  theme(axis.text.y = element_text(face = "italic"))
fig <- p_counts + p_recurrent + plot_layout(widths = c(1.2, 1))
save_svg(fig, "ec_secreted_degs.svg", width = 13, height = 6, dir = path_figures)
fig

secreted_hits |> filter(celltype == VOLCANO_CT) |>
  transmute(comparison = CONTRAST_LABEL[as.character(contrast)], gene, effect_scale = ifelse(contrast == "Interaction", "interaction_log2", "avg_log2FC"), effect = round(avg_log2FC, 3),
            p_bonferroni = format_pvalue(p_val_adj), `pct.1 (proportion)` = pct.1, `pct.2 (proportion)` = pct.2,
            animals_detected = paste(n_animals_detected_i1, n_animals_detected_i2, sep = " vs "),
            `direction_concordance (proportion)` = direction_concordance, single_animal_driven, low_detection, background_suspect, uniprot) |>
  arrange(comparison, effect) |>
  knitr::kable(caption = sprintf("%s secreted-protein genes passing the reported screen.", VOLCANO_CT))

TARGET_GENES <- readLines(file.path(PROJECT_ROOT, "config", "candidate_genes.txt"))
target_status <- function(fitted, in_assay, n_detecting, screened, p_val)
  case_when(!fitted ~ "comparison skipped", !in_assay ~ "absent from assay", !(n_detecting > 0) ~ "not detected",
            !screened ~ "below 5% floor, not tested", !is.finite(p_val) ~ "not testable", TRUE ~ "tested")
target_detected <- GetAssayData(ec_obj, layer = "data")[intersect(TARGET_GENES, rownames(ec_obj)), , drop = FALSE] > 0
block_groups    <- c(setNames(lapply(PAIRWISE, function(contrast_def) c(contrast_def$i1, contrast_def$i2)), vapply(PAIRWISE, `[[`, "", "name")),
                     list(Interaction = GROUP_ORDER))
n_detecting_in_block <- function(subtype_name, contrast_name, gene)
  sum(target_detected[intersect(gene, rownames(target_detected)), ec_obj$EC_subtype == subtype_name & ec_obj$disease %in% block_groups[[contrast_name]], drop = FALSE])
screen_rows <- de_all |> select(celltype, contrast, gene, p_val, avg_log2FC, pct.1, pct.2, p_val_adj) |>
  mutate(across(c(celltype, contrast), as.character)) |> filter(gene %in% TARGET_GENES)
target_all <- expand_grid(celltype = SUBTYPE_ORDER, contrast = vapply(CONTRASTS, `[[`, "", "name"), gene = TARGET_GENES) |>
  left_join(screen_rows, by = c("celltype", "contrast", "gene")) |>
  mutate(status = target_status(paste(celltype, contrast) %in% paste(summary_all$celltype, summary_all$contrast), gene %in% rownames(ec_obj),
                                mapply(n_detecting_in_block, celltype, contrast, gene, USE.NAMES = FALSE),
                                paste(celltype, contrast, gene) %in% paste(screen_rows$celltype, screen_rows$contrast, screen_rows$gene), p_val))
readr::write_csv(target_all, file.path(path_out, "ec_DE_target_genes.csv"))
target_all <- add_support(classify(target_all))
readr::write_csv(target_all, file.path(path_out, "ec_DE_target_genes_annotated.csv"))
animal_expr_target <- animal_stats(VOLCANO_CT, TARGET_GENES) |> left_join(animal_meta |> select(individual, disease, qc_flag), by = "individual")

target_plot_df <- target_all |>
  mutate(gene = factor(gene, levels = TARGET_GENES),
         mark = paste0(ifelse(sig, as.character(fdr_stars(p_val_adj)), ""), flag_mark(sig, single_animal_driven, low_detection, background_suspect)),
         lfc_capped = pmax(pmin(avg_log2FC, 1.2), -1.2))
fig <- ggplot(target_plot_df, aes(contrast, celltype, fill = lfc_capped)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = mark, colour = abs(lfc_capped) > .6), fontface = "bold", size = 3.6, show.legend = FALSE) +
  facet_wrap(~ gene, nrow = 1) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, limits = c(-1.2, 1.2),
                       na.value = "grey85", breaks = c(-1.2, -0.6, 0, 0.6, 1.2), name = "effect\navg_log2FC |\ninteraction_log2") +
  scale_colour_manual(values = c(`FALSE` = "grey15", `TRUE` = "white")) +
  scale_x_discrete(labels = CONTRAST_LABEL_WRAP) +
  scale_y_discrete(labels = short_subtype, limits = rev(SUBTYPE_ORDER)) +
  labs(x = NULL, y = NULL, title = "Target-gene expression contrasts in endothelial subtypes",
       subtitle = sprintf("*/**/*** Bonferroni p < 0.05/0.01/0.001 with |effect| > %s; \u2020 single-animal-driven; \u00b0 low detection; \u2021 background suspect; grey = unavailable.\nPairwise fill = avg_log2FC (Seurat log2 ratio of group means); Age x HFpEF fill = interaction_log2 (MAST expected log expression, difference of differences / ln 2; positive = larger HFpEF effect in old)", LFC_CUT)) +
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.ticks = element_blank())
save_svg(fig, "ec_target_genes_heatmap.svg", width = 12, height = 3.4, dir = path_figures)
fig

target_plot_df |> filter(celltype == VOLCANO_CT) |>
  transmute(gene, comparison = CONTRAST_LABEL[as.character(contrast)], effect_scale = ifelse(contrast == "Interaction", "interaction_log2", "avg_log2FC"), effect = round(avg_log2FC, 3),
            p_bonferroni = format_pvalue(p_val_adj), `pct.1 (proportion)` = pct.1, `pct.2 (proportion)` = pct.2,
            animals_detected = paste(n_animals_detected_i1, n_animals_detected_i2, sep = " vs "),
            `direction_concordance (proportion)` = direction_concordance, `max_animal_share (proportion)` = round(max_animal_share, 2), single_animal_driven, low_detection, background_suspect,
            passes_screen = sig, status) |>
  arrange(gene, comparison) |>
  knitr::kable(caption = sprintf("%s target results with per-animal support; pairwise effects = avg_log2FC (FindMarkers log2 ratio of group means), Age x HFpEF effect = interaction_log2 (MAST expected log expression, difference of differences / ln 2).", VOLCANO_CT))

target_plot_df |>
  transmute(gene, subtype = short_subtype(as.character(celltype)), comparison = CONTRAST_LABEL[as.character(contrast)],
            effect_scale = ifelse(contrast == "Interaction", "interaction_log2", "avg_log2FC"), effect = round(avg_log2FC, 3), p_bonferroni = format_pvalue(p_val_adj),
            animals_detected = paste(n_animals_detected_i1, n_animals_detected_i2, sep = " vs "),
            `direction_concordance (proportion)` = direction_concordance, single_animal_driven, low_detection, background_suspect, passes_screen = sig, status) |>
  arrange(gene, subtype, comparison)

fig <- animal_dotplot(TARGET_GENES, VOLCANO_CT, sprintf("Target genes, %s: per-animal expression", VOLCANO_CT), expr_df = animal_expr_target)
save_svg(fig, "ec_target_genes_per_animal.svg", width = 11, height = 4.2, dir = path_figures)
fig

target_plot_df |> group_by(gene) |>
  summarise(tested_blocks = sum(status == "tested"), passing_blocks = sum(sig, na.rm = TRUE),
            passing_single_animal_driven = sum(sig & single_animal_driven %in% TRUE),
            passing_low_detection = sum(sig & low_detection %in% TRUE),
            passing_background_suspect = sum(sig & background_suspect %in% TRUE), .groups = "drop") |>
  knitr::kable(caption = sprintf("Target screening results across %d blocks per gene.",
                                 n_distinct(target_plot_df$celltype) * n_distinct(target_plot_df$contrast)))

session_footer()
