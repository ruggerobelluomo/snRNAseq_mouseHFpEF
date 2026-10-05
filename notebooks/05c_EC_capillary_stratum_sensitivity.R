# 05c_EC_capillary_stratum_sensitivity.R
# Code-only export of the analysis pipeline. Data are not included; see README.md.

params <- list(rds_in = "results/03_ec_annotation/Endothelial_subcluster_annotated_clean.rds", 
    mast_csv = "results/05_ec_de/ec_DE_all_MAST_annotated.csv", gsea_rds = "results/06_ec_gsea/ec_gsea_mast.rds", 
    out_dir = "results/05c_ec_capillary_sensitivity")


knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = TRUE,
                      fig.align = "center", dpi = 150)
options(width = 110)

suppressPackageStartupMessages({
  library(Seurat); library(MAST); library(Matrix); library(tidyr); library(patchwork); library(data.table); library(fgsea)
})
source("_common.R")
knitr::opts_chunk$set(cache = TRUE, autodep = TRUE, cache.lazy = FALSE, cache.extra = input_md5(params))
path_rds_in   <- resolve_path(params$rds_in)
path_mast_csv <- resolve_path(params$mast_csv)
path_gsea_rds <- resolve_path(params$gsea_rds)
path_out      <- resolve_path(params$out_dir)
path_figures  <- file.path(path_out, "figures"); dir.create(path_figures, recursive = TRUE, showWarnings = FALSE)

nb_chunk <- function(notebook, label) {
  rmd   <- readLines(file.path(PROJECT_ROOT, "notebooks", notebook))
  start <- grep(sprintf("^```\\{r %s[,}]", label), rmd)
  end   <- start + match(TRUE, grepl("^```\\s*$", rmd[-seq_len(start)]))
  rmd[(start + 1):(end - 1)]
}
knitr::knit_code$set(`fit-function` = nb_chunk("05_EC_DE_analysis.Rmd", "fit-function"))
nb06_params <- rmarkdown::yaml_front_matter(file.path(PROJECT_ROOT, "notebooks", "06_EC_MAST_GSEA_analysis.Rmd"))$params
nb06_env    <- list2env(list(params = nb06_params), parent = environment())
for (fn_def in Filter(function(e) is.call(e) && identical(e[[1]], as.name("<-")) && deparse(e[[2]]) %in% c("make_rank", "run_fgsea"),
                      parse(text = nb_chunk("06_EC_MAST_GSEA_analysis.Rmd", "run-gsea")))) eval(fn_def, nb06_env)
make_rank <- nb06_env$make_rank; run_fgsea <- nb06_env$run_fgsea

PADJ_CUT     <- 0.05
LFC_CUT      <- 0.2
MIN_CELLS    <- 30
MIN_ANIMALS  <- 3
MAST_COVARIATES <- c(DESIGN_COVARIATES, "cngeneson", "exon_prop")
GSEA_MIN_PCT <- 0.01

FDR_CUT  <- 0.05
MIN_SIZE <- 10; MAX_SIZE <- 500
SEED     <- 1

CAPILLARY <- "Capillary EC"
STRATA    <- c(low_content = "Capillary EC low-content stratum", cytoplasm_rich = "Capillary EC cytoplasm-rich stratum")
REFIT_LABEL <- c(no_strata = "without both strata", no_cytoplasm_rich = "without the cytoplasm-rich stratum")
REFIT_STRIP <- c(no_strata = "without\nboth strata", no_cytoplasm_rich = "without the\ncytoplasm-rich\nstratum")
FIT_LABEL   <- c(full = "all nuclei (notebooks 05/06)", REFIT_LABEL)
FIT_FILL    <- c(full = "grey45", no_strata = "#0072B2", no_cytoplasm_rich = "#E69F00")
TARGET_GENES <- readLines(file.path(PROJECT_ROOT, "config", "candidate_genes.txt"))
N_LEADING_SETS <- 5
LFC_CAP <- 2
span <- function(x, fmt = "%.1f") paste(unique(sprintf(fmt, range(x))), collapse = "–")   #

print_versions(c("ggplot2", "dplyr", "Seurat", "MAST", "Matrix", "tidyr", "patchwork", "data.table", "fgsea"))

ec_obj <- readRDS(path_rds_in)
DefaultAssay(ec_obj) <- "RNA"
cap_obj <- subset(ec_obj, cells = colnames(ec_obj)[ec_obj$EC_subtype == CAPILLARY])
rm(ec_obj)
cap_obj$qc_flag <- unname(coalesce(QC_FLAGGED_ANIMALS[as.character(cap_obj$sample)], "pass"))
PAIRWISE  <- CONTRASTS
CONTRASTS <- c(PAIRWISE, list(list(name = "Interaction", i1 = NA, i2 = NA, meaning = "age x HFpEF interaction, (OH - OC) - (YH - YC)")))
stratum_subcluster <- vapply(STRATA, function(stratum) paste(unique(as.character(cap_obj$merged_cluster[cap_obj$EC_label == stratum])), collapse = "+"), "")
strata_tab <- cap_obj@meta.data |>
  mutate(disease = factor(disease, levels = GROUP_ORDER)) |>
  group_by(disease, individual, qc_flag) |>
  summarise(capillary_nuclei = n(), low_content = sum(EC_label == STRATA[["low_content"]]),
            cytoplasm_rich = sum(EC_label == STRATA[["cytoplasm_rich"]]), .groups = "drop") |>
  mutate(strata_pct = round(100 * (low_content + cytoplasm_rich) / capillary_nuclei, 1),
         cytoplasm_rich_pct = round(100 * cytoplasm_rich / capillary_nuclei, 1))
readr::write_csv(strata_tab, file.path(path_out, "ec_capillary_strata_per_animal.csv"))
strata_group <- strata_tab |> group_by(disease) |>
  summarise(capillary_nuclei = sum(capillary_nuclei), strata_pct = 100 * sum(low_content + cytoplasm_rich) / capillary_nuclei,
            cytoplasm_rich_pct = 100 * sum(cytoplasm_rich) / capillary_nuclei, .groups = "drop")
knitr::kable(strata_tab |> rename(group = disease, `strata (%)` = strata_pct, `cytoplasm-rich (%)` = cytoplasm_rich_pct),
             caption = "Capillary EC nuclei per animal and nuclei of the two technical strata (percentages of the animal's Capillary EC nuclei).")

animal_weight <- strata_tab |> group_by(disease) |>
  mutate(weight_full_pct = 100 * capillary_nuclei / sum(capillary_nuclei),
         weight_no_strata_pct = 100 * (capillary_nuclei - low_content - cytoplasm_rich) / sum(capillary_nuclei - low_content - cytoplasm_rich),
         weight_no_cytoplasm_rich_pct = 100 * (capillary_nuclei - cytoplasm_rich) / sum(capillary_nuclei - cytoplasm_rich)) |>
  ungroup() |> mutate(across(starts_with("weight"), ~ round(.x, 1)))
readr::write_csv(animal_weight |> select(disease, individual, starts_with("weight")), file.path(path_out, "ec_capillary_animal_weight.csv"))
knitr::kable(animal_weight |> transmute(group = disease, individual, `all nuclei (%)` = weight_full_pct, `without both strata (%)` = weight_no_strata_pct,
                                        `without the cytoplasm-rich stratum (%)` = weight_no_cytoplasm_rich_pct),
             caption = "Share of each animal in its group's Capillary EC nuclei (%) in each fit.")
weight_change <- animal_weight |> mutate(change = weight_no_strata_pct - weight_full_pct) |> slice_max(abs(change), n = 1)

fit_capillary <- function(excluded) {
  sub <- subset(cap_obj, cells = colnames(cap_obj)[!cap_obj$EC_label %in% excluded])
  block_fits <- setNames(lapply(CONTRASTS, function(contrast_def) {
    cat(format(Sys.time(), "%H:%M"), "fitting", CAPILLARY, "without", paste(excluded, collapse = " + "), "|", contrast_def$name, "\n", file = stderr())
    if (contrast_def$name == "Interaction") fit_interaction(sub) else fit_one(sub, contrast_def)
  }), vapply(CONTRASTS, `[[`, "", "name"))
  collect <- function(part) bind_rows(lapply(names(block_fits), function(contrast_name)
    cbind(celltype = CAPILLARY, contrast = contrast_name, block_fits[[contrast_name]][[part]])))
  list(adj = collect("adj"), uni = collect("uni"), meta = collect("meta"))
}
write_refit <- function(refit, name) {
  readr::write_csv(refit$adj,  file.path(path_out, sprintf("ec_capillary_DE_MAST_%s.csv", name)))
  readr::write_csv(refit$uni,  file.path(path_out, sprintf("ec_capillary_DE_MAST_gsea_universe_%s.csv", name)))
  readr::write_csv(refit$meta, file.path(path_out, sprintf("ec_capillary_DE_summary_%s.csv", name)))
}

t_start <- Sys.time()
refit_no_strata <- fit_capillary(unname(STRATA))
minutes_no_strata <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
write_refit(refit_no_strata, "no_strata")

t_start <- Sys.time()
refit_no_cytoplasm_rich <- fit_capillary(STRATA[["cytoplasm_rich"]])
minutes_no_cytoplasm_rich <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
write_refit(refit_no_cytoplasm_rich, "no_cytoplasm_rich")

REFITS <- list(no_strata = refit_no_strata, no_cytoplasm_rich = refit_no_cytoplasm_rich)
full_de <- readr::read_csv(path_mast_csv, show_col_types = FALSE) |> filter(celltype == CAPILLARY)
fit_de  <- bind_rows(c(list(full = full_de |> select(names(refit_no_strata$adj))), lapply(REFITS, `[[`, "adj")), .id = "fit") |>
  mutate(reported = is.finite(p_val_adj) & p_val_adj < PADJ_CUT & abs(avg_log2FC) > LFC_CUT & !is_technical_gene(gene) & !is_sex_gene(gene),
         fit = factor(fit, levels = names(FIT_LABEL)), contrast = factor(contrast, levels = CONTRAST_ORDER)) |>
  select(fit, contrast, gene, p_val, avg_log2FC, pct.1, pct.2, p_val_adj, reported)
block_nuclei <- bind_rows(lapply(REFITS, `[[`, "meta"), .id = "fit") |> transmute(fit, contrast, nuclei = n_ident1 + n_ident2, n_hits_raw, n_hits_adj)
full_nuclei  <- lapply(CONTRASTS, function(contrast_def) {
  groups <- if (contrast_def$name == "Interaction") GROUP_ORDER else c(contrast_def$i1, contrast_def$i2)
  data.frame(contrast = contrast_def$name, nuclei_full = sum(cap_obj$disease %in% groups))
}) |> bind_rows()
refit_summary <- fit_de |> count(fit, contrast, reported) |> filter(reported) |> select(fit, contrast, n_hits_reported = n) |>
  mutate(across(c(fit, contrast), as.character)) |>
  left_join(block_nuclei, by = c("fit", "contrast")) |> left_join(full_nuclei, by = "contrast") |>
  mutate(nuclei = coalesce(nuclei, nuclei_full), nuclei_kept_pct = round(100 * nuclei / nuclei_full, 1))
readr::write_csv(refit_summary, file.path(path_out, "ec_capillary_refit_summary.csv"))
refit_summary |>
  transmute(fit = FIT_LABEL[fit], contrast = CONTRAST_LABEL[contrast], nuclei, `nuclei kept (%)` = nuclei_kept_pct,
            n_hits_raw, n_hits_adj, n_hits_reported) |>
  knitr::kable(caption = "Nuclei and hits per fit and contrast. Hits: Bonferroni p < 0.05 and |effect| > 0.2 in the unadjusted and the covariate-adjusted fit; reported = adjusted hits without technical and sex-chromosome genes. The unadjusted and adjusted counts of all nuclei are in notebook 05.")

nlp <- function(p) -log10(pmax(p, .Machine$double.xmin))
paired_de <- bind_rows(lapply(names(REFITS), function(refit_name)
  full_join(fit_de |> filter(fit == "full") |> select(-fit), fit_de |> filter(fit == refit_name) |> select(-fit),
            by = c("contrast", "gene"), suffix = c("_full", "_refit")) |> mutate(refit = refit_name))) |>
  mutate(hit_full = reported_full %in% TRUE, hit_refit = reported_refit %in% TRUE,
         refit = factor(refit, levels = names(REFIT_LABEL)))
hit_cmp <- paired_de |> group_by(refit, contrast) |>
  summarise(hits_full = sum(hit_full), hits_refit = sum(hit_refit), shared = sum(hit_full & hit_refit),
            full_only = sum(hit_full & !hit_refit), refit_only = sum(!hit_full & hit_refit),
            full_hits_tested = sum(hit_full & is.finite(p_val_refit)),
            direction_agree_pct = 100 * mean((sign(avg_log2FC_refit) == sign(avg_log2FC_full))[hit_full & is.finite(p_val_refit)]),
            r_effect = cor(avg_log2FC_full[hit_full], avg_log2FC_refit[hit_full], use = "complete.obs"),
            rho_nlp = cor(nlp(p_val_full[hit_full]), nlp(p_val_refit[hit_full]), method = "spearman", use = "complete.obs"),
            nlp_ratio = median((nlp(p_val_refit) / nlp(p_val_full))[hit_full], na.rm = TRUE), .groups = "drop") |>
  mutate(shared_pct = 100 * shared / hits_full) |>
  left_join(refit_summary |> select(refit = fit, contrast, nuclei_kept_pct) |> mutate(refit = factor(refit, levels = names(REFIT_LABEL)),
                                                                                     contrast = factor(contrast, levels = CONTRAST_ORDER)),
            by = c("refit", "contrast"))
readr::write_csv(hit_cmp, file.path(path_out, "ec_capillary_hit_agreement.csv"))
readr::write_csv(paired_de, file.path(path_out, "ec_capillary_paired_statistics.csv"))
hit_cmp |>
  transmute(refit = REFIT_LABEL[as.character(refit)], contrast = CONTRAST_LABEL[as.character(contrast)], hits_full, hits_refit, shared,
            `shared (% of full)` = round(shared_pct, 1), full_only, refit_only,
            `direction agreement (%)` = round(direction_agree_pct, 1), r_effect = round(r_effect, 3),
            rho_nlp = round(rho_nlp, 3), nlp_ratio = round(nlp_ratio, 2), `nuclei kept (%)` = nuclei_kept_pct) |>
  knitr::kable(caption = "Reported hits of all nuclei (full) and of each refit per contrast; for the full-data hits: direction agreement (% of hits tested in the refit), Pearson correlation of the effect (r_effect), Spearman correlation of -log10 p (rho_nlp) and median ratio of -log10 p, refit / full (nlp_ratio).")

hc <- function(refit_name, column) hit_cmp[[column]][hit_cmp$refit == refit_name]
hcv <- function(refit_name, contrast_name, column) hit_cmp[[column]][hit_cmp$refit == refit_name & hit_cmp$contrast == contrast_name]
pairwise_names <- vapply(PAIRWISE, `[[`, "", "name")
n_flipped <- function(refit_name) sum(round(hc(refit_name, "full_hits_tested") * (1 - hc(refit_name, "direction_agree_pct") / 100)))
lost_same_sign <- paired_de |> filter(hit_full, !hit_refit, is.finite(p_val_refit)) |>
  group_by(refit) |> summarise(n = n(), same_sign = sum(sign(avg_log2FC_refit) == sign(avg_log2FC_full)),
                               below_lfc = sum(abs(avg_log2FC_refit) <= LFC_CUT), .groups = "drop")
lst <- function(refit_name, column) lost_same_sign[[column]][lost_same_sign$refit == refit_name]

hit_bars <- hit_cmp |> select(refit, contrast, shared, full_only, refit_only) |>
  pivot_longer(c(shared, full_only, refit_only), names_to = "class", values_to = "n") |>
  mutate(n = ifelse(class == "full_only", -n, n),
         class = factor(class, levels = c("full_only", "refit_only", "shared"),
                        labels = c("all nuclei only (drawn below 0)", "refit only", "both")))
fig <- ggplot(hit_bars, aes(contrast, n, fill = class)) +
  geom_col(colour = "grey30", linewidth = .2, width = .7) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = .3) +
  facet_wrap(~ refit, nrow = 1, labeller = as_labeller(REFIT_LABEL)) +
  scale_fill_manual(values = c("grey70", "#E69F00", "#0072B2"), name = NULL) +
  scale_x_discrete(labels = CONTRAST_LABEL_WRAP) +
  scale_y_continuous(labels = abs) +
  labs(x = NULL, y = "reported hits", title = sprintf("%s: reported MAST hits of all nuclei and of each refit", CAPILLARY),
       subtitle = sprintf("Hit: Bonferroni p < %s, |effect| > %s, technical and sex-chromosome genes removed", PADJ_CUT, LFC_CUT)) +
  theme(legend.position = "bottom")
save_svg(fig, "ec_capillary_hit_overlap.svg", width = 12, height = 4, dir = path_figures)
fig

scatter_df <- paired_de |> filter(hit_full, is.finite(p_val_refit)) |>
  mutate(status = factor(ifelse(hit_refit, "reported in the refit", "not reported in the refit"),
                         levels = c("reported in the refit", "not reported in the refit")))
scatter_text <- hit_cmp |> transmute(refit, contrast, lab = sprintf("r = %.2f\nrho(-log10 p) = %.2f", r_effect, rho_nlp))
refit_facets <- facet_grid(refit ~ contrast, labeller = labeller(refit = REFIT_STRIP, contrast = CONTRAST_LABEL), scales = "free")
p_effect <- ggplot(scatter_df, aes(avg_log2FC_full, avg_log2FC_refit, colour = status)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60", linetype = "dashed", linewidth = .3) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = .3) + geom_vline(xintercept = 0, colour = "grey80", linewidth = .3) +
  geom_point(size = .5, alpha = .6) +
  geom_text(data = scatter_text, aes(x = -Inf, y = Inf, label = lab), inherit.aes = FALSE, hjust = -.05, vjust = 1.2, size = 2.6) +
  refit_facets +
  scale_colour_manual(values = c("#0072B2", "#D55E00"), name = NULL) +
  labs(x = "effect, all nuclei", y = "effect, refit",
       title = sprintf("%s: effects of the full-data hits in each refit", CAPILLARY),
       subtitle = "Effect = avg_log2FC of notebook 05 (pairwise: log2 ratio of group means; Age x HFpEF: difference of the HFpEF effects on the MAST hurdle scale); dashed line = identity") +
  theme(legend.position = "bottom")
save_svg(p_effect, "ec_capillary_effect_agreement.svg", width = 14, height = 6.5, dir = path_figures)
p_effect

nlp_ref <- hit_cmp |> transmute(refit, contrast, slope = nuclei_kept_pct / 100)
p_nlp <- ggplot(scatter_df, aes(nlp(p_val_full), nlp(p_val_refit), colour = status)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60", linetype = "dashed", linewidth = .3) +
  geom_abline(data = nlp_ref, aes(slope = slope, intercept = 0), colour = "grey30", linetype = "dotted", linewidth = .4) +
  geom_point(size = .5, alpha = .6) +
  refit_facets +
  scale_colour_manual(values = c("#0072B2", "#D55E00"), name = NULL) +
  labs(x = expression(-log[10]~"p, all nuclei"), y = expression(-log[10]~"p, refit"),
       title = sprintf("%s: evidence for the full-data hits in each refit", CAPILLARY),
       subtitle = "Dashed line = identity; dotted line = slope equal to the proportion of nuclei kept (expected -log10 p at an unchanged effect)") +
  theme(legend.position = "bottom")
save_svg(p_nlp, "ec_capillary_pvalue_agreement.svg", width = 14, height = 6.5, dir = path_figures)
p_nlp

leading <- full_de |> filter(reported) |> group_by(contrast) |> slice_min(p_val_adj, n = 3) |> ungroup() |>
  transmute(gene, leading_contrast = contrast)
gene_order <- c(unique(leading$gene), setdiff(TARGET_GENES, leading$gene))
gene_status <- expand_grid(fit = factor(names(FIT_LABEL), levels = names(FIT_LABEL)), contrast = factor(CONTRAST_ORDER, levels = CONTRAST_ORDER), gene = gene_order) |>
  left_join(fit_de, by = c("fit", "contrast", "gene")) |>
  mutate(status = ifelse(is.finite(p_val), ifelse(reported, "reported", "tested, not reported"), "below the 5% floor, not tested"),
         gene_class = ifelse(gene %in% leading$gene, "leading gene", "candidate gene"))
readr::write_csv(gene_status, file.path(path_out, "ec_capillary_leading_and_candidate_genes.csv"))
lead_wide <- leading |> left_join(gene_status, by = c("gene", "leading_contrast" = "contrast")) |>
  transmute(gene, contrast = CONTRAST_LABEL[as.character(leading_contrast)], fit,
            value = sprintf("%.2f (%s)%s", avg_log2FC, format_pvalue(p_val_adj), ifelse(reported, "", " not reported"))) |>
  pivot_wider(names_from = fit, values_from = value)
knitr::kable(lead_wide |> rename_with(~ FIT_LABEL[.x], any_of(names(FIT_LABEL))),
             caption = "Leading capillary genes in the contrast where they lead: effect (Bonferroni p) per fit.")
cand_wide <- gene_status |> filter(gene %in% TARGET_GENES) |>
  transmute(gene, contrast = CONTRAST_LABEL[as.character(contrast)], fit,
            value = ifelse(is.finite(p_val), sprintf("%.2f (%s)%s", avg_log2FC, format_pvalue(p_val_adj), ifelse(reported, " *", "")), "not tested")) |>
  pivot_wider(names_from = fit, values_from = value) |> arrange(gene)
knitr::kable(cand_wide |> rename_with(~ FIT_LABEL[.x], any_of(names(FIT_LABEL))),
             caption = "Candidate genes in every contrast: effect (Bonferroni p) per fit; * = reported hit.")

lead_status <- leading |> left_join(gene_status |> filter(fit != "full"), by = c("gene", "leading_contrast" = "contrast"))
lead_kept   <- tapply(lead_status$reported %in% TRUE, lead_status$fit, sum)[names(REFIT_LABEL)]
cand_status <- gene_status |> filter(gene %in% TARGET_GENES) |> select(fit, contrast, gene, reported, avg_log2FC) |>
  pivot_wider(names_from = fit, values_from = c(reported, avg_log2FC))
cand_full_hits <- cand_status |> filter(reported_full %in% TRUE)
cand_kept <- vapply(names(REFIT_LABEL), function(refit_name) sum(cand_full_hits[[paste0("reported_", refit_name)]] %in% TRUE), 0L)
cand_new  <- vapply(names(REFIT_LABEL), function(refit_name) sum(cand_status[[paste0("reported_", refit_name)]] %in% TRUE & !cand_status$reported_full %in% TRUE), 0L)
cand_sign <- vapply(names(REFIT_LABEL), function(refit_name) sum(sign(cand_full_hits[[paste0("avg_log2FC_", refit_name)]]) == sign(cand_full_hits$avg_log2FC_full), na.rm = TRUE), 0L)

heat_df <- gene_status |>
  mutate(gene = factor(gene, levels = rev(gene_order)), effect_capped = pmax(pmin(avg_log2FC, LFC_CAP), -LFC_CAP),
         mark = ifelse(reported %in% TRUE, as.character(fdr_stars(p_val_adj)), ""))
lead_cells <- leading |> mutate(gene = factor(gene, levels = rev(gene_order)), contrast = factor(leading_contrast, levels = CONTRAST_ORDER), gene_class = "leading gene") |>
  expand_grid(fit = factor(names(FIT_LABEL), levels = names(FIT_LABEL)))
fig <- ggplot(heat_df, aes(contrast, gene, fill = effect_capped)) +
  geom_tile(colour = "white", linewidth = .8) +
  geom_tile(data = lead_cells, fill = NA, colour = "black", linewidth = .5) +
  geom_text(aes(label = mark, colour = abs(effect_capped) > LFC_CAP / 2), size = 3, fontface = "bold", show.legend = FALSE) +
  facet_grid(gene_class ~ fit, scales = "free_y", space = "free_y", labeller = labeller(fit = FIT_LABEL)) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, limits = c(-LFC_CAP, LFC_CAP),
                       na.value = "grey85", name = "effect") +
  scale_colour_manual(values = c(`FALSE` = "grey15", `TRUE` = "white")) +
  scale_x_discrete(labels = CONTRAST_LABEL_WRAP) +
  labs(x = NULL, y = NULL, title = sprintf("%s: leading and candidate genes in each fit", CAPILLARY),
       subtitle = sprintf("*/**/*** reported hit at Bonferroni p < 0.05/0.01/0.001; black outline = contrast in which the gene leads; grey = not tested; |effect| > %s drawn at the limit", LFC_CAP)) +
  theme(panel.grid = element_blank(), axis.text.y = element_text(face = "italic"))
save_svg(fig, "ec_capillary_leading_genes_heatmap.svg", width = 12, height = 6.5, dir = path_figures)
fig

make_rank
run_fgsea

gs_snapshot <- readRDS(file.path(PATHS$resources, "ec_gsea_gene_sets.rds"))
GS_LIST <- gs_snapshot$sets
rank_table <- function(universe_df) {
  mast_tab <- as.data.table(universe_df)[!is.na(gene) & nzchar(gene) & is.finite(avg_log2FC) & is.finite(p_val) & p_val >= 0 & p_val <= 1]
  mast_tab[, saturated := p_val <= .Machine$double.xmin]
  mast_tab[, rankstat := sign(avg_log2FC) * (-log10(pmax(p_val, .Machine$double.xmin)) + saturated * abs(avg_log2FC))]
  mast_tab
}
gsea_refit <- rbindlist(lapply(names(REFITS), function(refit_name) {
  set.seed(SEED)
  mast_tab <- rank_table(REFITS[[refit_name]]$uni)
  rbindlist(lapply(CONTRAST_ORDER, function(contrast_name) {
    rank_stats <- make_rank(mast_tab[contrast == contrast_name])
    rbindlist(lapply(names(GS_LIST), function(gsn)
      run_fgsea(rank_stats, GS_LIST[[gsn]])[, `:=`(fit = refit_name, contrast = contrast_name, geneset = gsn, universe_n = length(rank_stats))]))
  }))
}))
gsea_full <- readRDS(path_gsea_rds)[celltype == CAPILLARY][, fit := "full"]
gsea_fits <- rbindlist(list(gsea_full, gsea_refit), fill = TRUE)[, sig := is.finite(padj) & padj < FDR_CUT]
fwrite(gsea_refit[, .(fit, contrast, geneset, pathway, pval, padj, log2err, ES, NES, size, universe_n,
                      leadingEdge = vapply(leadingEdge, paste, "", collapse = ";"))],
       file.path(path_out, "ec_capillary_gsea_refits.csv"))

paired_gsea <- rbindlist(lapply(names(REFITS), function(refit_name)
  merge(gsea_fits[fit == "full", .(contrast, geneset, pathway, NES_full = NES, padj_full = padj, sig_full = sig)],
        gsea_fits[fit == refit_name, .(contrast, geneset, pathway, NES_refit = NES, padj_refit = padj, sig_refit = sig)],
        by = c("contrast", "geneset", "pathway"), all = TRUE)[, refit := refit_name]))
gsea_cmp <- paired_gsea[, .(sets_tested_both = sum(is.finite(NES_full) & is.finite(NES_refit)),
                            sig_full = sum(sig_full %in% TRUE), sig_refit = sum(sig_refit %in% TRUE),
                            shared = sum(sig_full %in% TRUE & sig_refit %in% TRUE),
                            sign_agree_pct = 100 * mean((sign(NES_refit) == sign(NES_full))[sig_full %in% TRUE & is.finite(NES_refit)]),
                            r_NES = cor(NES_full, NES_refit, use = "complete.obs"),
                            r_NES_sig = cor(NES_full[sig_full %in% TRUE], NES_refit[sig_full %in% TRUE], use = "complete.obs")),
                        by = .(refit, contrast)][, `:=`(shared_pct = 100 * shared / sig_full, jaccard = shared / (sig_full + sig_refit - shared))]
gsea_cmp <- gsea_cmp[order(match(refit, names(REFIT_LABEL)), match(contrast, CONTRAST_ORDER))]
fwrite(gsea_cmp, file.path(path_out, "ec_capillary_gsea_agreement.csv"))
knitr::kable(gsea_cmp[, .(refit = REFIT_LABEL[refit], contrast = CONTRAST_LABEL[contrast], sets_tested_both, sig_full, sig_refit, shared,
                          `shared (% of full)` = round(shared_pct, 1), jaccard = round(jaccard, 2), `NES sign agreement (%)` = round(sign_agree_pct, 1),
                          r_NES = round(r_NES, 3), r_NES_sig = round(r_NES_sig, 3))],
             caption = sprintf("Significant sets (FDR < %s, all three databases) of all nuclei (full) and of each refit per contrast; NES sign agreement: %% of the full-data significant sets tested in the refit; r_NES: Pearson correlation of NES over all sets tested in both fits, r_NES_sig over the full-data significant sets.", FDR_CUT))

gc <- function(refit_name, column) gsea_cmp[refit == refit_name][[column]]
flipped_sets <- paired_gsea[is.finite(NES_full) & is.finite(NES_refit) & sign(NES_full) != sign(NES_refit),
                            .(n = .N, n_sig = sum(sig_full %in% TRUE | sig_refit %in% TRUE)), by = refit]
gsea_low <- gsea_cmp[refit == "no_strata"][which.min(shared_pct)]

nes_df <- paired_gsea[is.finite(NES_full) & is.finite(NES_refit)][, status := factor(fcase(
  sig_full & sig_refit, "significant in both", sig_full, "all nuclei only", sig_refit, "refit only", default = "neither"),
  levels = c("neither", "all nuclei only", "refit only", "significant in both"))][order(status)]
nes_df[, `:=`(refit = factor(refit, levels = names(REFIT_LABEL)), contrast = factor(contrast, levels = CONTRAST_ORDER))]
nes_text <- gsea_cmp[, .(refit = factor(refit, levels = names(REFIT_LABEL)), contrast = factor(contrast, levels = CONTRAST_ORDER),
                         lab = sprintf("r = %.2f\nshared %d / %d", r_NES, shared, sig_full))]
fig <- ggplot(nes_df, aes(NES_full, NES_refit, colour = status)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60", linetype = "dashed", linewidth = .3) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = .3) + geom_vline(xintercept = 0, colour = "grey80", linewidth = .3) +
  geom_point(size = .5, alpha = .6) +
  geom_text(data = nes_text, aes(x = -Inf, y = Inf, label = lab), inherit.aes = FALSE, hjust = -.05, vjust = 1.2, size = 2.6) +
  facet_grid(refit ~ contrast, labeller = labeller(refit = REFIT_STRIP, contrast = CONTRAST_LABEL)) +
  scale_colour_manual(values = c(neither = "grey80", `all nuclei only` = "#D55E00", `refit only` = "#E69F00", `significant in both` = "#0072B2"), name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  labs(x = "NES, all nuclei (notebook 06)", y = "NES, refit",
       title = sprintf("%s: GSEA normalised enrichment scores of all nuclei and of each refit", CAPILLARY),
       subtitle = sprintf("Hallmark, Reactome and GO:BP sets tested in both fits; significant = FDR < %s (BH within database x contrast); shared = significant in both / significant with all nuclei", FDR_CUT)) +
  theme(legend.position = "bottom")
save_svg(fig, "ec_capillary_gsea_nes_agreement.svg", width = 14, height = 6.5, dir = path_figures)
fig

leading_sets <- gsea_fits[fit == "full" & sig][order(padj, -abs(NES))][, head(.SD, N_LEADING_SETS), by = contrast][, .(contrast, geneset, pathway)]
leading_sets_tab <- merge(leading_sets, gsea_fits[, .(fit, contrast, geneset, pathway, NES, padj, sig)], by = c("contrast", "geneset", "pathway"), all.x = TRUE)
fwrite(leading_sets_tab, file.path(path_out, "ec_capillary_leading_pathways.csv"))
dcast(leading_sets_tab[, .(contrast = factor(contrast, levels = CONTRAST_ORDER), geneset, pathway, fit = factor(fit, levels = names(FIT_LABEL)),
                           value = sprintf("%.2f (%s)%s", NES, formatC(padj, format = "e", digits = 1), ifelse(sig, "", " n.s.")))],
      contrast + geneset + pathway ~ fit, value.var = "value", fill = "not tested")[, contrast := CONTRAST_LABEL[as.character(contrast)]] |>
  setnames(names(FIT_LABEL), FIT_LABEL, skip_absent = TRUE) |>
  knitr::kable(caption = sprintf("The %d most significant sets per contrast of all nuclei (notebook 06): NES (FDR) per fit; n.s. = FDR >= %s.", N_LEADING_SETS, FDR_CUT))
lead_set_kept <- merge(leading_sets_tab[fit != "full"], leading_sets_tab[fit == "full", .(contrast, geneset, pathway, NES_full = NES)], by = c("contrast", "geneset", "pathway"))[
  , .(kept = sum(sig), same_sign = sum(sign(NES) == sign(NES_full)), max_nes_change = max(abs(NES - NES_full)),
      lost = paste(sprintf("%s (%s)", pathway[!sig], CONTRAST_LABEL[contrast[!sig]]), collapse = ", ")), by = fit]

session_footer()
