# 05b_EC_pseudobulk_sensitivity.R
# Code-only export of the analysis pipeline. Data are not included; see README.md.

params <- list(rds_in = "results/03_ec_annotation/Endothelial_subcluster_annotated_clean.rds", 
    mast_csv = "results/05_ec_de/ec_DE_all_MAST_annotated.csv", out_dir = "results/05b_ec_pseudobulk")


knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = TRUE,
                      fig.align = "center", dpi = 150)
options(width = 110)

suppressPackageStartupMessages({
  library(Seurat); library(edgeR); library(limma); library(DESeq2); library(tidyr); library(patchwork)
})
source("_common.R")
knitr::opts_chunk$set(cache = TRUE, autodep = TRUE, cache.lazy = FALSE, cache.extra = input_md5(params))
path_rds_in   <- resolve_path(params$rds_in)
path_mast_csv <- resolve_path(params$mast_csv)
path_out  <- resolve_path(params$out_dir)
path_figures <- file.path(path_out, "figures"); dir.create(path_figures, recursive = TRUE, showWarnings = FALSE)

FDR_CUT   <- 0.05
FDR_LOOSE <- 0.10
VOLCANO_CT <- "Capillary EC"
MIN_PB_NUCLEI <- 25

print_versions(c("ggplot2", "dplyr", "Seurat", "edgeR", "limma", "DESeq2", "tidyr", "patchwork"))

ec_obj <- readRDS(path_rds_in)
DefaultAssay(ec_obj) <- "RNA"
stopifnot(all(c("EC_subtype", "sample", "age", "condition", "batch", DESIGN_COVARIATES) %in% colnames(ec_obj@meta.data)))
SUBTYPE_ORDER <- levels(droplevels(factor(ec_obj$EC_subtype)))
ec_obj$EC_subtype <- factor(ec_obj$EC_subtype, levels = SUBTYPE_ORDER)
ec_obj$qc_flag <- unname(coalesce(QC_FLAGGED_ANIMALS[as.character(ec_obj$sample)], "pass"))
ec_obj$note    <- unname(coalesce(ANIMAL_NOTES[as.character(ec_obj$sample)], ""))
animal_meta <- ec_obj@meta.data |> distinct(sample, age, condition, sex, seq, batch, qc_flag, note) |>
  mutate(across(where(is.factor), as.character), group = make_group(age, condition))
if (anyDuplicated(animal_meta$sample)) stop("Expected one metadata row per animal")
mast_tab <- readr::read_csv(path_mast_csv, show_col_types = FALSE)
knitr::kable(animal_meta |> arrange(match(group, GROUP_ORDER), sample), row.names = FALSE, caption = "Animals entering the pseudobulk model (none excluded).")
sex_tab   <- as.data.frame.matrix(table(factor(as.character(animal_meta$group), levels = GROUP_ORDER), animal_meta$sex))
unbalanced_groups     <- rownames(sex_tab)[sex_tab$M != sex_tab$F]

ec_obj$pb_key <- paste(ec_obj$sample, ec_obj$EC_subtype, sep = "|")
pseudobulk_counts <- AggregateExpression(ec_obj, assays = "RNA", group.by = "pb_key")$RNA
pseudobulk_key_map <- ec_obj@meta.data |> count(sample, EC_subtype, pb_key, name = "n_nuclei") |>
  mutate(col = gsub("_", "-", pb_key), across(c(sample, EC_subtype), as.character))
stopifnot(setequal(colnames(pseudobulk_counts), pseudobulk_key_map$col))
colnames(pseudobulk_counts) <- pseudobulk_key_map$pb_key[match(colnames(pseudobulk_counts), pseudobulk_key_map$col)]
pseudobulk_samples <- pseudobulk_key_map |> select(sample, celltype = EC_subtype, pb_key, n_nuclei) |>
  mutate(lib_size = colSums(pseudobulk_counts)[pb_key], in_model = n_nuclei >= MIN_PB_NUCLEI) |> left_join(animal_meta, by = "sample")
readr::write_csv(pseudobulk_samples, file.path(path_out, "ec_pseudobulk_samples.csv"))
pseudobulk_samples |> transmute(sample, celltype, nuclei = paste0(n_nuclei, ifelse(in_model, "", " (below minimum)"))) |>
  pivot_wider(names_from = celltype, values_from = nuclei, values_fill = "0") |>
  knitr::kable(caption = sprintf("Nuclei per animal x subtype pseudobulk profile; profiles below %d nuclei do not enter that subtype's model.", MIN_PB_NUCLEI))

fit_pseudobulk <- function(subtype_name) {
  all_samples <- pseudobulk_samples |> filter(celltype == subtype_name)
  subtype_samples <- filter(all_samples, in_model)
  subtype_counts <- round(as.matrix(pseudobulk_counts[, subtype_samples$pb_key]))
  design_matrix <- model.matrix(DESIGN_FORMULA, data = subtype_samples); colnames(design_matrix) <- sub("^group", "", colnames(design_matrix))
  subtype_counts <- subtype_counts[filterByExpr(DGEList(subtype_counts), group = subtype_samples$group), ]
  dds <- DESeq(DESeqDataSetFromMatrix(subtype_counts, colData = data.frame(subtype_samples, row.names = subtype_samples$pb_key),
                                      design = design_matrix), quiet = TRUE)
  contrast_matrix <- make_contrast_matrix(design_matrix)[resultsNames(dds), ]
  contrast_res <- bind_rows(lapply(colnames(contrast_matrix), function(contrast_name)
    as.data.frame(results(dds, contrast = contrast_matrix[, contrast_name])) |>
      transmute(gene = rownames(dds), contrast = contrast_name, logFC = log2FoldChange, baseMean, stat, P.Value = pvalue, FDR = padj))) |>
    mutate(celltype = subtype_name, n_samples = nrow(subtype_samples), n_genes_tested = nrow(dds))
  dds_all <- estimateSizeFactors(DESeqDataSetFromMatrix(round(as.matrix(pseudobulk_counts[rownames(dds), all_samples$pb_key])),
                                                        colData = data.frame(all_samples, row.names = all_samples$pb_key), design = ~ 1))
  dispersionFunction(dds_all) <- dispersionFunction(dds)
  list(res = contrast_res, vst = assay(varianceStabilizingTransformation(dds_all, blind = FALSE)), samples = all_samples)
}
pseudobulk_fits <- lapply(setNames(SUBTYPE_ORDER, SUBTYPE_ORDER), fit_pseudobulk)
pseudobulk_de <- bind_rows(lapply(pseudobulk_fits, `[[`, "res")) |>
  mutate(celltype = factor(celltype, levels = SUBTYPE_ORDER), contrast = factor(contrast, levels = CONTRAST_ORDER),
         gene_flag = case_when(is_technical_gene(gene) ~ "technical", is_sex_gene(gene) ~ "sex", TRUE ~ ""),
         hit = coalesce(FDR < FDR_CUT, FALSE), hit_loose = coalesce(FDR < FDR_LOOSE, FALSE))
readr::write_csv(pseudobulk_de, file.path(path_out, "ec_pseudobulk_DE.csv"))
pseudobulk_summary <- pseudobulk_de |> group_by(celltype, contrast) |>
  summarise(n_samples = n_samples[1], genes_tested = n(), hits_FDR05 = sum(hit), hits_FDR10 = sum(hit_loose),
            hits_FDR05_technical_or_sex = sum(hit & gene_flag != ""), .groups = "drop")
pseudobulk_summary |> mutate(subtype = short_subtype(as.character(celltype)), contrast = CONTRAST_LABEL[as.character(contrast)]) |>
  select(subtype, contrast, everything(), -celltype) |>
  knitr::kable(caption = "Animal-level pseudobulk hits per subtype and contrast (BH within contrast).")

mds_df <- bind_rows(lapply(SUBTYPE_ORDER, function(subtype_name) {
  mds_fit <- plotMDS(pseudobulk_fits[[subtype_name]]$vst, plot = FALSE)
  pseudobulk_fits[[subtype_name]]$samples |> mutate(dim1 = mds_fit$x, dim2 = mds_fit$y, celltype = subtype_name)
})) |> mutate(celltype = factor(celltype, levels = SUBTYPE_ORDER))
fig <- ggplot(mds_df, aes(dim1, dim2, fill = group, shape = qc_flag == "pass", alpha = in_model)) +
  geom_point(size = 3, colour = "grey25", stroke = .4) +
  ggrepel::geom_text_repel(aes(label = ifelse(nzchar(note), paste0(sample, "*"), sample)), size = 2, colour = "grey30",
                           max.overlaps = Inf, min.segment.length = 0, segment.colour = "grey70", segment.size = .2, seed = 1) +
  facet_wrap(~ celltype, nrow = 1, scales = "free", labeller = as_labeller(short_subtype)) +
  scale_fill_manual(values = GROUP_FILL, labels = GROUP_LABEL, name = NULL) +
  scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 24), labels = c(`TRUE` = "qc pass", `FALSE` = "qc flagged"), name = NULL) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = .3), labels = c(`TRUE` = "in model", `FALSE` = sprintf("< %d nuclei, not in model", MIN_PB_NUCLEI)), name = NULL) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  scale_y_continuous(expand = expansion(mult = .2)) + scale_x_continuous(expand = expansion(mult = .1)) +
  labs(x = "MDS 1 (leading logFC)", y = "MDS 2", title = "Pseudobulk profiles per animal and subtype",
       subtitle = sprintf("* animal with a biological note; triangles = animals with a library-quality flag; faded = profile below %d nuclei", MIN_PB_NUCLEI)) +
  theme(legend.position = "bottom")
save_svg(fig, "ec_pseudobulk_mds.svg", width = 13, height = 3.4, dir = path_figures)
fig

ia_hits <- pseudobulk_de |> filter(contrast == "Interaction", FDR < FDR_LOOSE) |> arrange(FDR)
if (nrow(ia_hits)) ia_hits |>
  transmute(subtype = short_subtype(as.character(celltype)), gene, logFC = round(logFC, 2), stat = round(stat, 2),
            P = signif(P.Value, 2), FDR = signif(FDR, 2), gene_flag) |>
  knitr::kable(row.names = FALSE, caption = sprintf("Age x HFpEF interaction, (OH - OC) - (YH - YC): genes with FDR < %s in any subtype.", FDR_LOOSE))
interaction_min <- pseudobulk_de |> filter(contrast == "Interaction") |> slice_min(FDR, n = 1, with_ties = FALSE)

shared_genes <- mast_tab |> filter(contrast %in% CONTRAST_ORDER, is.finite(avg_log2FC)) |>
  transmute(celltype, contrast, gene, avg_log2FC, p_val, p_val_adj, reported, gene_flag, single_animal_driven,
            mast_stat = sign(avg_log2FC) * -log10(pmax(p_val, 1e-300))) |>
  inner_join(pseudobulk_de |> mutate(across(c(celltype, contrast), as.character)) |>
               select(celltype, contrast, gene, logFC, stat, P.Value, FDR, hit, hit_loose), by = c("celltype", "contrast", "gene")) |>
  mutate(celltype = factor(celltype, levels = SUBTYPE_ORDER), contrast = factor(contrast, levels = CONTRAST_ORDER),
         same_direction = sign(avg_log2FC) == sign(logFC))
concordance <- shared_genes |> group_by(celltype, contrast) |>
  summarise(genes_shared = n(), mast_hits = sum(reported),
            direction_agreement_mast_hits = round(mean(same_direction[reported]), 2),
            direction_agreement_single_animal = round(mean(same_direction[reported & single_animal_driven %in% TRUE]), 2),
            spearman_stat = round(cor(mast_stat, stat, method = "spearman", use = "complete.obs"), 2),
            pb_hits_FDR05 = sum(hit), pb_hits_FDR10 = sum(hit_loose),
            mast_hits_pb_FDR05 = sum(reported & hit), mast_hits_pb_FDR10 = sum(reported & hit_loose),
            mast_hits_pb_p05_same_dir = sum(reported & P.Value < 0.05 & same_direction), .groups = "drop")
readr::write_csv(dplyr::rename_with(concordance, ~ paste0(.x, "_prop"), starts_with("direction_agreement")), file.path(path_out, "ec_pseudobulk_concordance.csv"))
readr::write_csv(shared_genes, file.path(path_out, "ec_pseudobulk_vs_mast_genes.csv"))
concordance |> mutate(subtype = short_subtype(as.character(celltype)), contrast = CONTRAST_LABEL[as.character(contrast)],
               across(starts_with("direction_agreement"), \(x) ifelse(is.nan(x), "-", x))) |>
  select(subtype, contrast, everything(), -celltype) |>
  dplyr::rename_with(~ paste(.x, "(proportion)"), starts_with("direction_agreement")) |>
  knitr::kable(caption = "MAST (notebook 05, reported hits) vs animal-level pseudobulk, per subtype and contrast.")

CLASS_COL <- c("not a MAST hit" = "grey80", "MAST hit" = "#0072B2", "MAST hit, single-animal driven" = "#E69F00")
volcano_subtype_df <- shared_genes |> filter(celltype == VOLCANO_CT) |>
  mutate(class = factor(case_when(!reported ~ "not a MAST hit", single_animal_driven %in% TRUE ~ "MAST hit, single-animal driven",
                                  TRUE ~ "MAST hit"), levels = names(CLASS_COL)))
fig <- ggplot(volcano_subtype_df, aes(mast_stat, stat, colour = class)) +
  geom_hline(yintercept = 0, colour = "grey70") + geom_vline(xintercept = 0, colour = "grey70") +
  geom_point(size = .8, alpha = .8) +
  geom_point(data = volcano_subtype_df |> filter(hit), shape = 21, size = 2.4, colour = "black", fill = NA, stroke = .6) +
  facet_wrap(~ contrast, nrow = 1, labeller = as_labeller(CONTRAST_LABEL)) +
  scale_colour_manual(values = CLASS_COL, name = NULL) +
  labs(x = expression("MAST signed" ~ -log[10] ~ p ~ "(per nucleus)"), y = "pseudobulk DESeq2 Wald statistic (per animal)",
       title = sprintf("%s: per-nucleus MAST vs animal-level pseudobulk", VOLCANO_CT),
       subtitle = sprintf("Black rings: pseudobulk FDR < %s", FDR_CUT)) +
  theme(legend.position = "bottom")
save_svg(fig, "ec_pseudobulk_vs_mast_scatter.svg", width = 13, height = 4.2, dir = path_figures)
fig

ma_df <- pseudobulk_de |> filter(celltype == VOLCANO_CT) |>
  mutate(class = factor(case_when(hit ~ sprintf("FDR < %s", FDR_CUT), hit_loose ~ sprintf("FDR < %s", FDR_LOOSE), TRUE ~ "ns"),
                        levels = c("ns", sprintf("FDR < %s", FDR_LOOSE), sprintf("FDR < %s", FDR_CUT))))
fig <- ggplot(ma_df, aes(baseMean, logFC, colour = class)) +
  scale_x_log10() +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_point(size = .7, alpha = .8) +
  ggrepel::geom_text_repel(data = ma_df |> filter(hit_loose), aes(label = gene), size = 2.4, fontface = "italic", max.overlaps = 25, show.legend = FALSE) +
  facet_wrap(~ contrast, nrow = 1, labeller = as_labeller(CONTRAST_LABEL)) +
  scale_colour_manual(values = setNames(c("grey80", "#56B4E9", "#D55E00"), levels(ma_df$class)), name = NULL) +
  labs(x = "mean normalised count (DESeq2 baseMean, log scale)", y = expression(log[2]*"FC (pseudobulk)"),
       title = sprintf("%s: pseudobulk MA plots, four contrasts and the interaction", VOLCANO_CT)) +
  theme(legend.position = "bottom")
save_svg(fig, "ec_pseudobulk_MA_capillary.svg", width = 13, height = 4.2, dir = path_figures)
fig

hits_plot <- pseudobulk_summary |> pivot_longer(c(hits_FDR05, hits_FDR10), names_to = "threshold", values_to = "hits") |>
  mutate(threshold = factor(threshold, levels = c("hits_FDR10", "hits_FDR05"), labels = sprintf("FDR < %s", c(FDR_LOOSE, FDR_CUT))))
fig <- ggplot(hits_plot, aes(contrast, hits, fill = threshold)) +
  geom_col(position = position_dodge(width = .7), width = .65, colour = "grey30", linewidth = .2) +
  facet_wrap(~ celltype, nrow = 1, labeller = as_labeller(short_subtype)) +
  scale_fill_manual(values = c("#56B4E9", "#D55E00"), name = NULL) +
  scale_x_discrete(labels = CONTRAST_LABEL_WRAP) +
  labs(x = NULL, y = "genes", title = sprintf("Animal-level pseudobulk hits per contrast (profiles with >= %d nuclei)", MIN_PB_NUCLEI)) +
  theme(legend.position = "bottom", axis.text.x = element_text(size = 7))
save_svg(fig, "ec_pseudobulk_hits.svg", width = 12, height = 3.6, dir = path_figures)
fig

reported_in_pseudobulk   <- shared_genes |> filter(reported)
pseudobulk_hits  <- pseudobulk_de |> filter(hit) |> arrange(celltype, contrast, FDR) |>
  left_join(shared_genes |> select(celltype, contrast, gene, reported), by = c("celltype", "contrast", "gene"))
hits_by_contrast <- pseudobulk_hits |> group_by(contrast, .drop = FALSE) |> summarise(n = n(), mast = sum(reported %in% TRUE), .groups = "drop")
hit_text  <- paste(sprintf("%d in %s (%d of them MAST hits)", hits_by_contrast$n, CONTRAST_LABEL[as.character(hits_by_contrast$contrast)],
                           hits_by_contrast$mast), collapse = "; ")
agreement <- concordance$direction_agreement_mast_hits[concordance$mast_hits > 0]
mast_fdr_by_contrast <- reported_in_pseudobulk |> group_by(contrast, .drop = FALSE) |> summarise(n = n(), fdr = sum(hit), .groups = "drop") |>
  mutate(contrast = as.character(contrast))
mast_fdr_age <- sum(mast_fdr_by_contrast$fdr[mast_fdr_by_contrast$contrast %in% c("OC_vs_YC", "OH_vs_YH")])
mds_ext  <- mds_df |> group_by(celltype) |>
  summarise(ext1 = sample[which.max(abs(dim1 - median(dim1)))], ext2 = sample[which.max(abs(dim2 - median(dim2)))], .groups = "drop")
mds_extreme1_tab <- sort(table(mds_ext$ext1), decreasing = TRUE)
mds_extreme2_text <- paste(sprintf("%s (%s)", mds_ext$ext2, short_subtype(as.character(mds_ext$celltype))), collapse = ", ")
min_prof <- pseudobulk_samples |> filter(in_model) |> slice_min(n_nuclei, n = 1, with_ties = FALSE)
below_min <- filter(pseudobulk_samples, !in_model)
mds_extreme_groups  <- unique(as.character(animal_meta$group[animal_meta$sample %in% c(names(mds_extreme1_tab)[1], mds_ext$ext2)]))

session_footer()
