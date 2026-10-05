# 07_CCC_analysis.R
# Code-only export of the analysis pipeline. Data are not included; see README.md.

params <- list(rds_in = "data/objects/integrated_LV_snRNAseq_seurat_clean.rds", cleaning_tab = "results/02_subclustering/atlas_cleaning_by_cell_type.csv", 
    clean_rds = list(Endothelial = "results/03_ec_annotation/Endothelial_subcluster_annotated_clean.rds", 
        Cardiomyocyte = "results/04_cm_annotation/Cardiomyocyte_subcluster_annotated_clean.rds"), 
    lr_resource = "data/resources/mouseconsensus_LR.csv", out_dir = "results/07_ccc", 
    n_perms = 100L, group_col = "disease")


knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = TRUE, fig.align = "center", dpi = 150)
options(width = 110)
source("_common.R")
knitr::opts_chunk$set(cache = TRUE, autodep = TRUE, cache.lazy = FALSE, cache.extra = input_md5(params))

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(limma); library(CellChat)
  library(tidyr); library(patchwork); library(readr); library(tibble)
})

MIN_CELLS    <- 10
EXPR_PROP    <- 0.05
EXPR_PROP_WC <- 0.10
DET_ANIMAL   <- 0.05
MIN_ANIMALS  <- 3
FDR_CUT      <- 0.05
CC_TYPE      <- "truncatedMean"; CC_TRIM <- EXPR_PROP
CC_PVAL      <- 0.05
N_PERMS      <- params$n_perms
SEED         <- 1
TOP_N        <- 40
TOP_PATHWAYS <- 30
FOCUS_TYPES        <- c("Endothelial", "Cardiomyocyte")
CANONICAL_LIGANDS <- c("Nrg1", "Dll4", "Jag1", "Igf1", "Bmp6")
format_count <- function(x) formatC(x, big.mark = ",", format = "d")
set.seed(SEED)

path_out     <- resolve_path(params$out_dir)
path_tables <- file.path(path_out, "tables"); path_figures <- file.path(path_out, "figures")
invisible(lapply(c(path_tables, path_figures), dir.create, showWarnings = FALSE, recursive = TRUE))

print_versions(c("ggplot2", "dplyr", "Seurat", "Matrix", "limma", "CellChat", "tidyr", "patchwork", "readr", "tibble"))

atlas_obj <- readRDS(resolve_path(params$rds_in))
DefaultAssay(atlas_obj) <- "RNA"
fine_kept    <- lapply(params$clean_rds, function(p) colnames(readRDS(resolve_path(p))))
fine_removed <- sapply(names(fine_kept), function(ct) sum(atlas_obj$cell_type == ct & !colnames(atlas_obj) %in% fine_kept[[ct]]))
keep_cells   <- !atlas_obj$cell_type %in% names(fine_kept) | colnames(atlas_obj) %in% unlist(fine_kept)
atlas_obj    <- subset(atlas_obj, cells = colnames(atlas_obj)[keep_cells])
required <- c("cell_type", "sample", "individual", "age", "condition", DESIGN_COVARIATES, "batch", "exon_prop", params$group_col)
if (!all(required %in% colnames(atlas_obj@meta.data)) || anyNA(atlas_obj@meta.data[, required])) stop("Required metadata missing or incomplete")
if (!all(atlas_obj$cell_type %in% CELL_TYPE_ORDER)) stop("Update CELL_TYPE_ORDER in _common.R to match atlas cell types")
if (!all(atlas_obj@meta.data[[params$group_col]] %in% GROUP_ORDER)) stop("Unexpected condition codes")
CT_PRESENT <- CELL_TYPE_ORDER[CELL_TYPE_ORDER %in% atlas_obj$cell_type]
CC_TYPES   <- setdiff(CT_PRESENT, "Proliferating")
atlas_obj$cell_type  <- factor(as.character(atlas_obj$cell_type), levels = CT_PRESENT)
atlas_obj$disease    <- factor(as.character(atlas_obj@meta.data[[params$group_col]]), levels = GROUP_ORDER)
atlas_obj$sample     <- as.character(atlas_obj$sample)
atlas_obj$individual <- as.character(atlas_obj$individual)
atlas_obj$qc_flag     <- dplyr::coalesce(sub(":.*", "", unname(QC_FLAGGED_ANIMALS[atlas_obj$sample])), "pass")
atlas_obj$animal_note <- dplyr::coalesce(sub(" \\(.*", "", unname(ANIMAL_NOTES[atlas_obj$sample])), "")
animal_meta <- atlas_obj@meta.data |>
  distinct(sample, individual, batch, seq, age, condition, disease, sex, qc_flag, animal_note) |>
  arrange(disease, sample)
rownames(animal_meta) <- NULL
if (anyDuplicated(animal_meta$sample) || anyDuplicated(animal_meta$individual)) stop("Expected one sample and one metadata row per animal")
knitr::kable(animal_meta, caption = "Animal-level experimental metadata, library-quality flag and biological note.")
knitr::kable(data.frame(animal = c(names(QC_FLAGGED_ANIMALS), names(ANIMAL_NOTES)), kind = rep(c("qc_flag", "animal_note"), c(length(QC_FLAGGED_ANIMALS), length(ANIMAL_NOTES))),
                        meaning = c(QC_FLAGGED_ANIMALS, ANIMAL_NOTES), row.names = NULL),
             caption = "Meaning of the animal annotations (from _common.R). Nucleus counts in the flag text refer to the atlas after notebook-01 QC; the per-animal QC table below counts the nuclei analysed here. Flags annotate; they never exclude.")
cleaning_by_type <- readr::read_csv(resolve_path(params$cleaning_tab), show_col_types = FALSE)
knitr::kable(dplyr::rename(cleaning_by_type, `removed (%)` = pct_removed),
             caption = sprintf("Nuclei removed per compartment by the notebook-02 cleaning (input atlas = retained column; the Macrophage row includes the %s nuclei notebook 02 relabels Dendritic_cell).", format_count(sum(atlas_obj$cell_type == "Dendritic_cell"))))
path_lr_resource   <- resolve_path(params$lr_resource)

design <- atlas_obj@meta.data |>
  group_by(disease) |>
  summarise(nuclei = n(), animals = n_distinct(individual),
            females = n_distinct(individual[sex == "F"]), males = n_distinct(individual[sex == "M"]), .groups = "drop")
knitr::kable(design, caption = "Nuclei per group and their split into animals and sexes.", format.args = list(big.mark = ","))

animal_qc <- atlas_obj@meta.data |>
  group_by(sample) |>
  summarise(nuclei = n(), exon_prop_median = median(exon_prop), exon_prop_ec = median(exon_prop[cell_type == "Endothelial"]),
            exon_prop_cm = median(exon_prop[cell_type == "Cardiomyocyte"]), .groups = "drop") |>
  left_join(atlas_obj@meta.data |> dplyr::count(sample, cell_type) |>
              pivot_wider(names_from = cell_type, values_from = n, values_fill = 0), by = "sample") |>
  left_join(animal_meta |> dplyr::select(sample, disease, sex, seq, batch, qc_flag, animal_note), by = "sample") |>
  arrange(disease, sample) |>
  relocate(sample, disease, sex, seq, batch, nuclei, exon_prop_median, exon_prop_ec, exon_prop_cm, all_of(CT_PRESENT), qc_flag, animal_note)
readr::write_csv(dplyr::rename(animal_qc, median_exon_prop = exon_prop_median, ec_median_exon_prop = exon_prop_ec, cm_median_exon_prop = exon_prop_cm), file.path(path_tables, "animal_qc.csv"))
knitr::kable(dplyr::rename(animal_qc, `exon_prop_median (proportion)` = exon_prop_median, `exon_prop_ec (proportion)` = exon_prop_ec, `exon_prop_cm (proportion)` = exon_prop_cm), digits = 4, caption = "Per-animal QC: nuclei per cell type, median exonic read fraction over all nuclei, over the EC nuclei and over the CM nuclei (diagnostics, not model terms), library-quality flag and biological note. Every animal is retained.")

animal_meta_fit <- animal_meta |> mutate(group = make_group(age, condition), sex = factor(sex), seq = factor(seq))
DESIGN <- model.matrix(DESIGN_FORMULA, data = animal_meta_fit); colnames(DESIGN) <- sub("^group", "", colnames(DESIGN))
RESIDUAL_DF    <- nrow(DESIGN) - qr(DESIGN)$rank
sex_tab   <- table(animal_meta$disease, animal_meta$sex)
sex_text   <- sprintf("old controls are %d male / %d female (%s in the other groups)", sex_tab["OC", "M"], sex_tab["OC", "F"],
                     paste(unique(sprintf("%d / %d", sex_tab[rownames(sex_tab) != "OC", "M"], sex_tab[rownames(sex_tab) != "OC", "F"])), collapse = ", "))

path_cellchat <- file.path(path_out, "cellchat_list.rds")
CELLCHAT_DB   <- subsetDB(CellChatDB.mouse, search = c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"), key = "annotation")

run_one_cellchat <- function(atlas_subset_obj, group_id) {
  cells <- colnames(atlas_subset_obj)[atlas_subset_obj$disease == group_id & atlas_subset_obj$cell_type %in% CC_TYPES]
  expr_data   <- GetAssayData(atlas_subset_obj, assay = "RNA", layer = "data")[, cells, drop = FALSE]
  meta  <- data.frame(labels  = droplevels(atlas_subset_obj$cell_type[match(cells, colnames(atlas_subset_obj))]),
                      samples = factor(atlas_subset_obj$sample[match(cells, colnames(atlas_subset_obj))]), row.names = cells)
  cellchat_obj <- createCellChat(object = expr_data, meta = meta, group.by = "labels")
  cellchat_obj@DB <- CELLCHAT_DB
  cellchat_obj <- subsetData(cellchat_obj)
  cellchat_obj <- identifyOverExpressedGenes(cellchat_obj)
  cellchat_obj <- identifyOverExpressedInteractions(cellchat_obj)
  cellchat_obj <- computeCommunProb(cellchat_obj, type = CC_TYPE, trim = CC_TRIM, population.size = FALSE, nboot = N_PERMS, seed.use = SEED)
  cellchat_obj <- filterCommunication(cellchat_obj, min.cells = MIN_CELLS)
  cellchat_obj <- computeCommunProbPathway(cellchat_obj, thresh = CC_PVAL)
  aggregateNet(cellchat_obj, thresh = CC_PVAL)
}

future::plan("sequential"); options(future.globals.maxSize = 6 * 1024^3)
cellchat_list <- setNames(lapply(GROUP_ORDER, function(g) { cat("  CellChat:", g, "\n"); run_one_cellchat(atlas_obj, g) }), GROUP_ORDER)
saveRDS(cellchat_list, path_cellchat)

invisible(capture.output(suppressMessages({
  cellchat_list <- lapply(cellchat_list, liftCellChat, group.new = CC_TYPES)
  cellchat_merged  <- suppressWarnings(mergeCellChat(cellchat_list, add.names = GROUP_ORDER, cell.prefix = TRUE))
})))
group_index <- setNames(seq_along(GROUP_ORDER), GROUP_ORDER)
inferred_prob <- function(x) x@net$prob * (x@net$pval < CC_PVAL)

interaction_totals <- data.frame(
  group        = GROUP_ORDER,
  nuclei       = as.integer(table(atlas_obj$disease[atlas_obj$cell_type %in% CC_TYPES])[GROUP_ORDER]),
  n_inferred   = sapply(cellchat_list, function(x) sum(x@net$count)),
  n_lr_pairs   = sapply(cellchat_list, function(x) sum(apply(inferred_prob(x) > 0, 3, any))),
  n_pathways   = sapply(cellchat_list, function(x) length(x@netP$pathways)),
  EC_to_CM     = sapply(cellchat_list, function(x) sum(inferred_prob(x)["Endothelial", "Cardiomyocyte", ] > 0)),
  CM_to_EC     = sapply(cellchat_list, function(x) sum(inferred_prob(x)["Cardiomyocyte", "Endothelial", ] > 0)),
  total_weight = sapply(cellchat_list, function(x) sum(x@net$weight)), row.names = NULL)
knitr::kable(interaction_totals, digits = 3, format.args = list(big.mark = ","),
             caption = "Nuclei, inferred interactions (sender-receiver-LR-pair combinations), distinct LR pairs, active pathways, inferred EC->CM and CM->EC LR pairs, and summed interaction weight per group.")

cellchat_genes  <- intersect(extractGene(CELLCHAT_DB), rownames(atlas_obj))
detection_focus <- sapply(FOCUS_TYPES, function(cell_type_name) Matrix::rowMeans(GetAssayData(atlas_obj, assay = "RNA", layer = "data")[cellchat_genes, atlas_obj$cell_type == cell_type_name] > 0))
n_detected     <- sapply(c(triMean = 0.25, truncatedMean = CC_TRIM), function(th) colSums(detection_focus > th))

ec_cm_lr <- bind_rows(lapply(GROUP_ORDER, function(g) {
  prob_mat <- inferred_prob(cellchat_list[[g]])
  bind_rows(lapply(list(FOCUS_TYPES, rev(FOCUS_TYPES)), function(cell_pair) data.frame(group = g, source = cell_pair[1], target = cell_pair[2],
                                                                     interaction_name = dimnames(prob_mat)[[3]][prob_mat[cell_pair[1], cell_pair[2], ] > 0])))
}))
ecm_dirs <- ec_cm_lr |> group_by(source, target) |> summarise(n_lr = n_distinct(interaction_name), .groups = "drop")
BUBBLE_H <- 1.6 + 0.2 * max(ecm_dirs$n_lr); BUBBLE_W <- 1.5 + 3.5 * nrow(ecm_dirs)

cellchat_lr_db <- CellChatDB.mouse$interaction
path_flow <- bind_rows(lapply(GROUP_ORDER, function(g) {
  pathway_prob    <- cellchat_list[[g]]@netP$prob
  n_interactions_of <- apply(inferred_prob(cellchat_list[[g]]) > 0, 3, sum)
  n_interactions_of <- tapply(n_interactions_of, cellchat_lr_db$pathway_name[match(names(n_interactions_of), rownames(cellchat_lr_db))], sum)
  data.frame(group = g, pathway = dimnames(pathway_prob)[[3]], flow = apply(pathway_prob, 3, sum), n_interactions = as.integer(n_interactions_of[dimnames(pathway_prob)[[3]]]))
})) |> group_by(group) |> mutate(rel_flow = flow / sum(flow)) |> ungroup()
readr::write_csv(dplyr::rename(path_flow, rel_flow_prop = rel_flow), file.path(path_tables, "cellchat_pathway_flow.csv"))
n_union <- sapply(CONTRASTS, function(contrast_def) n_distinct(path_flow$pathway[path_flow$group %in% c(contrast_def$i1, contrast_def$i2)]))

bar_panel <- function(m) compareInteractions(cellchat_merged, show.legend = FALSE, group = seq_along(GROUP_ORDER), measure = m, color.use = unname(GROUP_FILL[GROUP_ORDER])) +
  scale_x_discrete(labels = sub(" ", "\n", GROUP_LABEL)) + theme_house() + theme(legend.position = "none")
fig <- bar_panel("count") + bar_panel("weight") +
  plot_annotation(tag_levels = "A",
                  caption = sprintf("Pooled-group inferences; nuclei per group: %s",
                                    paste(interaction_totals$group, format_count(interaction_totals$nuclei), sep = " = ", collapse = ", ")))
save_svg(fig, "ccc_total_interactions.svg", width = 8, height = 4, dir = path_figures)
fig

cellchat_merged_short <- cellchat_merged
cellchat_merged_short@net <- lapply(cellchat_merged@net, function(n) { dimnames(n$weight) <- rep(list(unname(CELL_TYPE_SHORT[rownames(n$weight)])), 2); n })
circle_file <- file.path(path_figures, "ccc_diff_interaction_circles.svg")
svg(circle_file, width = 10, height = 10)
par(mfrow = c(2, 2), xpd = TRUE, oma = c(3.5, 0, 0, 0), mar = c(1, 1, 2, 1))
for (contrast_def in CONTRASTS)

  netVisual_diffInteraction(cellchat_merged_short, weight.scale = TRUE, measure = "weight", comparison = c(group_index[[contrast_def$i2]], group_index[[contrast_def$i1]]),
                            color.use = unname(CELL_TYPE_FILL[CC_TYPES]), vertex.label.cex = 0.8, margin = 0.15,
                            title.name = sprintf("%s  (red = up in %s)", CONTRAST_LABEL[[contrast_def$name]], contrast_def$i1))
mtext(c("Edge colour: red = higher summed communication probability in the first group, blue = in the second; edge width scaled within each panel",
        sapply(split(sprintf("%s = %s", CELL_TYPE_SHORT[CC_TYPES], CELL_TYPE_LABEL[CC_TYPES]), ceiling(seq_along(CC_TYPES) / 6)), paste, collapse = ", ")),
      side = 1, line = c(-0.4, 0.6, 1.6), outer = TRUE, cex = 0.7, col = "grey30")
invisible(dev.off())
knitr::include_graphics(circle_file, rel_path = FALSE)

complex_db  <- CellChatDB.mouse$complex
lr_pair_key <- function(x) { m <- match(x, rownames(complex_db)); x[!is.na(m)] <- apply(complex_db[m[!is.na(m)], , drop = FALSE], 1, function(r) paste(sort(r[nzchar(r)]), collapse = "_")); x }
prob_long <- function(g) as.data.frame.table(inferred_prob(cellchat_list[[g]]), responseName = "prob", stringsAsFactors = FALSE) |>
  setNames(c("source", "target", "interaction_name", "prob"))

edge_decomposition <- bind_rows(lapply(CONTRASTS, function(contrast_def)
  full_join(prob_long(contrast_def$i1) |> rename(prob_group1 = prob), prob_long(contrast_def$i2) |> rename(prob_group2 = prob), by = c("source", "target", "interaction_name")) |>
    mutate(across(c(prob_group1, prob_group2), ~ replace_na(.x, 0)), contrast = contrast_def$name, group1 = contrast_def$i1, group2 = contrast_def$i2) |>
    filter(prob_group1 > 0 | prob_group2 > 0))) |>
  mutate(delta = prob_group1 - prob_group2, pathway = cellchat_lr_db$pathway_name[match(interaction_name, rownames(cellchat_lr_db))],
         lk = lr_pair_key(cellchat_lr_db$ligand[match(interaction_name, rownames(cellchat_lr_db))]), rk = lr_pair_key(cellchat_lr_db$receptor[match(interaction_name, rownames(cellchat_lr_db))])) |>
  group_by(contrast, source, target) |>
  mutate(edge_delta_total = sum(delta), share_of_edge = if_else(edge_delta_total == 0, NA_real_, delta / edge_delta_total)) |> ungroup() |>
  dplyr::select(contrast, group1, group2, source, target, pathway, interaction_name, prob_group1, prob_group2, delta, edge_delta_total, share_of_edge, lk, rk)

max_disc <- max(sapply(CONTRASTS, function(contrast_def) {
  weight_delta   <- cellchat_list[[contrast_def$i1]]@net$weight - cellchat_list[[contrast_def$i2]]@net$weight
  edge_totals <- edge_decomposition |> filter(contrast == contrast_def$name) |> distinct(source, target, edge_delta_total)
  m   <- matrix(0, nrow(weight_delta), ncol(weight_delta), dimnames = dimnames(weight_delta)); m[cbind(edge_totals$source, edge_totals$target)] <- edge_totals$edge_delta_total
  max(abs(weight_delta - m)) }))
cat(sprintf("Decomposition check (circle plot measure = weight): max |plotted edge difference - sum of pair deltas| = %.1e over %d edges\n",
            max_disc, n_distinct(edge_decomposition[, c("contrast", "source", "target")])))

top_edges <- edge_decomposition |> distinct(contrast, source, target, edge_delta_total) |> group_by(contrast) |> slice_max(abs(edge_delta_total), n = 6, with_ties = FALSE) |> ungroup()
bar_df <- edge_decomposition |> semi_join(top_edges, by = c("contrast", "source", "target")) |>
  group_by(contrast, source, target) |> slice_max(abs(delta), n = 8, with_ties = FALSE) |> ungroup() |>
  mutate(edge  = sprintf("%s \u2192 %s (total \u0394 = %+.3f)", CELL_TYPE_SHORT[source], CELL_TYPE_SHORT[target], edge_delta_total),
         row   = paste(contrast, source, target, interaction_name, sep = "||"))
bar_df$row  <- factor(bar_df$row, levels = bar_df$row[order(bar_df$delta)])
bar_df$edge <- factor(bar_df$edge, levels = unique(bar_df$edge[order(-abs(bar_df$edge_delta_total))]))
decomposition_plot <- function(contrast_def) {
  weight_delta <- bar_df |> filter(contrast == contrast_def$name)
  ggplot(weight_delta, aes(delta, row, fill = delta > 0)) +
    geom_col(width = .7) +
    geom_text(aes(label = pathway, hjust = ifelse(delta > 0, -0.1, 1.1)), size = 2.2, colour = "grey35") +
    geom_vline(xintercept = 0, colour = "grey50", linewidth = .3) +
    facet_wrap(~ edge, nrow = 2, scales = "free") +
    scale_fill_manual(values = c(`TRUE` = "#b2182b", `FALSE` = "#2166ac"), guide = "none") +
    scale_y_discrete(labels = function(x) sub("^.*\\|\\|", "", x)) +
    scale_x_continuous(expand = expansion(mult = c(0.7, 0.7))) +
    labs(x = sprintf("\u0394 communication probability (%s \u2212 %s)", contrast_def$i1, contrast_def$i2), y = NULL,
         title = sprintf("%s: red = higher in %s, blue = higher in %s", CONTRAST_LABEL[[contrast_def$name]], contrast_def$i1, contrast_def$i2),
         subtitle = "six edges with the largest |total \u0394|; up to eight pairs each; label = pathway") +
    theme(axis.text.y = element_text(size = 7), strip.text = element_text(size = 8))
}
fig <- wrap_plots(lapply(CONTRASTS, decomposition_plot), ncol = 1)
save_svg(fig, "ccc_diff_edge_decomposition.svg", width = 13, height = 22, dir = path_figures)
fig

signaling_role <- bind_rows(lapply(GROUP_ORDER, function(g) {
  weight_mat <- cellchat_list[[g]]@net$weight
  data.frame(group = g, cell_type = rownames(weight_mat), outgoing = rowSums(weight_mat), incoming = colSums(weight_mat))
})) |>
  mutate(group = factor(group, levels = GROUP_ORDER), cell_type = factor(cell_type, levels = CC_TYPES)) |>
  pivot_longer(c(outgoing, incoming), names_to = "direction", values_to = "strength")

fig <- ggplot(signaling_role, aes(group, strength, fill = group)) +
  geom_col(width = .7, colour = "grey20", linewidth = .3) +
  facet_grid(direction ~ cell_type, scales = "free_y", labeller = labeller(cell_type = CELL_TYPE_LABEL)) +
  scale_fill_manual(values = GROUP_FILL, labels = GROUP_LABEL, name = NULL) +
  labs(x = NULL, y = "summed communication probability",
       title = "Outgoing (sender) and incoming (receiver) signalling strength per cell type") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom")
save_svg(fig, "ccc_signaling_roles.svg", width = 13, height = 5, dir = path_figures)
fig

top_pathways <- function(i1, i2, n = TOP_PATHWAYS) {
  path_flow |> filter(group %in% c(i1, i2)) |> dplyr::select(pathway, group, rel_flow) |>
    pivot_wider(names_from = group, values_from = rel_flow, values_fill = 0) |>
    mutate(delta = .data[[i1]] - .data[[i2]]) |> slice_max(abs(delta), n = n, with_ties = FALSE)
}
rank_plots <- lapply(CONTRASTS, function(contrast_def)
  rankNet(cellchat_merged, mode = "comparison", measure = "weight", comparison = c(group_index[[contrast_def$i2]], group_index[[contrast_def$i1]]),
          signaling = top_pathways(contrast_def$i1, contrast_def$i2)$pathway,
          stacked = TRUE, do.stat = FALSE, color.use = GROUP_FILL[c(contrast_def$i2, contrast_def$i1)]) +
    labs(title = CONTRAST_LABEL[[contrast_def$name]]) + theme(plot.title = element_text(size = 9, face = "bold")))
fig <- wrap_plots(rank_plots, nrow = 2)
save_svg(fig, "ccc_information_flow.svg", width = 13, height = 9, dir = path_figures)
fig

bubble_plot <- function(source_type, target_type)
  netVisual_bubble(cellchat_merged, sources.use = source_type, targets.use = target_type, comparison = seq_along(GROUP_ORDER), thresh = CC_PVAL,
                   angle.x = 0, remove.isolate = TRUE, color.text = unname(GROUP_FILL[GROUP_ORDER])) +
    scale_x_discrete(labels = function(x) sub("^.*\\((.*)\\)$", "\\1", x)) +
    labs(title = sprintf("%s \u2192 %s", CELL_TYPE_LABEL[[source_type]], CELL_TYPE_LABEL[[target_type]]))
fig <- wrap_plots(Map(bubble_plot, ecm_dirs$source, ecm_dirs$target), nrow = 1)
save_svg(fig, "ccc_ecm_bubble.svg", width = BUBBLE_W, height = BUBBLE_H, dir = path_figures)
fig

path_lr_scores <- file.path(path_out, "lr_scores_by_sample.csv")
ligand_receptor <- readr::read_csv(path_lr_resource, show_col_types = FALSE) |> filter(!is.na(ligand), !is.na(receptor)) |> distinct(ligand, receptor)
expr_data <- GetAssayData(atlas_obj, assay = "RNA", layer = "data")
has_subunits <- function(complex_ids) vapply(strsplit(complex_ids, "_"), function(x) all(x %in% rownames(expr_data)), logical(1))
present <- has_subunits(ligand_receptor$ligand) & has_subunits(ligand_receptor$receptor)
readr::write_csv(ligand_receptor[!present, ], file.path(path_tables, "lr_missing_subunits.csv"))
ligand_receptor <- ligand_receptor[present, ]

in_scope   <- which(atlas_obj$cell_type %in% CC_TYPES)
idx_animal <- split(in_scope, paste(atlas_obj$sample,  atlas_obj$cell_type, sep = "||")[in_scope])
idx_animal <- idx_animal[lengths(idx_animal) >= MIN_CELLS]
idx_group  <- split(in_scope, paste(atlas_obj$disease, atlas_obj$cell_type, sep = "||")[in_scope])
mean_expr     <- sapply(idx_animal, function(i) Matrix::rowMeans(expr_data[, i, drop = FALSE]))
detect_frac     <- sapply(idx_animal, function(i) Matrix::rowMeans(expr_data[, i, drop = FALSE] > 0))
detect_frac_group <- sapply(idx_group,  function(i) Matrix::rowMeans(expr_data[, i, drop = FALSE] > 0))

reduce_complex <- function(value_mat, complex_ids) {
  m <- do.call(rbind, lapply(strsplit(complex_ids, "_"), function(su) apply(value_mat[su, , drop = FALSE], 2, min)))
  dimnames(m) <- list(complex_ids, colnames(value_mat)); m
}
ligands <- unique(ligand_receptor$ligand); receptors <- unique(ligand_receptor$receptor)
ligand_mean <- reduce_complex(mean_expr, ligands); ligand_detect <- reduce_complex(detect_frac, ligands); ligand_detect_group <- reduce_complex(detect_frac_group, ligands)
receptor_mean <- reduce_complex(mean_expr, receptors); receptor_detect <- reduce_complex(detect_frac, receptors); receptor_detect_group <- reduce_complex(detect_frac_group, receptors)

score_direction <- function(source_types, target_types) {
  group_detect_max <- do.call(pmax, lapply(GROUP_ORDER, function(g)
    pmin(ligand_detect_group[ligand_receptor$ligand, paste(g, source_types, sep = "||")], receptor_detect_group[ligand_receptor$receptor, paste(g, target_types, sep = "||")])))
  group_eligible  <- group_detect_max >= EXPR_PROP
  animals <- animal_meta$sample[paste(animal_meta$sample, source_types, sep = "||") %in% colnames(mean_expr) & paste(animal_meta$sample, target_types, sep = "||") %in% colnames(mean_expr)]
  if (!any(group_eligible) || !length(animals)) return(NULL)
  ligand_ok <- ligand_receptor$ligand[group_eligible]; receptor_ok <- ligand_receptor$receptor[group_eligible]
  source_keys  <- paste(animals, source_types, sep = "||"); target_keys <- paste(animals, target_types, sep = "||")
  data.frame(sample = rep(animals, each = length(ligand_ok)), source = source_types, target = target_types, ligand_complex = ligand_ok, receptor_complex = receptor_ok,
             group_detection = rep(group_detect_max[group_eligible], times = length(animals)),
             expr_prod = as.vector(ligand_mean[ligand_ok, source_keys, drop = FALSE] * receptor_mean[receptor_ok, target_keys, drop = FALSE]),
             detected  = as.vector(ligand_detect[ligand_ok, source_keys, drop = FALSE] >= DET_ANIMAL & receptor_detect[receptor_ok, target_keys, drop = FALSE] >= DET_ANIMAL))
}
cell_type_pairs  <- expand.grid(source = CC_TYPES, target = CC_TYPES, stringsAsFactors = FALSE)
lr_scores <- bind_rows(Map(score_direction, cell_type_pairs$source, cell_type_pairs$target)) |>
  left_join(animal_meta |> dplyr::select(sample, individual, disease, age, condition, sex, qc_flag, animal_note), by = "sample")
readr::write_csv(lr_scores, path_lr_scores)
lr_scores <- lr_scores |>
  mutate(disease = factor(disease, levels = GROUP_ORDER), ccc_pair = paste0(source, "->", target),
         lr = paste0(ligand_complex, " \u2192 ", receptor_complex))
KEY <- c("source", "target", "ligand_complex", "receptor_complex")
ecm_elig <- lr_scores |> filter(source %in% FOCUS_TYPES, target %in% FOCUS_TYPES, source != target) |> distinct(across(all_of(KEY)), group_detection)
cat(sprintf("Expression scores: %s rows, %s eligible pair-directions, %d animals\n",
            format_count(nrow(lr_scores)), format_count(n_distinct(lr_scores[, KEY])), n_distinct(lr_scores$sample)))

lr_resource_tab <- readr::read_csv(path_lr_resource, show_col_types = FALSE)
axis_lr  <- lr_resource_tab |> filter(ligand %in% CANONICAL_LIGANDS) |> distinct(ligand, receptor)
axis_genes <- intersect(unique(unlist(strsplit(c(axis_lr$ligand, axis_lr$receptor), "_"))), rownames(atlas_obj))
axis_cells <- split(which(atlas_obj$cell_type %in% FOCUS_TYPES), paste(atlas_obj$disease, atlas_obj$cell_type, sep = "||")[atlas_obj$cell_type %in% FOCUS_TYPES])
axis_detection   <- sapply(axis_cells, function(i) Matrix::rowMeans(GetAssayData(atlas_obj, assay = "RNA", layer = "data")[axis_genes, i, drop = FALSE] > 0))
complex_rate  <- function(complex_ids, key) vapply(strsplit(complex_ids, "_"), function(su) if (all(su %in% axis_genes)) min(axis_detection[su, key]) else 0, numeric(1))
canonical_axes <- bind_rows(lapply(list(FOCUS_TYPES, rev(FOCUS_TYPES)), function(cell_pair) {
  ligand_rate <- sapply(GROUP_ORDER, function(g) complex_rate(axis_lr$ligand,   paste(g, cell_pair[1], sep = "||")))
  receptor_rate <- sapply(GROUP_ORDER, function(g) complex_rate(axis_lr$receptor, paste(g, cell_pair[2], sep = "||")))
  axis_lr |> mutate(direction = sprintf("%s\u2192%s", CELL_TYPE_SHORT[cell_pair[1]], CELL_TYPE_SHORT[cell_pair[2]]),
                  ligand_det = apply(matrix(ligand_rate, nrow(axis_lr)), 1, max), receptor_det = apply(matrix(receptor_rate, nrow(axis_lr)), 1, max),
                  eligible = apply(matrix(pmin(ligand_rate, receptor_rate), nrow(axis_lr)), 1, max) >= EXPR_PROP)
})) |>
  group_by(ligand = factor(ligand, levels = CANONICAL_LIGANDS), direction) |>
  summarise(pairs_in_resource = n(), max_ligand_detection = round(max(ligand_det), 3),
            max_receptor_detection = round(max(receptor_det), 3), eligible = sum(eligible), .groups = "drop")
nrg_ligands <- sort(unique(grep("^Nrg", lr_resource_tab$ligand, value = TRUE)))
knitr::kable(dplyr::rename(canonical_axes, `max_ligand_detection (proportion)` = max_ligand_detection, `max_receptor_detection (proportion)` = max_receptor_detection), caption = sprintf("Literature EC-CM axes in the LIANA resource: pairs per ligand and direction, the highest group-level detection of the ligand in the sender and of any of its receptors in the receiver (all four groups), and pair-directions reaching the %s %% rule. Not in the resource: %s.",
                                     EXPR_PROP * 100, paste(setdiff(CANONICAL_LIGANDS, lr_resource_tab$ligand), collapse = ", ")))

per_sample <- lr_scores |>
  group_by(sample, source, target) |>
  summarise(strength = sum(expr_prod), mean_score = mean(expr_prod), n_pairs = n(), .groups = "drop")
n_pairs_df <- per_sample |> distinct(source, target, n_pairs)
score_by_group <- expand_grid(sample = animal_meta$sample, source = CC_TYPES, target = CC_TYPES) |>
  left_join(per_sample, by = c("sample", "source", "target")) |>
  left_join(animal_meta |> dplyr::select(sample, disease), by = "sample") |>
  group_by(disease, source, target) |>
  summarise(n_animals = sum(is.finite(mean_score)),
            mean_score = if (n_animals >= MIN_ANIMALS) mean(mean_score, na.rm = TRUE) else NA_real_, .groups = "drop") |>
  left_join(n_pairs_df, by = c("source", "target")) |>
  mutate(source = factor(source, levels = CC_TYPES), target = factor(target, levels = CC_TYPES))
readr::write_csv(score_by_group, file.path(path_tables, "lr_group_overview.csv"))

ecm_total <- per_sample |> filter(source %in% FOCUS_TYPES, target %in% FOCUS_TYPES, source != target) |>
  group_by(sample) |> summarise(ecm_score = sum(strength), .groups = "drop") |>
  left_join(animal_qc |> dplyr::select(sample, disease, exon_prop_ec, exon_prop_cm, qc_flag), by = "sample") |>
  mutate(score_rank = rank(-ecm_score), exon_prop_ec_rank = rank(-exon_prop_ec), exon_prop_cm_rank = rank(-exon_prop_cm))
cor_exon_ec <- cor.test(ecm_total$ecm_score, ecm_total$exon_prop_ec, method = "spearman")
cor_exon_cm <- cor.test(ecm_total$ecm_score, ecm_total$exon_prop_cm, method = "spearman")
knitr::kable(ecm_total |> arrange(disease, sample) |> dplyr::rename(`exon_prop_ec (proportion)` = exon_prop_ec, `exon_prop_cm (proportion)` = exon_prop_cm), digits = c(NA, 2, NA, 4, 4, NA, 0, 0, 0),
             caption = sprintf("Summed EC<->CM score per animal vs. median exonic read fraction of its EC nuclei and of its CM nuclei (ranks 1 = highest; Spearman rho = %.2f, p = %.3f with the EC fraction and rho = %.2f, p = %.3f with the CM fraction, n = %d).",
                               cor_exon_ec$estimate, cor_exon_ec$p.value, cor_exon_cm$estimate, cor_exon_cm$p.value, nrow(ecm_total)))

tile_axes <- list(scale_y_discrete(limits = rev(CC_TYPES), labels = CELL_TYPE_SHORT), labs(y = "sender (source)"),
                  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7), axis.text.y = element_text(size = 7), panel.grid = element_blank()))
p_pairs <- ggplot(n_pairs_df, aes(factor(target, levels = CC_TYPES), factor(source, levels = CC_TYPES))) +
  geom_tile(fill = "grey95", colour = "white", linewidth = .4) + geom_text(aes(label = n_pairs), size = 2.2, colour = "grey30") +
  scale_x_discrete(limits = CC_TYPES, labels = CELL_TYPE_SHORT) + tile_axes + labs(x = "receiver (target)", title = "eligible LR pairs")
p_scores <- ggplot(score_by_group, aes(target, source, fill = mean_score)) +
  geom_tile(colour = "white", linewidth = .4) +
  facet_wrap(~ disease, nrow = 1, labeller = as_labeller(GROUP_LABEL)) +
  scale_fill_gradient(low = "#f7fbff", high = "#08306b", name = "mean score\nper eligible\nLR pair", na.value = "grey85") +
  scale_x_discrete(limits = CC_TYPES, labels = CELL_TYPE_SHORT) + tile_axes +
  labs(x = "receiver (target)", title = "Mean expression-product score per eligible LR pair, by group",
       subtitle = sprintf("grey = fewer than %d animals with both populations", MIN_ANIMALS))
fig <- p_pairs + p_scores + plot_layout(widths = c(1, 4))
save_svg(fig, "ccc_expression_score_overview.svg", width = 13, height = 4, dir = path_figures)
fig

ecm_scores <- lr_scores |> filter(source %in% FOCUS_TYPES, target %in% FOCUS_TYPES, source != target) |>
  mutate(sample = factor(sample, levels = animal_meta_fit$sample))
score_wide <- ecm_scores |> dplyr::select(all_of(KEY), sample, expr_prod) |>
  pivot_wider(names_from = sample, values_from = expr_prod, names_expand = TRUE)
detected_wide <- ecm_scores |> dplyr::select(all_of(KEY), sample, detected) |>
  pivot_wider(names_from = sample, values_from = detected, names_expand = TRUE, values_fill = FALSE)
score_mat <- as.matrix(score_wide[, animal_meta_fit$sample]); detected_mat <- as.matrix(detected_wide[, animal_meta_fit$sample])
log_score_mat <- log1p(score_mat)
contrast_groups <- c(setNames(lapply(CONTRASTS, function(x) c(x$i1, x$i2)), sapply(CONTRASTS, `[[`, "name")), list(Interaction = GROUP_ORDER))

fit_groups  <- lmFit(log_score_mat, DESIGN)
fit   <- eBayes(contrasts.fit(fit_groups, make_contrast_matrix(DESIGN)), robust = TRUE)
n_detected_group  <- sapply(GROUP_ORDER, function(g) rowSums(detected_mat[, animal_meta_fit$sample[animal_meta_fit$group == g], drop = FALSE]))
mean_score_group <- sapply(GROUP_ORDER, function(g) rowMeans(score_mat[, animal_meta_fit$sample[animal_meta_fit$group == g], drop = FALSE], na.rm = TRUE))
colnames(n_detected_group) <- paste0("ndet_", GROUP_ORDER); colnames(mean_score_group) <- paste0("mean_", GROUP_ORDER)
diff_all <- bind_rows(lapply(colnames(fit$coefficients), function(contrast_name) data.frame(
  score_wide[KEY], contrast = contrast_name, logFC = fit$coefficients[, contrast_name], t_mod = fit$t[, contrast_name], p_value = fit$p.value[, contrast_name], n_detected_group, mean_score_group,
  eligible = apply(n_detected_group[, paste0("ndet_", contrast_groups[[contrast_name]]), drop = FALSE], 1, max) >= MIN_ANIMALS))) |>
  filter(eligible) |> dplyr::select(-eligible) |>
  group_by(contrast) |> mutate(FDR = p.adjust(p_value, "BH")) |> ungroup() |>
  mutate(FDR_all_contrasts = p.adjust(p_value, "BH"), passes = FDR < FDR_CUT, contrast = factor(contrast, levels = CONTRAST_ORDER),
         ccc_pair = paste0(source, "->", target))
readr::write_csv(diff_all, file.path(path_tables, "lr_diff_all.csv"))

ndet_str <- function(YC, OC, YH, OH) sprintf("%d/%d/%d/%d", YC, OC, YH, OH)
passing_hits <- diff_all |> filter(passes) |>
  transmute(comparison = CONTRAST_LABEL[as.character(contrast)], direction = ccc_pair,
            ligand = ligand_complex, receptor = receptor_complex,
            logFC = round(logFC, 3), t_mod = round(t_mod, 2), FDR = signif(FDR, 3),
            `detected YC/OC/YH/OH` = ndet_str(ndet_YC, ndet_OC, ndet_YH, ndet_OH)) |>
  arrange(comparison, FDR)
knitr::kable(passing_hits, caption = sprintf("EC-CM pair-directions with BH FDR < %s within contrast (%d-animal limma model)%s.", FDR_CUT, nrow(animal_meta),
                                     if (nrow(passing_hits) == 0) ": none" else ""))

sort_key <- function(x) vapply(strsplit(x, "_"), function(g) paste(sort(g), collapse = "_"), character(1))
edge_decomposition <- edge_decomposition |>
  left_join(diff_all |> transmute(contrast = as.character(contrast), source, target, lk = sort_key(ligand_complex), rk = sort_key(receptor_complex), limma_t = t_mod, limma_FDR = FDR),
            by = c("contrast", "source", "target", "lk", "rk"))
readr::write_csv(edge_decomposition |> dplyr::select(-lk, -rk), file.path(path_tables, "ccc_diff_edge_decomposition.csv"))
cellchat_ec_cm <- edge_decomposition |> filter(source %in% FOCUS_TYPES, target %in% FOCUS_TYPES, source != target) |>
  arrange(factor(contrast, levels = CONTRAST_ORDER), desc(abs(delta))) |>
  transmute(contrast = CONTRAST_LABEL[contrast], direction = sprintf("%s\u2192%s", CELL_TYPE_SHORT[source], CELL_TYPE_SHORT[target]), interaction_name, pathway,
            delta = round(delta, 4), share_of_edge = round(share_of_edge, 2), limma_t = round(limma_t, 2), limma_FDR = signif(limma_FDR, 3))
knitr::kable(cellchat_ec_cm, caption = sprintf("EC-CM pairs of the CellChat edge decomposition with their per-animal moderated t and FDR (NA = no matching pair-direction eligible in that contrast; share_of_edge = pair delta / edge total, outside [-1, 1] where pair deltas of opposite sign cancel). %d of %d entries matched; %d at FDR < %s.",
                                       sum(!is.na(cellchat_ec_cm$limma_t)), nrow(cellchat_ec_cm), sum(cellchat_ec_cm$limma_FDR < FDR_CUT, na.rm = TRUE), FDR_CUT))

ec_cm_diff <- diff_all |> mutate(lr_pair = paste(ccc_pair, ligand_complex, receptor_complex))
kept_pairs_tab <- ec_cm_diff |> group_by(lr_pair) |> summarise(max_t = max(abs(t_mod)), .groups = "drop") |>
  slice_max(max_t, n = TOP_N, with_ties = FALSE)
kept_pairs <- kept_pairs_tab$lr_pair
heat_df <- ec_cm_diff |> filter(lr_pair %in% kept_pairs) |>
  dplyr::select(lr_pair, source, target, ligand_complex, receptor_complex, contrast, logFC, passes) |>
  complete(nesting(lr_pair, source, target, ligand_complex, receptor_complex), contrast) |>
  left_join(kept_pairs_tab, by = "lr_pair") |>
  mutate(mark  = ifelse(passes %in% TRUE, "*", ""),
         label = sprintf("%s: %s\u2192%s", CELL_TYPE_SHORT[source], ligand_complex, receptor_complex))
heat_df$label <- factor(heat_df$label, levels = unique(heat_df$label[order(heat_df$max_t)]))
lfc_limit <- max(abs(heat_df$logFC), na.rm = TRUE)

fig <- ggplot(heat_df, aes(contrast, label, fill = logFC)) +
  geom_tile(colour = "white", linewidth = .8) +
  geom_text(aes(label = mark), fontface = "bold", size = 4, vjust = .72) +
  facet_grid(target ~ ., scales = "free_y", space = "free_y", labeller = labeller(target = function(x) paste("to", CELL_TYPE_SHORT[x]))) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, limits = c(-lfc_limit, lfc_limit),
                       name = "difference of\nmean log1p score", na.value = "grey85") +
  scale_x_discrete(drop = FALSE, labels = function(x) sub(" \\(", "\n(", CONTRAST_LABEL[x])) +
  labs(x = NULL, y = NULL, title = "EC \u2194 CM interaction changes across aging and HFpEF",
       subtitle = sprintf("rows: top %d by |moderated t|; * BH FDR < %s within contrast; grey = not eligible for that contrast", TOP_N, FDR_CUT)) +
  theme(panel.grid = element_blank(), axis.text.y = element_text(size = 7))
save_svg(fig, "ccc_ecm_differential_heatmap.svg", width = 12, height = 7, dir = path_figures)
fig

lead_pairs <- ec_cm_diff |> filter(lr_pair %in% head(kept_pairs, 6)) |> distinct(across(all_of(KEY)))
set.seed(SEED)
dot_df <- lr_scores |> inner_join(lead_pairs, by = KEY) |>
  mutate(panel = factor(sprintf("%s\u2192%s\n%s\u2192%s", CELL_TYPE_SHORT[source], CELL_TYPE_SHORT[target], ligand_complex, receptor_complex),
                        levels = unique(sprintf("%s\u2192%s\n%s\u2192%s", CELL_TYPE_SHORT[lead_pairs$source], CELL_TYPE_SHORT[lead_pairs$target], lead_pairs$ligand_complex, lead_pairs$receptor_complex))),
         group = factor(disease, levels = GROUP_ORDER), x_jit = as.numeric(group) + runif(n(), -.12, .12))
oh_top   <- dot_df |> filter(disease == "OH") |> group_by(panel) |> summarise(top_animal = sample[which.max(expr_prod)], .groups = "drop")
label_df <- dot_df |> filter(qc_flag != "pass" | animal_note != "") |>
  mutate(nudge = if_else(group == last(GROUP_ORDER), -.35, .35), hjust = if_else(group == last(GROUP_ORDER), 1, 0))
fig <- ggplot(dot_df, aes(group, expr_prod, fill = group)) +
  stat_summary(fun = mean, geom = "col", width = .7, colour = "grey20", linewidth = .3) +
  geom_point(aes(x = x_jit, shape = sex), size = 1.8, alpha = .9) +
  ggrepel::geom_text_repel(data = label_df, aes(x = x_jit, label = sample), nudge_x = label_df$nudge, hjust = label_df$hjust,
                           size = 2.2, direction = "y", min.segment.length = 0, segment.size = .2, segment.colour = "grey50", seed = SEED) +
  facet_wrap(~ panel, scales = "free_y", nrow = 1) +
  scale_x_discrete(drop = FALSE, labels = function(x) sub(" ", "\n", GROUP_LABEL[x])) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  scale_fill_manual(values = GROUP_FILL, guide = "none") +
  scale_shape_manual(values = c(F = 21, M = 24), name = "sex") +
  labs(x = NULL, y = "expression-product score", title = "Leading EC\u2013CM candidates by |moderated t|: per-animal scores") +
  theme(axis.text.x = element_text(size = 7), legend.position = "bottom")
save_svg(fig, "ccc_ecm_animal_dots.svg", width = 13, height = 4.7, dir = path_figures)
fig

synthesis <- diff_all |> group_by(contrast) |>
  summarise(tested = n(), passing = sum(passes), min_FDR = signif(min(FDR), 2), .groups = "drop") |>
  mutate(contrast = CONTRAST_LABEL[as.character(contrast)])
knitr::kable(synthesis, caption = sprintf("Eligible EC-CM pair-directions per contrast, hits at BH FDR < %s and the smallest FDR.", FDR_CUT))

top_share <- lr_scores |> inner_join(lead_pairs, by = KEY) |> group_by(across(all_of(KEY))) |>
  summarise(top_animal = sample[which.max(expr_prod)], top_animal_share = round(max(expr_prod) / sum(expr_prod), 2), .groups = "drop") |>
  mutate(top_animal_note = dplyr::coalesce(sub(":.*", "", unname(QC_FLAGGED_ANIMALS[top_animal])), sub(" \\(.*", "", unname(ANIMAL_NOTES[top_animal])), ""))

lead_table <- ec_cm_diff |> filter(lr_pair %in% head(kept_pairs, 6)) |> group_by(lr_pair) |> slice_max(abs(t_mod), n = 1, with_ties = FALSE) |> ungroup() |>
  arrange(desc(abs(t_mod))) |>
  transmute(direction = sprintf("%s\u2192%s", CELL_TYPE_SHORT[source], CELL_TYPE_SHORT[target]), interaction = paste(ligand_complex, "\u2192", receptor_complex),
            strongest_contrast = CONTRAST_LABEL[as.character(contrast)], logFC = round(logFC, 3), t_mod = round(t_mod, 2),
            FDR = signif(FDR, 3), `detected YC/OC/YH/OH` = ndet_str(ndet_YC, ndet_OC, ndet_YH, ndet_OH),
            groups_supporting = (ndet_YC >= MIN_ANIMALS) + (ndet_OC >= MIN_ANIMALS) + (ndet_YH >= MIN_ANIMALS) + (ndet_OH >= MIN_ANIMALS),
            across(all_of(KEY))) |>
  left_join(top_share, by = KEY) |> dplyr::select(-all_of(KEY))
knitr::kable(dplyr::rename(lead_table, `top_animal_share (proportion)` = top_animal_share), caption = sprintf("Leading EC-CM candidates: strongest contrast per pair-direction, per-animal detection support (groups_supporting = groups with >= %d detected animals), and the animal contributing most to the summed score (share of the 16-animal total; 1/16 = 0.06 would be an even contribution).", MIN_ANIMALS))
noted_leads <- lead_table |> filter(top_animal_note != "")
note_text <- noted_leads |> group_by(top_animal, top_animal_note) |>
  summarise(txt = sprintf(" %s (%s) contributes the largest single-animal share to %s (share %s of the summed score, proportion).", top_animal[1], top_animal_note[1],
                          paste(sprintf("%s %s", direction, interaction), collapse = ", "), paste(unique(range(top_animal_share)), collapse = "\u2013")), .groups = "drop") |>
  pull(txt) |> paste(collapse = "")
flagged_leads <- lead_table |> filter(top_animal %in% names(QC_FLAGGED_ANIMALS))
exon_groups   <- animal_qc |> group_by(disease) |> summarise(ec = mean(exon_prop_ec), cm = mean(exon_prop_cm), .groups = "drop")
single_group_leads  <- lead_table |> filter(groups_supporting == 1)
single_group_text  <- if (nrow(single_group_leads)) sprintf(" %s %s eligible through the detection of a single group only.", paste(sprintf("%s %s", single_group_leads$direction, single_group_leads$interaction), collapse = " and "),
                                        if (nrow(single_group_leads) == 1) "is" else "are") else ""
lead_text <- with(lead_table, paste(sprintf("%s %s (%s, t = %.1f, FDR = %.2g, detected %s)", direction, interaction, strongest_contrast, t_mod, FDR, `detected YC/OC/YH/OH`), collapse = "; "))
hits_text <- if (any(diff_all$passes)) sprintf("%d pair-direction x contrast combinations pass FDR < %s.", sum(diff_all$passes), FDR_CUT) else
  sprintf("No pair-direction passes FDR < %s in any contrast (smallest FDR %s).", FDR_CUT, signif(min(diff_all$FDR), 2))

n_interactions_of <- function(g, pathway_id) sum(path_flow$n_interactions[path_flow$group == g & path_flow$pathway == pathway_id])
flow_lead <- bind_rows(lapply(CONTRASTS, function(contrast_def) top_pathways(contrast_def$i1, contrast_def$i2, n = 5) |>
  transmute(contrast = factor(CONTRAST_LABEL[[contrast_def$name]], levels = CONTRAST_LABEL), pathway,
            rel_flow_grp1 = round(.data[[contrast_def$i1]], 3), rel_flow_grp2 = round(.data[[contrast_def$i2]], 3), delta = round(delta, 3),
            n_int_grp1 = vapply(pathway, function(pathway_id) n_interactions_of(contrast_def$i1, pathway_id), integer(1)),
            n_int_grp2 = vapply(pathway, function(pathway_id) n_interactions_of(contrast_def$i2, pathway_id), integer(1)))))
knitr::kable(dplyr::rename(flow_lead, `rel_flow_grp1 (proportion)` = rel_flow_grp1, `rel_flow_grp2 (proportion)` = rel_flow_grp2, `delta (proportion)` = delta), row.names = FALSE, caption = "CellChat: five pathways with the largest change in relative information flow per contrast (grp1 = first group of the contrast; n_int = inferred interactions of the pathway in that group).")
pathway_text <- flow_lead |> group_by(contrast) |> summarise(txt = sprintf("%s: %s", contrast[1], paste(sprintf("%s (%+.3f)", pathway, delta), collapse = ", ")), .groups = "drop") |> pull(txt) |> paste(collapse = "; ")
n_active_groups  <- c(flow_lead$n_int_grp1, flow_lead$n_int_grp2); n_active_groups <- n_active_groups[n_active_groups > 0]
one_sided <- flow_lead |> filter(n_int_grp1 == 0 | n_int_grp2 == 0)
pathway_support_text <- sprintf("%d\u2013%d inferred interactions per group where active%s", min(n_active_groups), max(n_active_groups),
                  if (nrow(one_sided)) sprintf("; %s %s active in one group of the contrast only, a presence/absence call rather than a graded change",
                                               paste(sprintf("%s (%s)", one_sided$pathway, one_sided$contrast), collapse = ", "), if (nrow(one_sided) == 1) "is" else "are") else "")

session_footer()
