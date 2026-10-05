# 02_subclustering.R
# Code-only export of the analysis pipeline. Data are not included; see README.md.

params <- list(rds_in = "data/objects/integrated_LV_snRNAseq_seurat.rds", out_dir = "results/02_subclustering")


knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = TRUE, fig.width = 12, fig.height = 7)
source("_common.R")
knitr::opts_chunk$set(cache = TRUE, autodep = TRUE, cache.lazy = FALSE, cache.extra = input_md5(params))

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(harmony); library(patchwork)
})
path_rds_in    <- resolve_path(params$rds_in)
path_out   <- resolve_path(params$out_dir)
clean_rds <- sub("\\.rds$", "_clean.rds", path_rds_in)
path_figures   <- file.path(path_out, "figures"); dir.create(path_figures, showWarnings = FALSE, recursive = TRUE)

SEED <- 1; set.seed(SEED)
N_HVG          <- 2000
N_PCS          <- 20
N_NEIGHBORS    <- 15
CYTO_PC_COR    <- 0.3
RES_GRID       <- round(seq(0.1, 1.0, by = 0.1), 1)
MAX_MARKERLESS <- 1
AUC_STOP       <- 0.65
RES_OVERRIDE   <- c(Endothelial = 0.4)
OVERRIDE_GENES <- list(Endothelial = c("Isg15", "Ifit1", "Ifit2", "Ifit3", "Rsad2", "Oasl2", "Herc6", "Ifi44", "Irf7", "Usp18"))
MIN_CLUSTER_SIZE    <- 20
MIN_SUPPORT_ANIMALS <- 3
MIN_SUPPORT_NUCLEI  <- 5
MIN_LINEAGE_NUCLEI  <- 100
N_WORKERS           <- 3
CLEAN_MIN_MARGIN  <- 0.1
CLEAN_DBL_RATIO   <- 2
ARTEFACT_FAMILY   <- c("Gm10800", "Gm10801", "Gm10718", "Gm21738")
NEAR_MARGIN       <- 0.2
NEAR_DBL_RATIO    <- 1.8
CLEAN_INCOH_REL   <- 0.25
DC_PANEL          <- c("Flt3", "Clec9a", "Xcr1", "Wdfy4", "Btla", "Cd209a", "Itgax", "Ccr7", "Cadm1")
EPI_MARKERS       <- c("Muc16", "Upk1b", "Upk3b", "Bnc1", "Msln", "Lrrn4")
EPI_PANEL         <- c("Msln", "Upk3b", "Wt1", "Krt19", "Upk1b", "Lrrn4", "Bnc1")
format_count <- function(x) formatC(x, big.mark = ",", format = "d")

print_versions(c("ggplot2", "dplyr", "Seurat", "Matrix", "harmony", "patchwork"))

atlas_obj <- readRDS(path_rds_in)
DefaultAssay(atlas_obj) <- "RNA"
atlas_obj$qc_flag <- unname(coalesce(sub(":.*", "", QC_FLAGGED_ANIMALS[as.character(atlas_obj$sample)]), "pass"))
atlas_obj$group   <- make_group(atlas_obj$age, atlas_obj$condition)
atlas_obj <- AddModuleScore(atlas_obj, features = list(intersect(DC_PANEL, rownames(atlas_obj))), name = "score_Dendritic_cell",
                            ctrl = 50, nbin = 25, seed = SEED)
atlas_obj$score_Dendritic_cell <- atlas_obj$score_Dendritic_cell1; atlas_obj$score_Dendritic_cell1 <- NULL
animal_meta <- distinct(atlas_obj[[]], sample, group, sex, batch, qc_flag)
write.csv(as.data.frame(table(sample = atlas_obj$sample, cell_type = atlas_obj$cell_type)),
          file.path(path_out, "input_cell_counts_by_sample.csv"), row.names = FALSE)
knitr::kable(count(animal_meta, qc_flag, name = "animals"), caption = "Library-quality flags (every animal is retained)")

cluster_palette <- function(levels)
  setNames(colorRampPalette(OKABE_ITO)(max(length(levels), length(OKABE_ITO)))[seq_along(levels)], levels)

animal_fill <- function(nucleus_meta) {
  animals_by_group <- split(as.character(nucleus_meta$individual), nucleus_meta$group)
  unlist(unname(Map(function(animal_ids, group_id) setNames(colorRampPalette(c("grey85", GROUP_FILL[[group_id]]))(length(animal_ids) + 1)[-1], animal_ids),
                    animals_by_group, names(animals_by_group))))
}

roc_markers <- function(lineage_obj)
  FindAllMarkers(lineage_obj, only.pos = TRUE, test.use = "roc", logfc.threshold = 0, min.pct = 0.1,
                 return.thresh = 0, verbose = FALSE)

cluster_table <- function(lineage_obj, roc_all, lineage) {
  cluster_ids <- Idents(lineage_obj); cluster_levels <- levels(cluster_ids)
  cluster_median <- function(x) as.numeric(tapply(x, cluster_ids, median)[cluster_levels])
  by_animal <- table(cluster_ids, lineage_obj$individual)[cluster_levels, , drop = FALSE]
  support   <- by_animal >= MIN_SUPPORT_NUCLEI
  best <- roc_all |> group_by(cluster) |> slice_max(myAUC, n = 1, with_ties = FALSE) |> ungroup()
  score_means_df    <- aggregate(lineage_obj[[]][, grep("^score_", colnames(lineage_obj[[]]))], list(cluster = cluster_ids), mean)
  score_means <- score_means_df[match(cluster_levels, score_means_df$cluster), -1]
  other_scores  <- score_means[, colnames(score_means) != paste0("score_", lineage)]
  data.frame(cluster = cluster_levels, n = as.integer(table(cluster_ids)[cluster_levels]),
             n_animals_supporting = as.integer(rowSums(support)),
             max_animal_share = round(apply(by_animal / rowSums(by_animal), 1, max), 3),
             top_animal = colnames(by_animal)[apply(by_animal, 1, which.max)],
             median_doublet_score = round(cluster_median(lineage_obj$cellbender_doublet_scores), 3),
             median_percent_mito  = round(cluster_median(lineage_obj$percent_mito), 3),
             median_exon_prop     = round(cluster_median(lineage_obj$exon_prop), 4),
             round(score_means, 3),
             top_other_lineage       = sub("^score_", "", colnames(other_scores)[apply(other_scores, 1, which.max)]),
             top_other_lineage_score = round(apply(other_scores, 1, max), 3),
             best_marker = best$gene[match(cluster_levels, best$cluster)],
             best_AUC    = round(best$myAUC[match(cluster_levels, best$cluster)], 3), row.names = NULL) |>
    mutate(resolved    = coalesce(best_AUC, 0) >= AUC_STOP,
           too_small   = n < MIN_CLUSTER_SIZE,
           unsupported = n_animals_supporting < MIN_SUPPORT_ANIMALS)
}

subcluster_one <- function(atlas_obj, lineage, path_out) {
  message("=== ", lineage, " ===")
  rds_path <- file.path(path_out, paste0(lineage, "_subcluster.rds"))

  lineage_obj <- subset(atlas_obj, subset = cell_type == lineage)
  if (ncol(lineage_obj) < MIN_LINEAGE_NUCLEI) { message("  <", MIN_LINEAGE_NUCLEI, " nuclei, skipping"); return(invisible(NULL)) }
  message(sprintf("  %d nuclei", ncol(lineage_obj)))

  lineage_obj <- FindVariableFeatures(lineage_obj, selection.method = "vst", nfeatures = N_HVG, verbose = FALSE)
  hvg_excluded <- VariableFeatures(lineage_obj)[is_technical_gene(VariableFeatures(lineage_obj))]
  VariableFeatures(lineage_obj) <- setdiff(VariableFeatures(lineage_obj), hvg_excluded)
  lineage_obj <- ScaleData(lineage_obj, verbose = FALSE)
  lineage_obj <- RunPCA(lineage_obj, npcs = N_PCS, verbose = FALSE)
  pc_cor   <- cor(Embeddings(lineage_obj, "pca"), lineage_obj$exon_prop)[, 1]
  dims_use <- which(abs(pc_cor) <= CYTO_PC_COR)
  dropped  <- setdiff(seq_len(N_PCS), dims_use)

  set.seed(SEED)
  lineage_obj <- RunHarmony(lineage_obj, group.by.vars = "individual", dims.use = dims_use, verbose = FALSE)
  lineage_obj <- FindNeighbors(lineage_obj, reduction = "harmony", dims = seq_along(dims_use), k.param = N_NEIGHBORS, verbose = FALSE)

  partitions <- list(); sweep <- data.frame()
  for (res in RES_GRID) {
    lineage_obj <- FindClusters(lineage_obj, resolution = res, algorithm = 1, random.seed = SEED, verbose = FALSE)
    partition <- as.character(Idents(lineage_obj))
    hit  <- Position(function(x) identical(x$part, partition), partitions)
    if (is.na(hit)) {
      message(sprintf("  res %.1f: %d clusters, computing ROC markers", res, nlevels(Idents(lineage_obj))))
      roc_all <- roc_markers(lineage_obj)
      partitions[[length(partitions) + 1]] <- list(part = partition, roc = roc_all, tab = cluster_table(lineage_obj, roc_all, lineage))
      hit <- length(partitions)
    }
    cluster_tab <- partitions[[hit]]$tab
    baseline <- cluster_tab$cluster[which.max(cluster_tab$n)]
    sweep <- rbind(sweep, data.frame(
      resolution = res, n_clusters = nrow(cluster_tab),
      n_markerless = sum(!cluster_tab$resolved & cluster_tab$cluster != baseline),
      min_cluster_size = min(cluster_tab$n), n_unsupported = sum(cluster_tab$unsupported), partition_id = hit))

    eligible <- sweep$n_markerless < MAX_MARKERLESS & sweep$min_cluster_size >= MIN_CLUSTER_SIZE & sweep$n_unsupported == 0
    rule_reached <- if (eligible[1]) !all(eligible) else any(eligible)
    if (rule_reached && res >= max(RES_OVERRIDE[lineage], na.rm = TRUE, -Inf)) break
  }
  sweep <- sweep |> mutate(eligible = n_markerless < MAX_MARKERLESS & min_cluster_size >= MIN_CLUSTER_SIZE & n_unsupported == 0)
  first_bad  <- match(FALSE, sweep$eligible)
  fallback   <- identical(first_bad, 1L)
  rule_res   <- sweep$resolution[if (is.na(first_bad)) nrow(sweep) else max(first_bad - 1L, 1L)]
  chosen_res <- if (lineage %in% names(RES_OVERRIDE)) RES_OVERRIDE[[lineage]] else rule_res
  if (fallback) warning(lineage, ": the lowest resolution is already ineligible; it is retained (", chosen_res, ")")
  write.csv(sweep, file.path(path_out, paste0(lineage, "_resolution_sweep.csv")), row.names = FALSE)

  Idents(lineage_obj) <- lineage_obj$seurat_clusters <- lineage_obj[[paste0("RNA_snn_res.", chosen_res)]][, 1]
  lineage_obj <- RunUMAP(lineage_obj, reduction = "harmony", dims = seq_along(dims_use), seed.use = SEED, verbose = FALSE)
  chosen <- partitions[[sweep$partition_id[sweep$resolution == chosen_res]]]
  roc_all <- chosen$roc
  roc_top <- roc_all |> group_by(cluster) |> slice_max(myAUC, n = 15, with_ties = FALSE) |> ungroup()
  wilcox_markers <- FindAllMarkers(lineage_obj, only.pos = TRUE, test.use = "wilcox", logfc.threshold = 0.25, min.pct = 0.1, verbose = FALSE)
  wilcox_top <- wilcox_markers |> group_by(cluster) |> slice_max(avg_log2FC, n = 15, with_ties = FALSE) |> ungroup()
  cluster_summary <- chosen$tab |> mutate(resolution = chosen_res, .after = cluster)
  r2_exon_prop <- summary(lm(lineage_obj$exon_prop ~ Idents(lineage_obj)))$r.squared

  lineage_obj$subcluster <- paste0(lineage, "-", as.character(Idents(lineage_obj)))
  lineage_obj@misc$subclustering <- list(
    cell_type = lineage, status = "processed", resolution = chosen_res, rule_resolution = rule_res, fallback = fallback, n_clusters = nrow(cluster_summary),
    n_pcs_computed = N_PCS, n_pcs_used = length(dims_use), dims_use = dims_use, dropped_pcs = dropped,
    pc_exon_prop_cor = round(pc_cor, 3), hvg_excluded = hvg_excluded, n_hvg = length(VariableFeatures(lineage_obj)),
    harmony_covariate = "individual", cluster_algorithm = "Louvain (FindClusters algorithm = 1)",
    k_param = N_NEIGHBORS, seed = SEED, auc_stop = AUC_STOP, max_markerless = MAX_MARKERLESS,
    min_cluster_size = MIN_CLUSTER_SIZE, min_support = c(animals = MIN_SUPPORT_ANIMALS, nuclei = MIN_SUPPORT_NUCLEI),
    r2_exon_prop_cluster = r2_exon_prop, sweep = sweep, clusters = cluster_summary, date = format(Sys.Date()))

  dot_genes <- roc_top |> filter(!is_technical_gene(gene)) |> group_by(cluster) |>
    slice_max(myAUC, n = 3, with_ties = FALSE) |> pull(gene) |> unique()
  p_dot <- DotPlot(lineage_obj, features = dot_genes) + coord_flip() + theme_house() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
  save_svg(p_dot, paste0(lineage, "_dotplot.svg"), width = max(7, 0.6 * nrow(cluster_summary) + 4),
           height = max(5, 0.25 * length(dot_genes) + 3), dir = path_figures)

  barcode_labels <- data.frame(barcode = colnames(lineage_obj), subcluster = lineage_obj$subcluster,
                    has_specific_marker = as.character(Idents(lineage_obj)) %in% cluster_summary$cluster[cluster_summary$resolved])
  write.csv(barcode_labels,        file.path(path_out, paste0(lineage, "_subclusters.csv")), row.names = FALSE)
  write.csv(roc_top,    file.path(path_out, paste0(lineage, "_subcluster_markers.csv")), row.names = FALSE)
  write.csv(wilcox_top,    file.path(path_out, paste0(lineage, "_subcluster_markers_wilcox.csv")), row.names = FALSE)
  write.csv(rename(cluster_summary, max_animal_share_prop = max_animal_share, median_mito_pct = median_percent_mito),
            file.path(path_out, paste0(lineage, "_subcluster_summary.csv")), row.names = FALSE)
  write.csv(as.data.frame.matrix(table(cluster = Idents(lineage_obj), individual = lineage_obj$individual)),
            file.path(path_out, paste0(lineage, "_by_individual.csv")))
  saveRDS(lineage_obj, rds_path)
  invisible(lineage_obj)
}

umap_plot <- function(plot_df, palette, title, subtitle = NULL, caption = NULL, labels = setNames(names(palette), names(palette))) {
  label_bin  <- diff(range(plot_df$UMAP_1)) / 40
  barcode_labels <- plot_df |> mutate(bx = round(UMAP_1 / label_bin), by = round(UMAP_2 / label_bin)) |> count(cluster, bx, by) |>
    group_by(cluster) |> slice_max(n, n = 1, with_ties = FALSE) |> ungroup() |>
    mutate(UMAP_1 = bx * label_bin, UMAP_2 = by * label_bin, label = labels[as.character(cluster)])
  ggplot(plot_df, aes(UMAP_1, UMAP_2, colour = cluster)) +
    scattermore::geom_scattermore(pointsize = if (nrow(plot_df) > 1e4) 2 else 4) +
    ggrepel::geom_text_repel(data = barcode_labels, aes(label = label), colour = "black", size = 3, seed = SEED) +
    scale_colour_manual(values = palette, guide = "none") +
    scale_x_continuous(expand = expansion(mult = 0.08)) +
    labs(title = title, subtitle = subtitle, caption = caption, x = "UMAP 1", y = "UMAP 2") +
    theme_house() + theme(axis.text = element_blank(), axis.ticks = element_blank(), plot.caption = element_text(hjust = 0))
}

lineage_figures <- function(lineage_obj, cells) {
  run_record <- lineage_obj@misc$subclustering; lineage <- run_record$cell_type
  palette <- cluster_palette(levels(cells$cluster))
  resolution_note <- if (run_record$resolution != run_record$rule_resolution) sprintf(", override; rule %.1f", run_record$rule_resolution) else if (run_record$fallback) ", fallback" else ""
  p_umap <- umap_plot(cells, palette, paste(CELL_TYPE_LABEL[[lineage]], "sub-clusters"),
                      sprintf("res %.1f%s | %d clusters | %d/%d PCs (dropped: %s) | Harmony on individual | %s nuclei",
                              run_record$resolution, resolution_note, run_record$n_clusters, run_record$n_pcs_used, run_record$n_pcs_computed,
                              if (length(run_record$dropped_pcs)) paste(run_record$dropped_pcs, collapse = ",") else "none", format_count(nrow(cells))))
  save_svg(p_umap, paste0(lineage, "_UMAP.svg"), width = 8, height = 6.5, dir = path_figures)

  animals <- distinct(lineage_obj[[]][, c("individual", "group")]) |> arrange(group, individual)
  composition <- as.data.frame(table(cluster = cells$cluster, individual = cells$individual)) |>
    mutate(individual = factor(individual, levels = animals$individual))
  p_comp <- ggplot(composition, aes(cluster, Freq, fill = individual)) +
    geom_col(position = "fill", colour = "white", linewidth = 0.1) +
    scale_fill_manual(values = animal_fill(animals), name = "animal (shaded by group)") +
    scale_y_continuous(labels = scales::percent) +
    labs(title = "Per-animal composition", y = "share of nuclei (%)", x = "sub-cluster")
  pc_index  <- seq_len(run_record$n_pcs_computed)
  p_pc <- ggplot(data.frame(pc = factor(pc_index), r = run_record$pc_exon_prop_cor, dropped = pc_index %in% run_record$dropped_pcs), aes(pc, r, fill = dropped)) +
    geom_col() + geom_hline(yintercept = c(-1, 1) * CYTO_PC_COR, linetype = 2, colour = "grey40") +
    scale_fill_manual(values = c(`FALSE` = "grey60", `TRUE` = OKABE_ITO[6])) +
    labs(title = "PC correlation with exonic read fraction", x = "PC", y = "Pearson r")
  p_sweep <- run_record$sweep |> select(resolution, n_clusters, n_markerless, eligible) |>
    tidyr::pivot_longer(c(n_clusters, n_markerless), names_to = "metric") |>
    ggplot(aes(resolution, value, colour = metric)) +
    geom_line() + geom_point(aes(shape = eligible), size = 2.5) +
    geom_vline(xintercept = run_record$resolution, linetype = 2) +
    scale_colour_manual(values = c(n_clusters = OKABE_ITO[5], n_markerless = OKABE_ITO[6]),
                        labels = c(n_clusters = "clusters", n_markerless = "clusters without marker"), name = NULL) +
    scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16)) + scale_x_continuous(breaks = RES_GRID) +
    labs(title = sprintf("Resolution sweep (selected %.1f%s)", run_record$resolution, resolution_note), y = "number of clusters")
  save_svg((p_umap | p_comp) / (p_pc | p_sweep), paste0(lineage, "_diagnostics.svg"), width = 15, height = 11, dir = path_figures)
  p_umap
}

lineages <- intersect(CELL_TYPE_ORDER, unique(as.character(atlas_obj$cell_type)))
lineages_by_size <- names(sort(table(as.character(atlas_obj$cell_type))[lineages], decreasing = TRUE))
processed <- setNames(unlist(parallel::mclapply(lineages_by_size, function(lineage) !is.null(subcluster_one(atlas_obj, lineage, path_out)),
                                                mc.cores = N_WORKERS, mc.preschedule = FALSE)), lineages_by_size)

run_rows <- list(); cluster_tabs <- list(); sweeps <- list(); cluster_cells <- list(); pc_sd <- list(); override_info <- list()
for (lineage in lineages) {
  lineage_obj <- if (processed[[lineage]]) readRDS(file.path(path_out, paste0(lineage, "_subcluster.rds")))
  if (is.null(lineage_obj)) {
    run_rows[[lineage]] <- data.frame(cell_type = lineage, n_nuclei = sum(atlas_obj$cell_type == lineage), status = "skipped_small")
    next
  }
  run_record <- lineage_obj@misc$subclustering
  cluster_tabs[[lineage]]  <- run_record$clusters
  sweeps[[lineage]]   <- run_record$sweep
  cluster_cells[[lineage]] <- data.frame(barcode = colnames(lineage_obj), individual = lineage_obj$individual,
                               cluster = factor(as.character(Idents(lineage_obj)), levels = levels(Idents(lineage_obj))),
                               setNames(as.data.frame(Embeddings(lineage_obj, "umap")[, 1:2]), c("UMAP_1", "UMAP_2")))
  pc_sd[[lineage]] <- lineage_obj[["pca"]]@stdev
  if (lineage %in% names(RES_OVERRIDE)) {
    override_counts <- GetAssayData(lineage_obj, layer = "counts")[OVERRIDE_GENES[[lineage]], ]
    override_info[[lineage]] <- list(
      xtab = table(rule = lineage_obj[[paste0("RNA_snn_res.", run_record$rule_resolution)]][, 1], chosen = Idents(lineage_obj)),
      det  = sapply(split(seq_len(ncol(lineage_obj)), Idents(lineage_obj)), function(i) Matrix::rowMeans(override_counts[, i, drop = FALSE] > 0)),
      det3 = tapply(Matrix::colSums(override_counts > 0) >= 3, Idents(lineage_obj), mean))
  }
  run_rows[[lineage]] <- data.frame(
    cell_type = lineage, n_nuclei = ncol(lineage_obj), status = run_record$status,
    resolution = run_record$resolution, rule_resolution = run_record$rule_resolution, fallback = run_record$fallback, top_of_grid = run_record$resolution == max(RES_GRID),
    n_clusters = run_record$n_clusters, n_resolved = sum(run_record$clusters$resolved),
    n_pcs_used = run_record$n_pcs_used, dropped_pcs = paste(run_record$dropped_pcs, collapse = ";"),
    n_hvg_excluded = length(run_record$hvg_excluded), n_hvg_used = run_record$n_hvg, r2_exon_prop_cluster = round(run_record$r2_exon_prop_cluster, 3), date = run_record$date)

  cat(sprintf("\n\n## %s (res %.1f, %d clusters%s)\n\n", CELL_TYPE_LABEL[[lineage]], run_record$resolution, run_record$n_clusters,
              if (run_record$resolution != run_record$rule_resolution) sprintf(", override; rule %.1f", run_record$rule_resolution) else if (run_record$fallback) ", fallback" else ""))
  cat(sprintf("%s nuclei. %d technical genes removed from the variable features (%d used). PCs with |r| > %.1f against exon_prop dropped: %s (%d of %d PCs used).\n\n",
              format_count(ncol(lineage_obj)), length(run_record$hvg_excluded), run_record$n_hvg, CYTO_PC_COR,
              if (length(run_record$dropped_pcs)) paste(run_record$dropped_pcs, collapse = ", ") else "none", run_record$n_pcs_used, run_record$n_pcs_computed))
  cat(knitr::kable(select(run_record$sweep, -partition_id), caption = "Resolution sweep"), sep = "\n"); cat("\n\n")
  cat(knitr::kable(run_record$clusters |> select(cluster, n, n_animals_supporting, `max_animal_share (proportion)` = max_animal_share, top_animal,
                                       median_doublet_score, `median_exon_prop (proportion)` = median_exon_prop, top_other_lineage,
                                       top_other_lineage_score, best_marker, best_AUC, resolved),
                   caption = sprintf("Per-cluster diagnostics (R^2 exon_prop ~ cluster = %.3f)", run_record$r2_exon_prop_cluster)), sep = "\n")
  cat("\n\n")
  print(lineage_figures(lineage_obj, cluster_cells[[lineage]]))
  cat("\n\n")
}
rm(lineage_obj)
run_summary <- bind_rows(run_rows)
write.csv(run_summary, file.path(path_out, "run_summary.csv"), row.names = FALSE)
cat(knitr::kable(run_summary, caption = "Run summary"), sep = "\n")

processed_runs <- filter(run_summary, status == "processed")
and_list <- function(x) if (length(x) < 2) x else paste(paste(head(x, -1), collapse = ", "), "and", tail(x, 1))

override_tab <- bind_rows(lapply(names(override_info), function(lineage) {
  override <- override_info[[lineage]]
  rule_of <- apply(override$xtab, 2, function(x) rownames(override$xtab)[which.max(x)])
  cluster_tabs[[lineage]] |> transmute(cell_type = lineage, cluster, n, rule_cluster = rule_of[cluster], best_marker, best_AUC, resolved,
                             top_other_lineage, top_other_lineage_score, `>= 3 genes (%)` = as.numeric(round(100 * override$det3[cluster]))) |>
    bind_cols(as.data.frame(t(round(100 * override$det)))[cluster_tabs[[lineage]]$cluster, , drop = FALSE] |> rename_with(~ paste0(.x, " (%)")))
}))
knitr::kable(override_tab, row.names = FALSE, caption = "Overridden lineages: selected partition against the rule partition")
ec_override <- filter(override_tab, cell_type == "Endothelial")
ifn_cluster <- slice_max(ec_override, `>= 3 genes (%)`, n = 1, with_ties = FALSE)
ifn_markers <- read.csv(file.path(path_out, "Endothelial_subcluster_markers.csv")) |> filter(cluster == ifn_cluster$cluster, !is_technical_gene(gene)) |> pull(gene) |> head(5)
split_halves <- filter(ec_override, rule_cluster == rule_cluster[!resolved & n < max(n)][1])

stratum <- function(lineage) { t <- cluster_tabs[[lineage]]; stratum_row <- t[which.max(t$median_exon_prop), ]; baseline_row <- t[which.max(t$n), ]
  data.frame(stratum_cluster = stratum_row$cluster, stratum_median_exon_prop = stratum_row$median_exon_prop, baseline_cluster = baseline_row$cluster,
             baseline_median_exon_prop = baseline_row$median_exon_prop, stratum_best_marker = stratum_row$best_marker) }
strata_tab <- bind_rows(lapply(setNames(processed_runs$cell_type, processed_runs$cell_type), stratum), .id = "cell_type")
knitr::kable(rename_with(strata_tab, ~ paste(.x, "(proportion)"), ends_with("exon_prop")), row.names = FALSE, caption = "Cluster with the highest median exon_prop per lineage, against the baseline cluster")
r2_text <- processed_runs |> arrange(desc(r2_exon_prop_cluster)) |> with(and_list(sprintf("%s %.3f", CELL_TYPE_LABEL[cell_type], r2_exon_prop_cluster)))

sd_ratio_pc1_pc15  <- sapply(pc_sd, function(x) x[1] / x[15])
sd_drop_pc15_last <- sapply(pc_sd, function(x) 1 - x[N_PCS] / x[15])

fallback_lineages <- processed_runs$cell_type[processed_runs$fallback]
fallback_tab <- bind_rows(lapply(fallback_lineages, function(lineage) {
  lineage_sweep <- sweeps[[lineage]]
  data.frame(cell_type = lineage, nuclei = processed_runs$n_nuclei[processed_runs$cell_type == lineage], lineage_sweep[1, c("resolution", "n_clusters", "n_markerless", "min_cluster_size", "n_unsupported")],
             first_eligible_resolution = lineage_sweep$resolution[match(TRUE, lineage_sweep$eligible)])
}))
presence <- bind_rows(lapply(fallback_lineages, function(lineage) count(distinct(cluster_cells[[lineage]], cluster, individual), cluster, name = "animals_present") |>
                               mutate(cell_type = lineage, cluster = as.character(cluster))))
unsupported_clusters <- bind_rows(lapply(fallback_lineages, function(lineage) mutate(filter(cluster_tabs[[lineage]], unsupported), cell_type = lineage))) |>
  left_join(presence, by = c("cell_type", "cluster"))
n_animals <- n_distinct(atlas_obj$individual)
knitr::kable(fallback_tab, row.names = FALSE, caption = "Fallback lineages: the lowest resolution of the grid and the first eligible resolution (NA = none)")

prolif_tab <- cluster_tabs$Proliferating
prolif_fail <- filter(prolif_tab, score_Proliferating - top_other_lineage_score < CLEAN_MIN_MARGIN)
prolif_text <- and_list(sprintf("%s %s (%s panel %.2f, cell-cycle panel %.2f)", prolif_tab$cluster, prolif_tab$best_marker,
                               coalesce(CELL_TYPE_LABEL[prolif_tab$top_other_lineage], prolif_tab$top_other_lineage), prolif_tab$top_other_lineage_score, prolif_tab$score_Proliferating))

clean_tab <- bind_rows(lapply(names(cluster_tabs), function(lineage) {
  lineage_tab <- cluster_tabs[[lineage]]
  score_mat <- as.matrix(lineage_tab[, grep("^score_", names(lineage_tab), value = TRUE)])
  best_identity <- sub("^score_", "", colnames(score_mat)[max.col(score_mat, ties.method = "first")])
  identity <- ifelse(lineage == "Macrophage" & best_identity == "Dendritic_cell", "Dendritic_cell", lineage)
  own_col <- cbind(seq_along(identity), match(paste0("score_", identity), colnames(score_mat)))
  other_scores <- replace(score_mat, own_col, -Inf)
  outside_scores <- replace(other_scores, cbind(seq_along(identity), match(paste0("score_", lineage), colnames(score_mat))), -Inf)
  data.frame(cell_type = lineage, cluster = as.character(lineage_tab$cluster), identity = identity, n = lineage_tab$n,
             own_score = score_mat[own_col], best_other_score = apply(other_scores, 1, max),
             best_other_lineage = sub("^score_", "", colnames(other_scores)[max.col(other_scores, ties.method = "first")]),
             best_outside_score = apply(outside_scores, 1, max),
             median_doublet_score = lineage_tab$median_doublet_score, best_marker = lineage_tab$best_marker)
})) |>
  group_by(cell_type) |>
  mutate(margin        = own_score - best_other_score,
         doublet_ratio = median_doublet_score / median(rep(median_doublet_score[identity == cell_type], n[identity == cell_type])),
         mixed_lineage = identity != cell_type & best_outside_score > CLEAN_INCOH_REL * own_score,
         reason = case_when(best_marker %in% ARTEFACT_FAMILY    ~ "artefact transcripts",
                            margin < CLEAN_MIN_MARGIN           ~ "lineage identity",
                            doublet_ratio >= CLEAN_DBL_RATIO & (identity == cell_type | mixed_lineage) ~ "doublet score",
                            mixed_lineage                       ~ "mixed lineage"),
         removed = !is.na(reason),
         decision = case_when(removed ~ "removed", identity != cell_type ~ "relabelled", TRUE ~ "retained"),
         near_threshold = !removed & (margin < NEAR_MARGIN | doublet_ratio >= NEAR_DBL_RATIO)) |>
  ungroup()

cells <- bind_rows(cluster_cells, .id = "cell_type") |>
  mutate(cluster = as.character(cluster), sample = as.character(atlas_obj$sample)[match(barcode, colnames(atlas_obj))]) |>
  left_join(select(clean_tab, cell_type, cluster, removed, reason, decision), by = c("cell_type", "cluster"))
stopifnot(!anyNA(cells$removed))
removed_bc  <- cells$barcode[cells$removed]
input_label <- sub("^Unassigned$", "Epicardial", as.character(atlas_obj$cell_type))
clean_label <- replace(input_label, colnames(atlas_obj) %in% cells$barcode[cells$decision == "relabelled"], "Dendritic_cell")
atlas_clean <- subset(atlas_obj, cells = setdiff(colnames(atlas_obj), removed_bc))
atlas_clean$cell_type <- factor(clean_label[match(colnames(atlas_clean), colnames(atlas_obj))], levels = intersect(CELL_TYPE_ORDER, clean_label))
stopifnot(!anyNA(atlas_clean$cell_type))
atlas_clean@misc$cleaning <- list(
  source = basename(path_rds_in), n_input = ncol(atlas_obj), n_removed = length(removed_bc), n_retained = ncol(atlas_clean),
  min_margin = CLEAN_MIN_MARGIN, doublet_ratio = CLEAN_DBL_RATIO, incoherence_ratio = CLEAN_INCOH_REL, doublet_score = "cellbender_doublet_scores",
  artefact_family = ARTEFACT_FAMILY, removed_clusters = select(filter(clean_tab, removed), cell_type, cluster, n, reason),
  dc_panel = DC_PANEL, relabelled_clusters = select(filter(clean_tab, decision == "relabelled"), cell_type, cluster, n, identity),
  epicardial = list(from = "Unassigned", n = sum(input_label == "Epicardial"), markers = EPI_MARKERS, panel = EPI_PANEL),
  date = format(Sys.Date()))
saveRDS(atlas_clean, clean_rds)

nucleus_meta <- atlas_obj[[]] |> transmute(sample = as.character(sample), group, qc_flag, cell_type = input_label,
                             doublet_score = cellbender_doublet_scores, removed = colnames(atlas_obj) %in% removed_bc)
removal_counts <- function(lineage_cells) summarise(lineage_cells, nuclei = n(), removed = sum(removed), retained = nuclei - removed,
                               pct_removed = round(100 * removed / nuclei, 1), .groups = "drop")
by_lineage     <- nucleus_meta |> group_by(cell_type) |> removal_counts() |> arrange(match(cell_type, CELL_TYPE_ORDER))
by_animal <- nucleus_meta |> group_by(sample, group, qc_flag) |> removal_counts() |> arrange(desc(pct_removed))
by_group  <- bind_rows(nucleus_meta, mutate(nucleus_meta, cell_type = "all")) |> group_by(group, cell_type) |> removal_counts()
write.csv(clean_tab, file.path(path_out, "atlas_cleaning_clusters.csv"), row.names = FALSE)
write.csv(filter(cells, removed) |> select(barcode, cell_type, cluster, reason),
          file.path(path_out, "atlas_cleaning_removed_barcodes.csv"), row.names = FALSE)
write.csv(by_lineage, file.path(path_out, "atlas_cleaning_by_cell_type.csv"), row.names = FALSE)
write.csv(rename(by_animal, removed_pct = pct_removed), file.path(path_out, "atlas_cleaning_by_animal.csv"), row.names = FALSE)
write.csv(rename(by_group, removed_pct = pct_removed), file.path(path_out, "atlas_cleaning_by_group.csv"), row.names = FALSE)

cat(knitr::kable(filter(clean_tab, removed) |> select(cell_type, cluster, n, own_score, best_other_lineage, best_other_score,
                                                      margin, doublet_ratio, best_marker, reason),
                 digits = 2, caption = "Removed sub-clusters"), sep = "\n")
cat("\n\n"); cat(knitr::kable(filter(clean_tab, decision == "relabelled") |> select(cell_type, cluster, identity, n, own_score, best_other_lineage,
                                                      best_other_score, margin, doublet_ratio, best_outside_score, best_marker),
                 digits = 2, caption = "Sub-clusters relabelled Dendritic_cell"), sep = "\n")
epi_detect <- GetAssayData(atlas_obj, layer = "counts")[union(EPI_MARKERS, EPI_PANEL), ] > 0
epi_tab <- bind_rows(lapply(split(seq_len(ncol(atlas_obj)), input_label), function(i)
  data.frame(nuclei = length(i), t(100 * Matrix::rowMeans(epi_detect[EPI_MARKERS, i, drop = FALSE])),
             panel_2plus = 100 * mean(Matrix::colSums(epi_detect[EPI_PANEL, i, drop = FALSE]) >= 2),
             Wt1 = 100 * mean(epi_detect["Wt1", i]))), .id = "cell_type") |> arrange(match(cell_type, CELL_TYPE_ORDER))
cat("\n\n"); cat(knitr::kable(rename_with(epi_tab, ~ paste(.x, "(%)"), c(all_of(EPI_MARKERS), Wt1)) |> rename(`>= 2 EPI_PANEL genes (%)` = panel_2plus),
                              digits = 1, caption = "Epicardial markers: % of nuclei detecting each gene per broad compartment (Notebook-1 label; Epicardial = nuclei without a Notebook-1 broad label)"), sep = "\n")
cat("\n\n"); cat(knitr::kable(rename(by_lineage, `removed (%)` = pct_removed), caption = "Nuclei removed per compartment"), sep = "\n")
cat("\n\n"); cat(knitr::kable(rename(by_animal, `removed (%)` = pct_removed), caption = "Nuclei removed per animal (no animal is excluded)"), sep = "\n")
group_wide <- by_group |> mutate(cell_type = factor(cell_type, levels = c(CELL_TYPE_ORDER, "all"), labels = c(CELL_TYPE_LABEL[CELL_TYPE_ORDER], "all"))) |>
  select(group, cell_type, pct_removed) |> tidyr::pivot_wider(names_from = cell_type, values_from = pct_removed, names_sort = TRUE) |>
  arrange(match(group, GROUP_ORDER)) |> mutate(group = GROUP_LABEL[as.character(group)])
cat("\n\n"); cat(knitr::kable(group_wide, caption = "Nuclei removed per group and compartment (% of the group's nuclei in the compartment)"), sep = "\n")
cat("\n\n")

n_removed_total   <- length(removed_bc)
pct_removed_of <- function(lineage) by_lineage$pct_removed[by_lineage$cell_type == lineage]
top2_animals   <- head(by_animal, 2)
dominated <- filter(cells, removed) |> group_by(cell_type, cluster) |>
  summarise(n = n(), share = mean(sample %in% top2_animals$sample), .groups = "drop") |> filter(share > 0.5) |> arrange(desc(n))
other_animals   <- filter(by_animal, !sample %in% top2_animals$sample)
doublet_median <- nucleus_meta |> group_by(sample, qc_flag) |> summarise(med = median(doublet_score), .groups = "drop")
group_removal <- filter(by_group, cell_type == "all") |> arrange(match(group, GROUP_ORDER))
prolif_rm <- filter(clean_tab, removed, reason == "lineage identity", best_other_lineage == "Proliferating")
adipocyte_tab     <- filter(clean_tab, best_other_lineage == "Adipocyte")
near_threshold_tab      <- filter(clean_tab, near_threshold)
dc_tab    <- filter(clean_tab, identity == "Dendritic_cell")
epi_meta  <- atlas_obj[[]][input_label == "Epicardial", ]
epi_row   <- filter(epi_tab, cell_type == "Epicardial"); epi_rest <- filter(epi_tab, cell_type != "Epicardial")
epi_panel_median <- sort(sapply(epi_meta[, setdiff(grep("^score_", colnames(epi_meta), value = TRUE), "score_Dendritic_cell")], median), decreasing = TRUE)

umap_note <- list(Endothelial = "Clusters labelled Capillary EC in notebook 3")
umap_clean <- function(lineage) {
  lineage_cells    <- filter(cells, cell_type == lineage)
  cluster_levels   <- levels(cluster_cells[[lineage]]$cluster)
  retained_cells <- filter(lineage_cells, !removed) |> mutate(cluster = factor(cluster, levels = cluster_levels))
  removed_clusters <- filter(clean_tab, cell_type == lineage, removed)$cluster
  umap_plot(retained_cells, cluster_palette(cluster_levels), paste(CELL_TYPE_LABEL[[lineage]], "after cleaning"),
            sprintf("%s of %s nuclei retained | removed clusters: %s", format_count(nrow(retained_cells)), format_count(nrow(lineage_cells)),
                    if (length(removed_clusters)) paste(removed_clusters, collapse = ", ") else "none"),
            caption = umap_note[[lineage]])
}
p_clean <- lapply(setNames(names(cluster_cells), names(cluster_cells)), umap_clean)
invisible(lapply(names(p_clean), function(lineage) save_svg(p_clean[[lineage]], paste0(lineage, "_UMAP_clean.svg"), width = 7, height = 6, dir = path_figures)))
lineages_present <- intersect(CELL_TYPE_ORDER, unique(as.character(atlas_clean$cell_type)))
atlas_df <- data.frame(cluster = factor(as.character(atlas_clean$cell_type), levels = lineages_present),
                       setNames(as.data.frame(Embeddings(atlas_clean, "umap")[, 1:2]), c("UMAP_1", "UMAP_2")))
p_atlas <- umap_plot(atlas_df, CELL_TYPE_FILL[lineages_present], "Cleaned atlas", sprintf("%s nuclei", format_count(ncol(atlas_clean))), labels = CELL_TYPE_LABEL)
save_svg(p_atlas, "atlas_UMAP_clean.svg", width = 8, height = 7, dir = path_figures)
print(patchwork::wrap_plots(p_clean, ncol = 3))
print(p_atlas)

session_footer()
