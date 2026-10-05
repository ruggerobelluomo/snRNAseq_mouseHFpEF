# 07b_EC_CM_subtype_CCC.R
# Code-only export of the analysis pipeline. Data are not included; see README.md.

params <- list(ec_rds = "results/03_ec_annotation/Endothelial_subcluster_annotated_clean.rds", 
    cm_rds = "results/04_cm_annotation/Cardiomyocyte_subcluster_annotated_clean.rds", 
    broad_diff = "results/07_ccc/tables/lr_diff_all.csv", lr_resource = "data/resources/mouseconsensus_LR.csv", 
    out_dir = "results/07b_ec_cm_subtype_ccc", n_perms = 100L, group_col = "disease")


knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = TRUE, fig.align = "center", dpi = 150)
options(width = 110)
source("_common.R")
knitr::opts_chunk$set(cache = TRUE, autodep = TRUE, cache.lazy = FALSE, cache.extra = input_md5(params))

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(limma); library(CellChat)
  library(tidyr); library(patchwork); library(readr); library(tibble)
})

MIN_CELLS     <- 10
EXPR_PROP     <- 0.05
DET_ANIMAL    <- 0.05
MIN_ANIMALS   <- 3
DOMINANCE_MAX <- 0.5
FDR_CUT       <- 0.05
CC_TYPE       <- "truncatedMean"; CC_TRIM <- EXPR_PROP
CC_PVAL       <- 0.05
N_PERMS       <- params$n_perms
SEED          <- 1
TOP_N         <- 40
TOP_LR_EDGE   <- 3
CANONICAL_LIGANDS <- c("Nrg1", "Dll4", "Jag1", "Igf1", "Bmp6")
format_count <- function(x) formatC(x, big.mark = ",", format = "d")
set.seed(SEED)

path_out     <- resolve_path(params$out_dir)
path_tables <- file.path(path_out, "tables"); path_figures <- file.path(path_out, "figures")
invisible(lapply(c(path_tables, path_figures), dir.create, showWarnings = FALSE, recursive = TRUE))

print_versions(c("ggplot2", "dplyr", "Seurat", "Matrix", "limma", "CellChat", "tidyr", "patchwork", "readr", "tibble"))

read_compartment <- function(path) {
  compartment_obj <- readRDS(resolve_path(path))
  subtype_col     <- grep("_subtype$", colnames(compartment_obj@meta.data), value = TRUE)
  compartment_obj$subtype <- as.character(compartment_obj@meta.data[[subtype_col]])
  list(obj = DietSeurat(compartment_obj, assays = "RNA", layers = "data"), column = subtype_col,
       levels = levels(factor(compartment_obj@meta.data[[subtype_col]])))
}
ec_part <- read_compartment(params$ec_rds); cm_part <- read_compartment(params$cm_rds)
EC_TYPES <- ec_part$levels; CM_TYPES <- cm_part$levels
SUBTYPE_ORDER <- c(EC_TYPES, CM_TYPES)
sub_obj <- JoinLayers(merge(ec_part$obj, cm_part$obj))
rm(ec_part, cm_part); invisible(gc())
DefaultAssay(sub_obj) <- "RNA"
required <- c("subtype", "cell_type", "sample", "individual", "age", "condition", DESIGN_COVARIATES, "batch", params$group_col)
if (!all(required %in% colnames(sub_obj@meta.data)) || anyNA(sub_obj@meta.data[, required])) stop("Required metadata missing or incomplete")
if (anyDuplicated(colnames(sub_obj))) stop("Barcodes shared between the EC and CM objects")
sub_obj$subtype     <- factor(sub_obj$subtype, levels = SUBTYPE_ORDER)
sub_obj$lineage     <- ifelse(sub_obj$subtype %in% EC_TYPES, "Endothelial", "Cardiomyocyte")
sub_obj$disease     <- factor(as.character(sub_obj@meta.data[[params$group_col]]), levels = GROUP_ORDER)
sub_obj$sample      <- as.character(sub_obj$sample)
sub_obj$individual  <- as.character(sub_obj$individual)
sub_obj$qc_flag     <- dplyr::coalesce(sub(":.*", "", unname(QC_FLAGGED_ANIMALS[sub_obj$sample])), "pass")
sub_obj$animal_note <- dplyr::coalesce(sub(" \\(.*", "", unname(ANIMAL_NOTES[sub_obj$sample])), "")
SUBTYPE_COLS  <- subtype_palette(SUBTYPE_ORDER)
SUBTYPE_SHORT <- setNames(c(short_subtype(EC_TYPES), sub(" CM$", "", CM_TYPES)), SUBTYPE_ORDER)
animal_meta <- sub_obj@meta.data |>
  distinct(sample, individual, batch, seq, age, condition, disease, sex, qc_flag, animal_note) |>
  arrange(disease, sample)
rownames(animal_meta) <- NULL
if (anyDuplicated(animal_meta$sample)) stop("Expected one metadata row per animal")
knitr::kable(animal_meta, caption = "Animal-level experimental metadata, library-quality flag and biological note.")
broad_diff <- readr::read_csv(resolve_path(params$broad_diff), show_col_types = FALSE)
path_lr_resource <- resolve_path(params$lr_resource)

design <- sub_obj@meta.data |>
  group_by(disease) |>
  summarise(nuclei = n(), EC = sum(lineage == "Endothelial"), CM = sum(lineage == "Cardiomyocyte"), animals = n_distinct(individual),
            females = n_distinct(individual[sex == "F"]), males = n_distinct(individual[sex == "M"]), .groups = "drop")
knitr::kable(design, caption = "Nuclei per group (EC and CM), and their split into animals and sexes.", format.args = list(big.mark = ","))

subtype_counts  <- sub_obj@meta.data |> dplyr::count(subtype, disease, sample, name = "nuclei")
subtype_support <- subtype_counts |> group_by(subtype) |>
  summarise(nuclei_total = sum(nuclei), top_animal = sample[which.max(nuclei)], top_animal_share = max(nuclei) / sum(nuclei), .groups = "drop") |>
  left_join(subtype_counts |> group_by(subtype, disease) |> summarise(animals = sum(nuclei >= MIN_CELLS), .groups = "drop") |>
              pivot_wider(names_from = disease, values_from = animals, names_prefix = "animals_", values_fill = 0), by = "subtype") |>
  mutate(group_tested = top_animal_share < DOMINANCE_MAX)
ONE_ANIMAL_TYPES <- as.character(subtype_support$subtype[!subtype_support$group_tested])
NPPA_CM       <- grep("^Nppa-high", CM_TYPES, value = TRUE)
NPPA_CM_GROUP <- GROUP_LABEL[[as.character(animal_meta$disease[animal_meta$sample == subtype_support$top_animal[subtype_support$subtype == NPPA_CM]])]]
readr::write_csv(dplyr::rename(subtype_support, top_animal_share_prop = top_animal_share), file.path(path_tables, "subtype_animal_support.csv"))
knitr::kable(dplyr::rename(subtype_support, `top_animal_share (proportion)` = top_animal_share), digits = 3, format.args = list(big.mark = ","),
             caption = sprintf("Nuclei per subtype, the animal contributing most nuclei and its share (proportion), animals with >= %d nuclei per group (animals_YC/OC/YH/OH), and whether the subtype is group-tested (top-animal share < %s).", MIN_CELLS, DOMINANCE_MAX))

animal_qc <- sub_obj@meta.data |> dplyr::count(sample, subtype) |>
  pivot_wider(names_from = subtype, values_from = n, values_fill = 0) |>
  left_join(animal_meta |> dplyr::select(sample, disease, sex, seq, qc_flag, animal_note), by = "sample") |>
  arrange(disease, sample) |> relocate(sample, disease, sex, seq, any_of(SUBTYPE_ORDER), qc_flag, animal_note)
readr::write_csv(animal_qc, file.path(path_tables, "animal_subtype_nuclei.csv"))
knitr::kable(animal_qc, format.args = list(big.mark = ","), caption = "Nuclei per animal and subtype. Every animal is retained.")

animal_meta_fit <- animal_meta |> mutate(group = make_group(age, condition), sex = factor(sex), seq = factor(seq))
DESIGN <- model.matrix(DESIGN_FORMULA, data = animal_meta_fit); colnames(DESIGN) <- sub("^group", "", colnames(DESIGN))
RESIDUAL_DF <- nrow(DESIGN) - qr(DESIGN)$rank

path_cellchat <- file.path(path_out, "cellchat_list.rds")
CELLCHAT_DB   <- subsetDB(CellChatDB.mouse, search = c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"), key = "annotation")

run_one_cellchat <- function(subtype_obj, group_id) {
  cells <- colnames(subtype_obj)[subtype_obj$disease == group_id]
  expr_data <- GetAssayData(subtype_obj, assay = "RNA", layer = "data")[, cells, drop = FALSE]
  meta <- data.frame(labels  = droplevels(subtype_obj$subtype[match(cells, colnames(subtype_obj))]),
                     samples = factor(subtype_obj$sample[match(cells, colnames(subtype_obj))]), row.names = cells)
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
cellchat_list <- setNames(lapply(GROUP_ORDER, function(g) { cat("  CellChat:", g, "\n"); run_one_cellchat(sub_obj, g) }), GROUP_ORDER)
saveRDS(cellchat_list, path_cellchat)
invisible(capture.output(suppressMessages(cellchat_list <- lapply(cellchat_list, liftCellChat, group.new = SUBTYPE_ORDER))))

cellchat_lr_db <- CellChatDB.mouse$interaction
is_cross <- function(source, target) (source %in% EC_TYPES & target %in% CM_TYPES) | (source %in% CM_TYPES & target %in% EC_TYPES)
direction_of <- function(source) if_else(source %in% EC_TYPES, "EC\u2192CM", "CM\u2192EC")

cross_long <- bind_rows(lapply(GROUP_ORDER, function(g) {
  net <- cellchat_list[[g]]@net
  as.data.frame.table(net$prob, responseName = "prob", stringsAsFactors = FALSE) |> setNames(c("source", "target", "interaction_name", "prob")) |>
    mutate(pval = as.vector(net$pval), group = g)
})) |>
  filter(is_cross(source, target), prob > 0) |>
  mutate(inferred = pval < CC_PVAL, direction = direction_of(source),
         pathway  = cellchat_lr_db$pathway_name[match(interaction_name, rownames(cellchat_lr_db))],
         ligand   = cellchat_lr_db$ligand[match(interaction_name, rownames(cellchat_lr_db))],
         lr_label = cellchat_lr_db$interaction_name_2[match(interaction_name, rownames(cellchat_lr_db))])
readr::write_csv(cross_long, file.path(path_tables, "cellchat_ec_cm_subtype_interactions.csv"))

edge_net <- bind_rows(lapply(GROUP_ORDER, function(g) {
  net <- cellchat_list[[g]]@net
  as.data.frame.table(net$count, responseName = "count", stringsAsFactors = FALSE) |> setNames(c("source", "target", "count")) |>
    mutate(weight = as.vector(net$weight), group = g)
})) |>
  filter(is_cross(source, target)) |>
  mutate(direction = direction_of(source), group = factor(group, levels = GROUP_ORDER))
readr::write_csv(edge_net, file.path(path_tables, "cellchat_ec_cm_subtype_edges.csv"))

interaction_totals <- edge_net |> group_by(group, direction) |> summarise(inferred = sum(count), weight = sum(weight), .groups = "drop") |>
  left_join(cross_long |> filter(inferred) |> group_by(group = factor(group, levels = GROUP_ORDER), direction) |>
              summarise(lr_pairs = n_distinct(interaction_name), .groups = "drop"), by = c("group", "direction")) |>
  pivot_wider(names_from = direction, values_from = c(inferred, lr_pairs, weight), names_glue = "{direction}_{.value}")
interaction_totals <- tibble(group = factor(GROUP_ORDER, levels = GROUP_ORDER), nuclei = as.integer(table(sub_obj$disease)[GROUP_ORDER])) |>
  left_join(interaction_totals, by = "group")
knitr::kable(interaction_totals, digits = 3, format.args = list(big.mark = ","),
             caption = "Nuclei, inferred EC->CM and CM->EC interactions (sender subtype x receiver subtype x LR pair), distinct LR pairs and summed interaction weight per group.")
top_edges <- edge_net |> group_by(direction, source, target) |> summarise(weight = mean(weight), count = mean(count), .groups = "drop") |>
  group_by(direction) |> slice_max(weight, n = 3, with_ties = FALSE) |> ungroup()

tick <- function(x) ifelse(x %in% ONE_ANIMAL_TYPES, paste0(SUBTYPE_SHORT[x], " \u2020"), SUBTYPE_SHORT[x])
heat_panel <- function(dir_label, sender_types, receiver_types) {
  ggplot(edge_net |> filter(direction == dir_label), aes(factor(target, levels = receiver_types), factor(source, levels = rev(sender_types)), fill = weight)) +
    geom_tile(colour = "white", linewidth = .5) +
    geom_text(aes(label = count, colour = weight > mean(range(weight))), size = 2.6) +
    scale_colour_manual(values = c(`FALSE` = "grey15", `TRUE` = "white"), guide = "none") +
    facet_wrap(~ group, nrow = 1, labeller = as_labeller(GROUP_LABEL)) +
    scale_fill_gradient(low = "#f7fbff", high = "#08306b", name = "summed\ncommunication\nprobability") +
    scale_x_discrete(labels = tick) + scale_y_discrete(labels = tick) +
    labs(x = "receiver subtype", y = "sender subtype", title = dir_label) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 7), axis.text.y = element_text(size = 7), panel.grid = element_blank())
}
fig <- heat_panel("EC\u2192CM", EC_TYPES, CM_TYPES) / heat_panel("CM\u2192EC", CM_TYPES, EC_TYPES) +
  plot_annotation(caption = sprintf("Tile number = inferred interactions; pooled-group CellChat networks; \u2020 = described, not group-tested (>= %d %% of nuclei from one animal)", DOMINANCE_MAX * 100))
save_svg(fig, "ccc_subtype_edge_heatmap.svg", width = 13, height = 7.5, dir = path_figures)
fig

cross_mask   <- outer(SUBTYPE_ORDER, SUBTYPE_ORDER, is_cross)
cross_weight <- lapply(cellchat_list, function(x) x@net$weight[SUBTYPE_ORDER, SUBTYPE_ORDER] * cross_mask)
circle_contrasts <- c(lapply(CONTRASTS, function(x) list(label = CONTRAST_LABEL[[x$name]], red = x$i1, diff = cross_weight[[x$i1]] - cross_weight[[x$i2]])),
                      list(list(label = CONTRAST_LABEL[["Interaction"]], red = "OH \u2212 OC > YH \u2212 YC",
                                diff = (cross_weight$OH - cross_weight$OC) - (cross_weight$YH - cross_weight$YC))))

CIRCLE_COLS <- replace(SUBTYPE_COLS, c("Interferon-response EC", "Conduction-like CM"), c("#8C564B", "#BCBD22"))
CIRCLE_COLS[ONE_ANIMAL_TYPES] <- "grey70"
circle_angle <- 2 * pi * (seq_along(SUBTYPE_ORDER) - 1) / length(SUBTYPE_ORDER)
circle_label <- unname(SUBTYPE_SHORT[SUBTYPE_ORDER])
draw_diff_circle <- function(diff_w, title) {
  g <- igraph::graph_from_adjacency_matrix(diff_w, mode = "directed", weighted = TRUE)
  ends_g <- igraph::ends(g, igraph::E(g))
  edge_col <- ifelse(ends_g[, 1] %in% ONE_ANIMAL_TYPES | ends_g[, 2] %in% ONE_ANIMAL_TYPES, "grey75",
                     ifelse(igraph::E(g)$weight > 0, "#b2182b", "#2166ac"))
  plot(g, layout = cbind(cos(circle_angle), sin(circle_angle)), margin = 0.15,
       vertex.size = 20, vertex.color = unname(CIRCLE_COLS[SUBTYPE_ORDER]), vertex.frame.color = unname(CIRCLE_COLS[SUBTYPE_ORDER]),
       vertex.label = circle_label, vertex.label.cex = 0.8, vertex.label.color = "black", vertex.label.family = "Helvetica",
       vertex.label.degree = -circle_angle, vertex.label.dist = 2 + 0.22 * nchar(circle_label) * abs(cos(circle_angle)),
       edge.color = grDevices::adjustcolor(edge_col, 0.6), edge.width = 0.3 + 8 * abs(igraph::E(g)$weight) / max(abs(igraph::E(g)$weight)),
       edge.curved = 0.2, edge.arrow.size = 0.2, edge.arrow.width = 1)
  title(title, font.main = 1, cex.main = 1.1)
}

one_animal_mask <- outer(SUBTYPE_ORDER, SUBTYPE_ORDER, function(s, t) s %in% ONE_ANIMAL_TYPES | t %in% ONE_ANIMAL_TYPES)
circle_summary  <- bind_rows(lapply(circle_contrasts, function(circle_def) {
  tested_diff <- circle_def$diff * (cross_mask & !one_animal_mask)
  top <- which.max(abs(tested_diff))
  data.frame(contrast = circle_def$label, up = sum(tested_diff > 0), down = sum(tested_diff < 0),
             top_edge = sprintf("%s \u2192 %s", SUBTYPE_ORDER[row(tested_diff)[top]], SUBTYPE_ORDER[col(tested_diff)[top]]), top_diff = tested_diff[top],
             one_animal_share = sum(abs(circle_def$diff[one_animal_mask])) / sum(abs(circle_def$diff)))
}))
circle_file <- file.path(path_figures, "ccc_subtype_diff_interaction_circles.svg")
svg(circle_file, width = 15, height = 10.5)
par(mfrow = c(2, 3), xpd = TRUE, oma = c(4.5, 0, 1, 0), mar = c(1, 1, 3, 1))
for (circle_def in circle_contrasts) draw_diff_circle(circle_def$diff, sprintf("%s  (red = up in %s)", circle_def$label, circle_def$red))
mtext(c("Edge colour: red = higher summed communication probability in the first group, blue = in the second; EC \u2192 CM and CM \u2192 EC edges only; edge width scaled within each panel",
        sprintf("Grey nodes and edges: %s, subtypes with >= %d %% of their nuclei from one animal; their edges reflect essentially that one animal",
                paste(sprintf("%s (%.0f %% %s)", circle_label[match(ONE_ANIMAL_TYPES, SUBTYPE_ORDER)], 100 * subtype_support$top_animal_share[match(ONE_ANIMAL_TYPES, subtype_support$subtype)],
                              subtype_support$top_animal[match(ONE_ANIMAL_TYPES, subtype_support$subtype)]), collapse = " and "), DOMINANCE_MAX * 100),
        sprintf("EC subtypes: %s", paste(SUBTYPE_SHORT[EC_TYPES], collapse = ", ")), sprintf("CM subtypes: %s", paste(SUBTYPE_SHORT[CM_TYPES], collapse = ", "))),
      side = 1, line = c(-0.4, 0.6, 1.6, 2.6), outer = TRUE, cex = 0.7, col = "grey30")
invisible(dev.off())
knitr::include_graphics(circle_file, rel_path = FALSE)

top_lr <- cross_long |> filter(inferred) |> group_by(source, target, interaction_name) |> summarise(max_prob = max(prob), .groups = "drop") |>
  group_by(source, target) |> slice_max(max_prob, n = TOP_LR_EDGE, with_ties = FALSE) |> ungroup()
bubble_df <- cross_long |> semi_join(top_lr, by = c("source", "target", "interaction_name")) |>
  mutate(group = factor(group, levels = GROUP_ORDER),
         p_cat = cut(pval, c(-Inf, 0.01, CC_PVAL, Inf), labels = c("p < 0.01", sprintf("0.01 \u2264 p < %s", CC_PVAL), sprintf("p \u2265 %s", CC_PVAL)), right = FALSE),
         edge  = factor(sprintf("%s \u2192\n%s", SUBTYPE_SHORT[source], SUBTYPE_SHORT[target]),
                        levels = unique(sprintf("%s \u2192\n%s", SUBTYPE_SHORT[c(rep(EC_TYPES, each = length(CM_TYPES)), rep(CM_TYPES, each = length(EC_TYPES)))],
                                                SUBTYPE_SHORT[c(rep(CM_TYPES, length(EC_TYPES)), rep(EC_TYPES, length(CM_TYPES)))]))))
readr::write_csv(bubble_df, file.path(path_tables, "cellchat_subtype_top_lr.csv"))
bubble_panel <- function(dir_label, n_receivers)
  ggplot(bubble_df |> filter(direction == dir_label), aes(group, lr_label, colour = prob, size = p_cat)) +
    geom_point() +
    facet_wrap(~ edge, scales = "free_y", ncol = n_receivers, drop = TRUE) +
    scale_colour_gradient(low = "#fee0d2", high = "#a50f15", name = "communication\nprobability") +
    scale_size_manual(values = c(3, 2, 0.8), drop = FALSE, name = "permutation p") +
    scale_x_discrete(drop = FALSE) +
    labs(x = NULL, y = NULL, title = sprintf("%s: top %d LR pairs per subtype edge", dir_label, TOP_LR_EDGE)) +
    theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 6), strip.text = element_text(size = 7))
fig <- bubble_panel("EC\u2192CM", length(CM_TYPES)) / bubble_panel("CM\u2192EC", length(EC_TYPES)) + plot_layout(guides = "collect")
save_svg(fig, "ccc_subtype_top_lr_bubble.svg", width = 14, height = 16, dir = path_figures)
fig

path_lr_scores  <- file.path(path_out, "lr_scores_by_sample.csv")
ligand_receptor_all <- readr::read_csv(path_lr_resource, show_col_types = FALSE) |> filter(!is.na(ligand), !is.na(receptor)) |> distinct(ligand, receptor)
expr_data <- GetAssayData(sub_obj, assay = "RNA", layer = "data")
has_subunits <- function(complex_ids) vapply(strsplit(complex_ids, "_"), function(x) all(x %in% rownames(expr_data)), logical(1))
ligand_receptor <- ligand_receptor_all[has_subunits(ligand_receptor_all$ligand) & has_subunits(ligand_receptor_all$receptor), ]

idx_animal <- split(seq_len(ncol(sub_obj)), paste(sub_obj$sample,  sub_obj$subtype, sep = "||"))
idx_animal <- idx_animal[lengths(idx_animal) >= MIN_CELLS]
idx_group  <- split(seq_len(ncol(sub_obj)), paste(sub_obj$disease, sub_obj$subtype, sep = "||"))
mean_expr         <- sapply(idx_animal, function(i) Matrix::rowMeans(expr_data[, i, drop = FALSE]))
detect_frac       <- sapply(idx_animal, function(i) Matrix::rowMeans(expr_data[, i, drop = FALSE] > 0))
detect_frac_group <- sapply(idx_group,  function(i) Matrix::rowMeans(expr_data[, i, drop = FALSE] > 0))

reduce_complex <- function(value_mat, complex_ids) {
  m <- do.call(rbind, lapply(strsplit(complex_ids, "_"), function(su) apply(value_mat[su, , drop = FALSE], 2, min)))
  dimnames(m) <- list(complex_ids, colnames(value_mat)); m
}
ligands <- unique(ligand_receptor$ligand); receptors <- unique(ligand_receptor$receptor)
ligand_mean <- reduce_complex(mean_expr, ligands); ligand_detect <- reduce_complex(detect_frac, ligands); ligand_detect_group <- reduce_complex(detect_frac_group, ligands)
receptor_mean <- reduce_complex(mean_expr, receptors); receptor_detect <- reduce_complex(detect_frac, receptors); receptor_detect_group <- reduce_complex(detect_frac_group, receptors)
group_cols <- function(g, subtype_id) intersect(paste(g, subtype_id, sep = "||"), colnames(detect_frac_group))

score_direction <- function(source_types, target_types) {
  group_detect_max <- do.call(pmax, c(lapply(GROUP_ORDER, function(g) {
    src <- group_cols(g, source_types); tgt <- group_cols(g, target_types)
    if (length(src) && length(tgt)) pmin(ligand_detect_group[ligand_receptor$ligand, src], receptor_detect_group[ligand_receptor$receptor, tgt]) else 0 }), list(na.rm = TRUE)))
  group_eligible <- group_detect_max >= EXPR_PROP
  animals <- animal_meta$sample[paste(animal_meta$sample, source_types, sep = "||") %in% colnames(mean_expr) & paste(animal_meta$sample, target_types, sep = "||") %in% colnames(mean_expr)]
  if (!any(group_eligible) || !length(animals)) return(NULL)
  ligand_ok <- ligand_receptor$ligand[group_eligible]; receptor_ok <- ligand_receptor$receptor[group_eligible]
  source_keys <- paste(animals, source_types, sep = "||"); target_keys <- paste(animals, target_types, sep = "||")
  data.frame(sample = rep(animals, each = length(ligand_ok)), source = source_types, target = target_types, ligand_complex = ligand_ok, receptor_complex = receptor_ok,
             group_detection = rep(group_detect_max[group_eligible], times = length(animals)),
             expr_prod = as.vector(ligand_mean[ligand_ok, source_keys, drop = FALSE] * receptor_mean[receptor_ok, target_keys, drop = FALSE]),
             detected  = as.vector(ligand_detect[ligand_ok, source_keys, drop = FALSE] >= DET_ANIMAL & receptor_detect[receptor_ok, target_keys, drop = FALSE] >= DET_ANIMAL))
}
subtype_pairs <- bind_rows(expand.grid(source = EC_TYPES, target = CM_TYPES, stringsAsFactors = FALSE), expand.grid(source = CM_TYPES, target = EC_TYPES, stringsAsFactors = FALSE))
lr_scores <- bind_rows(Map(score_direction, subtype_pairs$source, subtype_pairs$target)) |>
  left_join(animal_meta |> dplyr::select(sample, individual, disease, age, condition, sex, qc_flag, animal_note), by = "sample") |>
  mutate(group_tested = !(source %in% ONE_ANIMAL_TYPES | target %in% ONE_ANIMAL_TYPES))
readr::write_csv(lr_scores, path_lr_scores)
lr_scores <- lr_scores |> mutate(disease = factor(disease, levels = GROUP_ORDER), direction = direction_of(source), ccc_pair = paste0(source, "->", target))
KEY <- c("source", "target", "ligand_complex", "receptor_complex")
scored_edges <- lr_scores |> distinct(across(all_of(KEY)), direction, group_tested)
cat(sprintf("Expression scores: %s rows, %s scored pair-directions (%s on group-tested subtypes), %d animals\n",
            format_count(nrow(lr_scores)), format_count(nrow(scored_edges)), format_count(sum(scored_edges$group_tested)), n_distinct(lr_scores$sample)))

cc_receptor_genes <- strsplit(gsub(" ", "", cellchat_lr_db$receptor.symbol), ",")
best_rate <- function(genes, subtype_ids) {
  cols  <- colnames(detect_frac_group)[sub("^.*\\|\\|", "", colnames(detect_frac_group)) %in% subtype_ids]
  genes <- intersect(genes, rownames(detect_frac_group))
  rate  <- if (length(genes)) apply(detect_frac_group[genes, cols, drop = FALSE], 2, min) else setNames(rep(0, length(cols)), cols)
  list(rate = max(rate), where = sub("^.*\\|\\|", "", names(which.max(rate))))
}
canonical_axes <- bind_rows(lapply(CANONICAL_LIGANDS, function(ligand_id) bind_rows(lapply(list(list("EC\u2192CM", EC_TYPES, CM_TYPES), list("CM\u2192EC", CM_TYPES, EC_TYPES)), function(d) {
  inferred_ligand <- cross_long |> filter(inferred, direction == d[[1]], ligand == ligand_id)
  inferred_edges  <- inferred_ligand |> distinct(source, target)
  receptor_best  <- lapply(cc_receptor_genes[cellchat_lr_db$ligand == ligand_id], best_rate, subtype_ids = d[[3]])
  ligand_best    <- best_rate(ligand_id, d[[2]])
  tibble(ligand = ligand_id, direction = d[[1]],
         cellchat_pairs = sum(cellchat_lr_db$ligand == ligand_id), liana_pairs = sum(ligand_receptor_all$ligand == ligand_id),
         max_ligand_detection = round(ligand_best$rate, 3), ligand_top_sender = ligand_best$where,
         max_receptor_detection = round(max(c(0, sapply(receptor_best, `[[`, "rate"))), 3),
         cellchat_inferred_edges = nrow(inferred_edges), cellchat_groups = n_distinct(inferred_ligand$group),
         cellchat_lr = paste(sort(unique(inferred_ligand$lr_label)), collapse = "; "),
         cellchat_edges = paste(sprintf("%s\u2192%s", SUBTYPE_SHORT[inferred_edges$source], SUBTYPE_SHORT[inferred_edges$target]), collapse = "; "),
         scored_pair_directions = sum(scored_edges$ligand_complex == ligand_id & scored_edges$direction == d[[1]]))
}))))
readr::write_csv(canonical_axes, file.path(path_tables, "canonical_axes_subtype.csv"))
knitr::kable(dplyr::rename(canonical_axes, `max_ligand_detection (proportion)` = max_ligand_detection, `max_receptor_detection (proportion)` = max_receptor_detection),
             caption = "Literature EC-CM axes at the subtype level: LR pairs per resource, highest group-level detection of the ligand in a sender subtype (and which) and of the best receptor in a receiver subtype (CellChat receptor definitions), CellChat-inferred subtype edges in any group (with the groups and LR pairs), and per-animal LIANA pair-directions scored.")
axes_inferred <- canonical_axes |> filter(cellchat_inferred_edges > 0)
axes_absent   <- canonical_axes |> filter(cellchat_inferred_edges == 0, scored_pair_directions == 0)
nrg_axis      <- canonical_axes |> filter(ligand == "Nrg1", direction == "EC\u2192CM")

per_sample <- lr_scores |> group_by(sample, source, target, direction) |>
  summarise(mean_score = mean(expr_prod), n_pairs = n(), .groups = "drop")
score_by_group <- per_sample |> left_join(animal_meta |> dplyr::select(sample, disease), by = "sample") |>
  group_by(disease, source, target, direction) |>
  summarise(n_animals = n(), mean_score = if (n_animals >= MIN_ANIMALS) mean(mean_score) else NA_real_, n_pairs = n_pairs[1], .groups = "drop")
readr::write_csv(score_by_group, file.path(path_tables, "lr_subtype_group_overview.csv"))
score_panel <- function(dir_label, sender_types, receiver_types)
  ggplot(score_by_group |> filter(direction == dir_label) |> complete(disease, source, target, fill = list(n_animals = 0L)), aes(factor(target, levels = receiver_types), factor(source, levels = rev(sender_types)), fill = mean_score)) +
    geom_tile(colour = "white", linewidth = .5) +
    geom_text(aes(label = n_animals, colour = !is.na(mean_score) & mean_score > mean(range(mean_score, na.rm = TRUE))), size = 2.4) +
    scale_colour_manual(values = c(`FALSE` = "grey30", `TRUE` = "white"), guide = "none") +
    facet_wrap(~ disease, nrow = 1, labeller = as_labeller(GROUP_LABEL)) +
    scale_fill_gradient(low = "#f7fbff", high = "#08306b", name = "mean score\nper scored\nLR pair", na.value = "grey85") +
    scale_x_discrete(labels = tick) + scale_y_discrete(labels = tick) +
    labs(x = "receiver subtype", y = "sender subtype", title = dir_label) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 7), axis.text.y = element_text(size = 7), panel.grid = element_blank())
fig <- score_panel("EC\u2192CM", EC_TYPES, CM_TYPES) / score_panel("CM\u2192EC", CM_TYPES, EC_TYPES) +
  plot_annotation(caption = sprintf("Tile number = animals with both subtype profiles (>= %d nuclei each); grey = fewer than %d; \u2020 = described, not group-tested", MIN_CELLS, MIN_ANIMALS))
save_svg(fig, "ccc_subtype_score_overview.svg", width = 13, height = 7.5, dir = path_figures)
fig

ecm_scores <- lr_scores |> filter(group_tested) |> mutate(sample = factor(sample, levels = animal_meta_fit$sample))
score_wide <- ecm_scores |> dplyr::select(all_of(KEY), sample, expr_prod) |>
  pivot_wider(names_from = sample, values_from = expr_prod, names_expand = TRUE)
detected_wide <- ecm_scores |> dplyr::select(all_of(KEY), sample, detected) |>
  pivot_wider(names_from = sample, values_from = detected, names_expand = TRUE, values_fill = FALSE)
score_mat <- as.matrix(score_wide[, animal_meta_fit$sample]); detected_mat <- as.matrix(detected_wide[, animal_meta_fit$sample])
log_score_mat <- log1p(score_mat)
contrast_groups <- c(setNames(lapply(CONTRASTS, function(x) c(x$i1, x$i2)), sapply(CONTRASTS, `[[`, "name")), list(Interaction = GROUP_ORDER))

fit_groups <- lmFit(log_score_mat, DESIGN)
fit <- eBayes(contrasts.fit(fit_groups, make_contrast_matrix(DESIGN)), robust = TRUE)
n_detected_group <- sapply(GROUP_ORDER, function(g) rowSums(detected_mat[, animal_meta_fit$sample[animal_meta_fit$group == g], drop = FALSE]))
n_scored_group   <- sapply(GROUP_ORDER, function(g) rowSums(is.finite(score_mat[, animal_meta_fit$sample[animal_meta_fit$group == g], drop = FALSE])))
group_mean <- sapply(GROUP_ORDER, function(g) rowMeans(score_mat[, animal_meta_fit$sample[animal_meta_fit$group == g], drop = FALSE], na.rm = TRUE))
colnames(n_detected_group) <- paste0("ndet_", GROUP_ORDER); colnames(n_scored_group) <- paste0("nscored_", GROUP_ORDER); colnames(group_mean) <- paste0("mean_", GROUP_ORDER)
diff_all <- bind_rows(lapply(colnames(fit$coefficients), function(contrast_name) {
  se <- fit$stdev.unscaled[, contrast_name] * sqrt(fit$s2.post); half_width <- qt(0.975, fit$df.total) * se
  data.frame(score_wide[KEY], contrast = contrast_name, logFC = fit$coefficients[, contrast_name],
             CI_low = fit$coefficients[, contrast_name] - half_width, CI_high = fit$coefficients[, contrast_name] + half_width,
             t_mod = fit$t[, contrast_name], p_value = fit$p.value[, contrast_name], n_detected_group, n_scored_group, group_mean,
             eligible = apply(n_detected_group[, paste0("ndet_", contrast_groups[[contrast_name]]), drop = FALSE], 1, max) >= MIN_ANIMALS &
                        is.finite(fit$p.value[, contrast_name]))
})) |>
  filter(eligible) |> dplyr::select(-eligible) |>
  group_by(contrast) |> mutate(FDR = p.adjust(p_value, "BH")) |> ungroup() |>
  mutate(passes = FDR < FDR_CUT, contrast = factor(contrast, levels = CONTRAST_ORDER), direction = direction_of(source), ccc_pair = paste0(source, "->", target))
readr::write_csv(diff_all, file.path(path_tables, "lr_subtype_diff_all.csv"))
test_counts <- diff_all |> group_by(contrast) |>
  summarise(tested = n(), `EC->CM` = sum(direction == "EC\u2192CM"), `CM->EC` = sum(direction == "CM\u2192EC"), subtype_edges = n_distinct(ccc_pair),
            passing = sum(passes), min_FDR = signif(min(FDR), 2), .groups = "drop") |>
  mutate(contrast = CONTRAST_LABEL[as.character(contrast)])
knitr::kable(test_counts, caption = sprintf("Tests per contrast (subtype edge x LR pair), their split by direction, the subtype edges involved, hits at BH FDR < %s and the smallest FDR.", FDR_CUT))

ndet_str <- function(YC, OC, YH, OH) sprintf("%d/%d/%d/%d", YC, OC, YH, OH)
passing_hits <- diff_all |> filter(passes) |>
  transmute(comparison = CONTRAST_LABEL[as.character(contrast)], sender = source, receiver = target,
            ligand = ligand_complex, receptor = receptor_complex,
            logFC = round(logFC, 3), CI_low = round(CI_low, 3), CI_high = round(CI_high, 3), t_mod = round(t_mod, 2),
            p_value = signif(p_value, 3), FDR = signif(FDR, 3),
            `detected YC/OC/YH/OH` = ndet_str(ndet_YC, ndet_OC, ndet_YH, ndet_OH), `scored YC/OC/YH/OH` = ndet_str(nscored_YC, nscored_OC, nscored_YH, nscored_OH)) |>
  arrange(comparison, FDR)
readr::write_csv(passing_hits, file.path(path_tables, "lr_subtype_hits.csv"))
knitr::kable(passing_hits, caption = sprintf("Subtype pair-directions with BH FDR < %s within contrast (logFC = difference of mean log1p score, 95 %% CI)%s.", FDR_CUT,
                                     if (nrow(passing_hits) == 0) ": none" else ""))

ec_cm_diff <- diff_all |> mutate(lr_pair = paste(ccc_pair, ligand_complex, receptor_complex))
kept_pairs_tab <- ec_cm_diff |> group_by(lr_pair) |> summarise(max_t = max(abs(t_mod)), .groups = "drop") |> slice_max(max_t, n = TOP_N, with_ties = FALSE)
kept_pairs <- kept_pairs_tab$lr_pair
heat_df <- ec_cm_diff |> filter(lr_pair %in% kept_pairs) |>
  dplyr::select(lr_pair, source, target, direction, ligand_complex, receptor_complex, contrast, logFC, passes) |>
  complete(nesting(lr_pair, source, target, direction, ligand_complex, receptor_complex), contrast) |>
  left_join(kept_pairs_tab, by = "lr_pair") |>
  mutate(mark = ifelse(passes %in% TRUE, "*", ""),
         label = sprintf("%s \u2192 %s: %s\u2192%s", SUBTYPE_SHORT[source], SUBTYPE_SHORT[target], ligand_complex, receptor_complex))
heat_df$label <- factor(heat_df$label, levels = unique(heat_df$label[order(heat_df$max_t)]))
lfc_limit <- max(abs(heat_df$logFC), na.rm = TRUE)
fig <- ggplot(heat_df, aes(contrast, label, fill = logFC)) +
  geom_tile(colour = "white", linewidth = .8) +
  geom_text(aes(label = mark), fontface = "bold", size = 4, vjust = .72) +
  facet_grid(direction ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, limits = c(-lfc_limit, lfc_limit),
                       name = "difference of\nmean log1p score", na.value = "grey85") +
  scale_x_discrete(drop = FALSE, labels = function(x) sub(" \\(", "\n(", CONTRAST_LABEL[x])) +
  labs(x = NULL, y = NULL, title = "EC \u2194 CM subtype interaction changes across aging and HFpEF",
       subtitle = sprintf("rows: top %d by |moderated t|; * BH FDR < %s within contrast; grey = not tested in that contrast", TOP_N, FDR_CUT)) +
  theme(panel.grid = element_blank(), axis.text.y = element_text(size = 6.5))
save_svg(fig, "ccc_subtype_differential_heatmap.svg", width = 12, height = 8, dir = path_figures)
fig

lead_pairs <- ec_cm_diff |> filter(passes | lr_pair %in% head(kept_pairs, 6)) |> distinct(across(all_of(KEY)))
set.seed(SEED)
panel_label <- function(s, t, l, r) sprintf("%s \u2192\n%s\n%s\u2192%s", SUBTYPE_SHORT[s], SUBTYPE_SHORT[t], l, r)
dot_df <- lr_scores |> inner_join(lead_pairs, by = KEY) |>
  mutate(panel = factor(panel_label(source, target, ligand_complex, receptor_complex),
                        levels = unique(panel_label(lead_pairs$source, lead_pairs$target, lead_pairs$ligand_complex, lead_pairs$receptor_complex))),
         group = factor(disease, levels = GROUP_ORDER), x_jit = as.numeric(group) + runif(n(), -.12, .12))
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
  labs(x = NULL, y = "expression-product score", title = if (any(diff_all$passes)) sprintf("Subtype EC\u2013CM pairs at FDR < %s and leading candidates by |moderated t|: per-animal scores", FDR_CUT) else
         sprintf("Leading subtype EC\u2013CM candidates by |moderated t| (none at FDR < %s): per-animal scores", FDR_CUT)) +
  theme(axis.text.x = element_text(size = 7), legend.position = "bottom")
save_svg(fig, "ccc_subtype_animal_dots.svg", width = 13, height = 4.6, dir = path_figures)
fig

top_share <- lr_scores |> inner_join(lead_pairs, by = KEY) |> group_by(across(all_of(KEY))) |>
  summarise(top_animal = sample[which.max(expr_prod)], top_animal_share = round(max(expr_prod) / sum(expr_prod), 2), n_animals = n(), .groups = "drop") |>
  mutate(top_animal_note = dplyr::coalesce(sub(":.*", "", unname(QC_FLAGGED_ANIMALS[top_animal])), sub(" \\(.*", "", unname(ANIMAL_NOTES[top_animal])), ""))
lead_table <- ec_cm_diff |> semi_join(lead_pairs, by = KEY) |> group_by(lr_pair) |> slice_max(abs(t_mod), n = 1, with_ties = FALSE) |> ungroup() |>
  arrange(desc(abs(t_mod))) |>
  transmute(edge = sprintf("%s \u2192 %s", source, target), interaction = paste(ligand_complex, "\u2192", receptor_complex),
            strongest_contrast = CONTRAST_LABEL[as.character(contrast)], logFC = round(logFC, 3), CI_low = round(CI_low, 3), CI_high = round(CI_high, 3),
            t_mod = round(t_mod, 2), FDR = signif(FDR, 3), `detected YC/OC/YH/OH` = ndet_str(ndet_YC, ndet_OC, ndet_YH, ndet_OH), across(all_of(KEY))) |>
  left_join(top_share, by = KEY) |> dplyr::select(-all_of(KEY))
knitr::kable(dplyr::rename(lead_table, `top_animal_share (proportion)` = top_animal_share),
             caption = "Hits and leading subtype candidates: strongest contrast, logFC with 95 % CI, per-animal detection, and the animal contributing most to the summed score (share of the total over the animals with both profiles, proportion).")
lead_text <- with(lead_table, paste(sprintf("%s %s (%s, logFC %.3f [%.3f, %.3f], FDR = %.2g)", edge, interaction, strongest_contrast, logFC, CI_low, CI_high, FDR), collapse = "; "))
lead_subtypes <- table(unlist(lead_pairs[c("source", "target")]))
lead_subtype  <- names(which.max(lead_subtypes))
lead_support  <- subtype_support |> filter(subtype == lead_subtype)
flagged_leads <- lead_table |> filter(top_animal %in% names(QC_FLAGGED_ANIMALS))
hits_text <- if (any(diff_all$passes)) sprintf("%d subtype pair-direction x contrast combinations pass FDR < %s.", sum(diff_all$passes), FDR_CUT) else
  sprintf("No subtype pair-direction passes FDR < %s in any contrast (smallest FDR %s).", FDR_CUT, signif(min(diff_all$FDR), 2))

session_footer()
