
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

if (!exists("PROJECT_ROOT")) PROJECT_ROOT <- normalizePath(file.path(getwd(), ".."), mustWork = FALSE)
PATHS <- list(
  raw        = file.path(PROJECT_ROOT, "data", "raw"),
  objects    = file.path(PROJECT_ROOT, "data", "objects"),
  resources  = file.path(PROJECT_ROOT, "data", "resources"),
  results    = file.path(PROJECT_ROOT, "results"),
  figures    = file.path(PROJECT_ROOT, "figures")
)
resolve_path <- function(p) if (grepl("^(/|[A-Za-z]:)", p)) p else file.path(PROJECT_ROOT, p)
split_csv_param <- function(x) if (nzchar(x)) trimws(strsplit(x, ",")[[1]]) else NULL

input_md5 <- function(params) {
  files <- unlist(params)
  files <- vapply(files[grepl("\\.(rds|csv)$", files)], resolve_path, "")
  unname(tools::md5sum(c(files, file.path(PROJECT_ROOT, "notebooks", "_common.R"))))
}
results_dir <- function(notebook_dir) {
  dir_path <- file.path(PATHS$results, notebook_dir); dir.create(dir_path, showWarnings = FALSE, recursive = TRUE); dir_path
}

GROUP_ORDER <- c("YC", "OC", "YH", "OH")
GROUP_LABEL <- c(YC = "young control", OC = "old control", YH = "young HFpEF", OH = "old HFpEF")
GROUP_FILL  <- c(YC = "#9ecae1", OC = "#1f78b4", YH = "#fdbf6f", OH = "#e31a1c")
GROUP_LONG_ORDER <- unname(GROUP_LABEL[GROUP_ORDER])

CONTRAST_LABEL <- c(OH_vs_OC    = "HFpEF (old)",
                    YH_vs_YC    = "HFpEF (young)",
                    OC_vs_YC    = "Age (control)",
                    OH_vs_YH    = "Age (HFpEF)",
                    Interaction = "Age x HFpEF")
CONTRAST_ORDER  <- names(CONTRAST_LABEL)
CONTRAST_FILL <- c(OH_vs_OC = "#D55E00", YH_vs_YC = "#E69F00", OC_vs_YC = "#0072B2", OH_vs_YH = "#56B4E9", Interaction = "#009E73")
CONTRAST_LABEL_WRAP <- c(OH_vs_OC = "HFpEF\n(old)", YH_vs_YC = "HFpEF\n(young)", OC_vs_YC = "Age\n(control)",
                         OH_vs_YH = "Age\n(HFpEF)", Interaction = "Age x\nHFpEF")
CONTRAST_CODE <- c(OH_vs_OC = "HFo", YH_vs_YC = "HFy", OC_vs_YC = "Ac", OH_vs_YH = "Ah", Interaction = "AxH")
CONTRASTS <- list(
  list(name = "OH_vs_OC", i1 = "OH", i2 = "OC", meaning = "HFpEF effect in old"),
  list(name = "YH_vs_YC", i1 = "YH", i2 = "YC", meaning = "HFpEF effect in young"),
  list(name = "OC_vs_YC", i1 = "OC", i2 = "YC", meaning = "age effect in control"),
  list(name = "OH_vs_YH", i1 = "OH", i2 = "YH", meaning = "age effect in HFpEF")
)

make_contrast_matrix <- function(design) {
  limma::makeContrasts(
    YH_vs_YC    = YH - YC,
    OH_vs_OC    = OH - OC,
    OC_vs_YC    = OC - YC,
    OH_vs_YH    = OH - YH,
    Interaction = (OH - OC) - (YH - YC),
    levels = design)
}
make_group <- function(age, condition) {
  g <- dplyr::case_when(age == "young" & condition == "control" ~ "YC",
                        age == "young" & condition == "HFpEF"   ~ "YH",
                        age == "old"   & condition == "control" ~ "OC",
                        age == "old"   & condition == "HFpEF"   ~ "OH",
                        TRUE ~ NA_character_)
  if (anyNA(g)) stop("Unexpected age/condition labels while constructing the four groups.")
  factor(g, levels = c("YC", "YH", "OC", "OH"))
}
DESIGN_COVARIATES <- c("sex", "seq")
DESIGN_FORMULA    <- reformulate(c("0", "group", DESIGN_COVARIATES))

CELL_TYPE_ORDER <- c("Cardiomyocyte", "Endothelial", "Fibroblast", "Pericyte", "VSMC",
              "Macrophage", "Dendritic_cell", "Lymphocyte", "Lymphatic_EC", "Glia", "Epicardial", "Proliferating", "B_cell")

CELL_TYPE_SHORT <- c(Cardiomyocyte = "CM", Endothelial = "EC", Fibroblast = "FB", Pericyte = "PC", VSMC = "SMC", Macrophage = "Mac",
                     Dendritic_cell = "DC", Lymphocyte = "Lym", Lymphatic_EC = "LEC", Glia = "Glia", Epicardial = "Epi", Proliferating = "Prolif", B_cell = "B")

read_animal_notes <- function(kind) {
  f <- file.path(PROJECT_ROOT, "config", "animal_notes.csv")
  if (!file.exists(f)) return(setNames(character(), character()))
  d <- utils::read.csv(f, stringsAsFactors = FALSE)
  d <- d[d$kind == kind, ]
  setNames(d$note, d$animal)
}
QC_FLAGGED_ANIMALS <- read_animal_notes("qc_flag")

ANIMAL_NOTES <- read_animal_notes("biological_note")

is_mt_gene    <- function(g) grepl("^mt-", g)
is_ribo_gene  <- function(g) grepl("^Rp[ls]\\d|^Rp[ls]p\\d|^Rplp", g)
is_rrna_overlap_gene <- function(g) g %in% c("Gm42418", "AY036118", "Gm26917", "Gm47283")
is_technical_gene <- function(g) is_mt_gene(g) | is_ribo_gene(g) | is_rrna_overlap_gene(g)
is_sex_gene <- function(g) g %in% c("Xist", "Tsix", "Ddx3y", "Eif2s3y", "Kdm5d", "Uty")

WARD_MERGE_FRACTION <- 0.25
ward_merge <- function(obj, cluster_col, features, exclude = character()) {
  cluster_ids   <- as.character(obj[[cluster_col, drop = TRUE]])
  kept_clusters <- setdiff(unique(cluster_ids), exclude)
  kept_clusters <- kept_clusters[order(as.numeric(kept_clusters))]
  expr <- Seurat::GetAssayData(obj, layer = "data")[features, , drop = FALSE]
  profile_dist <- dist(t(sapply(kept_clusters, function(k) Matrix::rowMeans(expr[, cluster_ids == k, drop = FALSE]))))
  tree  <- hclust(profile_dist, method = "ward.D2")
  cut_h <- WARD_MERGE_FRACTION * max(profile_dist)
  merge_group   <- cutree(tree, h = cut_h)
  merged_ids <- unname(tapply(names(merge_group), merge_group, paste, collapse = "+")[as.character(merge_group)])
  merge_of <- c(setNames(merged_ids, names(merge_group)), setNames(exclude, exclude))
  list(tree = tree, cut = cut_h, max_dist = max(profile_dist), merge_of = merge_of[order(as.numeric(names(merge_of)))])
}
plot_ward_merge <- function(wm, title) {
  plot(wm$tree, hang = -1, main = title, sub = "", xlab = "Notebook 2 sub-cluster",
       ylab = "Ward.D2 height (Euclidean, mean log-expression)", cex = 0.9)
  abline(h = wm$cut, lty = 2, col = "firebrick")
  mtext(sprintf("cut = %.2f x largest profile distance (%.2f)", WARD_MERGE_FRACTION, wm$max_dist), side = 3, line = 0.2, cex = 0.8, col = "firebrick")
}

theme_house <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(panel.grid.major = element_line(colour = "grey92", linewidth = 0.25),
          panel.grid.minor = element_blank(),
          panel.border     = element_rect(colour = "grey70", fill = NA, linewidth = 0.3),
          axis.ticks       = element_line(colour = "grey70", linewidth = 0.3),
          axis.line        = element_blank(),
          strip.background = element_rect(fill = "grey95", colour = NA),
          strip.text = element_text(face = "bold", size = base_size - 1),
          plot.subtitle = element_text(size = base_size - 2, colour = "grey35"))
}
theme_set(theme_house())

OKABE_ITO <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
               "#999999", "#000000", "#8C564B", "#17BECF", "#BCBD22", "#7F7F7F")
subtype_palette <- function(levels) {
  stopifnot(length(levels) <= length(OKABE_ITO))
  setNames(OKABE_ITO[seq_along(levels)], levels)
}
TECH_GREY <- "grey75"
CELL_TYPE_FILL   <- setNames(c("#D55E00", "#0072B2", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#882255",
                        "#F0E442", "#8C564B", "#17BECF", "#999933", "#BCBD22", "#000000"), CELL_TYPE_ORDER)
CELL_TYPE_LABEL  <- setNames(gsub("_", " ", CELL_TYPE_ORDER), CELL_TYPE_ORDER)

save_svg <- function(plot, filename, width, height, dir = NULL, ...) {
  if (!is.null(dir)) { dir.create(dir, showWarnings = FALSE, recursive = TRUE); filename <- file.path(dir, filename) }
  ggsave(filename, plot, width = width, height = height, device = "svg", bg = "white", ...)
  invisible(filename)
}
fdr_stars <- function(fdr) cut(fdr, c(-Inf, .001, .01, .05, Inf), labels = c("***", "**", "*", ""))
short_subtype <- function(x) sub("Interferon-response", "IFN", sub(" EC", "", x))
format_pvalue <- function(p) ifelse(p < 1e-100, "<1e-100", formatC(p, format = "e", digits = 2))

fit_composition <- function(subtype_props, animal_meta, covariates = NULL) {
  animal_meta_ordered <- animal_meta[match(colnames(subtype_props$TransformedProps), as.character(animal_meta$sample)), ]
  stopifnot(identical(as.character(animal_meta_ordered$sample), colnames(subtype_props$TransformedProps)))
  animal_meta_ordered$group <- make_group(animal_meta_ordered$age, animal_meta_ordered$condition)
  design_matrix <- model.matrix(reformulate(c("0", "group", covariates)), data = animal_meta_ordered)
  colnames(design_matrix)[seq_along(levels(animal_meta_ordered$group))] <- levels(animal_meta_ordered$group)
  fit <- limma::lmFit(subtype_props$TransformedProps, design_matrix)
  fit <- limma::contrasts.fit(fit, make_contrast_matrix(design_matrix))
  fit <- limma::eBayes(fit, robust = TRUE, trend = FALSE)
  t_crit <- stats::qt(0.975, df = fit$df.total)
  composition_res <- dplyr::bind_rows(lapply(seq_len(ncol(fit$coefficients)), function(i) {
    std_error <- fit$stdev.unscaled[, i] * sqrt(fit$s2.post)
    data.frame(subtype = rownames(fit$coefficients), contrast = colnames(fit$coefficients)[i],
               logit_effect = fit$coefficients[, i], SE = std_error,
               CI_low = fit$coefficients[, i] - t_crit * std_error, CI_high = fit$coefficients[, i] + t_crit * std_error,
               p_value = fit$p.value[, i], n_animals = ncol(subtype_props$TransformedProps),
               stringsAsFactors = FALSE)
  })) |>
    dplyr::group_by(contrast) |> dplyr::mutate(FDR = p.adjust(p_value, "BH")) |> dplyr::ungroup() |>
    dplyr::mutate(FDR_all_contrasts = p.adjust(p_value, "BH"),
                  contrast = factor(contrast, levels = CONTRAST_ORDER))
  composition_res
}

print_versions <- function(pkgs) cat(R.version.string, sprintf("%s %s", pkgs, vapply(pkgs, function(p) as.character(packageVersion(p)), "")), sep = "\n")
session_footer <- function() print(sessionInfo())
