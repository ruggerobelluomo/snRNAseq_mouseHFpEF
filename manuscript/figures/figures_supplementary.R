# figures_supplementary.R
# Code-only export of the analysis pipeline. Data are not included; see README.md.

if (!exists("FIG_DIR")) FIG_DIR <- normalizePath(getwd())
if (!exists("save_figure")) source(file.path(FIG_DIR, "code", "figure_style.R"))
suppressPackageStartupMessages({ library(data.table); library(ggrepel); library(scales); library(stringr) })
PANELS_DIR <- file.path(FIG_DIR, "panels"); TABLES_DIR <- file.path(FIG_DIR, "tables")
dir.create(PANELS_DIR, showWarnings = FALSE); dir.create(TABLES_DIR, showWarnings = FALSE)
DATA <- function(...) file.path(FIG_ROOT, "data", ...)
T01  <- function(f) RES("01_qc_atlas", "tables", f)
PY_EXE <- Sys.getenv("RETICULATE_PYTHON", "python3")

GREY_QC <- "#B8BDC4"; KEEP_QC <- "#2E6FB0"; REMOVE_QC <- "#C0392B"
FLAG_RED <- "#C0392B"
FLAGGED  <- names(QC_FLAGGED_ANIMALS)
NOTED    <- names(ANIMAL_NOTES)
flag_short <- function(s) unname(ifelse(s %in% FLAGGED, sub(":.*$", "", QC_FLAGGED_ANIMALS[s]), "pass"))
animals <- read.csv(RES("03_ec_annotation", "tables", "ec_sample_metadata.csv"), stringsAsFactors = FALSE) |>
  dplyr::mutate(group = factor(group, levels = GROUP_ORDER)) |> dplyr::arrange(group, sample)
SAMPLE_ORDER <- animals$sample
sample_col   <- function(s) ifelse(s %in% FLAGGED, FLAG_RED, "black")
sample_lab   <- function(s) ifelse(s %in% NOTED, paste0(s, "*"), s)
sfmt_n  <- function(x) formatC(round(x), big.mark = ",", format = "d")
fmt_f   <- function(x, d = 2) formatC(x, format = "f", digits = d)
PT6 <- 6 / .pt; PT7 <- 7 / .pt
RASTER_DPI <- 600
px_for <- function(w_mm, h_mm = w_mm) round(c(w_mm, h_mm) / 25.4 * RASTER_DPI)
theme_supp <- function() theme_fig() + theme(axis.text = element_text(size = 6), legend.text = element_text(size = 6),
                                              legend.title = element_text(size = 6), legend.key.size = unit(2.5, "mm"),
                                              legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(1, "mm"),
                                              strip.text = element_text(size = BASE_PT, face = "bold"),
                                              panel.grid.major = element_blank())
no_axes <- theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())
wrap <- function(p) wrap_elements(full = p)

supp_note <- function(fig, ...) list(...)

svg_min_font <- function(f) { s <- readLines(f, warn = FALSE); v <- as.numeric(unlist(regmatches(s, gregexpr("(?<=font-size: )[0-9.]+(?=px)", s, perl = TRUE)))); min(v) }
supp_save <- function(fig, name, width_mm, height_mm) {
  stopifnot(height_mm <= 230, width_mm %in% c(WIDTH_1COL, WIDTH_2COL))
  fl <- save_figure(fig, name, width_mm, height_mm)
  pdf_file <- file.path(FIG_OUT, paste0(name, ".pdf"))
  if (!file.exists(pdf_file) || file.mtime(pdf_file) < file.mtime(file.path(FIG_OUT, paste0(name, ".svg"))) - 60) {
    grDevices::cairo_pdf(pdf_file, width = width_mm * MM, height = height_mm * MM, family = "Arial")
    print(fig + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = LETTER_PT, face = "bold", family = FONT)))
    invisible(dev.off())
  }
  min_pt <- svg_min_font(file.path(FIG_OUT, paste0(name, ".svg")))
  if (min_pt < 6 - 1e-6) warning(sprintf("%s: smallest text %.2f pt < 6 pt", name, min_pt))
  invisible(if (is.data.frame(fl)) fl else data.frame(figure = name, width_mm = width_mm, height_mm = height_mm, min_font_pt = min_pt))
}
write_table <- function(df, name) {
  data.table::fwrite(df, file.path(TABLES_DIR, paste0(name, ".csv")))
  if (requireNamespace("openxlsx", quietly = TRUE)) openxlsx::write.xlsx(df, file.path(TABLES_DIR, paste0(name, ".xlsx")))
  invisible(nrow(df))
}

nb_const <- function(rmd, name, env = globalenv()) {
  lines <- readLines(file.path(FIG_ROOT, "notebooks", rmd), warn = FALSE)
  fence <- grepl("^```", lines); open_r <- grepl("^```\\{r", lines); inside <- cumsum(open_r) - cumsum(fence & !open_r) > 0 & !fence
  exprs <- parse(text = lines[inside], keep.source = FALSE)
  hit <- Filter(function(e) is.call(e) && identical(e[[1]], as.name("<-")) && identical(e[[2]], as.name(name)), as.list(exprs))
  stopifnot(length(hit) >= 1)
  eval(hit[[1]][[3]], envir = env)
}
Sys.setenv(RETICULATE_PYTHON = PY_EXE, RETICULATE_CONDA = file.path(dirname(dirname(dirname(dirname(PY_EXE)))), "bin", "micromamba"))
invisible(reticulate::py_config())

t0_supp <- Sys.time()
rv   <- read.csv(DATA("raw", "recount_validation.csv"), stringsAsFactors = FALSE) |> dplyr::rename(sample = library)
rvn  <- fread(DATA("raw", "recount_validation_per_nucleus.csv.gz")) |> dplyr::rename(sample = library)
lib  <- read.csv(T01("qc_validation_library_metrics.csv"), stringsAsFactors = FALSE)
sqf  <- read.csv(T01("sample_qc_flags.csv"), stringsAsFactors = FALSE)
qc   <- fread(T01("cellQC_results.txt.gz"), select = c("sample", "exon_prop", "cellranger_ncount", "outlier"))
cbm  <- dplyr::bind_rows(lapply(SAMPLE_ORDER, function(s) {
  m <- read.csv(DATA("raw", s, "cellbender_metrics.csv"), header = FALSE, col.names = c("metric", "value"))
  data.frame(sample = s, t(setNames(m$value, m$metric)), check.names = FALSE) }))
brk  <- fread(file.path(PANELS_DIR, "supp_barcode_rank.csv.gz"))
brs  <- read.csv(file.path(PANELS_DIR, "supp_barcode_rank_summary.csv"), stringsAsFactors = FALSE)
amb  <- read.csv(file.path(PANELS_DIR, "supp_ambient_detection.csv"), stringsAsFactors = FALSE)
stopifnot(setequal(rv$sample, SAMPLE_ORDER), setequal(lib$sample, SAMPLE_ORDER), nrow(rvn) == sum(rv$cellranger_nuclei),
          all(brs$umis_at_rank_25000 == sqf$umi_last_cellbender_droplet[match(brs$sample, sqf$sample)]),
          all(brs$n_called_cellbender == lib$cb_found_cells[match(brs$sample, lib$sample)]))
exon_by_lib <- qc[, .(median_exon_prop = median(exon_prop, na.rm = TRUE), q25 = quantile(exon_prop, .25, na.rm = TRUE),
                      q75 = quantile(exon_prop, .75, na.rm = TRUE), n = .N), by = sample]
time_supp_s1_data <- as.numeric(difftime(Sys.time(), t0_supp, units = "secs"))

t0_supp <- Sys.time()
yl <- function() scale_y_discrete(limits = rev(SAMPLE_ORDER), labels = sample_lab)
ycol <- theme(axis.text.y = element_text(colour = rev(sample_col(SAMPLE_ORDER)), size = 6))

pA <- ggplot(rvn, aes(umis_cellranger, umis_exon_only)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey70", linewidth = .25, linetype = "22") +
  scattermore::geom_scattermore(pointsize = 4, alpha = .3, colour = KEEP_QC, pixels = px_for(52)) +
  scale_x_log10(labels = label_number(scale_cut = cut_short_scale())) + scale_y_log10(labels = label_number(scale_cut = cut_short_scale())) +
  annotate("text", x = min(rvn$umis_cellranger), y = max(rvn$umis_exon_only), hjust = 0, vjust = 1, size = PT6, lineheight = .9,
           label = sprintf("%s nuclei, 16 libraries\nmedian ratio %s-%s\nPearson r (log1p) %s", sfmt_n(nrow(rvn)),
                           fmt_f(min(rv$median_ratio_exon_only_to_cellranger), 3), fmt_f(max(rv$median_ratio_exon_only_to_cellranger), 3),
                           fmt_f(min(rv$pearson_r_nucleus_log1p_umis), 2))) +
  labs(x = "UMIs per nucleus, Cell Ranger", y = "UMIs per nucleus,\nexon-only recount") + coord_equal() + theme_supp()

rvn[, gain := umis_intron_inclusive / umis_cellranger]
pB <- ggplot(rvn, aes(gain, sample)) +
  geom_boxplot(outlier.shape = NA, width = .6, linewidth = .25, fill = "grey92", colour = "grey30") +
  geom_point(data = rv, aes(median_ratio_intron_inclusive_to_cellranger, sample), size = .9, colour = KEEP_QC) +
  geom_text(data = rv, aes(x = 6.4, y = sample, label = fmt_f(median_ratio_intron_inclusive_to_cellranger, 2)), size = PT6, hjust = 1) +
  yl() + scale_x_continuous(limits = c(1, 6.5), breaks = 1:6) +
  labs(x = "intron-inclusive / Cell Ranger\nUMIs per nucleus (ratio)", y = NULL) + theme_supp() + ycol

pC <- ggplot(qc, aes(exon_prop, sample)) +
  geom_violin(fill = "grey85", colour = NA, scale = "width", width = .85) +
  geom_point(data = exon_by_lib, aes(median_exon_prop, sample), size = .9, colour = KEEP_QC) +
  geom_point(data = sqf, aes(x = 1 - intronic_reads_pct / 100, y = sample), shape = 4, size = 1, stroke = .4, colour = "black") +
  yl() + scale_x_continuous(limits = c(0, 1), breaks = c(0, .5, 1), labels = c("0", "0.5", "1")) +
  labs(x = "exonic fraction per nucleus\n(exon_prop, proportion)", y = NULL) + theme_supp() + ycol +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

fate_col <- c(retained = KEEP_QC, `called, removed by QC` = REMOVE_QC, `not called` = GREY_QC)
brk_p <- brk |> dplyr::mutate(sample = factor(sample, levels = SAMPLE_ORDER), fate = factor(fate, levels = rev(names(fate_col)))) |> dplyr::arrange(fate)
lines_df <- lib |> dplyr::transmute(sample = factor(sample, levels = SAMPLE_ORDER), starsolo = n_cellranger)
lab_df <- brs |> dplyr::mutate(sample = factor(sample, levels = SAMPLE_ORDER), lab = sprintf("rank 25,000:\n%s UMIs", sfmt_n(umis_at_rank_25000)))
strip_lab <- setNames(paste0(sample_lab(SAMPLE_ORDER), ifelse(SAMPLE_ORDER %in% FLAGGED, " (flagged)", "")), SAMPLE_ORDER)
pD <- ggplot(brk_p, aes(rank, umis, colour = fate)) +
  scattermore::geom_scattermore(pointsize = 3.2, pixels = px_for(42, 24)) +
  geom_vline(data = lines_df, aes(xintercept = starsolo, linetype = "STARsolo cell calls"), linewidth = .3, colour = "black") +
  geom_vline(aes(xintercept = 25000, linetype = "CellBender window (25,000 droplets)"), linewidth = .3, colour = "black") +
  geom_hline(aes(yintercept = 200, linetype = "200-UMI floor"), linewidth = .3, colour = "black") +
  geom_text(data = lab_df, aes(x = 1.5, y = 1.3, label = lab), inherit.aes = FALSE, hjust = 0, vjust = 0, size = PT6, lineheight = .9,
            colour = ifelse(lab_df$umis_at_rank_25000 >= 200, FLAG_RED, "grey25")) +
  facet_wrap(~ sample, ncol = 4, labeller = as_labeller(strip_lab)) +
  scale_x_log10(breaks = 10^(0:6), labels = label_number(scale_cut = cut_short_scale())) +
  scale_y_log10(breaks = 10^(0:5), labels = label_number(scale_cut = cut_short_scale())) +
  scale_colour_manual(values = fate_col, breaks = names(fate_col), name = NULL,
                      labels = c(retained = "retained by QC", `called, removed by QC` = "CellBender cell, removed by QC", `not called` = "not called (20,000 shown)")) +
  scale_linetype_manual(values = c(`STARsolo cell calls` = "22", `CellBender window (25,000 droplets)` = "42", `200-UMI floor` = "11"), name = NULL) +
  guides(colour = guide_legend(override.aes = list(pointsize = 5), order = 1, nrow = 1), linetype = guide_legend(order = 2, nrow = 1)) +
  labs(x = "barcode rank", y = "raw (intron-inclusive) UMIs per barcode") + theme_supp() +
  theme(legend.position = "bottom", legend.box = "vertical", legend.spacing.y = unit(0, "mm"), strip.text = element_text(size = 6, face = "bold"), panel.spacing = unit(1.5, "mm"),
        strip.text.x = element_text(colour = "black"))

calls <- lib |> dplyr::transmute(sample, `STARsolo cell calls` = n_cellranger, `CellBender droplets entering QC` = n_cellbender, `retained by QC` = n_retained) |>
  pivot_longer(-sample, names_to = "set", values_to = "n") |> dplyr::mutate(set = factor(set, levels = c("STARsolo cell calls", "CellBender droplets entering QC", "retained by QC")))
pE1 <- ggplot(calls, aes(n, sample, colour = set, shape = set)) +
  geom_line(aes(group = sample), colour = "grey80", linewidth = .3) + geom_point(size = 1.3, stroke = .5) + yl() +
  scale_colour_manual(values = c("black", "#8C8C8C", KEEP_QC), name = NULL) + scale_shape_manual(values = c(4, 1, 16), name = NULL) +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale()), limits = c(0, NA), expand = expansion(mult = c(0, .04))) +
  guides(colour = guide_legend(nrow = 2)) +
  labs(x = "barcodes per library", y = NULL) + theme_supp() + ycol + theme(legend.position = "bottom", legend.justification = "left")
pE2 <- ggplot(cbm, aes(100 * as.numeric(fraction_counts_removed), sample)) +
  geom_segment(aes(x = 0, xend = 100 * as.numeric(fraction_counts_removed), yend = sample), colour = "grey75", linewidth = .3) +
  geom_point(size = 1, colour = REMOVE_QC) + yl() +
  scale_x_continuous(limits = c(0, 26), breaks = c(0, 10, 20), expand = c(0, 0)) +
  labs(x = "counts removed by\nCellBender (%)", y = NULL) + theme_supp() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
pE <- wrap(pE1 + pE2 + plot_layout(widths = c(2.3, 1)))

amb_l <- amb |> pivot_longer(c(raw_pct, cellbender_pct), names_to = "layer", values_to = "pct") |>
  dplyr::mutate(layer = factor(layer, levels = c("raw_pct", "cellbender_pct"), labels = c("raw counts", "CellBender counts")), gene = factor(gene, levels = amb$gene))
pF <- ggplot(amb_l, aes(gene, pct, fill = layer)) +
  geom_col(position = position_dodge(width = .8), width = .75) +
  scale_fill_manual(values = c(GREY_QC, KEEP_QC), name = NULL) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, .02))) +
  labs(x = NULL, y = "atlas nuclei with\ndetected counts (%)") + theme_supp() +
  theme(axis.text.x = element_text(face = "italic", size = 6), legend.position = "top")

fig_s1 <- wrap_plots(list(pA, pB, pC), nrow = 1, widths = c(1.05, 1, .8))
fig_s1 <- fig_s1 / wrap(pD) / (pE | wrap(pF)) + plot_layout(heights = c(55, 100, 58))
figlog_S1 <- supp_save(fig_s1, "Figure_S1", WIDTH_2COL, 228)
LOG_Figure_S1 <- supp_note("Figure_S1", n_recount_nuclei = nrow(rvn), ratio_min = min(rv$median_ratio_exon_only_to_cellranger), ratio_max = max(rv$median_ratio_exon_only_to_cellranger),
          r_nuc_min = min(rv$pearson_r_nucleus_log1p_umis), r_gene_min = min(rv$pearson_r_gene_log1p_totals), r_gene_max = max(rv$pearson_r_gene_log1p_totals),
          gain_min = min(rv$median_ratio_intron_inclusive_to_cellranger), gain_max = max(rv$median_ratio_intron_inclusive_to_cellranger),
          intronic_min = min(sqf$intronic_reads_pct), intronic_max = max(sqf$intronic_reads_pct),
          exon_prop_med_min = min(exon_by_lib$median_exon_prop), exon_prop_med_max = max(exon_by_lib$median_exon_prop), n_qc_input = nrow(qc),
          umi25k = setNames(brs$umis_at_rank_25000, brs$sample), cb_removed = setNames(100 * as.numeric(cbm$fraction_counts_removed), cbm$sample),
          n_starsolo = sum(lib$n_cellranger), n_cellbender = sum(lib$n_cellbender), n_cb_found = sum(lib$cb_found_cells), n_retained = sum(lib$n_retained), amb = amb)
time_supp_fig_s1 <- as.numeric(difftime(Sys.time(), t0_supp, units = "secs"))

t0_supp <- Sys.time()
sep <- read.csv(T01("qc_validation_scrublet_separation.csv"), stringsAsFactors = FALSE)
ctrl <- read.csv(T01("qc_validation_scrublet_controls.csv"), stringsAsFactors = FALSE)
stopifnot(setequal(unique(sep$sample), SAMPLE_ORDER))
sep <- sep |> dplyr::mutate(fixed = grepl("fixed", method), src = factor(source, levels = c("cellranger", "cellbender"), labels = c("raw counts", "CellBender counts")))
fixed_cut <- sep |> dplyr::filter(fixed) |> dplyr::distinct(src, threshold_used)
stopifnot(nrow(fixed_cut) == 2)
FIXED_LIBS <- sort(unique(sep$sample[sep$fixed & sep$source == "cellranger"]))
scat <- function(yv, yl) {
  ggplot(sep, aes(auroc_sim_gt_obs, .data[[yv]])) +
    geom_vline(xintercept = .5, colour = GREY_QC, linewidth = .3) +
    geom_point(aes(shape = src, fill = fixed), size = 1.5, stroke = .25, colour = "black") +
    geom_text_repel(data = dplyr::filter(sep, fixed, source == "cellranger"), aes(label = sample), colour = FLAG_RED, size = PT6, seed = 1,
                    box.padding = .25, min.segment.length = 0, segment.size = .2) +
    scale_shape_manual(values = c(21, 22), name = NULL) +
    scale_fill_manual(values = c(`FALSE` = KEEP_QC, `TRUE` = REMOVE_QC), labels = c(`FALSE` = "automatic threshold", `TRUE` = "fixed fallback threshold"), name = NULL) +
    guides(fill = guide_legend(override.aes = list(shape = 21), order = 1), shape = guide_legend(order = 2)) +
    scale_x_continuous(limits = c(.5, 1)) + labs(x = "AUROC (simulated > observed)", y = yl) + theme_supp()
}
pA <- scat("threshold_auto", "automatic Scrublet threshold") +
  geom_hline(data = fixed_cut, aes(yintercept = threshold_used, linetype = src), linewidth = .3, colour = "grey30") +
  scale_linetype_manual(values = c("11", "22"), name = "fixed fallback") + theme(legend.position = "none")
pB <- scat("detectable_prop", "detectable doublet fraction\n(proportion)") + theme(legend.position = "right")
hist_df <- dplyr::bind_rows(lapply(SAMPLE_ORDER, function(s) {
  t <- fread(T01(sprintf("%s_cellranger_scrublet_scores.csv.gz", s)))
  dplyr::bind_rows(lapply(c("observed", "simulated"), function(k) { v <- t[[k]][!is.na(t[[k]])]; h <- hist(v, breaks = seq(0, 1, length.out = 31), plot = FALSE)
    data.frame(sample = s, kind = k, x0 = head(h$breaks, -1), x1 = h$breaks[-1], density = h$density) })) }))
lab_df <- sep |> dplyr::filter(source == "cellranger") |>
  dplyr::mutate(strip = sprintf("%s%s\nAUROC %.2f, thr %.2f%s", sample_lab(sample), ifelse(sample %in% FLAGGED, " (flagged)", ""), auroc_sim_gt_obs, threshold_used, ifelse(fixed, " fixed", "")))
hist_df <- hist_df |> dplyr::left_join(dplyr::select(lab_df, sample, strip, threshold_used), by = "sample") |>
  dplyr::mutate(strip = factor(strip, levels = lab_df$strip[match(SAMPLE_ORDER, lab_df$sample)]))
sim_step <- hist_df |> dplyr::filter(kind == "simulated") |> dplyr::group_by(strip) |> dplyr::reframe(x = c(rbind(x0, x1)), y = rep(density, each = 2))
pC <- ggplot() +
  geom_rect(data = dplyr::filter(hist_df, kind == "observed"), aes(xmin = x0, xmax = x1, ymin = 0, ymax = density, fill = "observed nuclei"), alpha = .8) +
  geom_line(data = sim_step, aes(x, y, colour = "simulated doublets"), linewidth = .3) +
  geom_vline(data = dplyr::distinct(hist_df, strip, threshold_used), aes(xintercept = threshold_used, linetype = "threshold used"), colour = REMOVE_QC, linewidth = .3) +
  facet_wrap(~ strip, ncol = 4, scales = "free_y") +
  scale_fill_manual(values = KEEP_QC, name = NULL) + scale_colour_manual(values = "black", name = NULL) + scale_linetype_manual(values = "22", name = NULL) +
  scale_x_continuous(breaks = c(0, .5, 1), labels = c("0", "0.5", "1"), limits = c(0, 1), expand = c(0, 0)) + scale_y_continuous(expand = expansion(mult = c(0, .05))) +
  labs(x = "Scrublet doublet score (raw counts; score used for QC)", y = "density") + theme_supp() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), strip.text = element_text(size = 6, face = "plain", lineheight = .9),
        legend.position = "top", panel.spacing.y = unit(1, "mm"), panel.spacing.x = unit(2.5, "mm"))
kcol <- c(reference = KEEP_QC, depth = "#4C9F70", size = "#E0A03C", input = REMOVE_QC)
ctrl <- ctrl |> dplyr::mutate(control = gsub("CellRanger cell calls", "STARsolo cell calls", control),
                       control = factor(control, levels = rev(control)), kind = factor(kind, levels = names(kcol)),
                       lab = sprintf("thr %.2f, called %.1f %%", threshold, 100 * predicted_prop))
pD <- ggplot(ctrl, aes(auroc_sim_gt_obs, control, fill = kind)) +
  geom_col(width = .7) + geom_vline(xintercept = .5, colour = GREY_QC, linewidth = .3) +
  geom_text(aes(label = lab), hjust = -.05, size = PT6) +
  scale_fill_manual(values = kcol, name = NULL) + coord_cartesian(xlim = c(.5, 1.12), expand = FALSE) +
  scale_x_continuous(breaks = seq(.5, 1, .1)) +
  labs(x = "AUROC (simulated > observed), automatic Scrublet call", y = NULL) + theme_supp() + theme(legend.position = "right")
fig_s2 <- (pA | pB) / wrap(pC) / wrap(pD) + plot_layout(heights = c(48, 122, 45))
figlog_S2 <- supp_save(fig_s2, "Figure_S2", WIDTH_2COL, 228)
sep_raw <- dplyr::filter(sep, source == "cellranger")
LOG_Figure_S2 <- supp_note("Figure_S2", fixed_raw = fixed_cut$threshold_used[fixed_cut$src == "raw counts"], fixed_cb = fixed_cut$threshold_used[fixed_cut$src == "CellBender counts"],
          fixed_libs = FIXED_LIBS, auroc_min = min(sep_raw$auroc_sim_gt_obs), auroc_min_lib = sep_raw$sample[which.min(sep_raw$auroc_sim_gt_obs)],
          auroc_max = max(sep_raw$auroc_sim_gt_obs), auroc_fixed = setNames(sep_raw$auroc_sim_gt_obs[sep_raw$fixed], sep_raw$sample[sep_raw$fixed]),
          ctrl = ctrl |> dplyr::select(control, kind, auroc_sim_gt_obs, threshold, predicted_prop, n))
time_supp_fig_s2 <- as.numeric(difftime(Sys.time(), t0_supp, units = "secs"))

t0_supp <- Sys.time()
qc_all <- fread(T01("cellQC_results.txt.gz"), select = c("V1", "sample", "cellbender_ngenes", "cellbender_ncount", "percent_mito", "exon_prop",
                                                       "cellranger_doublet_scores", "leiden_overcluster", "outlier"))
meta01 <- fread(T01("cell_metadata.csv"), select = c("V1", "sample", "cell_type"))
clean_by_ct <- read.csv(RES("02_subclustering", "atlas_cleaning_by_cell_type.csv"), stringsAsFactors = FALSE)
reason <- read.csv(T01("qc_validation_primary_reason_by_sample.csv"), check.names = FALSE, stringsAsFactors = FALSE)
prof <- read.csv(T01("preqc_cluster_profiles.csv"), stringsAsFactors = FALSE)
pre  <- fread(file.path(PANELS_DIR, "supp_preqc_umap.csv.gz"))
stopifnot(nrow(pre) == nrow(qc_all), nrow(meta01) == sum(clean_by_ct$nuclei), sum(qc_all$outlier) + sum(!qc_all$outlier) == nrow(qc_all))
stopifnot(sum(!qc_all$outlier) == nrow(meta01))
funnel <- data.frame(stage = c("QC input\n(CellBender)", "after QC\n(atlas)", "after atlas\ncleaning"),
                     n = c(nrow(qc_all), nrow(meta01), sum(clean_by_ct$retained)))
time_supp_s3_data <- as.numeric(difftime(Sys.time(), t0_supp, units = "secs"))

t0_supp <- Sys.time()
yl <- function() scale_y_discrete(limits = rev(SAMPLE_ORDER), labels = sample_lab)
ycol <- theme(axis.text.y = element_text(colour = rev(sample_col(SAMPLE_ORDER)), size = 6))
animal_fill <- function() { grp <- split(animals$sample, animals$group)
  unlist(unname(Map(function(a, g) setNames(colorRampPalette(c("grey85", GROUP_FILL[[g]]))(length(a) + 1)[-1], a), grp, names(grp)))) }
ANIMAL_COLS <- animal_fill()
pA <- ggplot(funnel, aes(n, factor(stage, levels = rev(stage)))) +
  geom_col(width = .65, fill = rev(c(GREY_QC, KEEP_QC, "#1F4E79"))[3:1]) +
  geom_text(aes(label = sprintf("%s\n(%.0f %%)", sfmt_n(n), 100 * n / n[1])), hjust = -.08, size = PT6, lineheight = .9) +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale()), limits = c(0, max(funnel$n) * 2.1), breaks = c(0, 1e5), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(x = "nuclei", y = NULL) + theme_supp() + theme(axis.text.y = element_text(size = 6, lineheight = .9))
set.seed(1); pre_o <- pre[sample(nrow(pre))]
tech <- as.character(prof$leiden_overcluster[prof$is_technical %in% c(TRUE, "True")])
pB <- ggplot(pre_o, aes(UMAP_1, UMAP_2, colour = sample)) + scattermore::geom_scattermore(pointsize = 1, pixels = px_for(55)) +
  scale_colour_manual(values = ANIMAL_COLS, breaks = SAMPLE_ORDER, name = NULL) + coord_equal() +
  guides(colour = guide_legend(override.aes = list(pointsize = 4), ncol = 2, keyheight = unit(2.2, "mm"))) +
  labs(x = "UMAP 1", y = "UMAP 2") + theme_supp() + no_axes + theme(legend.text = element_text(size = 6), legend.key.spacing.y = unit(0, "mm"))
pre_o[, tech := as.character(leiden_overcluster) %in% tech]
pC <- ggplot(pre_o[order(tech)], aes(UMAP_1, UMAP_2, colour = tech)) + scattermore::geom_scattermore(pointsize = 1, pixels = px_for(55)) +
  scale_colour_manual(values = c(`FALSE` = GREY_QC, `TRUE` = REMOVE_QC), labels = c(`FALSE` = "other clusters", `TRUE` = "technical clusters"), name = NULL) +
  guides(colour = guide_legend(override.aes = list(pointsize = 4))) + coord_equal() +
  labs(x = "UMAP 1", y = "UMAP 2") + theme_supp() + no_axes + theme(legend.position = "right", legend.direction = "vertical")
reasons <- sub("_prop$", "", setdiff(colnames(reason), "sample"))
rpal <- setNames(c("#0072B2", "#56B4E9", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#F0E442", "#999999")[seq_along(reasons)], reasons)
rl <- reason |> pivot_longer(-sample, names_to = "reason", values_to = "prop") |> dplyr::mutate(reason = factor(sub("_prop$", "", reason), levels = rev(reasons)))
pD <- ggplot(rl, aes(100 * prop, sample, fill = reason)) + geom_col(width = .75) + yl() +
  scale_fill_manual(values = rpal, breaks = reasons, name = NULL) + scale_x_continuous(expand = expansion(mult = c(0, .03))) +
  guides(fill = guide_legend(ncol = 2)) + labs(x = "droplets removed, by primary reason (%)", y = NULL) + theme_supp() + ycol + theme(legend.position = "top")
pE <- ggplot(qc_all, aes(pmax(cellbender_ngenes, 1), sample)) + geom_violin(fill = "grey80", colour = NA, scale = "width", width = .85) +
  stat_summary(fun = median, geom = "point", size = .7, colour = KEEP_QC) +
  geom_vline(xintercept = 150, linetype = "22", linewidth = .3) + yl() +
  scale_x_log10(breaks = c(10, 100, 1000, 10000), labels = label_number(scale_cut = cut_short_scale())) +
  labs(x = "genes per nucleus,\nQC input", y = NULL) + theme_supp() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
comp <- meta01[, .N, by = .(sample, cell_type)][, prop := N / sum(N), by = sample]
comp[, cell_type := factor(cell_type, levels = rev(c(CELL_TYPE_ORDER, "Unassigned")))]
ct_fill <- c(CELL_TYPE_FILL, Unassigned = CELL_TYPE_FILL[["Epicardial"]])
ntot <- meta01[, .N, by = sample]
pF <- ggplot(comp, aes(prop, sample, fill = cell_type)) + geom_col(width = .8) + yl() +
  geom_text(data = ntot, aes(x = 1.02, y = sample, label = sfmt_n(N)), inherit.aes = FALSE, hjust = 0, size = PT6) +
  scale_fill_manual(values = ct_fill, breaks = c(intersect(CELL_TYPE_ORDER, comp$cell_type), "Unassigned"),
                    labels = function(x) ifelse(x == "Unassigned", "Unassigned (Epicardial)", CELL_TYPE_LABEL[x]), name = NULL) +
  scale_x_continuous(breaks = c(0, .5, 1), labels = c("0", "0.5", "1"), expand = c(0, 0)) + coord_cartesian(xlim = c(0, 1.22), clip = "off") +
  labs(x = "proportion of atlas nuclei\n(n at right)", y = NULL) + theme_supp() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.key.size = unit(2.2, "mm"))
post <- qc_all[outlier == FALSE]
mets <- c(cellbender_ngenes = "genes per nucleus", cellbender_ncount = "UMIs per nucleus", percent_mito = "mitochondrial\nUMIs (%)",
          exon_prop = "exonic fraction\n(proportion)", cellranger_doublet_scores = "Scrublet\nscore")
post_l <- melt(post[, c("sample", names(mets)), with = FALSE], id.vars = "sample")[, metric := factor(mets[as.character(variable)], levels = mets)]
post_l[metric %in% mets[1:2], value := log10(pmax(value, 1))]
refs <- post_l[!sample %in% FLAGGED, .(ref = median(value)), by = metric]
pG <- ggplot(post_l, aes(value, sample)) + geom_violin(fill = "grey80", colour = NA, scale = "width", width = .85) +
  stat_summary(fun = median, geom = "point", size = .6, colour = KEEP_QC) +
  geom_vline(data = refs, aes(xintercept = ref), linetype = "11", linewidth = .3) +
  facet_wrap(~ metric, nrow = 1, scales = "free_x", labeller = as_labeller(function(x) ifelse(x %in% mets[1:2], paste0(x, "\n(log10)"), x))) + yl() + scale_x_continuous(n.breaks = 4) +
  labs(x = NULL, y = NULL) + theme_supp() + ycol + theme(strip.text = element_text(size = 6, face = "bold"), panel.spacing.x = unit(2.5, "mm"))
fig_s3 <- wrap_plots(list(pA, pB, pC), nrow = 1, widths = c(.8, 1.3, 1)) / wrap_plots(list(pD, pE, pF), nrow = 1, widths = c(1.5, .55, 1)) / pG +
  plot_layout(heights = c(55, 70, 52))
figlog_S3 <- supp_save(fig_s3, "Figure_S3", WIDTH_2COL, 205)
rem <- read.csv(T01("qc_validation_library_metrics.csv"), stringsAsFactors = FALSE)
LOG_Figure_S3 <- supp_note("Figure_S3", funnel = funnel, n_tech = length(tech), tech = tech, n_overclusters = nrow(prof),
          removed_flagged = setNames(100 * rem$removed_prop[rem$sample %in% FLAGGED], rem$sample[rem$sample %in% FLAGGED]),
          removed_other_range = range(100 * rem$removed_prop[!rem$sample %in% FLAGGED]),
          lowgenes_flagged = setNames(100 * reason[["low genes_prop"]][reason$sample %in% FLAGGED], reason$sample[reason$sample %in% FLAGGED]),
          atlas_n_range = range(ntot$N), atlas_n_min_lib = ntot$sample[which.min(ntot$N)])
time_supp_fig_s3 <- as.numeric(difftime(Sys.time(), t0_supp, units = "secs"))

t0_supp <- Sys.time()
atl  <- fread(file.path(PANELS_DIR, "supp_atlas_obs_umap.csv.gz"))
lisi <- fread(T01("ilisi_per_nucleus.csv.gz")); lisi_s <- read.csv(T01("ilisi_by_sample.csv"))
sexm <- read.csv(T01("sex_check_by_sample.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(atl) == nrow(meta01), setequal(sexm$sample, SAMPLE_ORDER))
time_supp_s4_data <- as.numeric(difftime(Sys.time(), t0_supp, units = "secs"))

t0_supp <- Sys.time()
set.seed(1); atl_o <- atl[sample(nrow(atl))]
umap_cat <- function(col, pal, ncol_leg = 1, title = NULL) ggplot(atl_o, aes(UMAP_1, UMAP_2, colour = .data[[col]])) +
  scattermore::geom_scattermore(pointsize = 1, pixels = px_for(40)) + scale_colour_manual(values = pal, name = title) + coord_equal() +
  guides(colour = guide_legend(override.aes = list(pointsize = 4), ncol = ncol_leg, keyheight = unit(2, "mm"))) +
  labs(x = "UMAP 1", y = "UMAP 2") + theme_supp() + no_axes + theme(legend.position = "bottom", legend.key.spacing.y = unit(0, "mm"))
bat <- sort(unique(atl$batch)); sq <- sort(unique(atl$seq))
pA1 <- umap_cat("sample", ANIMAL_COLS[SAMPLE_ORDER], 4, "animal") + theme(legend.position = "none")
pA2 <- umap_cat("batch", setNames(OKABE_ITO[seq_along(bat)], bat), 4, "batch")
pA3 <- umap_cat("seq", setNames(c("#0072B2", "#E69F00")[seq_along(sq)], sq), 2, "sequencing run")
li <- melt(lisi[, .(pre_harmony, post_harmony)], measure.vars = 1:2)[, variable := factor(variable, labels = c("PCA", "Harmony"))]
n_samp <- length(SAMPLE_ORDER)
pB <- ggplot(li, aes(variable, value, fill = variable)) + geom_boxplot(outlier.shape = NA, width = .6, linewidth = .25) +
  geom_hline(yintercept = n_samp, linetype = "11", linewidth = .3, colour = "grey40") +
  annotate("text", x = .5, y = n_samp - .3, label = sprintf("%d = even mixing", n_samp), hjust = 0, vjust = 1, size = PT6) +
  stat_summary(fun = median, geom = "text", aes(label = fmt_f(after_stat(y), 1)), hjust = -1.4, size = PT6) +
  scale_fill_manual(values = c(GREY_QC, KEEP_QC), guide = "none") + scale_y_continuous(limits = c(0, n_samp + .5)) +
  labs(x = NULL, y = "iLISI (effective number of\nanimals per neighbourhood)") + theme_supp()
qcm <- list(cellbender_ncount = "UMIs per nucleus (log10)", cellbender_ngenes = "genes per nucleus (log10)", percent_mito = "mitochondrial UMIs (%)",
            exon_prop = "exonic fraction (proportion)", cellranger_doublet_scores = "Scrublet score")
qc_umap <- function(col) {
  v <- atl[[col]]; if (grepl("ncount|ngenes", col)) v <- log10(pmax(v, 1))
  lim <- quantile(v, c(.01, .99), na.rm = TRUE); d <- data.frame(UMAP_1 = atl$UMAP_1, UMAP_2 = atl$UMAP_2, v = v)[order(v), ]
  ggplot(d, aes(UMAP_1, UMAP_2, colour = pmin(pmax(v, lim[1]), lim[2]))) + scattermore::geom_scattermore(pointsize = 1, pixels = px_for(34)) +
    scale_colour_viridis_c(name = qcm[[col]], n.breaks = 3, labels = function(x) formatC(x, format = "g", digits = 2)) + coord_equal() +
    guides(colour = guide_colourbar(title.position = "top", barheight = unit(1.6, "mm"), barwidth = unit(18, "mm"))) +
    labs(x = NULL, y = NULL) + theme_supp() + no_axes + theme(legend.position = "bottom", legend.title = element_text(size = 6))
}
pC <- wrap(wrap_plots(lapply(names(qcm), qc_umap), nrow = 1))
sx <- sexm |> pivot_longer(c(Y, Xist), names_to = "gene_set", values_to = "expr") |>
  dplyr::mutate(gene_set = factor(gene_set, levels = c("Y", "Xist"), labels = c("Y genes (Ddx3y, Uty, Eif2s3y, Kdm5d; mean)", "Xist")),
         sex = factor(sex, levels = c("M", "F"), labels = c("recorded male", "recorded female")), sample = factor(sample, levels = rev(SAMPLE_ORDER)))
pD <- ggplot(sx, aes(expr, sample, colour = gene_set, shape = gene_set)) + geom_point(size = 1.3, stroke = .5) +
  facet_grid(sex ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = c(KEEP_QC, REMOVE_QC), name = NULL) + scale_shape_manual(values = c(16, 4), name = NULL) +
  scale_y_discrete(labels = sample_lab) +
  labs(x = "mean log-normalised expression per nucleus", y = NULL) + theme_supp() +
  theme(legend.position = "top", legend.direction = "vertical", strip.text.y = element_text(angle = 0, size = 6))
ls_l <- lisi_s |> pivot_longer(-sample, names_to = "space", values_to = "ilisi") |>
  dplyr::mutate(space = factor(space, levels = c("pre_harmony", "post_harmony"), labels = c("PCA", "Harmony")))
pE <- ggplot(ls_l, aes(ilisi, sample)) + geom_line(aes(group = sample), colour = "grey70", linewidth = .3) +
  geom_point(aes(colour = space), size = 1.2) + scale_y_discrete(limits = rev(SAMPLE_ORDER), labels = sample_lab) +
  scale_colour_manual(values = c(GREY_QC, KEEP_QC), name = NULL) +
  labs(x = "median iLISI of the animal's nuclei", y = NULL) + theme_supp() +
  theme(axis.text.y = element_text(colour = rev(sample_col(SAMPLE_ORDER))), legend.position = "top")
fig_s4 <- wrap(wrap_plots(list(pA1, pA2, pA3), nrow = 1)) + pB + plot_layout(widths = c(3, .75))
fig_s4 <- fig_s4 / pC / (pD | pE) + plot_layout(heights = c(68, 45, 62))
figlog_S4 <- supp_save(fig_s4, "Figure_S4", WIDTH_2COL, 190)
LOG_Figure_S4 <- supp_note("Figure_S4", n_atlas = nrow(atl), n_lisi = nrow(lisi), lisi_pre = median(lisi$pre_harmony), lisi_post = median(lisi$post_harmony),
          lisi_post_range = range(lisi_s$post_harmony), lisi_post_min_lib = lisi_s$sample[which.min(lisi_s$post_harmony)], lisi_pre_min_lib = lisi_s$sample[which.min(lisi_s$pre_harmony)],
          n_male = sum(sexm$sex == "M"), n_female = sum(sexm$sex == "F"),
          y_male = range(sexm$Y[sexm$sex == "M"]), xist_male_max = max(sexm$Xist[sexm$sex == "M"]),
          y_female_max = max(sexm$Y[sexm$sex == "F"]), xist_female = range(sexm$Xist[sexm$sex == "F"]), n_batch = length(bat), n_seq = length(sq))
time_supp_fig_s4 <- as.numeric(difftime(Sys.time(), t0_supp, units = "secs"))

t0_supp <- Sys.time()
DC_PANEL    <- nb_const("02_subclustering.Rmd", "DC_PANEL")
EPI_MARKERS <- nb_const("02_subclustering.Rmd", "EPI_MARKERS")
EPI_PANEL   <- nb_const("02_subclustering.Rmd", "EPI_PANEL")
mk <- fread(T01("broad_cluster_markers.txt.gz"))
mk[cell_type == "Unassigned", cell_type := "Epicardial"]
clean_cl <- read.csv(RES("02_subclustering", "atlas_cleaning_clusters.csv"), stringsAsFactors = FALSE)
dc_cl <- clean_cl |> dplyr::filter(identity == "Dendritic_cell")
stopifnot(nrow(dc_cl) == 1)
dc_mk <- read.csv(RES("02_subclustering", paste0(dc_cl$cell_type, "_subcluster_markers.csv")), stringsAsFactors = FALSE) |>
  dplyr::filter(as.character(cluster) == as.character(dc_cl$cluster)) |> dplyr::arrange(dplyr::desc(myAUC)) |> dplyr::slice_head(n = 3) |>
  dplyr::transmute(names = gene, cell_type = "Dendritic_cell", AUC = myAUC)
topg <- dplyr::bind_rows(as.data.frame(mk)[, c("names", "cell_type", "AUC")], dc_mk) |>
  dplyr::mutate(cell_type = factor(cell_type, levels = CELL_TYPE_ORDER)) |> dplyr::filter(!is.na(cell_type)) |>
  dplyr::arrange(cell_type, dplyr::desc(AUC)) |> dplyr::group_by(cell_type) |> dplyr::slice_head(n = 3) |> dplyr::ungroup()
marker_genes <- unique(topg$names)
evid_genes <- unique(c(DC_PANEL, EPI_PANEL, EPI_MARKERS, "Adgre1", "C1qa", "Dcn"))
so <- readRDS(DATA("objects", "integrated_LV_snRNAseq_seurat_clean.rds"))
genes_all <- intersect(unique(c(marker_genes, evid_genes)), rownames(so))
X <- Seurat::GetAssayData(so, assay = "RNA", layer = "data")[genes_all, ]
ct_clean <- factor(as.character(so$cell_type), levels = CELL_TYPE_ORDER)
grp_clean <- as.character(so$disease)
rm(so); invisible(gc())
stopifnot(!anyNA(ct_clean), all(marker_genes %in% genes_all))
idx <- split(seq_along(ct_clean), droplevels(ct_clean))
mean_e <- sapply(idx, function(i) Matrix::rowMeans(X[, i, drop = FALSE]))
frac_e <- sapply(idx, function(i) Matrix::rowMeans(X[, i, drop = FALSE] > 0))
n_ct <- table(droplevels(ct_clean))
epi_pos2 <- sapply(idx, function(i) mean(Matrix::colSums(X[intersect(EPI_PANEL, genes_all), i, drop = FALSE] > 0) >= 2))
rm(X); invisible(gc())
by_grp <- read.csv(RES("02_subclustering", "atlas_cleaning_by_group.csv"), stringsAsFactors = FALSE)
by_ct  <- read.csv(RES("02_subclustering", "atlas_cleaning_by_cell_type.csv"), stringsAsFactors = FALSE)
stopifnot(sum(n_ct) == sum(by_ct$retained))
time_supp_s5_data <- as.numeric(difftime(Sys.time(), t0_supp, units = "secs"))

t0_supp <- Sys.time()
ct_lv <- names(n_ct)
dot_df <- function(genes) {
  sc <- t(apply(mean_e[genes, , drop = FALSE], 1, function(v) (v - min(v)) / max(max(v) - min(v), 1e-9)))
  data.frame(gene = rep(genes, ncol(sc)), cell_type = rep(colnames(sc), each = length(genes)), scaled = as.vector(sc),
             pct = 100 * as.vector(frac_e[genes, colnames(sc), drop = FALSE])) |>
    dplyr::mutate(gene = factor(gene, levels = genes), cell_type = factor(cell_type, levels = rev(ct_lv)))
}
ylab_ct <- function(x) paste0(CELL_TYPE_LABEL[x], " (", sfmt_n(as.integer(n_ct[x])), ")")
dot_scales <- list(scale_colour_gradient(low = "grey90", high = "#08306b", name = "mean expression\n(scaled per gene)", breaks = c(0, .5, 1)),
                   scale_size(range = c(0, 2.6), name = "nuclei detecting\nthe gene (%)", limits = c(0, 100), breaks = c(25, 50, 75, 100)),
                   scale_y_discrete(labels = ylab_ct))
dA <- dot_df(marker_genes)
grp_x <- topg |> dplyr::distinct(names, .keep_all = TRUE) |> dplyr::mutate(x = match(names, marker_genes)) |> dplyr::group_by(cell_type) |>
  dplyr::summarise(x0 = min(x), x1 = max(x), .groups = "drop")
pA <- ggplot(dA, aes(gene, cell_type)) + geom_point(aes(size = pct, colour = scaled)) +
  annotate("segment", x = grp_x$x0 - .35, xend = grp_x$x1 + .35, y = length(ct_lv) + .75, yend = length(ct_lv) + .75, linewidth = .3, colour = "grey40") +
  annotate("text", x = (grp_x$x0 + grp_x$x1) / 2, y = length(ct_lv) + 1.0, label = CELL_TYPE_SHORT[as.character(grp_x$cell_type)], size = PT6, vjust = 0) +
  dot_scales + coord_cartesian(ylim = c(.5, length(ct_lv) + 1.6), clip = "off") + labs(x = NULL, y = NULL) + theme_supp() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, face = "italic"), legend.position = "right")
bgp <- by_grp |> dplyr::filter(cell_type != "all") |> dplyr::mutate(cell_type = factor(cell_type, levels = CELL_TYPE_ORDER), group = factor(group, levels = GROUP_ORDER))
bcl <- by_ct |> dplyr::mutate(cell_type = factor(cell_type, levels = CELL_TYPE_ORDER), lab = sprintf("%.1f %%\n(%s)", pct_removed, sfmt_n(removed)))
ymax <- max(bgp$removed_pct)
pB <- ggplot(bgp, aes(cell_type, removed_pct, fill = group)) +
  geom_col(position = position_dodge(width = .8), width = .75) +
  geom_text(data = bcl, aes(cell_type, ymax * 1.06, label = lab), inherit.aes = FALSE, size = PT6, lineheight = .85, vjust = 0) +
  scale_fill_manual(values = GROUP_FILL, labels = GROUP_LABEL, name = NULL) + scale_x_discrete(labels = CELL_TYPE_SHORT) +
  scale_y_continuous(expand = expansion(mult = c(0, .02)), limits = c(0, ymax * 1.32)) +
  labs(x = "sub-clustered lineage (notebook 02)", y = "nuclei removed in\natlas cleaning (%)") + theme_supp() + theme(legend.position = "right")
ev_genes <- intersect(c(DC_PANEL, "Adgre1", "C1qa", unique(c(EPI_MARKERS, EPI_PANEL)), "Dcn"), rownames(mean_e))
dC <- dot_df(ev_genes)
epi_lab <- data.frame(cell_type = factor(ct_lv, levels = rev(ct_lv)), lab = sprintf("%.0f %%", 100 * epi_pos2[ct_lv]))
pC <- ggplot(dC, aes(gene, cell_type)) + geom_point(aes(size = pct, colour = scaled)) +
  geom_text(data = epi_lab, aes(x = length(ev_genes) + 1.2, y = cell_type, label = lab), inherit.aes = FALSE, size = PT6, hjust = 0) +
  annotate("text", x = length(ev_genes) + 1.2, y = length(ct_lv) + .9, label = "\u2265 2 epicardial\npanel genes", size = PT6, hjust = 0, vjust = 0, lineheight = .85) +
  annotate("segment", x = c(.65, length(DC_PANEL) + 2.65), xend = c(length(DC_PANEL) + .35, length(DC_PANEL) + 2 + length(unique(c(EPI_MARKERS, EPI_PANEL))) + .35),
           y = length(ct_lv) + .75, yend = length(ct_lv) + .75, linewidth = .3, colour = "grey40") +
  annotate("text", x = c((1 + length(DC_PANEL)) / 2, length(DC_PANEL) + 2 + (1 + length(unique(c(EPI_MARKERS, EPI_PANEL)))) / 2), y = length(ct_lv) + 1,
           label = c("dendritic-cell panel", "epicardial markers"), size = PT6, vjust = 0) +
  dot_scales + coord_cartesian(xlim = c(.5, length(ev_genes) + 3.2), ylim = c(.5, length(ct_lv) + 2), clip = "off") +
  labs(x = NULL, y = NULL) + theme_supp() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, face = "italic"), legend.position = "none")
fig_s5 <- pA / pB / pC + plot_layout(heights = c(62, 38, 62))
figlog_S5 <- supp_save(fig_s5, "Figure_S5", WIDTH_2COL, 200)
det <- function(g, ct) 100 * frac_e[g, ct]
LOG_Figure_S5 <- supp_note("Figure_S5", n_clean = sum(n_ct), n_ct = n_ct, n_marker_genes = length(marker_genes), dc_cluster = dc_cl, dc_markers = dc_mk$names,
          removed_total = sum(by_ct$removed), nuclei_total = sum(by_ct$nuclei), by_ct = by_ct,
          epi_det = sapply(EPI_MARKERS, det, ct = "Epicardial"), epi_det_other_max = sapply(EPI_MARKERS, function(g) max(frac_e[g, setdiff(ct_lv, "Epicardial")]) * 100),
          epi_pos2 = 100 * epi_pos2, dc_det = sapply(intersect(DC_PANEL, rownames(frac_e)), det, ct = "Dendritic_cell"),
          dc_det_mac = sapply(intersect(DC_PANEL, rownames(frac_e)), det, ct = "Macrophage"))
time_supp_fig_s5 <- as.numeric(difftime(Sys.time(), t0_supp, units = "secs"))

t0_supp <- Sys.time()
ifn <- fread(file.path(PANELS_DIR, "supp_ifn_score.csv.gz"))
chk <- read.csv(T01("single_animal_cluster_check.csv"), stringsAsFactors = FALSE)
cab <- read.csv(T01("cluster_annotation_before_exclusions.csv"), stringsAsFactors = FALSE)
sa  <- cab |> dplyr::filter(cell_type == "Cardiomyocyte", largest_sample_fraction > .5)
FOCUS <- names(which.max(table(sa$largest_sample))); SINGLE <- as.character(sa$leiden[sa$largest_sample == FOCUS])
SINGLE_OTHER <- sa[sa$largest_sample != FOCUS, c("leiden", "largest_sample", "largest_sample_fraction")]
stopifnot(FOCUS %in% NOTED)
support <- read.csv(RES("07b_ec_cm_subtype_ccc", "tables", "subtype_animal_support.csv"), stringsAsFactors = FALSE, check.names = FALSE)
cmo <- atl[cell_type == "Cardiomyocyte"][, leiden := as.character(leiden)]
cm_cl <- as.character(chk$leiden)
cmo <- cmo[leiden %in% cm_cl][, leiden := factor(leiden, levels = cm_cl)]
lowc <- 0.5 * median(atl[cell_type == "Cardiomyocyte", cellbender_ngenes])
foc_col <- GROUP_FILL[[as.character(animals$group[animals$sample == FOCUS])]]
ccol <- setNames(ifelse(cm_cl %in% SINGLE, foc_col, ifelse(cm_cl %in% as.character(SINGLE_OTHER$leiden), "grey45", GREY_QC)), cm_cl)
sh <- cmo[, .(focus = mean(individual == FOCUS), n = .N), by = leiden]
pA <- ggplot(sh, aes(leiden, focus)) + geom_col(aes(y = 1), fill = GREY_QC, width = .75) + geom_col(fill = foc_col, width = .75) +
  geom_hline(yintercept = 1 / 16, linetype = "11", linewidth = .3) +
  scale_x_discrete(labels = function(x) paste0(x, "\n", sfmt_n(sh$n[match(x, sh$leiden)]))) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0), breaks = c(0, 1 / 16, .5, 1), labels = c("0", "1/16", "0.5", "1")) +
  labs(x = "cardiomyocyte cluster (notebook 01)\nand nuclei", y = sprintf("proportion of nuclei\nfrom %s", FOCUS)) + theme_supp()
qm <- melt(cmo[, .(leiden, cellbender_ngenes, cellranger_doublet_scores)], id.vars = "leiden")[, variable := factor(variable, labels = c("genes per nucleus", "Scrublet score"))]
pB <- ggplot(qm, aes(leiden, value, fill = leiden)) + geom_boxplot(outlier.shape = NA, width = .7, linewidth = .25) +
  geom_hline(data = data.frame(variable = factor("genes per nucleus", levels = levels(qm$variable)), y = lowc), aes(yintercept = y), linetype = "11", linewidth = .3) +
  facet_wrap(~ variable, scales = "free_y", nrow = 1) + scale_fill_manual(values = ccol, guide = "none") +
  labs(x = "cardiomyocyte cluster", y = NULL) + theme_supp()
sc_ifn <- merge(atl[cell_type %in% c("Cardiomyocyte", "Endothelial"), .(barcode, individual, cell_type, leiden)], ifn, by = "barcode")
stopifnot(nrow(sc_ifn) == nrow(ifn))
an <- sc_ifn[, .(score = mean(score_IFN)), by = .(individual, cell_type)][, group := animals$group[match(individual, animals$sample)]]
ord <- an[cell_type == "Cardiomyocyte"][order(score), individual]
pC <- ggplot(an, aes(factor(individual, levels = ord), score, shape = cell_type, fill = group)) +
  geom_point(size = 1.6, stroke = .3, colour = "grey25", position = position_dodge(width = .5)) +
  scale_shape_manual(values = c(Cardiomyocyte = 21, Endothelial = 22), labels = c("cardiomyocytes", "endothelial cells"), name = NULL) +
  scale_fill_manual(values = GROUP_FILL, labels = GROUP_LABEL, name = NULL) +
  guides(fill = guide_legend(override.aes = list(shape = 21), nrow = 2), shape = guide_legend(nrow = 2)) +
  scale_x_discrete(labels = sample_lab) + labs(x = NULL, y = "mean IFN-\u03b3 GTPase\nscore per animal") + theme_supp() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, colour = sample_col(ord), face = ifelse(ord == FOCUS, "bold", "plain")),
        legend.position = "top", legend.box = "horizontal")
b5 <- sc_ifn[individual == FOCUS & cell_type == "Cardiomyocyte"]
b5[, set := ifelse(leiden %in% SINGLE, paste("cluster", leiden), "other clusters")]
b5_lv <- c(paste("cluster", SINGLE), "other clusters"); b5[, set := factor(set, levels = b5_lv)]
ref_med <- median(sc_ifn[cell_type == "Cardiomyocyte" & individual != FOCUS, score_IFN])
pD <- ggplot(b5, aes(set, score_IFN, fill = set != "other clusters")) + geom_violin(colour = NA, scale = "width", width = .8) +
  stat_summary(fun = median, geom = "point", size = 1, colour = "black") +
  geom_hline(yintercept = ref_med, linetype = "11", linewidth = .3) +
  scale_fill_manual(values = c(`TRUE` = foc_col, `FALSE` = GREY_QC), guide = "none") +
  scale_x_discrete(labels = function(x) paste0(x, "\n(n = ", sfmt_n(table(b5$set)[x]), ")")) +
  labs(x = NULL, y = sprintf("IFN-\u03b3 GTPase score,\n%s cardiomyocyte nuclei", FOCUS)) + theme_supp()
sup <- support |> dplyr::mutate(lineage = ifelse(grepl(" EC$", subtype), "EC subtypes", "CM subtypes"), subtype = factor(subtype, levels = rev(subtype)),
                         tested = ifelse(group_tested %in% c(TRUE, "TRUE"), "tested", "one_animal"))
pE <- ggplot(sup, aes(100 * top_animal_share_prop, subtype, fill = tested)) + geom_col(width = .7) +
  geom_text(aes(label = sample_lab(top_animal)), hjust = -.1, size = PT6, colour = ifelse(sup$top_animal %in% FLAGGED, FLAG_RED, "black")) +
  geom_vline(xintercept = 50, linetype = "22", linewidth = .3) + facet_grid(lineage ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c(tested = KEEP_QC, one_animal = GREY_QC), labels = c(tested = "group-tested", one_animal = "one animal \u2265 50 %"), name = NULL) +
  scale_x_continuous(limits = c(0, 115), breaks = c(0, 25, 50, 75, 100), expand = c(0, 0)) +
  labs(x = "nuclei from the animal contributing most (%)", y = NULL) + theme_supp() +
  theme(legend.position = "top", strip.text.y = element_text(angle = 0, size = 6))
fig_s6 <- (wrap(pA) + wrap(pB) + plot_layout(widths = c(1.25, 1))) / (wrap(pC) + wrap(pD) + plot_layout(widths = c(1.4, 1))) /
  (wrap(pE) + plot_spacer() + plot_layout(widths = c(1.6, 1))) + plot_layout(heights = c(50, 52, 50))
figlog_S6 <- supp_save(fig_s6, "Figure_S6", WIDTH_2COL, 175)
LOG_Figure_S6 <- supp_note("Figure_S6", focus = FOCUS, single = SINGLE, single_other = SINGLE_OTHER, share_single = setNames(sh$focus[match(SINGLE, sh$leiden)], SINGLE),
          n_single = setNames(sh$n[match(SINGLE, sh$leiden)], SINGLE), lowc = lowc, an = an, ref_med = ref_med,
          b5_n = table(b5$set), b5_med = tapply(b5$score_IFN, b5$set, median), support = support, n_cm_clusters = length(cm_cl))
time_supp_fig_s6 <- as.numeric(difftime(Sys.time(), t0_supp, units = "secs"))

t0_supp <- Sys.time()
suppressPackageStartupMessages(library(Seurat))
run_sum  <- read.csv(RES("02_subclustering", "run_summary.csv"), stringsAsFactors = FALSE)
LINEAGES <- intersect(CELL_TYPE_ORDER, run_sum$cell_type[run_sum$status == "processed"])
removed_bc <- fread(RES("02_subclustering", "atlas_cleaning_removed_barcodes.csv"))
clean_cl2  <- read.csv(RES("02_subclustering", "atlas_cleaning_clusters.csv"), stringsAsFactors = FALSE)
GM10800_FAMILY <- nb_const("04_Cardiomyocyte_fine_annotation.Rmd", "GM10800_FAMILY")
DOT_PANELS <- list(Endothelial = nb_const("03_Endothelial_fine_annotation.Rmd", "dot_panels"),
                   Cardiomyocyte = nb_const("04_Cardiomyocyte_fine_annotation.Rmd", "dot_panels"))
lineage <- list()
for (ct in LINEAGES) {
  so <- readRDS(RES("02_subclustering", paste0(ct, "_subcluster.rds")))
  emb <- Embeddings(so, "umap")[, 1:2]
  d <- data.table(barcode = colnames(so), cluster = as.character(Idents(so)), individual = as.character(so$individual),
                  UMAP_1 = emb[, 1], UMAP_2 = emb[, 2])
  out <- list(cells = d, levels = levels(Idents(so)))
  if (ct %in% names(DOT_PANELS)) {
    s <- so@misc$subclustering
    out$misc <- s[c("resolution", "rule_resolution", "fallback", "n_clusters", "n_pcs_computed", "n_pcs_used", "dropped_pcs", "pc_exon_prop_cor", "sweep")]
    panels <- lapply(DOT_PANELS[[ct]], intersect, y = rownames(so))
    out$dot <- DotPlot(so, features = panels)$data
    rm_cl <- as.character(clean_cl2$cluster[clean_cl2$cell_type == ct & clean_cl2$removed %in% c(TRUE, "TRUE")])
    out$ward <- ward_merge(so, "seurat_clusters", VariableFeatures(so), exclude = rm_cl)
    out$removed_clusters <- rm_cl
  }
  lineage[[ct]] <- out
  rm(so); invisible(gc())
}
time_supp_lineage_data <- as.numeric(difftime(Sys.time(), t0_supp, units = "secs"))

t0_supp <- Sys.time()
CYTO_PC_COR <- nb_const("02_subclustering.Rmd", "CYTO_PC_COR")
cluster_palette <- function(lv) setNames(colorRampPalette(OKABE_ITO)(max(length(lv), length(OKABE_ITO)))[seq_along(lv)], lv)
cl_labels <- function(d) { bw <- diff(range(d$UMAP_1)) / 40
  d[, .(bx = round(UMAP_1 / bw), by = round(UMAP_2 / bw)), by = cluster][, .N, by = .(cluster, bx, by)][order(-N)][, .SD[1], by = cluster][, `:=`(UMAP_1 = bx * bw, UMAP_2 = by * bw)] }
umap_cl <- function(d, pal, w_mm = 40, psize = 1.2) ggplot(d, aes(UMAP_1, UMAP_2, colour = cluster)) +
  scattermore::geom_scattermore(pointsize = psize, pixels = px_for(w_mm)) +
  geom_text_repel(data = cl_labels(d), aes(label = cluster), colour = "black", size = PT6, seed = 1, box.padding = .1, min.segment.length = Inf) +
  scale_colour_manual(values = pal, guide = "none") + labs(x = "UMAP 1", y = "UMAP 2") + theme_supp() + no_axes
diag_row <- function(ct) {
  L <- lineage[[ct]]; s <- L$misc; d <- copy(L$cells)[, cluster := factor(cluster, levels = L$levels)]; pal <- cluster_palette(L$levels)
  p_umap <- umap_cl(d, pal) + coord_equal()
  comp <- d[, .N, by = .(cluster, individual)][, individual := factor(individual, levels = SAMPLE_ORDER)]
  p_comp <- ggplot(comp, aes(cluster, N, fill = individual)) + geom_col(position = "fill", width = .85) +
    scale_fill_manual(values = ANIMAL_COLS, breaks = SAMPLE_ORDER, name = NULL) + scale_y_continuous(expand = c(0, 0), breaks = c(0, .5, 1), labels = c("0", "0.5", "1")) +
    labs(x = "sub-cluster", y = "proportion of nuclei") + theme_supp() + theme(legend.position = "none")
  pcs <- seq_len(s$n_pcs_computed)
  p_pc <- ggplot(data.frame(pc = pcs, r = unname(s$pc_exon_prop_cor), dropped = pcs %in% s$dropped_pcs), aes(pc, r, fill = dropped)) +
    geom_col(width = .8) + geom_hline(yintercept = c(-1, 1) * CYTO_PC_COR, linetype = "22", colour = "grey40", linewidth = .25) +
    geom_hline(yintercept = 0, colour = "grey30", linewidth = .25) +
    scale_fill_manual(values = c(`FALSE` = "grey60", `TRUE` = OKABE_ITO[6]), labels = c(`FALSE` = "PC used", `TRUE` = "PC dropped"), name = NULL) +
    scale_x_continuous(breaks = c(1, seq(5, max(pcs), 5))) + labs(x = "principal component", y = "Pearson r with\nexonic fraction") + theme_supp() +
    theme(legend.position = "bottom")
  sw <- as.data.frame(s$sweep) |> dplyr::select(resolution, n_clusters, n_markerless, eligible) |> pivot_longer(c(n_clusters, n_markerless), names_to = "metric")
  sel <- data.frame(x = c(s$rule_resolution, s$resolution), what = c("rule", "used"))
  p_sw <- ggplot(sw, aes(resolution, value, colour = metric)) +
    geom_vline(data = sel[!duplicated(sel$x) | sel$what == "used", ], aes(xintercept = x, linetype = what), linewidth = .3, colour = "grey30") +
    geom_line(linewidth = .4) + geom_point(aes(shape = eligible), size = 1.3) +
    scale_colour_manual(values = c(n_clusters = OKABE_ITO[5], n_markerless = OKABE_ITO[6]), labels = c(n_clusters = "clusters", n_markerless = "without marker"), name = NULL) +
    scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), labels = c(`FALSE` = "ineligible", `TRUE` = "eligible"), name = NULL) +
    scale_linetype_manual(values = c(rule = "11", used = "22"), labels = c(rule = "rule resolution", used = "resolution used"), name = NULL) +
    scale_x_continuous(breaks = s$sweep$resolution) + scale_y_continuous(limits = c(0, NA)) +
    guides(colour = guide_legend(order = 1, ncol = 1), shape = guide_legend(order = 2, ncol = 1), linetype = guide_legend(order = 3, ncol = 1)) +
    labs(x = "Louvain resolution", y = "number of clusters") + theme_supp() + theme(legend.position = "right")
  list(p_umap, p_comp, p_pc, p_sw)
}
rows <- lapply(c("Endothelial", "Cardiomyocyte"), diag_row)
fig_s7 <- wrap_plots(c(rows[[1]], rows[[2]]), nrow = 2, widths = c(1.25, 1, .95, 1.3))
figlog_S7 <- supp_save(fig_s7, "Figure_S7", WIDTH_2COL, 130)
LOG_Figure_S7 <- supp_note("Figure_S7", misc = lapply(lineage[c("Endothelial", "Cardiomyocyte")], `[[`, "misc"), n = sapply(lineage[c("Endothelial", "Cardiomyocyte")], function(L) nrow(L$cells)),
          cyto_pc_cor = CYTO_PC_COR)
time_supp_fig_s7 <- as.numeric(difftime(Sys.time(), t0_supp, units = "secs"))

t0_supp <- Sys.time()
by_ct8 <- read.csv(RES("02_subclustering", "atlas_cleaning_by_cell_type.csv"), stringsAsFactors = FALSE)
relab <- clean_cl2 |> dplyr::filter(decision == "relabelled")
s8_rows <- list()
panel8 <- function(ct) {
  L <- lineage[[ct]]; d <- copy(L$cells)
  keep <- d[!barcode %in% removed_bc$barcode[removed_bc$cell_type == ct]][, cluster := factor(cluster, levels = L$levels)]
  gone <- sort(as.integer(unique(removed_bc$cluster[removed_bc$cell_type == ct])))
  stopifnot(nrow(keep) == by_ct8$retained[by_ct8$cell_type == ct], nrow(d) == by_ct8$nuclei[by_ct8$cell_type == ct])
  rl <- relab[relab$cell_type == ct, ]
  s8_rows[[ct]] <<- data.frame(cell_type = ct, nuclei = nrow(d), retained = nrow(keep), removed_clusters = paste(gone, collapse = ", "),
                               relabelled = if (nrow(rl)) sprintf("%s (cluster %s)", rl$identity, rl$cluster) else "")
  strip <- sprintf("%s\n%s of %s retained\nremoved: %s", CELL_TYPE_LABEL[[ct]], sfmt_n(nrow(keep)), sfmt_n(nrow(d)), if (length(gone)) paste(gone, collapse = ", ") else "none")
  lab <- cl_labels(keep)
  if (nrow(rl)) lab[cluster %in% rl$cluster, cluster := paste0(cluster, " (", CELL_TYPE_SHORT[rl$identity], ")")]
  keep$strip <- strip
  ggplot(keep, aes(UMAP_1, UMAP_2, colour = cluster)) +
    scattermore::geom_scattermore(pointsize = if (nrow(keep) > 1e4) 1.2 else 3, pixels = px_for(42)) +
    geom_text_repel(data = lab, aes(label = cluster), colour = "black", size = PT6, seed = 1, box.padding = .1, min.segment.length = Inf) +
    facet_wrap(~ strip) + scale_colour_manual(values = cluster_palette(L$levels), guide = "none") +
    labs(x = NULL, y = NULL) + theme_supp() + no_axes + theme(strip.text = element_text(size = 6, face = "plain", lineheight = .95), strip.clip = "off")
}
p8 <- lapply(LINEAGES, panel8)
fig_s8 <- wrap(wrap_plots(p8, ncol = 4))
figlog_S8 <- supp_save(fig_s8, "Figure_S8", WIDTH_2COL, 165)
LOG_Figure_S8 <- supp_note("Figure_S8", table = dplyr::bind_rows(s8_rows), n_lineages = length(LINEAGES))
time_supp_fig_s8 <- as.numeric(difftime(Sys.time(), t0_supp, units = "secs"))

hc_segments <- function(hc) {
  n <- length(hc$order); xl <- setNames(match(seq_len(n), hc$order), seq_len(n)); x <- numeric(n - 1); y <- hc$height; seg <- list()
  pos <- function(k) if (k < 0) xl[[as.character(-k)]] else x[k]
  hgt <- function(k) if (k < 0) 0 else y[k]
  for (i in seq_len(n - 1)) { a <- hc$merge[i, 1]; b <- hc$merge[i, 2]; x[i] <- (pos(a) + pos(b)) / 2
    seg[[i]] <- data.frame(x = c(pos(a), pos(b), pos(a)), xend = c(pos(a), pos(b), pos(b)), y = c(hgt(a), hgt(b), y[i]), yend = c(y[i], y[i], y[i])) }
  list(seg = dplyr::bind_rows(seg), leaves = data.frame(x = seq_len(n), label = hc$labels[hc$order]))
}
annotation_figure <- function(ct, prefix, qc_file, module_file, ward_file, prop_file, subtype_col, fig_name, heights, height_mm) {
  L <- lineage[[ct]]
  cq <- read.csv(RES(qc_file), stringsAsFactors = FALSE) |> dplyr::mutate(cluster = as.character(cluster))
  subtypes <- unique(read.csv(RES(prop_file), stringsAsFactors = FALSE)[[subtype_col]])
  stopifnot(setequal(subtypes, unique(cq$label[cq$kind == "subtype"])))
  st_col <- subtype_palette(subtypes)
  parent <- vapply(cq$label, function(l) { p <- subtypes[startsWith(l, subtypes)]; if (length(p)) p[1] else NA_character_ }, "")
  col_of <- setNames(ifelse(cq$kind == "subtype", st_col[cq$label], ifelse(cq$kind == "technical" & !is.na(parent), colorspace::lighten(st_col[parent], .35, space = "HLS"), "grey50")), cq$cluster)
  col_of[is.na(col_of)] <- "grey50"
  txt_col <- setNames(colorspace::darken(col_of, .25), names(col_of))
  cl_lab <- setNames(sprintf("%s  %s", cq$cluster, cq$label), cq$cluster)

  ms <- read.csv(RES(module_file), check.names = FALSE, stringsAsFactors = FALSE)
  ev <- ms |> pivot_longer(-module, names_to = "cluster", values_to = "score") |> dplyr::mutate(cluster = sub("^c", "", cluster),
          kind = factor(dplyr::case_when(grepl("_contam$", module) ~ "contamination", grepl("_tech$", module) ~ "technical", TRUE ~ "signature"), levels = c("signature", "technical", "contamination")),
          module = factor(module, levels = rev(ms$module)), cluster = factor(cluster, levels = cq$cluster))
  cap <- 1.2
  pA <- ggplot(ev, aes(cluster, module, fill = score)) + geom_tile(colour = "white", linewidth = .3) +
    geom_text(aes(label = sprintf("%.2f", round(score, 2) + 0), colour = abs(score) > cap / 2), size = PT6) +
    facet_grid(kind ~ ., scales = "free_y", space = "free_y") +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-cap, cap), oob = squish, name = sprintf("mean module score\n(capped at \u00b1%.1f)", cap)) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey20"), guide = "none") +
    scale_x_discrete(labels = function(x) str_wrap(cl_lab[x], 18)) + scale_y_discrete(labels = function(x) gsub("_", " ", x)) +
    labs(x = NULL, y = NULL) + theme_supp() +
    theme(axis.text.x = element_text(colour = txt_col[cq$cluster], lineheight = .85, angle = 35, hjust = 1), axis.line = element_blank(), axis.ticks = element_blank(),
          strip.text.y = element_text(angle = 0, size = 6), legend.key.height = unit(5, "mm"))

  dd <- L$dot |> dplyr::mutate(id = factor(as.character(id), levels = rev(cq$cluster)), feature.groups = factor(feature.groups, levels = names(DOT_PANELS[[ct]]), labels = sub("^Atrial$", "Atrial genes", names(DOT_PANELS[[ct]]))))
  pB <- ggplot(dd, aes(features.plot, id, size = pct.exp, colour = avg.exp.scaled)) + geom_point() +
    facet_grid(. ~ feature.groups, scales = "free_x", space = "free_x") +
    scale_colour_gradient(low = "grey88", high = "firebrick", name = "mean expression\n(scaled)") +
    scale_size(range = c(0, 2.4), name = "nuclei detecting\nthe gene (%)", limits = c(0, 100)) +
    scale_y_discrete(labels = function(x) cl_lab[x]) + labs(x = NULL, y = NULL) + theme_supp() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, face = "italic"), axis.text.y = element_text(colour = rev(txt_col[cq$cluster])),
          strip.text.x = element_text(size = 6, angle = 90, hjust = 0), strip.clip = "off", panel.spacing.x = unit(.6, "mm"), legend.position = "bottom")

  wm <- L$ward; saved <- read.csv(RES(ward_file), stringsAsFactors = FALSE)
  stopifnot(all(unname(wm$merge_of[as.character(saved$subcluster)]) == as.character(saved$merged_cluster)))
  hs <- hc_segments(wm$tree)
  pC <- ggplot(hs$seg) + geom_segment(aes(x = x, xend = xend, y = y, yend = yend), linewidth = .3, colour = "grey30") +
    geom_hline(yintercept = wm$cut, linetype = "22", colour = FLAG_RED, linewidth = .3) +
    annotate("text", x = .6, y = wm$cut, label = sprintf("cut = %.2f x largest distance", WARD_MERGE_FRACTION), hjust = 0, vjust = -.4, size = PT6, colour = FLAG_RED) +
    scale_x_continuous(breaks = hs$leaves$x, labels = str_wrap(cl_lab[hs$leaves$label], 16), expand = expansion(add = .5)) +
    scale_y_continuous(expand = expansion(mult = c(0, .05))) +
    labs(x = NULL, y = "Ward.D2 height\n(Euclidean, mean log-expression)") + theme_supp() +
    theme(axis.text.x = element_text(colour = txt_col[hs$leaves$label], angle = 35, hjust = 1, lineheight = .85))
  fig <- pA / wrap(pB) / pC + plot_layout(heights = heights)
  font <- supp_save(fig, fig_name, WIDTH_2COL, height_mm)
  c(list(font = font), supp_note(fig_name, cluster_qc = cq, n_nuclei = nrow(L$cells), ward_cut = wm$cut, ward_max = wm$max_dist, excluded = L$removed_clusters,
            merged = length(unique(wm$merge_of)) < length(wm$merge_of), dot_genes = sum(lengths(DOT_PANELS[[ct]]))))
}

t0_supp <- Sys.time()
LOG_Figure_S9 <- annotation_figure("Endothelial", "ec", "03_ec_annotation/tables/ec_cluster_qc.csv", "03_ec_annotation/tables/ec_module_scores_by_cluster.csv",
                  "03_ec_annotation/tables/ec_ward_merge.csv", "03_ec_annotation/tables/ec_subtype_proportions_per_sample.csv", "EC_subtype",
                  "Figure_S9", heights = c(80, 62, 45), height_mm = 215)
figlog_S9 <- LOG_Figure_S9$font
time_supp_fig_s9 <- as.numeric(difftime(Sys.time(), t0_supp, units = "secs"))

t0_supp <- Sys.time()
LOG_Figure_S10 <- annotation_figure("Cardiomyocyte", "cm", "04_cm_annotation/tables/cm_cluster_qc.csv", "04_cm_annotation/tables/cm_module_scores_by_cluster.csv",
                  "04_cm_annotation/tables/cm_ward_merge.csv", "04_cm_annotation/tables/cm_subtype_proportions_per_sample.csv", "CM_subtype",
                  "Figure_S10", heights = c(90, 62, 42), height_mm = 222)
figlog_S10 <- LOG_Figure_S10$font
time_supp_fig_s10 <- as.numeric(difftime(Sys.time(), t0_supp, units = "secs"))

t0_supp <- Sys.time()

suppressPackageStartupMessages({ library(Seurat); library(DESeq2); library(edgeR); library(limma) })
PB <- function(f) RES("05b_ec_pseudobulk", f)
MIN_PB_NUCLEI <- nb_const("05b_EC_pseudobulk_sensitivity.Rmd", "MIN_PB_NUCLEI")
pb_de      <- fread(PB("ec_pseudobulk_DE.csv"))
pb_samples <- read.csv(PB("ec_pseudobulk_samples.csv"), stringsAsFactors = FALSE)
ec_obj <- readRDS(RES("03_ec_annotation", "Endothelial_subcluster_annotated_clean.rds")); DefaultAssay(ec_obj) <- "RNA"
PB_SUBTYPES <- levels(droplevels(factor(ec_obj$EC_subtype)))
ec_obj$pb_key <- paste(ec_obj$sample, ec_obj$EC_subtype, sep = "|")
pb_counts <- AggregateExpression(ec_obj, assays = "RNA", group.by = "pb_key")$RNA
key_map <- ec_obj@meta.data |> dplyr::count(sample, EC_subtype, pb_key) |> dplyr::mutate(col = gsub("_", "-", pb_key))
colnames(pb_counts) <- key_map$pb_key[match(colnames(pb_counts), key_map$col)]
rm(ec_obj); invisible(gc())
mds <- dplyr::bind_rows(lapply(PB_SUBTYPES, function(st) {
  all_s <- pb_samples |> dplyr::filter(celltype == st); fit_s <- dplyr::filter(all_s, in_model %in% c(TRUE, "TRUE"))
  cnt <- round(as.matrix(pb_counts[, fit_s$pb_key]))
  dm <- model.matrix(DESIGN_FORMULA, data = fit_s); colnames(dm) <- sub("^group", "", colnames(dm))
  cnt <- cnt[filterByExpr(DGEList(cnt), group = fit_s$group), ]
  dds <- DESeq(DESeqDataSetFromMatrix(cnt, colData = data.frame(fit_s, row.names = fit_s$pb_key), design = dm), quiet = TRUE)
  saved <- pb_de[celltype == st & contrast == "OC_vs_YC"]
  cm <- make_contrast_matrix(dm)[resultsNames(dds), ]
  chk <- as.data.frame(results(dds, contrast = cm[, "OC_vs_YC"]))
  stopifnot(nrow(dds) == unique(saved$n_genes_tested), isTRUE(all.equal(chk$log2FoldChange, saved$logFC[match(rownames(dds), saved$gene)], tolerance = 1e-6)))
  dds_all <- estimateSizeFactors(DESeqDataSetFromMatrix(round(as.matrix(pb_counts[rownames(dds), all_s$pb_key])), colData = data.frame(all_s, row.names = all_s$pb_key), design = ~ 1))
  dispersionFunction(dds_all) <- dispersionFunction(dds)
  m <- plotMDS(assay(varianceStabilizingTransformation(dds_all, blind = FALSE)), plot = FALSE)
  all_s |> dplyr::mutate(dim1 = m$x, dim2 = m$y)
}))
time_supp_s11_mds <- as.numeric(difftime(Sys.time(), t0_supp, units = "secs"))

t0_supp <- Sys.time()
TARGET_GENES <- nb_const("05_EC_DE_analysis.Rmd", "TARGET_GENES")
FDR_CUT <- nb_const("05b_EC_pseudobulk_sensitivity.Rmd", "FDR_CUT"); FDR_LOOSE <- nb_const("05b_EC_pseudobulk_sensitivity.Rmd", "FDR_LOOSE")
shared <- fread(PB("ec_pseudobulk_vs_mast_genes.csv")); conc <- read.csv(PB("ec_pseudobulk_concordance.csv"), stringsAsFactors = FALSE)
fct <- function(d) d[, `:=`(celltype = factor(celltype, levels = PB_SUBTYPES), contrast = factor(contrast, levels = CONTRAST_ORDER))]
fct(pb_de); fct(shared)
CAP_CT <- PB_SUBTYPES[1]
pt_theme <- theme(legend.position = "bottom", panel.spacing.x = unit(2, "mm"))
mds_p <- mds |> dplyr::mutate(celltype = factor(celltype, levels = PB_SUBTYPES), group = factor(group, levels = GROUP_ORDER),
                       flagged = sample %in% FLAGGED, in_model = in_model %in% c(TRUE, "TRUE"))
pA <- ggplot(mds_p, aes(dim1, dim2, fill = group, shape = flagged, alpha = in_model)) +
  geom_point(size = 1.5, colour = "grey25", stroke = .3) +
  geom_text_repel(data = function(d) dplyr::filter(d, flagged | sample %in% NOTED), aes(label = sample_lab(sample)), size = PT6, colour = "grey25",
                  max.overlaps = Inf, min.segment.length = 0, segment.colour = "grey70", segment.size = .15, seed = 1, box.padding = .12, show.legend = FALSE) +
  facet_wrap(~ celltype, nrow = 1, scales = "free", labeller = as_labeller(short_subtype)) +
  scale_fill_manual(values = GROUP_FILL, labels = GROUP_LABEL, breaks = GROUP_ORDER, name = NULL) +
  scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 24), labels = c(`FALSE` = "QC pass", `TRUE` = "QC flag"), name = NULL) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = .35), labels = c(`TRUE` = "in model", `FALSE` = sprintf("< %d nuclei, not in model", MIN_PB_NUCLEI)), name = NULL) +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 2)), shape = guide_legend(override.aes = list(size = 2))) +
  scale_y_continuous(expand = expansion(mult = .15), n.breaks = 3) + scale_x_continuous(expand = expansion(mult = .12), n.breaks = 3) +
  labs(x = "MDS 1 (leading log2FC)", y = "MDS 2") + theme_supp() + pt_theme
CLASS_COL <- c("not a MAST hit" = "grey80", "MAST hit" = "#0072B2", "MAST hit, single-animal driven" = "#E69F00")
sc_ <- shared[celltype == CAP_CT][, class := factor(fcase(!reported, "not a MAST hit", single_animal_driven %in% TRUE, "MAST hit, single-animal driven", default = "MAST hit"), levels = names(CLASS_COL))][order(class)]
pB <- ggplot(sc_, aes(mast_stat, stat, colour = class)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = .25) + geom_vline(xintercept = 0, colour = "grey70", linewidth = .25) +
  scattermore::geom_scattermore(pointsize = 3, alpha = .8, pixels = px_for(34)) +
  geom_point(data = sc_[hit %in% TRUE], shape = 21, size = 1.4, colour = "black", fill = NA, stroke = .35) +
  facet_wrap(~ contrast, nrow = 1, labeller = as_labeller(CONTRAST_LABEL_WRAP)) +
  scale_colour_manual(values = CLASS_COL, name = NULL) + guides(colour = guide_legend(override.aes = list(pointsize = 4))) +
  labs(x = "MAST signed \u2212log10 p (per nucleus)", y = "pseudobulk Wald\nstatistic (per animal)") + theme_supp() + pt_theme
pb_sum <- pb_de[, .(FDR05 = sum(FDR < FDR_CUT, na.rm = TRUE), FDR10 = sum(FDR < FDR_LOOSE, na.rm = TRUE)), by = .(celltype, contrast)]
hits_l <- melt(pb_sum, id.vars = c("celltype", "contrast"))[, variable := factor(variable, levels = c("FDR10", "FDR05"), labels = sprintf("FDR < %s", c(FDR_LOOSE, FDR_CUT)))]
pC <- ggplot(hits_l, aes(contrast, value, fill = variable)) + geom_col(position = position_dodge(width = .75), width = .7) +
  facet_wrap(~ celltype, nrow = 1, labeller = as_labeller(short_subtype)) + scale_fill_manual(values = c("#56B4E9", "#D55E00"), name = NULL) +
  scale_x_discrete(labels = CONTRAST_CODE) + scale_y_continuous(expand = expansion(mult = c(0, .08))) +
  labs(x = NULL, y = "genes (pseudobulk)") + theme_supp() + pt_theme
ma <- pb_de[celltype == CAP_CT][, class := factor(fcase(FDR < FDR_CUT, sprintf("FDR < %s", FDR_CUT), FDR < FDR_LOOSE, sprintf("FDR < %s", FDR_LOOSE), default = "not significant"),
                                                  levels = c("not significant", sprintf("FDR < %s", FDR_LOOSE), sprintf("FDR < %s", FDR_CUT)))][order(class)]
pD <- ggplot(ma, aes(baseMean, logFC, colour = class)) + geom_hline(yintercept = 0, colour = "grey60", linewidth = .25) +
  scattermore::geom_scattermore(data = ma[class == "not significant"], pointsize = 3, alpha = .8, pixels = px_for(34)) +
  geom_point(data = ma[class != "not significant"], size = .8) +
  geom_text_repel(data = ma[FDR < FDR_CUT][order(FDR)][, head(.SD, 4), by = contrast], aes(label = gene), size = PT6, fontface = "italic",
                  max.overlaps = 30, show.legend = FALSE, segment.size = .15, min.segment.length = 0, seed = 1) +
  facet_wrap(~ contrast, nrow = 1, labeller = as_labeller(CONTRAST_LABEL_WRAP)) + scale_x_log10(labels = label_number(scale_cut = cut_short_scale())) +
  scale_colour_manual(values = setNames(c("grey80", "#56B4E9", "#D55E00"), levels(ma$class)), name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  labs(x = "mean normalised count (baseMean)", y = "log2FC (pseudobulk)") + theme_supp() + pt_theme
tg <- shared[celltype == CAP_CT & gene %in% TARGET_GENES][, `:=`(gene = factor(gene, levels = TARGET_GENES), is_int = contrast == "Interaction")]
CO_FILL <- CONTRAST_FILL
lim <- max(abs(c(tg$avg_log2FC, tg$logFC)), na.rm = TRUE) * 1.1
pE <- ggplot(tg, aes(avg_log2FC, logFC)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey75", linewidth = .25, linetype = "22") +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = .25) + geom_vline(xintercept = 0, colour = "grey70", linewidth = .25) +
  geom_point(aes(fill = contrast, shape = reported), size = 1.7, colour = "grey20", stroke = .3) +
  facet_wrap(~ gene, nrow = 1) + scale_fill_manual(values = CO_FILL, labels = CONTRAST_LABEL, name = NULL) +
  scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 22), labels = c(`TRUE` = "reported MAST hit", `FALSE` = "not reported"), name = NULL) +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 2), nrow = 1), shape = guide_legend(override.aes = list(size = 2))) +
  coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  labs(x = "MAST effect (per nucleus): log2FC; Age x HFpEF:\ninteraction effect (MAST, log2 scale)", y = "pseudobulk log2FC") +
  theme_supp() + pt_theme + theme(strip.text = element_text(face = "bold.italic"), legend.box = "vertical", legend.spacing.y = unit(0, "mm"))
fig_s11 <- pA / pB / pC / pD / pE + plot_layout(heights = c(1, 1, .75, 1, 1.05))
figlog_S11 <- supp_save(fig_s11, "Figure_S11", WIDTH_2COL, 228)
LOG_Figure_S11 <- supp_note("Figure_S11", n_profiles = nrow(pb_samples), n_in_model = sum(pb_samples$in_model %in% c(TRUE, "TRUE")), below = pb_samples |> dplyr::filter(!in_model %in% c(TRUE, "TRUE")) |> dplyr::select(sample, celltype, n_nuclei),
          min_nuclei = MIN_PB_NUCLEI, n_genes = pb_de[, .(n = n_genes_tested[1]), by = celltype], pb_sum = pb_sum, conc = conc, targets = tg,
          fdr_cut = FDR_CUT, fdr_loose = FDR_LOOSE, cap_ct = CAP_CT, target_genes = TARGET_GENES,
          pb_target_hits = shared[gene %in% TARGET_GENES & hit %in% TRUE, .(celltype, contrast, gene, logFC, FDR)],
          pb_all = pb_de[, .(n_cut = sum(FDR < FDR_CUT, na.rm = TRUE), n_loose = sum(FDR < FDR_LOOSE, na.rm = TRUE), n_tests = .N)])
time_supp_fig_s11 <- as.numeric(difftime(Sys.time(), t0_supp, units = "secs"))

t0_supp <- Sys.time()
gsea <- fread(RES("06_ec_gsea", "ec_gsea_mast_results.csv"))
GSEA_FDR <- 0.05; NES_LIMITS <- c(-2, 2); SIZE_CAP <- 10; JACCARD_COLLAPSE <- 0.5
gsea[, `:=`(sig = is.finite(padj) & padj < GSEA_FDR, lfdr = -log10(pmax(padj, 1e-300)), le = strsplit(leadingEdge, ";", fixed = TRUE))]
GSEA_SUBTYPES <- intersect(PB_SUBTYPES, unique(gsea$celltype))
stopifnot(setequal(unique(gsea$celltype), GSEA_SUBTYPES), setequal(unique(gsea$contrast), CONTRAST_ORDER))
collapse_paths <- function(d) {
  best <- d[order(padj, -abs(NES)), .SD[1], by = pathway][order(padj, -abs(NES))]; kept <- list(); rep_of <- character()
  for (i in seq_len(nrow(best))) { le <- best$le[[i]]
    jac <- vapply(kept, function(k) length(intersect(le, k)) / length(union(le, k)), numeric(1))
    if (!length(kept) || max(jac) <= JACCARD_COLLAPSE) { kept[[best$pathway[i]]] <- le; rep_of[best$pathway[i]] <- best$pathway[i] } else rep_of[best$pathway[i]] <- names(kept)[which.max(jac)] }
  data.table(pathway = names(rep_of), representative = unname(rep_of))
}
time_supp_gsea_data <- as.numeric(difftime(Sys.time(), t0_supp, units = "secs"))

ACRONYMS <- c("ATP","RNA","DNA","NMD","SRP","MHC","I","II","III","IV","IGF","TNF","TNFA","IFN","IL","VEGF","PDGF","TGF","TGFB","EGF","EGFR","FGF","FGFR","GPCR",
              "NADH","NAD","ROS","ER","UV","HIV","EIF2AK4","GCN2","EIF2","ECM","LDL","HDL","PI3K","AKT","MTOR","MTORC1","MAPK","ERK","JAK","STAT","STAT3","STAT5",
              "E2F","G2M","KRAS","MYC","P53","TP53","WNT","TCA","UPR","SNARE","RHO","ABC","SLC","HSP","HSF1","NOTCH","TLR","NOD","CD4","CD8","NK","ADP","PPAR","PPARA",
              "ATF4","ATF6","XBP1","HDAC","FC","MAP","SUMO","ROBO","SLIT","PLC","PKA","PKC","CDK","G1","G2","NCAM1","MET","L1CAM","RAF","RAS","TCR","BCR","IRE1","PERK",
              "NOX","NOS","LPS","OXPHOS","PTEN","MITF","CDH1","GTP","GDP","AMP","CAMP","CGMP","RUNX1","RUNX2","RUNX3","SMAD","YAP1","TAZ","HIF","NFE2L2","KEAP1",
              "MIRO","RHOBTB3","SLC2A4","GLUT4","PD","L1","CD274","MLL3","MLL4","BRAF","RAF1","ALK","NGF","IFNA","IFNB","EPHB","CD28","MAP2K","IRF3","CD40","IGF1R","PDGFR","VEGFR","VEGFR2","EPH","NRAGE","P75NTR","NTRK1","BMP","FCERI","FCGR","DAP12","SUMOYLATION")
SPECIAL <- c(PPARALPHA = "PPARalpha", MRNA = "mRNA", MIRNA = "miRNA", RRNA = "rRNA", TRNA = "tRNA", NCRNA = "ncRNA", LNCRNA = "lncRNA", NFKB = "NF-kB", GTPASE = "GTPase", GTPASES = "GTPases")
pretty_lab <- function(x, db, width) {
  x <- gsub("_", " ", sub("^(HALLMARK|REACTOME|GOBP)_", "", x))
  if (db != "Hallmark") x <- vapply(strsplit(sub("^(\\w)", "\\U\\1", tolower(x), perl = TRUE), " "), function(v) { u <- toupper(v)
    paste(ifelse(u %in% names(SPECIAL), SPECIAL[u], ifelse(u %in% ACRONYMS, u, v)), collapse = " ") }, character(1))
  str_wrap(x, width)
}
gsea_figure <- function(db, max_n, fig_name, height_mm, wrap_w) {
  d <- gsea[geneset == db]; ds <- d[sig == TRUE]
  cl <- collapse_paths(ds)
  ranked <- ds[pathway %in% cl[pathway == representative]$pathway, .(n = .N, best = max(lfdr)), by = pathway][order(-n, -best)]
  paths <- head(ranked$pathway, max_n)
  saved <- fread(RES("06_ec_gsea", "ec_gsea_collapsed_sets.csv"))[geneset == db]
  collapse_ok <- setequal(paste(cl[pathway != representative]$pathway, cl[pathway != representative]$representative), paste(saved$pathway, saved$representative))
  grid <- CJ(pathway = paths, celltype = GSEA_SUBTYPES, contrast = CONTRAST_ORDER)
  pd <- merge(grid, d, by = c("pathway", "celltype", "contrast"), all.x = TRUE)
  labs_v <- make.unique(pretty_lab(paths, db, wrap_w))
  pd[, `:=`(pathway_lab = factor(labs_v[match(pathway, paths)], levels = rev(labs_v)), contrast = factor(contrast, levels = CONTRAST_ORDER),
            celltype = factor(celltype, levels = GSEA_SUBTYPES), tested = !is.na(NES), fill_nes = ifelse(technical_candidate %in% TRUE, NA_real_, NES))]
  p <- ggplot(pd[tested == TRUE], aes(celltype, pathway_lab)) +
    geom_point(aes(fill = fill_nes, size = pmin(lfdr, SIZE_CAP), alpha = sig, colour = sig), shape = 21, stroke = .3) +
    geom_point(data = pd[tested == TRUE & single_animal_candidate %in% TRUE], size = .45, colour = "black", shape = 16, show.legend = FALSE) +
    geom_text(data = pd[tested == FALSE], label = "\u2013", colour = "grey55", size = BASE_PT / .pt) +
    scale_fill_gradientn(colours = c("#2166AC", "#6587C0", "#F7F7F7", "#DB7E78", "#B2182B"), limits = NES_LIMITS, oob = squish, na.value = TECH_GREY, name = "NES") +
    scale_size_continuous(range = c(.5, 3.4), breaks = c(1, 2, 5, 10), limits = c(0, SIZE_CAP), labels = c("1", "2", "5", "\u2265 10"), name = "\u2212log10 FDR") +
    scale_alpha_manual(values = c(`FALSE` = .45, `TRUE` = 1), guide = "none") + scale_colour_manual(values = c(`FALSE` = "#BBBBBB", `TRUE` = "#333333"), guide = "none") +
    facet_wrap(~ contrast, nrow = 1, drop = FALSE, labeller = as_labeller(CONTRAST_LABEL_WRAP)) +
    scale_x_discrete(drop = FALSE, labels = short_subtype) + scale_y_discrete(drop = FALSE) + labs(x = NULL, y = NULL) +
    guides(fill = guide_colourbar(order = 1, barwidth = unit(18, "mm"), barheight = unit(2.2, "mm"), title.vjust = .9),
           size = guide_legend(order = 2, override.aes = list(fill = "grey60", colour = "#333333"))) +
    theme_supp() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), axis.text.y = element_text(lineheight = .85),
                         panel.grid.major = element_line(colour = "grey94", linewidth = .2), panel.spacing.x = unit(1, "mm"),
                         legend.position = "bottom", legend.box = "horizontal")
  font <- supp_save(p, fig_name, WIDTH_2COL, height_mm)
  c(list(font = font), supp_note(fig_name, db = db, n_sets_tested = uniqueN(d$pathway), n_sig_sets = uniqueN(ds$pathway), n_sig_cells = nrow(ds),
            n_sig_single = sum(ds$single_animal_candidate %in% TRUE), n_sig_tech = sum(ds$technical_candidate %in% TRUE), n_collapsed = nrow(cl[pathway != representative]),
            n_shown = length(paths), max_n = max_n, collapse_matches_saved = collapse_ok, n_not_tested_cells = sum(!pd$tested),
            universe = range(d$universe_n, na.rm = TRUE), by_contrast = ds[, .N, by = contrast]))
}

t0_supp <- Sys.time()
LOG_Figure_S12 <- gsea_figure("Hallmark", 60, "Figure_S12", 150, 45)
figlog_S12 <- LOG_Figure_S12$font
time_supp_fig_s12 <- as.numeric(difftime(Sys.time(), t0_supp, units = "secs"))

t0_supp <- Sys.time()
LOG_Figure_S13 <- gsea_figure("Reactome", 60, "Figure_S13", 228, 90)
figlog_S13 <- LOG_Figure_S13$font
time_supp_fig_s13 <- as.numeric(difftime(Sys.time(), t0_supp, units = "secs"))

t0_supp <- Sys.time()
LOG_Figure_S14 <- gsea_figure("GO:BP", 45, "Figure_S14", 222, 62)
figlog_S14 <- LOG_Figure_S14$font
time_supp_fig_s14 <- as.numeric(difftime(Sys.time(), t0_supp, units = "secs"))

t0_supp <- Sys.time()
suppressPackageStartupMessages(library(CellChat))
cc_list <- readRDS(RES("07_ccc", "cellchat_list.rds"))[GROUP_ORDER]
CC_TYPES <- intersect(CELL_TYPE_ORDER, unique(unlist(lapply(cc_list, function(x) levels(x@idents)))))
CC_PVAL <- nb_const("07_CCC_analysis.Rmd", "CC_PVAL"); TOP_PATHWAYS <- nb_const("07_CCC_analysis.Rmd", "TOP_PATHWAYS")
cc_weight <- lapply(cc_list, function(x) x@net$weight); cc_count <- lapply(cc_list, function(x) x@net$count)
FOCUS_CT <- c("Endothelial", "Cardiomyocyte")
ecm_lr <- dplyr::bind_rows(lapply(GROUP_ORDER, function(g) { x <- cc_list[[g]]; p <- x@net$prob; pv <- x@net$pval
  dplyr::bind_rows(lapply(list(FOCUS_CT, rev(FOCUS_CT)), function(d) data.frame(group = g, source = d[1], target = d[2], interaction_name = dimnames(p)[[3]],
                                                                         prob = p[d[1], d[2], ], pval = pv[d[1], d[2], ]))) })) |> dplyr::filter(prob > 0)
lr_pathway <- dplyr::bind_rows(lapply(cc_list, function(x) x@LR$LRsig[, c("interaction_name", "pathway_name")])) |> dplyr::distinct()
rm(cc_list); invisible(gc())
time_supp_ccc_data <- as.numeric(difftime(Sys.time(), t0_supp, units = "secs"))

t0_supp <- Sys.time()
role <- dplyr::bind_rows(lapply(GROUP_ORDER, function(g) { w <- cc_weight[[g]]; data.frame(group = g, cell_type = rownames(w), outgoing = rowSums(w), incoming = colSums(w)) })) |>
  dplyr::mutate(group = factor(group, levels = GROUP_ORDER), cell_type = factor(cell_type, levels = CC_TYPES)) |>
  pivot_longer(c(outgoing, incoming), names_to = "direction", values_to = "strength") |> dplyr::mutate(direction = factor(direction, levels = c("outgoing", "incoming")))
pA <- ggplot(role, aes(group, strength, fill = group)) + geom_col(width = .75) +
  facet_grid(direction ~ cell_type, scales = "free_y", labeller = labeller(cell_type = CELL_TYPE_SHORT)) +
  scale_fill_manual(values = GROUP_FILL, labels = GROUP_LABEL, name = NULL) + scale_y_continuous(expand = expansion(mult = c(0, .08)), n.breaks = 3) +
  labs(x = NULL, y = "summed communication\nprobability") + theme_supp() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom", panel.spacing = unit(1.2, "mm"), strip.text.y = element_text(size = 6))
flow <- read.csv(RES("07_ccc", "tables", "cellchat_pathway_flow.csv"), stringsAsFactors = FALSE)
flow_panel <- function(ctr) {
  f <- flow |> dplyr::filter(group %in% c(ctr$i1, ctr$i2)) |> dplyr::select(pathway, group, flow, rel_flow_prop) |>
    pivot_wider(names_from = group, values_from = c(flow, rel_flow_prop), values_fill = 0)
  f$delta <- f[[paste0("rel_flow_prop_", ctr$i1)]] - f[[paste0("rel_flow_prop_", ctr$i2)]]
  f <- f |> dplyr::slice_max(abs(delta), n = TOP_PATHWAYS, with_ties = FALSE) |>
    dplyr::mutate(share1 = .data[[paste0("flow_", ctr$i1)]] / (.data[[paste0("flow_", ctr$i1)]] + .data[[paste0("flow_", ctr$i2)]]))
  f$pathway <- factor(f$pathway, levels = f$pathway[order(f$share1)])
  fl <- dplyr::bind_rows(data.frame(pathway = f$pathway, group = ctr$i1, share = f$share1), data.frame(pathway = f$pathway, group = ctr$i2, share = 1 - f$share1)) |>
    dplyr::mutate(group = factor(group, levels = c(ctr$i2, ctr$i1)))
  ggplot(fl, aes(share, pathway, fill = group)) + geom_col(width = .8) + geom_vline(xintercept = .5, colour = "white", linewidth = .3) +
    scale_fill_manual(values = GROUP_FILL, guide = "none") + scale_x_continuous(breaks = c(0, .5, 1), labels = c("0", "0.5", "1"), expand = c(0, 0)) +
    labs(x = sprintf("%s: %s vs %s\n(proportion of flow)", CONTRAST_LABEL[[ctr$name]], ctr$i1, ctr$i2), y = NULL) + theme_supp() +
    theme(axis.title.x = element_text(size = 6))
}
pB <- wrap(wrap_plots(lapply(CONTRASTS, flow_panel), nrow = 1))
fig_s15 <- pA / pB + plot_layout(heights = c(50, 120))
figlog_S15 <- supp_save(fig_s15, "Figure_S15", WIDTH_2COL, 185)
LOG_Figure_S15 <- supp_note("Figure_S15", cc_types = CC_TYPES, n_pathways = tapply(flow$pathway, flow$group, function(x) length(unique(x))), n_pathways_union = length(unique(flow$pathway)),
          top = TOP_PATHWAYS, cc_pval = CC_PVAL, total_count = sapply(cc_count, sum), total_weight = sapply(cc_weight, sum))
time_supp_fig_s15 <- as.numeric(difftime(Sys.time(), t0_supp, units = "secs"))

t0_supp <- Sys.time()
dec <- fread(RES("07_ccc", "tables", "ccc_diff_edge_decomposition.csv"))
N_EDGES <- 6; N_PAIRS <- 8
top_edges <- unique(dec[, .(contrast, source, target, edge_delta_total)])[order(-abs(edge_delta_total))][, head(.SD, N_EDGES), by = contrast]
bars <- dec[top_edges, on = .(contrast, source, target, edge_delta_total)][order(-abs(delta))][, head(.SD, N_PAIRS), by = .(contrast, source, target)]
bars[, edge := sprintf("%s\u2192%s  \u0394 %+.3f", CELL_TYPE_SHORT[source], CELL_TYPE_SHORT[target], edge_delta_total)]
bars[, row := paste(contrast, source, target, interaction_name, sep = "||")]
dec_plot <- function(ctr) {
  d <- bars[contrast == ctr$name]
  d[, row := factor(row, levels = row[order(delta)])]; d[, edge := factor(edge, levels = unique(edge[order(-abs(edge_delta_total))]))]
  ggplot(d, aes(delta, row, fill = delta > 0)) + geom_col(width = .7) + geom_vline(xintercept = 0, colour = "grey50", linewidth = .25) +
    facet_wrap(~ edge, nrow = 2, scales = "free") + scale_fill_manual(values = c(`TRUE` = "#b2182b", `FALSE` = "#2166ac"), guide = "none") +
    scale_y_discrete(labels = function(x) sub("^.*\\|\\|", "", x)) + scale_x_continuous(n.breaks = 3) +
    labs(x = sprintf("\u0394 communication probability, %s \u2212 %s (%s)", ctr$i1, ctr$i2, CONTRAST_LABEL[[ctr$name]]), y = NULL) + theme_supp() +
    theme(axis.title.x = element_text(size = BASE_PT), strip.text = element_text(size = 6, face = "bold"), panel.spacing.x = unit(1.5, "mm"), panel.spacing.y = unit(.8, "mm"))
}
fig_s16 <- wrap(wrap_plots(lapply(CONTRASTS, dec_plot), ncol = 1))
figlog_S16 <- supp_save(fig_s16, "Figure_S16", WIDTH_2COL, 228)

strength_total <- unique(rbind(dec[, .(group = group1, source, target, pathway, interaction_name, prob = prob_group1)],
                               dec[, .(group = group2, source, target, pathway, interaction_name, prob = prob_group2)]))[, .(total = sum(prob)), by = group]
strength_total <- setNames(strength_total$total, strength_total$group)[GROUP_ORDER]
LOG_Figure_S16 <- supp_note("Figure_S16", n_rows = nrow(dec), n_edges = N_EDGES, n_pairs = N_PAIRS, top_edges = top_edges, strength_total = strength_total)
time_supp_fig_s16 <- as.numeric(difftime(Sys.time(), t0_supp, units = "secs"))

t0_supp <- Sys.time()
MAX_LR_SHOWN <- 30
inferred <- ecm_lr |> dplyr::group_by(source, target, interaction_name) |> dplyr::filter(any(pval < CC_PVAL)) |> dplyr::ungroup()
n_inferred <- inferred |> dplyr::distinct(source, target, interaction_name) |> dplyr::count(source, target, name = "n_lr")
shown <- inferred |> dplyr::group_by(source, target, interaction_name) |> dplyr::summarise(maxp = max(prob[pval < CC_PVAL]), .groups = "drop") |>
  dplyr::group_by(source, target) |> dplyr::slice_max(maxp, n = MAX_LR_SHOWN, with_ties = FALSE) |> dplyr::ungroup()
bub <- inferred |> dplyr::semi_join(shown, by = c("source", "target", "interaction_name")) |>
  dplyr::mutate(dir = sprintf("%s \u2192 %s", CELL_TYPE_SHORT[source], CELL_TYPE_SHORT[target]), group = factor(group, levels = GROUP_ORDER),
         pcat = cut(pval, c(-Inf, .01, .05, Inf), labels = c("lt01", "lt05", "ge05")))
bub$lr <- factor(bub$interaction_name, levels = rev(unique(shown$interaction_name[order(shown$source, -shown$maxp)])))
pA <- ggplot(bub, aes(group, lr)) + geom_point(aes(colour = prob, size = pcat)) +
  facet_wrap(~ dir, nrow = 1, scales = "free_y") +
  scale_colour_gradientn(colours = c("#fee8c8", "#fc8d59", "#b30000"), name = "communication\nprobability", n.breaks = 3) +
  scale_size_manual(values = c(lt01 = 2.4, lt05 = 1.6, ge05 = .7), labels = c(lt01 = "p < 0.01", lt05 = "0.01 \u2264 p < 0.05", ge05 = "p \u2265 0.05"), name = "permutation p", drop = FALSE) +
  scale_x_discrete(labels = GROUP_ORDER) + scale_y_discrete(labels = function(x) gsub("_", " ", x)) +
  labs(x = NULL, y = NULL) + theme_supp() + theme(axis.text.y = element_text(size = 6), panel.grid.major = element_line(colour = "grey94", linewidth = .2),
                                                  axis.text.x = element_text(colour = colorspace::darken(GROUP_FILL[GROUP_ORDER], .3), face = "bold"))
st <- read.csv(RES("07_ccc", "tables", "lr_group_overview.csv"), stringsAsFactors = FALSE) |>
  dplyr::mutate(source = factor(source, levels = CC_TYPES), target = factor(target, levels = CC_TYPES), disease = factor(disease, levels = GROUP_ORDER))
npairs <- st |> dplyr::distinct(source, target, n_pairs)
tile_axes <- list(scale_y_discrete(limits = rev(CC_TYPES), labels = CELL_TYPE_SHORT), scale_x_discrete(limits = CC_TYPES, labels = CELL_TYPE_SHORT),
                  labs(x = "receiver", y = "sender"), coord_fixed(), theme_supp(),
                  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5), axis.line = element_blank(), axis.ticks = element_blank()))
p_n <- ggplot(npairs, aes(target, source)) + geom_tile(fill = "grey95", colour = "white", linewidth = .3) +
  geom_text(aes(label = n_pairs), size = PT6, colour = "grey20") + tile_axes
p_s <- ggplot(st, aes(target, source, fill = mean_score)) + geom_tile(colour = "white", linewidth = .3) +
  facet_wrap(~ disease, nrow = 2, labeller = as_labeller(GROUP_LABEL)) +
  scale_fill_gradient(low = "#f7fbff", high = "#08306b", name = "mean score per\neligible LR pair", na.value = "grey85") + tile_axes +
  theme(legend.key.height = unit(5, "mm"))
fig_s17 <- pA / (wrap(p_n) + wrap(p_s) + plot_layout(widths = c(1, 1.15))) + plot_layout(heights = c(112, 100))
figlog_S17 <- supp_save(fig_s17, "Figure_S17", WIDTH_2COL, 225)
LOG_Figure_S17 <- supp_note("Figure_S17", n_inferred = n_inferred, max_shown = MAX_LR_SHOWN, npairs_range = range(npairs$n_pairs), n_animals_range = range(st$n_animals))
time_supp_fig_s17 <- as.numeric(difftime(Sys.time(), t0_supp, units = "secs"))

t0_supp <- Sys.time()
edges <- read.csv(RES("07b_ec_cm_subtype_ccc", "tables", "cellchat_ec_cm_subtype_edges.csv"), stringsAsFactors = FALSE)
sup18 <- read.csv(RES("07b_ec_cm_subtype_ccc", "tables", "subtype_animal_support.csv"), stringsAsFactors = FALSE, check.names = FALSE)
ov18 <- read.csv(RES("07b_ec_cm_subtype_ccc", "tables", "lr_subtype_group_overview.csv"), stringsAsFactors = FALSE)
ONE_ANIMAL_TYPES <- sup18$subtype[!sup18$group_tested %in% c(TRUE, "TRUE")]
ST_ORDER <- sup18$subtype; EC_ST <- ST_ORDER[grepl(" EC$", ST_ORDER)]; CM_ST <- ST_ORDER[grepl(" CM$", ST_ORDER)]
circle_order <- c(EC_ST, rev(CM_ST))
node <- data.frame(subtype = circle_order, ang = pi / 2 - 2 * pi * (seq_along(circle_order) - 1) / length(circle_order)) |>
  dplyr::mutate(x = cos(ang), y = sin(ang), lineage = ifelse(subtype %in% EC_ST, "EC", "CM"), one = subtype %in% ONE_ANIMAL_TYPES,
         lab = sub("Interferon-response", "IFN-response", sub(" (EC|CM)$", "", subtype)))
node$fill <- ifelse(node$one, "grey75", ifelse(node$lineage == "EC", CELL_TYPE_FILL[["Endothelial"]], CELL_TYPE_FILL[["Cardiomyocyte"]]))
W <- edges |> dplyr::select(source, target, group, weight) |> pivot_wider(names_from = group, values_from = weight, values_fill = 0)
diffs <- c(lapply(CONTRASTS, function(ctr) list(name = ctr$name, d = W[[ctr$i1]] - W[[ctr$i2]])),
           list(list(name = "Interaction", d = (W$OH - W$OC) - (W$YH - W$YC))))
R_NODE <- .13
circle_plot <- function(dd) {
  e <- W |> dplyr::select(source, target) |> dplyr::mutate(d = dd$d) |> dplyr::filter(d != 0) |>
    dplyr::left_join(node |> dplyr::select(source = subtype, x0 = x, y0 = y), by = "source") |> dplyr::left_join(node |> dplyr::select(target = subtype, x1 = x, y1 = y), by = "target") |>
    dplyr::mutate(len = sqrt((x1 - x0)^2 + (y1 - y0)^2), xs = x0 + (x1 - x0) * R_NODE / len, ys = y0 + (y1 - y0) * R_NODE / len,
           xe = x1 - (x1 - x0) * R_NODE / len, ye = y1 - (y1 - y0) * R_NODE / len, one = source %in% ONE_ANIMAL_TYPES | target %in% ONE_ANIMAL_TYPES,
           col = ifelse(one, "one-animal subtype", ifelse(d > 0, "higher in first group", "higher in second group")), lw = .15 + 1.6 * abs(d) / max(abs(d)))
  ggplot() + geom_curve(data = e, aes(x = xs, y = ys, xend = xe, yend = ye, colour = col, linewidth = lw), curvature = .2, alpha = .75,
                        arrow = arrow(length = unit(1.1, "mm"), type = "closed")) +
    geom_point(data = node, aes(x, y), shape = 21, size = 4.2, fill = node$fill, colour = "white", stroke = .3) +
    geom_text(data = node, aes(1.32 * x, 1.2 * y, label = lab, hjust = ifelse(abs(x) < .2, .5, ifelse(x > 0, 0, 1))), size = PT6) +
    scale_colour_manual(values = c(`higher in first group` = "#b2182b", `higher in second group` = "#2166ac", `one-animal subtype` = "grey70"), name = NULL, drop = FALSE) +
    scale_linewidth_identity() + coord_equal(xlim = c(-2.1, 2.1), ylim = c(-1.35, 1.35), clip = "off") +
    labs(title = NULL, x = NULL, y = NULL, subtitle = NULL) +
    annotate("text", x = 0, y = 1.62, label = if (dd$name == "Interaction") "Age x HFpEF: (OH \u2212 OC) \u2212 (YH \u2212 YC)" else
      sprintf("%s: %s \u2212 %s", CONTRAST_LABEL[[dd$name]], CONTRASTS[[match(dd$name, CONTRAST_ORDER)]]$i1, CONTRASTS[[match(dd$name, CONTRAST_ORDER)]]$i2), size = PT7, fontface = "bold") +
    theme_void(base_size = BASE_PT, base_family = FONT) + theme(legend.position = "none", plot.margin = margin(4, 1, 1, 1, "mm"))
}
circ <- lapply(diffs, circle_plot)
leg_df <- data.frame(x = 1:3, y = 1, col = c("higher in first group", "higher in second group", "one-animal subtype"))
circ_leg <- ggplot(leg_df, aes(x, y, colour = col)) + geom_point(alpha = 0) +
  scale_colour_manual(values = c(`higher in first group` = "#b2182b", `higher in second group` = "#2166ac", `one-animal subtype` = "grey70"), name = NULL) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, shape = 15, size = 3), ncol = 1)) + theme_void(base_size = BASE_PT, base_family = FONT) +
  theme(legend.text = element_text(size = 6), legend.position = "inside", legend.position.inside = c(.5, .5))
pA <- wrap(wrap_plots(c(circ, list(circ_leg)), nrow = 2))
ed <- edges |> dplyr::mutate(edge = sprintf("%s \u2192 %s", sub(" (EC|CM)$", "", source), sub(" (EC|CM)$", "", target)), one = source %in% ONE_ANIMAL_TYPES | target %in% ONE_ANIMAL_TYPES,
                      group = factor(group, levels = GROUP_ORDER))
edge_lv <- ed |> dplyr::distinct(direction, source, target, edge) |> dplyr::arrange(direction, match(source, ST_ORDER), match(target, ST_ORDER)) |> dplyr::pull(edge)
ed$edge <- factor(ed$edge, levels = rev(edge_lv))
one_edges <- unique(as.character(ed$edge[ed$one]))
pB <- ggplot(ed, aes(group, edge, fill = weight)) + geom_tile(colour = "white", linewidth = .3) +
  geom_text(aes(label = sprintf("%.2f", weight)), size = PT6, colour = ifelse(ed$weight > .6 * max(ed$weight), "white", "grey20")) +
  facet_grid(direction ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient(low = "#f7fbff", high = "#08306b", name = "summed\ncommunication\nprobability") +
  labs(x = NULL, y = NULL) + theme_supp() +
  theme(axis.text.y = element_text(colour = ifelse(levels(ed$edge) %in% one_edges, "grey55", "black")), axis.line = element_blank(), axis.ticks = element_blank(),
        strip.text.y = element_text(angle = 0, size = 6), legend.key.height = unit(5, "mm"))
ov <- ov18 |> dplyr::mutate(edge = factor(sprintf("%s \u2192 %s", sub(" (EC|CM)$", "", source), sub(" (EC|CM)$", "", target)), levels = levels(ed$edge)),
                     disease = factor(disease, levels = GROUP_ORDER))
pC <- ggplot(ov, aes(disease, edge, fill = mean_score)) + geom_tile(colour = "white", linewidth = .3) +
  geom_text(aes(label = n_animals), size = PT6, colour = "grey30") +
  facet_grid(direction ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient(low = "#fff7ec", high = "#b30000", name = "mean animal-level\nscore per eligible\nLR pair", na.value = "grey88") +
  labs(x = NULL, y = NULL) + theme_supp() +
  theme(axis.text.y = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(), strip.text.y = element_text(angle = 0, size = 6), legend.key.height = unit(5, "mm"))
fig_s18 <- pA / (pB + pC + plot_layout(widths = c(1, 1))) + plot_layout(heights = c(95, 120))
figlog_S18 <- supp_save(fig_s18, "Figure_S18", WIDTH_2COL, 228)
LOG_Figure_S18 <- supp_note("Figure_S18", one_animal = sup18[!sup18$group_tested %in% c(TRUE, "TRUE"), ], n_edges = length(edge_lv), ec = EC_ST, cm = CM_ST,
          n_animals_range = range(ov18$n_animals), npairs_range = range(ov18$n_pairs))
time_supp_fig_s18 <- as.numeric(difftime(Sys.time(), t0_supp, units = "secs"))

t0_supp <- Sys.time()
C5 <- function(f) RES("05c_ec_capillary_sensitivity", f)
REFIT_LAB <- c(no_cytoplasm_rich = "without the cytoplasm-rich stratum", no_strata = "without both strata")
FIT_LAB   <- c(full = "all nuclei", REFIT_LAB)
FIT_COL   <- c(full = "grey45", no_cytoplasm_rich = "#E69F00", no_strata = "#0072B2")
strata  <- read.csv(C5("ec_capillary_strata_per_animal.csv"), stringsAsFactors = FALSE)
aw      <- read.csv(C5("ec_capillary_animal_weight.csv"), stringsAsFactors = FALSE)
hits    <- read.csv(C5("ec_capillary_hit_agreement.csv"), stringsAsFactors = FALSE)
gagr    <- read.csv(C5("ec_capillary_gsea_agreement.csv"), stringsAsFactors = FALSE)
refsum  <- read.csv(C5("ec_capillary_refit_summary.csv"), stringsAsFactors = FALSE)
genes5c <- read.csv(C5("ec_capillary_leading_and_candidate_genes.csv"), stringsAsFactors = FALSE)
paired  <- fread(C5("ec_capillary_paired_statistics.csv"))
grefit  <- fread(C5("ec_capillary_gsea_refits.csv"))
stopifnot(setequal(strata$individual, SAMPLE_ORDER))

sa <- strata |> dplyr::transmute(sample = individual, low_content = 100 * low_content / capillary_nuclei, cytoplasm_rich = 100 * cytoplasm_rich / capillary_nuclei, capillary_nuclei) |>
  tidyr::pivot_longer(c(low_content, cytoplasm_rich), names_to = "stratum", values_to = "pct") |>
  dplyr::mutate(stratum = factor(stratum, levels = c("cytoplasm_rich", "low_content"), labels = c("cytoplasm-rich stratum", "low-content stratum")))
yl19 <- scale_y_discrete(limits = rev(SAMPLE_ORDER), labels = sample_lab)
pA <- ggplot(sa, aes(pct, sample, fill = stratum)) + geom_col(width = .75) + yl19 +
  geom_text(data = strata, aes(x = strata_pct + .8, y = individual, label = sfmt_n(capillary_nuclei)), inherit.aes = FALSE, hjust = 0, size = PT6, colour = "grey35") +
  scale_fill_manual(values = c(`cytoplasm-rich stratum` = FIT_COL[["no_cytoplasm_rich"]], `low-content stratum` = FIT_COL[["no_strata"]]), name = NULL) +
  scale_x_continuous(limits = c(0, max(strata$strata_pct) * 1.35), expand = c(0, 0)) + guides(fill = guide_legend(ncol = 1)) +
  labs(x = "capillary nuclei in the strata (%)\n(n capillary nuclei at right)", y = NULL) + theme_supp() +
  theme(axis.text.y = element_text(colour = rev(sample_col(SAMPLE_ORDER))), legend.position = "top")

hb <- hits |> dplyr::select(refit, contrast, shared, full_only, refit_only, nuclei_kept_pct) |>
  tidyr::pivot_longer(c(shared, full_only, refit_only), names_to = "set", values_to = "n") |>
  dplyr::mutate(set = factor(set, levels = c("refit_only", "full_only", "shared"), labels = c("refit only", "all nuclei only", "in both")),
                contrast = factor(contrast, levels = rev(CONTRAST_ORDER)), refit = factor(refit, levels = names(REFIT_LAB)))
hb_lab <- hits |> dplyr::mutate(contrast = factor(contrast, levels = rev(CONTRAST_ORDER)), refit = factor(refit, levels = names(REFIT_LAB)),
                                lab = sprintf("%.0f %%", nuclei_kept_pct), x = shared + full_only + refit_only)
pB <- ggplot(hb, aes(n, contrast, fill = set)) + geom_col(width = .7) +
  geom_text(data = hb_lab, aes(x = x, y = contrast, label = lab), inherit.aes = FALSE, hjust = -.06, size = PT6, colour = "grey30") +
  facet_wrap(~ refit, nrow = 1, labeller = as_labeller(REFIT_LAB)) +
  scale_fill_manual(values = c(`in both` = "grey45", `all nuclei only` = "grey75", `refit only` = "#56B4E9"), breaks = c("in both", "all nuclei only", "refit only"), name = NULL) +
  scale_y_discrete(labels = CONTRAST_LABEL) + scale_x_continuous(expand = expansion(mult = c(0, .25)), labels = label_number(scale_cut = cut_short_scale())) +
  labs(x = "reported MAST hits (Capillary EC); label, nuclei kept", y = NULL) + theme_supp() + theme(legend.position = "top", strip.text = element_text(size = 6))

pc <- paired[hit_full == TRUE & is.finite(avg_log2FC_refit)][, `:=`(refit = factor(refit, levels = names(REFIT_LAB)), contrast = factor(contrast, levels = CONTRAST_ORDER),
                                                                    kept = ifelse(hit_refit %in% TRUE, "reported in the refit", "not reported in the refit"))][order(-rank(kept))]
rtxt <- hits |> dplyr::mutate(refit = factor(refit, levels = names(REFIT_LAB)), contrast = factor(contrast, levels = CONTRAST_ORDER), lab = sprintf("r = %.2f", r_effect))
lim_c <- quantile(abs(c(pc$avg_log2FC_full, pc$avg_log2FC_refit)), .999)
pC <- ggplot(pc, aes(avg_log2FC_full, avg_log2FC_refit, colour = kept)) +
  geom_abline(slope = 1, intercept = 0, linetype = "22", colour = "grey60", linewidth = .25) +
  scattermore::geom_scattermore(pointsize = 3, pixels = px_for(34, 26)) +
  geom_text(data = rtxt, aes(x = -lim_c, y = lim_c, label = lab), inherit.aes = FALSE, hjust = 0, vjust = 1, size = PT6) +
  facet_grid(refit ~ contrast, labeller = labeller(refit = as_labeller(c(no_cytoplasm_rich = "without\ncytoplasm-rich", no_strata = "without\nboth strata")), contrast = as_labeller(CONTRAST_LABEL))) +
  scale_colour_manual(values = c(`reported in the refit` = "#0072B2", `not reported in the refit` = "#D55E00"), name = NULL) +
  guides(colour = guide_legend(override.aes = list(pointsize = 5))) + coord_cartesian(xlim = c(-lim_c, lim_c), ylim = c(-lim_c, lim_c)) +
  labs(x = "effect, all nuclei (log2FC; Age x HFpEF: interaction effect, MAST, log2 scale)", y = "effect, refit") + theme_supp() +
  theme(legend.position = "right", strip.text = element_text(size = 6), panel.spacing = unit(1.5, "mm"))

TARGET_GENES <- nb_const("05_EC_DE_analysis.Rmd", "TARGET_GENES")
gd <- genes5c |> dplyr::mutate(fit = factor(fit, levels = names(FIT_LAB)), contrast = factor(contrast, levels = CONTRAST_ORDER),
                               eff = pmax(pmin(avg_log2FC, 2), -2), tested = is.finite(p_val), reported = reported %in% TRUE,
                               gene_class = factor(ifelse(gene_class == "leading gene", "leading genes", "candidate genes"), levels = c("leading genes", "candidate genes")))
g_order <- unique(genes5c$gene); gd$gene <- factor(gd$gene, levels = rev(g_order))
pD <- ggplot(gd, aes(fit, gene)) + geom_tile(aes(fill = ifelse(tested, eff, NA)), colour = "white", linewidth = .3) +
  geom_point(data = dplyr::filter(gd, reported), size = .5, colour = "black") +
  geom_text(data = dplyr::filter(gd, !tested), label = "\u2013", colour = "grey45", size = PT7) +
  facet_grid(gene_class ~ contrast, scales = "free_y", space = "free_y", labeller = labeller(contrast = as_labeller(CONTRAST_LABEL))) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-2, 2), na.value = "grey92", name = "effect\n(log2 scale;\ncapped at \u00b12)") +
  scale_x_discrete(labels = c(full = "all", no_cytoplasm_rich = "\u2212cyto", no_strata = "\u2212both")) +
  labs(x = NULL, y = NULL) + theme_supp() +
  theme(axis.text.y = element_text(face = "italic"), axis.line = element_blank(), axis.ticks = element_blank(), strip.text = element_text(size = 6),
        strip.text.y = element_text(angle = 0), panel.spacing.x = unit(1.5, "mm"), legend.key.height = unit(4, "mm"))

gfull <- gsea[celltype == "Capillary EC", .(contrast, geneset, pathway, NES_full = NES, sig_full = sig)]
ge <- merge(grefit[, .(fit, contrast, geneset, pathway, NES_refit = NES, sig_refit = is.finite(padj) & padj < GSEA_FDR)], gfull, by = c("contrast", "geneset", "pathway"))
ge[, `:=`(refit = factor(fit, levels = names(REFIT_LAB)), contrast = factor(contrast, levels = CONTRAST_ORDER),
          cls = factor(fcase(sig_full & sig_refit, "FDR < 0.05 in both", sig_full | sig_refit, "FDR < 0.05 in one", default = "not significant"),
                       levels = c("not significant", "FDR < 0.05 in one", "FDR < 0.05 in both")))]
setorder(ge, cls)
r_mine <- ge[, .(r_mine = cor(NES_full, NES_refit, use = "complete.obs")), by = .(refit, contrast)]
etxt <- gagr |> dplyr::mutate(refit = factor(refit, levels = names(REFIT_LAB)), contrast = factor(contrast, levels = CONTRAST_ORDER), lab = sprintf("r = %.2f", r_NES))
pE <- ggplot(ge, aes(NES_full, NES_refit, colour = cls)) +
  geom_abline(slope = 1, intercept = 0, linetype = "22", colour = "grey60", linewidth = .25) +
  scattermore::geom_scattermore(data = ge[cls == "not significant"], pointsize = 2.5, pixels = px_for(34, 26)) +
  geom_point(data = ge[cls != "not significant"], size = .6) +
  geom_text(data = etxt, aes(x = -3.2, y = 3.2, label = lab), inherit.aes = FALSE, hjust = 0, vjust = 1, size = PT6) +
  facet_grid(refit ~ contrast, labeller = labeller(refit = as_labeller(c(no_cytoplasm_rich = "without\ncytoplasm-rich", no_strata = "without\nboth strata")), contrast = as_labeller(CONTRAST_LABEL))) +
  scale_colour_manual(values = c(`not significant` = "grey80", `FDR < 0.05 in one` = "#E69F00", `FDR < 0.05 in both` = "#0072B2"), name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) + coord_cartesian(xlim = c(-3.3, 3.3), ylim = c(-3.3, 3.3)) +
  labs(x = "NES, all nuclei", y = "NES, refit") + theme_supp() + theme(legend.position = "right", strip.text = element_text(size = 6), panel.spacing = unit(1.5, "mm"))
fig_s19 <- (wrap(pA) + wrap(pB) + plot_layout(widths = c(.8, 1.4))) / pC / wrap(pD) / pE + plot_layout(heights = c(52, 44, 62, 44))
figlog_S19 <- supp_save(fig_s19, "Figure_S19", WIDTH_2COL, 228)
LOG_Figure_S19 <- supp_note("Figure_S19", strata = strata, aw = aw, hits = hits, gagr = gagr, refsum = refsum, genes = genes5c, r_mine = r_mine,
                            n_lead = length(unique(genes5c$gene[genes5c$gene_class == "leading gene"])), target_genes = TARGET_GENES,
                            lead_targets = intersect(TARGET_GENES, genes5c$gene[genes5c$gene_class == "leading gene"]),
                            nes_sets = ge[, .N, by = refit], nes_flip_sig = ge[(sig_full | sig_refit) & sign(NES_full) != sign(NES_refit), .N, by = refit],
                            nes_flip_all = ge[sign(NES_full) != sign(NES_refit), .N, by = refit],
                            flipped = paired[hit_full == TRUE & is.finite(avg_log2FC_refit), .(n = sum(sign(avg_log2FC_full) != sign(avg_log2FC_refit)), tested = .N), by = refit])
time_supp_fig_s19 <- as.numeric(difftime(Sys.time(), t0_supp, units = "secs"))

SUPP_FLAG_SHORT  <- c(low_yield_library = "low yield", manual_doublet_threshold = "manual doublet thr.")
SUPP_ANIMAL_KEY  <- c(setNames(sprintf("%s (animal note)", names(ANIMAL_NOTES)), names(ANIMAL_NOTES)),
                      setNames(sprintf("%s (%s)", names(FLAG_OF), SUPP_FLAG_SHORT[FLAG_OF]), names(FLAG_OF)), other = "other animals")
SUPP_ANIMAL_SHAPE <- setNames(c(rep(17, length(ANIMAL_NOTES)), c(1, 0, 5)[seq_along(FLAG_OF)], 16), SUPP_ANIMAL_KEY)
supp_animal_class <- function(s) factor(unname(ifelse(s %in% names(SUPP_ANIMAL_KEY), SUPP_ANIMAL_KEY[s], SUPP_ANIMAL_KEY[["other"]])), levels = SUPP_ANIMAL_KEY)
short_ec <- function(x) sub(" EC$", "", x)

per_animal_plot <- function(d, ylab, nrow = 1, ncol = NULL) {
  d <- d |> dplyr::mutate(group = factor(group, levels = GRP_ORDER), x = as.numeric(group), animal = supp_animal_class(sample))
  d$x_j <- d$x + withr::with_seed(1, runif(nrow(d), -.16, .16))
  mn <- d |> dplyr::group_by(facet, group, x) |> dplyr::summarise(m = mean(value), .groups = "drop")
  ggplot(d) + geom_segment(data = mn, aes(x = x - .32, xend = x + .32, y = m, yend = m, colour = group), linewidth = .8) +
    geom_point(aes(x_j, value, shape = animal), size = .9, stroke = .4, colour = "grey10") +
    facet_wrap(~ facet, scales = "free_y", nrow = nrow, ncol = ncol) +
    scale_x_continuous(breaks = seq_along(GRP_ORDER), labels = GRP_ORDER, expand = expansion(add = .5)) +
    scale_colour_manual(values = GRP_FILL, labels = GRP_LABEL, name = NULL) +
    scale_shape_manual(values = SUPP_ANIMAL_SHAPE, name = NULL, drop = FALSE) +
    labs(x = NULL, y = ylab) + theme_supp() + theme(strip.text = element_text(size = 6.5, face = "bold"), panel.spacing = unit(2, "mm"))
}

mast_heatmap <- function(d, cap = 1.2, ncol = NULL, nrow = NULL) {
  d <- d |> dplyr::mutate(eff = pmax(pmin(avg_log2FC, cap), -cap), tested = is.finite(avg_log2FC),
                          lab = ifelse(tested, mark, "\u2013"), lab = ifelse(is.na(lab), "", lab))
  ggplot(d, aes(contrast, celltype)) + geom_tile(aes(fill = ifelse(tested, eff, NA)), colour = "white", linewidth = .4) +
    geom_text(aes(label = lab, colour = tested), size = PT6, show.legend = FALSE) +
    facet_wrap(~ gene, ncol = ncol, nrow = nrow) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-cap, cap), na.value = "grey90",
                         name = sprintf("effect\n(log2 scale,\ncapped at \u00b1%s)", cap)) +
    scale_colour_manual(values = c(`TRUE` = "grey10", `FALSE` = "grey50")) +
    scale_x_discrete(labels = CONTRAST_CODE) + scale_y_discrete(labels = short_ec) + labs(x = NULL, y = NULL) + theme_supp() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1),
          strip.text = element_text(size = 6.5, face = "bold.italic"), panel.spacing = unit(1.2, "mm"), legend.key.height = unit(4, "mm"))
}

nb_py_src <- function(ipynb) { nb <- jsonlite::fromJSON(file.path(FIG_ROOT, "notebooks", ipynb), simplifyVector = FALSE)
  paste(unlist(lapply(Filter(function(c) c$cell_type == "code", nb$cells), function(c) paste(unlist(c$source), collapse = ""))), collapse = "\n") }
nb_py_value <- function(src, name) { m <- regmatches(src, regexpr(sprintf("(?m)^%s\\s*=\\s*.+$", name), src, perl = TRUE)); stopifnot(length(m) == 1)
  trimws(sub("^[^=]+=\\s*", "", m)) }
py_strings <- function(x) gsub("^'|'$", "", regmatches(x, gregexpr("'[^']*'", x))[[1]])

t0_supp <- Sys.time()
B01 <- function(f) RES("01b_integration_benchmark", "tables", f)
src01b  <- nb_py_src("01b_integration_benchmark.ipynb")
METHODS <- py_strings(nb_py_value(src01b, "METHODS"))
METHOD_COL <- setNames(py_strings(sub(".*METHODS,\\s*", "", nb_py_value(src01b, "METHOD_PALETTE"))), METHODS)
read_scib <- function(f) { x <- fread(B01(f)); list(type = unlist(x[method == "Metric Type"][, -1]), res = x[method != "Metric Type"][, lapply(.SD, as.numeric), by = method]) }
s_sample <- read_scib("scib_results_batch_sample_unscaled.csv"); s_seq <- read_scib("scib_results_batch_seq_unscaled.csv")
AGG <- c("Bio conservation", "Batch correction", "Total")
metric_type <- s_sample$type[setdiff(names(s_sample$type), AGG)]
stopifnot(setequal(s_sample$res$method, METHODS), all(metric_type %in% AGG[1:2]))
ranked <- fread(B01("scib_ranked_batch_sample.csv"), skip = 3, header = FALSE)[, .(method = V1, rank_unscaled = V8)]
lf  <- fread(B01("label_free_check.csv")); setnames(lf, 1, "method")
pur <- fread(B01("knn_panel_purity_by_panel_confident.csv")); setnames(pur, 1, "method")
pan <- fread(B01("panel_label_summary.csv"))
silb <- fread(B01("silhouette_batch_sample.csv")); setnames(silb, c("method", "silhouette_batch"))
umap_meta <- fread(file.path(PANELS_DIR, "supp_umap_by_method_meta.csv"))

w_total <- round(coef(lm(Total ~ 0 + `Bio conservation` + `Batch correction`, data = s_sample$res)), 2)
label_free_batch <- sub(" \\(sample\\)$", "", grep(" \\(sample\\)$", names(lf), value = TRUE))
label_based <- setdiff(names(metric_type), label_free_batch)
meth_order <- s_sample$res[order(-Total), method]

GRP_A <- c(`Bio conservation` = "bio conservation", `Batch correction` = "batch correction", aggregate = "aggregate scores", rank = " ", extra = "  ")
ta <- melt(s_sample$res, id.vars = "method", variable.name = "metric", value.name = "value", variable.factor = FALSE)
ta[, grp := ifelse(metric %in% AGG, "aggregate", metric_type[metric])]
ta <- rbind(ta, ranked[, .(method, metric = "rank (total)", value = NA_real_, grp = "rank", txt = as.character(rank_unscaled))],
            silb[, .(method, metric = "silhouette batch", value = silhouette_batch, grp = "extra")], fill = TRUE)
ta[is.na(txt), txt := fmt_f(value, 2)]
ta[, `:=`(method = factor(method, levels = rev(meth_order)), grp = factor(GRP_A[grp], levels = GRP_A),
          metric_lab = ifelse(metric %in% label_based, paste0(metric, " \u00a7"), metric))]
ta[, metric_lab := factor(metric_lab, levels = unique(metric_lab[order(grp, match(metric, c(names(metric_type), AGG, "rank (total)", "silhouette batch")))]))]
pA <- ggplot(ta, aes(metric_lab, method)) + geom_tile(aes(fill = value), colour = "white", linewidth = .5) +
  geom_text(aes(label = txt, colour = !is.na(value) & value > .72), size = PT6, show.legend = FALSE) +
  facet_grid(~ grp, scales = "free_x", space = "free_x", labeller = label_wrap_gen(14)) +
  scale_fill_gradient(low = "#F2F6FB", high = "#2E6FB0", limits = c(0, 1), na.value = "white", name = "score\n(0-1)") +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey10")) + scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL) + theme_supp() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x.top = element_text(angle = 40, hjust = 0, vjust = 0),
        axis.text.y = element_text(colour = METHOD_COL[rev(meth_order)], face = "bold"), strip.text = element_text(size = 6.5, face = "bold"),
        strip.placement = "outside", panel.spacing.x = unit(1.2, "mm"), legend.key.height = unit(4, "mm"))

n_runs <- length(unique(animals$seq))
KEY_LAB <- c(sample = sprintf("batch = sample (%d libraries)", length(SAMPLE_ORDER)), seq = sprintf("batch = sequencing run (%d runs)", n_runs))
tb <- rbind(s_sample$res[, .(method, bio = `Bio conservation`, batch = `Batch correction`, total = Total, key = KEY_LAB[["sample"]])],
            s_seq$res[, .(method, bio = `Bio conservation`, batch = `Batch correction`, total = Total, key = KEY_LAB[["seq"]])])
tb[, key := factor(key, levels = KEY_LAB)]
harm <- grep("^Harmony", METHODS, value = TRUE)
tb[, lab := fifelse(method %in% harm, fifelse(method == harm[1], paste(harm, collapse = " /\n"), ""), method)]
iso <- data.table(total = seq(.55, .8, by = .05))[, .(intercept = total / w_total[[1]], slope = -w_total[[2]] / w_total[[1]], total)]
pB <- ggplot(tb, aes(batch, bio, colour = method)) +
  geom_abline(data = iso, aes(intercept = intercept, slope = slope), colour = "grey88", linewidth = .25) +
  geom_point(size = 1.3) + ggrepel::geom_text_repel(aes(label = lab), size = PT6, lineheight = .85, max.overlaps = Inf, segment.size = .2, min.segment.length = 0,
                                                     box.padding = .15, point.padding = .1, seed = 1, show.legend = FALSE) +
  facet_wrap(~ key, ncol = 1) + scale_colour_manual(values = METHOD_COL, guide = "none") +
  labs(x = "batch correction (unscaled)", y = "bio conservation (unscaled)") + theme_supp() +
  theme(strip.text = element_text(size = 6.5, face = "bold"))

lf_long <- melt(lf[, .(method, all = `kNN panel purity, all nuclei`, confident = `kNN panel purity, confident nuclei`)], id.vars = "method", variable.name = "set")
lf_long[, `:=`(method = factor(method, levels = rev(meth_order)), set = factor(set, levels = c("all", "confident"), labels = c("all nuclei", "confident nuclei")))]
pC1 <- ggplot(lf_long, aes(value, method)) + geom_line(aes(group = method), colour = "grey75", linewidth = .4) +
  geom_point(aes(shape = set), size = 1.3, colour = "grey15") + scale_shape_manual(values = c(16, 1), name = NULL) +
  scale_x_continuous(breaks = scales::breaks_pretty(3)) +
  labs(x = "purity (proportion)", y = NULL) + theme_supp() + guides(shape = guide_legend(ncol = 1)) +
  theme(axis.text.y = element_text(colour = METHOD_COL[rev(meth_order)], face = "bold"), legend.position = "top")
panel_order <- pan[order(-`confident nuclei`), panel_label]
tc <- melt(pur, id.vars = "method", variable.name = "panel", value.name = "purity", variable.factor = FALSE)
tc[, `:=`(method = factor(method, levels = rev(meth_order)), panel = factor(sprintf("%s (%s)", panel, sfmt_n(pan$`confident nuclei`[match(panel, pan$panel_label)])),
                                                                            levels = sprintf("%s (%s)", panel_order, sfmt_n(pan[match(panel_order, panel_label), `confident nuclei`]))))]
pC2 <- ggplot(tc, aes(panel, method, fill = purity)) + geom_tile(colour = "white", linewidth = .4) +
  geom_text(aes(label = sub("^0", "", fmt_f(purity, 2)), colour = purity > .65), size = PT6, show.legend = FALSE) +
  scale_fill_gradient(low = "#F3FAF4", high = "#1B7837", limits = c(0, 1), name = "purity\n(proportion)") +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey10")) + labs(x = NULL, y = NULL) + theme_supp() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_text(angle = 40, hjust = 1),
        legend.key.height = unit(4, "mm"))

um <- fread(file.path(PANELS_DIR, "supp_umap_by_method.csv.gz"))
ct_levels <- c(intersect(CELL_TYPE_ORDER, unique(um$cell_type)), setdiff(unique(um$cell_type), CELL_TYPE_ORDER))
ct_col <- c(CELL_TYPE_FILL, setNames(rep("grey60", length(ct_levels)), ct_levels))[ct_levels]
um <- um[withr::with_seed(0, sample(.N))][, `:=`(method = factor(method, levels = METHODS), cell_type = factor(cell_type, levels = ct_levels))]
um_lab <- umap_meta[, .(method = factor(method, levels = METHODS), lab = ifelse(umap == "recomputed", "", "stored atlas UMAP"))]
pD <- ggplot(um, aes(x, y, colour = cell_type)) + scattermore::geom_scattermore(pointsize = 1.2, pixels = px_for(40, 38)) +
  geom_text(data = um_lab, aes(x = -Inf, y = -Inf, label = lab), inherit.aes = FALSE, hjust = -.05, vjust = -.4, size = PT6, colour = "grey40") +
  facet_wrap(~ method, nrow = 2) + scale_colour_manual(values = ct_col, labels = function(x) gsub("_", " ", x), name = "cell type\n(notebook 01)") +
  guides(colour = guide_legend(override.aes = list(pointsize = 4), ncol = 2)) + theme_supp() + no_axes + labs(x = NULL, y = NULL) +
  theme(strip.text = element_text(size = 6.5, face = "bold"), legend.position = "inside", legend.position.inside = c(.88, .25), aspect.ratio = 1)
fig_s20 <- wrap(pA) / (wrap(pB) + wrap(pC1 + pC2 + plot_layout(widths = c(.4, 1))) + plot_layout(widths = c(.9, 2.1))) / wrap(pD) + plot_layout(heights = c(60, 66, 92))
figlog_S20 <- supp_save(fig_s20, "Figure_S20", WIDTH_2COL, 228)
LOG_Figure_S20 <- supp_note("Figure_S20", sample = s_sample$res, seq = s_seq$res, ranked = ranked, lf = lf, silb = silb, pan = pan, umap_meta = umap_meta,
                            w_total = w_total, label_based = label_based, label_free_batch = label_free_batch, n_nuclei = nrow(um) / length(METHODS),
                            metric_type = metric_type, n_runs = n_runs, n_types = length(ct_levels))
rm(um); invisible(gc())
time_supp_fig_s20 <- as.numeric(difftime(Sys.time(), t0_supp, units = "secs"))

t0_supp <- Sys.time()
K08 <- function(f) RES("08_ec_klf", f)
NB08 <- "08_EC_KLF2_KLF4_shear.Rmd"
SHEAR_PROG <- nb_const(NB08, "SHEAR_PROG"); KLF_GENES <- nb_const(NB08, "KLF_GENES")
LFC_CUT08 <- nb_const(NB08, "LFC_CUT"); PADJ_CUT08 <- nb_const(NB08, "PADJ_CUT"); FDR_CUT08 <- nb_const(NB08, "FDR_CUT"); MIN_PB08 <- nb_const(NB08, "MIN_PB_NUCLEI")
prof <- fread(K08("ec_klf_profile_readouts.csv")); modm <- fread(K08("ec_klf_module_animal_model.csv"))
mg   <- fread(K08("ec_klf_mast_module_genes.csv")); act <- fread(K08("ec_klf_activity_pseudobulk.csv"))
scr  <- fread(K08("ec_klf_tf_activity_screen.csv"))
SUBTYPES08 <- unique(mg$celltype); POOLED <- setdiff(unique(prof$level), SUBTYPES08)
CAP <- grep("^Capillary", SUBTYPES08, value = TRUE); LEVELS_SHOWN <- c(CAP, POOLED)
AGE_CONTRASTS <- grep("^O._vs_Y.$", CONTRAST_ORDER, value = TRUE)
stopifnot(length(CAP) == 1, length(POOLED) == 1, all(prof$sample %in% SAMPLE_ORDER), all(KLF_GENES %in% scr$TF))
module_genes <- unique(mg$gene[mg$status != "not tested (< 1% detection)"])

pa <- prof[level %in% LEVELS_SHOWN, .(sample, group, value = module_score, facet = factor(level, levels = LEVELS_SHOWN), in_model)]
pA <- per_animal_plot(pa, "shear-module score\n(per-animal mean)") + theme(legend.position = "none")

tb <- modm[, .(level = factor(level, levels = rev(c(SUBTYPES08, POOLED))), contrast = factor(contrast, levels = CONTRAST_ORDER), effect, FDR,
               lab = paste0(sprintf("%+.2f", effect), ifelse(FDR < FDR_CUT08, as.character(fdr_stars(FDR)), "")))]
pB <- ggplot(tb, aes(contrast, level, fill = effect)) + geom_tile(colour = "white", linewidth = .5) + geom_text(aes(label = lab), size = PT6) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-1, 1) * max(abs(tb$effect)), name = "module-score\ndifference") +
  scale_x_discrete(labels = CONTRAST_CODE) + scale_y_discrete(labels = short_ec) + labs(x = NULL, y = NULL) + theme_supp() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), legend.key.height = unit(4, "mm"))

tc <- mg[, .(celltype = factor(celltype, levels = rev(SUBTYPES08)), contrast = factor(contrast, levels = CONTRAST_ORDER),
             gene = factor(gene, levels = SHEAR_PROG), avg_log2FC, mark, status)]
pC <- mast_heatmap(tc, nrow = 2)

pd <- melt(prof[level %in% LEVELS_SHOWN, c("sample", "group", "level", paste0(KLF_GENES, "_activity_pseudobulk")), with = FALSE],
           id.vars = c("sample", "group", "level"), variable.name = "tf", value.name = "value")
pd[, facet := factor(sprintf("%s, %s", toupper(sub("_activity_pseudobulk$", "", tf)), level), levels = as.vector(outer(toupper(KLF_GENES), LEVELS_SHOWN, paste, sep = ", ")))]
pD <- per_animal_plot(pd, "regulon activity\n(ULM score, pseudobulk)", nrow = 1) + theme(legend.position = "bottom", legend.box = "vertical", legend.spacing.y = unit(0, "mm"))

te <- scr[level %in% LEVELS_SHOWN & contrast %in% AGE_CONTRASTS][, `:=`(
  facet = factor(sprintf("%s, %s", CONTRAST_LABEL[contrast], level), levels = as.vector(t(outer(CONTRAST_LABEL[AGE_CONTRASTS], LEVELS_SHOWN, paste, sep = ", ")))),
  cls = factor(fifelse(TF %in% KLF_GENES, "KLF2 / KLF4", fifelse(FDR < FDR_CUT08, "other TF, FDR < 0.05", "other TF")), levels = c("other TF", "other TF, FDR < 0.05", "KLF2 / KLF4")))]
setorder(te, cls)
te_lab <- te[TF %in% KLF_GENES, .(facet, effect, p_value, lab = sprintf("%s %d/%d", toupper(TF), rank_p, n_TF))]
pE <- ggplot(te, aes(effect, -log10(p_value), colour = cls)) + geom_point(size = .45) +
  ggrepel::geom_text_repel(data = te_lab, mapping = aes(effect, -log10(p_value), label = lab), inherit.aes = FALSE, size = PT6,
                           colour = "#D55E00", segment.size = .2, min.segment.length = 0, box.padding = .3, seed = 1) +
  facet_wrap(~ facet, nrow = 1) + scale_colour_manual(values = c(`other TF` = "grey80", `other TF, FDR < 0.05` = "grey35", `KLF2 / KLF4` = "#D55E00"), name = NULL) +
  labs(x = "activity difference (ULM score)", y = "\u2212log10 p") + theme_supp() +
  theme(legend.position = "bottom", strip.text = element_text(size = 6.5, face = "bold"))
fig_s21 <- (wrap(pA) + wrap(pB) + plot_layout(widths = c(1, 1.15))) / wrap(pC) / wrap(pD) / wrap(pE) + plot_layout(heights = c(44, 62, 50, 50))
figlog_S21 <- supp_save(fig_s21, "Figure_S21", WIDTH_2COL, 228)
LOG_Figure_S21 <- supp_note("Figure_S21", shear = SHEAR_PROG, klf = KLF_GENES, module_genes = module_genes, lfc = LFC_CUT08, padj = PADJ_CUT08, modm = modm, mg = mg, act = act,
                            prof = prof, cap = CAP, pooled = POOLED, age = AGE_CONTRASTS, subtypes = SUBTYPES08,
                            fdr = FDR_CUT08, min_pb = MIN_PB08, screen_sig = scr[FDR < FDR_CUT08, .N, by = .(level, contrast)], n_tf = scr[, .(n_TF = n_TF[1]), by = level])
time_supp_fig_s21 <- as.numeric(difftime(Sys.time(), t0_supp, units = "secs"))

t0_supp <- Sys.time()
V09 <- function(f) RES("09_hfpef_validation", f)
NB09 <- "09_HFpEF_model_validation.Rmd"
FDR_CUT09 <- nb_const(NB09, "FDR_CUT"); LFC_CUT09 <- nb_const(NB09, "LFC_CUT"); PADJ_CUT09 <- nb_const(NB09, "PADJ_CUT")
conc <- fread(V09("hfpef_concordance_table.csv")); am <- fread(V09("hfpef_animal_model.csv"))
pav  <- fread(V09("hfpef_per_animal_values.csv")); ecm <- fread(V09("hfpef_ec_activation_mast.csv"))
ecs  <- fread(V09("hfpef_ec_activation_mast_summary.csv")); gss <- fread(V09("hfpef_gene_set_sizes.csv")); labf <- fread(V09("hfpef_lab_findings_reported.csv"))
stopifnot(all(pav$sample %in% SAMPLE_ORDER))

cr <- conc[!is.na(readout)][, row := .I]
cells <- rbind(am[, .(readout = label, contrast, effect, sig = FDR < FDR_CUT09, n_lower = NA_integer_, n_higher = NA_integer_)],
               ecs[, .(readout = sprintf("%s (EC, MAST)", gene), contrast, effect = median_effect, sig = NA, n_lower, n_higher)])
stopifnot(all(cr$readout %in% cells$readout))
ta <- merge(cr[, .(row, phenotype, readout, basis, lab_direction, lab_contrasts)], cells, by = "readout", allow.cartesian = TRUE)
ta[, `:=`(lab_contrast = mapply(function(cs, cn) !is.na(cs) && cn %in% strsplit(cs, ";")[[1]], lab_contrasts, contrast),
          lab_sign = c(up = 1, down = -1)[lab_direction], is_mast = !is.na(n_lower))]
ta[, side := fifelse(is_mast & lab_contrast, lab_sign, fifelse(is_mast & n_higher > 0 & n_lower == 0, 1, fifelse(is_mast & n_lower > 0 & n_higher == 0, -1, sign(effect))))]
ta[is_mast == TRUE, sig := fifelse(side > 0, n_higher > 0, n_lower > 0)]
ta[, opposite := lab_contrast & ((is_mast & fifelse(lab_sign > 0, n_lower > 0, n_higher > 0) %in% TRUE) | (!is_mast & sign(effect) == -lab_sign & sig %in% TRUE))]
ta[, `:=`(fill_dir = factor(fifelse(sign(effect) > 0, "higher", "lower"), levels = c("higher", "lower")),
          arrow = fifelse(lab_contrast, fifelse(lab_sign > 0, "\u2191", "\u2193"), ""),
          mark = fifelse(opposite, "\u00d7", fifelse(sig %in% TRUE, "\u2022", "")))]
ta[, `:=`(txt = paste0(arrow, mark), contrast = factor(contrast, levels = CONTRAST_ORDER),
          readout = factor(readout, levels = rev(cr$readout)), phenotype = factor(phenotype, levels = unique(cr$phenotype)))]

frac <- function(x) { v <- do.call(rbind, strsplit(x[x != "n/a"], "/")); c(k = sum(as.integer(v[, 1])), n = sum(as.integer(v[, 2]))) }
conc_tot <- rbind(direction = frac(cr$concordant_direction), significant = frac(cr$significant), opposite = frac(cr$opposite_significant))
mine <- ta[lab_contrast == TRUE, .(direction = sum(sign(effect) == lab_sign), significant = sum(sig %in% TRUE & (is_mast | sign(effect) == lab_sign)), opposite = sum(opposite), n = .N)]
stopifnot(mine$n == conc_tot["direction", "n"], mine$direction == conc_tot["direction", "k"], mine$significant == conc_tot["significant", "k"], mine$opposite == conc_tot["opposite", "k"])
SIDE_COLS <- c(dir = "dir.", thr = "thr.", opp = "opp.")
side_txt <- melt(cr[, .(readout = factor(readout, levels = rev(cr$readout)), phenotype = factor(phenotype, levels = unique(cr$phenotype)),
                        dir = concordant_direction, thr = significant, opp = opposite_significant)], id.vars = c("readout", "phenotype"), variable.name = "col", value.name = "lab")
side_txt[, col := factor(SIDE_COLS[as.character(col)], levels = c(CONTRAST_ORDER, SIDE_COLS))]
pA <- ggplot(ta, aes(contrast, readout)) + geom_tile(aes(fill = fill_dir), colour = "white", linewidth = .4) +
  geom_text(aes(label = txt), size = PT7, lineheight = .8) +
  geom_text(data = side_txt, aes(x = col, y = readout, label = lab), size = PT6, colour = "grey25") +
  facet_grid(phenotype ~ ., scales = "free_y", space = "free_y", labeller = label_wrap_gen(14)) +
  scale_fill_manual(values = c(higher = "#F4C7B8", lower = "#C3DBEC"), name = "snRNA effect", drop = FALSE) +
  scale_x_discrete(limits = c(CONTRAST_ORDER, SIDE_COLS), labels = c(CONTRAST_CODE, setNames(SIDE_COLS, SIDE_COLS))) +
  labs(x = NULL, y = NULL) + theme_supp() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), strip.text.y = element_text(angle = 0, hjust = 0, size = 6, face = "bold"),
        panel.spacing.y = unit(.8, "mm"), legend.position = "bottom")

MODULES_B <- c("hypertrophy", "calcium", "metabolic", "fb_activation", "senescence")
pb <- pav[(kind == "module" & feature %in% MODULES_B) | (kind == "gene" & feature == "Atp2a2" & lineage == "Cardiomyocyte")]
pb_order <- unique(pb[order(match(feature, c(MODULES_B[1:2], "Atp2a2", MODULES_B[-(1:2)])), lineage), readout])
pb[, facet := factor(readout, levels = pb_order)]
pB <- per_animal_plot(pb[, .(sample, group, value, facet)], "per-animal mean (module score;\nAtp2a2: log-normalised expression)", nrow = 2) +
  facet_wrap(~ facet, scales = "free_y", nrow = 2, labeller = label_wrap_gen(16)) +
  theme(legend.position = "right", legend.box = "vertical", strip.clip = "off")

tc <- ecm[, .(celltype = factor(celltype, levels = rev(unique(ecm$celltype))), contrast = factor(contrast, levels = CONTRAST_ORDER),
              gene = factor(gene, levels = unique(ecm$gene)), avg_log2FC, mark, status)]
pC <- mast_heatmap(tc, ncol = 1)
fig_s22 <- (wrap(pA) + wrap(pC) + plot_layout(widths = c(1.75, 1))) / wrap(pB) + plot_layout(heights = c(140, 80))
figlog_S22 <- supp_save(fig_s22, "Figure_S22", WIDTH_2COL, 228)
LOG_Figure_S22 <- supp_note("Figure_S22", conc = conc, cr = cr, conc_tot = conc_tot, mine = mine, ta = ta, am = am, ecm = ecm, ecs = ecs, gss = gss, labf = labf,
                            fdr = FDR_CUT09, lfc = LFC_CUT09, padj = PADJ_CUT09, pb_order = pb_order,
                            no_counterpart = conc[is.na(readout), lab_finding])
time_supp_fig_s22 <- as.numeric(difftime(Sys.time(), t0_supp, units = "secs"))

t0_supp <- Sys.time()

ret_m <- read.csv(T01("qc_validation_retained_metrics_by_sample.csv"), stringsAsFactors = FALSE)
clean_an <- read.csv(RES("02_subclustering", "atlas_cleaning_by_animal.csv"), stringsAsFactors = FALSE)
scr <- read.csv(T01("scrublet_summary.csv"), stringsAsFactors = FALSE) |> dplyr::filter(source == "cellranger")
sep_r <- read.csv(T01("qc_validation_scrublet_separation.csv"), stringsAsFactors = FALSE) |> dplyr::filter(source == "cellranger")
t1 <- animals |> dplyr::transmute(library = sample, group = as.character(group), age, condition, sex, sequencing_run = seq, batch,
                           qc_flag = flag_short(sample), animal_note = unname(ifelse(sample %in% NOTED, ANIMAL_NOTES[sample], ""))) |>
  dplyr::left_join(sqf |> dplyr::transmute(library = sample, cellranger_estimated_cells = estimated_cells, cellranger_reads = reads, cellranger_reads_in_cells_pct = reads_in_cells_pct,
                             cellranger_intronic_reads_pct = intronic_reads_pct, cellranger_median_genes_per_cell = median_genes_per_cell,
                             recorded_sex, inferred_sex, empty_droplet_flag, umis_at_droplet_rank_25000 = umi_last_cellbender_droplet), by = "library") |>
  dplyr::left_join(rv |> dplyr::transmute(library = sample, recount_cellranger_nuclei = cellranger_nuclei, recount_nuclei_matched = nuclei_in_recount,
                            median_ratio_exon_only_to_cellranger, pearson_r_nucleus_log1p_umis, pearson_r_gene_log1p_totals,
                            median_ratio_intron_inclusive_to_cellranger, median_umis_cellranger, median_umis_intron_inclusive), by = "library") |>
  dplyr::left_join(lib |> dplyr::transmute(library = sample, starsolo_cell_calls = n_cellranger, cellbender_cells = cb_found_cells, cellbender_droplets_nonempty_prop = cb_droplets_nonempty_prop,
                             cellbender_counts_removed_pct = 100 * cb_counts_removed_prop, qc_input_nuclei = n_cellbender, qc_removed_nuclei = n_removed,
                             qc_removed_pct = 100 * removed_prop, qc_retained_nuclei = n_retained), by = "library") |>
  dplyr::left_join(exon_by_lib |> as.data.frame() |> dplyr::transmute(library = sample, median_exon_prop_qc_input = median_exon_prop), by = "library") |>
  dplyr::left_join(sep_r |> dplyr::transmute(library = sample, scrublet_auroc_sim_gt_obs = auroc_sim_gt_obs, scrublet_threshold_auto = threshold_auto,
                               scrublet_threshold_used = threshold_used, scrublet_threshold_method = method, scrublet_called_doublets_pct = 100 * predicted_prop), by = "library") |>
  dplyr::left_join(ret_m |> dplyr::transmute(library = sample, retained_median_genes = median_genes, retained_median_umis = median_umi, retained_median_mito_pct = median_mito_pct,
                               atlas_nuclei = n_final_atlas), by = "library") |>
  dplyr::left_join(clean_an |> dplyr::transmute(library = sample, cleaned_atlas_nuclei = retained, atlas_cleaning_removed_pct = removed_pct), by = "library") |>
  dplyr::left_join(lisi_s |> dplyr::transmute(library = sample, ilisi_median_pca = pre_harmony, ilisi_median_harmony = post_harmony), by = "library") |>
  dplyr::left_join(sexm |> dplyr::transmute(library = sample, y_genes_mean_expression = Y, xist_mean_expression = Xist), by = "library")
stopifnot(nrow(t1) == 16, !anyNA(t1$atlas_nuclei), sum(t1$cleaned_atlas_nuclei) == sum(clean_by_ct$retained))
write_table(t1, "Table_S1_library_qc_recount")

bg2 <- read.csv(RES("02_subclustering", "atlas_cleaning_by_group.csv"), stringsAsFactors = FALSE)
bc2 <- read.csv(RES("02_subclustering", "atlas_cleaning_by_cell_type.csv"), stringsAsFactors = FALSE)
relab2 <- read.csv(RES("02_subclustering", "atlas_cleaning_clusters.csv"), stringsAsFactors = FALSE) |> dplyr::filter(decision == "relabelled")
t2 <- dplyr::bind_rows(bg2 |> dplyr::filter(cell_type != "all") |> dplyr::mutate(group_label = unname(GROUP_LABEL[group])),
                bc2 |> dplyr::transmute(group = "all", cell_type, nuclei, removed, retained, removed_pct = pct_removed, group_label = "all groups")) |>
  dplyr::transmute(lineage = cell_type, group, group_label, nuclei, removed_nuclei = removed, retained_nuclei = retained, removed_pct = round(100 * removed / nuclei, 1),
            note = dplyr::case_when(lineage %in% relab2$cell_type ~ sprintf("retained nuclei include %s nuclei of sub-cluster %s relabelled %s", sfmt_n(relab2$n[match(lineage, relab2$cell_type)]),
                                                                     relab2$cluster[match(lineage, relab2$cell_type)], relab2$identity[match(lineage, relab2$cell_type)]),
                             lineage == "Epicardial" ~ "not sub-clustered (157-nucleus notebook-01 cluster relabelled Epicardial)", TRUE ~ "")) |>
  dplyr::mutate(lineage = factor(lineage, levels = CELL_TYPE_ORDER), group = factor(group, levels = c(GROUP_ORDER, "all"))) |> dplyr::arrange(lineage, group)
stopifnot(all(abs(t2$removed_pct[t2$group == "all"] - bc2$pct_removed[match(t2$lineage[t2$group == "all"], bc2$cell_type)]) < .06))
write_table(t2, "Table_S2_atlas_cleaning")

mast <- fread(RES("05_ec_de", "ec_DE_all_MAST_annotated.csv"))
secretome <- fread(DATA("resources", "mouse_secretome_swissprot.csv"))
sec_col <- intersect(c("gene", "Gene", "gene_symbol", "symbol"), names(secretome))[1]
t3 <- mast[reported == TRUE][, .(ec_subtype = celltype, contrast, contrast_label = unname(CONTRAST_LABEL[contrast]), gene,
                                 effect_type = ifelse(contrast == "Interaction", "interaction effect (MAST, log2 scale)", "log2FC (Seurat avg_log2FC)"),
                                 effect_log2 = avg_log2FC, p_value = p_val, p_adj_bonferroni = p_val_adj, detected_group1_prop = pct.1, detected_group2_prop = pct.2,
                                 n_animals_detected_group1 = n_animals_detected_i1, n_animals_detected_group2 = n_animals_detected_i2,
                                 direction_concordance_prop = direction_concordance, max_animal_share_prop = max_animal_share, max_share_animal,
                                 single_animal_driven, low_detection, background_suspect, secreted = gene %in% secretome[[sec_col]])]
t3 <- t3[order(match(ec_subtype, PB_SUBTYPES), match(contrast, CONTRAST_ORDER), p_value)]
write_table(t3, "Table_S3_MAST_reported_hits")

t4a <- fread(PB("ec_pseudobulk_DE.csv"))[, .(ec_subtype = celltype, contrast, contrast_label = unname(CONTRAST_LABEL[contrast]), gene, log2FC = logFC, baseMean,
                                              wald_stat = stat, p_value = P.Value, FDR, n_animal_profiles = n_samples, n_genes_tested, gene_flag, FDR_below_0.05 = hit, FDR_below_0.10 = hit_loose)]
write_table(t4a, "Table_S4a_pseudobulk_DESeq2")
t4b <- read.csv(PB("ec_pseudobulk_concordance.csv"), stringsAsFactors = FALSE) |> dplyr::mutate(contrast_label = unname(CONTRAST_LABEL[contrast]), .after = contrast)
write_table(t4b, "Table_S4b_pseudobulk_MAST_concordance")

t5a <- read.csv(RES("06_ec_gsea", "ec_gsea_pathway_panel.csv"), stringsAsFactors = FALSE) |> dplyr::mutate(contrast_label = unname(CONTRAST_LABEL[contrast]), .after = contrast)
write_table(t5a, "Table_S5a_GSEA_pathway_panel")
t5b <- gsea[sig == TRUE, .(database = geneset, ec_subtype = celltype, contrast, contrast_label = unname(CONTRAST_LABEL[contrast]), pathway, NES, p_value = pval, FDR = padj,
                           FDR_all_blocks = padj_all, set_size = size, ranking_universe_genes = universe_n, technical_share_prop = tech_prop, technical_candidate,
                           single_animal_share_prop = sad_share_prop, single_animal = sad_animal, single_animal_candidate, leading_edge = leadingEdge)]
write_table(t5b[order(database, ec_subtype, contrast, FDR)], "Table_S5b_GSEA_significant")

lr07  <- fread(RES("07_ccc", "tables", "lr_diff_all.csv"))[, layer := "broad cell types (07)"]
lr07b <- fread(RES("07b_ec_cm_subtype_ccc", "tables", "lr_subtype_diff_all.csv"))[, layer := "EC and CM subtypes (07b)"]
t6 <- rbindlist(list(lr07, lr07b), fill = TRUE)
setnames(t6, c("t_mod", "p_value"), c("moderated_t", "p_value"), skip_absent = TRUE)
t6[, contrast_label := unname(CONTRAST_LABEL[contrast])]
setcolorder(t6, c("layer", "source", "target", "ligand_complex", "receptor_complex", "contrast", "contrast_label"))
write_table(t6, "Table_S6_LR_animal_tests")
TABLE_SUMMARY <- data.frame(table = c("Table_S1_library_qc_recount", "Table_S2_atlas_cleaning", "Table_S3_MAST_reported_hits", "Table_S4a_pseudobulk_DESeq2",
                                      "Table_S4b_pseudobulk_MAST_concordance", "Table_S5a_GSEA_pathway_panel", "Table_S5b_GSEA_significant", "Table_S6_LR_animal_tests"),
                            rows = c(nrow(t1), nrow(t2), nrow(t3), nrow(t4a), nrow(t4b), nrow(t5a), nrow(t5b), nrow(t6)))
LOG_Tables <- supp_note("Tables", summary = TABLE_SUMMARY, n_mast_reported = nrow(t3), n_single = sum(t3$single_animal_driven), n_low = sum(t3$low_detection),
          n_bkg = sum(t3$background_suspect), n_lr07 = nrow(lr07), n_lr07b = nrow(lr07b), min_fdr07 = min(lr07$FDR, na.rm = TRUE), min_fdr07b = min(lr07b$FDR, na.rm = TRUE),
          n07_fdr05 = sum(lr07$FDR < .05, na.rm = TRUE), n07b_fdr05 = sum(lr07b$FDR < .05, na.rm = TRUE))
knitr::kable(TABLE_SUMMARY, caption = "Supplementary tables written to tables/ (CSV; XLSX when openxlsx is installed).")
time_supp_tables <- as.numeric(difftime(Sys.time(), t0_supp, units = "secs"))

# [Removed from the public code export: generation of the supplementary figure/table legend text. Available from the authors.]
