# figures_main.R
# Code-only export of the analysis pipeline. Data are not included; see README.md.


IN <- list(
  shared = c(RES("03_ec_annotation", "tables", "ec_sample_metadata.csv"), RES("02_subclustering", "atlas_cleaning_removed_barcodes.csv")),
  fig1   = c(OBJ("integrated_LV_snRNAseq_seurat_clean.rds"), RES("01_qc_atlas", "tables", "qc_retention_by_sample.csv"),
             RES("02_subclustering", "atlas_cleaning_by_cell_type.csv")),
  table1 = c(RES("01_qc_atlas", "tables", "sample_qc_flags.csv"), RES("04_cm_annotation", "tables", "cm_sample_metadata.csv"),
             RES("07_ccc", "tables", "animal_qc.csv")),
  fig2   = c(RES("02_subclustering", "Endothelial_subcluster.rds"), RES("02_subclustering", "Cardiomyocyte_subcluster.rds"),
             RES("03_ec_annotation", "tables", c("Endothelial_subcluster_barcode_annotations.csv", "ec_subtype_proportions_per_sample.csv", "ec_propeller_results.csv")),
             RES("04_cm_annotation", "tables", c("Cardiomyocyte_subcluster_barcode_annotations.csv", "cm_subtype_proportions_per_sample.csv", "cm_propeller_results.csv"))),
  fig3   = RES("05_ec_de", c("ec_DE_all_MAST_annotated.csv", "ec_DE_target_genes_annotated.csv", "ec_DE_all_summary_annotated.csv", "ec_DE_per_animal_expression.csv")),
  fig4   = c(RESOURCE("mouse_secretome_swissprot.csv"), RES("05b_ec_pseudobulk", c("ec_pseudobulk_DE.csv", "ec_pseudobulk_vs_mast_genes.csv"))),
  fig5   = RES("06_ec_gsea", "ec_gsea_pathway_panel.csv"),
  fig6   = c(RES("07_ccc", "cellchat_list.rds"), RES("07_ccc", "tables", "ccc_diff_edge_decomposition.csv")),
  fig7   = c(RES("07_ccc", "tables", "lr_diff_all.csv"), RES("07_ccc", "lr_scores_by_sample.csv"),
             RES("07b_ec_cm_subtype_ccc", "tables", c("lr_subtype_diff_all.csv", "subtype_animal_support.csv"))))
invisible(lapply(IN, input_stamp))

animals <- read.csv(RES("03_ec_annotation", "tables", "ec_sample_metadata.csv")) |>
  mutate(group = factor(group, levels = GRP_ORDER)) |> arrange(group, sample)
stopifnot(nrow(animals) == 16, !anyDuplicated(animals$sample))
removed_barcodes <- data.table::fread(RES("02_subclustering", "atlas_cleaning_removed_barcodes.csv"))
FLAG_COL  <- c(none = "grey25", `single-animal driven` = "#E69F00", `low detection` = "#CC79A7", both = "#D55E00")
HIT_CLASS <- c("not significant" = "grey85", "technical gene" = "grey40", "sex-chromosome gene" = "#8C564B",
               "low detection" = FLAG_COL[["low detection"]], "single-animal driven" = FLAG_COL[["single-animal driven"]],
               "reported hit" = "#0072B2")
flag_mark <- function(sig, sad, low) paste0(ifelse(sig & sad %in% TRUE, "\u2020", ""), ifelse(sig & low %in% TRUE, "\u00b0", ""))
short_cm  <- function(x) sub(" CM$", "", x)
wrap_lab  <- function(w) function(x) vapply(x, function(s) paste(strwrap(s, w), collapse = "\n"), character(1))

FLAG_SHORT <- c(low_yield_library = "low yield", manual_doublet_threshold = "manual doublet thr.")
ANIMAL_KEY <- c(setNames(sprintf("%s (animal note)", names(ANIMAL_NOTES)), names(ANIMAL_NOTES)),
                setNames(sprintf("%s (%s)", names(FLAG_OF), FLAG_SHORT[FLAG_OF]), names(FLAG_OF)), other = "other animals")
ANIMAL_SHAPE <- setNames(c(17, c(1, 0, 5)[seq_along(FLAG_OF)], 16), ANIMAL_KEY)
animal_class <- function(s) factor(unname(ifelse(s %in% names(ANIMAL_KEY), ANIMAL_KEY[s], ANIMAL_KEY[["other"]])), levels = ANIMAL_KEY)
CONTRAST_STRIP <- CONTRAST_LABEL; CONTRAST_STRIP[["Interaction"]] <- "Age x HFpEF\n(interaction)"
INTERACTION_NOTE <- "Age x HFpEF: interaction effect (MAST, log2 scale)"

t0 <- Sys.time()
atlas <- readRDS(OBJ("integrated_LV_snRNAseq_seurat_clean.rds"))
meta1 <- atlas@meta.data |> mutate(cell_type = as.character(cell_type))
ct1 <- intersect(CT_ORDER, unique(meta1$cell_type))
stopifnot(setequal(ct1, unique(meta1$cell_type)))
um1 <- data.frame(cell_type = factor(meta1$cell_type, levels = ct1),
                  setNames(as.data.frame(Embeddings(atlas, "umap")[, 1:2]), c("UMAP_1", "UMAP_2")))

MARKERS_F1 <- list(Cardiomyocyte = c("Mhrt", "Atp2a2", "Nppa"), Endothelial = c("Flt1", "Kdr", "Npr3"),
                   Fibroblast = c("Abca8a", "Col3a1", "Lama2"), Pericyte = c("Notch3", "Vtn", "Higd1b"),
                   VSMC = c("Flna", "Lmod1", "Mylk"), Macrophage = c("Ptprc", "Lyz2", "Adgre1"),
                   Dendritic_cell = c("Flt3", "Xcr1", "Clec9a"), Lymphocyte = c("Cd247", "Bcl11b", "Prkcq"),
                   Lymphatic_EC = c("Ccl21a", "Lyve1", "Pdpn"), Glia = c("Sox10", "Adam23", "Gfra3"),
                   Epicardial = c("Muc16", "Upk1b", "Upk3b"), Proliferating = c("Hjurp", "Knl1", "Cenpf"),
                   B_cell = c("Ighm", "Bank1", "Fcmr"))
stopifnot(setequal(names(MARKERS_F1), ct1))
MARKERS_F1 <- lapply(MARKERS_F1[ct1], intersect, y = rownames(atlas))
genes1 <- unlist(MARKERS_F1, use.names = FALSE)
X1  <- GetAssayData(atlas, assay = "RNA", layer = "data")[genes1, , drop = FALSE]
Fm1 <- fac2sparse(factor(meta1$cell_type, levels = ct1))
n_ct1 <- rowSums(Fm1)
avg1 <- as.matrix(Fm1 %*% t(X1)) / n_ct1
pct1 <- as.matrix(Fm1 %*% t(X1 > 0)) / n_ct1
rm(atlas, X1, Fm1); invisible(gc())
dot1 <- as.data.frame(as.table(avg1), stringsAsFactors = FALSE) |> setNames(c("cell_type", "gene", "mean")) |>
  left_join(as.data.frame(as.table(pct1), stringsAsFactors = FALSE) |> setNames(c("cell_type", "gene", "frac")), by = c("cell_type", "gene")) |>
  group_by(gene) |> mutate(scaled = (mean - min(mean)) / (max(mean) - min(mean))) |> ungroup() |>
  mutate(panel = factor(rep(names(MARKERS_F1), lengths(MARKERS_F1))[match(gene, genes1)], levels = names(MARKERS_F1)),
         gene = factor(gene, levels = genes1), cell_type = factor(cell_type, levels = rev(ct1)))

comp1 <- meta1 |> count(individual, disease, qc_flag, cell_type) |>
  complete(nesting(individual, disease, qc_flag), cell_type = ct1, fill = list(n = 0)) |>
  group_by(individual) |> mutate(prop = n / sum(n)) |> ungroup() |>
  mutate(cell_type = factor(cell_type, levels = ct1), group = factor(disease, levels = GRP_ORDER))
design1 <- meta1 |> distinct(individual, disease, sex) |> count(disease, sex) |> pivot_wider(names_from = sex, values_from = n, values_fill = 0)
nuclei_clean_by_animal <- meta1 |> count(individual, name = "nuclei_after_cleaning")
n_ct_tab <- table(factor(meta1$cell_type, levels = ct1))
n_clean  <- nrow(meta1)
qc_ret   <- read.csv(RES("01_qc_atlas", "tables", "qc_retention_by_sample.csv"))
clean_by_ct <- read.csv(RES("02_subclustering", "atlas_cleaning_by_cell_type.csv"))
n_preqc <- sum(qc_ret$n_preqc); n_qc <- sum(qc_ret$n_retained)
stopifnot(sum(clean_by_ct$nuclei) == n_qc, sum(clean_by_ct$retained) == n_clean)
rm(meta1); invisible(gc())

des <- data.frame(group = GRP_ORDER) |>
  mutate(age = ifelse(grepl("^Y", group), "young", "old"), cond = ifelse(grepl("C$", group), "control", "HFpEF")) |>
  left_join(design1, by = c(group = "disease")) |>
  mutate(x = ifelse(age == "young", 1, 2), y = ifelse(cond == "control", 2, 1), n = F + M,
         lab = sprintf("%s  n = %d\n%d F / %d M", group, n, F, M), txt = ifelse(group %in% c("OC", "OH"), "white", "black"))
stopifnot(sum(des$n) == 16)
steps <- c(sprintf("Cell Ranger BAMs: %d libraries, 1 per animal", nrow(qc_ret)),
           "STARsolo recount, intron-inclusive (GeneFull)",
           sprintf("CellBender (ambient RNA): %s droplets", fmt_n(n_preqc)),
           sprintf("Nucleus QC + Scrublet: %s nuclei", fmt_n(n_qc)),
           "Harmony integration; broad cell types",
           sprintf("Sub-cluster cleaning: %s nuclei", fmt_n(n_clean)),
           "(dendritic and epicardial cells relabelled)",
           "EC / CM subtypes; animal-level composition",
           "MAST (nuclei); DESeq2 pseudobulk (animals)",
           "GSEA (MAST ranking)",
           "CellChat (pooled) + per-animal LR model")
arrow_after <- c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE)
wf <- data.frame(step = steps, y = rev(seq_along(steps)), arrow = arrow_after, indent = c(rep(0, 6), 1, rep(0, 4)))
pA_des <- ggplot(des) +
  geom_tile(aes(x, y, fill = group), width = .95, height = .9, colour = NA) +
  geom_text(aes(x, y, label = lab, colour = txt), size = pt_text(), lineheight = .9) +
  annotate("text", x = c(1, 2), y = 2.66, label = c("young (Y)", "old (O)"), size = pt_text(), fontface = "bold") +
  annotate("text", x = .47, y = c(2, 1), label = c("control (C)", "HFpEF (H)"), size = pt_text(), fontface = "bold", hjust = 1) +
  scale_fill_manual(values = GRP_FILL, guide = "none") + scale_colour_identity() +
  coord_cartesian(xlim = c(.5, 2.5), ylim = c(.5, 2.82), expand = FALSE, clip = "off") +
  theme_void(base_size = BASE_PT) + theme(plot.margin = margin(1, 8, 1, 22, "mm"))
pA_wf <- ggplot(wf, aes(0, y)) +
  geom_segment(data = wf[wf$arrow & wf$y > 1, ], aes(x = 0.15, xend = 0.15, y = y - .2, yend = y - .8),
               arrow = arrow(length = unit(0.9, "mm"), type = "closed"), linewidth = .25, colour = "grey45") +
  geom_text(aes(x = 0.45 + indent * 0.3, label = step), hjust = 0, size = pt_text()) +
  coord_cartesian(xlim = c(0, 10), ylim = c(.5, nrow(wf) + .5), expand = FALSE, clip = "off") +
  theme_void(base_size = BASE_PT) + theme(plot.margin = margin(1, 1, 1, 1, "mm"))
panelA1 <- wrap_elements(full = (wrap_elements(full = pA_des) / wrap_elements(full = pA_wf)) + plot_layout(heights = c(27, 55)))

bw  <- diff(range(um1$UMAP_1)) / 40
lab1 <- um1 |> mutate(bx = round(UMAP_1 / bw), by = round(UMAP_2 / bw)) |> count(cell_type, bx, by) |>
  group_by(cell_type) |> slice_max(n, n = 1, with_ties = FALSE) |> ungroup() |>
  mutate(UMAP_1 = bx * bw, UMAP_2 = by * bw, label = CT_LABEL[as.character(cell_type)])
pB1 <- ggplot(um1, aes(UMAP_1, UMAP_2, colour = cell_type)) +
  scattermore::geom_scattermore(pointsize = 1.6, pixels = raster_px(95, 80)) +
  ggrepel::geom_text_repel(data = lab1, aes(label = label), colour = "black", size = pt_text(), seed = 1,
                           box.padding = .3, min.segment.length = 0, segment.size = .25, segment.colour = "grey30") +
  scale_colour_manual(values = CT_FILL[ct1], guide = "none") +
  scale_x_continuous(expand = expansion(mult = 0.04)) + scale_y_continuous(expand = expansion(mult = 0.04)) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_text(size = BASE_PT))

STRIP_C1 <- CT_SHORT
pC1 <- ggplot(dot1, aes(gene, cell_type)) +
  geom_point(aes(size = 100 * frac, colour = scaled)) +
  facet_grid(. ~ panel, scales = "free_x", space = "free_x", labeller = as_labeller(STRIP_C1)) +
  scale_size_area(max_size = 2.8, breaks = c(25, 50, 75, 100), name = "nuclei\ndetecting (%)") +
  scale_colour_distiller(palette = "Reds", direction = 1, limits = c(0, 1), breaks = c(0, .5, 1), name = "mean expr.\n(scaled 0-1)") +
  scale_y_discrete(labels = CT_LABEL) +
  labs(x = NULL, y = NULL) +
  guides(colour = guide_colourbar(barwidth = unit(2, "mm"), barheight = unit(10, "mm"), order = 1)) +
  theme_fig(grid = "both") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, face = "italic"),
        panel.spacing.x = unit(.5, "mm"), strip.text = element_text(size = BASE_PT, face = "bold", vjust = 0),
        strip.clip = "off", legend.box = "horizontal", legend.spacing.x = unit(1, "mm"))

w1 <- .19
cd1 <- comp1 |> mutate(xc = as.numeric(cell_type) + (as.numeric(group) - 2.5) * w1,
                       xj = xc + withr::with_seed(1, runif(n(), -.045, .045)),
                       flagged = ifelse(qc_flag == "pass", "pass", "library-quality flag"))
med1 <- cd1 |> group_by(cell_type, group, xc) |> summarise(m = median(prop), .groups = "drop")
zero1 <- cd1 |> filter(prop == 0)
ymin1 <- min(cd1$prop[cd1$prop > 0]) * .7
pD1 <- ggplot(cd1 |> filter(prop > 0), aes(xj, 100 * prop)) +
  geom_point(aes(fill = group, colour = flagged), shape = 21, size = 1.2, stroke = .35) +
  geom_segment(data = med1, aes(x = xc - w1 * .4, xend = xc + w1 * .4, y = 100 * m, yend = 100 * m), inherit.aes = FALSE, linewidth = .45) +
  scale_y_log10(limits = c(100 * ymin1, 100), breaks = c(.01, .1, 1, 10, 100), labels = c("0.01", "0.1", "1", "10", "100")) +
  scale_x_continuous(breaks = seq_along(ct1), labels = CT_LABEL[ct1], expand = expansion(add = .4)) +
  scale_fill_manual(values = GRP_FILL, labels = GRP_LABEL[GRP_ORDER], name = NULL) +
  scale_colour_manual(values = c(pass = "transparent", `library-quality flag` = "black"), breaks = "library-quality flag", name = NULL) +
  guides(fill = guide_legend(override.aes = list(size = 1.8, colour = NA), nrow = 1, order = 1),
         colour = guide_legend(override.aes = list(fill = "white", size = 1.8))) +
  labs(x = NULL, y = "nuclei of the animal (%)") +
  theme_fig(grid = "y") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "top", legend.justification = "right",
        legend.margin = margin(0, 0, -2, 0, "mm"), legend.key.width = unit(2.5, "mm"))

top1 <- (panelA1 | pB1) + plot_layout(widths = c(85, 95))
fig1 <- top1 / pC1 / pD1 + plot_layout(heights = c(82, 50, 50))
figlog_1 <- save_figure(fig1, "Figure_1", WIDTH_2COL, 196)
rm(fig1, top1, pB1, um1); invisible(gc())
L1 <- list(n_preqc = n_preqc, n_qc = n_qc, n_clean = n_clean, n_ct = n_ct_tab, design = design1, n_zero = nrow(zero1),
           min_prop = min(cd1$prop[cd1$prop > 0]), n_types = length(ct1), markers = MARKERS_F1,
           removed_pct = 100 * (1 - n_clean / n_qc))
time_main_fig1 <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

qc_flags_01 <- read.csv(RES("01_qc_atlas", "tables", "sample_qc_flags.csv"))
cm_meta <- read.csv(RES("04_cm_annotation", "tables", "cm_sample_metadata.csv"))
animal_qc7 <- read.csv(RES("07_ccc", "tables", "animal_qc.csv"))
table1 <- animals |>
  transmute(animal = sample, group = as.character(group), age, condition, sex, batch, sequencing_run = seq) |>
  left_join(qc_ret |> select(animal = sample, droplets_cellbender = n_preqc, nuclei_after_qc = n_retained), by = "animal") |>
  left_join(nuclei_clean_by_animal |> rename(animal = individual), by = "animal") |>
  left_join(animals |> select(animal = sample, EC_nuclei = n_nuclei), by = "animal") |>
  left_join(cm_meta |> select(animal = sample, CM_nuclei = n_nuclei), by = "animal") |>
  left_join(animal_qc7 |> select(animal = sample, EC_median_exon_prop = ec_median_exon_prop, CM_median_exon_prop = cm_median_exon_prop), by = "animal") |>
  mutate(across(ends_with("exon_prop"), ~ round(.x, 3))) |>
  left_join(qc_flags_01 |> select(animal = sample, qc_flag), by = "animal") |>
  mutate(animal_note = ifelse(animal %in% names(ANIMAL_NOTES), sub(";.*$", "", ANIMAL_NOTES[animal]), ""),
         removed_at_qc_pct = round(100 * (1 - nuclei_after_qc / droplets_cellbender), 1),
         removed_at_cleaning_pct = round(100 * (1 - nuclei_after_cleaning / nuclei_after_qc), 1)) |>
  relocate(removed_at_qc_pct, .after = nuclei_after_qc) |> relocate(removed_at_cleaning_pct, .after = nuclei_after_cleaning)
stopifnot(nrow(table1) == 16, !anyNA(table1$nuclei_after_cleaning), sum(table1$nuclei_after_cleaning) == L1$n_clean,
          identical(table1$qc_flag, animals$qc_flag[match(table1$animal, animals$sample)]))
write_table(table1, "Table_1_animals", "animals")
knitr::kable(table1, caption = "Table 1 (tables/Table_1_animals.csv)")

t0 <- Sys.time()
LINEAGES <- list(
  EC = list(rds = RES("02_subclustering", "Endothelial_subcluster.rds"),
            ann = RES("03_ec_annotation", "tables", "Endothelial_subcluster_barcode_annotations.csv"),
            props = RES("03_ec_annotation", "tables", "ec_subtype_proportions_per_sample.csv"),
            res = RES("03_ec_annotation", "tables", "ec_propeller_results.csv"),
            label_col = "EC_label", kind_col = "EC_kind", subtype_col = "EC_subtype", short = short_subtype, lineage = "EC"),
  CM = list(rds = RES("02_subclustering", "Cardiomyocyte_subcluster.rds"),
            ann = RES("04_cm_annotation", "tables", "Cardiomyocyte_subcluster_barcode_annotations.csv"),
            props = RES("04_cm_annotation", "tables", "cm_subtype_proportions_per_sample.csv"),
            res = RES("04_cm_annotation", "tables", "cm_propeller_results.csv"),
            label_col = "CM_label", kind_col = "CM_kind", subtype_col = "CM_subtype", short = short_cm, lineage = "CM"))

build_lineage <- function(L) {
  res <- read.csv(L$res)
  subtypes <- unique(res$subtype)
  cols  <- subtype_palette(subtypes)
  tcols <- setNames(colorspace::darken(cols, 0.25), subtypes)
  ann <- read.csv(L$ann)
  ann$label <- ann[[L$label_col]]; ann$kind <- ann[[L$kind_col]]
  ann$removed <- ann$barcode %in% removed_barcodes$barcode
  stopifnot(all(ann$label[!ann$removed & ann$kind == "subtype"] %in% subtypes))
  parent_of <- function(l) vapply(l, function(x) { p <- subtypes[startsWith(x, subtypes)]; if (length(p)) p[1] else NA_character_ }, "")
  ann <- ann |> mutate(parent = ifelse(kind == "subtype", label, parent_of(label)),
                       stratum = case_when(removed ~ "removed", kind == "subtype" ~ label, !is.na(parent) ~ paste0("merged:", parent), TRUE ~ "removed"))
  stopifnot(!any(ann$stratum == "removed" & !ann$removed))
  cl_ids <- function(sel) paste(sort(unique(ann$cluster[sel])), collapse = ", ")
  merged_parents <- unique(ann$parent[grepl("^merged:", ann$stratum)])
  keys <- c()
  for (p in merged_parents) keys[paste0("merged:", p)] <- sprintf("merged into %s (cluster%s %s)", p,
                                                                  ifelse(length(unique(ann$cluster[ann$stratum == paste0("merged:", p)])) > 1, "s", ""),
                                                                  cl_ids(ann$stratum == paste0("merged:", p)))
  if (any(ann$removed)) keys["removed"] <- sprintf("removed in atlas cleaning (cluster%s %s)",
                                                   ifelse(length(unique(ann$cluster[ann$removed])) > 1, "s", ""), cl_ids(ann$removed))
  key_cols <- c(setNames(colorspace::lighten(cols[merged_parents], 0.6), paste0("merged:", merged_parents)), removed = TECH_GREY)[names(keys)]

  obj <- readRDS(L$rds)
  emb <- as.data.frame(Embeddings(obj, "umap")[, 1:2]) |> setNames(c("UMAP_1", "UMAP_2"))
  rm(obj); invisible(gc())
  stopifnot(all(ann$barcode %in% rownames(emb)))
  um <- cbind(emb[ann$barcode, ], ann[, c("label", "stratum", "cluster")])
  um <- um[order(um$stratum %in% subtypes), ]
  um_lab <- um |> filter(stratum %in% subtypes) |> group_by(stratum) |>
    summarise(across(c(UMAP_1, UMAP_2), median), .groups = "drop") |> mutate(label = L$short(stratum))
  p_umap <- ggplot(um, aes(UMAP_1, UMAP_2, colour = stratum)) +
    scattermore::geom_scattermore(pointsize = 1.8, pixels = raster_px(70, 60)) +
    ggrepel::geom_text_repel(data = um_lab, aes(label = label), colour = "black", size = pt_text(), seed = 1,
                             min.segment.length = Inf, box.padding = .15) +
    scale_colour_manual(values = c(cols, key_cols), breaks = names(keys), labels = unname(keys), name = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 1.8, pointsize = 3))) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_text(size = BASE_PT),
          legend.position = "top", legend.justification = "right", legend.direction = "vertical",
          legend.key.size = unit(2.5, "mm"), legend.margin = margin(0, 0, 0, 0, "mm"))

  one_animal <- res |> distinct(subtype, max_animal, max_animal_share_prop) |> mutate(one_animal = max_animal_share_prop > 0.5)
  strip <- setNames(paste0(L$short(subtypes), ifelse(one_animal$one_animal[match(subtypes, one_animal$subtype)], " \u2020", "")), subtypes)
  ps <- read.csv(L$props) |> rename(subtype = !!L$subtype_col) |>
    mutate(subtype = factor(subtype, levels = subtypes), group = factor(group, levels = GRP_ORDER), animal = animal_class(sample),
           x = as.numeric(group), x_j = x + withr::with_seed(1, runif(n(), -.14, .14)))
  stopifnot(!anyNA(ps$subtype))
  p_comp <- ggplot(ps, aes(y = 100 * proportion)) +
    stat_summary(aes(x = x, fill = group), fun = mean, geom = "col", width = .72) +
    geom_point(aes(x = x_j, shape = animal), size = .9, stroke = .4, colour = "grey10") +
    facet_wrap(~ subtype, scales = "free_y", nrow = 1, labeller = as_labeller(wrap_lab(16)(strip))) +
    scale_x_continuous(breaks = seq_along(GRP_ORDER), labels = GRP_ORDER) +
    scale_y_continuous(expand = expansion(mult = c(0, .06))) +
    scale_fill_manual(values = GRP_FILL, labels = GRP_LABEL, name = NULL) +
    scale_shape_manual(values = ANIMAL_SHAPE, name = NULL, drop = FALSE) +
    labs(x = NULL, y = sprintf("%s nuclei of the animal (%%)", L$lineage)) +
    theme_fig(grid = "y") +
    theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1), panel.spacing.x = unit(1.5, "mm"),
          strip.text = element_text(size = BASE_PT, face = "bold", lineheight = .9), strip.clip = "off")

  eff <- res |> mutate(subtype = factor(subtype, levels = rev(subtypes)),
                       contrast = factor(CONTRAST_LABEL[contrast], levels = unname(CONTRAST_LABEL)), star = star_fdr(FDR))
  p_eff <- ggplot(eff, aes(logit_effect, subtype, colour = subtype)) +
    geom_vline(xintercept = 0, colour = "grey60", linewidth = .3) +
    geom_errorbar(aes(xmin = CI_low, xmax = CI_high), width = .3, orientation = "y", linewidth = .35) +
    geom_point(size = 1.2) +
    geom_text(aes(label = star), colour = "black", nudge_y = .36, fontface = "bold", size = pt_text(TITLE_PT)) +
    facet_wrap(~ contrast, nrow = 1) +
    scale_colour_manual(values = tcols, guide = "none") +
    scale_y_discrete(labels = function(x) strip[x]) +
    scale_x_continuous(breaks = function(lim) pretty(lim, n = 4)) +
    labs(x = "logit-proportion effect (95% CI)", y = NULL) +
    theme_fig(grid = "x") + theme(panel.spacing.x = unit(2.5, "mm"))
  list(umap = p_umap, comp = p_comp, eff = p_eff, res = res, one_animal = one_animal, ps = ps, subtypes = subtypes,
       n_label = table(ann$label), n_stratum = table(ann$stratum), keys = keys, n_total = nrow(ann),
       n_retained = sum(!ann$removed), n_subtype = table(ann$parent[!ann$removed]))
}
ec2 <- build_lineage(LINEAGES$EC)
cm2 <- build_lineage(LINEAGES$CM)

leg_theme <- theme(legend.position = "bottom", legend.box = "vertical", legend.spacing.y = unit(0, "mm"),
                   legend.margin = margin(0, 0, 0, 0), legend.key.size = unit(2.6, "mm"), legend.box.just = "left")
pB2 <- ec2$comp + guides(fill = guide_legend(nrow = 1, order = 1),
                         shape = guide_legend(ncol = 2, byrow = TRUE, order = 2, override.aes = list(size = 1.4))) + leg_theme
pE2 <- cm2$comp + theme(legend.position = "none")
row1 <- (wrap_elements(full = ec2$umap) | pB2) + plot_layout(widths = c(70, 110))
row3 <- (wrap_elements(full = cm2$umap) | pE2) + plot_layout(widths = c(70, 110))
fig2 <- row1 / ec2$eff / row3 / cm2$eff + plot_layout(heights = c(78, 34, 66, 32))
figlog_2 <- save_figure(fig2, "Figure_2", WIDTH_2COL, 226)
rm(fig2, row1, row3, pB2, pE2); ec2$umap <- NULL; cm2$umap <- NULL; invisible(gc())
time_main_fig2 <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

table2 <- bind_rows(lapply(list(EC = ec2, CM = cm2), function(x) x$res), .id = "lineage") |>
  mutate(contrast_label = unname(CONTRAST_LABEL[contrast]),
         one_animal_subtype = max_animal_share_prop > 0.5,
         nuclei = as.integer(c(ec2$n_subtype, cm2$n_subtype)[subtype])) |>
  transmute(lineage, subtype, nuclei, contrast, contrast_label, logit_effect, SE, CI_low_95 = CI_low, CI_high_95 = CI_high,
            p_value, FDR_within_contrast = FDR, FDR_across_contrasts = FDR_all_contrasts, n_animals, n_animals_present,
            max_animal, max_animal_share_prop, one_animal_subtype) |>
  mutate(across(c(logit_effect, SE, CI_low_95, CI_high_95, p_value, FDR_within_contrast, FDR_across_contrasts, max_animal_share_prop), ~ signif(.x, 4)))
stopifnot(nrow(table2) == nrow(ec2$res) + nrow(cm2$res), !anyNA(table2$nuclei))
write_table(table2, "Table_2_composition_tests", "composition")
knitr::kable(table2 |> filter(FDR_within_contrast < 0.1) |> mutate(across(where(is.numeric), ~ signif(.x, 3))),
             caption = "Table 2 rows with FDR within contrast < 0.1 (full table: tables/Table_2_composition_tests.csv)")

t0 <- Sys.time()
PADJ_CUT <- 0.05; LFC_CUT <- 0.2; VOLCANO_CT <- "Capillary EC"
de_all  <- readr::read_csv(RES("05_ec_de", "ec_DE_all_MAST_annotated.csv"), show_col_types = FALSE) |> mutate(gene_flag = coalesce(gene_flag, ""))
tgt_all <- readr::read_csv(RES("05_ec_de", "ec_DE_target_genes_annotated.csv"), show_col_types = FALSE) |> mutate(gene_flag = coalesce(gene_flag, ""))
summ5   <- readr::read_csv(RES("05_ec_de", "ec_DE_all_summary_annotated.csv"), show_col_types = FALSE)
ae5     <- readr::read_csv(RES("05_ec_de", "ec_DE_per_animal_expression.csv"), show_col_types = FALSE)
EC_SUBTYPES <- ec2$subtypes
CANDIDATES  <- unique(tgt_all$gene)
N_LABEL3 <- 5; N_TOP3 <- 2
stopifnot(setequal(unique(de_all$celltype), EC_SUBTYPES))
code_key <- paste(sprintf("%s = %s", CONTRAST_CODE, CONTRAST_LABEL), collapse = ";  ")

HIT_SEG <- c("single-animal driven", "low detection", "other reported hits")
hits3 <- de_all |> filter(reported) |>
  mutate(class = case_when(single_animal_driven ~ HIT_SEG[1], low_detection ~ HIT_SEG[2], TRUE ~ HIT_SEG[3])) |>
  count(celltype, contrast, class) |>
  complete(celltype = EC_SUBTYPES, contrast = CO_ORDER, class = HIT_SEG, fill = list(n = 0))
tot3 <- hits3 |> group_by(celltype, contrast) |> summarise(n = sum(n), .groups = "drop")
chk3 <- tot3 |> inner_join(summ5 |> select(celltype, contrast, n_hits_reported), by = c("celltype", "contrast"))
stopifnot(nrow(chk3) == 25, all(chk3$n == chk3$n_hits_reported))
hits3 <- hits3 |> mutate(celltype = factor(celltype, levels = EC_SUBTYPES), contrast = factor(contrast, levels = rev(CO_ORDER)),
                         class = factor(class, levels = rev(HIT_SEG)))
tot3 <- tot3 |> mutate(celltype = factor(celltype, levels = EC_SUBTYPES), contrast = factor(contrast, levels = rev(CO_ORDER)))
pA3 <- ggplot(hits3, aes(n, contrast)) +
  geom_col(aes(fill = class), width = .72) +
  geom_text(data = tot3, aes(label = fmt_n(n)), hjust = -.15, size = pt_text()) +
  facet_wrap(~ celltype, nrow = 1, scales = "free_x", labeller = as_labeller(function(x) short_subtype(x))) +
  scale_fill_manual(values = c(`other reported hits` = HIT_CLASS[["reported hit"]], `low detection` = HIT_CLASS[["low detection"]],
                               `single-animal driven` = HIT_CLASS[["single-animal driven"]]), breaks = rev(HIT_SEG), name = NULL) +
  scale_y_discrete(labels = CONTRAST_LABEL) +
  scale_x_continuous(expand = expansion(mult = c(0, .38)), breaks = scales::breaks_extended(3)) +
  coord_cartesian(clip = "off") +
  labs(x = "reported MAST hits (genes)", y = NULL) +
  theme_fig(grid = "x") +
  theme(legend.position = "top", legend.justification = "right", legend.margin = margin(0, 0, -2, 0, "mm"),
        panel.spacing.x = unit(3, "mm"))

volc3 <- de_all |> filter(celltype == VOLCANO_CT) |>
  mutate(class = factor(case_when(!sig ~ "not significant", gene_flag == "technical" ~ "technical gene",
                                  gene_flag == "sex" ~ "sex-chromosome gene", single_animal_driven ~ "single-animal driven",
                                  low_detection ~ "low detection", TRUE ~ "reported hit"), levels = names(HIT_CLASS)),
         nlp = -log10(pmax(p_val_adj, 1e-300)), contrast = factor(contrast, levels = CO_ORDER)) |>
  arrange(class)
labs3 <- volc3 |> filter(reported, !gene %in% CANDIDATES) |> group_by(contrast) |> slice_min(p_val_adj, n = N_LABEL3, with_ties = FALSE) |> ungroup()

Y_MAX3 <- max(volc3$nlp); Y_BASE3 <- 1.06 * Y_MAX3; Y_STEP3 <- 0.085 * Y_MAX3
xr3 <- volc3 |> group_by(contrast) |> summarise(xmin = min(avg_log2FC), xmax = max(avg_log2FC), .groups = "drop")
lab3 <- labs3 |> left_join(xr3, by = "contrast") |> mutate(side = ifelse(avg_log2FC < 0, "left", "right")) |>
  arrange(contrast, desc(side == "left"), ifelse(side == "left", avg_log2FC, -avg_log2FC)) |>
  group_by(contrast) |>
  mutate(y_lab = Y_BASE3 + (n() - row_number()) * Y_STEP3, x_lab = ifelse(side == "left", xmin, xmax),
         hjust = ifelse(side == "left", 0, 1)) |> ungroup()
pB3 <- ggplot(volc3, aes(avg_log2FC, nlp, colour = class)) +
  geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), colour = "grey60", linetype = "dotted", linewidth = .25) +
  geom_hline(yintercept = -log10(PADJ_CUT), colour = "grey60", linetype = "dotted", linewidth = .25) +
  scattermore::geom_scattermore(data = filter(volc3, class == "not significant"), pointsize = 1.6, pixels = raster_px(34, 45)) +
  geom_point(data = filter(volc3, class != "not significant"), size = .45, stroke = 0) +
  geom_segment(data = lab3, aes(x = x_lab + ifelse(side == "left", 1, -1) * 0.03 * (xmax - xmin), xend = avg_log2FC,
                                y = y_lab, yend = nlp), colour = "grey60", linewidth = .2, inherit.aes = FALSE) +
  geom_label(data = lab3, aes(x = x_lab, y = y_lab, label = gene, hjust = hjust), size = pt_text(MIN_PT), fontface = "italic",
             colour = "black", fill = "white", border.colour = NA, label.padding = unit(0.25, "mm"), vjust = 0.3,
             inherit.aes = FALSE) +
  scale_y_continuous(breaks = pretty(c(0, Y_MAX3)), expand = expansion(mult = c(.02, .04))) +
  facet_wrap(~ contrast, nrow = 1, labeller = as_labeller(CONTRAST_STRIP), scales = "free_x") +
  scale_x_continuous(breaks = function(lim) { b <- pretty(lim, n = 3); b[b >= lim[1] & b <= lim[2]] }) +
  scale_colour_manual(values = HIT_CLASS, name = NULL, drop = FALSE,
                      labels = c("not significant", "technical", "sex chromosome", "low detection", "single-animal driven", "reported hit")) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 1.6))) +
  labs(x = "effect: avg_log2FC (pairwise contrasts) | interaction_log2 (Age x HFpEF)", y = "\u2212log10 Bonferroni p") +
  theme(legend.position = "bottom", legend.margin = margin(-2, 0, 0, 0, "mm"), panel.spacing.x = unit(2.5, "mm"),
        strip.clip = "off", strip.text = element_text(size = BASE_PT, face = "bold", lineheight = .9))

top3 <- de_all |> filter(celltype == VOLCANO_CT, reported, !gene %in% CANDIDATES) |> mutate(contrast = factor(contrast, levels = CO_ORDER)) |>
  group_by(contrast) |> slice_min(p_val_adj, n = N_TOP3, with_ties = FALSE) |> ungroup() |> distinct(gene) |> pull(gene)
aeC3 <- ae5 |> filter(celltype == VOLCANO_CT, gene %in% top3) |>
  mutate(gene = factor(gene, levels = top3), disease = factor(disease, levels = GRP_ORDER), animal = animal_class(individual))
stopifnot(setequal(unique(as.character(aeC3$gene)), top3), all(table(aeC3$gene) == 16))
codes3 <- de_all |> filter(celltype == VOLCANO_CT, gene %in% top3, sig) |>
  mutate(contrast = factor(contrast, levels = CO_ORDER)) |> arrange(contrast) |> group_by(gene) |>
  summarise(lab = paste(strwrap(paste0(CONTRAST_CODE[as.character(contrast)], star_fdr(p_val_adj),
                                       flag_mark(sig, single_animal_driven, low_detection), collapse = " "), width = 16), collapse = "\n"))
strip3 <- setNames(top3, top3); strip3[codes3$gene] <- paste0(codes3$gene, "\n", codes3$lab)
pC3 <- ggplot(aeC3, aes(disease, mean_expr)) +
  geom_point(aes(fill = disease, shape = animal, size = n_nuclei), colour = "grey20", stroke = .3,
             position = position_jitter(width = .12, height = 0, seed = 1)) +
  facet_wrap(~ gene, scales = "free_y", nrow = 2, labeller = as_labeller(strip3)) +
  scale_x_discrete(labels = function(x) x) +
  scale_fill_manual(values = GRP_FILL, labels = GRP_LABEL, name = NULL) +
  scale_shape_manual(values = c(setNames(c(24, 1, 0, 5, 21), ANIMAL_KEY)), name = NULL, drop = FALSE) +
  scale_size_continuous(range = c(.7, 2.3), name = "nuclei", breaks = c(500, 2000, 4000)) +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 2), order = 1, nrow = 1),
         shape = guide_legend(order = 2, nrow = 2, byrow = TRUE, override.aes = list(size = 1.6, fill = "grey60")), size = guide_legend(order = 3, nrow = 1)) +
  labs(x = NULL, y = "mean log-normalised expression\nper animal", caption = code_key) +
  theme_fig(grid = "y") +
  theme(strip.text = element_text(size = BASE_PT, face = "plain", hjust = 0, lineheight = .9),
        axis.text.x = element_text(size = MIN_PT), legend.position = "bottom", legend.box = "vertical", legend.box.just = "left",
        legend.key.size = unit(2.6, "mm"), legend.spacing.y = unit(0.6, "mm"),
        legend.margin = margin(0, 0, 0, 0, "mm"), legend.key.spacing.x = unit(0.3, "mm"),
        plot.caption.position = "plot", panel.spacing.y = unit(1, "mm"), panel.spacing.x = unit(3, "mm"))

fig3 <- pA3 / pB3 / pC3 + plot_layout(heights = c(40, 70, 86))
figlog_3 <- save_figure(fig3, "Figure_3", WIDTH_2COL, 230)
rm(fig3, pB3); invisible(gc())
rep3 <- de_all |> filter(reported)
L3 <- list(tot = tot3 |> mutate(contrast = factor(as.character(contrast), levels = CO_ORDER)) |> arrange(celltype, contrast),
           n_rep = nrow(rep3), n_sad = sum(rep3$single_animal_driven), n_low = sum(rep3$low_detection),
           n_both = sum(rep3$single_animal_driven & rep3$low_detection), n_bkg = sum(rep3$background_suspect),
           sad_top_animal = rep3 |> filter(single_animal_driven) |> count(max_share_animal, sort = TRUE),
           n_tested_cap = table(volc3$contrast), top = top3, summ = summ5,
           n_cap_nuclei = range(aeC3$n_nuclei), n_floor = sum(volc3$p_val_adj < 1e-300))
time_main_fig3 <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

t0 <- Sys.time()
secretome <- readr::read_csv(RESOURCE("mouse_secretome_swissprot.csv"), show_col_types = FALSE)
pb5 <- readr::read_csv(RES("05b_ec_pseudobulk", "ec_pseudobulk_DE.csv"), show_col_types = FALSE)
vs5 <- readr::read_csv(RES("05b_ec_pseudobulk", "ec_pseudobulk_vs_mast_genes.csv"), show_col_types = FALSE)
TARGET_GENES <- unique(tgt_all$gene)
SEC_TARGETS  <- TARGET_GENES[TARGET_GENES %in% secretome$gene]
OTHER_TARGETS <- setdiff(TARGET_GENES, SEC_TARGETS)
sig_sec <- de_all |> filter(reported, gene %in% secretome$gene)
MIN_BLOCKS <- 4

cnt4 <- expand_grid(celltype = EC_SUBTYPES, contrast = CO_ORDER) |>
  left_join(sig_sec |> count(celltype, contrast, name = "n"), by = c("celltype", "contrast")) |>
  mutate(n = coalesce(n, 0L), celltype = factor(celltype, levels = rev(EC_SUBTYPES)), contrast = factor(contrast, levels = CO_ORDER))
pA4 <- ggplot(cnt4, aes(contrast, celltype, fill = n)) +
  geom_tile(colour = "white", linewidth = .5) +
  geom_text(aes(label = n, colour = n > max(n) * .55), size = pt_text()) +
  scale_fill_gradient(low = "grey96", high = "#54278f", name = "secreted\ngene hits", breaks = function(l) pretty(l, n = 2)) +
  scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "white"), guide = "none") +
  scale_x_discrete(labels = CONTRAST_CODE) + scale_y_discrete(labels = short_subtype) +
  guides(fill = guide_colourbar(barwidth = unit(2, "mm"), barheight = unit(12, "mm"))) +
  labs(x = NULL, y = NULL) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        legend.position = "bottom", legend.title = element_text(vjust = .9)) +
  guides(fill = guide_colourbar(barwidth = unit(12, "mm"), barheight = unit(2, "mm"), title = "hits", title.vjust = .9))

rec4 <- sig_sec |> count(gene, name = "n_blocks") |> arrange(desc(n_blocks), gene)
recB <- rec4 |> filter(n_blocks >= MIN_BLOCKS)
dB4 <- sig_sec |> semi_join(recB, by = "gene") |>
  mutate(gene = factor(gene, levels = rev(recB$gene)), celltype = factor(celltype, levels = EC_SUBTYPES),
         contrast = factor(contrast, levels = CO_ORDER), lfc = pmax(pmin(avg_log2FC, 1.5), -1.5),
         mark = flag_mark(TRUE, single_animal_driven, low_detection))
nB4 <- recB |> mutate(gene = factor(gene, levels = rev(recB$gene)))
pB4 <- ggplot(dB4, aes(contrast, gene, fill = lfc)) +
  geom_tile(colour = "white", linewidth = .4) +
  geom_text(aes(label = mark), size = pt_text(), colour = "black") +
  facet_grid(. ~ celltype, labeller = as_labeller(short_subtype)) +
  scale_fill_gradient2(low = "#2166ac", mid = "grey97", high = "#b2182b", midpoint = 0, limits = c(-1.5, 1.5), breaks = c(-1.5, 0, 1.5),
                       name = "effect\navg_log2FC |\ninteraction_log2\n(AxH; capped \u00b11.5)") +
  geom_vline(xintercept = length(CO_ORDER) - .5, colour = "grey45", linewidth = .25, linetype = "dashed") +
  scale_x_discrete(drop = FALSE, labels = CONTRAST_CODE) +
  scale_y_discrete(drop = FALSE, labels = function(g) sprintf("%s  (%d)", g, recB$n_blocks[match(g, recB$gene)])) +
  guides(fill = guide_colourbar(barwidth = unit(2, "mm"), barheight = unit(12, "mm"))) +
  labs(x = NULL, y = NULL) +
  theme_fig(grid = "both") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.y = element_text(face = "italic"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5), panel.spacing.x = unit(1.5, "mm"),
        strip.text = element_text(size = BASE_PT, face = "bold"), strip.clip = "off")

tgt4 <- tgt_all |>
  mutate(gene = factor(gene, levels = c(SEC_TARGETS[order(-rec4$n_blocks[match(SEC_TARGETS, rec4$gene)])], OTHER_TARGETS)),
         celltype = factor(celltype, levels = rev(EC_SUBTYPES)), contrast = factor(contrast, levels = CO_ORDER),
         tested = status == "tested",
         mark = ifelse(!tested, "\u2013", paste0(ifelse(sig %in% TRUE, star_fdr(p_val_adj), ""), flag_mark(sig %in% TRUE, single_animal_driven, low_detection))),
         lfc_capped = pmax(pmin(avg_log2FC, 1.2), -1.2))
gene_strip4 <- setNames(ifelse(levels(tgt4$gene) %in% SEC_TARGETS, paste0(levels(tgt4$gene), "\n "), paste0(levels(tgt4$gene), "\n(not in secretome)")), levels(tgt4$gene))
pC4 <- ggplot(tgt4, aes(contrast, celltype, fill = ifelse(tested, lfc_capped, NA))) +
  geom_tile(colour = "white", linewidth = .6) +
  geom_text(aes(label = mark, colour = ifelse(!tested, "na", ifelse(abs(lfc_capped) > .6, "w", "k"))), fontface = "bold", size = pt_text(), show.legend = FALSE) +
  facet_wrap(~ gene, nrow = 1, labeller = as_labeller(gene_strip4)) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, limits = c(-1.2, 1.2),
                       na.value = "grey90", breaks = c(-1.2, -0.6, 0, 0.6, 1.2), name = "effect\navg_log2FC |\ninteraction_log2\n(AxH; capped \u00b11.2)") +
  geom_vline(xintercept = length(CO_ORDER) - .5, colour = "grey30", linewidth = .35, linetype = "dashed") +
  scale_colour_manual(values = c(k = "grey15", w = "white", na = "grey55")) +
  scale_x_discrete(labels = c(CONTRAST_CODE[-length(CONTRAST_CODE)], Interaction = "AxH\n(int.)")) + scale_y_discrete(labels = short_subtype) +
  guides(fill = guide_colourbar(barwidth = unit(2, "mm"), barheight = unit(14, "mm"))) +
  labs(x = NULL, y = NULL) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        strip.text = element_text(face = "bold.italic", size = BASE_PT, lineheight = .9), strip.clip = "off", panel.spacing.x = unit(3, "mm"))

aeD4 <- ae5 |> filter(celltype == VOLCANO_CT, gene %in% TARGET_GENES) |>
  mutate(gene = factor(gene, levels = levels(tgt4$gene)), disease = factor(disease, levels = GRP_ORDER), animal = animal_class(individual))
stopifnot(all(table(aeD4$gene) == 16))
codes4 <- tgt_all |> filter(celltype == VOLCANO_CT, sig %in% TRUE) |>
  mutate(contrast = factor(contrast, levels = CO_ORDER)) |> arrange(contrast) |> group_by(gene) |>
  summarise(lab = paste(strwrap(paste0(CONTRAST_CODE[as.character(contrast)], star_fdr(p_val_adj),
                                       flag_mark(sig, single_animal_driven, low_detection), collapse = " "), width = 9), collapse = "\n"))
strip4 <- setNames(levels(tgt4$gene), levels(tgt4$gene)); strip4[codes4$gene] <- paste0(codes4$gene, "\n", codes4$lab)
pD4 <- ggplot(aeD4, aes(disease, mean_expr)) +
  geom_point(aes(fill = disease, shape = animal, size = n_nuclei), colour = "grey20", stroke = .3,
             position = position_jitter(width = .12, height = 0, seed = 1)) +
  facet_wrap(~ gene, scales = "free_y", nrow = 1, labeller = as_labeller(strip4)) +
  scale_fill_manual(values = GRP_FILL, labels = GRP_LABEL, name = NULL) +
  scale_shape_manual(values = c(setNames(c(24, 1, 0, 5, 21), ANIMAL_KEY)), name = NULL, drop = FALSE) +
  scale_size_continuous(range = c(.7, 2.3), name = "nuclei", breaks = c(500, 2000, 4000)) +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 2), order = 1, nrow = 1),
         shape = guide_legend(order = 2, ncol = 2, byrow = TRUE, override.aes = list(size = 1.6, fill = "grey60")), size = guide_legend(order = 3, nrow = 1)) +
  labs(x = NULL, y = "mean log-normalised\nexpression per animal") +
  theme_fig(grid = "y") +
  theme(strip.text = element_text(size = BASE_PT, face = "plain", hjust = 0, lineheight = .9),
        axis.text.x = element_text(size = MIN_PT, angle = 90, vjust = .5, hjust = 1), legend.position = "bottom", legend.box = "vertical", legend.box.just = "left",
        legend.key.size = unit(2.4, "mm"), legend.spacing.y = unit(0.6, "mm"), legend.margin = margin(0, 0, 0, 0, "mm"),
        panel.spacing.x = unit(3, "mm"))

agree4 <- tgt_all |> filter(status == "tested") |>
  select(celltype, contrast, gene, avg_log2FC, sig, reported) |>
  left_join(pb5 |> select(celltype, contrast, gene, logFC, FDR), by = c("celltype", "contrast", "gene")) |>
  mutate(gene = factor(gene, levels = levels(tgt4$gene)), contrast = factor(contrast, levels = CO_ORDER),
         pb_class = case_when(is.na(logFC) ~ "not in pseudobulk", !is.na(FDR) & FDR < 0.05 ~ "pseudobulk FDR < 0.05",
                              !is.na(FDR) & FDR < 0.1 ~ "pseudobulk FDR < 0.1", TRUE ~ "pseudobulk FDR \u2265 0.1"))
ag_lab <- agree4 |> filter(pb_class %in% c("pseudobulk FDR < 0.05", "pseudobulk FDR < 0.1")) |>
  mutate(lab = sprintf("%s %s", short_subtype(celltype), CONTRAST_CODE[as.character(contrast)]))

agr_ok <- filter(agree4, !is.na(logFC))
ag_lab05 <- ag_lab |> filter(pb_class == "pseudobulk FDR < 0.05") |> group_by(gene) |>
  mutate(x_lab = min(agr_ok$avg_log2FC), y_lab = max(agr_ok$logFC) - (row_number() - 1) * 0.12 * diff(range(agr_ok$logFC))) |> ungroup()
pE4 <- ggplot(filter(agree4, !is.na(logFC)), aes(avg_log2FC, logFC)) +
  geom_hline(yintercept = 0, colour = "grey75", linewidth = .25) + geom_vline(xintercept = 0, colour = "grey75", linewidth = .25) +
  geom_point(aes(fill = contrast, shape = sig %in% TRUE), size = 1.4, stroke = .3, colour = "grey20") +
  geom_point(data = filter(agree4, pb_class == "pseudobulk FDR < 0.05"), shape = 21, size = 3, stroke = .5, colour = "black", fill = NA) +
  geom_point(data = filter(agree4, pb_class == "pseudobulk FDR < 0.1"), shape = 21, size = 3, stroke = .35, colour = "grey45", fill = NA, linetype = "dashed") +
  geom_segment(data = ag_lab05, aes(x = x_lab + 0.1, y = y_lab - 0.08 * diff(range(agr_ok$logFC)), xend = avg_log2FC, yend = logFC),
               colour = "grey55", linewidth = .2) +
  geom_text(data = ag_lab05, aes(x = x_lab, y = y_lab, label = lab), size = pt_text(MIN_PT), hjust = 0, vjust = 1) +
  facet_wrap(~ gene, nrow = 2) +
  scale_x_continuous(breaks = function(lim) pretty(lim, n = 3)) + scale_y_continuous(breaks = function(lim) pretty(lim, n = 3)) +
  scale_fill_manual(values = CONTRAST_FILL, labels = CONTRAST_CODE, name = NULL) +
  scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 22), labels = c(`TRUE` = "MAST hit", `FALSE` = "not a MAST hit"), name = NULL) +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 1.8), nrow = 1, order = 1), shape = guide_legend(nrow = 1, order = 2, override.aes = list(size = 1.6))) +
  labs(x = "MAST avg_log2FC | interaction_log2 (AxH)", y = "pseudobulk log2FC (DESeq2)") +
  theme_fig(grid = "both") +
  theme(strip.text = element_text(face = "bold.italic", size = BASE_PT), legend.position = "bottom", legend.box = "vertical",
        legend.key.size = unit(2.4, "mm"), legend.spacing.y = unit(0.6, "mm"), legend.margin = margin(0, 0, 0, 0, "mm"),
        legend.key.spacing.x = unit(0.3, "mm"), panel.spacing = unit(2, "mm"))

row1_4 <- (wrap_elements(full = pA4) | wrap_elements(full = pB4)) + plot_layout(widths = c(44, 136))
row3_4 <- (wrap_elements(full = pD4) | wrap_elements(full = pE4)) + plot_layout(widths = c(104, 76))
fig4 <- (row1_4 / pC4 / row3_4 + plot_layout(heights = c(66, 40, 62))) +
  plot_annotation(caption = code_key, theme = theme(plot.caption = element_text(size = BASE_PT, hjust = 0, colour = "grey20")))
figlog_4 <- save_figure(fig4, "Figure_4", WIDTH_2COL, 212)
rm(fig4); invisible(gc())
L4 <- list(n_secretome = nrow(secretome), retrieved = as.character(secretome$retrieved_on[1]), n_sec_hits = nrow(sig_sec),
           n_sec_genes = n_distinct(sig_sec$gene), n_recurrent2 = sum(rec4$n_blocks >= 2), recB = recB, rec = rec4,
           sec_targets = SEC_TARGETS, other_targets = OTHER_TARGETS, cnt = cnt4,
           tgt_pass = tgt_all |> group_by(gene) |> summarise(tested = sum(status == "tested"), passing = sum(sig %in% TRUE), sad = sum(sig %in% TRUE & single_animal_driven %in% TRUE), .groups = "drop"),
           tgt_cap = tgt_all |> filter(celltype == VOLCANO_CT) |> select(gene, contrast, avg_log2FC, p_val_adj, sig, single_animal_driven, low_detection),
           agree = agree4, pb_n = pb5 |> distinct(celltype, n_samples, n_genes_tested),
           n_sec_hits_sad = sum(sig_sec$single_animal_driven),
           pb_mast = vs5 |> filter(reported) |> summarise(n = n(), n_fdr05 = sum(FDR < 0.05, na.rm = TRUE),
                                                          n_same = sum(same_direction, na.rm = TRUE), n_dir = sum(!is.na(same_direction))),
           pb_mast_by = vs5 |> filter(reported) |> group_by(contrast) |> summarise(n = n(), n_fdr05 = sum(FDR < 0.05, na.rm = TRUE), .groups = "drop"))
time_main_fig4 <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

t0 <- Sys.time()
FDR_CUT5 <- 0.05; SIZE_CAP5 <- 10
pp5 <- readr::read_csv(RES("06_ec_gsea", "ec_gsea_pathway_panel.csv"), show_col_types = FALSE)
theme_order5 <- unique(pp5$theme); set_order5 <- unique(pp5$set)
lab_contrast <- setNames(names(CONTRAST_LABEL), CONTRAST_LABEL)
pp5 <- pp5 |>
  mutate(tested = status == "tested", sig = tested & FDR < FDR_CUT5,
         contrast = factor(lab_contrast[contrast], levels = CO_ORDER), celltype = factor(subtype, levels = EC_SUBTYPES),
         theme = factor(theme, levels = theme_order5), set = factor(set, levels = rev(set_order5)),
         sac = tested & single_animal_candidate %in% TRUE)
stopifnot(!anyNA(pp5$contrast), !anyNA(pp5$celltype), nrow(pp5) == length(set_order5) * 25)
nes_lim5 <- ceiling(max(abs(pp5$NES), na.rm = TRUE) * 2) / 2
p5 <- ggplot(filter(pp5, tested), aes(contrast, set)) +
  geom_point(aes(fill = NES, size = pmin(-log10(FDR), SIZE_CAP5), colour = sig, stroke = ifelse(sig, .6, .25)), shape = 21) +
  geom_point(data = filter(pp5, sac), aes(shape = "single-animal candidate"), size = .5, colour = "black") +
  geom_point(data = filter(pp5, !tested), aes(shape = "not tested"), size = 3, colour = "grey40") +
  scale_shape_manual(values = c(`single-animal candidate` = 16, `not tested` = 45), name = NULL) +
  facet_grid(theme ~ celltype, scales = "free_y", space = "free_y", switch = "y",
             labeller = labeller(celltype = short_subtype, theme = wrap_lab(13))) +
  scale_x_discrete(labels = CONTRAST_LABEL) +
  scale_y_discrete(labels = wrap_lab(36), position = "right") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, limits = c(-nes_lim5, nes_lim5), breaks = function(l) pretty(l, n = 3), name = "NES") +
  scale_size_continuous(range = c(.5, 3.2), limits = c(0, SIZE_CAP5), breaks = c(1, 2, 5, 10), labels = c("1", "2", "5", "\u2265 10"),
                        name = "\u2212log10 FDR") +
  scale_colour_manual(values = c(`FALSE` = "grey65", `TRUE` = "black"),
                      labels = c(`FALSE` = sprintf("FDR \u2265 %s", FDR_CUT5), `TRUE` = sprintf("FDR < %s", FDR_CUT5)), name = NULL) +
  guides(fill = guide_colourbar(barwidth = unit(20, "mm"), barheight = unit(2, "mm"), order = 1, title.vjust = .9),
         size = guide_legend(order = 2, nrow = 1), colour = guide_legend(order = 3, nrow = 2, override.aes = list(size = 2.2, fill = "white", stroke = .6)),
         shape = guide_legend(order = 4, nrow = 2, override.aes = list(size = c(3, 1.2)))) +
  labs(x = NULL, y = NULL) +
  theme_fig(grid = "both") +
  theme(axis.line = element_blank(), axis.ticks = element_line(linewidth = .2, colour = "grey60"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5), axis.text.y.right = element_text(lineheight = .8, hjust = 0),
        strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1, size = BASE_PT, face = "bold", lineheight = .85),
        strip.text.x = element_text(size = TITLE_PT, face = "bold"), strip.clip = "off",
        panel.spacing.x = unit(1.5, "mm"), panel.spacing.y = unit(1, "mm"),
        legend.position = "bottom", legend.box.spacing = unit(1, "mm"), legend.spacing.x = unit(4, "mm"))
figlog_5 <- save_figure(wrap_plots(p5), "Figure_5", WIDTH_2COL, 225, tags = FALSE)
age5 <- pp5 |> filter(sig, contrast %in% c("OC_vs_YC", "OH_vs_YH"), NES < 0) |> group_by(set, contrast) |>
  summarise(n_lower = n_distinct(celltype), .groups = "drop") |> filter(n_lower >= 3)
L5 <- list(n_sets = length(set_order5), n_themes = length(theme_order5), n_tested = sum(pp5$tested), n_not = sum(!pp5$tested),
           not_tested = pp5 |> filter(!tested) |> count(set), sig_by = pp5 |> filter(sig) |> count(contrast, .drop = FALSE),
           n_sig = sum(pp5$sig), n_sig_sac = sum(pp5$sig & pp5$sac), n_sac = sum(pp5$sac),
           sac_top = pp5 |> filter(sig, sac) |> nrow(), age = age5,
           sets_by_db = pp5 |> distinct(set, pathway) |> mutate(db = sub("_.*$", "", pathway)) |> count(db))
time_main_fig5 <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

t0 <- Sys.time()
suppressPackageStartupMessages(library(CellChat))
cc_list <- readRDS(RES("07_ccc", "cellchat_list.rds"))
stopifnot(identical(names(cc_list), GRP_ORDER))
CC_TYPES <- levels(cc_list[[1]]@idents)
gi <- setNames(seq_along(GRP_ORDER), GRP_ORDER)
tot6 <- data.frame(group = factor(GRP_ORDER, levels = GRP_ORDER),
                   count = sapply(cc_list, function(x) sum(x@net$count)),
                   weight = sapply(cc_list, function(x) sum(x@net$weight))) |>
  pivot_longer(c(count, weight), names_to = "measure") |>
  mutate(measure = factor(measure, levels = c("count", "weight"), labels = c("inferred interactions", "interaction strength")),
         lab = ifelse(measure == "inferred interactions", fmt_n(round(value)), sprintf("%.1f", value)))
pA6 <- ggplot(tot6, aes(group, value, fill = group)) +
  geom_col(width = .7) +
  geom_text(aes(label = lab), vjust = -.35, size = pt_text(MIN_PT)) +
  facet_wrap(~ measure, scales = "free_y", nrow = 1, strip.position = "left") +
  scale_fill_manual(values = GRP_FILL, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, .14)), labels = scales::label_comma()) +
  labs(x = NULL, y = NULL) +
  theme_fig(grid = "y") +
  theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 90, face = "plain", size = BASE_PT), panel.spacing.x = unit(5, "mm"))

invisible(capture.output(suppressMessages({
  cc_l   <- lapply(cc_list, liftCellChat, group.new = CC_TYPES)
  merged <- suppressWarnings(mergeCellChat(cc_l, add.names = GRP_ORDER, cell.prefix = TRUE))
})))
rm(cc_l, cc_list); invisible(gc())
merged@net <- lapply(merged@net, function(n) { dimnames(n$weight) <- rep(list(unname(CT_SHORT[rownames(n$weight)])), 2); n })
draw_circle <- function(ctr) {
  par(xpd = NA, mar = c(0.2, 0.8, 1.1, 0.8), family = FONT)

  netVisual_diffInteraction(merged, weight.scale = TRUE, measure = "weight", comparison = c(gi[[ctr$i2]], gi[[ctr$i1]]),
                            color.use = unname(CT_FILL[CC_TYPES]), vertex.label.cex = BASE_PT / 12, margin = 0.22,
                            edge.width.max = 4, arrow.size = 0.1, vertex.size.max = 9, title.name = "")
  mtext(sprintf("%s: %s vs %s", CONTRAST_LABEL[[ctr$name]], ctr$i1, ctr$i2), side = 3, line = 0, cex = BASE_PT / 12)
}
circle_one <- function(ctr) wrap_elements(full = local({ ctr <- ctr; ~ draw_circle(ctr) }))
pB6 <- wrap_elements(full = wrap_plots(lapply(CONTRASTS, circle_one), ncol = 2)) + theme(plot.margin = margin(3, 0, 0, 0, "mm"))

dec6 <- readr::read_csv(RES("07_ccc", "tables", "ccc_diff_edge_decomposition.csv"), show_col_types = FALSE)
cmec6 <- dec6 |> filter(source == "Cardiomyocyte", target == "Endothelial") |>
  group_by(contrast) |> slice_max(abs(delta), n = 8, with_ties = FALSE) |> ungroup() |>
  mutate(contrast = factor(contrast, levels = sapply(CONTRASTS, `[[`, "name")), row = paste(contrast, interaction_name, sep = "||"))
edge_tot6 <- dec6 |> filter(source == "Cardiomyocyte", target == "Endothelial") |> distinct(contrast, edge_delta_total)
strip_C6 <- setNames(sprintf("%s\ntotal \u0394 = %+.3f", CONTRAST_LABEL[edge_tot6$contrast], edge_tot6$edge_delta_total), edge_tot6$contrast)
cmec6$row <- factor(cmec6$row, levels = cmec6$row[order(cmec6$contrast, cmec6$delta)])
pC6 <- ggplot(cmec6, aes(delta, row, fill = delta > 0)) +
  geom_col(width = .7) +
  geom_vline(xintercept = 0, colour = "grey50", linewidth = .3) +
  facet_wrap(~ contrast, nrow = 1, scales = "free", labeller = as_labeller(strip_C6)) +
  scale_fill_manual(values = c(`TRUE` = "#b2182b", `FALSE` = "#2166ac"), guide = "none") +
  scale_y_discrete(labels = function(x) gsub("_", " ", sub("^.*\\|\\|", "", x))) +
  scale_x_continuous(expand = expansion(mult = c(.12, .12)), breaks = function(lim) pretty(lim, n = 2)) +
  labs(x = "\u0394 communication probability, CM \u2192 EC (first \u2212 second group)", y = NULL) +
  theme_fig(grid = "x") +
  theme(axis.text.y = element_text(size = MIN_PT), panel.spacing.x = unit(2, "mm"), strip.clip = "off",
        strip.text = element_text(size = BASE_PT, face = "plain", lineheight = .9))

kd6 <- data.frame(ct = CC_TYPES, i = seq_along(CC_TYPES)) |>
  mutate(col = (i - 1) %/% 6, row = (i - 1) %% 6, lab = sprintf("%s = %s", CT_SHORT[ct], CT_LABEL[ct]))
pKey6 <- ggplot(kd6, aes(col * 1.5, -row)) +
  geom_point(aes(colour = ct), size = 2) +
  geom_text(aes(x = col * 1.5 + .1, label = lab), hjust = 0, size = pt_text()) +
  annotate("segment", x = c(0, 1.5), xend = c(.25, 1.75), y = -6.5, yend = -6.5, colour = c("#b2182b", "#2166ac"), linewidth = 1) +
  annotate("text", x = c(.32, 1.82), y = -6.5, hjust = 0, size = pt_text(), label = c("higher in first group", "higher in second group")) +
  scale_colour_manual(values = CT_FILL[CC_TYPES], guide = "none") +
  coord_cartesian(xlim = c(-.1, 3.1), ylim = c(-7, .5), expand = FALSE, clip = "off") +
  theme_void(base_size = BASE_PT)
left6 <- (pA6 / pKey6) + plot_layout(heights = c(56, 36))
top6  <- (wrap_elements(full = left6) | pB6) + plot_layout(widths = c(80, 100))
fig6  <- top6 / pC6 + plot_layout(heights = c(100, 55))
figlog_6 <- save_figure(fig6, "Figure_6", WIDTH_2COL, 165)
rm(fig6, top6, pB6, merged); invisible(gc())
L6 <- list(tot = tot6 |> select(group, measure, value), edge_tot = edge_tot6, cmec = cmec6 |> select(contrast, interaction_name, pathway, delta, share_of_edge),
           types = CC_TYPES)
time_main_fig6 <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

t0 <- Sys.time()
FDR_CUT7 <- 0.05; TOP_N7 <- 40; SEED7 <- 1
KEY7 <- c("source", "target", "ligand_complex", "receptor_complex")
diff7  <- readr::read_csv(RES("07_ccc", "tables", "lr_diff_all.csv"), show_col_types = FALSE) |> mutate(contrast = factor(contrast, levels = CO_ORDER))
score7 <- readr::read_csv(RES("07_ccc", "lr_scores_by_sample.csv"), show_col_types = FALSE)
diff7b <- readr::read_csv(RES("07b_ec_cm_subtype_ccc", "tables", "lr_subtype_diff_all.csv"), show_col_types = FALSE) |> mutate(contrast = factor(contrast, levels = CO_ORDER))
support7b <- read.csv(RES("07b_ec_cm_subtype_ccc", "tables", "subtype_animal_support.csv"))

ecm7 <- diff7 |> mutate(lr_pair = paste(ccc_pair, ligand_complex, receptor_complex))
keep7 <- ecm7 |> group_by(lr_pair) |> summarise(max_t = max(abs(t_mod)), .groups = "drop") |> slice_max(max_t, n = TOP_N7, with_ties = FALSE)
hm7 <- ecm7 |> filter(lr_pair %in% keep7$lr_pair) |>
  select(lr_pair, source, target, ligand_complex, receptor_complex, contrast, logFC, passes) |>
  complete(nesting(lr_pair, source, target, ligand_complex, receptor_complex), contrast) |>
  left_join(keep7, by = "lr_pair") |>
  mutate(label = sprintf("%s: %s \u2192 %s", CT_SHORT[source], gsub("_", " ", ligand_complex), gsub("_", " ", receptor_complex)),
         target_lab = paste("to", CT_SHORT[target]))
hm7$label <- factor(hm7$label, levels = unique(hm7$label[order(hm7$max_t)]))
stopifnot(!any(hm7$passes %in% TRUE))
lim7 <- max(abs(hm7$logFC), na.rm = TRUE)
pA7 <- ggplot(hm7, aes(contrast, label, fill = logFC)) +
  geom_tile(colour = "white", linewidth = .4) +
  facet_grid(target_lab ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, limits = c(-lim7, lim7),
                       name = "difference of\nmean log1p\nscore", na.value = "grey88") +
  scale_x_discrete(drop = FALSE, labels = CONTRAST_LABEL) +
  guides(fill = guide_colourbar(barwidth = unit(2, "mm"), barheight = unit(16, "mm"))) +
  labs(x = NULL, y = NULL) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.y = element_text(size = MIN_PT),
        axis.text.x = element_text(angle = 45, hjust = 1), strip.text.y = element_text(angle = -90, face = "bold", size = BASE_PT),
        panel.spacing.y = unit(1, "mm"), legend.position = "right", legend.justification = "top")

hist7 <- bind_rows(diff7 |> transmute(layer = "broad types", contrast, p_value, FDR),
                   diff7b |> transmute(layer = "subtypes", contrast, p_value, FDR)) |>
  mutate(layer = factor(layer, levels = c("broad types", "subtypes")))
summ7 <- hist7 |> group_by(layer, contrast) |>
  summarise(n = n(), n_p05 = sum(p_value < 0.05), expected = n() * 0.05, min_fdr = min(FDR), n_fdr = sum(FDR < FDR_CUT7), .groups = "drop") |>
  mutate(lab = sprintf("p < 0.05: %s / %s\nexpected %s\nmin FDR %.2f", fmt_n(n_p05), fmt_n(n), fmt_n(round(expected)), min_fdr))
pB7 <- ggplot(hist7, aes(p_value)) +
  geom_histogram(aes(y = after_stat(density)), breaks = seq(0, 1, .05), fill = "grey60", colour = "white", linewidth = .15) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = .3, colour = "grey20") +
  geom_text(data = summ7, aes(x = .5, y = Inf, label = lab), hjust = .5, vjust = 1.15, size = pt_text(MIN_PT), lineheight = .9) +
  facet_grid(contrast ~ layer, labeller = labeller(contrast = as_labeller(CONTRAST_CODE)), switch = "y") +
  scale_x_continuous(breaks = c(0, .5, 1), labels = c("0", "0.5", "1"), expand = expansion(mult = .01)) +
  scale_y_continuous(expand = expansion(mult = c(0, 1.25)), breaks = c(0, 1)) +
  labs(x = "nominal p-value", y = "density (uniform = 1)") +
  theme(strip.text.y.left = element_text(angle = 0, size = BASE_PT, face = "plain", lineheight = .9, hjust = 1), strip.placement = "outside",
        strip.text.x = element_text(size = BASE_PT), strip.clip = "off", panel.spacing.y = unit(1.2, "mm"), panel.spacing.x = unit(2, "mm"))

lead7 <- ecm7 |> filter(lr_pair %in% head(keep7$lr_pair, 6)) |> distinct(across(all_of(c("lr_pair", KEY7))))
lead7 <- lead7[match(head(keep7$lr_pair, 6), lead7$lr_pair), ]
panel_lab7 <- function(s, t, l, r) sprintf("%s \u2192 %s\n%s \u2192\n%s", CT_SHORT[s], CT_SHORT[t], gsub("_", " ", l), gsub("_", " ", r))
dots7 <- score7 |> inner_join(lead7, by = KEY7) |>
  mutate(panel = factor(panel_lab7(source, target, ligand_complex, receptor_complex),
                        levels = panel_lab7(lead7$source, lead7$target, lead7$ligand_complex, lead7$receptor_complex)),
         group = factor(disease, levels = GRP_ORDER), x_jit = as.numeric(group) + withr::with_seed(SEED7, runif(n(), -.12, .12)))
stopifnot(all(table(dots7$panel) == 16))
FLAG_SAMPLE7 <- names(FLAG_OF)[FLAG_OF == "low_yield_library"]
ytop7 <- dots7 |> group_by(panel) |> summarise(ytop = max(expr_prod), .groups = "drop")
lab7 <- dots7 |> filter(sample %in% FLAG_SAMPLE7) |> left_join(ytop7, by = "panel") |>
  mutate(y_lab = ytop * 1.2, x_lab = pmin(x_jit, length(GRP_ORDER) - .15))
pC7 <- ggplot(dots7, aes(group, expr_prod)) +
  stat_summary(aes(fill = group), fun = mean, geom = "col", width = .7) +
  geom_point(aes(x = x_jit, shape = sex), size = 1.1, stroke = .3, fill = "grey85", colour = "grey10") +
  geom_segment(data = lab7, aes(x = x_lab, xend = x_jit, y = y_lab * .96, yend = expr_prod * 1.03), linewidth = .2, colour = "grey45") +
  geom_text(data = lab7, aes(x = x_lab, y = y_lab, label = sample), size = pt_text(MIN_PT), vjust = 0, hjust = .8) +
  facet_wrap(~ panel, scales = "free_y", nrow = 1) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
  scale_fill_manual(values = GRP_FILL, labels = GRP_LABEL, name = NULL) +
  scale_shape_manual(values = c(F = 21, M = 24), labels = c(F = "female", M = "male"), name = NULL) +
  guides(fill = guide_legend(nrow = 1, order = 1), shape = guide_legend(order = 2, override.aes = list(size = 1.6))) +
  labs(x = NULL, y = "expression-product score") +
  theme_fig(grid = "y") +
  theme(strip.text = element_text(size = BASE_PT, face = "plain", lineheight = .9), strip.clip = "off", legend.position = "bottom",
        legend.margin = margin(-2, 0, 0, 0, "mm"), panel.spacing.x = unit(2, "mm"))

edge7b <- diff7b |> group_by(direction, source, target, contrast) |>
  summarise(n = n(), n_p05 = sum(p_value < 0.05), min_fdr = min(FDR), .groups = "drop") |>
  mutate(pct_p05 = 100 * n_p05 / n,
         ec = ifelse(direction == unique(direction)[grepl("^EC", unique(direction))], source, target),
         cm = ifelse(ec == source, target, source),
         edge = sprintf("%s \u2192 %s", ifelse(ec == source, short_subtype(ec), short_cm(cm)), ifelse(ec == source, short_cm(cm), short_subtype(ec))))
edge_order <- edge7b |> distinct(direction, ec, cm, edge) |>
  mutate(ec = factor(ec, levels = EC_SUBTYPES), cm = factor(cm, levels = cm2$subtypes)) |> arrange(direction, cm, ec)
edge7b <- edge7b |> mutate(edge = factor(edge, levels = unique(edge_order$edge)), contrast = factor(contrast, levels = rev(CO_ORDER)))
pD7 <- ggplot(edge7b, aes(edge, contrast, fill = pct_p05)) +
  geom_tile(colour = "white", linewidth = .4) +
  geom_text(aes(label = sprintf("%.0f", pct_p05)), size = pt_text(MIN_PT)) +
  facet_grid(. ~ direction, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 5, limits = c(0, max(10, max(edge7b$pct_p05))),
                       name = "tests with\np < 0.05 (%)") +
  scale_y_discrete(labels = CONTRAST_LABEL) +
  guides(fill = guide_colourbar(barwidth = unit(2, "mm"), barheight = unit(12, "mm"))) +
  labs(x = NULL, y = NULL) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_text(angle = 40, hjust = 1, size = MIN_PT),
        strip.text = element_text(size = BASE_PT, face = "bold"), panel.spacing.x = unit(2, "mm"))

top7 <- (pA7 | wrap_elements(full = pB7)) + plot_layout(widths = c(100, 80))
fig7 <- top7 / wrap_elements(full = pC7) / pD7 + plot_layout(heights = c(112, 52, 46))
figlog_7 <- save_figure(fig7, "Figure_7", WIDTH_2COL, 225)
rm(fig7); invisible(gc())
L7 <- list(summ = summ7 |> select(-lab), n_tests = nrow(diff7), n_pairs = n_distinct(ecm7$lr_pair), min_fdr = min(diff7$FDR),
           min_fdr_all = min(diff7$FDR_all_contrasts), n_p05 = sum(diff7$p_value < 0.05), n_fdr = sum(diff7$FDR < FDR_CUT7),
           n_tests_b = nrow(diff7b), min_fdr_b = min(diff7b$FDR), n_p05_b = sum(diff7b$p_value < 0.05), n_fdr_b = sum(diff7b$FDR < FDR_CUT7),
           n_edges_b = n_distinct(paste(diff7b$source, diff7b$target)), lead = lead7,
           lead_best = ecm7 |> filter(lr_pair %in% lead7$lr_pair) |> group_by(lr_pair) |> slice_max(abs(t_mod), n = 1, with_ties = FALSE) |>
             ungroup() |> select(lr_pair, contrast, logFC, t_mod, p_value, FDR),
           flag_rank = dots7 |> group_by(panel, group) |> mutate(rank = rank(-expr_prod)) |> ungroup() |> filter(sample %in% FLAG_SAMPLE7) |> select(panel, group, expr_prod, rank),
           flag_sample = FLAG_SAMPLE7, n_animals = n_distinct(score7$sample), support = support7b,
           edge = edge7b, pct_range = range(edge7b$pct_p05))
time_main_fig7 <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

# [Removed from the public code export: main figure legend text (code/legends_main.R). Available from the authors.]
