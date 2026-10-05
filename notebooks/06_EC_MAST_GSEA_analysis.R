# 06_EC_MAST_GSEA_analysis.R
# Code-only export of the analysis pipeline. Data are not included; see README.md.

params <- list(universe_csv = "results/05_ec_de/ec_DE_all_MAST_gsea_universe.csv", animal_csv = "results/05_ec_de/ec_DE_per_animal_expression.csv", 
    out_dir = "results/06_ec_gsea", nperm_simple = 100000L)


knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = TRUE,
                      fig.align = "center", dpi = 150)
options(width = 110)

suppressPackageStartupMessages({
  library(data.table); library(fgsea); library(scales); library(svglite); library(DT); library(stringr)
})
source("_common.R")
knitr::opts_chunk$set(cache = TRUE, autodep = TRUE, cache.lazy = FALSE, cache.extra = c(input_md5(params), params$nperm_simple))
path_universe_csv    <- resolve_path(params$universe_csv)
path_animal_csv <- resolve_path(params$animal_csv)
path_out    <- resolve_path(params$out_dir)
path_figures <- file.path(path_out, "figures"); dir.create(path_figures, recursive = TRUE, showWarnings = FALSE)
path_gsea_rds   <- file.path(path_out, "ec_gsea_mast.rds")

saved_to <- function(path) { cat(sprintf("Saved: %s\n", sub(paste0(PROJECT_ROOT, "/"), "", normalizePath(path, winslash = "/"), fixed = TRUE))); invisible(path) }
dt_table <- function(df, caption) DT::datatable(df, rownames = FALSE, filter = "top", caption = caption, class = "compact stripe hover",
  options = list(pageLength = 15, lengthMenu = c(10, 15, 25, 50, 100), scrollX = TRUE, searchHighlight = TRUE, autoWidth = TRUE))

FDR_CUT    <- 0.05
MIN_SIZE   <- 10; MAX_SIZE <- 500
NES_LIMITS <- c(-2, 2)
SIZE_CAP   <- 10
TECH_FRAC  <- 0.5
MAX_ANIMAL_SHARE <- 0.5
JACCARD_COLLAPSE <- 0.5
DASH_HALF  <- 0.15
TERMS_PER_PAGE <- 25
TABLE_ROWS <- 200
PATHWAY_PANEL <- list(
  "Angiogenesis / VEGF"      = c("HALLMARK_ANGIOGENESIS", "REACTOME_SIGNALING_BY_VEGF", "GOBP_SPROUTING_ANGIOGENESIS", "GOBP_ENDOTHELIAL_CELL_PROLIFERATION"),
  "EndMT / TGF-beta"         = c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_TGF_BETA_SIGNALING"),
  "Inflammation"             = c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_INFLAMMATORY_RESPONSE",
                                 "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE"),
  "Leukocyte adhesion"       = c("REACTOME_CELL_SURFACE_INTERACTIONS_AT_THE_VASCULAR_WALL", "GOBP_LEUKOCYTE_ADHESION_TO_VASCULAR_ENDOTHELIAL_CELL",
                                 "GOBP_LEUKOCYTE_MIGRATION"),
  "Senescence"               = c("HALLMARK_P53_PATHWAY", "REACTOME_CELLULAR_SENESCENCE", "REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP"),
  "Oxidative stress"         = c("HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"),
  "Nitric oxide"             = c("REACTOME_METABOLISM_OF_NITRIC_OXIDE_NOS3_ACTIVATION_AND_REGULATION", "GOBP_NITRIC_OXIDE_MEDIATED_SIGNAL_TRANSDUCTION"),
  "Metabolism"               = c("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_ADIPOGENESIS",
                                 "HALLMARK_GLYCOLYSIS", "HALLMARK_HYPOXIA"),
  "Growth signalling"        = c("HALLMARK_MTORC1_SIGNALING", "HALLMARK_PI3K_AKT_MTOR_SIGNALING"),
  "Coagulation / complement" = c("HALLMARK_COAGULATION", "HALLMARK_COMPLEMENT"),
  "ECM / basement membrane"  = c("REACTOME_COLLAGEN_FORMATION", "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION", "REACTOME_LAMININ_INTERACTIONS"),
  "Apoptosis"                = c("HALLMARK_APOPTOSIS"))
SEED <- 1; set.seed(SEED)

theme_gsea <- function() theme_house(base_size = 10) +
  theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1, vjust = 1), axis.text.y = element_text(size = 7),
        axis.title = element_blank(), plot.caption = element_text(size = 7, colour = "grey40", hjust = 0),
        legend.title = element_text(size = 7), legend.text = element_text(size = 6), legend.key = element_blank())

print_versions(c("ggplot2", "dplyr", "data.table", "fgsea", "scales", "svglite", "DT", "stringr"))

mast_all <- fread(path_universe_csv)
mast_tab <- mast_all[!is.na(gene) & nzchar(gene) & is.finite(avg_log2FC) & is.finite(p_val) & p_val >= 0 & p_val <= 1]
not_estimable <- mast_all[contrast == "Interaction", .(N = sum(!is.finite(avg_log2FC))), by = celltype]
SUBTYPE_ORDER <- unique(mast_tab$celltype)
CONTRASTS_RUN        <- intersect(CONTRAST_ORDER, unique(mast_tab$contrast))
CONTRAST_PANELS        <- unname(CONTRAST_LABEL[CONTRASTS_RUN])
mast_tab[, saturated := p_val <= .Machine$double.xmin]
mast_tab[, rankstat := sign(avg_log2FC) * (-log10(pmax(p_val, .Machine$double.xmin)) + saturated * abs(avg_log2FC))]
mast_tab[, technical := is_technical_gene(gene)]
universe <- mast_tab[, .(universe_n = .N, technical_genes = sum(technical),
                     saturated = sum(saturated), tied_scores = .N - uniqueN(rankstat)), by = .(celltype, contrast)]

animal_expr     <- fread(path_animal_csv)
group_of_animal     <- unique(animal_expr[, .(individual, disease)])[, setNames(disease, individual)]
animal_expr_mat <- lapply(split(animal_expr, by = "celltype"), function(celltype_dt) {
  wide <- dcast(celltype_dt, gene ~ individual, value.var = "mean_expr"); expr_mat <- as.matrix(wide[, -1]); rownames(expr_mat) <- wide$gene; expr_mat })
UNIVERSE_TXT <- universe[, .(txt = sprintf("%s %s", short_subtype(celltype[1]), if (uniqueN(universe_n) == 1) format(universe_n[1], big.mark = ",") else
                                            sprintf("%s\u2013%s", format(min(universe_n), big.mark = ","), format(max(universe_n), big.mark = ",")))),
                         by = celltype][, paste(txt, collapse = "; ")]
knitr::kable(universe[order(match(celltype, SUBTYPE_ORDER), match(contrast, CONTRASTS_RUN))],
             caption = "Ranking universe per block (genes detected in >= 1% of the nuclei of either compared group; Interaction: of the subtype's pooled nuclei; MAST statistics from notebook 05), technical genes among them, saturated p-values and remaining tied scores.")

path_gene_sets <- file.path(PATHS$resources, "ec_gsea_gene_sets.rds")
gs_snapshot <- readRDS(path_gene_sets)
GS_LIST <- gs_snapshot$sets
cat(sprintf("MSigDB %s (msigdbr %s, retrieved %s): Hallmark %d | Reactome %d | GO:BP %d sets\n", gs_snapshot$db_version,
            gs_snapshot$msigdbr, gs_snapshot$retrieved_on, length(GS_LIST$Hallmark), length(GS_LIST$Reactome), length(GS_LIST[["GO:BP"]])))
n_sets_mt <- sum(vapply(unlist(GS_LIST, recursive = FALSE), function(set_genes) any(is_mt_gene(set_genes)), logical(1)))

make_rank <- function(rank_dt) {
  stopifnot(!anyDuplicated(rank_dt$gene))
  rank_dt <- rank_dt[order(-rankstat, gene)]
  setNames(rank_dt$rankstat, rank_dt$gene)
}
run_fgsea <- function(rank_stats, gene_sets) as.data.table(fgseaMultilevel(pathways = gene_sets, stats = rank_stats, minSize = MIN_SIZE, maxSize = MAX_SIZE,
                                                             eps = 0, nproc = 1, nPermSimple = params$nperm_simple))
pair_of <- setNames(lapply(CONTRASTS, \(cc) c(cc$i1, cc$i2)), vapply(CONTRASTS, `[[`, "", "name"))
set_share <- function(leading_edge, set_subtype, contrast_name, nes_value) {
  expr_mat <- animal_expr_mat[[set_subtype]]; le_expr <- expr_mat[intersect(leading_edge, rownames(expr_mat)), , drop = FALSE]; le_expr <- le_expr[rowSums(le_expr) > 0, , drop = FALSE]
  animal_group <- group_of_animal[colnames(le_expr)]
  if (contrast_name == "Interaction") {
    group_mean_score <- tapply(colMeans(le_expr / rowMeans(le_expr)), animal_group, mean)
    hfpef_pair <- if (abs(group_mean_score[["OH"]] - group_mean_score[["OC"]]) >= abs(group_mean_score[["YH"]] - group_mean_score[["YC"]])) c("OH", "OC") else c("YH", "YC")
    high_low <- if (group_mean_score[[hfpef_pair[1]]] >= group_mean_score[[hfpef_pair[2]]]) hfpef_pair else rev(hfpef_pair)
  } else high_low <- if (nes_value > 0) pair_of[[contrast_name]] else rev(pair_of[[contrast_name]])
  le_expr  <- le_expr[, animal_group %in% high_low, drop = FALSE]; le_expr <- le_expr[rowSums(le_expr) > 0, , drop = FALSE]; animal_group <- group_of_animal[colnames(le_expr)]
  set_score <- colMeans(le_expr / rowMeans(le_expr))
  excess_signal <- pmax(set_score[animal_group == high_low[1]] - median(set_score[animal_group == high_low[2]]), 0)
  list(share = if (sum(excess_signal) > 0) max(excess_signal) / sum(excess_signal) else NA_real_, animal = if (sum(excess_signal) > 0) names(which.max(excess_signal)) else "",
       expr_share = max(set_score[animal_group == high_low[1]]) / sum(set_score[animal_group == high_low[1]]))
}
skip_file <- file.path(path_out, "ec_gsea_skipped.csv")
set.seed(SEED); t_start <- Sys.time()
all_results <- list(); skipped <- list()
for (set_subtype in SUBTYPE_ORDER) for (contrast_name in CONTRASTS_RUN) {
  rank_dt <- mast_tab[celltype == set_subtype & contrast == contrast_name]
  rank_stats <- make_rank(rank_dt)
  if (length(rank_stats) < 50 || !any(rank_stats != 0)) {
    skipped[[paste(set_subtype, contrast_name)]] <- data.frame(celltype = set_subtype, contrast = contrast_name, reason = "fewer than 50 genes or no nonzero ranking scores"); next
  }
  for (gsn in names(GS_LIST)) {
    fgsea_res <- run_fgsea(rank_stats, GS_LIST[[gsn]])
    fgsea_res[, `:=`(celltype = set_subtype, contrast = contrast_name, panel = CONTRAST_LABEL[[contrast_name]], geneset = gsn, universe_n = length(rank_stats))]
    all_results[[paste(set_subtype, contrast_name, gsn)]] <- fgsea_res
  }
}
gsea_minutes <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
skip_table <- rbindlist(c(list(data.table(celltype = character(), contrast = character(), reason = character())), skipped))
fwrite(skip_table, skip_file)
gsea <- rbindlist(all_results)
gsea[, padj_all := p.adjust(pval, "BH")]
gsea[, tech_frac := vapply(leadingEdge, function(leading_edge) mean(is_technical_gene(leading_edge)), numeric(1))]
gsea[, technical_candidate := tech_frac >= TECH_FRAC]
gsea[, `:=`(sad_share = NA_real_, sad_animal = "")]
gsea[padj < FDR_CUT, c("sad_share", "sad_animal") :=
       rbindlist(Map(\(leading_edge, set_subtype, contrast_name, nes_value) set_share(leading_edge, set_subtype, contrast_name, nes_value)[c("share", "animal")], leadingEdge, celltype, contrast, NES))]
gsea[, single_animal_candidate := (sad_share > MAX_ANIMAL_SHARE) %in% TRUE]
saveRDS(gsea, path_gsea_rds)
gsea[, `:=`(sig = is.finite(padj) & padj < FDR_CUT, lfdr = -log10(pmax(padj, 1e-300)))]

if (nrow(skip_table)) knitr::kable(skip_table, caption = "Comparisons not run")
knitr::kable(gsea[, .(tested = .N, significant = sum(sig), technical_candidate = sum(sig & technical_candidate),
                      single_animal_candidate = sum(sig & single_animal_candidate)), by = .(geneset)],
             caption = sprintf("Tested sets and hits (FDR < %s) per database; hits whose leading edge is >= %d%% technical genes, or whose high-side signal comes > %d%% from one animal.", FDR_CUT, TECH_FRAC * 100, MAX_ANIMAL_SHARE * 100))
sad_sets <- gsea[sig & single_animal_candidate]
single_animal_tab  <- dcast(sad_sets[, .N, by = .(sad_animal, geneset)], sad_animal ~ geneset, value.var = "N", fill = 0)
knitr::kable(single_animal_tab[order(-rowSums(single_animal_tab[, -1]))], caption = "Significant single-animal candidate sets by animal and database.")
hallmark_lab <- function(x) gsub("_", " ", sub("^HALLMARK_", "", x))
sig_results  <- gsea[sig == TRUE]
n_families    <- uniqueN(gsea[, .(geneset, celltype, contrast)])
example_set   <- sig_results[geneset == "Hallmark" & single_animal_candidate][order(-sad_share)][1]
example_share <- set_share(example_set$leadingEdge[[1]], example_set$celltype, example_set$contrast, example_set$NES)
hall_sad <- sad_sets[geneset == "Hallmark", .(n = .N, animals = paste(names(sort(table(sad_animal), decreasing = TRUE)), collapse = "/")), by = pathway][order(-n, pathway)]
note_animal   <- names(ANIMAL_NOTES)[1]

PARALLEL <- list(HFpEF = c("OH_vs_OC", "YH_vs_YC"), Age = c("OC_vs_YC", "OH_vs_YH"))
reversals <- rbindlist(lapply(names(PARALLEL), function(effect) {
  pair <- PARALLEL[[effect]]
  merge(sig_results[contrast == pair[1], .(geneset, celltype, pathway, NES_1 = NES, sad_1 = sad_animal, sad_flag_1 = single_animal_candidate)],
        sig_results[contrast == pair[2], .(geneset, celltype, pathway, NES_2 = NES, sad_2 = sad_animal, sad_flag_2 = single_animal_candidate)],
        by = c("geneset", "celltype", "pathway"))[sign(NES_1) != sign(NES_2)][, effect := effect]
}))
reversals <- gsea[contrast == "Interaction", .(geneset, celltype, pathway, interaction_sig = sig)][reversals, on = c("geneset", "celltype", "pathway")]
fwrite(reversals[, .(effect, geneset, celltype, pathway, NES_1, sad_1, sad_flag_1, NES_2, sad_2, sad_flag_2, interaction_sig)],
       file.path(path_out, "ec_gsea_sign_reversals.csv"))
rev_hf  <- reversals[effect == "HFpEF"]; rev_age <- reversals[effect == "Age"]
top_old <- sort(table(rev_hf$sad_1), decreasing = TRUE); top_young <- sort(table(rev_hf$sad_2), decreasing = TRUE)
count_text <- function(x) paste(sprintf("%d %s", as.vector(table(short_subtype(x$celltype))), names(table(short_subtype(x$celltype)))), collapse = ", ")

ACRONYMS <- c("ATP", "RNA", "DNA", "CDNA", "NMD", "SRP", "MHC", "I", "II", "III", "IV", "IGF", "TNF", "TNFA", "IFN", "IL", "VEGF",
              "PDGF", "TGF", "TGFB", "EGF", "EGFR", "FGF", "FGFR", "GPCR", "NADH", "NAD", "FAD", "GTP", "GDP", "AMP", "CAMP",
              "CGMP", "ROS", "ER", "UV", "HIV", "SARS", "EIF2AK4", "GCN2", "EIF2", "EIF4", "ECM", "LDL", "HDL", "VLDL", "PI3K", "AKT",
              "MTOR", "MTORC1", "MAPK", "ERK", "JAK", "STAT", "STAT3", "STAT5", "E2F", "G2M", "KRAS", "MYC", "P53", "TP53", "WNT", "TCA",
              "UPR", "SNARE", "RHO", "RAC1", "CDC42", "ABC", "SLC", "HSP", "HSF1", "NOTCH", "TLR", "NOD", "RIG", "CD4", "CD8", "TH1", "TH2",
              "TH17", "NK", "ADP", "PPAR", "PPARA", "ATF4", "ATF6", "XBP1", "HDAC", "FC", "IGE", "IGG", "MAP", "SUMO", "NLRP3", "ROBO",
              "SLIT", "SEMA", "PLC", "PKA", "PKC", "AURKA", "PLK1", "CDK", "S", "G1", "G2", "M", "NCAM1", "MET", "L1CAM", "RAF", "RAS",
              "DAP12", "TNFR2", "TNFR1", "NTRK", "NTRKS", "TRKA", "TRKB", "TRKC", "CD28", "CTLA4", "PD1", "TCR", "BCR", "AU", "GLI",
              "SHH", "IRE1", "PERK", "MDM2", "NOX", "NOS", "LPS", "PAMP", "DAMP", "OXPHOS", "EIF2AK1", "HRI", "PTEN", "MAP2K", "MLL3", "MLL4", "MITF", "CDH1")
SPECIAL  <- c(MRNA = "mRNA", MIRNA = "miRNA", RRNA = "rRNA", TRNA = "tRNA", NCRNA = "ncRNA", SNRNA = "snRNA", SNORNA = "snoRNA",
              LNCRNA = "lncRNA", COV = "CoV", NFKB = "NF-kB", IGFBPS = "IGFBPs", NTRKS = "NTRKs", CDKS = "CDKs", GTPASE = "GTPase", GTPASES = "GTPases")
pretty_lab <- function(x, geneset_db, width = 45) {
  x <- gsub("_", " ", sub("^(HALLMARK|REACTOME|GOBP)_", "", x))
  if (geneset_db != "Hallmark") {
    words <- strsplit(sub("^(\\w)", "\\U\\1", tools::toTitleCase(tolower(x)), perl = TRUE), " ")
    x <- vapply(words, function(word) { upper_word <- toupper(word)
      paste(ifelse(upper_word %in% names(SPECIAL), SPECIAL[upper_word], ifelse(upper_word %in% ACRONYMS, upper_word, word)), collapse = " ") }, character(1))
    x <- gsub("\\bSARS CoV\\b( (\\d))?", "SARS-CoV\\2", gsub("\\bNf Kb\\b", "NF-kB", x))
    x <- gsub("SARS-CoV(\\d)", "SARS-CoV-\\1", gsub("\\bMITF M\\b", "MITF-M", x))
  }
  str_wrap(x, width)
}

collapse_paths <- function(gsea_dt) {
  best <- gsea_dt[order(padj, -abs(NES)), .SD[1], by = pathway][order(padj, -abs(NES))]
  kept_edges <- list(); representative_of <- character()
  for (i in seq_len(nrow(best))) {
    leading_edge <- best$leadingEdge[[i]]
    jaccard <- vapply(kept_edges, function(k) length(intersect(leading_edge, k)) / length(union(leading_edge, k)), numeric(1))
    if (!length(kept_edges) || max(jaccard) <= JACCARD_COLLAPSE) { kept_edges[[best$pathway[i]]] <- leading_edge; representative_of[best$pathway[i]] <- best$pathway[i] }
    else representative_of[best$pathway[i]] <- names(kept_edges)[which.max(jaccard)]
  }
  data.table(pathway = names(representative_of), representative = unname(representative_of))
}

pick_paths <- function(gsea, geneset_db, max_paths) {
  gsea_dt <- gsea[geneset == geneset_db & sig]
  if (!nrow(gsea_dt)) return(list(paths = character(), collapsed = data.table()))
  collapsed <- collapse_paths(gsea_dt)
  ranked_paths <- gsea_dt[pathway %in% collapsed[pathway == representative]$pathway, .(n = .N, best = max(lfdr)), by = pathway][order(-n, -best)]
  list(paths = head(ranked_paths$pathway, max_paths), collapsed = collapsed[pathway != representative])
}

make_plot_df <- function(gsea, geneset_db, paths) {
  gsea_dt <- gsea[geneset == geneset_db & pathway %in% paths]
  path_labels <- make.unique(pretty_lab(paths, geneset_db))
  gsea_dt[, pathway_lab := factor(path_labels[match(pathway, paths)], levels = rev(path_labels))]
  gsea_dt[, panel := factor(panel, levels = CONTRAST_PANELS)]
  gsea_dt[, celltype := factor(celltype, levels = SUBTYPE_ORDER)]
  gsea_dt[, fill_nes := ifelse(technical_candidate, NA_real_, NES)]
  gsea_dt[]
}

single_animal_dot <- function(dt) geom_point(data = dt[single_animal_candidate == TRUE], shape = 16, size = 0.9, colour = "black", alpha = 1, show.legend = FALSE)
not_tested_dash   <- function(dt, x, y) geom_segment(data = dt, aes(x = as.integer(.data[[x]]) - DASH_HALF, xend = as.integer(.data[[x]]) + DASH_HALF,
                                                                    y = .data[[y]], yend = .data[[y]]), inherit.aes = FALSE, colour = "grey55", linewidth = 0.5)

gsea_dotplot <- function(plot_dt, title, subtitle, caption) {
  not_tested <- CJ(celltype = levels(plot_dt$celltype), panel = levels(plot_dt$panel), pathway_lab = levels(plot_dt$pathway_lab))[
    !plot_dt[, .(celltype = as.character(celltype), panel = as.character(panel), pathway_lab = as.character(pathway_lab))], on = c("celltype", "panel", "pathway_lab")]
  not_tested[, `:=`(celltype = factor(celltype, levels = levels(plot_dt$celltype)), panel = factor(panel, levels = levels(plot_dt$panel)),
                    pathway_lab = factor(pathway_lab, levels = levels(plot_dt$pathway_lab)))]
  ggplot(plot_dt, aes(celltype, pathway_lab, fill = fill_nes, size = pmin(lfdr, SIZE_CAP), alpha = sig, colour = sig)) +
    geom_point(shape = 21, stroke = 0.4) +
    single_animal_dot(plot_dt) +
    not_tested_dash(not_tested, "celltype", "pathway_lab") +
    scale_fill_gradientn(colours = c("#2166AC", "#6587C0", "#F7F7F7", "#DB7E78", "#B2182B"), limits = NES_LIMITS,
                         oob = squish, na.value = TECH_GREY, name = "NES") +
    scale_size_continuous(range = c(0.8, 5.2), breaks = c(1, 2, 5, 10), limits = c(0, SIZE_CAP),
                          labels = c("1", "2", "5", "\u2265 10"), name = expression(-log[10]~FDR)) +
    scale_alpha_manual(values = c(`FALSE` = 0.45, `TRUE` = 1), guide = "none") +
    scale_colour_manual(values = c(`FALSE` = "#BBBBBB", `TRUE` = "#333333"), guide = "none") +
    facet_wrap(~ panel, nrow = 1, drop = FALSE) +
    scale_x_discrete(drop = FALSE, labels = short_subtype) +
    labs(title = title, subtitle = subtitle, caption = caption) +
    guides(fill = guide_colourbar(order = 1, barwidth = 0.6, barheight = 4), size = guide_legend(order = 2)) +
    theme_gsea()
}

results_wide <- function(gsea, geneset_db) {
  gsea_dt <- gsea[geneset == geneset_db]
  gsea_dt <- gsea_dt[pathway %in% gsea_dt[sig == TRUE, unique(pathway)]]
  gsea_dt[, cell := sprintf("%+.2f (FDR %.2g)%s%s%s", NES, padj, ifelse(sig, "*", ""), ifelse(technical_candidate, " T", ""), ifelse(single_animal_candidate, " S", ""))]
  gsea_dt[, panel := factor(panel, levels = CONTRAST_PANELS)]
  words <- dcast(gsea_dt, celltype + pathway ~ panel, value.var = "cell", fill = "")
  words[, pathway := pretty_lab(pathway, geneset_db, width = 200)]
  words[order(celltype, pathway)]
}

FIG_SPEC <- list(
  Hallmark = list(file = "ec_gsea_hallmark_bubble.svg", title = "Hallmark enrichment in endothelial subtypes across age and HFpEF contrasts", w = 12.5, h = 5.4, max_n = 60),
  Reactome = list(file = "ec_gsea_reactome_bubble.svg", title = "Reactome enrichment in endothelial subtypes", w = 13.5, h = 10.5, max_n = 60),
  `GO:BP`  = list(file = "ec_gsea_gobp_bubble.svg",     title = "GO Biological Process enrichment in endothelial subtypes", w = 14.0, h = 9.5, max_n = 40))
FIG_CAPTION <- paste0(
  "Fill = NES (red: enriched toward genes higher in the first group of the contrast); grey = technical candidate (leading edge >= ", TECH_FRAC * 100, "% mt-/ribosomal genes).\n",
  "Size = -log10 FDR within database x subtype x contrast, capped at ", SIZE_CAP, " (", n_families, " families); dark ring = FDR < ", FDR_CUT, ", faded = not significant;\nblack centre dot = single-animal candidate (one animal > ", MAX_ANIMAL_SHARE * 100, "% of the high-side set signal); grey dash = not tested (< ", MIN_SIZE, " genes in the ranking universe).\n",
  "Redundant sets collapsed at leading-edge Jaccard > ", JACCARD_COLLAPSE, ". Per-nucleus MAST rankings (notebook 05); animal-level replication is assessed in notebook 05b.\n",
  "Ranking universe (genes detected in >= 1% of the nuclei of either compared group; Age x HFpEF: of the subtype's pooled nuclei) per subtype, across contrasts: ", UNIVERSE_TXT, ".")
COLLAPSED <- list()
make_fig <- function(geneset_db) {
  fig_spec <- FIG_SPEC[[geneset_db]]
  picked <- pick_paths(gsea, geneset_db, fig_spec$max_n)
  COLLAPSED[[geneset_db]] <<- picked$collapsed
  if (!length(picked$paths)) return(ggplot() + annotate("text", x = 0, y = 0, label = "No pathways pass the FDR threshold") + theme_void() + labs(title = fig_spec$title))
  fig <- gsea_dotplot(make_plot_df(gsea, geneset_db, picked$paths), title = fig_spec$title,
                    subtitle = sprintf("%s GSEA on covariate-adjusted MAST -log10 p signed by the effect direction (notebook 05); %d of %d significant sets shown, %d collapsed",
                                       geneset_db, length(picked$paths), uniqueN(gsea[geneset == geneset_db & sig]$pathway), nrow(picked$collapsed)),
                    caption = FIG_CAPTION)
  saved_to(save_svg(fig, fig_spec$file, width = fig_spec$w, height = fig_spec$h, dir = path_figures))
  fig
}

dt_table(results_wide(gsea, "Hallmark"), "Hallmark GSEA: pathways significant in >= 1 subtype x contrast. Cell = NES (FDR); * FDR < 0.05; T technical candidate; S single-animal candidate.")

make_fig("Hallmark")

dt_table(head(results_wide(gsea, "Reactome"), TABLE_ROWS), sprintf("Reactome GSEA: first %d rows of pathways significant in >= 1 subtype x contrast (full export in CSV).", TABLE_ROWS))

make_fig("Reactome")

dt_table(head(results_wide(gsea, "GO:BP"), TABLE_ROWS), sprintf("GO:BP GSEA: first %d rows of pathways significant in >= 1 subtype x contrast (full export in CSV).", TABLE_ROWS))

make_fig("GO:BP")

panel_lab <- function(x) {
  words <- strsplit(gsub("_", " ", sub("^(HALLMARK|REACTOME|GOBP)_", "", x)), " ")
  vapply(words, function(word) {
    word <- ifelse(word %in% names(SPECIAL), SPECIAL[word], ifelse(word %in% c(ACRONYMS, "SASP", "NOS3") | grepl("^[A-Z]+[0-9]+[A-Z]*$", word), word, tolower(word)))
    sub("^(\\w)", "\\U\\1", paste(word, collapse = " "), perl = TRUE) }, character(1))
}
panel_def <- data.table(theme = rep(names(PATHWAY_PANEL), lengths(PATHWAY_PANEL)), pathway = unlist(PATHWAY_PANEL, use.names = FALSE))
panel_def[, set := panel_lab(pathway)]
panel_grid <- CJ(celltype = SUBTYPE_ORDER, contrast = CONTRASTS_RUN, pathway = panel_def$pathway, sorted = FALSE)[panel_def, on = "pathway"]
panel_grid <- gsea[, .(celltype, contrast, pathway, NES, padj, size, single_animal_candidate, leadingEdge)][panel_grid, on = c("celltype", "contrast", "pathway")]
panel_grid[, `:=`(tested = !is.na(NES), sig = padj < FDR_CUT & !is.na(padj),
          celltype = factor(celltype, levels = SUBTYPE_ORDER), contrast = factor(contrast, levels = CONTRASTS_RUN),
          theme = factor(theme, levels = names(PATHWAY_PANEL)), set = factor(set, levels = rev(panel_def$set)))]
nes_lim <- ceiling(max(abs(panel_grid$NES), na.rm = TRUE) * 2) / 2
fig <- ggplot(panel_grid[tested == TRUE], aes(contrast, set)) +
  geom_point(aes(fill = NES, size = pmin(-log10(padj), SIZE_CAP), colour = sig, stroke = ifelse(sig, .8, .3)), shape = 21) +
  single_animal_dot(panel_grid[tested == TRUE]) +
  not_tested_dash(panel_grid[tested == FALSE], "contrast", "set") +
  facet_grid(theme ~ celltype, scales = "free_y", space = "free_y", labeller = labeller(celltype = short_subtype)) +
  scale_x_discrete(labels = CONTRAST_LABEL_WRAP) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, limits = c(-nes_lim, nes_lim), name = "NES") +
  scale_size_continuous(range = c(1, 5), limits = c(0, SIZE_CAP), name = expression(-log[10]~FDR)) +
  scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = "black"), labels = c(`FALSE` = sprintf("FDR >= %s", FDR_CUT), `TRUE` = sprintf("FDR < %s", FDR_CUT)),
                      name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 3, fill = "white", stroke = .8))) +
  labs(x = NULL, y = NULL, title = "Pathways of interest in ageing and HFpEF (predefined panel)",
       subtitle = sprintf("NES and FDR from the full GSEA (BH within database x subtype x contrast); black ring = FDR < %s; black centre dot = single-animal candidate (one animal > %d%% of the high-side set signal);\ngrey dash = not tested (fewer than %d genes in the ranking universe); size capped at %d",
                          FDR_CUT, MAX_ANIMAL_SHARE * 100, MIN_SIZE, SIZE_CAP)) +
  theme_house(base_size = 9) +
  theme(axis.text.x = element_text(size = 6), panel.spacing.x = unit(8, "pt"), axis.text.y = element_text(size = 7.5), strip.text.y = element_text(angle = 0, hjust = 0, size = 7.5),
        panel.spacing.y = unit(2, "pt"), legend.position = "right")
save_svg(fig, "ec_gsea_pathway_panel.svg", width = 16, height = 10.5, dir = path_figures)
ggsave(file.path(path_figures, "ec_gsea_pathway_panel.png"), plot = fig, width = 16, height = 10.5, dpi = 150)
fig

panel_out <- panel_grid[order(celltype, contrast, theme, match(pathway, panel_def$pathway)),
                .(subtype = as.character(celltype), contrast = CONTRAST_LABEL[as.character(contrast)], theme = as.character(theme),
                  set = as.character(set), pathway, NES = round(NES, 2), FDR = signif(padj, 2), size,
                  leading_edge_top10 = vapply(leadingEdge, \(g) paste(head(g, 10), collapse = ", "), character(1)),
                  single_animal_candidate, status = ifelse(tested, "tested", "not tested"))]
fwrite(panel_out, file.path(path_out, "ec_gsea_pathway_panel.csv"))
dt_table(panel_out[, -"pathway"], "Predefined pathway panel: NES and FDR from the full GSEA; leading edge = first 10 genes in ranking order.")
age_contrasts <- c("OC_vs_YC", "OH_vs_YH"); hfpef_contrasts <- c("OH_vs_OC", "YH_vs_YC")
both_age <- panel_grid[contrast %in% age_contrasts & sig == TRUE, .(n = .N, dir = unique(sign(NES))[1], same = uniqueN(sign(NES)) == 1), by = .(celltype, set)][
  n == 2 & same == TRUE, .(subtypes = .N), by = .(set, dir)][subtypes >= 3][order(dir, -subtypes)]
sig_by_contrast <- panel_grid[, .(n = sum(sig)), by = contrast][, setNames(n, as.character(contrast))]
ifn_age_control    <- panel_grid[contrast == "OC_vs_YC" & sig == TRUE & grepl("INTERFERON", pathway)]
n_tested  <- panel_grid[, sum(tested)]
lower_first <- function(x) { x <- as.character(x); ifelse(grepl("^[A-Z][a-z]", x), paste0(tolower(substr(x, 1, 1)), substring(x, 2)), x) }

results_export <- gsea[, .(geneset, panel, contrast, celltype, pathway, NES, ES, pval, padj, padj_all, size, log2err, universe_n, tech_prop = tech_frac,
                technical_candidate, sad_share_prop = sad_share, sad_animal, single_animal_candidate, leadingEdge = vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1)))][
        order(geneset, match(contrast, CONTRAST_ORDER), match(celltype, SUBTYPE_ORDER), padj)]
csv_path <- file.path(path_out, "ec_gsea_mast_results.csv"); fwrite(results_export, csv_path); saved_to(csv_path)
collapsed_sets <- rbindlist(COLLAPSED, idcol = "geneset")
fwrite(collapsed_sets, file.path(path_out, "ec_gsea_collapsed_sets.csv"))
cat(sprintf("%d rows (%d significant at FDR < %.2f; %d of these technical candidates); %d sets collapsed into representatives in the figures\n",
            nrow(results_export), sum(results_export$padj < FDR_CUT, na.rm = TRUE), FDR_CUT, sum(results_export$padj < FDR_CUT & results_export$technical_candidate, na.rm = TRUE), nrow(collapsed_sets)))
dt_table(head(results_export[padj < FDR_CUT][order(padj)][, .(geneset, panel, celltype, pathway, NES = round(NES, 2), padj = signif(padj, 2), size,
                                                   technical_candidate)], TABLE_ROWS),
         sprintf("Top %d subtype x pathway enrichments (FDR < %s), all databases.", TABLE_ROWS, FDR_CUT))

page_plot <- function(page_dt, geneset_db, subtype_name, contrast_name, i, n_pages) {
  path_labels <- make.unique(pretty_lab(page_dt$pathway, geneset_db, width = 55))
  page_dt[, pathway_lab := factor(path_labels, levels = rev(path_labels))]
  page_dt[, class := fifelse(technical_candidate, "technical candidate", fifelse(NES > 0, "NES > 0", "NES < 0"))]
  ggplot(page_dt, aes(NES, pathway_lab, fill = class, size = pmin(lfdr, SIZE_CAP))) +
    geom_vline(xintercept = 0, colour = "grey75", linetype = "dashed") +
    geom_point(shape = 21, colour = "grey25", stroke = 0.3) +
    single_animal_dot(page_dt) +
    scale_fill_manual(values = c(`NES > 0` = "#B2182B", `NES < 0` = "#2166AC", `technical candidate` = TECH_GREY), name = NULL) +
    scale_size_continuous(range = c(2, 6), limits = c(0, SIZE_CAP), name = expression(-log[10]~FDR)) +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 3))) +
    labs(title = sprintf("%s \u2014 %s \u2014 %s", geneset_db, subtype_name, CONTRAST_LABEL[[contrast_name]]),
         subtitle = sprintf("FDR < %s | universe %s genes | page %d of %d | black centre dot = single-animal candidate", FDR_CUT, format(page_dt$universe_n[1], big.mark = ","), i, n_pages),
         x = "Normalized enrichment score (NES)", y = NULL) +
    theme_house() + theme(axis.text.y = element_text(size = 8))
}
pdf_files <- vapply(names(GS_LIST), function(geneset_db) {
  path_pdf <- file.path(path_figures, sprintf("ec_gsea_%s_per_comparison.pdf", tolower(gsub("[^A-Za-z]", "", geneset_db))))
  cairo_pdf(path_pdf, width = 11, height = 8.5, onefile = TRUE)
  for (subtype_name in SUBTYPE_ORDER) for (contrast_name in CONTRASTS_RUN) {
    page_dt <- gsea[geneset == geneset_db & celltype == subtype_name & contrast == contrast_name & sig & is.finite(NES)][order(padj, -abs(NES), pathway)]
    page_rows <- split(seq_len(nrow(page_dt)), ceiling(seq_len(nrow(page_dt)) / TERMS_PER_PAGE))
    for (i in seq_along(page_rows)) print(page_plot(page_dt[page_rows[[i]]], geneset_db, subtype_name, contrast_name, i, length(page_rows)))
  }
  dev.off(); path_pdf
}, character(1))
invisible(lapply(pdf_files, saved_to))

session_footer()
