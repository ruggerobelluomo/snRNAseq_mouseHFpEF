
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(svglite); library(dplyr); library(tidyr) })

if (!exists("FIG_DIR")) {
  FIG_DIR <- normalizePath(".")
  if (basename(FIG_DIR) == "code") FIG_DIR <- dirname(FIG_DIR)
}
FIG_DIR      <- normalizePath(FIG_DIR)
if (!exists("PROJECT_ROOT") || !nzchar(PROJECT_ROOT)) PROJECT_ROOT <- file.path(FIG_DIR, "..", "..")
PROJECT_ROOT <- normalizePath(PROJECT_ROOT)
FIG_ROOT     <- PROJECT_ROOT
stopifnot(file.exists(file.path(PROJECT_ROOT, "notebooks", "_common.R")))
source(file.path(PROJECT_ROOT, "notebooks", "_common.R"))
RES      <- function(...) file.path(PROJECT_ROOT, "results", ...)
OBJ      <- function(...) file.path(PROJECT_ROOT, "data", "objects", ...)
RESOURCE <- function(...) file.path(PROJECT_ROOT, "data", "resources", ...)
FIG_OUT     <- FIG_DIR
PREVIEW_DIR <- file.path(FIG_DIR, "previews")
TABLE_DIR   <- file.path(FIG_DIR, "tables")
PANEL_DIR   <- file.path(FIG_DIR, "panels")
for (d in c(PREVIEW_DIR, TABLE_DIR, PANEL_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

GRP_ORDER <- GROUP_ORDER; GRP_LABEL <- GROUP_LABEL; GRP_FILL <- GROUP_FILL
CO_ORDER  <- CONTRAST_ORDER
CT_ORDER  <- CELL_TYPE_ORDER; CT_FILL <- CELL_TYPE_FILL; CT_LABEL <- CELL_TYPE_LABEL; CT_SHORT <- CELL_TYPE_SHORT
star_fdr  <- function(p) as.character(fdr_stars(p))

MM         <- 1 / 25.4
WIDTH_1COL <- 85
WIDTH_2COL <- 180
MAX_HEIGHT <- 230
BASE_PT    <- 7
TITLE_PT   <- 8
LETTER_PT  <- 10
MIN_PT     <- 6
FONT       <- "sans"
RASTER_DPI <- 600
pt_text    <- function(pt = BASE_PT) pt / .pt

raster_px <- function(w_mm, h_mm, dpi = RASTER_DPI) round(c(w_mm, h_mm) * MM * dpi)

AXIS_GREY <- "grey45"
GRID_GREY <- "grey93"
theme_fig <- function(grid = c("none", "x", "y", "both")) {
  grid <- match.arg(grid)
  gl <- element_line(colour = GRID_GREY, linewidth = 0.2)
  theme_classic(base_size = BASE_PT, base_family = FONT) +
    theme(axis.line = element_line(linewidth = 0.25, colour = AXIS_GREY),
          axis.ticks = element_line(linewidth = 0.25, colour = AXIS_GREY), axis.ticks.length = unit(0.8, "mm"),
          axis.text = element_text(size = BASE_PT, colour = "black"), axis.title = element_text(size = TITLE_PT),
          panel.border = element_blank(), panel.background = element_blank(),
          panel.grid.major.x = if (grid %in% c("x", "both")) gl else element_blank(),
          panel.grid.major.y = if (grid %in% c("y", "both")) gl else element_blank(),
          strip.background = element_blank(), strip.text = element_text(size = TITLE_PT, face = "bold"),
          legend.title = element_text(size = BASE_PT), legend.text = element_text(size = BASE_PT),
          legend.key.size = unit(3, "mm"), legend.background = element_blank(),
          plot.title = element_blank(), plot.subtitle = element_blank(),
          plot.caption = element_text(size = BASE_PT, hjust = 0, colour = "grey20"),
          plot.margin = margin(1.5, 1.5, 1.5, 1.5, "mm"))
}
theme_set(theme_fig())

FLAG_OF    <- sub(":.*$", "", QC_FLAGGED_ANIMALS)
FLAG_LABEL <- c(low_yield_library = "low-yield library", manual_doublet_threshold = "manual doublet threshold")
animal_tag <- function(sample) {
  out <- rep("", length(sample))
  i <- sample %in% names(FLAG_OF); out[i] <- paste0(sample[i], " (", FLAG_LABEL[FLAG_OF[sample[i]]], ")")
  j <- sample %in% names(ANIMAL_NOTES); out[j] <- paste0(sample[j], " (animal note)")
  out
}
fmt_n   <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
fmt_p   <- function(p, digits = 2) ifelse(p < 0.001, formatC(p, format = "e", digits = 1), formatC(signif(p, digits), format = "fg", digits = digits))
fmt_pct <- function(x, digits = 1) formatC(100 * x, format = "f", digits = digits)

input_stamp <- function(...) {
  f <- unlist(list(...), use.names = FALSE)
  stopifnot(all(file.exists(f)))
  paste(f, file.size(f), format(file.mtime(f), "%Y-%m-%d %H:%M:%OS3"), sep = "|")
}

save_figure <- function(fig, name, width_mm = WIDTH_2COL, height_mm, tags = TRUE, preview_dpi = 200) {
  stopifnot(width_mm <= WIDTH_2COL + 0.01, height_mm <= MAX_HEIGHT + 0.01)
  if (tags) fig <- fig + plot_annotation(tag_levels = "A")
  fig <- fig & theme(plot.tag = element_text(size = LETTER_PT, face = "bold", family = FONT))
  svg_file <- file.path(FIG_OUT, paste0(name, ".svg"))
  svglite::svglite(svg_file, width = width_mm * MM, height = height_mm * MM)
  print(fig); invisible(dev.off())
  svg <- readLines(svg_file, warn = FALSE)
  writeLines(gsub("font-family: [^;']*", "font-family: Arial, Helvetica, sans-serif", svg), svg_file)
  grDevices::cairo_pdf(file.path(FIG_OUT, paste0(name, ".pdf")), width = width_mm * MM, height = height_mm * MM, family = "sans")
  print(fig); invisible(dev.off())
  ragg::agg_png(file.path(PREVIEW_DIR, paste0(name, ".png")), width = width_mm * MM, height = height_mm * MM, units = "in", res = preview_dpi)
  print(fig); invisible(dev.off())
  fs <- svg_font_sizes(svg_file)
  if (min(fs) < MIN_PT - 0.01) warning(sprintf("%s: smallest SVG text %.2f pt (< %d pt)", name, min(fs), MIN_PT))
  invisible(data.frame(figure = name, width_mm = width_mm, height_mm = height_mm,
                       min_font_pt = min(fs), max_font_pt = max(fs), svg_kb = round(file.size(svg_file) / 1024)))
}

svg_font_sizes <- function(svg_file) {
  svg <- paste(readLines(svg_file, warn = FALSE), collapse = "\n")
  as.numeric(regmatches(svg, gregexpr("(?<=font-size: )[0-9.]+(?=px)", svg, perl = TRUE))[[1]])
}

write_table <- function(df, name, sheet = name) {
  csv <- file.path(TABLE_DIR, paste0(name, ".csv"))
  utils::write.csv(df, csv, row.names = FALSE, na = "")
  if (requireNamespace("openxlsx", quietly = TRUE)) openxlsx::write.xlsx(df, file.path(TABLE_DIR, paste0(name, ".xlsx")), sheetName = substr(sheet, 1, 31), overwrite = TRUE)
  invisible(csv)
}
