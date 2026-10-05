# make_figures.R
# Code-only export of the analysis pipeline. Data are not included; see README.md.

params <- list(fig_dir = "manuscript/figures", project_root = "", run_main = TRUE, run_supp = TRUE)


t_knit_start <- Sys.time()
FIG_DIR <- normalizePath(params$fig_dir)
if (nzchar(params$project_root)) PROJECT_ROOT <- params$project_root
knitr::opts_knit$set(root.dir = FIG_DIR)
source(file.path(FIG_DIR, "code", "figure_style.R"))

if (!isTRUE(l10n_info()[["UTF-8"]])) stop("knit in a UTF-8 locale, e.g. LANG=C.UTF-8 (current LC_CTYPE: ", Sys.getlocale("LC_CTYPE"), ")")
STYLE_MD5 <- c(unname(tools::md5sum(c(file.path(FIG_DIR, "code", "figure_style.R"), file.path(PROJECT_ROOT, "notebooks", "_common.R")))),
               Sys.getlocale("LC_CTYPE"))
knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE, cache = TRUE, autodep = TRUE, cache.lazy = FALSE,
                      cache.extra = STYLE_MD5, fig.show = "hide", results = "hold")

suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(dplyr); library(tidyr); library(Seurat); library(Matrix) })
theme_set(theme_fig())

for (fn in c("count", "filter", "select", "rename", "slice", "slice_max", "slice_min", "desc", "first", "last", "between",
             "lag", "lead", "n", "n_distinct", "summarise", "mutate", "arrange", "distinct", "group_by", "ungroup", "coalesce"))
  assign(fn, getExportedValue("dplyr", fn), envir = knitr::knit_global())

svgs <- list.files(FIG_OUT, pattern = "^Figure_.*\\.svg$", full.names = TRUE)
font_check <- bind_rows(lapply(svgs, function(f) {
  fs <- svg_font_sizes(f)
  x <- tryCatch(xml2::read_xml(f), error = function(e) NULL)
  data.frame(file = basename(f), min_font_pt = min(fs), max_font_pt = max(fs), n_text = length(fs),
             xml_ok = !is.null(x), fonts = paste(unique(regmatches(paste(readLines(f, warn = FALSE), collapse = ""),
                                                           gregexpr("font-family: [^;']*", paste(readLines(f, warn = FALSE), collapse = "")))[[1]]), collapse = " | "),
             size_kb = round(file.size(f) / 1024))
}))
knitr::kable(font_check, caption = "SVG checks: text sizes at print size (pt), XML validity, declared fonts")
stopifnot(all(font_check$xml_ok), all(font_check$min_font_pt >= MIN_PT - 0.01))

t_total <- as.numeric(difftime(Sys.time(), t_knit_start, units = "secs"))
cat(sprintf("knit wall-clock time: %.0f s\n", t_total))
time_objs <- ls(pattern = "^time_")
if (length(time_objs)) knitr::kable(data.frame(chunk = sub("^time_", "", time_objs), seconds = round(unlist(mget(time_objs)), 1)),
                                    caption = "Chunk run times at their last uncached run (s)")
fig_objs <- ls(pattern = "^figlog_")
if (length(fig_objs)) knitr::kable(bind_rows(mget(fig_objs)), caption = "Saved figures: print size (mm), SVG text-size range (pt), SVG size (kB)")

print_versions(c("ggplot2", "patchwork", "svglite", "scattermore", "Seurat", "CellChat", "knitr", "rmarkdown"))

if (isTRUE(params$run_main)) source(file.path(params$fig_dir, "figures_main.R"))
if (isTRUE(params$run_supp)) source(file.path(params$fig_dir, "figures_supplementary.R"))
