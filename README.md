# snRNA-seq of the ageing and HFpEF mouse heart: analysis code

Analysis code for a single-nucleus RNA-seq study of the left ventricle in a 2 x 2 design (young/old x control/HFpEF mice).

This is a **code-only** release. The sequencing data, processed objects, results and figures are not included, and the code cannot be run without them. Data will be made available on publication; until then, contact the authors.

## How this export was made
- The analysis was written as R Markdown and Jupyter notebooks. Each one is provided here as a plain script:
  - R Markdown notebooks were converted with `knitr::purl(documentation = 0)`.
  - The Jupyter notebooks contain only their code cells, with `# %%` cell markers that VS Code and Jupytext can open as notebooks.
- Narrative text, inline results, comments and notebook outputs were removed.
- Some study-specific inputs that were typed into the original code are now read from private files in `config/` (see `config/README.md`). These are the per-animal notes, the candidate gene list and the lab phenotyping tables.
- Figure-legend text generation is omitted.

## Pipeline (run in this order)
| Step | Script | Purpose |
|---|---|---|
| 0 | `starsolo/` | Recount each library from its Cell Ranger BAM with STARsolo (intron-inclusive counts), then CellBender ambient-RNA removal; see `starsolo/README.md` |
| 0 | `data/resources/build_reference_resources.R` | Build the reference resources (TF regulons, ligand-receptor table, secretome list, gene sets) |
| 1 | `notebooks/01_snRNAseq_QC_integration_annotation.py` | Per-nucleus QC, Scrublet doublet detection, Harmony integration, broad cell-type annotation (scanpy) |
| 1b | `notebooks/01b_integration_benchmark.py` | Benchmark of integration methods (scib-metrics) |
| 1c | `notebooks/export_seurat.R` | Convert the annotated AnnData atlas to a Seurat object |
| 2 | `notebooks/02_subclustering.R` | Per-lineage sub-clustering and atlas cleaning |
| 3 | `notebooks/03_Endothelial_fine_annotation.R` | Endothelial subtype annotation and composition (propeller/limma) |
| 4 | `notebooks/04_Cardiomyocyte_fine_annotation.R` | Cardiomyocyte subtype annotation and composition |
| 5 | `notebooks/05_EC_DE_analysis.R` | Endothelial differential expression (MAST) |
| 5b | `notebooks/05b_EC_pseudobulk_sensitivity.R` | Pseudobulk (DESeq2/edgeR) sensitivity analysis of the DE results |
| 5c | `notebooks/05c_EC_capillary_stratum_sensitivity.R` | Sensitivity of DE/GSEA to capillary sub-strata |
| 6 | `notebooks/06_EC_MAST_GSEA_analysis.R` | Gene-set enrichment (fgsea) on the MAST rankings |
| 7 | `notebooks/07_CCC_analysis.R` | Cell-cell communication (CellChat) between broad cell types, with animal-level tests |
| 7b | `notebooks/07b_EC_CM_subtype_CCC.R` | Endothelial-cardiomyocyte subtype-level communication |
| 8 | `notebooks/08_EC_KLF2_KLF4_shear.R` | Endothelial transcription-factor programme analysis (pseudobulk + decoupleR) |
| 9 | `notebooks/09_HFpEF_model_validation.R` | Comparison of snRNA-seq readouts with the phenotyping of the model |
| Figures | `manuscript/figures/make_figures.R` | Builds the main and supplementary figures (`figures_main.R`, `figures_supplementary.R`, style in `code/figure_style.R`) |

Shared settings and helpers for the R scripts are in `notebooks/_common.R`. The scripts expect to run from `notebooks/`, with the project root one level up (`PROJECT_ROOT`). The shell scripts in `starsolo/` read `PROJECT_ROOT` from the environment.

## Software
- Python: `env/environment_python.yml` (QC/integration); `env/environment_integration_bench.yml` (benchmark).
- R: package versions in `env/R_packages.csv`; install with `env/install_R_packages.R`.
- STARsolo/CellBender: `starsolo/env/`.

## License
MIT; see `LICENSE`.
