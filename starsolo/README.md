# Intron-inclusive recount (STARsolo) and ambient-RNA removal (CellBender)

Scripts that turn the Cell Ranger BAM of each library into the count matrices read by notebook 01:
intron-inclusive counts, exon-only counts, spliced / unspliced / ambiguous layers, the exonic read fraction
per nucleus, and the CellBender-corrected matrix. Versions of every tool are in `env/software_versions.csv`;
the two conda environments are in `env/`.

## Files

| file | purpose |
|---|---|
| `build_reference.sh` | downloads the 10x reference refdata-gex-mm10-2020-A and builds the STAR index |
| `recount_library.sh` | one library: extracts its BAM from the zip, runs STARsolo, converts the output, deletes the BAM and scratch |
| `starsolo_to_10x.py` | converts STARsolo matrices into 10x HDF5 / h5ad files and computes `intron_metrics.csv` |
| `cellbender_library.sh` | one library: CellBender remove-background on the intron-inclusive raw matrix (GPU) |
| `run_zip.sh` | runs both steps for every library in one zip; restartable (finished libraries are skipped) |
| `library_map.csv` | Cell Ranger library ID -> library name used in the analysis (private; not included) |
| `env/environment_starsolo.yml`, `env/environment_cellbender.yml` | conda environments |
| `env/software_versions.csv` | versions of all software, reference and hardware |

## How to run

```bash
conda env create -f env/environment_starsolo.yml
conda env create -f env/environment_cellbender.yml
conda run -n starsolo bash build_reference.sh /path/to/ref          # once; ~40 GB disk, ~32 GB RAM
export STAR_INDEX=/path/to/ref/star_index SCRATCH=/path/to/scratch  # scratch needs space for one BAM
export ZIP_PASSWORD=...   # archive password, if the zip is encrypted (keep it out of files and logs)
bash run_zip.sh /path/to/AMC_S182.zip
```

Outputs per library in `data/raw/<library>/`:

| file | content |
|---|---|
| `raw_feature_bc_matrix.h5` | intron-inclusive UMI counts, all barcodes (10x HDF5) |
| `filtered_feature_bc_matrix.h5` | same, barcodes called by STARsolo EmptyDrops_CR |
| `exonic_raw_feature_bc_matrix.h5` | exon-only UMI counts, all barcodes |
| `spliced_unspliced.h5ad` | spliced / unspliced / ambiguous layers (Velocyto convention), called barcodes |
| `intron_metrics.csv` | per barcode: exon-only UMIs (`Gene`), total UMIs (`GeneFull_Ex50pAS`), `exon_prop` = exon-only / total, and spliced / unspliced / ambiguous molecules |
| `cellbender_filtered.h5`, `cellbender_cell_barcodes.csv`, `cellbender_metrics.csv`, `cellbender_report.html` | CellBender outputs |
| `cellranger/` | the original Cell Ranger 6.0.2 `metrics_summary.csv` and `web_summary.html` |
| `starsolo_logs/` | STAR `Log.final.out`, STARsolo `Summary.csv` per feature set, `Barcodes.stats`, CellBender log |

## Methods (text for the manuscript)

**Intron-inclusive quantification.** Libraries had originally been quantified with Cell Ranger 6.0.2 (refdata-gex-mm10-2020-A, introns not
included), which counts only reads compatible with exons; a large share of reads from nuclei maps to introns. Because most reads from nuclei originate from unspliced pre-mRNA, each library
was re-quantified from its Cell Ranger BAM with STARsolo (STAR 2.7.11b)^1,2^ against the same reference
(10x Genomics refdata-gex-mm10-2020-A: GRCm38 with the GENCODE vM23 / Ensembl 98 annotation, 32,285 genes;
index built with `--sjdbOverhang 90`). All primary and unmapped reads carrying a corrected cell barcode (CB tag)
and UMI (UB tag) were re-aligned (`--readFilesType SAM SE`; reads stored on the reverse strand are
reverse-complemented by STAR), using the Cell Ranger barcode and UMI corrections recorded in the BAM
(matched exactly to the corrected barcodes present in each BAM: `--soloCBwhitelist <BAM barcodes> --soloCBmatchWLtype Exact`, which
gives every feature the same barcode list), Cell Ranger-equivalent UMI collapsing and filtering (`--soloUMIdedup 1MM_CR`,
`--soloUMIfiltering MultiGeneUMI_CR`), Cell Ranger 4 adapter and poly(A) clipping (`--clipAdapterType CellRanger4`),
`--outFilterScoreMin 30` and sense-strand counting (`--soloStrand Forward`). Reads were counted as in Cell Ranger ≥ 7
with introns included (`--soloFeatures GeneFull_Ex50pAS`), as exon-only counts (`Gene`) and as spliced, unspliced and
ambiguous molecules (`Velocyto`)^3^. Non-empty droplets were called with STARsolo's EmptyDrops_CR
(`--soloCellFilter EmptyDrops_CR`)^4^. For each barcode the fraction of reads mapping exclusively to exons (`exon_prop` = exon-only UMIs / intron-inclusive
UMIs) was computed as a measure of cytoplasmic carry-over, following Simonson et al.^5^

**Ambient RNA removal.** Ambient RNA was removed from the intron-inclusive raw matrices with CellBender
remove-background 0.4.0^6^ on a GPU (PyTorch 2.8, CUDA 12.8) with 6,000 expected cells, 25,000 droplets included,
a target false-positive rate of 0.01, 150 epochs, a latent dimension of 100 and a learning rate of 1 × 10^-4^, following
the settings of Simonson et al.^5^

**References**
1. Dobin A, et al. STAR: ultrafast universal RNA-seq aligner. *Bioinformatics* 2013;29:15–21. doi:10.1093/bioinformatics/bts635
2. Kaminow B, Yunusov D, Dobin A. STARsolo: accurate, fast and versatile mapping/quantification of single-cell and single-nucleus RNA-seq data. *bioRxiv* 2021. doi:10.1101/2021.05.05.442755
3. La Manno G, et al. RNA velocity of single cells. *Nature* 2018;560:494–498. doi:10.1038/s41586-018-0414-6
4. Lun ATL, et al. EmptyDrops: distinguishing cells from empty droplets in droplet-based single-cell RNA sequencing data. *Genome Biol* 2019;20:63. doi:10.1186/s13059-019-1662-y
5. Simonson B, Chaffin M, et al. Single-nucleus RNA sequencing in ischemic cardiomyopathy reveals common transcriptional profile underlying end-stage heart failure. *Cell Rep* 2023;42:112086. doi:10.1016/j.celrep.2023.112086
6. Fleming SJ, et al. Unsupervised removal of systematic background noise from droplet-based single-cell experiments using CellBender. *Nat Methods* 2023;20:1323–1335. doi:10.1038/s41592-023-01943-7

## Validation against the original Cell Ranger counts
`python validate_recount.py` (env starsolo) compares, for every library, the STARsolo exon-only recount with the original Cell Ranger 6.0.2
filtered matrix (`data/raw/<library>/cellranger/filtered_feature_bc_matrix.h5`) and with the intron-inclusive counts. Outputs:
`data/raw/recount_validation.csv` (one row per library) and `data/raw/recount_validation_per_nucleus.csv.gz` (UMIs per Cell Ranger nucleus).

