"""Compare the STARsolo exon-only recount of every library with its original Cell Ranger 6.0.2 filtered matrix.

Inputs per library (data/raw/<library>/): cellranger/filtered_feature_bc_matrix.h5 (Cell Ranger 6.0.2, include-introns = false),
exonic_raw_feature_bc_matrix.h5 (STARsolo Gene, exon-only) and raw_feature_bc_matrix.h5 (STARsolo GeneFull_Ex50pAS, intron-inclusive).
Outputs (data/raw/): recount_validation.csv (one row per library) and recount_validation_per_nucleus.csv.gz (Cell Ranger nuclei:
UMIs in Cell Ranger, exon-only and intron-inclusive counts).
Usage: python validate_recount.py   (env starsolo)
"""
import h5py, numpy as np, pandas as pd, scipy.sparse as sp
import os
from pathlib import Path

RAW = Path(os.environ["PROJECT_ROOT"]) / "data" / "raw"

def read_10x(path):
    with h5py.File(path, "r") as h5:
        m = h5["matrix"]
        counts = sp.csc_matrix((m["data"][:], m["indices"][:], m["indptr"][:]), shape=m["shape"][:])
        return counts, m["barcodes"][:].astype(str), m["features"]["id"][:].astype(str)

rows, per_nucleus = [], []
for lib_dir in sorted(p for p in RAW.iterdir() if (p / "cellranger" / "filtered_feature_bc_matrix.h5").exists()):
    cr, cr_bc, cr_genes = read_10x(lib_dir / "cellranger" / "filtered_feature_bc_matrix.h5")
    ex, ex_bc, ex_genes = read_10x(lib_dir / "exonic_raw_feature_bc_matrix.h5")
    full, full_bc, full_genes = read_10x(lib_dir / "raw_feature_bc_matrix.h5")
    assert (cr_genes == ex_genes).all() and (ex_genes == full_genes).all()
    cells = np.intersect1d(cr_bc, ex_bc)
    cr_s, ex_s, full_s = (m[:, pd.Index(b).get_indexer(cells)] for m, b in ((cr, cr_bc), (ex, ex_bc), (full, full_bc)))
    umi_cr, umi_ex, umi_full = (np.asarray(m.sum(0)).ravel() for m in (cr_s, ex_s, full_s))
    gene_cr, gene_ex = (np.asarray(m.sum(1)).ravel() for m in (cr_s, ex_s))
    expressed = (gene_cr + gene_ex) >= 10
    rows.append({"library": lib_dir.name, "cellranger_nuclei": len(cr_bc), "nuclei_in_recount": len(cells),
                 "median_ratio_exon_only_to_cellranger": np.median(umi_ex / np.maximum(umi_cr, 1)),
                 "pearson_r_nucleus_log1p_umis": np.corrcoef(np.log1p(umi_cr), np.log1p(umi_ex))[0, 1],
                 "pearson_r_gene_log1p_totals": np.corrcoef(np.log1p(gene_cr[expressed]), np.log1p(gene_ex[expressed]))[0, 1],
                 "median_ratio_intron_inclusive_to_cellranger": np.median(umi_full / np.maximum(umi_cr, 1)),
                 "median_umis_cellranger": np.median(umi_cr), "median_umis_intron_inclusive": np.median(umi_full)})
    per_nucleus.append(pd.DataFrame({"library": lib_dir.name, "barcode": cells, "umis_cellranger": umi_cr,
                                     "umis_exon_only": umi_ex, "umis_intron_inclusive": umi_full}))
    print(lib_dir.name, len(cells), flush=True)
pd.DataFrame(rows).round(4).to_csv(RAW / "recount_validation.csv", index=False)
pd.concat(per_nucleus).to_csv(RAW / "recount_validation_per_nucleus.csv.gz", index=False)
