"""Convert one library's STARsolo output into the files the notebooks read.

Usage: python starsolo_to_10x.py <solo_out_dir> <library_out_dir>

Writes to <library_out_dir>:
  raw_feature_bc_matrix.h5        intron-inclusive counts (GeneFull_Ex50pAS), all barcodes, 10x HDF5 format
  filtered_feature_bc_matrix.h5   same, barcodes called by STARsolo EmptyDrops_CR
  exonic_raw_feature_bc_matrix.h5 exon-only counts (Gene), all barcodes
  spliced_unspliced.h5ad          Velocyto spliced / unspliced / ambiguous layers, called barcodes
  intron_metrics.csv              per barcode with >= 1 UMI: exon-only and total UMIs, exon_prop, Velocyto molecule counts
Barcodes carry the Cell Ranger suffix "-1" so they match the Cell Ranger barcodes of each library.
"""
import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.io

path_solo = Path(sys.argv[1])
path_out = Path(sys.argv[2])
path_out.mkdir(parents=True, exist_ok=True)
BARCODE_SUFFIX = "-1"
GENOME = "mm10"


def read_mtx_dir(path_dir):
    barcodes = (pd.read_csv(path_dir / "barcodes.tsv", header=None)[0].astype(str) + BARCODE_SUFFIX).values
    features = pd.read_csv(path_dir / "features.tsv", header=None, sep="\t")
    counts = scipy.io.mmread(path_dir / "matrix.mtx").tocsc().astype(np.int32)   # genes x barcodes
    return counts, barcodes, features


def write_10x_h5(counts, barcodes, features, path_h5):
    with h5py.File(path_h5, "w") as h5:
        matrix = h5.create_group("matrix")
        matrix.create_dataset("barcodes", data=np.array(barcodes, dtype="S"))
        matrix.create_dataset("data", data=counts.data, compression="gzip")
        matrix.create_dataset("indices", data=counts.indices.astype(np.int64), compression="gzip")
        matrix.create_dataset("indptr", data=counts.indptr.astype(np.int64))
        matrix.create_dataset("shape", data=np.array(counts.shape, dtype=np.int32))
        feature_group = matrix.create_group("features")
        feature_group.create_dataset("id", data=features[0].values.astype("S"))
        feature_group.create_dataset("name", data=features[1].values.astype("S"))
        feature_group.create_dataset("feature_type", data=np.array(["Gene Expression"] * len(features), dtype="S"))
        feature_group.create_dataset("genome", data=np.array([GENOME] * len(features), dtype="S"))
        feature_group.create_dataset("_all_tag_keys", data=np.array(["genome"], dtype="S"))


# Intron-inclusive counts (Cell Ranger >= 7 include-introns convention) and exon-only counts
counts_raw, barcodes_raw, features = read_mtx_dir(path_solo / "GeneFull_Ex50pAS" / "raw")
write_10x_h5(counts_raw, barcodes_raw, features, path_out / "raw_feature_bc_matrix.h5")
counts_called, barcodes_called, _ = read_mtx_dir(path_solo / "GeneFull_Ex50pAS" / "filtered")
write_10x_h5(counts_called, barcodes_called, features, path_out / "filtered_feature_bc_matrix.h5")
counts_exonic, barcodes_exonic, _ = read_mtx_dir(path_solo / "Gene" / "raw")
write_10x_h5(counts_exonic, barcodes_exonic, features, path_out / "exonic_raw_feature_bc_matrix.h5")

# Exonic read fraction (Simonson et al. 2023: proportion of reads mapping exclusively to exons) = exon-only UMIs (Gene) /
# intron-inclusive UMIs (GeneFull_Ex50pAS). With the exact whitelist all STARsolo features share one barcode list.
path_velocyto = path_solo / "Velocyto" / "raw"
velocyto_barcodes = (pd.read_csv(path_velocyto / "barcodes.tsv", header=None)[0].astype(str) + BARCODE_SUFFIX).values
assert (barcodes_exonic == barcodes_raw).all() and (velocyto_barcodes == barcodes_raw).all()
exonic_umis = np.asarray(counts_exonic.sum(axis=0)).ravel()
total_umis = np.asarray(counts_raw.sum(axis=0)).ravel()

# Spliced / unspliced / ambiguous molecules per barcode (Velocyto convention)
velocyto_features = pd.read_csv(path_velocyto / "features.tsv", header=None, sep="\t")
layers = {name: scipy.io.mmread(path_velocyto / f"{name}.mtx").T.tocsr().astype(np.int32)   # barcodes x genes
          for name in ["spliced", "unspliced", "ambiguous"]}
molecules = {name: np.asarray(matrix.sum(axis=1)).ravel() for name, matrix in layers.items()}
has_umis = total_umis > 0
pd.DataFrame({"barcode": barcodes_raw[has_umis],
              "exonic_umis": exonic_umis[has_umis],
              "total_umis": total_umis[has_umis],
              "exon_prop": exonic_umis[has_umis] / total_umis[has_umis],
              "spliced": molecules["spliced"][has_umis],
              "unspliced": molecules["unspliced"][has_umis],
              "ambiguous": molecules["ambiguous"][has_umis]}
             ).to_csv(path_out / "intron_metrics.csv", index=False)

is_called = np.isin(barcodes_raw, barcodes_called)
velocyto_obj = ad.AnnData(X=layers["spliced"][is_called],
                          obs=pd.DataFrame(index=barcodes_raw[is_called]),
                          var=pd.DataFrame({"gene_ids": velocyto_features[0].values}, index=velocyto_features[1].values),
                          layers={name: matrix[is_called] for name, matrix in layers.items()})
velocyto_obj.var_names_make_unique()
velocyto_obj.write_h5ad(path_out / "spliced_unspliced.h5ad", compression="gzip")

exon_prop_called = exonic_umis[is_called] / total_umis[is_called]
print(f"raw barcodes {counts_raw.shape[1]:,} | called {len(barcodes_called):,} | "
      f"median UMIs called {np.median(np.asarray(counts_called.sum(axis=0)).ravel()):.0f} | "
      f"median exon_prop called {np.median(exon_prop_called):.3f}")
