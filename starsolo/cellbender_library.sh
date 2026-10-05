#!/usr/bin/env bash
# CellBender remove-background on one library's intron-inclusive raw matrix (GPU).
# Usage: cellbender_library.sh <library_name>
# Environment: conda env "cellbender" (CellBender 0.4.0, PyTorch 2.8 / CUDA 12.8).
# Settings follow Simonson et al. 2023 (Cell Rep): 6,000 expected cells, 25,000 droplets, FPR 0.01, 150 epochs,
# latent dimension 100, learning rate 1e-4.
set -euo pipefail

library=$1
path_library=${PROJECT_ROOT:?set PROJECT_ROOT to the project folder}/data/raw/$library
path_scratch=${SCRATCH:?set SCRATCH to a scratch folder on the Linux disk}/cellbender_$library
mkdir -p "$path_scratch" && cd "$path_scratch"

echo "[$(date +%H:%M)] $library: CellBender"
cellbender remove-background --cuda \
  --input "$path_library/raw_feature_bc_matrix.h5" --output "$path_scratch/cellbender.h5" \
  --expected-cells 6000 --total-droplets-included 25000 --fpr 0.01 --epochs 150 \
  --z-dim 100 --learning-rate 1e-4 > "$path_scratch/cellbender.log" 2>&1

cp cellbender_filtered.h5 "$path_library/cellbender_filtered.h5"
cp cellbender_cell_barcodes.csv "$path_library/cellbender_cell_barcodes.csv"
cp cellbender_metrics.csv "$path_library/cellbender_metrics.csv"
cp cellbender_report.html "$path_library/cellbender_report.html" 2>/dev/null || true
cp cellbender.log "$path_library/starsolo_logs/cellbender.log"
cd / && rm -rf "$path_scratch"
echo "[$(date +%H:%M)] $library: CellBender done"
