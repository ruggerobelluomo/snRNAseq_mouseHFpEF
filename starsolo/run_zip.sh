#!/usr/bin/env bash
# Recount (STARsolo) and CellBender every library contained in one zip of Cell Ranger outputs.
# Usage: STAR_INDEX=<index> SCRATCH=<scratch folder> [ZIP_PASSWORD=<password>] run_zip.sh <zip_file>
# The archive password is passed only through the environment, never written to a file or log.
# Libraries are taken from the possorted_genome_bam.bam entries of the zip; their Cell Ranger library IDs
# (the folder name above outs/) are translated to library names with library_map.csv.
# Libraries whose outputs already exist are skipped, so the script can be restarted.
set -euo pipefail

zip_file=$1
path_scripts=$(cd "$(dirname "$0")" && pwd)
path_raw=${PROJECT_ROOT:?set PROJECT_ROOT to the project folder}/data/raw

bam_list=$( (7za l -slt ${ZIP_PASSWORD:+"-p$ZIP_PASSWORD"} "$zip_file" || true) | awk -F' = ' '/^Path = .*possorted_genome_bam\.bam$/ {print $2}')
echo "$bam_list" | while read -r bam_in_zip; do
  cellranger_id=$(echo "$bam_in_zip" | tr '\\' '/' | awk -F/ '{for (i = 1; i <= NF; i++) if ($i == "outs") print $(i-1)}')
  library=$(awk -F, -v id="$cellranger_id" '$1 == id {print $2}' "$path_scripts/library_map.csv")
  [ -n "$library" ] || { echo "no library name for '$cellranger_id' ($bam_in_zip)"; exit 1; }
  [ -f "$path_raw/$library/raw_feature_bc_matrix.h5" ] || conda run -n starsolo bash "$path_scripts/recount_library.sh" "$zip_file" "$bam_in_zip" "$library"
  [ -f "$path_raw/$library/cellbender_filtered.h5" ]   || conda run -n cellbender bash "$path_scripts/cellbender_library.sh" "$library"
done
