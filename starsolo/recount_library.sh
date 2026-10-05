#!/usr/bin/env bash
# Intron-inclusive recount of one library from its Cell Ranger BAM with STARsolo.
# Usage: recount_library.sh <zip_file> <path_of_bam_inside_zip> <library_name>
# Environment: conda env "starsolo" (STAR 2.7.11b, samtools 1.21, 7-Zip, python with anndata/h5py/scipy).
# Only the converted outputs are kept (data/raw/<library>/); the extracted BAM and STAR scratch are deleted.
set -euo pipefail

zip_file=$1; bam_in_zip=$2; library=$3
path_v3=${PROJECT_ROOT:?set PROJECT_ROOT to the project folder}
path_star_index=${STAR_INDEX:?set STAR_INDEX to the mm10-2020-A STAR index}
path_scratch=${SCRATCH:?set SCRATCH to a scratch folder on the Linux disk}/$library
path_library_out=$path_v3/data/raw/$library
threads=8

mkdir -p "$path_scratch" "$path_library_out/starsolo_logs"
echo "[$(date +%H:%M)] $library: extracting BAM"
# ZIP_PASSWORD: archive password, set in the environment only. The delivered archives carry ~100 kB of data after the
# end of the zip, which 7-Zip reports with exit code 2; the extracted BAM is therefore verified against its stored CRC32.
7za e -y -bd ${ZIP_PASSWORD:+"-p$ZIP_PASSWORD"} "$zip_file" "$bam_in_zip" -o"$path_scratch" > /dev/null || true
bam_file=$path_scratch/$(basename "$bam_in_zip")
crc_zip=$( (7za l -slt ${ZIP_PASSWORD:+"-p$ZIP_PASSWORD"} "$zip_file" "$bam_in_zip" || true) | awk -F' = ' '/^CRC = /{print $2}')
crc_bam=$(7za h -scrcCRC32 "$bam_file" | awk '/for data:/{print $NF}')
[ "$crc_zip" = "$crc_bam" ] || { echo "$library: extracted BAM fails its CRC32 check ($crc_bam vs $crc_zip)"; exit 1; }

# Primary and unmapped reads that carry a corrected cell barcode (CB) and UMI (UB), exactly the reads Cell Ranger counts from
cat > "$path_scratch/bam_reads.sh" <<'EOF'
#!/usr/bin/env bash
samtools view -@ 2 -F 0x900 -e '[CB] && [UB]' "$1"
EOF
chmod +x "$path_scratch/bam_reads.sh"

# Whitelist = the corrected cell barcodes present in the BAM, matched exactly. Cell Ranger has already corrected them, so the
# counts equal those without a whitelist, and every STARsolo feature (Gene, GeneFull_Ex50pAS, Velocyto) shares one barcode list
# (without a whitelist the Velocyto matrix columns do not follow its barcodes.tsv).
"$path_scratch/bam_reads.sh" "$bam_file" | grep -o 'CB:Z:[ACGT]\{16\}' | cut -c6- | sort -u -S 4G -T "$path_scratch" > "$path_scratch/whitelist.txt"

echo "[$(date +%H:%M)] $library: STARsolo"
# The barcode "read" is CB + UB concatenated (CCCCCCCCCCCCCCCC-1UUUUUUUUUUUU): CB = positions 1-16, UMI = 19-30.
# CB and UB are already corrected by Cell Ranger (exact whitelist match, no correction); UMI and cell-calling rules follow Cell Ranger.
STAR --runThreadN $threads --genomeDir "$path_star_index" \
     --readFilesIn "$bam_file" --readFilesType SAM SE --readFilesCommand "$path_scratch/bam_reads.sh" \
     --soloType CB_UMI_Simple --soloInputSAMattrBarcodeSeq CB UB --soloInputSAMattrBarcodeQual - \
     --soloCBwhitelist "$path_scratch/whitelist.txt" --soloCBmatchWLtype Exact --soloCBstart 1 --soloCBlen 16 --soloUMIstart 19 --soloUMIlen 12 --soloBarcodeReadLength 0 \
     --soloStrand Forward --soloFeatures Gene GeneFull_Ex50pAS Velocyto \
     --soloUMIdedup 1MM_CR --soloUMIfiltering MultiGeneUMI_CR --soloCellFilter EmptyDrops_CR \
     --clipAdapterType CellRanger4 --outFilterScoreMin 30 \
     --outSAMtype None --outFileNamePrefix "$path_scratch/star/" --outTmpDir "$path_scratch/star_tmp"

echo "[$(date +%H:%M)] $library: converting outputs"
python "$path_v3/starsolo/starsolo_to_10x.py" "$path_scratch/star/Solo.out" "$path_library_out"
cp "$path_scratch/star/Log.final.out" "$path_library_out/starsolo_logs/"
for feature_set in Gene GeneFull_Ex50pAS Velocyto; do
  cp "$path_scratch/star/Solo.out/$feature_set/Summary.csv" "$path_library_out/starsolo_logs/${feature_set}_Summary.csv" 2>/dev/null || true
done
cp "$path_scratch/star/Solo.out/Barcodes.stats" "$path_library_out/starsolo_logs/"

rm -rf "$path_scratch"
echo "[$(date +%H:%M)] $library: done"
