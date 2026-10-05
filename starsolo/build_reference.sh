#!/usr/bin/env bash
# Download the 10x Genomics mouse reference used by Cell Ranger (refdata-gex-mm10-2020-A: GRCm38 / Ensembl 98,
# GENCODE vM23 annotation, 32,285 genes) and build the STAR index used for the recount.
# Usage: build_reference.sh <output_folder>        (needs ~40 GB disk and ~32 GB RAM; ~20 min on 8 cores)
# Environment: conda env "starsolo" (env/environment_starsolo.yml).
set -euo pipefail

path_ref=$1
mkdir -p "$path_ref" && cd "$path_ref"
curl -L -o refdata-gex-mm10-2020-A.tar.gz https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz
tar -xzf refdata-gex-mm10-2020-A.tar.gz refdata-gex-mm10-2020-A/fasta refdata-gex-mm10-2020-A/genes refdata-gex-mm10-2020-A/reference.json
rm refdata-gex-mm10-2020-A.tar.gz
gunzip -f refdata-gex-mm10-2020-A/genes/genes.gtf.gz 2>/dev/null || true

# sjdbOverhang = read length - 1 (Chromium Single Cell 3' v3: 91-bp cDNA read)
mkdir -p star_index
STAR --runMode genomeGenerate --runThreadN 8 --genomeDir star_index \
     --genomeFastaFiles refdata-gex-mm10-2020-A/fasta/genome.fa \
     --sjdbGTFfile refdata-gex-mm10-2020-A/genes/genes.gtf --sjdbOverhang 90
