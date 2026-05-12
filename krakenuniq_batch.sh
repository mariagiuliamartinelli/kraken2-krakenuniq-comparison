#!/bin/bash
#SBATCH -A naiss2025-5-616
#SBATCH -p memory
#SBATCH --mem=650GB
#SBATCH -t 06:00:00
#SBATCH -J krakenuniq_batch

ml krakenuniq/1.0.4

MIC_NT="/cfs/klemming/projects/supr/archaeogenetics/databases/DBDIR_KrakenUniq_MicrobialNT"
INPUT_LIST="$1"  # Path to a text file with one FASTQ file per line
OUTPUT_DIR="/cfs/klemming/projects/supr/naiss2025-6-437/results/metaJAM/Fiskartjarn/krakenuniq/"

mkdir -p "$OUTPUT_DIR"

while IFS= read -r INPUT; do
    BASENAME=$(basename "$INPUT" _merged.dust.rmdup.fastq.gz)
    krakenuniq --preload \
      --db "$MIC_NT" \
      --fastq-input "$INPUT" \
      --threads 256 \
      --output "$OUTPUT_DIR/${BASENAME}.krakenuniq.sequences" \
      --report-file "$OUTPUT_DIR/${BASENAME}.krakenuniq.output" \
      --gzip-compressed \
      --only-classified-out \
      &> "$OUTPUT_DIR/${BASENAME}.krakenuniq.log"
done < "$INPUT_LIST"
