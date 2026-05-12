#!/bin/bash
#SBATCH -A naiss2025-5-616
#SBATCH -p memory
#SBATCH --mem=650GB
#SBATCH -t 15:00:00
#SBATCH -J krakenuniq

ml krakenuniq/1.0.4

MIC_NT="/cfs/klemming/projects/supr/archaeogenetics/databases/DBDIR_KrakenUniq_MicrobialNT"
INPUT="/cfs/klemming/projects/supr/naiss2025-6-437/results/metaJAM/new_metaJAM/metaJAM/FISK2006_deep/02_sga/FISK2006_merged.dust.rmdup.fastq.gz"
OUTPUT_DIR="/cfs/klemming/projects/supr/naiss2025-6-437/results/metaJAM/new_metaJAM/metaJAM/FISK2006_deep/krakenuniq/output"

mkdir -p "$OUTPUT_DIR"

krakenuniq --preload \
  --db "$MIC_NT" \
  --fastq-input "$INPUT" \
  --threads 256 \
  --output "$OUTPUT_DIR/krakenuniq.sequences" \
  --report-file "$OUTPUT_DIR/krakenuniq.output" \
  --gzip-compressed \
  --only-classified-out \
  &> "$OUTPUT_DIR/krakenuniq.log"