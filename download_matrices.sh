#!/bin/bash
# Scarica le due matrici (Kraken2 + KrakenUniq) da Dardel nella cartella CPG.
# Lanciare DAL TERMINAL DEL MAC (non da Dardel):
#     bash ~/Desktop/Università/Magistrale/CPG/download_matrices.sh

KEY="$HOME/.ssh/id-ed25519-pdc"
USER="mm1202"
HOST="dardel.pdc.kth.se"
REMOTE_BASE="/cfs/klemming/projects/supr/archaeogenetics/mm1202"
DEST="$HOME/Desktop/Università/Magistrale/CPG"

echo "[1/2] kraken2_family_matrix.tsv ..."
scp -i "$KEY" "${USER}@${HOST}:${REMOTE_BASE}/Kraken2_analysis/kraken2_family_matrix.tsv" "$DEST/"

echo "[2/2] krakenuniq_family_matrix.tsv ..."
scp -i "$KEY" "${USER}@${HOST}:${REMOTE_BASE}/KrakenUniq_analysis/krakenuniq_family_matrix.tsv" "$DEST/"

echo "Fatto. File salvati in: $DEST"
ls -lh "$DEST"/{kraken2,krakenuniq}_family_matrix.tsv 2>/dev/null
