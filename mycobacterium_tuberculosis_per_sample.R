library(readxl)
library(dplyr)

mt <- read.table("mycobacterium_tuberculosis_per_sample.tsv", sep="\t", header=FALSE,
                 col.names=c("sample_raw","clade_reads","direct_reads","unique_kmers","coverage"))
mt$NGI_ID <- gsub("_S[0-9]+_L[0-9]+.*", "", mt$sample_raw)

# Alcune righe sono doppie (sottospecie) — tieni il massimo per campione
mt <- mt |>
  group_by(NGI_ID) |>
  slice_max(unique_kmers, n=1, with_ties=FALSE) |>
  ungroup()

mt_filt <- mt[mt$NGI_ID %in% metadata$NGI_ID, ] |>
  merge(metadata[, c("NGI_ID","Site","Period","Material")], by="NGI_ID") |>
  arrange(desc(unique_kmers))

head(mt_filt, 10)