library(readxl)
library(dplyr)

tf <- read.table("tannerella_forsythia_per_sample.tsv", sep="\t", header=FALSE,
                 col.names=c("sample_raw","clade_reads","direct_reads","unique_kmers","coverage"))
tf$NGI_ID <- gsub("_S[0-9]+_L[0-9]+.*", "", tf$sample_raw)
tf_filt <- tf[tf$NGI_ID %in% metadata$NGI_ID, ] |>
  merge(metadata[, c("NGI_ID","Site","Period","Material")], by="NGI_ID") |>
  arrange(desc(unique_kmers))

# Stampa i top con soglia k-mers >= 1000
head(tf_filt[tf_filt$unique_kmers >= 1000, ], 10)