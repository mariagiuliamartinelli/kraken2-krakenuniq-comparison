library(readxl)
library(dplyr)

yp <- read.table("yersinia_pestis_per_sample.tsv", sep="\t", header=FALSE,
                 col.names=c("sample_raw","clade_reads","direct_reads","unique_kmers","coverage"))

metadata <- read_excel("metadataGOG.xlsx")

yp$NGI_ID <- gsub("_S[0-9]+_L[0-9]+.*", "", yp$sample_raw)

yp_filt <- yp[yp$NGI_ID %in% metadata$NGI_ID, ] |>
  merge(metadata[, c("NGI_ID","Site","Period","Material")], by="NGI_ID") |>
  arrange(desc(unique_kmers))

cat("Campioni nel metadata con Y. pestis reads:", nrow(yp_filt), "/", nrow(yp), "\n\n")
print(yp_filt[, c("NGI_ID","direct_reads","unique_kmers","Period","Site","Material")], row.names=FALSE)