library(readxl)
library(dplyr)

# --------------------------------------------------------------------------
# Pathogen Screen — Yersinia pestis (family: Yersiniaceae)
# --------------------------------------------------------------------------
# Approach: Kraken2 and KrakenUniq report reads at family level.
# We look for the family containing Y. pestis (Yersiniaceae) in both matrices
# and extract per-sample raw read counts.
#
# Other candidate families checked:
#   Tannerella forsythia   -> Tannerellaceae  (periodontal pathogen)
#   Brucella spp.          -> Brucellaceae    (brucellosis, common in pastoralist populations)
#   Mycobacterium tuberculosis -> Mycobacteriaceae (tuberculosis — already in top-20 abundant)
# --------------------------------------------------------------------------

path_k2   <- "kraken2_family_matrix.tsv"
path_ku   <- "krakenuniq_family_matrix.tsv"
path_meta <- "metadataGOG.xlsx"

k2_raw   <- read.table(path_k2,   sep="\t", header=TRUE, row.names=1, check.names=FALSE)
ku_raw   <- read.table(path_ku,   sep="\t", header=TRUE, row.names=1, check.names=FALSE)
metadata <- read_excel(path_meta)

clean_sample_names <- function(names) gsub("_S[0-9]+_L[0-9]+.*", "", names)

valid_ids <- metadata$NGI_ID

k2_filt <- k2_raw[clean_sample_names(rownames(k2_raw)) %in% valid_ids, ]
ku_filt <- ku_raw[clean_sample_names(rownames(ku_raw)) %in% valid_ids, ]

if ("Hominidae" %in% colnames(k2_filt)) k2_filt <- k2_filt[, colnames(k2_filt) != "Hominidae"]
if ("Hominidae" %in% colnames(ku_filt)) ku_filt <- ku_filt[, colnames(ku_filt) != "Hominidae"]

k2_filt <- k2_filt[rowSums(k2_filt) > 0, ]
ku_filt <- ku_filt[rowSums(ku_filt) > 0, ]

# --------------------------------------------------------------------------
# 1. Which pathogen-associated families are present in the matrices?
# --------------------------------------------------------------------------
pathogen_families <- c(
  "Yersiniaceae",       # Yersinia pestis (plague)
  "Tannerellaceae",     # Tannerella forsythia (periodontal)
  "Brucellaceae",       # Brucella spp. (brucellosis)
  "Mycobacteriaceae",   # M. tuberculosis (tuberculosis)
  "Pasteurellaceae",    # Haemophilus, Pasteurella
  "Rickettsiaceae",     # Rickettsia (typhus)
  "Francisellaceae"     # Francisella tularensis (tularemia)
)

cat("=== Pathogen family presence check ===\n")
cat(sprintf("%-22s  %-10s  %-10s\n", "Family", "Kraken2", "KrakenUniq"))
cat(strrep("-", 48), "\n")
for (fam in pathogen_families) {
  in_k2 <- fam %in% colnames(k2_filt)
  in_ku <- fam %in% colnames(ku_filt)
  cat(sprintf("%-22s  %-10s  %-10s\n", fam,
              ifelse(in_k2, "YES", "no"),
              ifelse(in_ku, "YES", "no")))
}
cat("\n")

# --------------------------------------------------------------------------
# 2. Per-sample read counts for each detected pathogen family
# --------------------------------------------------------------------------
extract_pathogen <- function(mat, fam_name, tool_name) {
  if (!fam_name %in% colnames(mat)) {
    cat(tool_name, ": family", fam_name, "not found in matrix.\n")
    return(NULL)
  }
  ids    <- clean_sample_names(rownames(mat))
  counts <- mat[, fam_name]
  df <- data.frame(
    NGI_ID = ids,
    reads  = counts,
    stringsAsFactors = FALSE
  )
  df <- df[df$reads > 0, ]  # keep only samples with at least 1 read
  df <- df[order(df$reads, decreasing = TRUE), ]
  df$tool <- tool_name
  df$family <- fam_name
  return(df)
}

for (fam in pathogen_families) {
  r_k2 <- extract_pathogen(k2_filt, fam, "Kraken2")
  r_ku <- extract_pathogen(ku_filt, fam, "KrakenUniq")

  if (is.null(r_k2) && is.null(r_ku)) next

  cat("=== ", fam, " ===\n", sep="")

  if (!is.null(r_k2) && nrow(r_k2) > 0) {
    cat("Kraken2: ", nrow(r_k2), " samples with reads\n", sep="")
    cat("  Total reads:", sum(r_k2$reads), "\n")
    cat("  Max reads in a single sample:", max(r_k2$reads),
        "(", r_k2$NGI_ID[1], ")\n")
    cat("  Distribution: min =", min(r_k2$reads),
        "| median =", median(r_k2$reads),
        "| max =", max(r_k2$reads), "\n")
    cat("  Top 10 samples:\n")
    print(head(r_k2[, c("NGI_ID","reads")], 10), row.names=FALSE)
  }

  if (!is.null(r_ku) && nrow(r_ku) > 0) {
    cat("KrakenUniq: ", nrow(r_ku), " samples with reads\n", sep="")
    cat("  Total reads:", sum(r_ku$reads), "\n")
    cat("  Max reads in a single sample:", max(r_ku$reads),
        "(", r_ku$NGI_ID[1], ")\n")
    cat("  Distribution: min =", min(r_ku$reads),
        "| median =", median(r_ku$reads),
        "| max =", max(r_ku$reads), "\n")
    cat("  Top 10 samples:\n")
    print(head(r_ku[, c("NGI_ID","reads")], 10), row.names=FALSE)
  }
  cat("\n")
}

# --------------------------------------------------------------------------
# 3. Focus: Yersinia pestis — detailed per-sample table with metadata
# --------------------------------------------------------------------------
cat("=== FOCUS: Yersiniaceae (Yersinia pestis) — detailed ===\n")

build_detail_table <- function(mat, fam, meta) {
  if (!fam %in% colnames(mat)) return(NULL)
  ids <- clean_sample_names(rownames(mat))
  df  <- data.frame(
    NGI_ID = ids,
    reads  = mat[, fam],
    stringsAsFactors = FALSE
  )
  df <- merge(df, meta[, c("NGI_ID","Site","Period","Material")], by="NGI_ID", all.x=TRUE)
  df <- df[order(df$reads, decreasing=TRUE), ]
  return(df)
}

yp_k2 <- build_detail_table(k2_filt, "Yersiniaceae", metadata)
yp_ku <- build_detail_table(ku_filt, "Yersiniaceae", metadata)

if (!is.null(yp_k2)) {
  positives_k2 <- yp_k2[yp_k2$reads > 0, ]
  cat("Kraken2 — samples with Yersiniaceae reads:", nrow(positives_k2), "/", nrow(yp_k2), "\n")
  if (nrow(positives_k2) > 0) {
    print(positives_k2, row.names=FALSE)
  }
  cat("\n")
}

if (!is.null(yp_ku)) {
  positives_ku <- yp_ku[yp_ku$reads > 0, ]
  cat("KrakenUniq — samples with Yersiniaceae reads:", nrow(positives_ku), "/", nrow(yp_ku), "\n")
  if (nrow(positives_ku) > 0) {
    print(positives_ku, row.names=FALSE)
  }
}

# --------------------------------------------------------------------------
# 4. Plot: reads per sample for Yersiniaceae (if present in both)
# --------------------------------------------------------------------------
if (!is.null(yp_k2) && !is.null(yp_ku) &&
    nrow(yp_k2[yp_k2$reads > 0,]) > 0) {

  par(mar=c(8, 5, 4, 2))

  # Merge for side-by-side comparison on shared samples
  both <- merge(
    yp_k2[, c("NGI_ID","reads")],
    yp_ku[, c("NGI_ID","reads")],
    by="NGI_ID", suffixes=c("_k2","_ku")
  )
  both <- both[order(both$reads_k2, decreasing=TRUE), ]
  both_pos <- both[both$reads_k2 > 0 | both$reads_ku > 0, ]

  if (nrow(both_pos) > 0) {
    mat_plot <- t(as.matrix(both_pos[, c("reads_k2","reads_ku")]))
    rownames(mat_plot) <- c("Kraken2","KrakenUniq")
    colnames(mat_plot) <- both_pos$NGI_ID

    barplot(mat_plot,
            beside=TRUE,
            col=c("steelblue","tomato"),
            las=2, cex.names=0.5,
            ylab="Read count (Yersiniaceae)",
            main=paste0("Yersiniaceae reads per sample\n(", nrow(both_pos), " samples with > 0 reads in either tool)"))
    legend("topright", fill=c("steelblue","tomato"),
           legend=c("Kraken2","KrakenUniq"), bty="n")
  }
}

cat("\nDone. Next step: if Yersiniaceae reads are present, validate with\n")
cat("species-level analysis (kraken2 --report at species level) or\n")
cat("ancient DNA damage patterns (mapDamage / PyDamage).\n")
