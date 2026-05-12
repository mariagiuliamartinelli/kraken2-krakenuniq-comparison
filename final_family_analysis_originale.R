library(dplyr)
library(vegan)
library(RColorBrewer)

# Functions
clean_names <- function(n) {
  n <- gsub("\\..*", "", n)
  n <- gsub("_", " ", n)
  return(n)
}

get_top_20 <- function(means) {
  top <- sort(means, decreasing=TRUE)[1:20]
  others <- sum(means) - sum(top)
  return(c(top, "Others" = others))
}

map_abund <- function(means, top_list) {
  v <- sapply(top_list, function(f) if(f %in% names(means)) means[f] else 0)
  c(v, "Others" = 100 - sum(v))
}

# Data Loading
path_k2 <- "kraken2_family_matrix.tsv"
path_ku <- "krakenuniq_family_matrix.tsv"

k2_raw <- read.table(path_k2, sep="\t", header=TRUE, row.names=1, check.names=FALSE)
ku_raw <- read.table(path_ku, sep="\t", header=TRUE, row.names=1, check.names=FALSE)

# Sample alignment: only samples present in BOTH tools
# avoids zero-padding artifacts in PCoA / Procrustes.
common_samples <- sort(intersect(rownames(k2_raw), rownames(ku_raw)))
cat("Sample condivisi tra i due tool:", length(common_samples), "\n")

k2_full <- k2_raw[common_samples, ]
ku_full <- ku_raw[common_samples, ]

# Remove Hominidae present only in KrakenUniq/MicrobialNT
k2_full <- k2_full[, !(colnames(k2_full) %in% c("Hominidae"))]
ku_full <- ku_full[, !(colnames(ku_full) %in% c("Hominidae"))]

# no row with sum 0 after Hominidae
n_zero_k2 <- sum(rowSums(k2_full) == 0)
n_zero_ku <- sum(rowSums(ku_full) == 0)
cat("Sample con somma 0 in Kraken2  :", n_zero_k2, "\n")
cat("Sample con somma 0 in KrakenUniq:", n_zero_ku, "\n")

# removing of 0 sample
keep <- rowSums(k2_full) > 0 & rowSums(ku_full) > 0
if (sum(!keep) > 0) {
  cat("Esclusi", sum(!keep), "sample con riga tutta zero in almeno un tool:",
      paste(rownames(k2_full)[!keep], collapse=", "), "\n")
  k2_full <- k2_full[keep, ]
  ku_full <- ku_full[keep, ]
}
cat("Sample finali in analisi:", nrow(k2_full), "\n")

common_f <- intersect(colnames(k2_full), colnames(ku_full))
cat("Famiglie condivise:", length(common_f), "\n\n")

m_k2 <- colMeans(k2_full)
m_ku <- colMeans(ku_full)

# Statistical Calculations
# Spearman correlation per shared family
rhos <- sapply(common_f, function(f) {
  x <- k2_full[,f]; y <- ku_full[,f]
  if(sum(x)>0 && sum(y)>0) cor(x, y, method="spearman") else NA
})

# PCoA & Procrustes on RELATIVE ABUNDANCES (rowSum = 1)
k2_norm <- sweep(k2_full, 1, rowSums(k2_full), "/")
ku_norm <- sweep(ku_full, 1, rowSums(ku_full), "/")

cat("Range rowSums dopo normalizzazione (Kraken2)   :", range(rowSums(k2_norm)), "\n")
cat("Range rowSums dopo normalizzazione (KrakenUniq):", range(rowSums(ku_norm)), "\n")

d_k2 <- vegdist(k2_norm, method="bray")
d_ku <- vegdist(ku_norm, method="bray")
pcoa_k2 <- cmdscale(d_k2, k=2); pcoa_ku <- cmdscale(d_ku, k=2)
prot <- protest(pcoa_k2, pcoa_ku)
pc1_test <- cor.test(pcoa_ku[,1], pcoa_k2[,1], method="spearman")

# Visualization
par(mfrow=c(1,1), mar=c(5,5,4,2))

# 1. PC1 Comparison
plot(x = pcoa_ku[,1], y = pcoa_k2[,1], xlab = "KrakenUniq PC1", ylab = "Kraken2 PC1",
     main = "PCoA PC1 Agreement", pch = 16, col = rgb(0,0,1,0.5))
abline(lm(pcoa_k2[,1] ~ pcoa_ku[,1]), col="red", lwd=2)
legend("topleft", legend=paste("r =", round(cor(pcoa_ku[,1], pcoa_k2[,1], method="spearman"), 3)), bty="n")

# 2. Mean Family Abundance
plot(x = log10(m_ku[common_f] + 1e-4), y = log10(m_k2[common_f] + 1e-4),
     xlab = "KrakenUniq log10(%)", ylab = "Kraken2 log10(%)",
     main = paste0("Mean Family Abundance"),
     pch = 16, col = rgb(0,0,0,0.4))
abline(0, 1, col="red", lty=2)

# 3. Venn Diagram
plot(0,0, type="n", axes=F, xlab="", ylab="", xlim=c(-1.5, 1.5), ylim=c(-1.5, 1.5), main="Shared Microbial Families")
symbols(x=c(-0.4, 0.4), y=c(0, 0), circles=c(0.9, 0.6), inches=F, add=T, bg=c(rgb(1,0,0,0.2), rgb(0,0,1,0.2)), fg="white")
text(-0.7, 0, length(setdiff(colnames(k2_full), colnames(ku_full)))); text(0.7, 0, length(setdiff(colnames(ku_full), colnames(k2_full)))); text(0, 0, length(common_f))

# 4. Composition Barplot (Normalized)
par(mar=c(5, 4, 4, 12))
m_k2_rel <- (m_k2 / sum(m_k2)) * 100; m_ku_rel <- (m_ku / sum(m_ku)) * 100
comp_k2 <- get_top_20(m_k2_rel); comp_ku <- get_top_20(m_ku_rel)
all_top_f <- unique(c(names(comp_k2), names(comp_ku)))
all_top_f <- all_top_f[all_top_f != "Others"]
barplot_data <- cbind("Kraken2" = map_abund(m_k2_rel, all_top_f), "KrakenUniq" = map_abund(m_ku_rel, all_top_f))
rownames(barplot_data) <- clean_names(rownames(barplot_data))
colors <- c(colorRampPalette(brewer.pal(12, "Paired"))(nrow(barplot_data)-1), "grey90")
barplot(barplot_data, col=colors, border=NA, main="Average Microbial Composition", ylab="Abundance (%)", xlim=c(0, 4), ylim=c(0, 100))
legend(x = 3.0, y = 100, legend = rownames(barplot_data), fill = colors, bty = "n", cex = 0.6, xpd = TRUE)

# 5. Distribution of per-family Spearman correlations
par(mar=c(5, 5, 4, 2))
rhos_clean <- na.omit(rhos)
rho_min <- floor(min(rhos_clean) / 0.05) * 0.05
rho_max <- ceiling(max(rhos_clean) / 0.05) * 0.05
breaks  <- seq(rho_min, rho_max, by = 0.05)
h_max   <- max(hist(rhos_clean, breaks = breaks, plot = FALSE)$counts)

hist(rhos_clean,
     breaks = breaks,
     xlim   = c(rho_min, rho_max),
     ylim   = c(0, h_max + 5),
     main   = "Distribution of per-family Spearman correlations",
     xlab   = "Spearman rho (per family, across samples)",
     ylab   = "Number of families",
     col    = "steelblue", border = "white")

abline(v = mean(rhos_clean),   col = "red",       lwd = 2, lty = 2)
abline(v = median(rhos_clean), col = "darkgreen", lwd = 2, lty = 2)

legend("topright",
       legend = c(sprintf("Mean = %.3f",   mean(rhos_clean)),
                  sprintf("Median = %.3f", median(rhos_clean)),
                  sprintf("N = %d valid families", length(rhos_clean))),
       col = c("red", "darkgreen", NA),
       lty = c(2, 2, NA),
       lwd = c(2, 2, NA),
       bty = "n",
       cex = 0.9,
       inset = 0.02)

# Console Output for Report
cat("\n--- KEY STATISTICS ---\n")
cat("Samples in analysis (intersection):", nrow(k2_full), "\n")
cat("Taxonomy (Families): K2 =", length(colnames(k2_full)), " | KU =", length(colnames(ku_full)), " | Shared =", length(common_f), "\n")
cat("Global Agreement (Procrustes r):", round(prot$t0, 4), " (p =", prot$signif, ")\n")
cat("PCoA PC1 Spearman r =", round(pc1_test$estimate, 4), " (p =", format.pval(pc1_test$p.value, digits=3, eps=2.2e-16), ")\n")
cat("Shared Families Mean   Spearman rho:", round(mean(rhos_clean), 3), "\n")
cat("Shared Families Median Spearman rho:", round(median(rhos_clean), 3), "\n")
cat("Top 5 Consistent Families:\n"); print(round(sort(rhos_clean, decreasing = TRUE)[1:5], 3))

cat("\n--- DISTRIBUTION OF PER-FAMILY SPEARMAN ---\n")
cat("Valid (non-NA) families :", length(rhos_clean), "/", length(rhos), "\n")
cat("Quartiles:\n"); print(round(quantile(rhos_clean, c(0.25, 0.5, 0.75, 0.95)), 3))
cat("Families with rho > 0.5 :", sum(rhos_clean > 0.5), " / ", length(rhos_clean),
    " (", round(100*sum(rhos_clean > 0.5)/length(rhos_clean), 1), "%)\n", sep="")
cat("Families with rho > 0.8 :", sum(rhos_clean > 0.8), " / ", length(rhos_clean),
    " (", round(100*sum(rhos_clean > 0.8)/length(rhos_clean), 1), "%)\n", sep="")
cat("Families with rho < 0   :", sum(rhos_clean < 0), " / ", length(rhos_clean),
    " (", round(100*sum(rhos_clean < 0)/length(rhos_clean), 1), "%)\n", sep="")
