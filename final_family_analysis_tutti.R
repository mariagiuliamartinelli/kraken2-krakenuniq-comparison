library(dplyr)
library(vegan)
library(RColorBrewer)
library(readxl)

# Functions
clean_names <- function(n) {
  n <- gsub("\\..*", "", n)
  n <- gsub("_", " ", n)
  return(n)
}

clean_sample_names <- function(names) {
  gsub("_S[0-9]+_L[0-9]+.*", "", names)
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
path_meta <- "metadataGOG.xlsx"

k2_raw <- read.table(path_k2, sep="\t", header=TRUE, row.names=1, check.names=FALSE)
ku_raw <- read.table(path_ku, sep="\t", header=TRUE, row.names=1, check.names=FALSE)
metadata <- read_excel(path_meta)

# filtro per metadata (tutti i campioni presenti nell'excel, indipendente per tool)
valid_ids <- metadata$NGI_ID

k2_cleaned <- clean_sample_names(rownames(k2_raw))
k2_full <- k2_raw[k2_cleaned %in% valid_ids, ]
cat("Kraken2: campioni dopo filtro metadata:", nrow(k2_full), "\n")

ku_cleaned <- clean_sample_names(rownames(ku_raw))
ku_full <- ku_raw[ku_cleaned %in% valid_ids, ]
cat("KrakenUniq: campioni dopo filtro metadata:", nrow(ku_full), "\n")

# rimozione Hominidae
k2_full <- k2_full[, !(colnames(k2_full) %in% c("Hominidae"))]
ku_full <- ku_full[, !(colnames(ku_full) %in% c("Hominidae"))]

# rimozione righe con somma 0
k2_full <- k2_full[rowSums(k2_full) > 0, ]
ku_full <- ku_full[rowSums(ku_full) > 0, ]
cat("Dopo rimozione vuoti: Kraken2:", nrow(k2_full), "| KrakenUniq:", nrow(ku_full), "\n")

# famiglie condivise per analisi Spearman
common_f <- intersect(colnames(k2_full), colnames(ku_full))
cat("Famiglie condivise:", length(common_f), "\n\n")

m_k2 <- colMeans(k2_full)
m_ku <- colMeans(ku_full)

# Spearman per famiglia condivisa — richiede stessi campioni in entrambi i tool.
# Usiamo l'intersezione dei campioni solo per questa analisi;
# le PCoA restano calcolate su tutti i campioni di ciascun tool.
shared_samples <- intersect(clean_sample_names(rownames(k2_full)),
                            clean_sample_names(rownames(ku_full)))
k2_shared <- k2_full[clean_sample_names(rownames(k2_full)) %in% shared_samples, ]
ku_shared <- ku_full[clean_sample_names(rownames(ku_full)) %in% shared_samples, ]
cat("Campioni usati per Spearman (intersezione):", length(shared_samples), "\n")

rhos <- sapply(common_f, function(f) {
  x <- k2_shared[, f]; y <- ku_shared[, f]
  if(sum(x) > 0 && sum(y) > 0) cor(x, y, method="spearman") else NA
})

# PCoA & Procrustes sui campioni condivisi (stessa logica dello Spearman)
k2_norm_sh <- sweep(k2_shared, 1, rowSums(k2_shared), "/")
ku_norm_sh <- sweep(ku_shared, 1, rowSums(ku_shared), "/")
d_k2_sh <- vegdist(k2_norm_sh, method="bray")
d_ku_sh <- vegdist(ku_norm_sh, method="bray")
pcoa_k2_sh <- cmdscale(d_k2_sh, k=2)
pcoa_ku_sh <- cmdscale(d_ku_sh, k=2)

# correzione segno arbitrario degli assi PCoA
if (cor(pcoa_ku_sh[,1], pcoa_k2_sh[,1]) < 0) pcoa_ku_sh[,1] <- -pcoa_ku_sh[,1]
if (cor(pcoa_ku_sh[,2], pcoa_k2_sh[,2]) < 0) pcoa_ku_sh[,2] <- -pcoa_ku_sh[,2]

prot     <- protest(pcoa_k2_sh, pcoa_ku_sh)
pc1_test <- cor.test(pcoa_ku_sh[,1], pcoa_k2_sh[,1], method="spearman")

par(mar=c(5,5,4,2))

# 1. PC1 Agreement (campioni condivisi)
plot(x = pcoa_ku_sh[,1], y = pcoa_k2_sh[,1],
     xlab = "KrakenUniq PC1", ylab = "Kraken2 PC1",
     main = paste0("PCoA PC1 Agreement (shared samples, n=", length(shared_samples), ")"),
     pch = 16, col = rgb(0,0,1,0.5))
abline(lm(pcoa_k2_sh[,1] ~ pcoa_ku_sh[,1]), col="red", lwd=2)
legend("topleft", legend=paste("r =", round(cor(pcoa_ku_sh[,1], pcoa_k2_sh[,1], method="spearman"), 3)), bty="n")

# 2. Mean Family Abundance
plot(x=log10(m_ku[common_f] + 1e-4), y=log10(m_k2[common_f] + 1e-4),
     xlab="KrakenUniq log10(%)", ylab="Kraken2 log10(%)",
     main="Mean Family Abundance",
     pch=16, col=rgb(0,0,0,0.4))
abline(0, 1, col="red", lty=2)

# 3. Venn Diagram
plot(0,0, type="n", axes=F, xlab="", ylab="", xlim=c(-1.5,1.5), ylim=c(-1.5,1.5), main="Shared Microbial Families")
symbols(x=c(-0.4,0.4), y=c(0,0), circles=c(0.9,0.6), inches=F, add=T,
        bg=c(rgb(1,0,0,0.2), rgb(0,0,1,0.2)), fg="white")
text(-0.7, 0, length(setdiff(colnames(k2_full), colnames(ku_full))))
text(0.7,  0, length(setdiff(colnames(ku_full), colnames(k2_full))))
text(0,    0, length(common_f))

# 4. Composition Barplot
par(mar=c(5,4,4,12))
m_k2_rel <- (m_k2 / sum(m_k2)) * 100
m_ku_rel <- (m_ku / sum(m_ku)) * 100
comp_k2 <- get_top_20(m_k2_rel); comp_ku <- get_top_20(m_ku_rel)
all_top_f <- unique(c(names(comp_k2), names(comp_ku)))
all_top_f <- all_top_f[all_top_f != "Others"]
barplot_data <- cbind("Kraken2"=map_abund(m_k2_rel, all_top_f), "KrakenUniq"=map_abund(m_ku_rel, all_top_f))
rownames(barplot_data) <- clean_names(rownames(barplot_data))
colors <- c(colorRampPalette(brewer.pal(12,"Paired"))(nrow(barplot_data)-1), "grey90")
barplot(barplot_data, col=colors, border=NA, main="Average Microbial Composition",
        ylab="Abundance (%)", xlim=c(0,4), ylim=c(0,100))
legend(x=3.0, y=100, legend=rownames(barplot_data), fill=colors, bty="n", cex=0.6, xpd=TRUE)

# 5. Histogram Spearman
par(mar=c(5,5,4,2))
rhos_clean <- na.omit(rhos)
rho_min <- floor(min(rhos_clean) / 0.05) * 0.05
rho_max <- ceiling(max(rhos_clean) / 0.05) * 0.05
breaks  <- seq(rho_min, rho_max, by = 0.05)
h_max   <- max(hist(rhos_clean, breaks = breaks, plot = FALSE)$counts)
hist(rhos_clean, breaks=breaks, xlim=c(rho_min, rho_max), ylim=c(0,h_max+5),
     main="Distribution of per-family Spearman correlations",
     xlab="Spearman rho (per family, across samples)", ylab="Number of families",
     col="steelblue", border="white")
abline(v=mean(rhos_clean),   col="red",       lwd=2, lty=2)
abline(v=median(rhos_clean), col="darkgreen", lwd=2, lty=2)
legend("topright",
       legend=c(sprintf("Mean = %.3f",   mean(rhos_clean)),
                sprintf("Median = %.3f", median(rhos_clean)),
                sprintf("N = %d valid families", length(rhos_clean))),
       col=c("red","darkgreen",NA), lty=c(2,2,NA), lwd=c(2,2,NA),
       bty="n", cex=0.9, inset=0.02)

# Console Output
cat("\n--- KEY STATISTICS (tutti i campioni; Procrustes su campioni condivisi) ---\n")
cat("Kraken2 campioni totali:", nrow(k2_full), "| KrakenUniq campioni totali:", nrow(ku_full), "\n")
cat("Campioni condivisi (usati per Procrustes e Spearman):", length(shared_samples), "\n")
cat("Famiglie: K2 =", ncol(k2_full), "| KU =", ncol(ku_full), "| Condivise =", length(common_f), "\n")
cat("Global Agreement (Procrustes r):", round(prot$t0, 4), " (p =", prot$signif, ")\n")
cat("PCoA PC1 Spearman r =", round(pc1_test$estimate, 4), " (p =", format.pval(pc1_test$p.value, digits=3, eps=2.2e-16), ")\n")
cat("Shared Families Mean   Spearman rho:", round(mean(rhos_clean), 3), "\n")
cat("Shared Families Median Spearman rho:", round(median(rhos_clean), 3), "\n")
cat("Top 5 Consistent Families:\n"); print(round(sort(rhos_clean, decreasing=TRUE)[1:5], 3))

cat("\n--- DISTRIBUTION OF PER-FAMILY SPEARMAN ---\n")
cat("Valid (non-NA) families:", length(rhos_clean), "/", length(rhos), "\n")
cat("Quartiles:\n"); print(round(quantile(rhos_clean, c(0.25,0.5,0.75,0.95)), 3))
cat("Families with rho > 0.5:", sum(rhos_clean > 0.5), "/", length(rhos_clean),
    "(", round(100*sum(rhos_clean>0.5)/length(rhos_clean),1), "%)\n")
cat("Families with rho > 0.8:", sum(rhos_clean > 0.8), "/", length(rhos_clean),
    "(", round(100*sum(rhos_clean>0.8)/length(rhos_clean),1), "%)\n")
