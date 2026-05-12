library(dplyr)
library(vegan)
library(readxl)

path_k2   <- "kraken2_family_matrix.tsv"
path_ku   <- "krakenuniq_family_matrix.tsv"
path_meta <- "metadataGOG.xlsx"

k2_raw   <- read.table(path_k2,   sep="\t", header=TRUE, row.names=1, check.names=FALSE)
ku_raw   <- read.table(path_ku,   sep="\t", header=TRUE, row.names=1, check.names=FALSE)
metadata <- read_excel(path_meta)

clean_sample_names <- function(names) gsub("_S[0-9]+_L[0-9]+.*", "", names)

valid_ids <- metadata$NGI_ID

# filtro metadata
k2_filtered <- k2_raw[clean_sample_names(rownames(k2_raw)) %in% valid_ids, ]
ku_filtered <- ku_raw[clean_sample_names(rownames(ku_raw)) %in% valid_ids, ]

# rimozione Hominidae e righe vuote
if ("Hominidae" %in% colnames(k2_filtered)) k2_filtered <- k2_filtered[, colnames(k2_filtered) != "Hominidae"]
if ("Hominidae" %in% colnames(ku_filtered)) ku_filtered <- ku_filtered[, colnames(ku_filtered) != "Hominidae"]
k2_filtered <- k2_filtered[rowSums(k2_filtered) > 0, ]
ku_filtered <- ku_filtered[rowSums(ku_filtered) > 0, ]

# intersezione campioni
common_ids  <- intersect(clean_sample_names(rownames(k2_filtered)),
                         clean_sample_names(rownames(ku_filtered)))
k2_filtered <- k2_filtered[clean_sample_names(rownames(k2_filtered)) %in% common_ids, ]
ku_filtered <- ku_filtered[clean_sample_names(rownames(ku_filtered)) %in% common_ids, ]
n <- length(common_ids)
cat("Campioni intersezione:", n, "\n")

# lookup metadata
make_lookup  <- function(meta, col) setNames(meta[[col]], meta$NGI_ID)
get_labels   <- function(ids, lookup) { lab <- lookup[ids]; lab[is.na(lab)] <- "Unknown"; lab }
make_palette <- function(l1, l2) { cats <- sort(unique(c(l1, l2))); setNames(rainbow(length(cats)), cats) }

k2_ids <- clean_sample_names(rownames(k2_filtered))
ku_ids <- clean_sample_names(rownames(ku_filtered))

k2_material <- get_labels(k2_ids, make_lookup(metadata, "Material"))
ku_material <- get_labels(ku_ids, make_lookup(metadata, "Material"))
pal_material <- make_palette(k2_material, ku_material)

k2_site <- get_labels(k2_ids, make_lookup(metadata, "Site"))
ku_site <- get_labels(ku_ids, make_lookup(metadata, "Site"))
pal_site <- make_palette(k2_site, ku_site)

k2_period <- get_labels(k2_ids, make_lookup(metadata, "Period"))
ku_period <- get_labels(ku_ids, make_lookup(metadata, "Period"))
pal_period <- make_palette(k2_period, ku_period)

all_mat_cats <- sort(unique(c(k2_material, ku_material)))
pch_values   <- c(21, 22, 23, 24, 25, 21, 22, 23)[seq_along(all_mat_cats)]
pal_pch      <- setNames(pch_values, all_mat_cats)

# normalizzazione e PCoA
k2_norm <- sweep(k2_filtered, 1, rowSums(k2_filtered), "/")
ku_norm <- sweep(ku_filtered, 1, rowSums(ku_filtered), "/")

dist_k2 <- vegdist(k2_norm, method="bray")
dist_ku <- vegdist(ku_norm, method="bray")
pcoa_k2 <- cmdscale(dist_k2, k=2, eig=TRUE)
pcoa_ku <- cmdscale(dist_ku, k=2, eig=TRUE)

var_exp_k2 <- round(pcoa_k2$eig / sum(pcoa_k2$eig) * 100, 1)
var_exp_ku <- round(pcoa_ku$eig / sum(pcoa_ku$eig) * 100, 1)

all_x  <- c(pcoa_k2$points[,1], pcoa_ku$points[,1])
all_y  <- c(pcoa_k2$points[,2], pcoa_ku$points[,2])
x_lims <- range(all_x) * 1.2
y_lims <- range(all_y) * 1.2

# --------------------------------------------------------------------------
# envfit-based biplot helper
#
# Method: envfit() fits a direction vector for each family such that
# the projected sample scores correlate maximally with that family's
# relative abundance across samples. Arrow length ∝ sqrt(r²) = |r|.
# Only families with permutation p ≤ p_cutoff are shown.
#
# Pre-filtering to top `pre_filter` families by mean abundance keeps
# the permutation test fast (avoids testing thousands of near-zero families).
# --------------------------------------------------------------------------
add_envfit_arrows <- function(data_norm, pcoa_points,
                               n_arrows     = 6,
                               pre_filter   = 80,
                               p_cutoff     = 0.05,
                               permutations = 999,
                               scale        = 0.72,
                               min_angle    = 25,    # min degrees between any two arrows
                               col_arrow    = "gray15",
                               col_lab      = "gray15",
                               cex_lab      = 0.55,
                               lwd_arrow    = 1.6) {

  # 1. Keep top families by mean relative abundance (speeds up envfit)
  mean_ab  <- colMeans(data_norm)
  top_fams <- names(sort(mean_ab, decreasing=TRUE))[1:min(pre_filter, length(mean_ab))]
  data_sub <- data_norm[, top_fams, drop=FALSE]

  # 2. envfit: fits correlation direction for each family onto the ordination
  ef    <- envfit(pcoa_points, data_sub, permutations=permutations, choices=1:2)
  dirs  <- ef$vectors$arrows   # unit-length direction vectors
  r2    <- ef$vectors$r        # r² (0–1)
  pvals <- ef$vectors$pvals

  # 3. Filter by significance, sort by r² (best first)
  sig <- which(pvals <= p_cutoff)
  if (length(sig) == 0) {
    message("No significant families at p <= ", p_cutoff)
    return(invisible(NULL))
  }
  candidates <- sig[order(r2[sig], decreasing=TRUE)]
  # take up to 3x n_arrows candidates to have enough after angle filtering
  candidates <- candidates[1:min(n_arrows * 3, length(candidates))]

  dirs_cand  <- dirs[candidates, , drop=FALSE]
  r2_cand    <- r2[candidates]
  pvals_cand <- pvals[candidates]

  # 4. Angular filter: if two arrows point within min_angle of each other,
  #    drop the one with lower r² (it conveys redundant information).
  #    Uses |cos(angle)| so opposite directions are also treated as similar.
  if (min_angle > 0 && nrow(dirs_cand) > 1) {
    min_cos <- cos(min_angle * pi / 180)
    keep    <- rep(TRUE, nrow(dirs_cand))
    for (i in 2:nrow(dirs_cand)) {
      if (!keep[i]) next
      for (j in 1:(i - 1)) {
        if (!keep[j]) next
        cos_ij <- abs(sum(dirs_cand[i, ] * dirs_cand[j, ]))
        if (cos_ij > min_cos) {
          keep[i] <- FALSE   # i has lower r² (list is sorted), remove it
          break
        }
      }
    }
    dirs_cand  <- dirs_cand[keep, , drop=FALSE]
    r2_cand    <- r2_cand[keep]
    pvals_cand <- pvals_cand[keep]
  }

  # 5. Keep final top N
  n_show     <- min(n_arrows, nrow(dirs_cand))
  dirs_top   <- dirs_cand[1:n_show, , drop=FALSE]
  r2_top     <- r2_cand[1:n_show]
  pvals_top  <- pvals_cand[1:n_show]

  # 6. Scale: arrow length ∝ |r| = sqrt(r²), longest arrow = scale * plot range
  plot_range   <- max(abs(pcoa_points)) * scale
  arrow_scaled <- dirs_top * outer(sqrt(r2_top), c(1, 1)) * plot_range

  # 7. Draw arrows
  arrows(0, 0, arrow_scaled[,1], arrow_scaled[,2],
         length=0.09, col=col_arrow, lwd=lwd_arrow)

  # 8. Smart label placement:
  #    - position 20% beyond arrow tip in the same direction
  #    - horizontal justification: right-of-tip if arrow goes right, left-of-tip if left
  #    - vertical justification: above-tip if arrow goes down, below-tip if up
  fam_names <- rownames(arrow_scaled)
  fam_short <- gsub("aceae$", ".", fam_names)

  for (i in seq_len(nrow(arrow_scaled))) {
    ax <- arrow_scaled[i, 1]
    ay <- arrow_scaled[i, 2]
    lx <- ax * 1.20
    ly <- ay * 1.20
    adj_x <- if (ax >= 0) 0 else 1          # left-align if pointing right, right-align if left
    adj_y <- if (ay >= 0) 0 else 1          # below if pointing up, above if pointing down
    text(lx, ly, labels=fam_short[i],
         cex=cex_lab, col=col_lab, font=2,
         adj=c(adj_x, adj_y), xpd=NA)
  }

  # 9. Console summary
  cat("\nenvfit biplot families shown:\n")
  print(data.frame(
    Family = fam_names,
    r2     = round(r2_top,    3),
    p      = round(pvals_top, 3)
  ), row.names=FALSE)
  invisible(list(arrows=arrow_scaled, r2=r2_top, p=pvals_top))
}

# --------------------------------------------------------------------------
# Plot functions
# --------------------------------------------------------------------------
plot_pcoa_biplot <- function(points, var_exp, title, data_norm,
                              colors, palette, legend_labels,
                              n_arrows=8, ...) {
  par(mar=c(5, 4, 4, 13))
  plot(points[,1], points[,2],
       pch=21, bg=colors, cex=1.5,
       xlim=x_lims, ylim=y_lims,
       xlab=paste0("PCoA1 (", var_exp[1], "%)"),
       ylab=paste0("PCoA2 (", var_exp[2], "%)"),
       main=title)
  grid()
  abline(h=0, v=0, lty=2, col="gray")
  add_envfit_arrows(data_norm, points, n_arrows=n_arrows, ...)
  if (!is.null(palette)) {
    usr <- par("usr")
    par(xpd=NA)
    legend(x=usr[2] + strwidth("MM"), y=usr[4],
           legend=legend_labels, pt.bg=palette[legend_labels],
           pch=21, pt.cex=1.0, cex=0.6, bty="n", y.intersp=0.85)
    par(xpd=FALSE)
  }
}

plot_pcoa_combined_biplot <- function(points, var_exp, title, data_norm,
                                       site_labels, mat_labels, n_arrows=8, ...) {
  par(mar=c(5, 4, 4, 15))
  plot(points[,1], points[,2],
       pch=pal_pch[mat_labels], bg=pal_site[site_labels],
       col="gray30", cex=1.5,
       xlim=x_lims, ylim=y_lims,
       xlab=paste0("PCoA1 (", var_exp[1], "%)"),
       ylab=paste0("PCoA2 (", var_exp[2], "%)"),
       main=title)
  grid()
  abline(h=0, v=0, lty=2, col="gray")
  add_envfit_arrows(data_norm, points, n_arrows=n_arrows, ...)
  usr   <- par("usr")
  x_leg <- usr[2] + strwidth("MM")
  par(xpd=NA)
  leg_site <- legend(x=x_leg, y=usr[4],
                     legend=names(pal_site), pt.bg=pal_site, pch=21,
                     pt.cex=1.0, cex=0.55, bty="n", y.intersp=0.85, title="Site")
  legend(x=x_leg, y=leg_site$rect$top - leg_site$rect$h - strheight("A") * 0.8,
         legend=all_mat_cats, pch=pal_pch[all_mat_cats],
         pt.bg="gray70", col="gray30",
         pt.cex=1.0, cex=0.65, bty="n", y.intersp=0.9, title="Material")
  par(xpd=FALSE)
}

# --------------------------------------------------------------------------
# Generate 10 biplot PCoA plots
# Note: envfit uses 999 permutations — may take ~1-2 min per plot.
#       Set permutations=99 for quick preview, 999 for final output.
# --------------------------------------------------------------------------
cat("\nGenerazione 10 biplot PCoA (envfit, 999 permutations)...\n")
cat("Tempo stimato: ~2-4 minuti totali.\n\n")

# 1-2: uniform
cat("--- Plot 1: Kraken2 uniform ---\n")
plot_pcoa_biplot(pcoa_k2$points, var_exp_k2,
                 paste0("PCoA Biplot - Kraken2 (GTDB) [n=", n, "]"),
                 k2_norm, "steelblue", NULL, NULL)

cat("--- Plot 2: KrakenUniq uniform ---\n")
plot_pcoa_biplot(pcoa_ku$points, var_exp_ku,
                 paste0("PCoA Biplot - KrakenUniq (MicrobialNT) [n=", n, "]"),
                 ku_norm, "tomato", NULL, NULL)

# 3-4: Material
cat("--- Plot 3: Kraken2 | Material ---\n")
plot_pcoa_biplot(pcoa_k2$points, var_exp_k2,
                 paste0("PCoA Biplot - Kraken2 | Material [n=", n, "]"),
                 k2_norm, pal_material[k2_material], pal_material, names(pal_material))

cat("--- Plot 4: KrakenUniq | Material ---\n")
plot_pcoa_biplot(pcoa_ku$points, var_exp_ku,
                 paste0("PCoA Biplot - KrakenUniq | Material [n=", n, "]"),
                 ku_norm, pal_material[ku_material], pal_material, names(pal_material))

# 5-6: Site
cat("--- Plot 5: Kraken2 | Site ---\n")
plot_pcoa_biplot(pcoa_k2$points, var_exp_k2,
                 paste0("PCoA Biplot - Kraken2 | Site [n=", n, "]"),
                 k2_norm, pal_site[k2_site], pal_site, names(pal_site))

cat("--- Plot 6: KrakenUniq | Site ---\n")
plot_pcoa_biplot(pcoa_ku$points, var_exp_ku,
                 paste0("PCoA Biplot - KrakenUniq | Site [n=", n, "]"),
                 ku_norm, pal_site[ku_site], pal_site, names(pal_site))

# 7-8: Period
cat("--- Plot 7: Kraken2 | Period ---\n")
plot_pcoa_biplot(pcoa_k2$points, var_exp_k2,
                 paste0("PCoA Biplot - Kraken2 | Period [n=", n, "]"),
                 k2_norm, pal_period[k2_period], pal_period, names(pal_period))

cat("--- Plot 8: KrakenUniq | Period ---\n")
plot_pcoa_biplot(pcoa_ku$points, var_exp_ku,
                 paste0("PCoA Biplot - KrakenUniq | Period [n=", n, "]"),
                 ku_norm, pal_period[ku_period], pal_period, names(pal_period))

# 9-10: Site + Material
cat("--- Plot 9: Kraken2 | Site + Material ---\n")
plot_pcoa_combined_biplot(pcoa_k2$points, var_exp_k2,
                           paste0("PCoA Biplot - Kraken2 | Site + Material [n=", n, "]"),
                           k2_norm, k2_site, k2_material)

cat("--- Plot 10: KrakenUniq | Site + Material ---\n")
plot_pcoa_combined_biplot(pcoa_ku$points, var_exp_ku,
                           paste0("PCoA Biplot - KrakenUniq | Site + Material [n=", n, "]"),
                           ku_norm, ku_site, ku_material)

cat("\nFatto! 10 biplot PCoA con envfit (p <= 0.05, top-8 per r^2).\n")
