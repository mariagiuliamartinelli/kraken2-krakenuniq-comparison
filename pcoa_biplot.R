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
make_lookup <- function(meta, col) setNames(meta[[col]], meta$NGI_ID)
get_labels  <- function(ids, lookup) { lab <- lookup[ids]; lab[is.na(lab)] <- "Unknown"; lab }
make_palette <- function(l1, l2) {
  cats <- sort(unique(c(l1, l2)))
  setNames(rainbow(length(cats)), cats)
}

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

dist_k2  <- vegdist(k2_norm, method="bray")
dist_ku  <- vegdist(ku_norm, method="bray")
pcoa_k2  <- cmdscale(dist_k2, k=2, eig=TRUE)
pcoa_ku  <- cmdscale(dist_ku, k=2, eig=TRUE)

var_exp_k2 <- round(pcoa_k2$eig / sum(pcoa_k2$eig) * 100, 1)
var_exp_ku <- round(pcoa_ku$eig / sum(pcoa_ku$eig) * 100, 1)

all_x  <- c(pcoa_k2$points[,1], pcoa_ku$points[,1])
all_y  <- c(pcoa_k2$points[,2], pcoa_ku$points[,2])
x_lims <- range(all_x) * 1.15   # slightly wider than usual to make room for arrow labels
y_lims <- range(all_y) * 1.15

# --------------------------------------------------------------------------
# Biplot helper: draws top-N family arrows on an existing PCoA plot.
#
# Uses vegan::wascores() — the weighted average of sample scores for each
# family, which is the standard ecometric biplot approach.
# Arrow length = correlation of that family with the ordination axes.
# --------------------------------------------------------------------------
add_biplot_arrows <- function(data_norm, pcoa_points, n_arrows = 10,
                               scale = 0.65, col = "gray20", cex_lab = 0.5) {
  # weighted average species scores
  sp      <- wascores(pcoa_points, data_norm)
  vec_len <- sqrt(sp[,1]^2 + sp[,2]^2)
  top_idx <- order(vec_len, decreasing = TRUE)[1:min(n_arrows, nrow(sp))]
  sp_top  <- sp[top_idx, , drop = FALSE]

  # scale: longest arrow reaches `scale` * max sample score
  plot_range <- max(abs(pcoa_points)) * scale
  max_len    <- max(sqrt(sp_top[,1]^2 + sp_top[,2]^2))
  if (max_len == 0) return(invisible(NULL))
  sp_sc <- sp_top * (plot_range / max_len)

  # shorten names: drop trailing 'aceae' -> '.' for legibility
  lab <- gsub("aceae$", ".", rownames(sp_sc))

  arrows(0, 0, sp_sc[,1], sp_sc[,2],
         length = 0.08, col = col, lwd = 1.5)
  text(sp_sc[,1] * 1.18, sp_sc[,2] * 1.18,
       labels = lab, cex = cex_lab, col = col, font = 2, xpd = NA)
}

# --------------------------------------------------------------------------
# Plot functions — biplot variants
# --------------------------------------------------------------------------

plot_pcoa_biplot <- function(points, var_exp, title, data_norm,
                              colors, palette, legend_labels, n_arrows = 10) {
  par(mar = c(5, 4, 4, 13))
  plot(points[,1], points[,2],
       pch = 21, bg = colors, cex = 1.5,
       xlim = x_lims, ylim = y_lims,
       xlab = paste0("PCoA1 (", var_exp[1], "%)"),
       ylab = paste0("PCoA2 (", var_exp[2], "%)"),
       main = title)
  grid()
  abline(h = 0, v = 0, lty = 2, col = "gray")

  # biplot arrows (drawn before legend so they sit under legend text)
  add_biplot_arrows(data_norm, points, n_arrows = n_arrows)

  if (!is.null(palette)) {
    usr <- par("usr")
    par(xpd = NA)
    legend(x = usr[2] + strwidth("MM"), y = usr[4],
           legend = legend_labels, pt.bg = palette[legend_labels],
           pch = 21, pt.cex = 1.0, cex = 0.6, bty = "n", y.intersp = 0.85)
    par(xpd = FALSE)
  }
}

plot_pcoa_combined_biplot <- function(points, var_exp, title, data_norm,
                                       site_labels, mat_labels, n_arrows = 10) {
  par(mar = c(5, 4, 4, 15))
  plot(points[,1], points[,2],
       pch = pal_pch[mat_labels], bg = pal_site[site_labels],
       col = "gray30", cex = 1.5,
       xlim = x_lims, ylim = y_lims,
       xlab = paste0("PCoA1 (", var_exp[1], "%)"),
       ylab = paste0("PCoA2 (", var_exp[2], "%)"),
       main = title)
  grid()
  abline(h = 0, v = 0, lty = 2, col = "gray")

  add_biplot_arrows(data_norm, points, n_arrows = n_arrows)

  usr   <- par("usr")
  x_leg <- usr[2] + strwidth("MM")
  par(xpd = NA)
  leg_site <- legend(x = x_leg, y = usr[4],
                     legend = names(pal_site), pt.bg = pal_site, pch = 21,
                     pt.cex = 1.0, cex = 0.55, bty = "n", y.intersp = 0.85,
                     title = "Site")
  legend(x = x_leg, y = leg_site$rect$top - leg_site$rect$h - strheight("A") * 0.8,
         legend = all_mat_cats, pch = pal_pch[all_mat_cats],
         pt.bg = "gray70", col = "gray30",
         pt.cex = 1.0, cex = 0.65, bty = "n", y.intersp = 0.9,
         title = "Material")
  par(xpd = FALSE)
}

# --------------------------------------------------------------------------
# Generate 10 biplot PCoA plots
# --------------------------------------------------------------------------
cat("Generazione 10 biplot PCoA...\n")

# 1-2: uniform + biplot
plot_pcoa_biplot(pcoa_k2$points, var_exp_k2,
                 paste0("PCoA Biplot - Kraken2 (GTDB) [n=", n, "]"),
                 k2_norm, "blue", NULL, NULL)

plot_pcoa_biplot(pcoa_ku$points, var_exp_ku,
                 paste0("PCoA Biplot - KrakenUniq (MicrobialNT) [n=", n, "]"),
                 ku_norm, "red", NULL, NULL)

# 3-4: Material
plot_pcoa_biplot(pcoa_k2$points, var_exp_k2,
                 paste0("PCoA Biplot - Kraken2 | Material [n=", n, "]"),
                 k2_norm, pal_material[k2_material], pal_material, names(pal_material))

plot_pcoa_biplot(pcoa_ku$points, var_exp_ku,
                 paste0("PCoA Biplot - KrakenUniq | Material [n=", n, "]"),
                 ku_norm, pal_material[ku_material], pal_material, names(pal_material))

# 5-6: Site
plot_pcoa_biplot(pcoa_k2$points, var_exp_k2,
                 paste0("PCoA Biplot - Kraken2 | Site [n=", n, "]"),
                 k2_norm, pal_site[k2_site], pal_site, names(pal_site))

plot_pcoa_biplot(pcoa_ku$points, var_exp_ku,
                 paste0("PCoA Biplot - KrakenUniq | Site [n=", n, "]"),
                 ku_norm, pal_site[ku_site], pal_site, names(pal_site))

# 7-8: Period
plot_pcoa_biplot(pcoa_k2$points, var_exp_k2,
                 paste0("PCoA Biplot - Kraken2 | Period [n=", n, "]"),
                 k2_norm, pal_period[k2_period], pal_period, names(pal_period))

plot_pcoa_biplot(pcoa_ku$points, var_exp_ku,
                 paste0("PCoA Biplot - KrakenUniq | Period [n=", n, "]"),
                 ku_norm, pal_period[ku_period], pal_period, names(pal_period))

# 9-10: Site + Material
plot_pcoa_combined_biplot(pcoa_k2$points, var_exp_k2,
                           paste0("PCoA Biplot - Kraken2 | Site + Material [n=", n, "]"),
                           k2_norm, k2_site, k2_material)

plot_pcoa_combined_biplot(pcoa_ku$points, var_exp_ku,
                           paste0("PCoA Biplot - KrakenUniq | Site + Material [n=", n, "]"),
                           ku_norm, ku_site, ku_material)

cat("Fatto! 10 biplot PCoA generati.\n")
cat("Nota: le frecce mostrano le top-10 famiglie per weighted average score (wascores).\n")
cat("I nomi sono abbreviati: 'aceae' -> '.' per leggibilita'.\n")
