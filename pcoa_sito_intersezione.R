library(dplyr)
library(vegan)
library(readxl)

path_k2 <- "kraken2_family_matrix.tsv"
path_ku <- "krakenuniq_family_matrix.tsv"
path_meta <- "metadataGOG.xlsx"

k2_raw <- read.table(path_k2, sep="\t", header=TRUE, row.names=1, check.names=FALSE)
ku_raw <- read.table(path_ku, sep="\t", header=TRUE, row.names=1, check.names=FALSE)
metadata <- read_excel(path_meta)

clean_sample_names <- function(names) {
  gsub("_S[0-9]+_L[0-9]+.*", "", names)
}

valid_ids <- metadata$NGI_ID

# filtro per metadata
k2_cleaned_names <- clean_sample_names(rownames(k2_raw))
k2_keep <- k2_cleaned_names %in% valid_ids
k2_filtered <- k2_raw[k2_keep, ]
k2_removed <- rownames(k2_raw)[!k2_keep]
cat("Kraken2: Campioni iniziali:", nrow(k2_raw), "-> Rimasti dopo filtro excel:", nrow(k2_filtered), "\n")
if (length(k2_removed) > 0) {
  cat("  Campioni rimossi (non nel metadata):\n")
  cat(paste0("    - ", k2_removed, "\n"), sep="")
}

ku_cleaned_names <- clean_sample_names(rownames(ku_raw))
ku_keep <- ku_cleaned_names %in% valid_ids
ku_filtered <- ku_raw[ku_keep, ]
ku_removed <- rownames(ku_raw)[!ku_keep]
cat("KrakenUniq: Campioni iniziali:", nrow(ku_raw), "-> Rimasti dopo filtro excel:", nrow(ku_filtered), "\n")
if (length(ku_removed) > 0) {
  cat("  Campioni rimossi (non nel metadata):\n")
  cat(paste0("    - ", ku_removed, "\n"), sep="")
}

# rimozione Hominidae
if ("Hominidae" %in% colnames(k2_filtered)) k2_filtered <- k2_filtered[, colnames(k2_filtered) != "Hominidae"]
if ("Hominidae" %in% colnames(ku_filtered)) ku_filtered <- ku_filtered[, colnames(ku_filtered) != "Hominidae"]

k2_filtered <- k2_filtered[rowSums(k2_filtered) > 0, ]
ku_filtered <- ku_filtered[rowSums(ku_filtered) > 0, ]

# intersezione: solo campioni presenti in entrambi i tool
common_ids <- intersect(clean_sample_names(rownames(k2_filtered)),
                        clean_sample_names(rownames(ku_filtered)))
k2_filtered <- k2_filtered[clean_sample_names(rownames(k2_filtered)) %in% common_ids, ]
ku_filtered <- ku_filtered[clean_sample_names(rownames(ku_filtered)) %in% common_ids, ]
cat("Campioni in comune (intersezione):", length(common_ids), "\n")

# lookup metadata
make_lookup <- function(meta, col) setNames(meta[[col]], meta$NGI_ID)
get_labels <- function(ids, lookup) { lab <- lookup[ids]; lab[is.na(lab)] <- "Unknown"; lab }
make_palette <- function(labels_k2, labels_ku) {
  cats <- sort(unique(c(labels_k2, labels_ku)))
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

cat("Material:", paste(names(pal_material), collapse=", "), "\n")
cat("Site:", paste(names(pal_site), collapse=", "), "\n")
cat("Period:", paste(names(pal_period), collapse=", "), "\n")

# normalizzazione e PCoA
k2_norm <- sweep(k2_filtered, 1, rowSums(k2_filtered), "/")
ku_norm <- sweep(ku_filtered, 1, rowSums(ku_filtered), "/")

cat("Calcolo distanze Bray-Curtis...\n")
dist_k2 <- vegdist(k2_norm, method="bray")
dist_ku <- vegdist(ku_norm, method="bray")

cat("Esecuzione PCoA...\n")
pcoa_k2 <- cmdscale(dist_k2, k=2, eig=TRUE)
pcoa_ku <- cmdscale(dist_ku, k=2, eig=TRUE)

var_exp_k2 <- round(pcoa_k2$eig / sum(pcoa_k2$eig) * 100, 1)
var_exp_ku <- round(pcoa_ku$eig / sum(pcoa_ku$eig) * 100, 1)

all_x <- c(pcoa_k2$points[,1], pcoa_ku$points[,1])
all_y <- c(pcoa_k2$points[,2], pcoa_ku$points[,2])
x_lims <- range(all_x) * 1.1
y_lims <- range(all_y) * 1.1

n <- length(common_ids)

all_mat_cats <- sort(unique(c(k2_material, ku_material)))
pch_values <- c(21, 22, 23, 24, 25, 21, 22, 23)[seq_along(all_mat_cats)]
pal_pch <- setNames(pch_values, all_mat_cats)

plot_pcoa <- function(points, var_exp, title, colors, palette, legend_labels) {
  par(mar=c(5, 4, 4, 13))
  plot(points[,1], points[,2], pch=21, bg=colors, cex=1.5,
       xlim=x_lims, ylim=y_lims,
       xlab=paste0("PCoA1 (", var_exp[1], "%)"),
       ylab=paste0("PCoA2 (", var_exp[2], "%)"), main=title)
  grid(); abline(h=0, v=0, lty=2, col="gray")
  if (!is.null(palette)) {
    usr <- par("usr")
    par(xpd=NA)
    legend(x = usr[2] + strwidth("MM"),
           y = usr[4],
           legend = legend_labels, pt.bg = palette[legend_labels],
           pch=21, pt.cex=1.0, cex=0.6, bty="n", y.intersp=0.85)
    par(xpd=FALSE)
  }
}

plot_pcoa_combined <- function(points, var_exp, title, site_labels, mat_labels) {
  par(mar=c(5, 4, 4, 15))
  plot(points[,1], points[,2], pch=pal_pch[mat_labels], bg=pal_site[site_labels],
       col="gray30", cex=1.5, xlim=x_lims, ylim=y_lims,
       xlab=paste0("PCoA1 (", var_exp[1], "%)"),
       ylab=paste0("PCoA2 (", var_exp[2], "%)"), main=title)
  grid(); abline(h=0, v=0, lty=2, col="gray")
  usr <- par("usr")
  x_leg <- usr[2] + strwidth("MM")
  par(xpd=NA)
  # Site in alto — legend() restituisce il bounding box usato per posizionare Material
  leg_site <- legend(x=x_leg, y=usr[4],
                     legend=names(pal_site), pt.bg=pal_site, pch=21,
                     pt.cex=1.0, cex=0.55, bty="n", y.intersp=0.85, title="Site")
  # Material posizionato esattamente sotto Site, con piccolo gap
  legend(x=x_leg, y=leg_site$rect$top - leg_site$rect$h - strheight("A") * 0.8,
         legend=all_mat_cats, pch=pal_pch[all_mat_cats],
         pt.bg="gray70", col="gray30",
         pt.cex=1.0, cex=0.65, bty="n", y.intersp=0.9, title="Material")
  par(xpd=FALSE)
}

cat("Generazione 10 grafici...\n")

plot_pcoa(pcoa_k2$points, var_exp_k2, paste0("PCoA - Kraken2 (GTDB) [intersezione, n=", n, "]"),           "blue", NULL, NULL)
plot_pcoa(pcoa_ku$points, var_exp_ku, paste0("PCoA - KrakenUniq (MicrobialNT) [intersezione, n=", n, "]"), "red",  NULL, NULL)

plot_pcoa(pcoa_k2$points, var_exp_k2, paste0("PCoA - Kraken2 | Material [n=", n, "]"),    pal_material[k2_material], pal_material, names(pal_material))
plot_pcoa(pcoa_ku$points, var_exp_ku, paste0("PCoA - KrakenUniq | Material [n=", n, "]"), pal_material[ku_material], pal_material, names(pal_material))

plot_pcoa(pcoa_k2$points, var_exp_k2, paste0("PCoA - Kraken2 | Site [n=", n, "]"),    pal_site[k2_site], pal_site, names(pal_site))
plot_pcoa(pcoa_ku$points, var_exp_ku, paste0("PCoA - KrakenUniq | Site [n=", n, "]"), pal_site[ku_site], pal_site, names(pal_site))

plot_pcoa(pcoa_k2$points, var_exp_k2, paste0("PCoA - Kraken2 | Period [n=", n, "]"),    pal_period[k2_period], pal_period, names(pal_period))
plot_pcoa(pcoa_ku$points, var_exp_ku, paste0("PCoA - KrakenUniq | Period [n=", n, "]"), pal_period[ku_period], pal_period, names(pal_period))

plot_pcoa_combined(pcoa_k2$points, var_exp_k2, paste0("PCoA - Kraken2 | Site + Material [n=", n, "]"),    k2_site, k2_material)
plot_pcoa_combined(pcoa_ku$points, var_exp_ku, paste0("PCoA - KrakenUniq | Site + Material [n=", n, "]"), ku_site, ku_material)

cat("Fatto! 10 grafici (intersezione).\n")
