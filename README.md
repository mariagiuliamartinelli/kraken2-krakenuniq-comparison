# Kraken2 vs KrakenUniq: Family-level Taxonomic Comparison

Analysis scripts for comparing two metagenomic classification pipelines applied to ancient DNA samples from archaeological necropolises in the Valencia area (Roman, Visigothic, Islamic, Medieval, and Christian periods).

## Project overview

Family-level abundance matrices were produced from Kraken2 (GTDB database) and KrakenUniq (MicrobialNT database) classification reports. The scripts perform:

- **Beta-diversity analysis** (PCoA, Bray-Curtis dissimilarity)  
- **Cross-tool agreement** (Procrustes analysis, Spearman correlations)  
- **Taxonomic composition comparison** (Venn diagram, stacked barplot)  
- **PCoA biplots** (family drivers of differentiation)  
- **Pathogen screening** (Yersiniaceae / *Yersinia pestis*, Mycobacteriaceae, and others)

## Scripts

| Script | Description |
|--------|-------------|
| `pcoa_sito_intersezione.R` | PCoA plots coloured by Site, Period, Material, and Site+Material. Uses the **intersection** of samples present in both tools (n = 211). **Main script for the report.** |
| `pcoa_sito_tutti.R` | Same PCoA plots, but each tool keeps all its own samples independently (Kraken2 ≈ 230, KrakenUniq ≈ 220). |
| `final_family_analysis_tutti.R` | PC1 Agreement, Procrustes, Spearman per-family correlations, Venn diagram, composition barplot. Procrustes uses the shared-sample intersection. |
| `pcoa_biplot.R` | PCoA plots with biplot arrows (top-10 family drivers, via `wascores()`). Intersection samples. |
| `pathogen_screen.R` | Screens both matrices for pathogen-associated families (Yersiniaceae, Mycobacteriaceae, Brucellaceae, etc.) and outputs per-sample read counts. |

## Data (not included — private until publication)

The following files are excluded via `.gitignore`:

- `kraken2_family_matrix.tsv` — Kraken2 family-level count matrix (samples × families)
- `krakenuniq_family_matrix.tsv` — KrakenUniq family-level count matrix
- `metadataGOG.xlsx` — sample metadata (NGI_ID, Site, Period, Material)

Scripts expect these files in the **working directory** (set with `setwd()` in RStudio, or run from the CPG folder).

## Dependencies (R)

```r
install.packages(c("dplyr", "vegan", "readxl", "RColorBrewer"))
```

## Usage

```r
# Set working directory to the CPG folder, then:
source("pcoa_sito_intersezione.R")   # main PCoA plots
source("final_family_analysis_tutti.R")  # agreement statistics
source("pcoa_biplot.R")              # biplot version
source("pathogen_screen.R")          # pathogen read counts
```

## Author

**Mariagiulia Martinelli**  
Sapienza University of Rome | Centre for Palaeogenetics, Stockholm University  
Supervisor: Gonzalo Orteo Garcia
