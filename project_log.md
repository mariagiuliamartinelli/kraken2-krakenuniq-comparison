# Project Log — Ancient DNA Metagenomic Classification (Kraken2 vs KrakenUniq)

**Author:** Mariagiulia Martinelli  
**Programme:** CPG (Centre for Palaeogenetics) — Erasmus exchange  
**Cluster:** Dardel (PDC, KTH)  
**Local environment:** macOS, R 4.5.2, Pages, Keynote  
**Period covered:** March – May 2026  
**AI assistance:** Claude (Anthropic) — used for confirmation of methodological choices, debugging and corrections on specific steps. All outputs verified by the author against intermediate data.

---

## Project overview

Two-stage metagenomic classification of 250 ancient DNA samples followed by a comparative analysis of two classifiers (Kraken2 + GTDB and KrakenUniq + MicrobialNT) at multiple taxonomic levels. The final output is a methodologically defensible comparison at family level, with a full suite of statistical analyses, two PDF reports, an R script, and a mini-presentation prepared for the supervisor.

**Final deliverables in `CPG/`:**
| File | Content |
|---|---|
| `Metagenomic Classification Results: Kraken2 vs KrakenUniq.pdf` | First report (April) — overall classification, top genera/species, % bacteria |
| `Family-level Abundance Matrix Construction (Kraken2 vs KrakenUniq).pdf` | Second report (April-end) — family-level matrices, statistical comparison, dual-reality conclusions |
| `Kraken2_vs_KrakenUniq_presentation.pptx` | 9-slide deck for the in-person mini-presentation to Gonzalo |
| `kraken2_family_matrix.tsv` | 250 × 5461 (samples × families, %) |
| `krakenuniq_family_matrix.tsv` | 248 × 1663 |
| `final_family_analysis.R` | End-to-end R analysis (cleaning, PCoA, Procrustes, plots) |
| `download_matrices.sh` | scp helper to refresh local copies |
| `family_correlations.tsv`, `pcoa_coordinates.tsv` | intermediate R outputs |
| `Rplot.pdf` … `Rplot04.pdf` | The 5 R-generated plots (PCoA, mean abundance, Venn, stacked barplot, histogram) |
| `metadataGOG.xlsx` | Sample metadata |
| `krakenuniq.sh`, `krakenuniq_batch.sh` | KrakenUniq SLURM submission scripts |
| `training/` | Onboarding material — handbooks, FastQC reports, papers, FileZilla |
| `application/` | Erasmus / programme application documents |

---

## Phase 0 — Onboarding & training (late March – early April)

Setup of access, tools, and basic bioinformatics literacy.

- Read CPG handbook (`training/CPG_Handbook_240115_EE.pdf`)
- NAISS user agreement signed (`training/naiss_user_agreement_2025-08-08_english.pdf`)
- Generated SSH key for Dardel (`training/Id25519.ppk`)
- Installed FileZilla locally for file transfer (`training/FileZilla_3.69.6_macos-arm64.app.tar.bz2`)
- Read Dardel/conda/aMeta connection guide (`training/How to Dardel, connect, conda (and install aMeta).pdf`)
- Studied Bioinformatics 101 material (`training/Bioinformatics101.txt`, `training/IntroBioInformatics101.pdf`)
- Read background literature: `s13059-019-1891-0.pdf` (the Kraken2 paper, Genome Biology 2019) and `s41596-022-00738-y.pdf` (Nature Protocols 2022 ancient DNA pipeline)

---

## Phase 1 — Read pre-processing (~April 9)

QC and adapter trimming on the raw FASTQ files of 250 ancient DNA samples.

- **Raw QC (FastQC)**: `P35761_1001_S1_L001_R1_001_fastqc.html`, `..._R2_001_fastqc.html` in `training/` were the per-sample QC outputs.
- **Adapter trimming (Cutadapt)** via SLURM scripts on Dardel: `Cutadapt_Job.sh`, `Cutadapt_Job2.sh`.
- **Post-trim QC**: `trimmed_R1_fastqc.html`, `trimmed_R2_fastqc.html` confirmed adapters were removed and quality acceptable.
- **Read merging** (paired-end → merged single-end FASTQ): `P35761_1001_merged_fastqc.html` shows the merged sample QC.

Output: 250 `*_fastp_merged.fq.gz` files on Dardel, used as input for both classifiers.

---

## Phase 2 — Kraken2 classification on Dardel (~April 13–18)

Ran Kraken2 against the GTDB database (bacteria + archaea only).

**SLURM job array submission** (`SLURM_Kraken2_Multi`, 251-job array). Each task ran:
```bash
kraken2 --db <GTDB_DB> --paired (or --single for merged) --threads <N> \
        --report Kraken2_rep_<sample>.report --output ... <input.fq.gz>
```
Output: 250 `Kraken2_rep_*.report` files (one job failed — sample with only 2 unclassified reads).

**First inspection of one report** (`training/P35761_1001_kraken_report.txt`):
```bash
head -20 Kraken2_rep_P32015_1001_S28_L005_fastp_merged.fq.gz.report
```

---

## Phase 3 — KrakenUniq classification on Dardel (~April 17–21)

Ran KrakenUniq against the MicrobialNT database (bacteria + archaea + fungi + viruses + eukaryotes), excluding *Homo sapiens* at species level.

**SLURM script (`krakenuniq.sh`)** — single-sample template:
```bash
#!/bin/bash
#SBATCH -A naiss2025-5-616
#SBATCH -p memory
#SBATCH --mem=650GB
#SBATCH -t 15:00:00
#SBATCH -J krakenuniq

ml krakenuniq/1.0.4
MIC_NT="/cfs/klemming/projects/supr/archaeogenetics/databases/DBDIR_KrakenUniq_MicrobialNT"
INPUT="..."
OUTPUT_DIR="..."
mkdir -p "$OUTPUT_DIR"
krakenuniq --preload --db "$MIC_NT" --fastq-input "$INPUT" --threads 256 \
  --output "$OUTPUT_DIR/krakenuniq.sequences" \
  --report-file "$OUTPUT_DIR/krakenuniq.output" \
  --gzip-compressed --only-classified-out \
  &> "$OUTPUT_DIR/krakenuniq.log"
```

**Batch version (`krakenuniq_batch.sh`)** to loop over a list of input FASTQs:
```bash
while IFS= read -r INPUT; do
    BASENAME=$(basename "$INPUT" _merged.dust.rmdup.fastq.gz)
    krakenuniq --preload --db "$MIC_NT" --fastq-input "$INPUT" --threads 256 \
      --output "$OUTPUT_DIR/${BASENAME}.krakenuniq.sequences" \
      --report-file "$OUTPUT_DIR/${BASENAME}.krakenuniq.output" \
      --gzip-compressed --only-classified-out \
      &> "$OUTPUT_DIR/${BASENAME}.krakenuniq.log"
done < "$INPUT_LIST"
```

Output: 250 `KrakenUniq_rep_*.report` files plus 1 anomalous file with empty sample ID.

---

## Phase 4 — First exploratory analysis & report (~April 18–21)

Initial exploration of the Kraken2/KrakenUniq output to summarise classification rates and most abundant taxa.

### Kraken2 — exploratory commands
```bash
# Overall classification rate
grep "unclassified" Kraken2_rep_*.report > summary_unclassified.txt
grep -w "root"      Kraken2_rep_*.report > summary_classified.txt

# Top 5 species per sample
for f in Kraken2_rep_*.report; do
    echo "=== $f ==="
    grep -w "S" $f | sort -rn | head -5
done > top_species.txt

# Top 20 species — most abundant overall
grep -h -w "S" Kraken2_rep_*.report | sort -k2 -rn | head -20 > top20_species.txt

# Top 20 species — sum of reads across all samples
grep -h -w "S" Kraken2_rep_*.report > all_species.txt
awk '{sum[$NF]+=$2} END {for (g in sum) print sum[g], g}' all_species.txt \
    | sort -rn | head -20 > top20_species_sum.txt

# Same logic for genera (rank "G")
grep -h -w "G" Kraken2_rep_*.report > all_genus.txt
awk '{sum[$NF]+=$2} END {for (g in sum) print sum[g], g}' all_genus.txt \
    | sort -rn | head -20 > top20_genus_sum.txt

# Top 20 genera with total reads + average % per sample
awk '{reads[$NF]+=$2; perc[$NF]+=$1; count[$NF]++} \
     END {for (g in reads) print reads[g], perc[g]/count[g], g}' all_genus.txt \
    | sort -rn | head -20 > top20_totalreads_%_.txt

# % Bacteria per sample (rank "D" + name "Bacteria")
grep -w "D" Kraken2_rep_*.report | grep "Bacteria" > bacteria_ranking.txt
sort -t: -k2 -n  bacteria_ranking.txt | head -10   # samples with FEWEST bacteria
sort -t: -k2 -rn bacteria_ranking.txt | head -10   # samples with MOST bacteria
```

### KrakenUniq — exploratory commands (different format: rank in col 8, with literal words)
```bash
grep "unclassified" KrakenUniq_rep_*.report > KrakenUniq_summary_unclassified.txt
grep -w "root"      KrakenUniq_rep_*.report > KrakenUniq_summary_classified.txt

# Top 5 species per sample
for f in KrakenUniq_rep_*.report; do
    echo "=== $f ==="
    grep -w "species" "$f" | sort -k2 -rn | head -5
done > KrakenUniq_top_species.txt

# Top 20 species, total reads + avg %, excluding Homo sapiens
grep -h -w "species" KrakenUniq_rep_*.report > KrakenUniq_all_species.txt
grep -v "Homo sapiens" KrakenUniq_all_species.txt | \
awk '{ name=""; for(i=9;i<=NF;i++) name=name (i>9?" ":"") $i;
       gsub(/^ +/, "", name);
       reads[name]+=$2; perc[name]+=$1; count[name]++ }
     END { for (g in reads) print reads[g], perc[g]/count[g], g }' \
    | sort -rn | head -20 > KrakenUniq_top20_species_noHomo.txt

# Top 20 genera, total reads + avg %, no Homo
grep -h -w "genus" KrakenUniq_rep_*.report > KrakenUniq_all_genus.txt
grep -v "Homo" KrakenUniq_all_genus.txt | \
awk '{reads[$NF]+=$2; perc[$NF]+=$1; count[$NF]++} \
     END { for (g in reads) print reads[g], perc[g]/count[g], g }' \
    | sort -rn | head -20 > KrakenUniq_top20_genus_noHomo.txt

# % Bacteria per sample (rank "superkingdom" + name "Bacteria")
grep -w "Bacteria" KrakenUniq_rep_*.report | grep "superkingdom" \
    > KrakenUniq_bacteria_ranking.txt
```

### First report PDF — *Metagenomic Classification Results: Kraken2 vs KrakenUniq*

Compiled in Pages with the bash command snippets, top-10 genera/species tables for both tools, and a comparison block.

**Key findings of the first report:**
- Both tools identified ***Burkholderia*** as the most abundant genus.
- Bacterial read percentages consistently higher in Kraken2 than in KrakenUniq:
  - Kraken2: 0.84% (P32015_1105) → 51.56% (P35761_1115)
  - KrakenUniq: 0.30% (P32015_1105) → 22.65% (P35761_1115)
- Same samples ranked extreme in both tools.
- Kraken2/GTDB assigned alphanumeric codes at species level (e.g. `sp018375725`); KrakenUniq used standard names (`Burkholderia cepacia complex`).

**Top-10 genera (Kraken2):** Burkholderia, Streptosporangium, Streptomyces, Spirillospora, Novosphingobium, Curvibacter, TA-21, Pseudomonas_E, Cupriavidus, Paraburkholderia.

**Top-10 genera (KrakenUniq):** Burkholderia, Streptomyces, Pseudomonas, Paraburkholderia, Streptosporangium, Cupriavidus, Actinomadura, Bradyrhizobium, Methylobacterium, Mesorhizobium.

This report was shared with Gonzalo, who proposed the next step: shift focus to **family level** and produce correlation/PCoA plots between the two tools.

---

## Phase 5 — Family-level matrix construction on Dardel (~April 25)

Following Gonzalo's instructions: build, for each tool, a per-sample TSV with families ordered by abundance, then assemble a `samples × families` matrix with percentages.

### Step 1 — Locate the report files
```bash
ls /cfs/klemming/.../Kraken2_results/*.report     | wc -l   # 250
ls /cfs/klemming/.../KrakenUniq_results/*.report  | wc -l   # 251 (1 anomalous)
```
Sample IDs identical between the two tools — no remapping needed.

### Step 2 — Single-sample test (Kraken2)
Family rank in column 6, exact match `F` (sub-ranks `F1`, `F2` excluded).
```bash
awk '$6=="F" {print $1"\t"$NF}' \
    Kraken2_rep_P32015_1001_S28_L005_fastp_merged.fq.gz.report \
  | sort -k1 -rn | head -20
```
Top hit: `Burkholderiaceae 0.47%`.

### Step 3 — Loop over all Kraken2 reports
```bash
mkdir -p /cfs/klemming/.../Kraken2_analysis/per_sample
cd /cfs/klemming/.../Kraken2_results
for f in Kraken2_rep_*.report; do
    sample=${f#Kraken2_rep_}
    sample=${sample%_fastp_merged.fq.gz.report}
    awk '$6=="F" {print $1"\t"$NF}' "$f" | sort -k1 -rn \
      > ../Kraken2_analysis/per_sample/${sample}_families.tsv
done
```
Output: 250 per-sample TSVs.

### Step 4 — Single-sample test (KrakenUniq)
KrakenUniq format: header lines start with `#`, rank in column 8 with literal `family`, taxon name may contain spaces.
```bash
awk '!/^#/ && $8=="family" {name=$9; for(i=10;i<=NF;i++) name=name" "$i; print $1"\t"name}' \
    KrakenUniq_rep_P32015_1001_S28_L005_fastp_merged.fq.gz.report \
  | sort -k1 -rn | head -20
```
Top hit: `Hominidae 42.94%` — *Homo sapiens* was excluded only at species level, so the family is still present.

### Step 5 — Loop over all KrakenUniq reports
Skip files where the sample ID does not start with `P<digit>` (catches the anomalous file).
```bash
for f in KrakenUniq_rep_*.report; do
    sample=${f#KrakenUniq_rep_}
    sample=${sample%_fastp_merged.fq.gz.report}
    [[ ! "$sample" =~ ^P[0-9] ]] && continue
    awk '!/^#/ && $8=="family" {name=$9; for(i=10;i<=NF;i++) name=name" "$i; print $1"\t"name}' "$f" \
      | sort -k1 -rn > ../KrakenUniq_analysis/per_sample/${sample}_families.tsv
done
```
Output: 250 per-sample TSVs (1 anomalous file skipped).

### Step 6 — Union of families per tool
```bash
cut -f2 Kraken2_analysis/per_sample/*_families.tsv    | sort -u \
    > Kraken2_analysis/all_families.txt        # 5460
cut -f2 KrakenUniq_analysis/per_sample/*_families.tsv | sort -u \
    > KrakenUniq_analysis/all_families.txt     # 1662
```
Kraken2 has 5460 unique family names — many are GTDB alphanumeric placeholders (`0-14-0-10-38-17`, `01-FULL-45-15b`).
KrakenUniq has 1662 standard names including eukaryotes (`Acanthamoebidae` etc.).

### Step 7 — Build the samples × families matrix (single awk one-liner)
```bash
awk -v famfile=<all_families.txt> '
  BEGIN{ FS="\t"; while((getline l<famfile)>0){F[++n]=l; ix[l]=n}; close(famfile);
         printf "Sample"; for(i=1;i<=n;i++) printf "\t%s",F[i]; print "" }
  FNR==1{ if(cur){ printf "%s",cur; for(i=1;i<=n;i++) printf "\t%s",(i in M?M[i]:"0"); print ""; delete M };
          cur=FILENAME; sub(/.*\//,"",cur); sub(/_families.tsv$/,"",cur) }
  { if($2 in ix) M[ix[$2]]=$1 }
  END{ if(cur){ printf "%s",cur; for(i=1;i<=n;i++) printf "\t%s",(i in M?M[i]:"0"); print "" } }
' <per_sample>/*_families.tsv > <output>.tsv
```

### Quality control — empty per-sample files
```bash
find <per_sample> -name "*_families.tsv" -size 0
```
- Kraken2: 1 empty → `P35761_1083_S83_L001` (only 2 reads, all unclassified — failed sample).
- KrakenUniq: 3 empty → `P35761_1101_S101_L001`, `P34558_1098_S253_L005`, `P34558_1100_S255_L005`.

### Output matrices
| Tool | Database | Rows (samples) | Cols (families) | Excluded |
|---|---|---|---|---|
| Kraken2 | GTDB (bacteria/archaea) | 249 | 5460 | P35761_1083 |
| KrakenUniq | MicrobialNT (incl. eukaryotes) | 247 | 1662 | P35761_1101, P34558_1098, P34558_1100 |

Saved on Dardel as `Kraken2_analysis/kraken2_family_matrix.tsv` and `KrakenUniq_analysis/krakenuniq_family_matrix.tsv`.

---

## Phase 6 — Local export via scp (~April 25)

Transferred matrices from Dardel to local CPG folder using the dedicated SSH key, from the **Mac terminal** (not from inside Dardel):

```bash
scp -i ~/.ssh/id-ed25519-pdc \
    mm1202@dardel.pdc.kth.se:/cfs/.../kraken2_family_matrix.tsv \
    ~/Desktop/Università/Magistrale/CPG/

scp -i ~/.ssh/id-ed25519-pdc \
    mm1202@dardel.pdc.kth.se:/cfs/.../krakenuniq_family_matrix.tsv \
    ~/Desktop/Università/Magistrale/CPG/
```
Wrapped in `download_matrices.sh` for repeated use. (Initial confusion with Pages auto-correcting `mm1202@...` into a mailto link and with running scp from inside Dardel by mistake — both resolved by using the helper script.)

---

## Phase 7 — R analysis (~April 25–27)

**File:** `final_family_analysis.R`. Libraries: `vegan`, `dplyr`, `RColorBrewer`.

### Cleaning steps
1. **Sample intersection** — only the 246 samples present in both tools (instead of the union of 250). Avoids zero-padding artifacts in multivariate analysis.
2. **Hominidae removal** — Kraken2/GTDB has no eukaryotes; Hominidae was 43% of KrakenUniq reads.
3. **Empty-row check** — one further sample (`P34558_1102_S257_L005`) had zero microbial reads in KrakenUniq after Hominidae removal → excluded. **Final: 245 samples.**
4. **Relative-abundance normalization** — each sample rescaled to `rowSum = 1` before Bray-Curtis (community-ecology standard).
5. **Family matching** — only families with the exact same name in both databases → **362 shared families**.

### Statistical analyses
**Per-family Spearman correlations** (245 samples × 362 shared families):
```r
rhos <- sapply(common_f, function(f) {
  x <- k2_full[,f]; y <- ku_full[,f]
  if(sum(x)>0 && sum(y)>0) cor(x, y, method="spearman") else NA
})
```

**Bray-Curtis → PCoA → Procrustes** on relative abundances:
```r
k2_norm  <- sweep(k2_full, 1, rowSums(k2_full), "/")
ku_norm  <- sweep(ku_full, 1, rowSums(ku_full), "/")
d_k2     <- vegdist(k2_norm, method="bray")
d_ku     <- vegdist(ku_norm, method="bray")
pcoa_k2  <- cmdscale(d_k2, k=2)
pcoa_ku  <- cmdscale(d_ku, k=2)
prot     <- protest(pcoa_k2, pcoa_ku)
pc1_test <- cor.test(pcoa_ku[,1], pcoa_k2[,1], method="spearman")
```

### Methodological revision rounds
The R script went through 3 rounds of revision driven by the goal of producing a defensible analysis:
1. **Round 1** — switched from `union` (with zero-padding) to `intersect` of sample sets. Procrustes r remained 0.96 → result robust.
2. **Round 2** — added relative-abundance normalization before Bray-Curtis. Procrustes r dropped to 0.90 → previous 0.96 was inflated by sample-size effects. The 0.90 number is methodologically defensible.
3. **Round 3** — added the histogram of per-family Spearman + summary stats (Q1/median/Q3/P95, % above thresholds) to expose the wide distribution behind the moderate mean. Also added `cor.test()` for the PC1 vs PC1 p-value.

### Key results
| Metric | Value |
|---|---|
| Samples in analysis | **245** |
| Shared families (exact name match) | **362** |
| Valid per-family rho | **201 / 362** (161 nominally shared but no data overlap) |
| Median per-family Spearman rho | **0.28** (mean = 0.31) |
| Families with rho > 0.5 | 31.8% (64 / 201) |
| Families with rho > 0.8 | 9.5% (19 / 201) |
| **PCoA PC1 vs PC1 correlation** | **r = 0.945, p < 0.001** |
| **Procrustes correlation** | **r = 0.9003, p < 0.001** |

**Top-correlated families (rho > 0.93):** Burkholderiaceae (0.99), Gemmatimonadaceae (0.96), Rhodanobacteraceae (0.94), Steroidobacteraceae (0.94), Micromonosporaceae (0.94).

### R-generated plots (in `Rplot.pdf` … `Rplot04.pdf`)
1. PCoA PC1 Agreement scatter (KrakenUniq X vs Kraken2 Y) — r = 0.945
2. Mean Family Abundance log10 scatter — 362 shared families, 1:1 line
3. Venn diagram of shared family names — 5098 / 362 / 1299
4. Stacked barplot of average microbial composition — top 20 families normalized to 100%
5. Histogram of per-family Spearman correlations — mean + median lines, 201 valid families

---

## Phase 8 — Second report (~April 27)

**File:** `Family-level Abundance Matrix Construction (Kraken2 vs KrakenUniq).pdf` (4 pages, written in Pages).

**Sections:**
- Disclosure of AI assistance + tool versions
- Steps 1-7 (bash pipeline) with `# Claude's help` inline markers on the four most complex blocks (Step 3 loop, Step 4 awk, Step 5 loop, Step 7 matrix one-liner)
- Quality control + output matrix table
- Step 8: scp export and R setup
- Step 9: Cleaning decisions in R
- Step 10: Statistical analysis (Spearman per-family + PCoA/Procrustes side by side)
- Step 11: Visualization — 5 plots with captions
- Final Conclusions: the **dual reality**

### Final Conclusions — "Dual reality"
- *Local taxonomic level*: significant nomenclature divergence (GTDB placeholders vs NCBI names). 161 of 362 "shared" families have no data overlap; only ~10% show strong per-family agreement (rho > 0.8).
- *Global ecological level*: highly consistent — PC1 r = 0.945, Procrustes r = 0.90 (p < 0.001).
- For downstream analyses (clustering, PCoA, group comparisons), both tools give the same biological conclusions. For taxon-specific work (e.g. pathogens), tool choice matters.

### Gonzalo's feedback
> *"Very nice. I guess there is quite a bit to unpack there but well done. I say, next time you can come you explain it to me as if it's a mock mini-presentation. Perhaps there is still a lot of complexity at the Family level, maybe I should have told you to do it at Phylum level but it's good. Next step is think of some pathogen. I was thinking* Yersinia pestis *because it's possible some Medieval samples carry it. But* Tannerella forsythia *too, that one we will find in the teeth with some luck."*

---

## Phase 9 — Mini-presentation (~May 1)

**File:** `Kraken2_vs_KrakenUniq_presentation.pptx` (9 slides, 16:9, palette Forest & Moss, with speaker notes).

| # | Slide | Visual |
|---|---|---|
| 1 | Title | Dark forest background |
| 2 | The question | "250" big stat + question prompt |
| 3 | Building the comparison | 3-step pipeline + 2 matrix tiles |
| 4 | Cleaning decisions | 2×2 grid (4 decisions) + outcome strip |
| 5 | Two questions, two analyses | 2-column comparison |
| 6 | Result A: Global agreement | Big stat r = 0.945 + Plot 1 |
| 7 | Result B: Per-family agreement | Stat tiles + Plot 5 + top-5 families |
| 8 | Why: a nomenclature problem | Plot 3 (Venn) + 3 numbers |
| 9 | The dual reality + closing question | Dark forest, 2 columns, accent strip |

Speaker notes (one paragraph per slide) provide the spoken script. Length target: 5-7 minutes.

---

## Phase 10 — PCoA con rimozione campioni non nel metadata (~May 2026)

**Files:** vedi File map (sezione dedicata Phase 10)

Obiettivo: produrre PCoA separati per Kraken2 e KrakenUniq filtrando i campioni sulla base di `metadataGOG.xlsx` e visualizzarli nella console di R (pannello Plots di RStudio) senza salvataggio automatico su disco.

### Cos'è il filtro metadata e perché serve

Le matrici TSV (`kraken2_family_matrix.tsv`, `krakenuniq_family_matrix.tsv`) contengono tutti i campioni classificati dal tool, compresi eventuali campioni di controllo, campioni mal classificati o campioni non pertinenti al progetto. Il file `metadataGOG.xlsx` è la lista ufficiale dei campioni del progetto GOG: contiene solo quelli che hanno superato i controlli di qualità e che sono biologicamente rilevanti per le analisi.

**Il filtro metadata fa esattamente questo: mantiene solo i campioni il cui `NGI_ID` è presente in `metadataGOG.xlsx`, rimuovendo tutti gli altri.**

Passaggi tecnici:
1. I nomi riga delle matrici contengono il suffisso `_S[0-9]+_L[0-9]+...` (aggiunto dalla pipeline NGI); viene rimosso con `gsub` per ottenere l'`NGI_ID` netto.
2. Solo le righe il cui `NGI_ID` corrisponde a un valore in `metadata$NGI_ID` vengono mantenute.
3. I campioni esclusi vengono stampati esplicitamente in console con il loro nome completo per tracciabilità.
4. Dopo la rimozione di Hominidae e dei campioni con somma = 0, i campioni rimanenti vengono normalizzati ad abbondanza relativa e si calcola la distanza Bray-Curtis → PCoA con `cmdscale`.

### Versioni `_tutti` vs `_intersezione`

Per ogni script esistono due varianti:

| Variante | Campioni usati | Procrustes |
|---|---|---|
| `_tutti` | Ogni tool mantiene i propri campioni post-filtro (possono differire tra K2 e KU) | Non eseguito (set diversi) |
| `_intersezione` | Solo campioni presenti **in entrambi** i tool post-filtro | Possibile (set identici) |

La variante `_originale` (`final_family_analysis_originale.R`) non applica nessun filtro metadata: usa l'intersezione grezza dei nomi riga così come escono dalle matrici (~246 campioni, r = 0.945).

La variante `_intersezione` con filtro metadata usa ~212 campioni e dà r = 0.953. La differenza (0.945 → 0.953) è attesa: si sta confrontando un sottoinsieme leggermente diverso di campioni, non c'è nessun problema metodologico.

> **Nota sul segno degli assi PCoA:** `cmdscale()` restituisce assi con segno arbitrario. Nella versione `_intersezione`, PC1 di KrakenUniq risultava anticorrelato rispetto a Kraken2 (r = −0.953 invece di +0.953). Lo script include un controllo automatico che inverte il segno se necessario:
> ```r
> if (cor(pcoa_ku[,1], pcoa_k2[,1]) < 0) pcoa_ku[,1] <- -pcoa_ku[,1]
> if (cor(pcoa_ku[,2], pcoa_k2[,2]) < 0) pcoa_ku[,2] <- -pcoa_ku[,2]
> ```
> Questo non altera la struttura della PCoA, solo l'orientamento visivo.

### Colorazione e visualizzazione
- Costruita una lookup `NGI_ID → colonna` dalla colonna corrispondente di `metadataGOG.xlsx` (Material, Site, Period).
- La palette è generata automaticamente con `rainbow()` sulle categorie presenti; campioni senza corrispondenza ricevono l'etichetta `"Unknown"`.
- La legenda è posizionata dentro il grafico (`"topright"`, `cex=0.55`) per non uscire dal pannello Plots.
- Le categorie trovate vengono stampate in console all'avvio per verifica immediata.

### Output — grafici nella cronologia Plots di RStudio

**Script `pcoa_rimozione_differenza_campioni_*.R`** — 4 grafici:
1. Kraken2 — colore uniforme (blu)
2. KrakenUniq — colore uniforme (rosso)
3. Kraken2 — colorato per `Material`
4. KrakenUniq — colorato per `Material`

**Script `pcoa_sito_*.R`** — 10 grafici:
1. Kraken2 — colore uniforme
2. KrakenUniq — colore uniforme
3. Kraken2 — per `Material`
4. KrakenUniq — per `Material`
5. Kraken2 — per `Site`
6. KrakenUniq — per `Site`
7. Kraken2 — per `Period`
8. KrakenUniq — per `Period`
9. Kraken2 — `Site` (colore) + `Material` (forma punto)
10. KrakenUniq — `Site` (colore) + `Material` (forma punto)

Assi comuni tra tutti i grafici per comparabilità diretta. Nessun file salvato automaticamente: l'export è manuale via "Export" in RStudio.

### Pulizia script (successiva)

Dopo una revisione, 4 script sono stati eliminati perché il loro contenuto era già interamente presente in `pcoa_sito_*.R`:

| Eliminato | Sostituito da |
|---|---|
| `pcoa_separate_intersezione.R` (2 plot) | plot 1-2 di `pcoa_sito_intersezione.R` |
| `pcoa_rimozione_differenza_campioni_intersezione.R` (4 plot) | plot 1-4 di `pcoa_sito_intersezione.R` |
| `pcoa_separate_tutti.R` (2 plot) | plot 1-2 di `pcoa_sito_tutti.R` |
| `pcoa_rimozione_differenza_campioni_tutti.R` (4 plot) | plot 1-4 di `pcoa_sito_tutti.R` |

Da `final_family_analysis_tutti.R` sono stati rimossi i 2 plot PCoA uniformi (già coperti da `pcoa_sito_tutti.R`); rimangono i 4 plot statistici unici (Mean Abundance, Venn, Barplot, Spearman).

**Struttura finale: 5 script** organizzati in `intersezione/` e `tutti/` + `final_family_analysis_originale.R` nella root.

---

## Open / next steps

Discussed with Gonzalo, not yet started:

1. **Pathogen analysis** — search the existing KrakenUniq reports for:
   - *Yersinia pestis* (plague, plausible in Medieval samples)
   - *Tannerella forsythia* (oral pathogen, expected in dental calculus / teeth samples)
   Use unique k-mer counts (KrakenUniq's strength) for low-abundance authentication.
2. **Open question to Gonzalo** — pathogens immediately, or first aggregate to Phylum to check whether the per-family disagreement was just a resolution effect.
3. **`P34558_1102_S257_L005`** — flagged: sample dominated by human DNA with negligible microbial biomass. Could be biologically interesting in the dental subset.

---

## Tools / software stack

| Layer | Where it runs | Tool |
|---|---|---|
| Read pre-processing (adapter trimming, QC, merging) | Dardel | Cutadapt, FastQC |
| Read classification | Dardel | Kraken2 + GTDB; KrakenUniq + MicrobialNT (via SLURM) |
| Data wrangling (extraction, matrices) | Dardel | bash + awk |
| File transfer | Mac terminal | scp |
| Statistical analysis | Mac | R 4.5.2 + vegan 2.7-3 + dplyr + RColorBrewer |
| Report writing | Mac | Pages → PDF |
| Presentation | Mac | Keynote (.pptx) |

---

## File map (CPG/)

```
CPG/
├── Phase 4 outputs
│   ├── Metagenomic Classification Results: Kraken2 vs KrakenUniq.pdf   (first report)
│   └── Metagenomic Classification Results: Kraken2 vs KrakenUniq.pages
├── Phase 5 outputs
│   ├── kraken2_family_matrix.tsv          (4.5 MB, 250×5461)
│   └── krakenuniq_family_matrix.tsv       (1.6 MB, 248×1663)
├── Phase 6 helpers
│   └── download_matrices.sh
├── Phase 7 outputs
│   ├── final_family_analysis.R
│   ├── family_correlations.tsv
│   ├── pcoa_coordinates.tsv
│   └── Rplot.pdf, Rplot01.pdf, Rplot02.pdf, Rplot03.pdf, Rplot04.pdf
├── Phase 8 outputs
│   ├── Family-level Abundance Matrix Construction (Kraken2 vs KrakenUniq).pdf
│   ├── Family-level Abundance Matrix Construction (Kraken2 vs KrakenUniq).pages
│   └── family_matrices_report.md          (early markdown draft)
├── Phase 9 outputs
│   └── Kraken2_vs_KrakenUniq_presentation.pptx
├── Phase 10 outputs
│   ├── final_family_analysis_originale.R      (no filtro metadata, ~246 campioni, r=0.945)
│   ├── intersezione/
│   │   ├── final_family_analysis_intersezione.R  (5 plot: PC1 agreement, Mean Abundance, Venn, Barplot, Spearman)
│   │   └── pcoa_sito_intersezione.R              (10 plot: uniform, Material, Site, Period, Site+Material)
│   └── tutti/
│       ├── final_family_analysis_tutti.R         (4 plot: Mean Abundance, Venn, Barplot, Spearman; no Procrustes)
│       └── pcoa_sito_tutti.R                     (10 plot: uniform, Material, Site, Period, Site+Material)
├── Phase 3 SLURM scripts
│   ├── krakenuniq.sh
│   └── krakenuniq_batch.sh
├── Sample metadata
│   └── metadataGOG.xlsx
├── This log
│   └── project_log.md
├── application/                           (Erasmus / programme docs)
└── training/                              (onboarding material, FastQC, papers, FileZilla)
```

