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

## Phase 11 — PCoA Biplot (envfit) (~May 2026)

**File:** `pcoa_biplot.R` (nuovo script; 10 grafici, stesso schema di `pcoa_sito_intersezione.R` + frecce biplot)

**Richiesta di Gonzalo:** *"I guess it could be also nice to add the biplot axis to see what drives the differentiation."*

### Metodo

Si usa `vegan::envfit()` che, per ogni famiglia, fitta la **direzione di massima correlazione** con gli assi PCoA tramite test di permutazione (999 permutazioni). È preferibile a `wascores()` (medie ponderate) perché:
- restituisce r² e p-value → mostra solo frecce statisticamente significative (p ≤ 0.05)
- lunghezza freccia proporzionale a |r| = √(r²) → più interpretabile
- le direzioni sono più distribuite nello spazio

Pre-filtro: solo le top-80 famiglie per abbondanza media vengono passate a `envfit()` (velocità). Filtro angolare aggiuntivo: se due famiglie formano un angolo < 25° tra loro, viene mantenuta solo quella con r² più alto (ridondanza). Default: max 6 frecce per grafico.

### Risultati

**Kraken2** (identici su tutti i 5 plot K2 — la PCoA è la stessa, cambiano solo i colori):
| Famiglia | r² | p |
|---|---|---|
| Streptosporangiaceae | 0.922 | 0.001 |
| Burkholderiaceae | 0.873 | 0.001 |
| Acetobacteraceae | 0.693 | 0.001 |

**KrakenUniq** (identici su tutti i 5 plot KU):
| Famiglia | r² | p |
|---|---|---|
| Burkholderiaceae | 0.968 | 0.001 |
| Rhodobacteraceae | 0.742 | 0.001 |
| Streptosporangiaceae | 0.659 | 0.001 |
| Oxalobacteraceae | 0.492 | 0.001 |

### Interpretazione biologica

Tutti i driver principali sono **batteri del suolo** (Burkholderiaceae, Streptosporangiaceae, Acetobacteraceae, Rhodobacteraceae, Oxalobacteraceae). Questo conferma che il segnale dominante della PCoA è la **contaminazione tafonomica**: il microbiota del sedimento funerario assorbito dalle ossa nel tempo, non differenze biologiche tra individui. La separazione degli ossicle (visibile nei plot Material) riflette la maggiore porosità dell'osso del martello/incudine/staffa rispetto all'osso corticale, che ne facilita l'infiltrazione microbica.

I driver secondari differiscono tra i due tool (Acetobacteraceae in K2, Rhodobacteraceae+Oxalobacteraceae in KU) — coerente con le differenze di database già documentate (GTDB vs MicrobialNT).

---

## Phase 12 — Yersinia pestis screen su Dardel (~May 2026)

**Obiettivo:** rispondere alla richiesta di Gonzalo — *"pick one [pathogen] and check how many reads each sample has for that."* Patogeno scelto: *Yersinia pestis* (famiglia Yersiniaceae, agente della peste bubbonica — plausibile in campioni medievali valencioti).

### Perché su Dardel e non in R locale

Le matrici family-level sul Mac aggregano i dati a livello di famiglia (Yersiniaceae), non di specie. Per distinguere *Y. pestis* da *Y. enterocolitica*, *Y. pseudotuberculosis* ecc., servono i **report grezzi KrakenUniq** a livello di specie — disponibili solo su Dardel in `KrakenUniq_results/KrakenUniq_rep_*.report`.

### Formato report KrakenUniq (colonne rilevanti)

```
col1: % reads nel clade
col2: reads nel clade (inclusi discendenti)
col3: reads dirette alla specie
col4: unique minimizers (k-mers unici) ← metrica chiave per autenticazione
col5: dup (hits medi per minimizer)
col6: coverage stimata
col7: taxID  col8: rank  col9+: nome
```

La colonna **unique k-mers** (col4) è la firma di KrakenUniq: k-mers che mappano *esclusivamente* su *Y. pestis* e non su altre Yersiniaceae. Più è alto, più il segnale è specifico e affidabile.

### Comandi eseguiti su Dardel

```bash
# Verifica: quanti report contengono almeno 1 riga Y. pestis
grep -l "Yersinia pestis" KrakenUniq_results/KrakenUniq_rep_*.report | wc -l
# → 126

# Estrazione per campione: clade_reads, direct_reads, unique_kmers, coverage
for f in KrakenUniq_results/KrakenUniq_rep_*.report; do
    sample=$(basename "$f" | sed 's/KrakenUniq_rep_//' | sed 's/_fastp_merged\.fq\.gz\.report//')
    awk -v s="$sample" '!/^#/ && $8=="species" && /Yersinia pestis/{
        printf "%s\t%s\t%s\t%s\t%s\n", s, $2, $3, $4, $6
    }' "$f"
done | sort -k4 -rn > KrakenUniq_analysis/yersinia_pestis_per_sample.tsv
# → 126 righe salvate
```

File scaricato sul Mac via scp: `CPG/yersinia_pestis_per_sample.tsv`

### Analisi del campione con più segnale grezzo (P35761_1115)

Grezzo (prima del filtro metadata): 354 clade reads, 352 direct reads, 28 unique k-mers. Report completo:
```
Yersiniaceae:                    2979 clade reads   1359 unique k-mers
 └─ Yersinia:                     472 reads           85 unique k-mers
     └─ Y. pseudotuberculosis complex:  373 reads     32 unique k-mers
         └─ Y. pestis:            352 direct reads    28 unique k-mers   coverage 3.96e-05
```
**Conclusione:** segnale presente ma non autentico — coverage 3.96e-05 (≈182 bp coperti su 4.6 Mb), dup=200 (reads concentrate in pochissime posizioni). *Y. pestis* è evolutivamente un clone di *Y. pseudotuberculosis* (97% genoma condiviso): molte reads classificate come *pestis* potrebbero essere *pseudotuberculosis* che cross-mappa sulle regioni uniche. **Questo campione è inoltre escluso dal metadata** → controllo, non campione del progetto.

### Filtro metadata + join in R

```r
yp <- read.table("yersinia_pestis_per_sample.tsv", sep="\t", header=FALSE,
                 col.names=c("sample_raw","clade_reads","direct_reads","unique_kmers","coverage"))
yp$NGI_ID <- gsub("_S[0-9]+_L[0-9]+.*", "", yp$sample_raw)
yp_filt <- yp[yp$NGI_ID %in% metadata$NGI_ID, ] |>
  merge(metadata[, c("NGI_ID","Site","Period","Material")], by="NGI_ID") |>
  arrange(desc(unique_kmers))
```

**Risultato:** 108 campioni del progetto con ≥1 read di *Y. pestis* (su 126 totali; 18 esclusi = controlli).

### Top candidati (campioni del progetto, unique k-mers ≥ 10)

| NGI_ID | Direct reads | Unique k-mers | Period | Site | Material |
|---|---|---|---|---|---|
| P32015_1038 | 1 | 20 | **Medieval** | Cementerio de San Lorenzo | petrous |
| P35503_1031 | 1 | 19 | Roman | Necropolis Occidental | petrous |
| P35503_1050 | 2 | 19 | Visigoth | Necropolis Visigoda, L'Almoina | ossicle |
| P35761_1091 | 52 | 18 | Roman | Legionarios Sertorianos, L'Almoina | other |
| P32015_1089 | 2 | 16 | **Medieval** | Cementerio de San Lorenzo | petrous |
| P35761_1087 | 154 | 16 | Roman | Legionarios Sertorianos, L'Almoina | other |
| P35503_1062 | 5 | 12 | Visigoth | Necropolis Visigoda, L'Almoina | ossicle |
| P35761_1100 | 42 | 12 | Islamic | La Rauda (Otros), L'Almoina | tooth |
| P32015_1020 | 0 | 11 | Islamic | Maqbara Arrabal de l'Alculdia | tooth |
| P35503_1042 | 3 | 11 | Visigoth | Necropolis Visigoda, L'Almoina | ossicle |
| P35503_1049 | 3 | 11 | Visigoth | Necropolis Visigoda, L'Almoina | ossicle |
| P35761_1112 | 4 | 11 | Christian | Masía de Arguinas (Raimundo Morelló) | phalanx |

### Pattern degni di nota

- **P35761_1087 e P35761_1091** provengono entrambi da **Legionarios Sertorianos, L'Almoina** (periodo Romano) e mostrano la combinazione più alta di reads dirette + unique k-mers tra i campioni del progetto. Due individui dello stesso sito romano con segnale Y. pestis — storicamente plausibile (prima pandemia di Giustiniano, VI sec. d.C., lambisce la Spagna).
- **P32015_1038 e P32015_1089** sono gli unici campioni **medievali** con segnale rilevante (Cementerio de San Lorenzo). La II pandemia (XIV sec., Morte Nera) colpisce duramente la Spagna — ma le reads dirette sono molto basse (1–2).
- La maggior parte dei 108 campioni ha ≤5 unique k-mers: **rumore di fondo** da cross-mapping con *Y. pseudotuberculosis* ambientale.

### Interpretazione e limiti

Nessun campione supera 20 unique k-mers. Per l'autenticazione di *Y. pestis* in aDNA le soglie tipiche della letteratura sono:
- ≥200 unique k-mers (molto stringente)
- ≥50 unique k-mers (moderata)
- ≥10 unique k-mers + profilo di danno aDNA (screening preliminare)

I candidati attuali rientrano nella categoria **screening preliminare**: segnale presente, non confermato. Il passo successivo è l'analisi di allineamento diretto al genoma di riferimento di *Y. pestis* CO92 con verifica dei pattern di danno (mapDamage/PyDamage: sostituzioni C→T al 5', G→A al 3').

---

## Phase 13 — Lettura articolo Oskolkov 2025 (biorxiv) e contestualizzazione

**Paper:** Oskolkov et al. 2025, biorxiv — *"Metagenomic Pathogen Detection in Ancient DNA: Benchmarking KrakenUniq and Kraken2 Filtering Strategies"* (titolo approssimativo; paper letto nella sua interezza, non ancora pubblicato su rivista peer-reviewed).

### Principali finding dell'articolo

1. **K-filter (unique k-mers) da solo è ottimale.** Combinare unique k-mers con altri filtri (reads dirette, coverage) non migliora la specificità — anzi può peggiorare la sensibilità. La soglia raccomandata per profiling conservativo è **≥1000 unique k-mers**.

2. **Soglia scalabile con la profondità di sequenziamento.** La relazione empirica è:
   > unique k-mers ≈ 200 × (reads_totali / 100 000)
   
   In pratica: per campioni con ~500 000 reads, ci si aspettano ~1000 k-mers unici se il patogeno è presente in quantità minima rilevabile. Per campioni più piccoli (50 000 reads), la soglia scende a ~100 k-mers.

3. **Kraken2 e KrakenUniq sono comparabili con il filtro ottimale.** La differenza principale è che KrakenUniq fornisce la colonna unique k-mers out of the box; con Kraken2 si può ottenere una stima analoga ma richiede passaggi aggiuntivi.

4. **≥200 reads dirette necessarie per mapDamage.** Profili di danno affidabili (C→T al 5', G→A al 3') richiedono almeno 200 reads direttamente assegnate alla specie target.

5. **≥50 reads** è la soglia minima assoluta per qualsiasi analisi downstream (non sufficiente per mapDamage ma abbastanza per un segnale indicativo).

### Implicazioni per i nostri risultati Y. pestis

| Metrica | Soglia Oskolkov 2025 | Nostro top candidato (P35761_1087) |
|---|---|---|
| Unique k-mers (conservativo) | ≥ 1000 | **16** |
| Unique k-mers (scaling, ~100k reads) | ≥ 200 | **16** |
| Direct reads per mapDamage | ≥ 200 | **154** (vicino, ma sotto soglia) |
| Direct reads minimo assoluto | ≥ 50 | 154 ✓ (P35761_1087) |

**Conclusione definitiva:** tutti i campioni del progetto valenciano con segnale *Y. pestis* si collocano nettamente al di sotto delle soglie validate nell'articolo. Il massimo osservato è 20 unique k-mers (P32015_1038, campione medievale) — 50× sotto la soglia conservativa di 1000. Questo è coerente con **rumore di fondo da cross-mapping** con *Y. pseudotuberculosis*, che condivide il 97% del genoma con *Y. pestis* e che è un batterio ambientale presente nel suolo funerario.

Il risultato non è negativo: lo screen ha permesso di escludere con metodologia validata la presenza di *Y. pestis* in quantità autentica nei campioni analizzati, e di identificare due candidati (P35761_1087, P35761_1091 — romani, L'Almoina) che, se si disponesse di FASTQ grezzi e infrastruttura BWA su Dardel, meriterebbero allineamento al CO92 come passo successivo.

### Note metodologiche per la comunicazione con Gonzalo

- Citare l'articolo Oskolkov 2025 quando si riferiscono le soglie (evita di sembrare soglie arbitrarie).
- Sottolineare che il K-filter è il criterio giusto da usare — non le reads senza unique k-mers.
- Lo stesso approccio (screen → unique k-mers → soglia scaling) è trasferibile a *Tannerella forsythia*.

---

## Phase 14 — Tannerella forsythia screen su Dardel (~May 2026)

**Obiettivo:** secondo patogeno richiesto da Gonzalo — *Tannerella forsythia*, agente della parodontite cronica, atteso soprattutto nei campioni dentali.

### Comandi eseguiti su Dardel

```bash
# Verifica: quanti report contengono almeno 1 riga T. forsythia
grep -l "Tannerella forsythia" KrakenUniq_results/KrakenUniq_rep_*.report | wc -l
# → 112

# Estrazione per campione
for f in KrakenUniq_results/KrakenUniq_rep_*.report; do
    sample=$(basename "$f" | sed 's/KrakenUniq_rep_//' | sed 's/_fastp_merged\.fq\.gz\.report//')
    awk -v s="$sample" '!/^#/ && $8=="species" && /Tannerella forsythia/{
        printf "%s\t%s\t%s\t%s\t%s\n", s, $2, $3, $4, $6
    }' "$f"
done | sort -k4 -rn > KrakenUniq_analysis/tannerella_forsythia_per_sample.tsv
# → 112 righe salvate
```

File scaricato sul Mac: `CPG/tannerella_forsythia_per_sample.tsv`

### Top candidati (prime 20 righe, ordinati per unique k-mers)

| Campione (raw) | Clade reads | Direct reads | Unique k-mers | Coverage |
|---|---|---|---|---|
| P35761_1114 | 8 908 | 7 129 | **151 159** | 0.03566 |
| P32015_1011 | 7 798 | 5 998 | **80 463** | 0.01898 |
| P35761_1113 | 3 064 | 2 439 | **42 707** | 0.01008 |
| P35503_1106 | 1 043 | 715 | **12 657** | 0.002986 |
| P32015_1009 | 766 | 562 | **8 372** | 0.001975 |
| P35761_1102 | 319 | 261 | **3 546** | 0.0008366 |
| P35503_1091 | 106 | 83 | **1 193** | 0.0002815 |
| P35503_1104 | 85 | 61 | 996 | 0.000235 |
| P35503_1056 | 83 | 66 | 910 | 0.0002147 |
| P35503_1094 | 28 | 22 | 331 | 7.809e-05 |
| P35761_1115 | 10 | 7 | 319 | 7.526e-05 |
| P35761_1092 | 11 | 8 | 315 | 7.432e-05 |
| P35761_1097 | 29 | 25 | 261 | 6.158e-05 |
| P35761_1111 | 34 | 28 | 255 | 6.016e-05 |
| P35503_1097 | 18 | 15 | 248 | 5.851e-05 |
| P35761_1087 | 18 | 15 | 241 | 5.686e-05 |
| P35761_1110 | 19 | 14 | 237 | 5.591e-05 |
| P32015_1021 | 35 | 29 | 235 | 5.544e-05 |
| P35761_1105 | 8 | 7 | 133 | 3.138e-05 |
| P35761_2001 | 5 | 5 | 132 | 3.114e-05 |

### Confronto con soglie Oskolkov 2025

| Soglia | Campioni che la superano |
|---|---|
| ≥ 1 000 unique k-mers (conservativa) | **7** (P35761_1114, P32015_1011, P35761_1113, P35503_1106, P32015_1009, P35761_1102, P35503_1091) |
| ≥ 200 direct reads (mapDamage) | **6** (tutti i primi 6 sopra) |
| Segnale rumore di fondo (< 1 000 k-mers) | 105 restanti |

**Contrasto netto con *Y. pestis*:** il top candidato *Y. pestis* aveva 20 unique k-mers; il top *T. forsythia* ne ha **151 159** — 7 500× di più. Il segnale è autentico e biologicamente atteso: *T. forsythia* è un patogeno parodontale che colonizza il solco gengivale e il tartaro dentale, e si conserva molto bene nell'aDNA dentale.

### Filtro metadata + join in R

```r
tf <- read.table("tannerella_forsythia_per_sample.tsv", sep="\t", header=FALSE,
                 col.names=c("sample_raw","clade_reads","direct_reads","unique_kmers","coverage"))
tf$NGI_ID <- gsub("_S[0-9]+_L[0-9]+.*", "", tf$sample_raw)
tf_filt <- tf[tf$NGI_ID %in% metadata$NGI_ID, ] |>
  merge(metadata[, c("NGI_ID","Site","Period","Material")], by="NGI_ID") |>
  arrange(desc(unique_kmers))

head(tf_filt[tf_filt$unique_kmers >= 1000, ], 10)
```

**Nota:** P35761_1114 (151 159 k-mers) e P32015_1011 (80 463 k-mers) — i due campioni con segnale più alto nel TSV grezzo — **non sono nel metadata** e vengono esclusi dal join. Sono probabilmente controlli di laboratorio.

### Campioni del progetto con ≥ 1 000 unique k-mers (soglia Oskolkov conservativa)

| NGI_ID | Direct reads | Unique k-mers | Period | Site | Material |
|---|---|---|---|---|---|
| P35761_1113 | 2 439 | **42 707** | Christian | Masía de Arguinas (Raimundo Morelló) | **phalanx** |
| P35503_1106 | 715 | **12 657** | Islamic | La Rauda (Otros), L'Almoina | **tooth** ✓ |
| P32015_1009 | 562 | **8 372** | Roman | Cementerio San Lorenzo | **petrous** |
| P35761_1102 | 261 | **3 546** | Islamic | La Rauda (Otros), L'Almoina | **tooth** ✓ |
| P35503_1091 | 83 | **1 193** | Islamic | La Rauda (Otros), L'Almoina | **petrous** |

### Interpretazione

**Campioni tooth (attesi):** P35503_1106 e P35761_1102 — entrambi islamici, entrambi da La Rauda, L'Almoina. Confermano l'ipotesi di Gonzalo: *T. forsythia* come patogeno parodontale è rilevabile nei campioni dentali con segnale forte e autentico (>1000 unique k-mers, >200 direct reads → eleggibili per mapDamage).

**Campioni non-dentali (inattesi):**
- **P35761_1113 (falange, Cristiano)** — il campione con più k-mers tra i campioni del progetto (42 707), materiale osseo non dentale. Possibili spiegazioni: (a) contaminazione post-mortem da denti vicini nella sepoltura; (b) batteriemia sistemica in individuo con parodontite grave (*T. forsythia* può entrare in circolo e depositarsi nelle ossa); (c) artefatto tassonomico (sequenze genomiche vicine nel database). Merita approfondimento.
- **P32015_1009 (petroso, Romano)** e **P35503_1091 (petroso, Islamico)** — stesso ragionamento; l'osso petroso è denso e meno soggetto a infiltrazione microbica rispetto alle ossa spugnose, il che rende la contaminazione ambientale meno probabile ma non impossibile.

**Confronto con *Y. pestis*:** il segnale *T. forsythia* è qualitativamente superiore — 7 campioni sopra la soglia conservativa vs. 0 per *Y. pestis*. Per i due campioni tooth islamici (≥200 direct reads), l'autenticazione con mapDamage è tecnicamente fattibile su Dardel.

---

## Phase 15 — Mycobacterium tuberculosis screen su Dardel (~May 2026)

**Razionale:** terzo patogeno aggiunto per completare il quadro — *M. tuberculosis* è documentato in popolazioni antiche europee inclusa la Spagna medievale e romana. Storicamente plausibile, metodologia identica alle fasi 12 e 14.

### Comandi eseguiti su Dardel

```bash
# Verifica presenza
grep -l "Mycobacterium tuberculosis" KrakenUniq_results/KrakenUniq_rep_*.report | wc -l
# → 188 file

# Estrazione per campione
for f in KrakenUniq_results/KrakenUniq_rep_*.report; do
    sample=$(basename "$f" | sed 's/KrakenUniq_rep_//' | sed 's/_fastp_merged\.fq\.gz\.report//')
    awk -v s="$sample" '!/^#/ && $8=="species" && /Mycobacterium tuberculosis/{
        printf "%s\t%s\t%s\t%s\t%s\n", s, $2, $3, $4, $6
    }' "$f"
done | sort -k4 -rn > KrakenUniq_analysis/mycobacterium_tuberculosis_per_sample.tsv
# → 312 righe (> 188 perché il pattern matcha anche sottospecie:
#   M. tuberculosis var. bovis, M. tuberculosis complex ecc. — alcune righe doppie per campione)
```

File scaricato sul Mac: `CPG/mycobacterium_tuberculosis_per_sample.tsv`

### Risultati — top 20 per unique k-mers

| Campione (raw) | Clade reads | Direct reads | Unique k-mers | Coverage |
|---|---|---|---|---|
| P35503_1031 | 64 | 6 | **334** | 5.933e-05 |
| P35503_1044 | 67 | 9 | 281 | 4.992e-05 |
| P32015_1091 | 77 | 2 | 248 | 4.406e-05 |
| P35503_1096 | 45 | 1 | 183 | 3.251e-05 |
| P35761_1087 | 202 | **0** | 162 | 2.878e-05 |
| P35761_1112 | 36 | 0 | 150 | 2.665e-05 |
| … | … | … | … | … |

**Massimo unique k-mers: 334** — nessun campione supera la soglia conservativa di 1 000 (Oskolkov 2025).

### Segnali di cross-mapping

- Molti campioni hanno **0 direct reads** pur avendo clade reads elevati (es. P35761_1087: 202 clade reads, 0 direct reads). Questo è il pattern classico di cross-mapping: le reads vengono assegnate a discendenti di *M. tuberculosis* nel database (Mycobacteriaceae ambientali del suolo) ma nessuna mappa direttamente alla specie bersaglio.
- Il numero di file con match (188) è il più alto tra i tre patogeni screenati — coerente con l'ubiquità delle Mycobacteriaceae nel suolo funerario.

### Filtro metadata + join in R

```r
mt <- read.table("mycobacterium_tuberculosis_per_sample.tsv", sep="\t", header=FALSE,
                 col.names=c("sample_raw","clade_reads","direct_reads","unique_kmers","coverage"))
mt$NGI_ID <- gsub("_S[0-9]+_L[0-9]+.*", "", mt$sample_raw)
mt <- mt |>
  group_by(NGI_ID) |>
  slice_max(unique_kmers, n=1, with_ties=FALSE) |>
  ungroup()
mt_filt <- mt[mt$NGI_ID %in% metadata$NGI_ID, ] |>
  merge(metadata[, c("NGI_ID","Site","Period","Material")], by="NGI_ID") |>
  arrange(desc(unique_kmers))
head(mt_filt, 10)
```

**Top 10 campioni del progetto per unique k-mers:**

| NGI_ID | Direct reads | Unique k-mers | Period | Site | Material |
|---|---|---|---|---|---|
| P35503_1031 | 6 | **334** | Roman | Necropolis Occidental | petrous |
| P35503_1044 | 9 | 281 | Visigoth | Necropolis Visigoda, L'Almoina | petrous |
| P32015_1091 | 2 | 248 | Islamic | Maqbara Arrabal de l'Alculdia | petrous |
| P35503_1096 | 1 | 183 | Islamic | La Rauda (Otros), L'Almoina | petrous |
| P35761_1087 | **0** | 162 | Roman | Legionarios Sertorianos, L'Almoina | other |
| P35761_1112 | 0 | 150 | Christian | Masía de Arguinas | phalanx |
| P35503_1057 | 0 | 145 | Visigoth | Necropolis Visigoda, L'Almoina | petrous |
| P35761_1104 | 0 | 144 | Islamic | La Rauda (Otros), L'Almoina | tooth |
| P35761_1109 | 0 | 137 | Islamic | La Rauda (Otros), L'Almoina | tooth |
| P35503_1046 | 0 | 129 | Visigoth | Necropolis Visigoda, L'Almoina | petrous |

**Osservazioni:**
- Massimo 334 unique k-mers — nessun campione vicino alla soglia di 1 000
- Molti campioni con **0 direct reads** nonostante clade reads elevati → cross-mapping con Mycobacteriaceae ambientali del suolo
- Materiali e periodi distribuiti casualmente — nessun pattern biologico riconoscibile (confronta con *T. forsythia* dove i top campioni erano concentrati su tooth islamici di La Rauda)

### Riepilogo comparativo dei tre patogeni

| Patogeno | File con match | Max unique k-mers | Soglia Oskolkov | Verdetto |
|---|---|---|---|---|
| *Yersinia pestis* | 126 | 20 | ≥ 1 000 | ✗ Rumore di fondo |
| *Mycobacterium tuberculosis* | 188 | 334 | ≥ 1 000 | ✗ Rumore di fondo |
| *Tannerella forsythia* | 112 | 151 159 (controllo) / 42 707 (progetto) | ≥ 1 000 | ✓ Segnale autentico |

**Conclusione:** nessun segnale autentico di TB rilevabile nel dataset di Valencia a questa profondità di sequenziamento. Il risultato è biologicamente coerente: le Mycobacteriaceae del suolo producono molto cross-mapping ma nessun unique k-mers abbastanza alto da superare la soglia di autenticazione. *T. forsythia* rimane l'unico patogeno con segnale genuino nei campioni del progetto.

---

## Open / next steps

1. **Comunicare risultati a Gonzalo** — sunto pronto (Y. pestis, T. forsythia, M. tuberculosis). Attendere feedback prima di procedere con ulteriori analisi.
2. **Validazione patogeni** — fuori scope attuale. Allineamento al genoma di riferimento + mapDamage richiederebbe istruzioni esplicite dal supervisore.
3. **GitHub push** — rimandato a fine progetto.
   - Allineamento al genoma di riferimento *Y. pestis* CO92 (BWA-MEM su Dardel)
   - Verifica danno aDNA con mapDamage o PyDamage (serve ≥200 direct reads → P35761_1087 è vicino con 154)
   - Verifica coverage sui loci plasmidici (pFra, pPla) per distinguere *pestis* da *pseudotuberculosis*
3. **`P34558_1102_S257_L005`** — flagged: campione dominato da DNA umano, microbiota quasi assente dopo rimozione Hominidae.
4. **GitHub** — repo locale inizializzato (`git init`, primo commit). Push da fare quando il progetto è completato.

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
├── Phase 11 outputs
│   └── pcoa_biplot.R                          (10 plot biplot: stessi di pcoa_sito_intersezione.R + frecce envfit)
├── Phase 12 outputs
│   └── yersinia_pestis_per_sample.tsv         (126 righe: NGI_ID | clade_reads | direct_reads | unique_kmers | coverage)
├── Phase 14 outputs
│   └── tannerella_forsythia_per_sample.tsv    (112 righe: NGI_ID | clade_reads | direct_reads | unique_kmers | coverage)
├── Phase 15 outputs
│   └── mycobacterium_tuberculosis_per_sample.tsv  (312 righe: NGI_ID | clade_reads | direct_reads | unique_kmers | coverage)
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

