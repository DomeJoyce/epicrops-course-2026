# Day 2 – Practical Exercise
## DNA Methylation Analysis: Decoding the Plant Methylome
### The Epi-Code – Florence Training School 2026
#### Instructors: Domenico Giosa & Letterio Giuffrè | University of Messina

---

## Learning Objectives

By the end of Day 2 you will be able to:

1. Assess quality of bisulfite-converted sequencing reads.
2. Align WGBS data to a reference genome using Bismark's bisulfite-aware aligner.
3. Extract per-cytosine methylation calls across CpG, CHG, and CHH contexts.
4. Identify differentially methylated cytosines (DMCs) and regions (DMRs) with methylKit.
5. Interpret context-specific methylation profiles in the context of plant stress epigenetics.
6. Connect DNA methylation patterns to the RNA-directed DNA methylation (RdDM) pathway
   and to the lncRNA transcriptional layer studied in Day 1.

---

## Biological Background

**DNA methylation in plants** occurs at three cytosine sequence contexts:

| Context | Sequence | Maintenance mechanism          | Biological role                          |
|---------|----------|-------------------------------|------------------------------------------|
| CpG     | CG       | MET1 (plant DNMT1 homolog)    | Gene body methylation, TE silencing      |
| CHG     | CHG      | CMT3                          | TE silencing, repetitive elements        |
| CHH     | CHH      | DRM2 (RdDM pathway) / CMT2   | RdDM, stress-responsive de novo meth.   |

**Long non-coding RNAs and DNA methylation are directly linked**: many lncRNAs act *in
cis*, recruiting chromatin- and DNA-methylation machinery to neighbouring loci, or are
themselves transcribed antisense to a differentially methylated region — exactly the
cis-regulatory relationship you tested for Day 1's DE lncRNAs. Small-RNA-guided
methylation (**RdDM**) is a second, well-established route to the same molecular outcome:
24 nt hc-siRNAs guide the methyltransferase DRM2 to silence transposable elements via CHH
methylation. Part 6 today tests the lncRNA↔methylation link directly, using the real DE
lncRNA coordinates produced on Day 1.

**Whole-Genome Bisulfite Sequencing (WGBS)** is the gold-standard method for
single-base resolution methylome profiling. Sodium bisulfite converts unmethylated
cytosines to uracil (read as thymine), while methylated cytosines remain unchanged.

Today we analyse **Arabidopsis thaliana Chr4** WGBS data (NCBI BioProject
**PRJNA700573**) comparing wild-type Col-0 control plants against **salt-stressed**
plants — **3 biological replicates per condition**. Salt stress triggers the same
ABA-mediated osmotic-stress signalling explored in Day 1, so the two days examine
two regulatory layers (lncRNA-based transcriptional regulation and DNA methylation) of
one stress response. The same study also profiled the *nrpd1* (Pol IV) RdDM mutant — a
natural reference point when we discuss the 24 nt hc-siRNA → CHH-methylation link.

---

## Part 0 – Start the Day 2 environment

Unlike Day 1, the Day 2 image has no browser UI bundled — you type commands in a plain
terminal attached to the container. You have two equivalent ways to do that; pick whichever
you're more comfortable with, both run the exact same commands inside the exact same
container:

**Option A — plain terminal (default, works everywhere: bash, zsh, or Windows PowerShell):**

```bash
docker compose run --rm course
```

This drops you straight into a shell inside the container. Every command in this
walkthrough goes here.

**Option B — JupyterLab-style browser terminal (optional, same experience as Day 1):**

```bash
docker compose up jupyter
```

Then open **http://localhost:8888/?token=epicode2026** in your browser, and use
**File → New → Terminal** — a terminal opens inside the container, and you type the exact
same commands there as in Option A.

Either way, results write to `results/day2_methylation/` **on your laptop** automatically
(bind-mounted) — unlike Day 1's self-contained image, there is no manual export step here.

> **Laptop note.** Commands below use **4** CPU threads. On a small laptop, change each `4`
> to `2`.

---

## Part 1 – Quality Control of Bisulfite-Seq Data

### 1.1 Start the Docker container

(Already done in Part 0 — skip if you're continuing in the same terminal.)

```bash
docker compose run --rm course
```

Verify tools are available:

```bash
bismark --version

trim_galore --version

samtools --version

R --version
```

### 1.2 Inspect raw FASTQ files

```bash
ls -lh /course/data/raw/wgbs/
```

WGBS reads are paired-end. Look at R1 and R2 of one sample:

```bash
seqtk seq -A /course/data/raw/wgbs/SRR13650161_1.fastq.gz | head -8

seqtk seq -A /course/data/raw/wgbs/SRR13650161_2.fastq.gz | head -8
```

> **Question 1:** Notice the extremely high T content compared to standard RNA-seq.
> Why does bisulfite conversion dramatically increase T frequency? What base
> composition pattern do you expect in the reverse complement strand (R2)?

### 1.3 FastQC quality control

First create the output folders once (a single plain command — nothing is hidden):

```bash
mkdir -p /course/results/day2_methylation/01_qc /course/results/day2_methylation/02_trimmed /course/results/day2_methylation/03_bismark_aligned /course/results/day2_methylation/04_methylation /course/results/day2_methylation/05_diffmeth
```

Run FastQC on all 12 FASTQ files (R1 + R2 of the 6 samples), then summarise with MultiQC:

```bash
fastqc --threads 4 --outdir /course/results/day2_methylation/01_qc /course/data/raw/wgbs/SRR13650161_1.fastq.gz /course/data/raw/wgbs/SRR13650161_2.fastq.gz /course/data/raw/wgbs/SRR13650163_1.fastq.gz /course/data/raw/wgbs/SRR13650163_2.fastq.gz /course/data/raw/wgbs/SRR13650165_1.fastq.gz /course/data/raw/wgbs/SRR13650165_2.fastq.gz /course/data/raw/wgbs/SRR13650194_1.fastq.gz /course/data/raw/wgbs/SRR13650194_2.fastq.gz /course/data/raw/wgbs/SRR13650196_1.fastq.gz /course/data/raw/wgbs/SRR13650196_2.fastq.gz /course/data/raw/wgbs/SRR13650198_1.fastq.gz /course/data/raw/wgbs/SRR13650198_2.fastq.gz

multiqc /course/results/day2_methylation/01_qc --outdir /course/results/day2_methylation/01_qc/multiqc --filename Day2_raw_QC
```

Open `./results/day2_methylation/01_qc/multiqc/Day2_raw_QC.html`.

> **Question 2:** The "Per base sequence content" module will almost certainly fail.
> Is this a problem? Why does bisulfite sequencing always produce non-uniform
> base composition at the library level?

---

## Part 2 – Adapter Trimming for WGBS

WGBS data requires additional trimming considerations:
- Paired-end reads require coordinated trimming of R1 and R2
- Additional end clipping (6 bp) removes residual random hexamer priming biases
- The `--paired` flag enforces concordant trimming

### 2.1 Trim each sample

Run Trim Galore **once per sample** — the same command six times, changing only the
sample accession (three control: `SRR13650161/163/165`, three salt: `SRR13650194/196/198`).

```bash
# --- SRR13650161 (control) ---
trim_galore --paired --cores 4 --quality 20 --stringency 5 --dont_gzip --clip_r1 6 --clip_r2 6 --three_prime_clip_R1 6 --three_prime_clip_R2 6 --output_dir /course/results/day2_methylation/02_trimmed /course/data/raw/wgbs/SRR13650161_1.fastq.gz /course/data/raw/wgbs/SRR13650161_2.fastq.gz
```

```bash
# --- SRR13650163 (control) ---
trim_galore --paired --cores 4 --quality 20 --stringency 5 --dont_gzip --clip_r1 6 --clip_r2 6 --three_prime_clip_R1 6 --three_prime_clip_R2 6 --output_dir /course/results/day2_methylation/02_trimmed /course/data/raw/wgbs/SRR13650163_1.fastq.gz /course/data/raw/wgbs/SRR13650163_2.fastq.gz
```

```bash
# --- SRR13650165 (control) ---
trim_galore --paired --cores 4 --quality 20 --stringency 5 --dont_gzip --clip_r1 6 --clip_r2 6 --three_prime_clip_R1 6 --three_prime_clip_R2 6 --output_dir /course/results/day2_methylation/02_trimmed /course/data/raw/wgbs/SRR13650165_1.fastq.gz /course/data/raw/wgbs/SRR13650165_2.fastq.gz
```

```bash
# --- SRR13650194 (salt) ---
trim_galore --paired --cores 4 --quality 20 --stringency 5 --dont_gzip --clip_r1 6 --clip_r2 6 --three_prime_clip_R1 6 --three_prime_clip_R2 6 --output_dir /course/results/day2_methylation/02_trimmed /course/data/raw/wgbs/SRR13650194_1.fastq.gz /course/data/raw/wgbs/SRR13650194_2.fastq.gz
```

```bash
# --- SRR13650196 (salt) ---
trim_galore --paired --cores 4 --quality 20 --stringency 5 --dont_gzip --clip_r1 6 --clip_r2 6 --three_prime_clip_R1 6 --three_prime_clip_R2 6 --output_dir /course/results/day2_methylation/02_trimmed /course/data/raw/wgbs/SRR13650196_1.fastq.gz /course/data/raw/wgbs/SRR13650196_2.fastq.gz
```

```bash
# --- SRR13650198 (salt) ---
trim_galore --paired --cores 4 --quality 20 --stringency 5 --dont_gzip --clip_r1 6 --clip_r2 6 --three_prime_clip_R1 6 --three_prime_clip_R2 6 --output_dir /course/results/day2_methylation/02_trimmed /course/data/raw/wgbs/SRR13650198_1.fastq.gz /course/data/raw/wgbs/SRR13650198_2.fastq.gz
```

### 2.2 Examine trimming report

```bash
grep -E "Reads|Pairs|Quality|Adapter" /course/results/day2_methylation/02_trimmed/SRR13650161_1.fastq.gz_trimming_report.txt
```

> **Question 3:** Why do we clip 6 bp from the 5′ end of each read (`--clip_r1 6`
> and `--clip_r2 6`)? What systematic bias does this correct?

---

## Part 3 – Bisulfite Alignment with Bismark

### 3.1 How Bismark works

Bismark performs an **in-silico conversion** of the genome before alignment:
it creates two converted versions (C→T and G→A) and aligns reads against both,
then determines the original methylation state from mismatches vs expected conversions.

```
Unmethylated C  →  Bisulfite converts C to U  →  reads as T after PCR
Methylated C    →  Bisulfite does NOT convert  →  reads as C (retained)
```

### 3.2 Verify the Bismark genome index

```bash
ls /course/data/reference/indices/bismark/Bisulfite_Genome/
```

You should see: `CT_conversion/` and `GA_conversion/` directories.

### 3.3 Run Bismark alignment

Align **one sample at a time**. Bismark takes the paired reads with `-1` and `-2`;
`--basename` gives each sample a clean output name, so Bismark writes `<SAMPLE>_pe.bam`.

```bash
# --- SRR13650161 (control) ---
bismark --genome /course/data/reference/indices/bismark --bowtie2 --multicore 4 --output_dir /course/results/day2_methylation/03_bismark_aligned --basename SRR13650161 -1 /course/results/day2_methylation/02_trimmed/SRR13650161_1_val_1.fq -2 /course/results/day2_methylation/02_trimmed/SRR13650161_2_val_2.fq
```

```bash
# --- SRR13650163 (control) ---
bismark --genome /course/data/reference/indices/bismark --bowtie2 --multicore 4 --output_dir /course/results/day2_methylation/03_bismark_aligned --basename SRR13650163 -1 /course/results/day2_methylation/02_trimmed/SRR13650163_1_val_1.fq -2 /course/results/day2_methylation/02_trimmed/SRR13650163_2_val_2.fq
```

```bash
# --- SRR13650165 (control) ---
bismark --genome /course/data/reference/indices/bismark --bowtie2 --multicore 4 --output_dir /course/results/day2_methylation/03_bismark_aligned --basename SRR13650165 -1 /course/results/day2_methylation/02_trimmed/SRR13650165_1_val_1.fq -2 /course/results/day2_methylation/02_trimmed/SRR13650165_2_val_2.fq
```

```bash
# --- SRR13650194 (salt) ---
bismark --genome /course/data/reference/indices/bismark --bowtie2 --multicore 4 --output_dir /course/results/day2_methylation/03_bismark_aligned --basename SRR13650194 -1 /course/results/day2_methylation/02_trimmed/SRR13650194_1_val_1.fq -2 /course/results/day2_methylation/02_trimmed/SRR13650194_2_val_2.fq
```

```bash
# --- SRR13650196 (salt) ---
bismark --genome /course/data/reference/indices/bismark --bowtie2 --multicore 4 --output_dir /course/results/day2_methylation/03_bismark_aligned --basename SRR13650196 -1 /course/results/day2_methylation/02_trimmed/SRR13650196_1_val_1.fq -2 /course/results/day2_methylation/02_trimmed/SRR13650196_2_val_2.fq
```

```bash
# --- SRR13650198 (salt) ---
bismark --genome /course/data/reference/indices/bismark --bowtie2 --multicore 4 --output_dir /course/results/day2_methylation/03_bismark_aligned --basename SRR13650198 -1 /course/results/day2_methylation/02_trimmed/SRR13650198_1_val_1.fq -2 /course/results/day2_methylation/02_trimmed/SRR13650198_2_val_2.fq
```

Each alignment takes a few minutes on this Chr4-sized data.
While waiting, read the Bismark paper (Krueger & Andrews, 2011, Bioinformatics).

### 3.4 Read the alignment report

```bash
cat /course/results/day2_methylation/03_bismark_aligned/SRR13650161_PE_report.txt
```

> **Question 4:** WGBS alignment rates are typically lower than standard RNA-seq
> (often 60–80% for plants). What are the main reasons for this lower efficiency?
> What is the bisulfite conversion efficiency of this sample?

### 3.5 Remove PCR duplicates

Remove duplicates **for each sample** (input `<SAMPLE>_pe.bam` → output
`<SAMPLE>_pe.deduplicated.bam`):

```bash
# --- SRR13650161 ---
deduplicate_bismark --paired --output_dir /course/results/day2_methylation/03_bismark_aligned /course/results/day2_methylation/03_bismark_aligned/SRR13650161_pe.bam
```

```bash
# --- SRR13650163 ---
deduplicate_bismark --paired --output_dir /course/results/day2_methylation/03_bismark_aligned /course/results/day2_methylation/03_bismark_aligned/SRR13650163_pe.bam
```

```bash
# --- SRR13650165 ---
deduplicate_bismark --paired --output_dir /course/results/day2_methylation/03_bismark_aligned /course/results/day2_methylation/03_bismark_aligned/SRR13650165_pe.bam
```

```bash
# --- SRR13650194 ---
deduplicate_bismark --paired --output_dir /course/results/day2_methylation/03_bismark_aligned /course/results/day2_methylation/03_bismark_aligned/SRR13650194_pe.bam
```

```bash
# --- SRR13650196 ---
deduplicate_bismark --paired --output_dir /course/results/day2_methylation/03_bismark_aligned /course/results/day2_methylation/03_bismark_aligned/SRR13650196_pe.bam
```

```bash
# --- SRR13650198 ---
deduplicate_bismark --paired --output_dir /course/results/day2_methylation/03_bismark_aligned /course/results/day2_methylation/03_bismark_aligned/SRR13650198_pe.bam
```

> **Question 5:** Why is deduplication especially important for WGBS data?
> What would happen to methylation level estimates if PCR duplicates were retained?

---

## Part 4 – Methylation Extraction

### 4.1 Run Bismark methylation extractor

Extract methylation **for each sample**. The extractor names its output after the input
BAM, so after each run we rename the report to a clean `<SAMPLE>.CX_report.txt` that
methylKit will pick up in Part 5.

> **Note on flags:** the correct flag is `--CX` (extracts CpG + CHG + CHH), *not*
> `--CX_context`, and `bismark_methylation_extractor` has **no** `--basename` option
> (that is an aligner-only flag). `--cytosine_report` produces the per-cytosine report.

```bash
# --- SRR13650161 ---
bismark_methylation_extractor --paired-end --comprehensive --cytosine_report --CX --genome_folder /course/data/reference/indices/bismark --parallel 4 --output_dir /course/results/day2_methylation/04_methylation /course/results/day2_methylation/03_bismark_aligned/SRR13650161_pe.deduplicated.bam

mv /course/results/day2_methylation/04_methylation/SRR13650161_pe.deduplicated.CX_report.txt /course/results/day2_methylation/04_methylation/SRR13650161.CX_report.txt
```

```bash
# --- SRR13650163 ---
bismark_methylation_extractor --paired-end --comprehensive --cytosine_report --CX --genome_folder /course/data/reference/indices/bismark --parallel 4 --output_dir /course/results/day2_methylation/04_methylation /course/results/day2_methylation/03_bismark_aligned/SRR13650163_pe.deduplicated.bam

mv /course/results/day2_methylation/04_methylation/SRR13650163_pe.deduplicated.CX_report.txt /course/results/day2_methylation/04_methylation/SRR13650163.CX_report.txt
```

```bash
# --- SRR13650165 ---
bismark_methylation_extractor --paired-end --comprehensive --cytosine_report --CX --genome_folder /course/data/reference/indices/bismark --parallel 4 --output_dir /course/results/day2_methylation/04_methylation /course/results/day2_methylation/03_bismark_aligned/SRR13650165_pe.deduplicated.bam

mv /course/results/day2_methylation/04_methylation/SRR13650165_pe.deduplicated.CX_report.txt /course/results/day2_methylation/04_methylation/SRR13650165.CX_report.txt
```

```bash
# --- SRR13650194 ---
bismark_methylation_extractor --paired-end --comprehensive --cytosine_report --CX --genome_folder /course/data/reference/indices/bismark --parallel 4 --output_dir /course/results/day2_methylation/04_methylation /course/results/day2_methylation/03_bismark_aligned/SRR13650194_pe.deduplicated.bam

mv /course/results/day2_methylation/04_methylation/SRR13650194_pe.deduplicated.CX_report.txt /course/results/day2_methylation/04_methylation/SRR13650194.CX_report.txt
```

```bash
# --- SRR13650196 ---
bismark_methylation_extractor --paired-end --comprehensive --cytosine_report --CX --genome_folder /course/data/reference/indices/bismark --parallel 4 --output_dir /course/results/day2_methylation/04_methylation /course/results/day2_methylation/03_bismark_aligned/SRR13650196_pe.deduplicated.bam

mv /course/results/day2_methylation/04_methylation/SRR13650196_pe.deduplicated.CX_report.txt /course/results/day2_methylation/04_methylation/SRR13650196.CX_report.txt
```

```bash
# --- SRR13650198 ---
bismark_methylation_extractor --paired-end --comprehensive --cytosine_report --CX --genome_folder /course/data/reference/indices/bismark --parallel 4 --output_dir /course/results/day2_methylation/04_methylation /course/results/day2_methylation/03_bismark_aligned/SRR13650198_pe.deduplicated.bam

mv /course/results/day2_methylation/04_methylation/SRR13650198_pe.deduplicated.CX_report.txt /course/results/day2_methylation/04_methylation/SRR13650198.CX_report.txt
```

### 4.2 Inspect methylation output

```bash
# Cytosine report format: chr  pos  strand  count_M  count_U  context  trinucleotide
head -20 /course/results/day2_methylation/04_methylation/SRR13650161.CX_report.txt

# Count cytosines per context
awk 'NR>1 {print $6}' /course/results/day2_methylation/04_methylation/SRR13650161.CX_report.txt | sort | uniq -c
```

> **Question 6:** The cytosine report has 7 columns. What information does each
> column encode? Calculate the mean CpG methylation level for this sample
> from the command line (hint: `awk` with context filter on column 6).

```bash
# Calculate mean CpG methylation
awk '$6 == "CG" && ($4+$5) > 0 {sum += $4/($4+$5); n++}

     END {printf "Mean CpG methylation: %.2f%%\n", sum/n*100}' /course/results/day2_methylation/04_methylation/SRR13650161.CX_report.txt
```

---

## Part 5 – Differential Methylation Analysis with methylKit (R)

### 5.1 Launch R and load libraries

```bash
R
```

```r
library(methylKit)
library(ggplot2)
library(dplyr)

# Set paths
meth_dir <- "/course/results/day2_methylation/04_methylation"
out_dir  <- "/course/results/day2_methylation/05_diffmeth"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
```

### 5.2 Load cytosine reports into methylKit

We load **all 6 samples (3 control + 3 salt)** from the sample sheet, so the
differential test has real biological replication.

```r
# Read the sample sheet (sample_id + condition for all 6 samples)
meta <- read.table("/course/data/raw/wgbs/samplesheet.tsv",
                   header = TRUE, sep = "\t")

file_list  <- as.list(file.path(meth_dir,
                paste0(meta$sample_id, ".CX_report.txt")))
sample_ids <- as.list(meta$sample_id)
treatment  <- ifelse(meta$condition == "control", 0L, 1L)   # 0 = control, 1 = salt

# CpG context
myobj <- methRead(
  file_list,
  sample.id  = sample_ids,
  assembly   = "TAIR10_Chr4",
  treatment  = treatment,
  context    = "CpG",
  mincov     = 1,     # keep every covered cytosine. This is a subsampled, low-
                      # depth teaching set — a high per-cytosine minimum leaves
                      # nothing shared across all 6 samples. We aggregate depth
                      # regionally (200 bp tiles) in step 5.6.
  pipeline   = "bismarkCytosineReport"
)

# Inspect the first sample's methylation distribution
getMethylationStats(myobj[[1]], plot = TRUE, both.strands = FALSE)
```

> **Question 7:** The methylation histogram shows a bimodal distribution
> (peaks near 0% and 100%). What biological meaning does each peak carry?
> Is this pattern expected in plants?

### 5.3 Filter by coverage and unite samples

```r
# Filter: min 3x, remove outlier high-coverage positions.
# (3x, not the textbook 10x — see the mincov note above: at this depth 10x in
#  every sample leaves no shared cytosine to unite.)
myobj_filtered <- filterByCoverage(myobj,
  lo.count = 3,
  hi.perc  = 99.9
)

# Unite: retain only positions covered in ALL 6 samples
meth <- unite(myobj_filtered, destrand = TRUE)
cat("CpG sites after filtering:", nrow(meth), "\n")   # ~510 on this teaching set
```

### 5.4 Sample-level QC

```r
# Correlation between samples
getCorrelation(meth, plot = TRUE)

# PCA
PCASamples(meth)

# Hierarchical clustering
clusterSamples(meth, dist = "correlation", method = "ward", plot = TRUE)
```

> **Question 8:** What Pearson correlation coefficient do you observe between
> control and salt-treated samples? Is a high (>0.95) or low (<0.80) correlation
> what you would expect, and why?

### 5.5 Calculate differential methylation

> We use `mc.cores = 1` (serial). methylKit's parallel backend (`mclapply`) can
> **deadlock** on some Docker setups — worse on many-core machines — so serial is the
> reliable choice; on this Chr4-sized data it is plenty fast.

```r
myDiff <- calculateDiffMeth(meth,
  mc.cores       = 1,
  test           = "F",
  overdispersion = "MN"    # matches the automated day2_methylkit.R
)

# Select DMCs: |Δmeth| ≥ 25%, q < 0.01
myDiff_25 <- getMethylDiff(myDiff,
  difference = 25,
  qvalue     = 0.01,
  type       = "all"
)

cat("Total DMCs  :", nrow(myDiff_25), "\n")
cat("Hypermethylated:", nrow(getMethylDiff(myDiff, 25, 0.01, "hyper")), "\n")
cat("Hypomethylated :", nrow(getMethylDiff(myDiff, 25, 0.01, "hypo")), "\n")

# Visualise the full per-cytosine differential (all tested Cs, not only the
# significant ones — on this shallow teaching set the single-cytosine test
# often returns 0 significant DMCs, so the regional analysis in 5.6 is the
# primary result).
plot(myDiff, chromosome = "Chr4",
     col = c("firebrick", "steelblue"), lwd = 1.5,
     main = "CpG differential methylation – salt vs control (Chr4)")
```

> **Note on depth.** Because this is a subsampled dataset, a strict single-
> cytosine test (|Δ| ≥ 25%, q < 0.01) typically yields **0 DMCs** — there simply
> is not enough per-base coverage. This is expected, and it is exactly why we
> move to **regional (tiled) DMRs** below, which pool coverage across each window.

### 5.6 Tiling window analysis (DMRs)

```r
# Tile the genome into 200 bp non-overlapping bins.
# IMPORTANT: tile the RAW per-sample object `myobj` (from 5.2) and unite the
# tiles AFTERWARDS. This pools every read inside a window before we ask for the
# window to be shared across samples. (Tiling the already-united single-cytosine
# object `meth` would inherit its sparsity and return no DMRs on this dataset.)
tiles      <- tileMethylCounts(myobj, win.size = 200, step.size = 200, cov.bases = 1)
meth_tiles <- unite(tiles, destrand = FALSE)

myDiff_tiles <- calculateDiffMeth(meth_tiles, mc.cores = 1, test = "F", overdispersion = "MN")

DMRs <- getMethylDiff(myDiff_tiles,
  difference = 20,
  qvalue     = 0.05,
  type       = "all"
)

cat("DMRs (200 bp tiles):", nrow(DMRs), "\n")

# Export as BED for IGV
DMR_df <- as.data.frame(DMRs)
write.table(
  data.frame(
    chr   = DMR_df$chr,
    start = DMR_df$start - 1,
    end   = DMR_df$end,
    name  = paste0("DMR_", seq_len(nrow(DMR_df))),
    score = round(DMR_df$meth.diff, 1),
    strand = "."
  ),
  file.path(out_dir, "DMRs_CpG_200bp.bed"),
  sep = "\t", quote = FALSE,
  row.names = FALSE, col.names = FALSE
)
cat("BED file written:", file.path(out_dir, "DMRs_CpG_200bp.bed"), "\n")
```

### 5.7 Context-level comparison (CpG vs CHG vs CHH)

```r
# Run the same analysis for CHG (all 6 samples, reusing file_list/sample_ids/treatment)
myobj_CHG <- methRead(file_list,
  sample.id = sample_ids,
  assembly  = "TAIR10_Chr4",
  treatment = treatment,
  context   = "CHG",
  mincov    = 1,     # low-depth teaching set — see the CpG note in 5.2
  pipeline  = "bismarkCytosineReport"
)
meth_CHG <- unite(filterByCoverage(myobj_CHG, lo.count = 3, hi.perc = 99.9),
  destrand = FALSE)
myDiff_CHG <- calculateDiffMeth(meth_CHG, mc.cores = 1, test = "F", overdispersion = "MN")
DMC_CHG <- getMethylDiff(myDiff_CHG, 25, 0.01, "all")

# And CHH
myobj_CHH <- methRead(file_list,
  sample.id = sample_ids,
  assembly  = "TAIR10_Chr4",
  treatment = treatment,
  context   = "CHH",
  mincov    = 1,     # low-depth teaching set — see the CpG note in 5.2
  pipeline  = "bismarkCytosineReport"
)
meth_CHH <- unite(filterByCoverage(myobj_CHH, lo.count = 3, hi.perc = 99.9),
  destrand = FALSE)
myDiff_CHH <- calculateDiffMeth(meth_CHH, mc.cores = 1, test = "F", overdispersion = "MN")
DMC_CHH <- getMethylDiff(myDiff_CHH, 25, 0.01, "all")

# Summary bar plot
summary_df <- data.frame(
  Context = c("CpG", "CHG", "CHH"),
  DMCs    = c(nrow(myDiff_25), nrow(DMC_CHG), nrow(DMC_CHH))
)
ggplot(summary_df, aes(x = Context, y = DMCs, fill = Context)) +
  geom_col() +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "DMCs per cytosine context – salt vs control",
       y = "Number of DMCs") +
  theme_classic()
```

> **Question 9:** In which cytosine context do you observe the most DMCs?
> Given the role of RdDM (24 nt hc-siRNAs → DRM2 → CHH methylation), propose a
> mechanistic model linking small-RNA changes to altered CHH methylation under salt
> stress. Separately — and this is the link Part 6 tests directly — could a *cis*-acting
> lncRNA, rather than a small RNA, recruit the same methylation machinery to a
> neighbouring locus instead? What evidence would distinguish the two mechanisms?

---

## Part 6 – Integrative Discussion

### Connecting Day 1 and Day 2 Results

On Day 1 we found **differentially expressed lncRNAs** on Chr4; today we found
**differentially methylated regions (DMRs)** on the same chromosome. A natural question is
whether these two regulatory marks tend to sit in the **same places** on the genome.

> **Heads-up — these are different experiments.** Day 1 used Tr-PET micro-nanoplastic
> stress and Day 2 uses salt stress, so this is a **positional co-localisation** exercise
> ("do the transcriptomic and epigenetic marks land in the same regions of Chr4?"), not a
> single joint differential test. It still teaches the mechanics of integrating two
> omics layers.

Day 1 wrote the coordinates of its DE lncRNAs to `10_cis/de_lnc.bed`. That file is visible
here as `/course/results/day1_lncRNA/10_cis/de_lnc.bed`, because both days share the same
`./results` folder on your host machine.

> **One catch — chromosome names.** Day 1's reference calls the chromosome `4`, while
> Day 2's calls it `Chr4`. If we compare them as-is, every overlap silently returns zero —
> so we rename Day 1's `4` to `Chr4` first.

```r
library(GenomicRanges)

# --- Day 1 output: differentially expressed lncRNA loci ---
lnc <- read.table("/course/results/day1_lncRNA/10_cis/de_lnc.bed", sep = "\t")
lnc_chr <- ifelse(lnc$V1 == "4", "Chr4", as.character(lnc$V1))   # harmonise chr name
lnc_GR  <- GRanges(seqnames = lnc_chr,
                   ranges   = IRanges(start = lnc$V2 + 1, end = lnc$V3))

# --- Day 2 output: differentially methylated regions (BED written in Part 5.6) ---
dmr <- read.table(file.path(out_dir, "DMRs_CpG_200bp.bed"), sep = "\t")
dmr_GR <- GRanges(seqnames = as.character(dmr$V1),
                  ranges   = IRanges(start = dmr$V2 + 1, end = dmr$V3))

# sanity check: both objects must use the same chromosome name, or overlaps will be empty
cat("Day 1 lncRNA contigs:", paste(unique(as.character(seqnames(lnc_GR))), collapse = ","), "\n")
cat("Day 2 DMR contigs   :", paste(unique(as.character(seqnames(dmr_GR))), collapse = ","), "\n")

# --- Do DE lncRNAs sit within 2 kb of a CpG DMR? ---
near <- findOverlaps(lnc_GR, dmr_GR, maxgap = 2000)
cat("DE lncRNAs within 2 kb of a CpG DMR:",
    length(unique(queryHits(near))), "of", length(lnc_GR), "\n")
```

> **Question 10:** How many DE lncRNAs sit within 2 kb of a differentially methylated
> region? Because the two datasets come from different stresses, what would a positive
> overlap actually tell you — and what single experiment (the *same* plants, assayed for
> both RNA and methylation) would you run to test whether lncRNA expression and local DNA
> methylation are *causally* linked?

---

## Shortcut – run the whole Day 2 pipeline at once

If you fall behind, the entire Day 2 pipeline (QC → trimming → Bismark alignment →
methylation extraction → methylKit) can be run non-interactively inside the course
container:

```bash
bash /course/scripts/run_day2.sh 4
```

---

## Final Questions for Group Discussion

1. **Stress memory:** Plants subjected to recurring stress episodes often show
   altered DNA methylation patterns at stress-responsive loci. Based on today's
   data, what genomic regions would you investigate further to test for
   "methylation memory" of salt stress?

2. **Heritability:** The WGBS data captures the methylome of somatic tissue.
   How would you design an experiment to determine whether the observed DMRs
   are transmitted to the next generation (transgenerational epigenetic inheritance)?

3. **Statistical power:** This course uses 3 biological replicates per condition —
   the minimum that lets methylKit estimate within-group variance and apply the
   logistic-regression test with overdispersion correction. Why do 3 replicates
   change the analysis compared with a single sample per group (where methylKit
   falls back to Fisher's exact test)? For a genome-wide study, how many replicates
   would you budget for, and why?

4. **Methylation and gene expression:** Propose how you would integrate the WGBS
   data from today with a standard RNA-seq dataset to test whether DMRs in gene
   promoters correlate with changes in gene expression.

---

## Closing Remarks

Congratulations on completing the Epi-Code practical course.

You have now implemented two complete bioinformatics workflows:
- **Day 1:** RNA-seq → QC → trim → HISAT2 → StringTie → gffcompare/CPC2 → featureCounts → edgeR
- **Day 2:** WGBS → QC → trim → Bismark align → methylation extraction → methylKit

These pipelines, containerised in Docker, are fully reproducible and can be
scaled to complete genomes on HPC systems by removing the Chr4 restriction.

All results are saved in `./results/` on your host machine.

Thank you for your participation and curiosity!

*D. Giosa & L. Giuffrè — University of Messina, ChiBioFarAm Dept.*
