# Commands Explained — what every step actually does

This guide walks through **every command** you run in the two practicals and explains
*what it does*, *why the specific options are there*, and *what the output means
biologically*. Read it alongside the walkthroughs
([Day 1](exercises/Day1_lncRNA_mRNA_Pipeline.md),
[Day 2](exercises/Day2_DNA_Methylation_Analysis.md)) — the goal is that by the end you can
read a bioinformatics pipeline and understand *why* each tool is where it is, not just
copy-paste it.

> Throughout, the analysis is deliberately kept small: *Arabidopsis thaliana*
> **Chromosome 4** only, reads subsampled, so everything finishes in minutes. The *logic*
> is exactly the same as a whole-genome study on a server.

---

## 0. Core concepts you need first

| Term | What it means |
|------|---------------|
| **Read** | One short DNA sequence produced by the sequencer (here ~100 bp). |
| **FASTQ** | Text file with, per read: an ID line, the bases, a `+`, and a per-base **quality** string (Phred score). Gzipped as `.fastq.gz`. |
| **Paired-end (R1/R2)** | The fragment is read from both ends → two files `_R1`/`_R2` (or `_1`/`_2`). The two mates belong together and must stay in sync. |
| **Adapter** | Short synthetic sequence ligated to every fragment during library prep. If the fragment is short, the sequencer reads *into* the adapter, so we trim it off. |
| **QC** | Quality control — looking at the reads *before* analysis to catch problems (low quality tails, adapter contamination, wrong base composition). |
| **Reference genome / annotation** | The known DNA sequence (`.fa`) and the catalogue of known genes on it (`.gtf`/`.gff`). We map our reads onto this coordinate system. |
| **Alignment (mapping)** | Finding where each read came from on the reference. Output is a **BAM** file (compressed, sorted, indexed table of aligned reads). |
| **Assembly** | Reconstructing the transcripts that were actually expressed *from* the aligned reads (StringTie) — needed to discover **new** transcripts the annotation doesn't list. |
| **Counts / quantification** | How many reads fall on each gene → a **count matrix** (genes × samples). The raw material for differential expression. |
| **Normalisation** | Correcting for the fact that libraries have different total read numbers and composition, so counts are comparable across samples. |
| **Replicates** | Independent biological samples per condition (here **3 vs 3**). They let statistics estimate *natural variation*, which is what makes a difference "significant". |
| **Multiple-testing correction (FDR / q-value)** | We test thousands of genes/positions at once; some look "significant" by chance. FDR/q-value control that false-positive rate. |
| **"Differential"** | The whole point: a feature (gene, cytosine, region) whose signal **changes between conditions** more than expected from noise. |

**The shape of both days** is the same classic pattern:

```
raw reads → QC → trim → align → summarise per feature → normalise → test between conditions → interpret
```

Day 1 summarises reads into **transcripts/genes** (expression); Day 2 summarises them into
**methylated cytosines** (epigenetic marks). Everything else is analogous.

---

# Day 1 — lncRNA & mRNA (RNA-seq)

**Biological question:** which genes — including previously *unknown* long non-coding RNAs
— change expression when *Arabidopsis* roots are exposed to PET micro-nanoplastics
(Tr-PET) vs control? Reproduces Galbo, Giosa et al. (2025).

**Pipeline:** `FastQC → Trimmomatic → HISAT2 → StringTie → gffcompare + CPC2 →
featureCounts → edgeR → visualisation → lncRNA cis-targets + GO`.

> **One idea that runs through Day 1: strandedness.** This is a *stranded* library, so we
> know which DNA strand each transcript came from. That is essential for lncRNAs, because a
> lncRNA is often transcribed **antisense** to a gene at the same locus — only strand info
> tells them apart. Every tool is told the library is reverse-stranded: HISAT2
> `--rna-strandness RF`, StringTie `--rf`, featureCounts `-s 2`.

### `fastqc` — read quality control
```bash
fastqc -t 4 -o results/01_qc <all 12 FASTQ files>
```
- Reads each FASTQ and produces an HTML report of per-base quality, adapter content, GC%,
  duplication, etc. `-t 4` = 4 files in parallel; `-o` = output folder.
- **Why:** you never trust reads you haven't looked at. For RNA-seq, some sequence
  *duplication is normal* (highly expressed genes are sequenced many times) — not a problem.

### `multiqc` — one report from many
```bash
multiqc -f -o results/01_qc results/01_qc
```
- Scans a folder for tool outputs (FastQC here) and merges them into **one** interactive
  report so you compare all samples on one page. `-f` overwrites an old report.

### `trimmomatic PE` — adapter & quality trimming
```bash
trimmomatic PE -threads 4 -phred33 R1.fastq.gz R2.fastq.gz \
  R1.paired.fq.gz R1.unpaired.fq.gz R2.paired.fq.gz R2.unpaired.fq.gz \
  ILLUMINACLIP:adapters.fa:2:30:10 HEADCROP:1 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:18 MINLEN:40
```
- `PE` = paired-end; it writes 4 files: for each mate a **paired** output (both mates
  survived) and an **unpaired** one (only this mate survived). We keep the **paired** files
  because the aligner needs both mates.
- `ILLUMINACLIP:adapters.fa:2:30:10` — remove Illumina adapters (allow 2 mismatches; the
  30/10 control how aggressively partial/palindromic adapters are clipped).
- `HEADCROP:1` drop the first base; `LEADING:3`/`TRAILING:3` trim ends below quality 3;
  `SLIDINGWINDOW:4:18` scan a 4-base window and cut once mean quality drops below 18;
  `MINLEN:40` discard reads shorter than 40 bp (too short to map uniquely).
- **Output line to note:** `Input Read Pairs … Both Surviving XX%` — the fraction kept.

### `hisat2` (+ `samtools sort`/`index`) — splice-aware alignment
```bash
hisat2 -p 4 -x reference/hisat2_index/chr4 -1 R1.paired.fq.gz -2 R2.paired.fq.gz \
  --rna-strandness RF --dta --known-splicesite-infile chr4.ss.txt \
  --new-summary --summary-file sample.hisat2.summary \
  | samtools sort -@ 4 -o sample.sorted.bam -
samtools index sample.sorted.bam
```
- **HISAT2** maps RNA-seq reads. RNA reads cross **exon–exon junctions** (introns are
  spliced out of mRNA but present in the genome), so a normal DNA aligner would fail at
  junctions — HISAT2 is **splice-aware** and maps a read across an intron.
- `-x` = the pre-built genome index; `-1/-2` = the paired reads.
- `--rna-strandness RF` = reverse-stranded library (see the strandedness box).
- `--dta` = "downstream transcriptome assembly" — makes alignments tidy enough for
  StringTie to assemble transcripts from them.
- `--known-splicesite-infile` gives HISAT2 the *known* junctions as hints (faster, more
  accurate); `--new-summary --summary-file` writes the alignment-rate report.
- **The pipe `| samtools sort`** turns HISAT2's stream directly into a **coordinate-sorted
  BAM** (no huge intermediate SAM); `samtools index` builds the `.bai` so tools can jump to
  any region instantly.
- **Read the result:** `grep 'Overall alignment rate' *.summary`. High here because the
  reads were pre-selected to Chr4.

### `stringtie` (per sample, then `--merge`) — transcript assembly
```bash
stringtie sample.sorted.bam --rf -p 4 -G reference/chr4.gtf -o sample.gtf -l sample
# then, across all 6:
stringtie --merge -p 4 -G reference/chr4.gtf -o merged.gtf sample1.gtf … sample6.gtf
```
- **StringTie** reconstructs the transcripts actually present in each sample from the
  aligned reads. `-G chr4.gtf` supplies the known annotation as a **guide** — known
  transcripts are recovered *and* genuinely new ones are still reported (this is how we
  find novel lncRNAs). `--rf` = reverse-stranded; `-l` sets a name prefix.
- **`--merge`** combines the six per-sample assemblies into **one non-redundant
  transcriptome**, so every sample is later quantified against the *same* set of
  transcripts. Novel transcripts get IDs like `MSTRG.*`.
- **Why guided, not de-novo?** Without `-G`, you would still find transcripts but lose the
  link to known gene names and get more fragmented assemblies.

### `gffcompare` — classify each transcript vs the annotation
```bash
gffcompare -r reference/chr4.gtf -o results/05_novel/gffcmp merged.gtf
```
- Compares every assembled transcript to the reference and tags it with a **class code**.
  Key codes for lncRNA discovery: **`u`** = intergenic (in a gene desert), **`i`** =
  intronic, **`x`** = antisense (overlaps a gene on the *opposite* strand), `=` = matches a
  known transcript, `j` = novel isoform of a known gene.
- We keep `u`, `i`, `x` (candidates that are **not** already annotated genes on the same
  strand) — the raw material for novel lncRNAs.

### `gtf_select.py` + `gffread` — pick candidates & get their sequence
```bash
python3 lib/gtf_select.py gffcmp.annotated.gtf --classes u,i,x --min-len 200 \
  --out-gtf candidates.gtf --out-tsv candidates.tsv
gffread -w candidates.fa -g reference/chr4.fa candidates.gtf
```
- The helper keeps class `u/i/x` transcripts **≥ 200 nt** (the length threshold in the
  definition of a *long* non-coding RNA).
- **`gffread -w`** extracts the **spliced transcript sequence** (FASTA) for each candidate
  from the genome — needed for the coding-potential test next.

### `CPC2` — is it really non-coding?
```bash
CPC2 -i candidates.fa -o results/06_lncrna/cpc2
```
- **Coding Potential Calculator 2** scores each transcript's likelihood of being
  translated (ORF length/quality, Fickett score, isoelectric point…) and labels it
  **coding** or **noncoding**.
- **Why needed** (a subtle but important point): "unannotated" does **not** mean
  "non-coding" — a candidate could be an unannotated protein-coding gene. CPC2 removes
  those, so what remains is a confident lncRNA set. We keep only the **noncoding** ones.

### `build_augmented.py` — one annotation with genes + new lncRNAs
```bash
python3 lib/build_augmented.py --ref-gtf chr4.gtf --lnc-gtf novel_lncRNA.gtf \
  --out-gtf augmented.gtf --out-tsv gene_biotype.tsv
```
- Merges known genes + confirmed novel lncRNAs into a single **augmented** annotation and
  writes a table labelling each gene `mRNA` / `lncRNA` / `other`. Everything downstream
  (counting, DE) uses this so mRNAs and lncRNAs are analysed **together**.

### `featureCounts` — build the count matrix
```bash
featureCounts -T 4 -p --countReadPairs -s 2 -t exon -g gene_id -a augmented.gtf \
  -o counts.txt  <all 6 sorted BAMs>
```
- Counts how many fragments land on each gene, for **all six samples at once** → one table,
  six count columns (the matrix DE needs).
- `-p --countReadPairs` = count **fragments** (read pairs), not individual mates.
- `-s 2` = **reverse-stranded** counting — critical so an antisense lncRNA's reads are not
  mis-assigned to its sense neighbour.
- `-t exon -g gene_id` = count reads over `exon` features, grouped by `gene_id`.
- The `.summary` file shows the **assignment rate** (what % of reads were counted).

### `edgeR` (in R) — differential expression
```r
y <- DGEList(counts = counts, group = group)
keep <- rowSums(cpm(y) > 1) >= 2          # drop genes barely expressed anywhere
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)                    # TMM normalisation
design <- model.matrix(~ group)
y   <- estimateDisp(y, design)             # estimate biological variability
fit <- glmFit(y, design); lrt <- glmLRT(fit, coef = 2)   # test TrPET vs control
res <- topTags(lrt, n = Inf)$table
```
- **edgeR** models RNA-seq counts with a **negative-binomial** distribution (counts are
  discrete and over-dispersed — variance > mean).
- **Filter first** (`cpm > 1 in ≥ 2 libraries`): genes with almost no reads carry no
  information and only cost you statistical power after multiple-testing correction.
- **`calcNormFactors` = TMM**: corrects for library size *and composition* so a few
  very-high genes don't make everything else look "down".
- `estimateDisp` learns how much genes naturally vary between replicates; `glmFit` +
  `glmLRT` run the **likelihood-ratio test** of Tr-PET vs control. `topTags` returns
  `logFC` (fold change) and `FDR` (corrected p-value).
- **Significance here:** `|log2FC| ≥ 0.5` and `FDR < 0.05` (the paper's thresholds). We
  then split hits into `mRNA` vs `lncRNA` using the biotype table.

### Visualisation — PCA / volcano / heatmap (ggplot2, pheatmap)
- **PCA** (`prcomp` on log-CPM): do the samples separate by condition? If control and
  Tr-PET split on PC1, the treatment has a genome-wide effect.
- **Volcano** (`logFC` vs `-log10 FDR`): each dot a gene; up/down-regulated significant
  genes fly out to the top corners. Colouring by lncRNA/mRNA shows whether lncRNAs are
  over-represented among strong responders.
- **Heatmap** (z-scored log-CPM of top DE genes): visual clustering of samples and genes.

### `bedtools closest` + `clusterProfiler::enrichGO` — lncRNA cis-targets & GO
```r
system("bedtools closest -a de_lnc.bed -b pc.bed -d -t first > closest.tsv")
ego <- enrichGO(gene = unique(cis$target_gene), OrgDb = org.At.tair.db,
                keyType = "TAIR", ont = "BP", pvalueCutoff = 0.05)
```
- lncRNAs often regulate their **neighbouring** protein-coding gene (*in cis*).
  **`bedtools closest`** finds, for each DE lncRNA, the nearest protein-coding gene and the
  distance; we keep pairs within 10 kb.
- **`enrichGO`** asks: are those neighbour genes enriched for particular **GO biological
  processes** (stress, transport, cell wall…)? This turns a positional prediction into a
  biological hypothesis. `ont="BP"` = biological process; `pvalueCutoff` filters terms.
- **Caveat taught here:** this is a *positional* prediction of regulation, not proof —
  testing causation needs a wet-lab experiment (e.g. knock down the lncRNA).

### `Ballgown` (bonus) — isoform-level analysis
- Re-quantifies with StringTie `-e -B` (estimate abundances only, output Ballgown tables)
  and tests differences at the level of **individual transcript isoforms** rather than
  whole genes — a finer-grained view of the same data.

---

# Day 2 — DNA methylation (WGBS)

**Biological question:** where does the DNA methylation pattern change when *Arabidopsis*
is salt-stressed vs control, across the three plant methylation contexts CpG/CHG/CHH?

**Pipeline:** `FastQC → Trim Galore → Bismark (align → deduplicate → extract) → methylKit`.

> **The key trick of WGBS:** sodium **bisulfite** converts **unmethylated C → U (→ T after
> PCR)**, while **methylated C stays C**. So after sequencing, a position that reads `C`
> was methylated and one that reads `T` was not. Every tool below exists to handle the fact
> that the reads no longer match the genome one-to-one.

### One-time setup: `00_setup_reference.sh` → `bismark_genome_preparation`
```bash
bash $SCRIPTS_DIR/00_setup_reference.sh 4
```
- Downloads Chr4 + annotation and builds the aligner indices. For WGBS the crucial one is
  the **Bismark bisulfite genome**: it creates two *in-silico converted* copies of the
  genome — a **C→T** version and a **G→A** version — because a bisulfite read must be
  compared against a converted reference, not the original.
- Produces `Bisulfite_Genome/CT_conversion/` and `GA_conversion/`.

### `download_data.sh` — get the reads (disk-safe)
```bash
bash $DATA_DIR/download_data.sh
```
- Streams the 6 WGBS libraries from ENA and keeps only the first N read-pairs (a full WGBS
  run is ~6 GB; we store a few hundred MB). Writes a `samplesheet.tsv` (sample → condition
  → files) that the R step reads.

### `fastqc` on bisulfite reads
- Same tool as Day 1, but **expect "Per base sequence content" to FAIL** — and that's
  normal. Bisulfite conversion makes the library cytosine-poor / thymine-rich, so the base
  composition is deliberately skewed. FastQC's pass/fail rules assume ordinary DNA.

### `trim_galore --paired … --clip_r1 6 --clip_r2 6 …`
```bash
trim_galore --paired --cores 4 --quality 20 --stringency 5 --dont_gzip \
  --clip_r1 6 --clip_r2 6 --three_prime_clip_R1 6 --three_prime_clip_R2 6 \
  --output_dir 02_trimmed  R1.fastq.gz R2.fastq.gz
```
- Trims adapters + low-quality bases (`--quality 20`), `--paired` keeps mates concordant.
- **The WGBS-specific bit:** `--clip_r1 6 --clip_r2 6` (and the 3′ clips) cut 6 bp off the
  read ends. WGBS libraries have a **methylation bias in the first few bases** from
  random-hexamer priming/end-repair (called *M-bias*); clipping removes that artefact so it
  doesn't distort methylation calls. Output files are named `*_val_1.fq` / `*_val_2.fq`.

### `bismark` — bisulfite-aware alignment
```bash
bismark --genome reference/indices/bismark --bowtie2 --multicore 4 \
  --output_dir 03_bismark_aligned --basename SAMPLE \
  -1 SAMPLE_1_val_1.fq -2 SAMPLE_2_val_2.fq
```
- Aligns the (bisulfite-converted) reads against both converted genomes using Bowtie2 under
  the hood, then infers each cytosine's original state from the pattern of matches. `-1/-2`
  = paired reads; `--multicore 4` = parallelism; `--basename SAMPLE` → output `SAMPLE_pe.bam`.
- **Read the report** `SAMPLE_PE_report.txt`: the **mapping efficiency** (WGBS is lower
  than RNA-seq — the converted genome has less sequence complexity) and the **bisulfite
  conversion rate** (should be ≥ ~99%; measured from non-CpG C's that should all be
  unmethylated).

### `deduplicate_bismark --paired`
```bash
deduplicate_bismark --paired --output_dir 03_bismark_aligned SAMPLE_pe.bam
```
- Removes **PCR duplicates** (read pairs mapping to the exact same position — copies of one
  original molecule). **Why it matters for WGBS:** a duplicate carries the *same* methylation
  call repeated; keeping duplicates would let a single molecule vote many times and bias the
  methylation percentage. Output: `SAMPLE_pe.deduplicated.bam`.

### `bismark_methylation_extractor … --cytosine_report --CX`
```bash
bismark_methylation_extractor --paired-end --comprehensive --cytosine_report --CX \
  --genome_folder reference/indices/bismark --parallel 4 \
  --output_dir 04_methylation  SAMPLE_pe.deduplicated.bam
# then rename SAMPLE_pe.deduplicated.CX_report.txt → SAMPLE.CX_report.txt
```
- Walks every aligned read and records, for each cytosine, how many reads said **methylated
  (C)** vs **unmethylated (T)**.
- `--CX` = report **all three contexts** CpG, CHG, CHH (plant methylation uses all three;
  animals mostly CpG). `--cytosine_report` = a genome-wide per-cytosine table
  (`chr pos strand countM countU context trinucleotide`). `--comprehensive` merges strands
  per context.
- We rename the report to `SAMPLE.CX_report.txt` so methylKit finds it by sample name.

### `methylKit` (in R) — differential methylation
```r
myobj <- methRead(file_list, sample.id, treatment, context="CpG", mincov=10,
                  pipeline="bismarkCytosineReport")
meth  <- unite(filterByCoverage(myobj, lo.count=10, hi.perc=99.9), destrand=TRUE)
myDiff <- calculateDiffMeth(meth, test="F", overdispersion="MN")
myDiff25 <- getMethylDiff(myDiff, difference=25, qvalue=0.01, type="all")
tiles  <- tileMethylCounts(meth, win.size=200, step.size=200)   # → DMRs
```
- **`methRead`** loads the per-cytosine reports; `treatment = 0/1` encodes control vs salt;
  `mincov=10` ignores positions with < 10× coverage (too little data to trust).
- **`filterByCoverage`** drops low-coverage and PCR-inflated (top 0.1%) positions;
  **`unite`** keeps only cytosines covered in **all** samples so every position is
  comparable. `destrand=TRUE` (CpG only) merges the two symmetric strands.
- **`calculateDiffMeth`** tests each cytosine for a methylation difference between
  conditions. With **3 replicates** it uses logistic regression with **overdispersion
  correction** (`test="F", overdispersion="MN"`) — the correct model when you have
  biological variation. (With a single sample per group it would silently fall back to
  Fisher's exact test — one reason replicates matter.)
- **`getMethylDiff`** selects **DMCs** — differentially methylated cytosines — at
  |Δmethylation| ≥ 25% and q < 0.01 (q = FDR-corrected p). `hyper`/`hypo` = gained/lost
  methylation in salt vs control.
- **`tileMethylCounts` → DMRs**: instead of single cytosines, sum methylation in 200 bp
  windows and test those — **differentially methylated regions**, more robust and easier to
  interpret than isolated cytosines. Exported as a **BED** file for genome browsers (IGV).
- **Context comparison (CpG/CHG/CHH):** the same test is run per context. In plants, the
  **CHH** context is the most stress-responsive because it is maintained *de novo* by the
  **RdDM** pathway (24 nt small RNAs → DRM2) — the mechanistic bridge back to Day 1's
  regulatory-RNA layer.

---

## Quick glossary of flags you'll see repeatedly

| Flag | Meaning |
|------|---------|
| `-t` / `-p` / `-T` / `--cores` / `--parallel` / `--multicore` | number of CPU threads (use 2 on a small laptop) |
| `-o` / `--output_dir` / `-out` | where to write results |
| `--rna-strandness RF`, `--rf`, `-s 2` | tell the tool the library is **reverse-stranded** |
| `-G` (StringTie) / `-a` (featureCounts) / `-r` (gffcompare) | supply the reference **annotation** (GTF) |
| `-x` (HISAT2) / `--genome` (Bismark) | the pre-built genome **index** |
| `FDR` / `q-value` | p-value **corrected** for testing many features at once |
| `logFC` / `meth.diff` | size of the change (expression fold-change / methylation % difference) |
| `mincov` / `MINLEN` / `cpm>1` | minimum-evidence thresholds that drop unreliable features |

---

*If a command ever fails, read its message top-to-bottom — the tools are usually explicit
about what they didn't like (a missing file, a wrong path, a strandedness mismatch). That
habit — read the error, fix the one thing, re-run — is most of practical bioinformatics.*
