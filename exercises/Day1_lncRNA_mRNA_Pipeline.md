# Day 1 – Practical Exercise
## lncRNA & mRNA Analysis: Identification and Differential Expression
### The Epi-Code – Florence Training School 2026
#### Instructors: Domenico Giosa & Letterio Giuffrè | University of Messina

---

## Learning Objectives

By the end of Day 1 you will be able to:

1. Assess the quality of strand-specific RNA-seq data with FastQC and MultiQC.
2. Trim adapters/low-quality bases with Trimmomatic and align RNA-seq reads to a plant genome with the splice-aware aligner HISAT2.
3. Reconstruct expressed transcripts with StringTie and discover **novel, unannotated long non-coding RNAs (lncRNAs)** with gffcompare.
4. Separate non-coding from coding transcripts using a coding-potential test (CPC2).
5. Build a gene-level count matrix with featureCounts and run **differential expression** for both mRNAs and lncRNAs with edgeR.
6. Relate differentially expressed lncRNAs to their neighbouring protein-coding genes (*cis*-targets) and interpret the result with GO enrichment.

This exercise reproduces, on a laptop-sized dataset, the workflow of
**Galbo, Giosa et al. (2025)** — *"lncRNA–mRNA–miRNA Networks in Arabidopsis thaliana
Exposed to Micro-Nanoplastics"*, Int. J. Plant Biol. 16, 70
([doi:10.3390/ijpb16020070](https://doi.org/10.3390/ijpb16020070)).

---

## Biological Background

Micro- and nano-plastics are an emerging environmental stress. Plants respond by
re-shaping their transcriptome, including **long non-coding RNAs** — transcripts longer
than 200 nt with little or no protein-coding potential that fine-tune the expression of
other genes, frequently their genomic neighbours.

Today's dataset (ArrayExpress **E-MTAB-13532** / ENA **PRJEB69526**) compares
*Arabidopsis thaliana* (Col-0) **roots** exposed to transparent PET micro-nanoplastics
(**Tr-PET**) against untreated **controls** — **3 biological replicates per condition**
(6 libraries), sequenced as **strand-specific** paired-end RNA-seq (TruSeq Stranded mRNA,
Illumina NovaSeq 2×100 bp).

The reference is ***A. thaliana* TAIR10, Chromosome 4** (~18 Mb, Araport11 annotation),
and the reads have been subsampled to Chr4, so the whole pipeline runs in minutes on a
laptop.

**Why "strand-specific" matters.** A stranded library records which DNA strand each transcript came from. This is essential for lncRNAs, because a lncRNA is frequently transcribed **antisense** to a protein-coding gene at the *same* locus — only  strandedness lets us tell the two apart. Throughout today we tell every tool the library is **reverse-stranded** (dUTP): HISAT2 `--rna-strandness RF`, StringTie `--rf`, featureCounts `-s 2`.

| RNA class            | Typical length | Role                                                  |
|----------------------|----------------|-------------------------------------------------------|
| mRNA (protein-coding)| variable       | Translated into protein                               |
| lincRNA (intergenic) | > 200 nt       | Regulation in *trans* / of neighbours (*cis*)         |
| antisense lncRNA     | > 200 nt       | Regulates the overlapping sense gene                  |
| intronic lncRNA      | > 200 nt       | Host-gene / splicing regulation                       |

---

## Part 0 – Start the Day 1 environment

Day 1 runs in a **self-contained image** that already includes every tool, the 6 Chr4
FASTQ libraries, the reference and the HISAT2 index — **nothing to download** beyond the
image itself.

On your **host machine** (a normal terminal, *not* inside a container):

```bash
docker pull leogiuffre/lncrna-mnps-workshop:1.0

docker run --rm -p 8888:8888 -e JAVA_TOOL_OPTIONS=-XX:TieredStopAtLevel=1 -e SHELL=/bin/bash leogiuffre/lncrna-mnps-workshop:1.0
```

**Why `-e JAVA_TOOL_OPTIONS=…`?** Trimmomatic (Part 2) is an old Java program and crashes under this image's newer Java compiler; the flag disables the faulty compiler tier (no measurable slowdown). Keep it on every `docker run` today.

**Why `-e SHELL=/bin/bash`?** Without it, JupyterLab's terminal defaults to `/bin/sh` (dash) — no arrow-key history, no Tab completion, no `history` builtin. This flag makes every terminal tab you open a real bash shell. (If a terminal is already open and stuck in dash, just type `exec bash` in it — no restart needed.)

Easier still: from the repo folder run **`docker compose up day1`** — both fixes are already baked into `docker-compose.yml`.

**Apple-Silicon Macs (M1/M2/M3):** add `--platform linux/amd64` to the commands.

Open **http://localhost:8888** in your browser. This is JupyterLab. We will **type the
commands ourselves in a terminal**, so open one:

**File → New → Terminal**

Everything you produce lands under `results/` in the JupyterLab file browser on the left
— open PNGs and tables by double-clicking them.

Prefer the command line only? On the host run
`docker run --rm -it -e JAVA_TOOL_OPTIONS=-XX:TieredStopAtLevel=1 leogiuffre/lncrna-mnps-workshop:1.0 bash` instead.

Prefer a fully guided click-through? The same steps are in the ready-made notebooks under `workshop/notebooks/` (run them top to bottom).

### 0.1 Create the output folders

The reads, the reference genome and the HISAT2 index are already inside the image under
`/home/student/`. Every step writes into a `results/` folder — create the folders once
with a single plain command (there is **no configuration file to source**, and we will
type every path in full so nothing is hidden):

```bash
mkdir -p /home/student/results/01_qc /home/student/results/02_trim /home/student/results/03_align /home/student/results/04_assembly /home/student/results/05_novel /home/student/results/06_lncrna /home/student/results/07_quant /home/student/results/08_de /home/student/results/09_figures /home/student/results/10_cis /home/student/results/11_ballgown
```

Our **6 libraries** are `control_1 control_2 control_3` (control) and
`trpet_1 trpet_2 trpet_3` (Tr-PET). We process them **one sample at a time**, then bring
them together only when a step genuinely needs all six.

### 0.2 Check the toolbox and tour the data

```bash
hisat2 --version | head -1

stringtie --version

CPC2 --help 2>&1 | head -1

featureCounts -v 2>&1 | tr -d '\n'; echo

Rscript -e 'cat(R.version.string); suppressMessages(library(edgeR)); cat(" | edgeR", as.character(packageVersion("edgeR")), "\n")'
```

```bash
column -t /home/student/data/reads/samples.tsv        # the 6 samples and their conditions

ls -lh /home/student/data/reads/*.fastq.gz            # the paired FASTQ (R1/R2)

samtools faidx /home/student/reference/chr4.fa

cat /home/student/reference/chr4.fa.fai               # Chromosome 4 length
```

**Question 1.** How many read pairs are in one library (hint: `zcat /home/student/data/reads/control_1_R1.fastq.gz | wc -l` then divide by 4)? In one sentence, why can we not discover *antisense* lncRNAs from an *unstranded* library?

**Laptop note.** Every command below uses **4** CPU threads. On a small laptop, change each `4` to `2` (i.e. `-threads 2`, `-p 2`, `-T 2`).

---

## Part 1 – Quality Control (FastQC + MultiQC)

Run FastQC on all 12 FASTQ files (R1 + R2 of the 6 samples), then summarise them with
MultiQC:

```bash
fastqc -t 4 -o /home/student/results/01_qc /home/student/data/reads/control_1_R1.fastq.gz /home/student/data/reads/control_1_R2.fastq.gz /home/student/data/reads/control_2_R1.fastq.gz /home/student/data/reads/control_2_R2.fastq.gz /home/student/data/reads/control_3_R1.fastq.gz /home/student/data/reads/control_3_R2.fastq.gz /home/student/data/reads/trpet_1_R1.fastq.gz /home/student/data/reads/trpet_1_R2.fastq.gz /home/student/data/reads/trpet_2_R1.fastq.gz /home/student/data/reads/trpet_2_R2.fastq.gz /home/student/data/reads/trpet_3_R1.fastq.gz /home/student/data/reads/trpet_3_R2.fastq.gz

multiqc -f -o /home/student/results/01_qc /home/student/results/01_qc
```

Open `results/01_qc/multiqc_report.html` in the JupyterLab file browser.

**Question 2.** Look at "Per base sequence quality" and "Adapter Content". Is there a 3′ quality drop? Which adapter is flagged? Why is some sequence duplication *expected* and not alarming in RNA-seq?

---

## Part 2 – Trimming (Trimmomatic)

We use **Trimmomatic** with the exact parameters from the paper. `ILLUMINACLIP` removes
TruSeq adapters; `SLIDINGWINDOW:4:18` trims once mean quality drops; `MINLEN:40` drops
reads too short to align uniquely. Paired-end mode writes a *paired* output (mate
survived) and an *unpaired* one — we keep the paired reads.

We run Trimmomatic **once per sample** — the same command six times, changing only the
sample name. Trimmomatic prints an `Input Read Pairs … Both Surviving` line to the screen
after each run (that is the number you need for Question 3).

```bash
# --- control_1 ---
trimmomatic PE -threads 4 -phred33 /home/student/data/reads/control_1_R1.fastq.gz /home/student/data/reads/control_1_R2.fastq.gz /home/student/results/02_trim/control_1_R1.paired.fq.gz /home/student/results/02_trim/control_1_R1.unpaired.fq.gz /home/student/results/02_trim/control_1_R2.paired.fq.gz /home/student/results/02_trim/control_1_R2.unpaired.fq.gz ILLUMINACLIP:/home/student/reference/adapters.fa:2:30:10 HEADCROP:1 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:18 MINLEN:40
```

```bash
# --- control_2 ---
trimmomatic PE -threads 4 -phred33 /home/student/data/reads/control_2_R1.fastq.gz /home/student/data/reads/control_2_R2.fastq.gz /home/student/results/02_trim/control_2_R1.paired.fq.gz /home/student/results/02_trim/control_2_R1.unpaired.fq.gz /home/student/results/02_trim/control_2_R2.paired.fq.gz /home/student/results/02_trim/control_2_R2.unpaired.fq.gz ILLUMINACLIP:/home/student/reference/adapters.fa:2:30:10 HEADCROP:1 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:18 MINLEN:40
```

```bash
# --- control_3 ---
trimmomatic PE -threads 4 -phred33 /home/student/data/reads/control_3_R1.fastq.gz /home/student/data/reads/control_3_R2.fastq.gz /home/student/results/02_trim/control_3_R1.paired.fq.gz /home/student/results/02_trim/control_3_R1.unpaired.fq.gz /home/student/results/02_trim/control_3_R2.paired.fq.gz /home/student/results/02_trim/control_3_R2.unpaired.fq.gz ILLUMINACLIP:/home/student/reference/adapters.fa:2:30:10 HEADCROP:1 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:18 MINLEN:40
```

```bash
# --- trpet_1 ---
trimmomatic PE -threads 4 -phred33 /home/student/data/reads/trpet_1_R1.fastq.gz /home/student/data/reads/trpet_1_R2.fastq.gz /home/student/results/02_trim/trpet_1_R1.paired.fq.gz /home/student/results/02_trim/trpet_1_R1.unpaired.fq.gz /home/student/results/02_trim/trpet_1_R2.paired.fq.gz /home/student/results/02_trim/trpet_1_R2.unpaired.fq.gz ILLUMINACLIP:/home/student/reference/adapters.fa:2:30:10 HEADCROP:1 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:18 MINLEN:40
```

```bash
# --- trpet_2 ---
trimmomatic PE -threads 4 -phred33 /home/student/data/reads/trpet_2_R1.fastq.gz /home/student/data/reads/trpet_2_R2.fastq.gz /home/student/results/02_trim/trpet_2_R1.paired.fq.gz /home/student/results/02_trim/trpet_2_R1.unpaired.fq.gz /home/student/results/02_trim/trpet_2_R2.paired.fq.gz /home/student/results/02_trim/trpet_2_R2.unpaired.fq.gz ILLUMINACLIP:/home/student/reference/adapters.fa:2:30:10 HEADCROP:1 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:18 MINLEN:40
```

```bash
# --- trpet_3 ---
trimmomatic PE -threads 4 -phred33 /home/student/data/reads/trpet_3_R1.fastq.gz /home/student/data/reads/trpet_3_R2.fastq.gz /home/student/results/02_trim/trpet_3_R1.paired.fq.gz /home/student/results/02_trim/trpet_3_R1.unpaired.fq.gz /home/student/results/02_trim/trpet_3_R2.paired.fq.gz /home/student/results/02_trim/trpet_3_R2.unpaired.fq.gz ILLUMINACLIP:/home/student/reference/adapters.fa:2:30:10 HEADCROP:1 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:18 MINLEN:40
```
**Question 3.** What fraction of read pairs survived trimming? Is it consistent across samples? Why do we keep only the **paired** output for a paired-end aligner?

---

## Part 3 – Splice-aware alignment (HISAT2)

RNA-seq reads span exon–exon junctions, so we use the splice-aware aligner **HISAT2**
against the Chr4 index. `--rna-strandness RF` records the transcribed strand;
`--dta` makes the alignments StringTie-friendly.

Align **one sample at a time**. Each command pipes HISAT2 straight into `samtools sort`
to make a sorted BAM, then indexes it. HISAT2 writes a small alignment summary next to it.

```bash
# --- control_1 ---
hisat2 -p 4 -x /home/student/reference/hisat2_index/chr4 -1 /home/student/results/02_trim/control_1_R1.paired.fq.gz -2 /home/student/results/02_trim/control_1_R2.paired.fq.gz --rna-strandness RF --dta --known-splicesite-infile /home/student/reference/chr4.ss.txt --new-summary --summary-file /home/student/results/03_align/control_1.hisat2.summary 2>/dev/null | samtools sort -@ 4 -o /home/student/results/03_align/control_1.sorted.bam -

samtools index /home/student/results/03_align/control_1.sorted.bam
```

```bash
# --- control_2 ---
hisat2 -p 4 -x /home/student/reference/hisat2_index/chr4 -1 /home/student/results/02_trim/control_2_R1.paired.fq.gz -2 /home/student/results/02_trim/control_2_R2.paired.fq.gz --rna-strandness RF --dta --known-splicesite-infile /home/student/reference/chr4.ss.txt --new-summary --summary-file /home/student/results/03_align/control_2.hisat2.summary 2>/dev/null | samtools sort -@ 4 -o /home/student/results/03_align/control_2.sorted.bam -

samtools index /home/student/results/03_align/control_2.sorted.bam
```

```bash
# --- control_3 ---
hisat2 -p 4 -x /home/student/reference/hisat2_index/chr4 -1 /home/student/results/02_trim/control_3_R1.paired.fq.gz -2 /home/student/results/02_trim/control_3_R2.paired.fq.gz --rna-strandness RF --dta --known-splicesite-infile /home/student/reference/chr4.ss.txt --new-summary --summary-file /home/student/results/03_align/control_3.hisat2.summary 2>/dev/null | samtools sort -@ 4 -o /home/student/results/03_align/control_3.sorted.bam -

samtools index /home/student/results/03_align/control_3.sorted.bam
```

```bash
# --- trpet_1 ---
hisat2 -p 4 -x /home/student/reference/hisat2_index/chr4 -1 /home/student/results/02_trim/trpet_1_R1.paired.fq.gz -2 /home/student/results/02_trim/trpet_1_R2.paired.fq.gz --rna-strandness RF --dta --known-splicesite-infile /home/student/reference/chr4.ss.txt --new-summary --summary-file /home/student/results/03_align/trpet_1.hisat2.summary 2>/dev/null | samtools sort -@ 4 -o /home/student/results/03_align/trpet_1.sorted.bam -

samtools index /home/student/results/03_align/trpet_1.sorted.bam
```

```bash
# --- trpet_2 ---
hisat2 -p 4 -x /home/student/reference/hisat2_index/chr4 -1 /home/student/results/02_trim/trpet_2_R1.paired.fq.gz -2 /home/student/results/02_trim/trpet_2_R2.paired.fq.gz --rna-strandness RF --dta --known-splicesite-infile /home/student/reference/chr4.ss.txt --new-summary --summary-file /home/student/results/03_align/trpet_2.hisat2.summary 2>/dev/null | samtools sort -@ 4 -o /home/student/results/03_align/trpet_2.sorted.bam -

samtools index /home/student/results/03_align/trpet_2.sorted.bam
```

```bash
# --- trpet_3 ---
hisat2 -p 4 -x /home/student/reference/hisat2_index/chr4 -1 /home/student/results/02_trim/trpet_3_R1.paired.fq.gz -2 /home/student/results/02_trim/trpet_3_R2.paired.fq.gz --rna-strandness RF --dta --known-splicesite-infile /home/student/reference/chr4.ss.txt --new-summary --summary-file /home/student/results/03_align/trpet_3.hisat2.summary 2>/dev/null | samtools sort -@ 4 -o /home/student/results/03_align/trpet_3.sorted.bam -

samtools index /home/student/results/03_align/trpet_3.sorted.bam
```

When all six are done, read the alignment rates and inspect one BAM:

```bash
grep 'Overall alignment rate' /home/student/results/03_align/*.hisat2.summary

samtools flagstat /home/student/results/03_align/control_1.sorted.bam
```

**Question 4.** What is the overall alignment rate (these libraries are Chr4-enriched, so it is high)? What does `--dta` change, and why does StringTie need it?

---

## Part 4 – Reference-guided transcript assembly (StringTie)

**StringTie** reconstructs the transcripts that were actually expressed. `-G` supplies
the known annotation as a *guide* (not a straitjacket): known transcripts are recovered
and genuinely new ones are still reported. We then **merge** the per-sample assemblies
into one non-redundant transcriptome so all samples are quantified against the same set.

First assemble transcripts **for each sample** (`--rf` = reverse-stranded; `-l` sets the
transcript name prefix):

```bash
# --- control_1 ---
stringtie /home/student/results/03_align/control_1.sorted.bam --rf -p 4 -G /home/student/reference/chr4.gtf -o /home/student/results/04_assembly/control_1.gtf -l control_1
```

```bash
# --- control_2 ---
stringtie /home/student/results/03_align/control_2.sorted.bam --rf -p 4 -G /home/student/reference/chr4.gtf -o /home/student/results/04_assembly/control_2.gtf -l control_2
```

```bash
# --- control_3 ---
stringtie /home/student/results/03_align/control_3.sorted.bam --rf -p 4 -G /home/student/reference/chr4.gtf -o /home/student/results/04_assembly/control_3.gtf -l control_3
```

```bash
# --- trpet_1 ---
stringtie /home/student/results/03_align/trpet_1.sorted.bam --rf -p 4 -G /home/student/reference/chr4.gtf -o /home/student/results/04_assembly/trpet_1.gtf -l trpet_1
```

```bash
# --- trpet_2 ---
stringtie /home/student/results/03_align/trpet_2.sorted.bam --rf -p 4 -G /home/student/reference/chr4.gtf -o /home/student/results/04_assembly/trpet_2.gtf -l trpet_2
```

```bash
# --- trpet_3 ---
stringtie /home/student/results/03_align/trpet_3.sorted.bam --rf -p 4 -G /home/student/reference/chr4.gtf -o /home/student/results/04_assembly/trpet_3.gtf -l trpet_3
```

**Now the first "bring everyone together" step:** merge the six per-sample assemblies into
one non-redundant transcriptome, so all samples are quantified against the same set. We
simply list the six GTFs after `--merge`:

```bash
stringtie --merge -p 4 -G /home/student/reference/chr4.gtf -o /home/student/results/04_assembly/merged.gtf /home/student/results/04_assembly/control_1.gtf /home/student/results/04_assembly/control_2.gtf /home/student/results/04_assembly/control_3.gtf /home/student/results/04_assembly/trpet_1.gtf /home/student/results/04_assembly/trpet_2.gtf /home/student/results/04_assembly/trpet_3.gtf

echo "merged transcripts: $(grep -c $'\ttranscript\t' /home/student/results/04_assembly/merged.gtf)"

echo "known transcripts : $(grep -c $'\ttranscript\t' /home/student/reference/chr4.gtf)"
```

**Question 5.** Does the merged assembly contain **more** transcripts than the annotation? Where do the extra ones come from? What would happen to novel-transcript discovery if we ran StringTie **without** `-G`?

---

## Part 5 – Novel lncRNA identification (gffcompare + CPC2)

`gffcompare` tags every assembled transcript with a **class code** relative to the
reference. We keep the ones that are **not** already annotated and do not overlap a known
gene on the same strand: **`u`** (intergenic), **`i`** (intronic), **`x`** (antisense) —
and require length **≥ 200 nt**.

```bash
gffcompare -r /home/student/reference/chr4.gtf -o /home/student/results/05_novel/gffcmp /home/student/results/04_assembly/merged.gtf

echo "== class-code distribution =="

grep $'\ttranscript\t' /home/student/results/05_novel/gffcmp.annotated.gtf | grep -oE 'class_code "."' | sort | uniq -c | sort -rn
```

Select the candidates (helper script) and extract their spliced sequences:

```bash
python3 /home/student/workshop/lib/gtf_select.py /home/student/results/05_novel/gffcmp.annotated.gtf --classes u,i,x --min-len 200 --out-gtf /home/student/results/05_novel/candidates.gtf --out-tsv /home/student/results/05_novel/candidates.tsv

gffread -w /home/student/results/05_novel/candidates.fa -g /home/student/reference/chr4.fa /home/student/results/05_novel/candidates.gtf

echo "candidate novel transcripts: $(( $(wc -l < /home/student/results/05_novel/candidates.tsv) - 1 ))"
```

A novel transcript is a **lncRNA** only if it is unlikely to be translated. Score coding
potential with **CPC2** and keep the **non-coding** transcripts:

```bash
CPC2 -i /home/student/results/05_novel/candidates.fa -o /home/student/results/06_lncrna/cpc2 >/dev/null

awk -F'\t' 'NR>1{print $1"\t"$NF}' /home/student/results/06_lncrna/cpc2.txt > /home/student/results/06_lncrna/cpc2.labels.tsv

awk -F'\t' '$2=="noncoding"{print $1}' /home/student/results/06_lncrna/cpc2.labels.tsv > /home/student/results/06_lncrna/lncRNA.ids.txt

echo "noncoding (lncRNA): $(wc -l < /home/student/results/06_lncrna/lncRNA.ids.txt)   coding: $(grep -c coding /home/student/results/06_lncrna/cpc2.labels.tsv)"

grep -F -w -f /home/student/results/06_lncrna/lncRNA.ids.txt /home/student/results/05_novel/candidates.gtf > /home/student/results/06_lncrna/novel_lncRNA.gtf
```

Finally, build the **augmented annotation** (known genes + novel lncRNAs) and a table
labelling every gene as `mRNA` / `lncRNA` / `other`:

```bash
python3 /home/student/workshop/lib/build_augmented.py --ref-gtf /home/student/reference/chr4.gtf --lnc-gtf /home/student/results/06_lncrna/novel_lncRNA.gtf --out-gtf /home/student/results/06_lncrna/augmented.gtf --out-tsv /home/student/results/06_lncrna/gene_biotype.tsv

echo "== gene classes in augmented annotation =="

tail -n +2 /home/student/results/06_lncrna/gene_biotype.tsv | cut -f3 | sort | uniq -c
```

**Question 6.** Which class code dominates your candidates — intergenic (`u`), intronic (`i`) or antisense (`x`)? What fraction of candidates did CPC2 call *coding* and discard? Why is a coding-potential filter necessary — aren't unannotated transcripts non-coding by definition?

---

## Part 6 – Gene-level quantification (featureCounts)

Count fragments per gene (mRNAs **and** lncRNAs) against the augmented annotation. Note
`-s 2`: **reverse-stranded** counting, so an antisense lncRNA is not confused with its
sense neighbour.

This is the **count-matrix "bring everyone together" step**: one `featureCounts` command
takes **all six sorted BAMs at once** and produces a single table with six count columns
(`-s 2` = reverse-stranded counting):

```bash
featureCounts -T 4 -p --countReadPairs -s 2 -t exon -g gene_id -a /home/student/results/06_lncrna/augmented.gtf -o /home/student/results/07_quant/counts.txt /home/student/results/03_align/control_1.sorted.bam /home/student/results/03_align/control_2.sorted.bam /home/student/results/03_align/control_3.sorted.bam /home/student/results/03_align/trpet_1.sorted.bam /home/student/results/03_align/trpet_2.sorted.bam /home/student/results/03_align/trpet_3.sorted.bam

# tidy the column names (strip paths and .sorted.bam)
sed -i '2s#[^\t]*/##g; 2s#\.sorted\.bam##g' /home/student/results/07_quant/counts.txt

column -t /home/student/results/07_quant/counts.txt.summary

grep -v '^#' /home/student/results/07_quant/counts.txt | cut -f1,7- | grep '^MSTRG' | head -5 | column -t
```

**Question 7.** What percentage of fragments were **assigned**? Do the novel lncRNAs (`MSTRG.*`) carry non-trivial counts? (They must, or the differential test will filter them out.)

---

## Part 7 – Differential expression with edgeR (R)

Now switch to **R**. Start it in the same terminal:

```bash
R
```

Inside R, load the paths and the count matrix. edgeR models raw counts with a
negative-binomial distribution; we filter lowly-expressed genes, normalise with TMM, and
test **Tr-PET vs control** with a Likelihood-Ratio Test. Significance: **|log₂FC| ≥ 0.5
and FDR < 0.05** (paper thresholds).

```r
source("/home/student/workshop/scripts/paths.R")   # sets QUANT_DIR, LNC_DIR, DE_DIR, samples, thresholds
suppressMessages(library(edgeR))

fc <- read.delim(file.path(QUANT_DIR, "counts.txt"), comment.char = "#", check.names = FALSE)
counts <- as.matrix(fc[, 7:ncol(fc)]); rownames(counts) <- fc$Geneid
counts <- counts[, samples$sample]                    # order columns as the sample sheet
group  <- factor(samples$condition, levels = c("control", "TrPET"))

bt <- read.delim(file.path(LNC_DIR, "gene_biotype.tsv"), stringsAsFactors = FALSE)
grp_of <- setNames(bt$group, bt$gene_id)

y <- DGEList(counts = counts, group = group)
keep <- rowSums(cpm(y) > 1) >= 2                      # keep genes CPM>1 in >=2 libraries
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)                               # TMM normalisation

design <- model.matrix(~ group)
y   <- estimateDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef = 2)                          # TrPET vs control

res <- topTags(lrt, n = Inf, sort.by = "PValue")$table
res$gene_id     <- rownames(res)
res$rna_class   <- ifelse(res$gene_id %in% names(grp_of), grp_of[res$gene_id], "other")
res$significant <- abs(res$logFC) >= LFC_THR & res$FDR < FDR_THR
res$direction   <- ifelse(!res$significant, "ns", ifelse(res$logFC > 0, "up", "down"))

# save for the plotting / cis-target steps
write.table(res, file.path(DE_DIR, "edger_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
logcpm <- cpm(y, log = TRUE, prior.count = 2)
saveRDS(list(res = res, logcpm = logcpm, group = group), file.path(DE_DIR, "de_objects.rds"))

# how many DE genes, by RNA class and direction?
sig <- res[res$significant, ]
addmargins(table(rna_class = sig$rna_class, direction = sig$direction))
```

**Question 8.** How many **lncRNAs** and how many **mRNAs** are differentially expressed, and in which direction? Why filter lowly-expressed genes *before* testing?

---

## Part 8 – Visualisation (PCA / MA / volcano / heatmap)

Still in R:

```r
suppressMessages({library(ggplot2); library(pheatmap)})
obj <- readRDS(file.path(DE_DIR, "de_objects.rds"))
res <- obj$res; logcpm <- obj$logcpm; group <- obj$group
res$sig_lab <- ifelse(res$significant, res$rna_class, "n.s.")
pal <- c(mRNA="#1b7837", lncRNA="#762a83", other="#999999", n.s.="#cccccc")

# PCA
pca <- prcomp(t(logcpm)); ve <- round(100*pca$sdev^2/sum(pca$sdev^2), 1)
dp <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], sample=colnames(logcpm), condition=group)
ggplot(dp, aes(PC1, PC2, colour=condition, label=sample)) +
  geom_point(size=4) + geom_text(vjust=-1, size=3) + theme_bw() +
  labs(x=sprintf("PC1 (%.1f%%)", ve[1]), y=sprintf("PC2 (%.1f%%)", ve[2]), title="PCA (log-CPM)")
ggsave(file.path(FIG_DIR, "PCA.png"), width=6, height=5, dpi=120)

# Volcano
res$neglog10FDR <- -log10(res$FDR)
ggplot(res, aes(logFC, neglog10FDR, colour=sig_lab)) +
  geom_point(size=0.8, alpha=0.6) +
  geom_vline(xintercept=c(-LFC_THR, LFC_THR), linetype=2) +
  geom_hline(yintercept=-log10(FDR_THR), linetype=2) +
  scale_colour_manual(values=pal) + theme_bw() +
  labs(title="Volcano (Tr-PET vs control)", x="log2 fold change", y="-log10 FDR", colour=NULL)
ggsave(file.path(FIG_DIR, "volcano.png"), width=6.5, height=5, dpi=120)

# Heatmap of the top DE genes
sig <- res[res$significant, ]; top <- head(sig[order(sig$FDR), "gene_id"], 40)
if (length(top) >= 2) {
  mat <- t(scale(t(logcpm[top, , drop=FALSE])))
  ac <- data.frame(condition=group); rownames(ac) <- colnames(mat)
  ar <- data.frame(class=res$rna_class[match(top, res$gene_id)]); rownames(ar) <- top
  pheatmap(mat, annotation_col=ac, annotation_row=ar, fontsize_row=6,
           main="Top DE genes (row z-score log-CPM)",
           filename=file.path(FIG_DIR, "heatmap_top_DE.png"), width=7, height=8)
}
```

Open `results/09_figures/*.png` in the file browser.

**Question 9.** Do control and Tr-PET separate on **PC1**? On the volcano, are lncRNAs (purple) over-represented among the strongly regulated genes?

---

## Part 9 – lncRNA *cis*-targets & GO enrichment (R)

A central finding of the paper is that these lncRNAs tend to regulate their
**neighbouring protein-coding genes**. Find the nearest protein-coding gene (±10 kb) of
each DE lncRNA, then ask which biological processes those neighbours are enriched for
(clusterProfiler + org.At.tair.db — an offline stand-in for the paper's gProfiler).

```r
suppressMessages({library(clusterProfiler); library(org.At.tair.db); library(ggplot2)})
WINDOW <- 10000
obj <- readRDS(file.path(DE_DIR, "de_objects.rds")); res <- obj$res
bt  <- read.delim(file.path(LNC_DIR, "gene_biotype.tsv")); grp <- setNames(bt$group, bt$gene_id)

read_gene_coords <- function(gtf) {
  x <- read.delim(gtf, header=FALSE, comment.char="#", quote="", stringsAsFactors=FALSE)
  x <- x[x$V3=="exon", ]; gid <- sub('.*gene_id "([^"]+)".*', "\\1", x$V9)
  data.frame(gene_id=tapply(gid,gid,`[`,1), chr=tapply(x$V1,gid,`[`,1),
             start=tapply(x$V4,gid,min), end=tapply(x$V5,gid,max),
             strand=tapply(x$V7,gid,`[`,1), stringsAsFactors=FALSE)
}
coords <- read_gene_coords(file.path(LNC_DIR, "augmented.gtf"))
de_lnc <- res$gene_id[res$significant & res$rna_class=="lncRNA"]
pc_all <- names(grp)[grp=="mRNA"]
cat(sprintf("DE lncRNAs: %d   protein-coding genes: %d\n", length(de_lnc), length(pc_all)))

to_bed <- function(ids){ d<-coords[coords$gene_id %in% ids & !is.na(coords$chr),]
  data.frame(chr=d$chr, start=pmax(0,d$start-1), end=d$end, name=d$gene_id, score=0, strand=d$strand) }
lb<-file.path(CIS_DIR,"de_lnc.bed"); pb<-file.path(CIS_DIR,"pc.bed")
write.table(to_bed(de_lnc), lb, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(to_bed(pc_all), pb, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
system(paste("sort -k1,1 -k2,2n", lb, "-o", lb)); system(paste("sort -k1,1 -k2,2n", pb, "-o", pb))
system(paste("bedtools closest -a", lb, "-b", pb, "-d -t first >", file.path(CIS_DIR,"closest.tsv")))

cl <- read.delim(file.path(CIS_DIR,"closest.tsv"), header=FALSE)
cl <- cl[cl$V13>=0 & cl$V13<=WINDOW, ]
cis <- data.frame(lncRNA=cl$V4, target_gene=cl$V10, distance=cl$V13)
write.table(cis, file.path(CIS_DIR,"cis_targets.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
cat(sprintf("cis lncRNA-target pairs within %d bp: %d\n", WINDOW, nrow(cis))); head(cis)

ego <- enrichGO(gene=unique(cis$target_gene), OrgDb=org.At.tair.db, keyType="TAIR",
                ont="BP", pvalueCutoff=0.05, qvalueCutoff=0.1)
if (!is.null(ego) && nrow(as.data.frame(ego))>0) {
  write.table(as.data.frame(ego), file.path(CIS_DIR,"go_enrichment.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
  ggsave(file.path(CIS_DIR,"go_dotplot.png"), dotplot(ego, showCategory=15), width=8, height=6, dpi=120)
  head(as.data.frame(ego)[, c("ID","Description","p.adjust","Count")], 10)
} else cat("No significant GO enrichment for this (small) Chr4 target set.\n")
```

Type `quit(save = "no")` to leave R.

**Question 10.** How many DE lncRNAs have a protein-coding neighbour within 10 kb? Which GO biological processes surface — do any relate to stress, transport or cell-wall responses? This is a *positional* prediction: what experiment would *test* that a lncRNA truly regulates its neighbour?

---

## Bonus – Isoform-level analysis (Ballgown)

Everything above was gene-level. To look at individual isoforms, re-estimate abundances
with StringTie `-e -B` and test with Ballgown:

```bash
# back at the shell (leave R first with quit(save="no"))
bash /home/student/workshop/scripts/11_ballgown.sh

Rscript /home/student/workshop/scripts/11_ballgown.R
```

---

## Shortcut – run the whole pipeline at once

If you fall behind, the entire Day 1 pipeline (Parts 1–9, plus the bonus) can be run
non-interactively:

```bash
bash /home/student/workshop/scripts/run_all.sh --with-ballgown
```

---

## Discussion Questions

1. Micro-nanoplastics are an osmotic-type abiotic stress. Would you expect stress-induced
   lncRNAs to act mostly in *cis* (on neighbours) or in *trans*? What evidence in your
   results supports either?
2. Day 1 (lncRNA/mRNA) and Day 2 (DNA methylation) probe the same stress-signalling
   biology. How could a *cis*-acting lncRNA and DNA methylation at the same locus be
   mechanistically linked (think RdDM)?
3. We restricted the analysis to Chromosome 4 and subsampled the reads. What biases does
   this introduce, and how would you scale to the whole genome on a server?
4. gffcompare class `x` (antisense) depends entirely on correct strandedness. What would
   happen to your antisense lncRNA calls if the library strand were mis-specified?

---

## End of Day 1 – Well Done!

Your tables and figures are under `results/` (visible in the JupyterLab file browser).
To keep them after the container stops, download the `results/` folder from the file
browser, or re-run with a bind-mount:

```bash
docker run --rm -p 8888:8888 -e JAVA_TOOL_OPTIONS=-XX:TieredStopAtLevel=1 -e SHELL=/bin/bash -v "$PWD/day1_results":/home/student/results leogiuffre/lncrna-mnps-workshop:1.0
```

Tomorrow (Day 2) we decode the **methylome** of the same organism with Whole-Genome
Bisulfite Sequencing — the second molecular layer of the same stress story.
