# EpiCrops Course — Student Handbook
## lncRNA transcriptomics & DNA methylation in plants
### Florence Training School, 15–16 July 2026 · University of Messina

Welcome! Over two days you will run **two complete bioinformatics analyses** on
*Arabidopsis thaliana* — start to finish, on your own laptop, entirely inside Docker
(no software to install by hand). Everything is subsampled to **Chromosome 4** so each
pipeline finishes in minutes.

| Day | You will… | Data | Walkthrough |
|-----|-----------|------|-------------|
| **1** | Identify **novel lncRNAs** and find genes that change expression under micro-nanoplastic stress (RNA-seq → HISAT2 → StringTie → CPC2 → edgeR) | control vs **Tr-PET** (3+3) | [exercises/Day1_lncRNA_mRNA_Pipeline.md](exercises/Day1_lncRNA_mRNA_Pipeline.md) |
| **2** | Map **DNA methylation** and find regions that change under salt stress (WGBS → Bismark → methylKit) | control vs **salt** (3+3) | [exercises/Day2_DNA_Methylation_Analysis.md](exercises/Day2_DNA_Methylation_Analysis.md) |

**First time working on a terminal?** Read [BASH_tutorial.md](BASH_tutorial.md) — it explains 
*the basic bash command commonly used on a linux terminal*, so you are able to do the basic command into your machine.

**New to the commands?** Read [COMMANDS_EXPLAINED.md](COMMANDS_EXPLAINED.md) — it explains
*what every command does and why*, so the pipelines make sense instead of being magic.

---

## Before the course (do this at home!)

Please complete **[SETUP_BEFORE_CLASS.md](SETUP_BEFORE_CLASS.md)** *before* you arrive:
install Docker, get the two images, and check they start. The downloads are large — you do
**not** want to be doing them on classroom Wi-Fi.

Minimum laptop: **Docker Desktop**, **8 GB RAM**, **4 cores**, **~20 GB free disk**.

---

## Get the code

```bash
git clone https://github.com/DomeJoyce/epicrops-course-2026.git
cd epicrops-course-2026
```

---

## Day 1 — lncRNA / mRNA (self-contained image)

Everything for Day 1 (tools, reference **and** the 6 RNA-seq libraries) is inside one
image — nothing to download or set up.

**Start the interactive environment (recommended):**
```bash
docker compose up day1        # then open http://localhost:8888  (JupyterLab)
```
In JupyterLab: **File → New → Terminal**, and type the commands from the
[Day 1 walkthrough](exercises/Day1_lncRNA_mRNA_Pipeline.md).

> The compose file already sets two fixes: `JAVA_TOOL_OPTIONS` (keeps Trimmomatic from
> crashing) and `SHELL=/bin/bash` (keeps the JupyterLab terminal from defaulting to a
> shell with no arrow-key history or Tab completion) — so prefer `docker compose up day1`.
> If you launch the image with a raw `docker run`, add
> `-e JAVA_TOOL_OPTIONS=-XX:TieredStopAtLevel=1 -e SHELL=/bin/bash` (the walkthrough shows
> this). Apple-Silicon Macs: also add `--platform linux/amd64`.

**Or run the whole pipeline in one go (no typing):**
```bash
bash scripts/run_day1.sh      # QC → trim → HISAT2 → StringTie → CPC2 → featureCounts → edgeR → …
```
Results land in `results/day1_lncRNA/` on your laptop (~1.5 GB; whole run ≈ 8–15 min).

---

## Day 2 — WGBS / DNA methylation (course image)

Day 2 uses a course image built from this repository.

**1. Get the image** — pull the pre-built image (fast, recommended):
```bash
docker pull djoyce86/epi-code-practical:2026
docker tag djoyce86/epi-code-practical:2026 epi-code-practical:latest
```
> *No internet at the venue, or you prefer building it yourself? It also works from
> source (~20–40 min the first time):* `docker compose build course`

**2. One-time setup + data (inside the container):**
```bash
docker compose run --rm course          # opens a shell in the container
bash $SCRIPTS_DIR/validate_env.sh        # should end with "ALL GOOD"
bash $SCRIPTS_DIR/00_setup_reference.sh 4   # build Chr4 indices (~1 min here; a few min on a laptop)
bash $DATA_DIR/download_data.sh          # download + subsample the 6 WGBS libraries
```
> *Prefer the browser, like Day 1?* Run **`docker compose up day2`** instead and open
> **http://localhost:8888/** → *File → New → Terminal*, then run the same commands there.
> (Stop Day 1 first — both use port 8888: `docker rm -f epi_code_day1`.)

**3. Run the analysis** — follow the
[Day 2 walkthrough](exercises/Day2_DNA_Methylation_Analysis.md) step by step, **or** run it
all at once:
```bash
bash $SCRIPTS_DIR/run_day2.sh 4          # use 2 instead of 4 on a small laptop
```
Results land in `results/day2_methylation/` on your laptop.

---

## Where are my results?

Both days write to `results/` **on your laptop** (bind-mounted), so they survive after the
container stops:
- `results/day1_lncRNA/` — QC, alignments, assembly, lncRNA set, DE tables, figures
- `results/day2_methylation/` — QC, methylation reports, DMC/DMR tables, figures

Open the `.png`/`.html` files by double-clicking them (JupyterLab file browser or your OS).

---

## Troubleshooting

| Symptom | Fix |
|---------|-----|
| Day 1 Trimmomatic prints "A fatal error … Java Runtime" | Use `docker compose up day1`, or add `-e JAVA_TOOL_OPTIONS=-XX:TieredStopAtLevel=1` to your `docker run`. |
| Day 1 terminal: no arrow-key history, no Tab completion, `history: not found` | You're in `dash`, not `bash` (JupyterLab's terminal has no `$SHELL` set by default). Immediate fix in an open terminal: type `exec bash`. Permanent fix: use `docker compose up day1` (already sets `SHELL=/bin/bash`), or add `-e SHELL=/bin/bash` to your `docker run`. |
| `no space left on device` | Free disk; you need ~20 GB. `docker system prune` reclaims space. |
| Day 2 tool "not found" | Run `bash $SCRIPTS_DIR/validate_env.sh` — it should end with "ALL GOOD". |
| `bowtie2 … Failed to launch x86-64-v3 version` | Harmless — it falls back automatically. |
| Apple-Silicon Mac, Day 1 won't start | Add `--platform linux/amd64`. |
| Slow laptop | Pass fewer threads: `run_day2.sh 2`; on Day 1 change each `4` to `2` in the commands. |

---

*The Epi-Code / EpiCrops — Florence 2026 · Domenico Giosa & Letterio Giuffrè · University of Messina, ChiBioFarAm Dept.*
