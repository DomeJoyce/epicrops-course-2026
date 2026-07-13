#!/usr/bin/env bash
# =============================================================================
# run_day2.sh  –  Day 2 Master Pipeline
# WGBS Analysis: QC → Trim → Bismark Align → Extract → methylKit (R)
#
# The Epi-Code – Florence Training School 2026
# =============================================================================

set -euo pipefail

THREADS="${1:-4}"
RAW_DIR="${DATA_DIR:-/course/data}/raw/wgbs"
REF_DIR="${REF_DIR:-/course/data/reference}"
OUT_DIR="${RESULTS_DIR:-/course/results}/day2_methylation"

echo ""
echo "╔════════════════════════════════════════════════════════════╗"
echo "║  Day 2 – Whole-Genome Bisulfite Sequencing Pipeline        ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# STEP 1 – Quality Control
# ─────────────────────────────────────────────────────────────────────────────
echo "[STEP 1/5] Quality Control with FastQC..."
mkdir -p "${OUT_DIR}/01_qc"

fastqc \
    --threads "${THREADS}" \
    --outdir "${OUT_DIR}/01_qc" \
    "${RAW_DIR}"/*.fastq.gz \
    2>&1 | tee "${OUT_DIR}/01_qc/fastqc.log"

multiqc \
    "${OUT_DIR}/01_qc" \
    --outdir "${OUT_DIR}/01_qc/multiqc" \
    --filename "Day2_raw_QC" \
    --quiet

echo "  → QC reports: ${OUT_DIR}/01_qc/multiqc/"

# ─────────────────────────────────────────────────────────────────────────────
# STEP 2 – Adapter Trimming (Trim Galore WGBS mode)
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "[STEP 2/5] Adapter and quality trimming with Trim Galore (WGBS mode)..."
mkdir -p "${OUT_DIR}/02_trimmed"

# Read sample sheet to get R1/R2 pairs
tail -n +2 "${RAW_DIR}/samplesheet.tsv" | while IFS=$'\t' read -r SAMPLE COND REP R1 R2; do
    echo "  Trimming: ${SAMPLE} (${COND})"
    trim_galore \
        --paired \
        --cores "${THREADS}" \
        --quality 20 \
        --stringency 5 \
        --dont_gzip \
        --clip_r1 6 \
        --clip_r2 6 \
        --three_prime_clip_R1 6 \
        --three_prime_clip_R2 6 \
        --output_dir "${OUT_DIR}/02_trimmed" \
        "${RAW_DIR}/${R1}" \
        "${RAW_DIR}/${R2}" \
        2>&1 | tee "${OUT_DIR}/02_trimmed/${SAMPLE}_trimming.log"
done

multiqc \
    "${OUT_DIR}/02_trimmed" \
    --outdir "${OUT_DIR}/02_trimmed/multiqc" \
    --filename "Day2_trimmed_QC" \
    --quiet

echo "  → Trimmed reads: ${OUT_DIR}/02_trimmed/"

# ─────────────────────────────────────────────────────────────────────────────
# STEP 3 – Bismark Alignment (bisulfite-aware)
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "[STEP 3/5] Bisulfite alignment with Bismark (Bowtie2 mode)..."
mkdir -p "${OUT_DIR}/03_bismark_aligned"

BISMARK_IDX="${REF_DIR}/indices/bismark"

tail -n +2 "${RAW_DIR}/samplesheet.tsv" | while IFS=$'\t' read -r SAMPLE COND REP R1 R2; do
    # Trim Galore output naming convention
    R1_TRIMMED="${OUT_DIR}/02_trimmed/$(basename ${R1} .fastq.gz)_val_1.fq"
    R2_TRIMMED="${OUT_DIR}/02_trimmed/$(basename ${R2} .fastq.gz)_val_2.fq"

    echo "  Aligning: ${SAMPLE}"
    # Bismark paired-end uses -1/-2 (not --paired-end); --multicore sets parallelism.
    # With --basename the output BAM is ${SAMPLE}_pe.bam.
    bismark \
        --genome "${BISMARK_IDX}" \
        -1 "${R1_TRIMMED}" \
        -2 "${R2_TRIMMED}" \
        --bowtie2 \
        --multicore "${THREADS}" \
        --output_dir "${OUT_DIR}/03_bismark_aligned" \
        --basename "${SAMPLE}" \
        2>&1 | tee "${OUT_DIR}/03_bismark_aligned/${SAMPLE}_bismark.log"

    echo "  Deduplicating: ${SAMPLE}"
    deduplicate_bismark \
        --paired \
        --output_dir "${OUT_DIR}/03_bismark_aligned" \
        "${OUT_DIR}/03_bismark_aligned/${SAMPLE}_pe.bam" \
        2>&1 | tee "${OUT_DIR}/03_bismark_aligned/${SAMPLE}_dedup.log"

    echo "  → BAM: ${OUT_DIR}/03_bismark_aligned/${SAMPLE}_pe.deduplicated.bam"
done

echo "  → Alignment stats: ${OUT_DIR}/03_bismark_aligned/*.log"

# ─────────────────────────────────────────────────────────────────────────────
# STEP 4 – Methylation Extraction
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "[STEP 4/5] Extracting methylation calls with Bismark..."
mkdir -p "${OUT_DIR}/04_methylation"

for BAM in "${OUT_DIR}/03_bismark_aligned"/*.deduplicated.bam; do
    SAMPLE=$(basename "${BAM}" _pe.deduplicated.bam)
    echo "  Extracting: ${SAMPLE}"

    # NOTE: bismark_methylation_extractor has NO --basename option (that is an
    # aligner-only flag); outputs are named after the input BAM. --CX makes the
    # genome-wide cytosine report cover CpG + CHG + CHH. The .bedGraph/.cov files
    # for IGV are produced automatically by --comprehensive.
    bismark_methylation_extractor \
        --paired-end \
        --comprehensive \
        --cytosine_report \
        --CX \
        --genome_folder "${BISMARK_IDX}" \
        --parallel "${THREADS}" \
        --output_dir "${OUT_DIR}/04_methylation" \
        "${BAM}" \
        2>&1 | tee "${OUT_DIR}/04_methylation/${SAMPLE}_extraction.log"

    # Rename the CX report to a clean per-sample name for methylKit
    # (from <SAMPLE>_pe.deduplicated.CX_report.txt -> <SAMPLE>.CX_report.txt)
    report=$(ls "${OUT_DIR}/04_methylation/${SAMPLE}"*.CX_report.txt 2>/dev/null | grep -v "^${OUT_DIR}/04_methylation/${SAMPLE}\.CX_report\.txt$" | head -1 || true)
    if [[ -n "${report}" ]]; then
        mv -f "${report}" "${OUT_DIR}/04_methylation/${SAMPLE}.CX_report.txt"
    fi
done

echo "  → Methylation reports: ${OUT_DIR}/04_methylation/"

# ─────────────────────────────────────────────────────────────────────────────
# STEP 5 – Differential Methylation Analysis with methylKit (R)
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "[STEP 5/5] Differential methylation analysis with methylKit (R)..."
mkdir -p "${OUT_DIR}/05_diffmeth"

Rscript "${SCRIPTS_DIR:-/course/scripts}/day2_methylkit.R" \
    "${OUT_DIR}/04_methylation" \
    "${RAW_DIR}/samplesheet.tsv" \
    "${OUT_DIR}/05_diffmeth" \
    2>&1 | tee "${OUT_DIR}/05_diffmeth/methylkit.log"

echo ""
echo "╔════════════════════════════════════════════════════════════╗"
echo "║  Day 2 Pipeline COMPLETE                                    ║"
echo "╠════════════════════════════════════════════════════════════╣"
echo "║  QC reports  : ${OUT_DIR}/01_qc/multiqc/    "
echo "║  Bismark BAMs: ${OUT_DIR}/03_bismark_aligned/"
echo "║  Methylation : ${OUT_DIR}/04_methylation/   "
echo "║  DM results  : ${OUT_DIR}/05_diffmeth/      "
echo "╚════════════════════════════════════════════════════════════╝"
echo ""
