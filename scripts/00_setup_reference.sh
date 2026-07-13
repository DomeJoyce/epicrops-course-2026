#!/usr/bin/env bash
# =============================================================================
# 00_setup_reference.sh
# Download and index the reference for DAY 2 (WGBS / methylation).
# Reference: Arabidopsis thaliana TAIR10, Chromosome 4 only (~18 Mb) — keeps the
# computational load suitable for student laptops.
#
# (DAY 1 – lncRNA/mRNA – needs no reference setup here: its Chr4 genome, GTF and
#  HISAT2 index are already bundled inside the Day 1 image.)
#
# Runtime estimate: ~10–15 minutes (CPU-dependent, run once before Day 2)
# =============================================================================

set -euo pipefail

REF_DIR="${REF_DIR:-/course/data/reference}"
GENOME_DIR="${REF_DIR}/genome"
ANNOT_DIR="${REF_DIR}/annotation"
IDX_BOWTIE2="${REF_DIR}/indices/bowtie2"
IDX_BISMARK="${REF_DIR}/indices/bismark"
THREADS="${1:-4}"

TAIR10_GENOME_URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.chromosome.4.fa.gz"
TAIR10_GFF_URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.59.chromosome.4.gff3.gz"

echo "============================================================"
echo " Setting up Day 2 reference data"
echo " Arabidopsis thaliana TAIR10 – Chromosome 4"
echo " Threads: ${THREADS}"
echo "============================================================"

# ── 1. Reference genome ──────────────────────────────────────────────────────
echo ""
echo "[1/4] Downloading TAIR10 Chr4 genome..."
mkdir -p "${GENOME_DIR}"

if [[ ! -f "${GENOME_DIR}/ath_chr4.fa" ]]; then
    wget -q --show-progress -O "${GENOME_DIR}/ath_chr4.fa.gz" "${TAIR10_GENOME_URL}"
    gunzip -k "${GENOME_DIR}/ath_chr4.fa.gz"
    # Rename chromosome header to a simple "Chr4" for tool compatibility
    sed -i 's/>4 .*/>Chr4/' "${GENOME_DIR}/ath_chr4.fa"
    echo "  → Genome: ${GENOME_DIR}/ath_chr4.fa"
else
    echo "  → Genome already present. Skipping."
fi

# ── 2. Gene annotation ───────────────────────────────────────────────────────
echo ""
echo "[2/4] Downloading TAIR10 Chr4 GFF3 annotation..."
mkdir -p "${ANNOT_DIR}"

if [[ ! -f "${ANNOT_DIR}/ath_chr4.gff3" ]]; then
    wget -q --show-progress -O "${ANNOT_DIR}/ath_chr4.gff3.gz" "${TAIR10_GFF_URL}"
    gunzip -k "${ANNOT_DIR}/ath_chr4.gff3.gz"
    # Patch chromosome name to match the genome file
    sed -i 's/^4\t/Chr4\t/' "${ANNOT_DIR}/ath_chr4.gff3"
    echo "  → Annotation: ${ANNOT_DIR}/ath_chr4.gff3"
else
    echo "  → Annotation already present. Skipping."
fi

# ── 3. Bowtie2 index (used by Bismark for bisulfite mapping) ─────────────────
echo ""
echo "[3/4] Building Bowtie2 index for Chr4..."
mkdir -p "${IDX_BOWTIE2}"

if [[ ! -f "${IDX_BOWTIE2}/ath_chr4.1.bt2" ]]; then
    bowtie2-build \
        --threads "${THREADS}" \
        "${GENOME_DIR}/ath_chr4.fa" \
        "${IDX_BOWTIE2}/ath_chr4" \
        > "${IDX_BOWTIE2}/bowtie2_build.log" 2>&1
    echo "  → Bowtie2 index: ${IDX_BOWTIE2}/ath_chr4.*"
else
    echo "  → Bowtie2 index already present. Skipping."
fi

# ── 4. Bismark index (for WGBS bisulfite alignment) ──────────────────────────
echo ""
echo "[4/4] Building Bismark bisulfite genome index for Chr4..."
mkdir -p "${IDX_BISMARK}"

if [[ ! -f "${IDX_BISMARK}/Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa" ]]; then
    cp "${GENOME_DIR}/ath_chr4.fa" "${IDX_BISMARK}/"
    bismark_genome_preparation \
        --bowtie2 \
        --parallel "${THREADS}" \
        "${IDX_BISMARK}/" \
        > "${IDX_BISMARK}/bismark_prepare.log" 2>&1
    # Keep the FASTA in the index folder: Bismark's aligner and the cytosine-report
    # extractor both read the reference genome from --genome / --genome_folder.
    echo "  → Bismark index: ${IDX_BISMARK}/Bisulfite_Genome/ (+ ath_chr4.fa)"
else
    echo "  → Bismark index already present. Skipping."
fi

# ── Summary ──────────────────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo " Day 2 reference setup COMPLETE"
echo "============================================================"
echo ""
echo "  Genome       : ${GENOME_DIR}/ath_chr4.fa"
echo "  Annotation   : ${ANNOT_DIR}/ath_chr4.gff3"
echo "  Bowtie2 idx  : ${IDX_BOWTIE2}/ath_chr4.*"
echo "  Bismark idx  : ${IDX_BISMARK}/Bisulfite_Genome/"
echo ""
echo "  You are now ready to start Day 2!"
echo "  → bash \$SCRIPTS_DIR/run_day2.sh"
echo ""
