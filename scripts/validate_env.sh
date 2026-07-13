#!/usr/bin/env bash
# =============================================================================
# validate_env.sh — quick sanity check that the container has every tool the
# course needs. Run this once after starting the container:
#     bash $SCRIPTS_DIR/validate_env.sh
# Exits non-zero if anything is missing, so it doubles as a CI smoke test.
# =============================================================================
set -uo pipefail

fail=0
check_cmd() {
    if command -v "$1" >/dev/null 2>&1; then
        printf "  [OK]   %-24s %s\n" "$1" "$("${@:2}" 2>&1 | head -n1)"
    else
        printf "  [MISS] %-24s NOT FOUND\n" "$1"; fail=1
    fi
}

echo "============================================================"
echo " Environment validation — The Epi-Code practical"
echo "============================================================"
echo "Command-line tools:"
check_cmd fastqc                 fastqc --version
check_cmd multiqc                multiqc --version
check_cmd trim_galore            trim_galore --version
check_cmd cutadapt               cutadapt --version
check_cmd bowtie2                bowtie2 --version
check_cmd samtools               samtools --version
check_cmd bedtools               bedtools --version
check_cmd bismark                bismark --version
check_cmd deduplicate_bismark    deduplicate_bismark --version
check_cmd bismark_methylation_extractor bismark_methylation_extractor --version
check_cmd seqtk                  seqtk
check_cmd R                      R --version

echo ""
echo "R / Bioconductor packages:"
Rscript -e '
pkgs <- c("methylKit","GenomicRanges","rtracklayer","EnhancedVolcano",
          "ggplot2","pheatmap","RColorBrewer","dplyr","tidyr","readr","tibble",
          "patchwork","viridis")
miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
for (p in setdiff(pkgs, miss)) cat(sprintf("  [OK]   %s\n", p))
if (length(miss)) { for (p in miss) cat(sprintf("  [MISS] %s\n", p)); quit(status = 1) }
' || fail=1

echo ""
if [[ "$fail" -eq 0 ]]; then
    echo "ALL GOOD — environment is ready. "
else
    echo "SOME COMPONENTS MISSING — see [MISS] lines above."
fi
exit "$fail"
