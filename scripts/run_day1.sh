#!/usr/bin/env bash
# =============================================================================
# run_day1.sh  –  Day 1 launcher (HOST-side)
# lncRNA & mRNA Analysis: QC → Trim → HISAT2 → StringTie → gffcompare → CPC2
#                         → featureCounts → edgeR → cis-targets/GO → Ballgown
#
# Day 1 runs in a SELF-CONTAINED image that already bundles every tool, the 6
# Chr4 RNA-seq libraries (E-MTAB-13532, control vs Tr-PET), the reference and
# the HISAT2 index. This script pulls that image and runs the whole pipeline
# non-interactively, writing all results to ./results/day1_lncRNA on the host.
#
# Run this on your HOST machine (not inside the course container):
#     bash scripts/run_day1.sh            # full pipeline incl. Ballgown
#
# The step-by-step teaching path is exercises/Day1_lncRNA_mRNA_Pipeline.md.
# The Epi-Code – Florence Training School 2026
# =============================================================================
set -euo pipefail

IMAGE="${DAY1_IMAGE:-leogiuffre/lncrna-mnps-workshop:1.0}"
OUT_DIR="${PWD}/results/day1_lncRNA"
PLATFORM_FLAG=()
# Apple Silicon / non-amd64 hosts: force the built (amd64) platform.
if [[ "$(uname -m)" != "x86_64" ]]; then PLATFORM_FLAG=(--platform linux/amd64); fi

echo ""
echo "╔════════════════════════════════════════════════════════════╗"
echo "║  Day 1 – lncRNA & mRNA Differential Expression Pipeline    ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo "  Image  : ${IMAGE}"
echo "  Results: ${OUT_DIR}"
echo ""

mkdir -p "${OUT_DIR}"
docker pull "${PLATFORM_FLAG[@]}" "${IMAGE}"

# Run the in-image pipeline; bind-mount results to the host so they persist.
# JAVA_TOOL_OPTIONS: Trimmomatic (old Java) SIGSEGVs under the image's OpenJDK 22
# C2 JIT; disabling the C2 tier is the reliable fix (validated on this dataset).
docker run --rm "${PLATFORM_FLAG[@]}" \
  -e JAVA_TOOL_OPTIONS=-XX:TieredStopAtLevel=1 \
  -v "${OUT_DIR}:/home/student/results" \
  "${IMAGE}" \
  bash /home/student/workshop/scripts/run_all.sh --with-ballgown

echo ""
echo "Day 1 complete. Tables & figures are under: ${OUT_DIR}"
echo "For the interactive, type-along version see:"
echo "  exercises/Day1_lncRNA_mRNA_Pipeline.md"
