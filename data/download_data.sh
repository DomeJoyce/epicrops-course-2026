#!/usr/bin/env bash
# =============================================================================
# download_data.sh
# Download the real, REPLICATED Arabidopsis thaliana WGBS dataset used on
# DAY 2 and subsample it to a laptop-friendly size.
#
#   Day 2 – WGBS (paired-end)               BioProject PRJNA700573
#     control  (Col-0)          : SRR13650161, SRR13650163, SRR13650165
#     salt-stress (Col-0 + NaCl): SRR13650194, SRR13650196, SRR13650198
#
# (DAY 1 – lncRNA/mRNA – needs no download: its 6 Chr4 RNA-seq libraries are
#  already bundled inside the Day 1 image, leogiuffre/lncrna-mnps-workshop.)
#
# DISK-SAFE DOWNLOAD
#   Full WGBS runs are large (~6 GB each). We never store a full run: we stream
#   the gzipped FASTQ straight from the ENA mirror and keep only the first N
#   reads (`zcat | head`), closing the pipe early. Each sample costs a few
#   hundred MB instead of gigabytes. First-N (not random) keeps R1/R2 in sync
#   and is fully deterministic/reproducible for teaching.
# =============================================================================

set -euo pipefail

DATA_DIR="${DATA_DIR:-/course/data}"
RAW_DIR="${DATA_DIR}/raw"

# Read PAIRS to keep per sample (override: download_data.sh <wgbs_pairs>)
WGBS_PAIRS="${1:-2000000}"        # read PAIRS (each mate truncated to this many)

ENA_API="https://www.ebi.ac.uk/ena/portal/api/filereport"

mkdir -p "${RAW_DIR}/wgbs"

echo "============================================================"
echo " Downloading Day 2 WGBS dataset (3 replicates / condition)"
echo " WGBS read-pairs/sample : ${WGBS_PAIRS}"
echo "============================================================"

# ── Resolve ENA FASTQ URL(s) for an accession (';'-separated) ────────────────
ena_fastq_urls() {
    local SRR=$1
    curl -s "${ENA_API}?accession=${SRR}&result=read_run&fields=fastq_ftp&format=tsv" \
        | tail -n +2 | cut -f2
}

# ── Stream a gzipped FASTQ and keep only the first N reads (4N lines) ─────────
# pipefail is OFF here: `head` closing the pipe makes curl/zcat exit non-zero
# by design (SIGPIPE); we validate the result by counting reads afterwards.
# Retries a few times, because the ENA mirror occasionally drops the first
# connection (a transient empty download that must NOT abort the whole run).
stream_head() {
    local URL=$1 OUT=$2 NREADS=$3
    local LINES=$(( NREADS * 4 ))
    local attempt got
    for attempt in 1 2 3 4; do
        ( set +e +o pipefail
          curl -s --retry 3 --retry-delay 3 --connect-timeout 30 "https://${URL}" \
            | zcat 2>/dev/null | head -n "${LINES}" | gzip > "${OUT}"
        )
        got=$(zcat "${OUT}" 2>/dev/null | head -n "${LINES}" | wc -l)
        if [[ "${got}" -ge 4 && $(( got % 4 )) -eq 0 ]]; then
            echo "    -> ${OUT}  ($(( got / 4 )) reads)"
            return 0
        fi
        echo "    .. attempt ${attempt} for $(basename "${OUT}") got ${got} lines; retrying..." >&2
        sleep 5
    done
    echo "    !! ERROR: ${OUT} still invalid after retries. Check your internet / the ENA mirror." >&2
    return 1
}

download_pe() {
    local SRR=$1 OUT_DIR=$2
    local O1="${OUT_DIR}/${SRR}_1.fastq.gz" O2="${OUT_DIR}/${SRR}_2.fastq.gz"
    if [[ -s "${O1}" && -s "${O2}" ]]; then echo "[=] ${SRR} present, skipping."; return 0; fi
    echo "[+] ${SRR} (paired-end)..."
    local urls u1 u2
    urls=$(ena_fastq_urls "${SRR}")
    u1=$(echo "${urls}" | cut -d';' -f1)
    u2=$(echo "${urls}" | cut -d';' -f2)
    [[ -n "${u1}" && -n "${u2}" && "${u1}" != "${u2}" ]] \
        || { echo "    !! ${SRR} is not paired on ENA (${urls})" >&2; return 1; }
    stream_head "${u1}" "${O1}" "${WGBS_PAIRS}"
    stream_head "${u2}" "${O2}" "${WGBS_PAIRS}"
}

# ── Day 2: WGBS ──────────────────────────────────────────────────────────────
echo ""; echo "--- Day 2: WGBS (control vs salt stress) ---"
for SRR in SRR13650161 SRR13650163 SRR13650165 SRR13650194 SRR13650196 SRR13650198; do
    download_pe "${SRR}" "${RAW_DIR}/wgbs"
done

# ── Sample sheet ─────────────────────────────────────────────────────────────
echo ""; echo "--- Writing sample sheet ---"

cat > "${RAW_DIR}/wgbs/samplesheet.tsv" << 'EOF'
sample_id	condition	replicate	fastq_r1	fastq_r2
SRR13650161	control	1	SRR13650161_1.fastq.gz	SRR13650161_2.fastq.gz
SRR13650163	control	2	SRR13650163_1.fastq.gz	SRR13650163_2.fastq.gz
SRR13650165	control	3	SRR13650165_1.fastq.gz	SRR13650165_2.fastq.gz
SRR13650194	salt	1	SRR13650194_1.fastq.gz	SRR13650194_2.fastq.gz
SRR13650196	salt	2	SRR13650196_1.fastq.gz	SRR13650196_2.fastq.gz
SRR13650198	salt	3	SRR13650198_1.fastq.gz	SRR13650198_2.fastq.gz
EOF

echo ""
echo "============================================================"
echo " Download COMPLETE"
echo "   WGBS : ${RAW_DIR}/wgbs/   (6 samples, control vs salt)"
echo "============================================================"
ls -lh "${RAW_DIR}/wgbs/"
