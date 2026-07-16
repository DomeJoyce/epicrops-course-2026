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

# control (Col-0) then salt-stress (Col-0 + NaCl) — used by both the download
# loop and the final completeness check below, so it lives in one place.
SAMPLES=(SRR13650161 SRR13650163 SRR13650165 SRR13650194 SRR13650196 SRR13650198)

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

# ── Download the first N reads of a remote gzipped FASTQ (robust/resumable) ───
# We only need the first N reads, not the whole multi-GB run. The old approach
# streamed the entire file through `zcat | head` and let `head` close the pipe.
# That is fragile: on a network that STALLS right after the TLS handshake — very
# common inside Docker (broken IPv6, an MTU mismatch on the bridge, or a
# proxy/firewall) — curl delivers only its first buffer (~1 KB, i.e. ~50
# decompressed lines) and every plain retry restarts from byte 0 and stalls
# again at the same point (the "got 50 lines" failure).
#
# Instead we fetch a BOUNDED, RESUMABLE byte-range prefix that already contains
# enough reads, then cut it to exactly N. Because each retry resumes from the
# bytes already on disk (Range: <have>-), a transient stall makes progress
# instead of looping forever on the same 1 KB. -4 avoids broken IPv6, -f fails
# on HTTP errors (so a proxy error page never masquerades as data), and the
# lenient --speed-limit/--speed-time only aborts a truly DEAD transfer.
stream_head() {
    local URL=$1 OUT=$2 NREADS=$3
    local LINES=$(( NREADS * 4 ))
    local tmp="${OUT}.part"
    # ~62 compressed bytes/read for this dataset; 100 gives ~1.6x headroom so a
    # worse-compressing run still yields >= N reads. Floor at 1 MB for tiny N.
    local want=$(( NREADS * 100 )); (( want < 1048576 )) && want=1048576
    rm -f "${tmp}"; : > "${tmp}"
    local have=0 stall=0 iter=0 new
    while (( have < want && stall < 6 && iter < 40 )); do
        iter=$(( iter + 1 ))
        # Resume: request the range starting at the bytes we already have and
        # APPEND. A byte-exact prefix of a gzip file is itself valid to decode.
        curl -s -4 -L -f --connect-timeout 30 --max-time 600 \
             --speed-limit 1024 --speed-time 30 --retry 3 --retry-delay 3 \
             -r "${have}-$(( want - 1 ))" "https://${URL}" >> "${tmp}" 2>/dev/null || true
        new=$(stat -c%s "${tmp}" 2>/dev/null || echo 0)
        if (( new > have )); then
            have=${new}; stall=0
        else
            stall=$(( stall + 1 ))
            echo "    .. $(basename "${OUT}"): no progress at $(( have/1024 )) KB (retry ${stall}/6)" >&2
            sleep 3
        fi
    done
    # Cut to exactly N whole reads, recompress, then drop the temp prefix.
    zcat "${tmp}" 2>/dev/null | head -n "${LINES}" | gzip > "${OUT}"
    rm -f "${tmp}"
    local got
    got=$(zcat "${OUT}" 2>/dev/null | wc -l)
    got=$(( got - got % 4 ))                       # ignore a partial trailing read
    # Success only if we actually got (essentially) all N reads. The prefix is
    # sized with ~1.6x headroom, so a healthy transfer yields exactly N; anything
    # much smaller means the transfer STALLED and was truncated. We must NOT
    # accept that silently — a handful of reads would poison the whole pipeline —
    # and we delete the partial file so a re-run re-fetches it (otherwise the
    # gzip-valid stub would be treated as "already downloaded").
    local need=$(( LINES * 95 / 100 ))
    if [[ "${got}" -ge "${need}" ]]; then
        echo "    -> ${OUT}  ($(( got / 4 )) reads)"
        return 0
    fi
    rm -f "${OUT}"
    echo "    !! ERROR: ${OUT} truncated — got only $(( got / 4 )) of ${NREADS} reads." >&2
    echo "       The download connected but the bulk transfer STALLED after the first" >&2
    echo "       kilobyte. This is the network, not the data (the ENA file is a valid" >&2
    echo "       ~3 GB gzip). The usual cause inside Docker is an MTU mismatch on the" >&2
    echo "       bridge network (full-size packets are silently dropped). Fixes:" >&2
    echo "         * set the Docker bridge MTU to 1400 (see README / daemon.json), or" >&2
    echo "         * run this download on the host, outside the container, or" >&2
    echo "         * switch to a network without a proxy/firewall in the way." >&2
    return 1
}

# A pair counts as "already done" only if BOTH mates exist AND are valid gzip.
# A download killed mid-stream leaves a truncated, non-empty file that `-s`
# would wrongly accept (feeding corrupt reads to the pipeline); `gzip -t`
# rejects it so the next run re-fetches it.
pair_ok() { [[ -s "$1" && -s "$2" ]] && gzip -t "$1" 2>/dev/null && gzip -t "$2" 2>/dev/null; }

download_pe() {
    local SRR=$1 OUT_DIR=$2
    local O1="${OUT_DIR}/${SRR}_1.fastq.gz" O2="${OUT_DIR}/${SRR}_2.fastq.gz"
    if pair_ok "${O1}" "${O2}"; then echo "[=] ${SRR} present, skipping."; return 0; fi
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
for SRR in "${SAMPLES[@]}"; do
    # One flaky sample must NOT abort the rest (set -e would). Log and carry on;
    # every file is resumable, so re-running the script fills in what's missing.
    download_pe "${SRR}" "${RAW_DIR}/wgbs" || echo "[!] ${SRR} incomplete — re-run to resume" >&2
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
missing=0
for SRR in "${SAMPLES[@]}"; do
    pair_ok "${RAW_DIR}/wgbs/${SRR}_1.fastq.gz" "${RAW_DIR}/wgbs/${SRR}_2.fastq.gz" \
        || { echo " [MISSING] ${SRR}"; missing=$(( missing + 1 )); }
done
if [[ ${missing} -eq 0 ]]; then
    echo " Download COMPLETE — all ${#SAMPLES[@]} WGBS samples present"
else
    echo " ${missing}/${#SAMPLES[@]} sample(s) still incomplete — re-run this script to resume"
fi
echo "   WGBS : ${RAW_DIR}/wgbs/   (control vs salt)"
echo "============================================================"
ls -lh "${RAW_DIR}/wgbs/"
[[ ${missing} -eq 0 ]]   # exit non-zero if anything is still missing (preflight-friendly)
