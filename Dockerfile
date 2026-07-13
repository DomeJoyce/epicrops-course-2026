# =============================================================================
# Dockerfile – "The Epi-Code" Practical Course
# Florence Training School, 14–17 July 2026
#
# Purpose : Reproducible bioinformatics environment for plant epigenomics
#           and ncRNA transcriptomics — Day 2 WGBS methylation (Day 1 lncRNA/mRNA
#           ships as its own image, leogiuffre/lncrna-mnps-workshop).
# Base OS  : Ubuntu 22.04 LTS (minimal footprint for laptop compatibility)
# Authors  : Domenico Giosa & Letterio Giuffrè
#            University of Messina, ChiBioFarAm Dept.
# =============================================================================

FROM ubuntu:22.04

LABEL maintainer="Domenico Giosa <dgiosa@unime.it>"
LABEL version="1.0"
LABEL description="Epigenomics & ncRNA transcriptomics practical environment"

# ---------------------------------------------------------------------------
# 1. System environment
# ---------------------------------------------------------------------------
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/Rome
ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8

# ---------------------------------------------------------------------------
# 2. Core system packages
# ---------------------------------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    wget \
    curl \
    git \
    unzip \
    bzip2 \
    gzip \
    pigz \
    ca-certificates \
    libncurses5-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    zlib1g-dev \
    python3 \
    python3-pip \
    python3-dev \
    perl \
    default-jre \
    tzdata \
    less \
    vim \
    nano \
    tree \
    htop \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# ---------------------------------------------------------------------------
# 3. Conda / Mamba (miniforge) – primary package manager
# ---------------------------------------------------------------------------
RUN wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh \
    -O /tmp/miniforge.sh && \
    bash /tmp/miniforge.sh -b -p /opt/conda && \
    rm /tmp/miniforge.sh

ENV PATH="/opt/conda/bin:$PATH"

RUN conda update -n base -c defaults conda -y && \
    conda install -n base -c conda-forge mamba -y

# ---------------------------------------------------------------------------
# 4. Bioinformatics tools via Mamba (pinned for reproducibility)
# ---------------------------------------------------------------------------
# ponytail: names only, not exact patch pins. Exact triples (e.g. samtools=1.19.2
# vs bismark's expected samtools) make the solver fail on a single conflict; letting
# conda pick one mutually-consistent set is what actually builds reliably on every
# machine. Reproducibility is still anchored by the pinned channels + build date.
# Day 2 (WGBS) tools only: Day 1's lncRNA/mRNA pipeline ships in its own image
# (leogiuffre/lncrna-mnps-workshop). bowtie2 is required by Bismark for bisulfite mapping.
# Version notes (verified by running Day 2 end-to-end):
#  - bismark: the exercise + run_day2.sh use the Bismark 3.x "Rust suite" CLI
#    (--genome / --paired-end / --CX / --basename), so '>=0.24' is left to resolve to 3.x.
#    Classic Bismark 0.24/0.25 uses <genome_folder> -1/-2 and rejects --paired-end.
#  - trim-galore: pinned to the classic 0.6.x line (<1) so it emits the *_val_1.fq
#    filenames the alignment step reads; cutadapt left unpinned (base-python compat).
RUN mamba install -y -c conda-forge -c bioconda \
    'fastqc>=0.12' \
    'multiqc>=1.19' \
    'trim-galore>=0.6.10,<1' \
    'cutadapt>=4.6' \
    'bowtie2>=2.5' \
    'samtools>=1.20' \
    'bedtools>=2.31' \
    'bismark>=0.24' \
    seqtk \
    sra-tools \
    parallel \
    && conda clean --all -y

# ---------------------------------------------------------------------------
# 5. Python packages (via conda, into the same env as the bio tools — the
#    system pip built pandas from source and failed; conda ships wheels/binaries
#    and keeps one consistent Python interpreter).
# ---------------------------------------------------------------------------
RUN mamba install -y -c conda-forge -c bioconda \
    pandas \
    numpy \
    matplotlib \
    seaborn \
    scipy \
    pysam \
    biopython \
    jupyter \
    notebook \
    ipykernel \
    && conda clean --all -y

# ---------------------------------------------------------------------------
# 6. R and Bioconductor packages
# ---------------------------------------------------------------------------
# ponytail: pin R minor version (4.3) so Bioconductor 3.18 packages resolve together;
# let the individual package patch versions float. clusterprofiler is lowercase in conda.
RUN mamba install -y -c conda-forge -c bioconda \
    'r-base>=4.3,<4.4' \
    r-essentials \
    bioconductor-genomeinfodbdata \
    bioconductor-methylkit \
    bioconductor-genomicranges \
    bioconductor-rtracklayer \
    bioconductor-enhancedvolcano \
    bioconductor-clusterprofiler \
    bioconductor-biocparallel \
    r-ggplot2 \
    r-pheatmap \
    r-rcolorbrewer \
    r-dplyr \
    r-tidyr \
    r-readr \
    r-tibble \
    r-ggrepel \
    r-patchwork \
    r-viridis \
    r-data.table \
    r-knitr \
    r-rmarkdown \
    && conda clean --all -y

# Two fixes bundled into one fast layer (keeps the big R solve above cached):
#  1. methylKit 1.28 imports `key<-` from data.table, which data.table >=1.15
#     no longer exports -> pin data.table below 1.15.
#  2. bioconda `genomeinfodbdata` downloads its data in a post-link hook that
#     fails silently on flaky networks, leaving GenomeInfoDb (and thus methylKit /
#     GenomicRanges) unloadable -> install it deterministically via Bioconductor.
# The final load test makes the build FAIL here rather than ship a broken image.
RUN mamba install -y -c conda-forge -c bioconda 'r-data.table<1.15' && conda clean --all -y && \
    Rscript -e 'if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org"); BiocManager::install("GenomeInfoDbData", update=FALSE, ask=FALSE, force=TRUE)' && \
    Rscript -e 'suppressMessages({library(methylKit); library(GenomicRanges); library(rtracklayer)}); cat("R stack OK\n")'

# ---------------------------------------------------------------------------
# 7. Working directory structure
# ---------------------------------------------------------------------------
# Day 1 (lncRNA/mRNA) runs in its own image and writes to ./results/day1_lncRNA on the
# host (mounted in via docker-compose). This course image only pre-creates Day 2's tree.
RUN mkdir -p \
    /course/data/raw \
    /course/data/reference/genome \
    /course/data/reference/annotation \
    /course/data/reference/indices/bowtie2 \
    /course/data/reference/indices/bismark \
    /course/results/day2_methylation/01_qc \
    /course/results/day2_methylation/02_trimmed \
    /course/results/day2_methylation/03_bismark_aligned \
    /course/results/day2_methylation/04_methylation \
    /course/results/day2_methylation/05_diffmeth \
    /course/scripts \
    /course/exercises \
    /course/notebooks

WORKDIR /course

# ---------------------------------------------------------------------------
# 8. Copy course scripts and exercises
# ---------------------------------------------------------------------------
COPY scripts/ /course/scripts/
COPY exercises/ /course/exercises/
COPY data/download_data.sh /course/data/

RUN chmod +x /course/scripts/*.sh /course/data/download_data.sh

# ---------------------------------------------------------------------------
# 9. Environment variables for course paths
# ---------------------------------------------------------------------------
ENV COURSE_DIR=/course
ENV DATA_DIR=/course/data
ENV REF_DIR=/course/data/reference
ENV RESULTS_DIR=/course/results
ENV SCRIPTS_DIR=/course/scripts

# ---------------------------------------------------------------------------
# 10. Expose Jupyter port (optional interactive use)
# ---------------------------------------------------------------------------
EXPOSE 8888

# ---------------------------------------------------------------------------
# 11. Default entrypoint: interactive bash with welcome message
# ---------------------------------------------------------------------------
COPY scripts/welcome.sh /usr/local/bin/welcome.sh
RUN chmod +x /usr/local/bin/welcome.sh

RUN echo 'source /usr/local/bin/welcome.sh' >> /root/.bashrc

CMD ["/bin/bash"]
