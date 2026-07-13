#!/usr/bin/env Rscript
# =============================================================================
# day2_methylkit.R
# Differential methylation analysis of WGBS data using methylKit.
# Covers CpG, CHG and CHH contexts; produces DMR tables and publication
# quality visualisations.
#
# The Epi-Code – Florence Training School 2026
# Day 2: Whole-Genome Bisulfite Sequencing Analysis
# =============================================================================

suppressPackageStartupMessages({
  library(methylKit)
  library(GenomicRanges)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(RColorBrewer)
  library(patchwork)
  library(viridis)
})

# ─────────────────────────────────────────────────────────────────────────────
# 0. Command-line arguments
# ─────────────────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)

meth_dir    <- ifelse(length(args) >= 1, args[1],
  "/course/results/day2_methylation/04_methylation")
samplesheet <- ifelse(length(args) >= 2, args[2],
  "/course/data/raw/wgbs/samplesheet.tsv")
out_dir     <- ifelse(length(args) >= 3, args[3],
  "/course/results/day2_methylation/05_diffmeth")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("============================================================\n")
cat("  methylKit – Differential Methylation Analysis\n")
cat("  Input  :", meth_dir, "\n")
cat("  Samples:", samplesheet, "\n")
cat("  Output :", out_dir, "\n")
cat("============================================================\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# 1. Load bismark cytosine reports
# ─────────────────────────────────────────────────────────────────────────────
cat("[1/7] Loading Bismark cytosine reports...\n")

meta <- read_tsv(samplesheet, show_col_types = FALSE)

# Build file list (Bismark cytosine report). Bismark writes ".CX_report.txt";
# accept a gzipped variant too so the script works regardless of --gzip.
pick_report <- function(id) {
  cand <- file.path(meth_dir, paste0(id, c(".CX_report.txt", ".CX_report.txt.gz")))
  hit  <- cand[file.exists(cand)]
  if (length(hit) == 0) stop("No CX report found for sample ", id,
                             " (looked for: ", paste(cand, collapse = ", "), ")")
  hit[1]
}
cyt_files <- vapply(meta$sample_id, pick_report, character(1))

cat("  Files to load:\n")
for (f in cyt_files) cat("   ", f, "\n")

# Treatment vector: force control = 0, everything else = 1 (reference-safe,
# independent of row order in the sample sheet).
treatment_v <- ifelse(meta$condition == "control", 0L, 1L)
conditions  <- c("control", setdiff(unique(meta$condition), "control"))

cat("  Treatment vector:", paste(treatment_v, collapse = " "), "\n")
cat("  (0 =", conditions[1], "| 1 =", conditions[2], ")\n\n")

# ── 1a. CpG context ──────────────────────────────────────────────────────────
cat("[1a] Loading CpG context...\n")

myobj_CpG <- methRead(
  as.list(cyt_files),
  sample.id  = as.list(meta$sample_id),
  assembly   = "TAIR10_Chr4",
  treatment  = treatment_v,
  context    = "CpG",
  mincov     = 1,          # keep every covered cytosine: the regional (tiled)
                           # DMR step below aggregates them per window. Shallow,
                           # subsampled WGBS has too little per-cytosine depth to
                           # survive a high mincov in all samples at once.
  pipeline   = "bismarkCytosineReport"
)

# ── 1b. CHG context ──────────────────────────────────────────────────────────
cat("[1b] Loading CHG context...\n")

myobj_CHG <- methRead(
  as.list(cyt_files),
  sample.id  = as.list(meta$sample_id),
  assembly   = "TAIR10_Chr4",
  treatment  = treatment_v,
  context    = "CHG",
  mincov     = 1,          # see CpG note above
  pipeline   = "bismarkCytosineReport"
)

# ── 1c. CHH context ──────────────────────────────────────────────────────────
cat("[1c] Loading CHH context...\n")

myobj_CHH <- methRead(
  as.list(cyt_files),
  sample.id  = as.list(meta$sample_id),
  assembly   = "TAIR10_Chr4",
  treatment  = treatment_v,
  context    = "CHH",
  mincov     = 1,          # see CpG note above
  pipeline   = "bismarkCytosineReport"
)

# ─────────────────────────────────────────────────────────────────────────────
# 2. Descriptive statistics and QC plots
# ─────────────────────────────────────────────────────────────────────────────
cat("\n[2/7] QC plots and descriptive statistics...\n")

## 2a. Methylation distribution histograms (CpG)
pdf(file.path(out_dir, "QC_methylation_histograms.pdf"), width = 10, height = 6)
getMethylationStats(myobj_CpG[[1]], plot = TRUE, both.strands = FALSE)
title(paste("CpG Methylation –", meta$sample_id[1]))
getMethylationStats(myobj_CpG[[2]], plot = TRUE, both.strands = FALSE)
title(paste("CpG Methylation –", meta$sample_id[2]))
dev.off()

## 2b. Coverage statistics
pdf(file.path(out_dir, "QC_coverage_histograms.pdf"), width = 10, height = 6)
getCoverageStats(myobj_CpG[[1]], plot = TRUE, both.strands = FALSE)
title(paste("Coverage –", meta$sample_id[1]))
getCoverageStats(myobj_CpG[[2]], plot = TRUE, both.strands = FALSE)
title(paste("Coverage –", meta$sample_id[2]))
dev.off()

cat("  → QC histograms saved\n")

# ─────────────────────────────────────────────────────────────────────────────
# 3. Filter by coverage, unite samples
# ─────────────────────────────────────────────────────────────────────────────
cat("\n[3/7] Filtering by coverage and uniting samples...\n")

run_analysis_for_context <- function(myobj, context_name) {

  # Filter: min coverage 3x, max 99.9th percentile (remove PCR duplicates).
  # 3x (rather than the textbook 10x) is deliberate: this is a subsampled, low-
  # depth WGBS teaching set, and requiring 10x in *every* sample leaves no shared
  # cytosine to unite. 3x keeps the "covered in all samples" logic while still
  # giving a usable number of single-cytosine sites.
  myobj_filtered <- filterByCoverage(myobj,
    lo.count    = 3,
    lo.perc     = NULL,
    hi.count    = NULL,
    hi.perc     = 99.9
  )

  # Unite: keep only positions covered in ALL samples
  meth_united <- unite(myobj_filtered, destrand = ifelse(context_name == "CpG", TRUE, FALSE))

  cat("  [", context_name, "] Sites after filtering:", nrow(meth_united), "\n")
  return(meth_united)
}

meth_CpG <- run_analysis_for_context(myobj_CpG, "CpG")
meth_CHG <- run_analysis_for_context(myobj_CHG, "CHG")
meth_CHH <- run_analysis_for_context(myobj_CHH, "CHH")

# ─────────────────────────────────────────────────────────────────────────────
# 4. Sample correlation and clustering
# ─────────────────────────────────────────────────────────────────────────────
cat("\n[4/7] Sample correlation and hierarchical clustering...\n")

pdf(file.path(out_dir, "Correlation_and_Clustering_CpG.pdf"), width = 8, height = 6)
getCorrelation(meth_CpG, plot = TRUE)
title("Sample Pearson correlation – CpG methylation")

clusterSamples(meth_CpG, dist = "correlation", method = "ward",
               plot = TRUE)
title("Hierarchical clustering – CpG methylation")

PCASamples(meth_CpG, screeplot = TRUE)
title("PCA screeplot – CpG methylation")

PCASamples(meth_CpG)
title("PCA – CpG methylation")
dev.off()

cat("  → Correlation and clustering plots saved\n")

# ─────────────────────────────────────────────────────────────────────────────
# 5. Differential methylation at cytosine level
# ─────────────────────────────────────────────────────────────────────────────
cat("\n[5/7] Testing differential methylation at cytosine level...\n")

diff_meth_cytosine <- function(meth_obj, context_name, q = 0.01, diff = 25) {

  myDiff <- calculateDiffMeth(meth_obj,
    mc.cores   = max(1, parallel::detectCores() - 1),
    test       = "F",
    overdispersion = "MN"
  )

  # Significant DMCs (|% diff| ≥ 25, q < 0.01)
  myDiff_sig <- getMethylDiff(myDiff,
    difference = diff,
    qvalue     = q,
    type       = "all"
  )

  cat("  [", context_name, "] Total DMCs:", nrow(myDiff_sig), "\n")
  cat("           Hypermethylated:", nrow(getMethylDiff(myDiff, diff, q, "hyper")), "\n")
  cat("           Hypomethylated :", nrow(getMethylDiff(myDiff, diff, q, "hypo")), "\n")

  # Save to file
  write.csv(as.data.frame(myDiff_sig),
    file.path(out_dir, paste0("DMC_", context_name, ".csv")),
    row.names = FALSE
  )

  return(myDiff)
}

myDiff_CpG <- diff_meth_cytosine(meth_CpG, "CpG")
myDiff_CHG <- diff_meth_cytosine(meth_CHG, "CHG")
myDiff_CHH <- diff_meth_cytosine(meth_CHH, "CHH")

# ─────────────────────────────────────────────────────────────────────────────
# 6. Tiling window analysis (DMRs, 200 bp tiles)
# ─────────────────────────────────────────────────────────────────────────────
cat("\n[6/7] Tiling window analysis (200 bp bins – DMRs)...\n")

tile_analysis <- function(myobj_raw, context_name, q = 0.05, diff = 20) {

  # Tile each sample's RAW per-cytosine counts into 200 bp non-overlapping
  # windows FIRST, then unite. This sums the reads within a window *before*
  # requiring the window to be shared across samples — the standard methylKit
  # approach for shallow / subsampled WGBS. (Uniting single cytosines first and
  # tiling afterwards, as one might expect, leaves almost no shared positions at
  # this coverage, so no DMR would ever be called.)
  tiles      <- tileMethylCounts(myobj_raw, win.size = 200, step.size = 200,
                                  cov.bases = 1)
  meth_tiles <- unite(tiles, destrand = FALSE)

  myDiff_tiles <- calculateDiffMeth(meth_tiles,
    mc.cores = max(1, parallel::detectCores() - 1),
    test     = "F",
    overdispersion = "MN"
  )

  myDiff_sig <- getMethylDiff(myDiff_tiles,
    difference = diff,
    qvalue     = q,
    type       = "all"
  )

  cat("  [", context_name, "] DMRs (200 bp tiles):", nrow(myDiff_sig), "\n")

  write.csv(as.data.frame(myDiff_sig),
    file.path(out_dir, paste0("DMR_200bp_", context_name, ".csv")),
    row.names = FALSE
  )

  return(myDiff_tiles)
}

# Pass the RAW methylRawList objects (not the united cytosine tables) so tiling
# aggregates coverage per window before the samples are intersected.
DMR_CpG <- tile_analysis(myobj_CpG, "CpG")
DMR_CHG <- tile_analysis(myobj_CHG, "CHG")
DMR_CHH <- tile_analysis(myobj_CHH, "CHH")

# ─────────────────────────────────────────────────────────────────────────────
# 7. Visualisations
# ─────────────────────────────────────────────────────────────────────────────
cat("\n[7/7] Generating visualisations...\n")

## 7a. Methylation level per context (bar plot) — mean over replicates by group
# Average the per-cytosine methylation fraction across the replicate columns of
# each condition (treatment 0 = control, 1 = treatment), for every context.
group_mean_meth <- function(mobj, tvec) {
  d <- getData(mobj)
  fr <- sapply(seq_along(tvec), function(i)
    d[[paste0("numCs", i)]] / d[[paste0("coverage", i)]])
  c(mean(fr[, tvec == 0, drop = FALSE], na.rm = TRUE),
    mean(fr[, tvec == 1, drop = FALSE], na.rm = TRUE)) * 100
}

grp_labels <- conditions[1:2]   # e.g. c("control", "salt")
ctx_long <- do.call(rbind, lapply(
  list(CpG = meth_CpG, CHG = meth_CHG, CHH = meth_CHH),
  function(m) group_mean_meth(m, treatment_v)
)) |>
  as.data.frame() |>
  setNames(grp_labels)
ctx_long$Context <- rownames(ctx_long)
ctx_long <- tidyr::pivot_longer(ctx_long, -Context,
  names_to = "Condition", values_to = "MeanMethylation")

p_context <- ggplot(ctx_long, aes(x = Context, y = MeanMethylation,
                                   fill = Condition)) +
  geom_col(position = "dodge", color = "white", linewidth = 0.3) +
  scale_fill_brewer(palette = "Set1") +
  labs(
    title   = "Mean DNA methylation level by cytosine context",
    y       = "Mean methylation (%)",
    x       = "Cytosine context",
    caption = "Arabidopsis thaliana Chr4 | Epi-Code 2026"
  ) +
  theme_classic(base_size = 13) +
  ylim(0, 100)

ggsave(file.path(out_dir, "Methylation_context_barplot.pdf"),
       p_context, width = 6, height = 5)
ggsave(file.path(out_dir, "Methylation_context_barplot.png"),
       p_context, width = 6, height = 5, dpi = 150)
cat("  → Context barplot saved\n")

## 7b. Chromosome methylation track (CpG, sliding window average)
pdf(file.path(out_dir, "Chr4_methylation_track_CpG.pdf"), width = 12, height = 4)
plot(myDiff_CpG, chromosome = "Chr4", col = c("firebrick", "steelblue"),
     lwd = 2, main = "CpG methylation track – Chr4")
dev.off()
cat("  → Chr4 methylation track saved\n")

## 7c. Diffmeth scatter plot: % difference vs coverage
diff_data_CpG <- as.data.frame(myDiff_CpG)
# A methylDiff object has no coverage columns — those live in the united object
# (meth_CpG), one per sample (coverage1..coverageN). Sum them for total per-position
# coverage; calculateDiffMeth keeps the same rows/order, so they line up.
.cov <- getData(meth_CpG)
.cc  <- grep("^coverage", colnames(.cov))
diff_data_CpG$total_cov <- if (nrow(.cov) == nrow(diff_data_CpG))
                             rowSums(.cov[, .cc, drop = FALSE]) else NA_integer_
p_scatter <- ggplot(diff_data_CpG,
    aes(x = total_cov, y = meth.diff,
        color = ifelse(qvalue < 0.01 & abs(meth.diff) >= 25,
                       "Significant", "Not significant"))) +
  geom_point(size = 0.6, alpha = 0.5) +
  scale_color_manual(values = c("Significant" = "firebrick",
                                 "Not significant" = "grey70")) +
  scale_x_log10() +
  labs(
    title  = "CpG differential methylation – coverage vs % difference",
    x      = "Total coverage (log10)",
    y      = "Methylation difference (%)",
    color  = NULL,
    caption = "Arabidopsis thaliana Chr4 | Epi-Code 2026"
  ) +
  theme_classic(base_size = 12) +
  geom_hline(yintercept = c(-25, 25), linetype = "dashed", color = "grey40")

ggsave(file.path(out_dir, "DiffMeth_scatter_CpG.pdf"),
       p_scatter, width = 8, height = 5)
ggsave(file.path(out_dir, "DiffMeth_scatter_CpG.png"),
       p_scatter, width = 8, height = 5, dpi = 150)
cat("  → Diffmeth scatter plot saved\n")

## Session info
sink(file.path(out_dir, "session_info.txt"))
sessionInfo()
sink()

cat("\n============================================================\n")
cat("  methylKit analysis COMPLETE\n")
cat("  Results saved to:", out_dir, "\n")
cat("============================================================\n\n")
