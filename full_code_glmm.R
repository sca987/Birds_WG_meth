#!/usr/bin/env Rscript
# ================================================================
# Title: Fledge Methylation Filtering Pipeline
# Author: Susan C. Anderson
# Date: 2025-11-04
# Repository: https://github.com/<your-username>/<your-repo>
# Description:
#   Filters CpG methylation data from chestnut-crowned babbler fledglings.
#   Sequential filters applied:
#     1. Remove low-quality individuals
#     2. Remove low-coverage sites (<10×)
#     3. Remove high-coverage outliers (>99.9th percentile per individual)
#     4. Keep sites present in ≥60% of individuals per carers group
#        (thresholds dynamically calculated)
#     5. Remove low-variance sites (bottom 10%)
# Output:
#   all_fledge_filtered_for_model_2025_v1.csv.gz
# ================================================================

# suppressPackageStartupMessages(library(data.table))

# -------------------- 1. Load Data --------------------
fledge_all <- fread("all_fledge_merged_long_20251019_v3_with_missing.csv.gz")
cat("✅ Loaded dataset with", nrow(fledge_all), "rows and",
    uniqueN(fledge_all$sample), "individuals\n")

# -------------------- 2. Remove Bad Individuals --------------------
bad_ids <- c(79596, 79686, 79690, 62603, 62644, 79566, 79692, 79688)
fledge_all <- fledge_all[!sample %in% bad_ids]
cat("Removed bad individuals. Remaining:", uniqueN(fledge_all$sample), "samples\n")

# -------------------- 3. Coverage ≥10 --------------------
fledge_10x <- fledge_all[!is.na(total_coverage) & total_coverage >= 10]
cat("After ≥10x filter:", nrow(fledge_10x), "rows\n")

# -------------------- 4. 99.9th Percentile per Individual --------------------
fledge_99 <- fledge_10x[, {
  thresh <- quantile(total_coverage, probs = 0.999, na.rm = TRUE)
  .SD[total_coverage <= thresh]
}, by = sample]
cat("After per-sample 99.9% filter:", nrow(fledge_99), "rows\n")

# -------------------- 5. ≥60% Presence per Carer Group (Dynamic thresholds) --------------------
# Count per-site sample presence by carers group
site_group_counts <- fledge_99[, .(n_individuals = uniqueN(sample)),
                               by = .(chr, pos, carers)]

# Total number of individuals per carers group
group_sizes <- fledge_99[, .(total_individuals = uniqueN(sample)), by = carers]
print(group_sizes)

# Dynamically calculate minimum count (60% of group size, rounded up)
min_thresholds <- group_sizes[, .(
  carers,
  min_individuals = ceiling(0.60 * total_individuals)
)]
print(min_thresholds)

# Merge and calculate proportion present
site_group_counts <- merge(site_group_counts, group_sizes, by = "carers")
site_group_counts[, prop_present := n_individuals / total_individuals]

# Merge in dynamic thresholds
site_group_counts <- merge(site_group_counts, min_thresholds, by = "carers")

# Keep sites meeting both 60% and minimum count thresholds
site_group_pass <- site_group_counts[
  prop_present >= 0.60 & n_individuals >= min_individuals
]

# Identify sites passing for both carers groups
site_group_pass_counts <- site_group_pass[, .N, by = .(chr, pos)]
sites_to_keep <- site_group_pass_counts[N == uniqueN(fledge_99$carers)]

# Subset dataset to retained sites
fledge_60 <- fledge_99[sites_to_keep, on = c("chr", "pos")]
cat("After 60% + min-individual filter:",
    uniqueN(fledge_60[, .(chr, pos)]), "unique sites retained\n")

# -------------------- 6. Remove Low-Variance Sites (bottom 10%) --------------------
site_var <- fledge_60[, .(var_meth = var(perc_meth, na.rm = TRUE)), by = .(chr, pos)]
threshold_var <- quantile(site_var$var_meth, probs = 0.10, na.rm = TRUE)
sites_high_var <- site_var[var_meth > threshold_var]
fledge_final <- fledge_60[sites_high_var, on = c("chr", "pos")]
cat("After variance filter:",
    uniqueN(fledge_final[, .(chr, pos)]), "unique sites retained\n")

# -------------------- 7. Save Final Dataset --------------------
out_file <- "all_fledge_filtered_for_model_2025_v1.csv.gz"
fwrite(fledge_final, out_file, compress = "gzip")

cat("\n✅ Final dataset saved:", out_file, "\n")
cat("Final unique sites:", uniqueN(fledge_final[, .(chr, pos)]), "\n")
cat("Final individuals:", uniqueN(fledge_final$sample), "\n")

# ================================================================
# End of Script
# ================================================================
#!/bin/bash
#SBATCH --job-name=glmm_beta_parallel
#SBATCH --output=logs/glmm_beta_parallel_%j.out
#SBATCH --error=logs/glmm_beta_parallel_%j.err
#SBATCH --partition=himem
#SBATCH --cpus-per-task=48
#SBATCH --mem=512G
#SBATCH --time=48:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=susan.c.anderson@coyotes.usd.edu

# ================================================================
# Script: run_glmm_beta_parallel.slurm
# Purpose:
#   Run the site-specific beta regression models (GLMM) in parallel
#   using glmmTMB across CpG sites from the fledge dataset.
# ================================================================

# --- 1. Load environment ---
module load R/4.4.2

# --- 2. Define environment variables used in glmm_beta_parallel.R ---
export FLEDGE_INPUT="all_fledge_filtered_for_model_2025_v1.csv.gz"  # input file
export OUT_DIR="glmm_beta_out"                                      # output directory
export OUT_ALL="glmm_results_beta_all.csv.gz"                       # full model output
export OUT_SIG="significant_glmm_sites_beta_FDR05.csv.gz"           # significant sites
export N_CORES=48                                                   # number of parallel cores
export FDR_ALPHA=0.05                                               # FDR threshold
export SITE_LIMIT=0                                                 # 0 = use all sites

# --- 3. Create logs directory if missing ---
mkdir -p logs

# --- 4. Run the analysis ---
echo "Job started on $(hostname) at $(date)"
echo "Running GLMM Beta model using $N_CORES cores and $SLURM_MEM_PER_NODE memory"

Rscript glmm_beta_parallel.R

echo "Job finished at $(date)"
# ================================================================
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(glmmTMB)
  library(parallel)
})

## ---------------- Config (env overrides) ----------------
infile     <- Sys.getenv("FLEDGE_INPUT", "all_fledge_filtered_for_model_2025_v1.csv.gz")
out_dir    <- Sys.getenv("OUT_DIR", "glmm_beta")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_all    <- file.path(out_dir, Sys.getenv("OUT_ALL", "glmm_results_beta_all.csv.gz"))
out_sig    <- file.path(out_dir, Sys.getenv("OUT_SIG", "significant_glmm_sites_beta_FDR05.csv.gz"))
parts_dir  <- file.path(out_dir, "parts")
dir.create(parts_dir, showWarnings = FALSE, recursive = TRUE)

n_cores    <- as.integer(Sys.getenv("N_CORES", "24"))
alpha_fdr  <- as.numeric(Sys.getenv("FDR_ALPHA", "0.05"))
site_limit <- as.integer(Sys.getenv("SITE_LIMIT", "0"))   # 0 = use all

## ---------------- Load data ----------------
dt <- fread(infile)

# make sure types are correct
dt[, carers       := as.factor(carers)]
if ("2" %in% levels(dt$carers)) dt[, carers := relevel(carers, ref = "2")]
dt[, natal_group  := as.factor(natal_group)]
dt[, age_scaled   := scale(age)]

# beta regression cannot have exactly 0 or 100
dt[, perc_meth_beta := fifelse(perc_meth <= 0,   0.0001,
                               fifelse(perc_meth >=100, 0.9999, perc_meth/100))]

# site id
dt[, site_id := paste(chr, pos, sep = "_")]

## ---------------- Sites & sharding ----------------
sites <- unique(dt$site_id)
if (site_limit > 0) sites <- sites[seq_len(min(site_limit, length(sites)))]
n_sites <- length(sites)
message("Total sites to process: ", n_sites)

shards <- split(sites, cut(seq_along(sites), breaks = n_cores, labels = FALSE))

## ---------------- Per-site fitter ----------------
fit_one_site <- function(s) {
  sub <- dt[site_id == s]
  sub <- sub[complete.cases(perc_meth_beta, carers, age_scaled, natal_group)]
  
  # must have both carers groups and some variance
  if (nrow(sub) < 10L || length(unique(sub$carers)) < 2L || var(sub$perc_meth_beta) == 0) return(NULL)
  
  fit <- try(
    glmmTMB(
      perc_meth_beta ~ carers + age_scaled + (1 | natal_group),
      data = sub,
      family = beta_family(link = "logit")
    ),
    silent = TRUE
  )
  if (inherits(fit, "try-error")) return(NULL)
  
  co <- try(summary(fit)$coefficients$cond, silent = TRUE)
  if (inherits(co, "try-error") || is.null(co)) return(NULL)
  
  data.table(
    site_id  = s,
    term     = rownames(co),
    estimate = co[, "Estimate"],
    se       = co[, "Std. Error"],
    z        = co[, "z value"],
    p        = co[, "Pr(>|z|)"],
    n_obs    = nrow(sub)
  )
}

## ---------------- Worker for one shard ----------------
run_shard <- function(idx) {
  svec <- shards[[idx]]
  part_file <- file.path(parts_dir, sprintf("part_%03d.csv.gz", idx))
  if (file.exists(part_file)) {
    message(sprintf("[Shard %d] exists, skipping.", idx))
    return(part_file)
  }
  message(sprintf("[Shard %d] sites=%d", idx, length(svec)))
  out_list <- vector("list", length(svec))
  for (i in seq_along(svec)) {
    if (i %% 1000 == 0) message(sprintf("[Shard %d] %d / %d", idx, i, length(svec)))
    out_list[[i]] <- fit_one_site(svec[[i]])
  }
  res <- rbindlist(out_list, use.names = TRUE, fill = TRUE)
  if (nrow(res)) fwrite(res, part_file, compress = "gzip")
  part_file
}

## ---------------- Run in parallel ----------------
if (.Platform$OS.type == "windows") {
  part_files <- lapply(seq_along(shards), run_shard)
} else {
  part_files <- mclapply(seq_along(shards), run_shard, mc.cores = n_cores)
}

## ---------------- Combine parts, apply FDR ----------------
parts <- list.files(parts_dir, pattern = "^part_\\d{3}\\.csv\\.gz$", full.names = TRUE)
if (!length(parts)) stop("No part files found.")

message("Combining ", length(parts), " shards ...")
all_res <- rbindlist(lapply(parts, fread), use.names = TRUE, fill = TRUE)

all_res <- all_res[!is.na(p) & !is.na(term) & !is.na(site_id)]

# 🆕 Save raw, unadjusted EWAS results for DMRFF
out_raw <- file.path(out_dir, "glmm_results_beta_raw_unadjusted.csv.gz")
fwrite(all_res, out_raw, compress = "gzip")
message("Saved raw model output (for DMRFF): ", out_raw)



# FDR corrections
all_res[, p_adj_by_term := p.adjust(p, method = "BH"), by = term]  # within effect
all_res[, p_adj_global  := p.adjust(p, method = "BH")]             # across all

# Save full
fwrite(all_res, out_all, compress = "gzip")

# Significant by *GLOBAL* FDR ≤ 0.01
sig_global <- all_res[p_adj_global <= 0.01 & term != "(Intercept)"]
fwrite(sig_global, out_sig, compress = "gzip")

# Significant by *TERM-SPECIFIC* FDR ≤ 0.01
out_sig_byterm <- file.path(out_dir, "significant_glmm_sites_beta_byterm_FDR01.csv.gz")
sig_byterm <- all_res[p_adj_by_term <= 0.01 & term != "(Intercept)"]
fwrite(sig_byterm, out_sig_byterm, compress = "gzip")

message("Done.\n  All results: ", out_all,
        "\n  Global FDR ≤ 0.01: ", out_sig,
        "\n  By-term FDR ≤ 0.01: ", out_sig_byterm)
#!/usr/bin/env Rscript
# ================================================================
# Title: Summarize GLMM Beta Regression Results (Interactive)
# Author: Susan C. Anderson
# Date: 2025-11-05
# Repository: https://github.com/<your-username>/<your-repo>
# Description:
#   Interactive workflow to summarize per-site beta regression
#   results (from glmmTMB) for methylation data in
#   chestnut-crowned babblers.
#   Steps:
#     1. Load complete model output
#     2. Deduplicate by keeping lowest-FDR per site-term
#     3. Pivot wide for intercept + carers3+
#     4. Compute predicted methylation and Δ% between groups
#     5. Filter by FDR ≤ 0.01 and |Δ%| ≥ 25
#     6. Save significant sites (~1,592 expected)
# ================================================================

# -------------------- Libraries --------------------
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
})

# -------------------- 1. Load Full Output --------------------
infile <- "glmm_beta/glmm_results_beta_all.csv.gz"
cat("📂 Loading:", infile, "\n")
full <- fread(infile)

cat("✅ Loaded", nrow(full), "rows across",
    uniqueN(full$site_id), "unique CpG sites\n\n")

# -------------------- 2. Deduplicate per Site-Term --------------------
cat("🧹 Cleaning and keeping lowest-FDR entry per site-term...\n")

full <- full[!is.na(site_id) & !is.na(term)]
before_n <- nrow(full)

setorder(full, site_id, term, p_adj_global)
full_min <- full[, .SD[1], by = .(site_id, term)]

after_n <- nrow(full_min)
cat(sprintf("✅ Reduced from %d → %d unique site-term pairs\n\n",
            before_n, after_n))

# -------------------- 3. Pivot to Wide Format --------------------
cat("🔄 Pivoting to wide format (Intercept + carers3+)...\n")

wide <- full_min %>%
  filter(term %in% c("(Intercept)", "carers3+")) %>%
  select(site_id, term, estimate, p_adj_global, p_adj_by_term) %>%
  pivot_wider(
    names_from = term,
    values_from = c(estimate, p_adj_global, p_adj_by_term),
    names_sep = "_"
  ) %>%
  rename(
    estimate_intercept       = `estimate_(Intercept)`,
    estimate_carers3plus     = `estimate_carers3+`,
    p_adj_global_carers3plus = `p_adj_global_carers3+`,
    p_adj_byterm_carers3plus = `p_adj_by_term_carers3+`
  )

cat("✅ Pivoted data to wide format — one row per CpG site\n\n")

# -------------------- 4. Compute Δβ and Δ% --------------------
cat("📈 Calculating Δβ and Δ% (3+ carers – 2 carers)...\n")

inv_logit <- function(x) exp(x) / (1 + exp(x))

wide <- wide %>%
  mutate(
    p_carers2     = inv_logit(estimate_intercept),
    p_carers3     = inv_logit(estimate_intercept + estimate_carers3plus),
    delta_beta    = p_carers3 - p_carers2,
    delta_percent = 100 * delta_beta
  )

cat("✅ Calculated predicted methylation change for all sites\n\n")

# -------------------- 5. Apply FDR and Effect Size Filters --------------------
cat("🎯 Applying filters (FDR ≤ 0.01 & |Δ%| ≥ 25)...\n")

summary_filtered <- wide %>%
  filter(
    p_adj_global_carers3plus <= 0.01,
    abs(delta_percent) >= 25
  )

cat("✅ Retained", nrow(summary_filtered),
    "sites with FDR ≤ 0.01 and |Δ%| ≥ 25%\n")
cat("Unique sites:", uniqueN(summary_filtered$site_id), "\n\n")

# -------------------- 6. Save Output --------------------
out_dir <- "glmm_beta_summary"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

outfile <- file.path(out_dir, "glmm_beta_significant_sites_FDR01_25pp.csv.gz")
fwrite(summary_filtered, outfile, compress = "gzip")

cat("💾 Saved filtered results to:", outfile, "\n")
cat("✅ Pipeline complete!\n")

# ================================================================
# End of Script
# ================================================================
#annotations for significant sites 

library(data.table)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)



## ---- 1. Load the significant CpG sites (correct file) ----
sig <- fread("glmm_beta_significant_sites_FDR01_25pp.csv.gz")

head(sig)

## ---- 2. Parse chr + start from site_id ----
sig[, c("chr","start") := tstrsplit(site_id, "_", fixed = TRUE)]
sig[, start := as.integer(start)]
sig[, end := start]

# Build CpG GRanges
cpg_gr <- GRanges(seqnames = sig$chr,
                  ranges   = IRanges(sig$start, sig$end),
                  strand   = "*")
mcols(cpg_gr)$site_id <- sig$site_id

## ---- 3. Load transcriptome / gene model ----
gff_file <- "ncbi_dataset/data/GCA_013400735.1/genomic.gff"
txdb <- makeTxDbFromGFF(gff_file, format = "gff")

## ---- 4. Extract gene TSS ----
gene_gr <- genes(txdb)

tss_pos <- ifelse(as.character(strand(gene_gr)) == "+",
                  start(gene_gr), end(gene_gr))

tss_gr <- GRanges(seqnames = seqnames(gene_gr),
                  ranges   = IRanges(tss_pos, tss_pos),
                  strand   = strand(gene_gr))

mcols(tss_gr)$gene_id <- names(gene_gr)

## ---- 5. Compute nearest TSS distance (strand-aware) ----

signed_dist <- function(cpg_pos, tss_pos, tss_strand) {
  ifelse(tss_strand == "+", cpg_pos - tss_pos, tss_pos - cpg_pos)
}
dir_label <- function(d) ifelse(d < 0, "upstream", ifelse(d > 0, "downstream", "at_TSS"))

hits <- GenomicRanges::distanceToNearest(cpg_gr, tss_gr, select="all", ignore.strand = FALSE)
cpg_idx <- queryHits(hits)
tss_idx <- subjectHits(hits)

dist_signed <- signed_dist(start(cpg_gr)[cpg_idx],
                           start(tss_gr)[tss_idx],
                           as.character(strand(tss_gr))[tss_idx])

nearest_dt <- data.table(
  site_id          = mcols(cpg_gr)$site_id[cpg_idx],
  gene_id          = mcols(tss_gr)$gene_id[tss_idx],
  distance_to_TSS  = as.integer(dist_signed),
  distance_direction = dir_label(dist_signed)
)

## ---- 6. Merge back into sig ----
sig_annot_final <- nearest_dt[sig, on="site_id", allow.cartesian = TRUE]

## ---- 7. Save ----
fwrite(sig_annot_final, "significant_fledge_sites_TSS_delta_beta4.csv.gz", compress="gzip")

fwrite(sig_annot_final, "significant_fledge_sites_TSS_annotated_delta_beta4.csv", compress = "none")


uniqueN(sig$site_id)






#.......................................................................................


library(GenomicRanges)
library(GenomicFeatures)
library(data.table)

# ---- 1. Build feature ranges from txdb ----
gene_gr   <- genes(txdb)
tx        <- transcripts(txdb)
utr5_gr   <- unlist(fiveUTRsByTranscript(txdb))
utr3_gr   <- unlist(threeUTRsByTranscript(txdb))
# Gene body defined as full gene span (robust for new GRanges versions)
exon_gr <- genes(txdb)

# TSS (single base per transcript)
tss_raw <- promoters(tx, upstream = 1, downstream = 1)
tss_gr  <- tss_raw
start(tss_gr) <- ifelse(strand(tss_raw)=="+", start(tss_raw), start(tss_raw))
end(tss_gr)   <- start(tss_gr)  # make them 1bp

# TSS window: 300bp upstream to 50bp downstream
tss_window <- promoters(tx, upstream = 300, downstream = 50)

# Promoter window: 2000bp upstream to 200bp downstream
promoter_window <- promoters(tx, upstream = 2000, downstream = 200)

# 10kb upstream of gene body
upstream_10k <- promoters(gene_gr, upstream = 10000, downstream = 0)

# 10kb downstream of gene body
downstream_10k <- flank(gene_gr, width = 10000, start = FALSE)

# ---- 2. Convert your CpGs to GRanges ----
cpg_gr <- GRanges(
  seqnames = sig_annot_final$chr,
  ranges   = IRanges(sig_annot_final$start, sig_annot_final$end),
  strand   = "*"
)
mcols(cpg_gr)$site_id <- sig_annot_final$site_id

# ---- 3. Annotate in PRIORITY order ----
annot <- rep("intergenic", length(cpg_gr))

annot_layer <- function(cpg, feat, label) {
  hits <- findOverlaps(cpg, feat)
  out <- rep(NA_character_, length(cpg))
  out[queryHits(hits)] <- label
  out
}

# 1) TSS
tss_hits <- annot_layer(cpg_gr, tss_window, "TSS")
annot[!is.na(tss_hits)] <- "TSS"

# 2) Promoter (only where not annotated yet)
prom_hits <- annot_layer(cpg_gr, promoter_window, "promoter")
annot[annot == "intergenic" & !is.na(prom_hits)] <- "promoter"

# 3) 5'UTR
utr5_hits <- annot_layer(cpg_gr, utr5_gr, "5UTR")
annot[annot == "intergenic" & !is.na(utr5_hits)] <- "5UTR"

# 4) 3'UTR
utr3_hits <- annot_layer(cpg_gr, utr3_gr, "3UTR")
annot[annot == "intergenic" & !is.na(utr3_hits)] <- "3UTR"

# 5) Gene body (union of exons)
gb_hits <- annot_layer(cpg_gr, exon_gr, "gene_body")
annot[annot == "intergenic" & !is.na(gb_hits)] <- "gene_body"

# 6) Upstream (within 10kb window)
up_hits <- annot_layer(cpg_gr, upstream_10k, "upstream10k")
annot[annot == "intergenic" & !is.na(up_hits)] <- "upstream10k"

# 7) Downstream (within 10kb window)
down_hits <- annot_layer(cpg_gr, downstream_10k, "downstream10k")
annot[annot == "intergenic" & !is.na(down_hits)] <- "downstream10k"

# Save to table
sig_annot_final$annotation <- annot

annot

fwrite(sig_annot_final, "significant_fledge_sites_annotated_delta_beta4.csv", compress="none")


#DMRS....................................................................
library(dmrff)
library(data.table)
library(dplyr)
library(GenomicRanges)
library(readr)


stats <- fread("glmm_results_beta_raw_unadjusted.csv.gz")
head(stats)


stats <- stats %>%
  separate(site_id, into = c("chr", "pos"), sep = "_", remove = FALSE) %>%
  mutate(
    pos = as.numeric(pos),
    chr = gsub("^chr", "", chr)  # remove "chr" prefix if present
  )


anno <- fread("significant_fledge_sites_annotated_beta.csv")

stats <- left_join(
  stats,
  dplyr::select(anno, site_id, chr, start, end, annotation),
  by = "site_id"
)

library(data.table)

# Read the methylation matrix without headers
meth_raw <- fread("dms_methylation_matrix_FDR01_ES25.csv", header = FALSE)

# Extract the first row (which contains sample IDs)
header_row <- as.character(unlist(meth_raw[1, ]))

# Assign that first row as the new column names
setnames(meth_raw, header_row)

# Drop the first row from the dataset
meth <- meth_raw[-1, ]

# Rename the first column to "site_id"
setnames(meth, 1, "site_id")

# Convert the rest of the columns to numeric (they are still characters)
num_cols <- names(meth)[-1]  # everything except site_id
meth[, (num_cols) := lapply(.SD, as.numeric), .SDcols = num_cols]

# Check the structure
str(meth)
head(meth)
fwrite(meth, "dms_methylation_matrix_FDR01_ES25_clean2.csv.gz")

stats <- stats %>%
  mutate(site_id = as.character(site_id)) %>%
  tidyr::separate(site_id, into = c("chr", "pos"), sep = "_", remove = FALSE) %>%
  mutate(pos = as.numeric(pos))


# Keep only CpGs in both datasets
common_sites <- intersect(stats$site_id, meth$site_id)
stats <- stats %>% filter(site_id %in% common_sites)
meth <- meth %>% filter(site_id %in% common_sites)

# Reorder methylation rows to match stats
meth <- meth[match(stats$site_id, meth$site_id), ]
stopifnot(identical(meth$site_id, stats$site_id))

meth_matrix <- as.matrix(meth[, -1, with = FALSE])
rownames(meth_matrix) <- meth$site_id

dmr_results <- dmrff(
  estimate    = stats$estimate,
  se          = stats$se,
  p.value     = stats$p,
  chr         = stats$chr,
  pos         = stats$pos,
  methylation = meth_matrix,
  maxgap      = 2000,     # merge CpGs ≤500 bp apart
  p.cutoff    = 0.05,    # include all CpGs with p ≤ 0.05
  minmem      = TRUE,
  verbose     = TRUE
)

dmr_filtered <- dmr_results %>%
  dplyr::filter(n >= 3, p.adjust < 0.05)

write.csv(dmr_filtered, "dmrff_fledge_final_DMRs2.csv", row.names = FALSE)



library(data.table)
library(GenomicRanges)
library(dplyr)
library(tidyr)

dmr <- fread("dmrff_fledge_final_DMRs2.csv")
head(dmr)

library(GenomicFeatures)
library(GenomicRanges)
library(data.table)
library(dplyr)

# Load the DMR table
dmr <- fread("dmrff_fledge_final_DMRs2.csv")

# Load or rebuild TxDb (only needs to be done once)
gff_file <- "ncbi_dataset/data/GCA_013400735.1/genomic.gff"
txdb <- makeTxDbFromGFF(gff_file, format = "gff")

# Extract genes and TSSs
gene_gr <- genes(txdb)
tx <- transcripts(txdb)

# Define promoter region: 2000 bp upstream to 200 bp downstream of TSS
promoter_gr <- promoters(tx, upstream = 2000, downstream = 200)

dmr_gr <- GRanges(
  seqnames = dmr$chr,
  ranges = IRanges(start = dmr$start, end = dmr$end),
  strand = "*"
)


# Gene body overlaps
hits_gene <- findOverlaps(dmr_gr, gene_gr)
gene_annot <- data.frame(
  chr = as.character(seqnames(dmr_gr)[queryHits(hits_gene)]),
  start = start(dmr_gr)[queryHits(hits_gene)],
  end = end(dmr_gr)[queryHits(hits_gene)],
  gene_id = mcols(gene_gr)$gene_id[subjectHits(hits_gene)],
  feature = "gene_body"
)

# Promoter overlaps
hits_prom <- findOverlaps(dmr_gr, promoter_gr)
prom_annot <- data.frame(
  chr = as.character(seqnames(dmr_gr)[queryHits(hits_prom)]),
  start = start(dmr_gr)[queryHits(hits_prom)],
  end = end(dmr_gr)[queryHits(hits_prom)],
  gene_id = mcols(promoter_gr)$tx_name[subjectHits(hits_prom)],
  feature = "promoter"
)

# Combine annotations
annot_combined <- bind_rows(gene_annot, prom_annot)
annot_combined <- distinct(annot_combined)


nearest_hits <- nearest(dmr_gr, gene_gr)
nearest_annot <- data.frame(
  chr = as.character(seqnames(dmr_gr)),
  start = start(dmr_gr),
  end = end(dmr_gr),
  nearest_gene = mcols(gene_gr)$gene_id[nearest_hits]
)

dmr_annotated <- dmr %>%
  left_join(annot_combined, by = c("chr", "start", "end")) %>%
  left_join(nearest_annot, by = c("chr", "start", "end")) %>%
  mutate(feature = ifelse(is.na(feature), "intergenic", feature))


dmr_annotated

fwrite(dmr_annotated, "dmrff_fledge_final_DMRs_annotated2.csv")
#===========================================================
#END SCRIPT 
