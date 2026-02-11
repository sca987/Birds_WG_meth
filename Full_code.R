

###############################################################################
# 1) filter_fledge_methylation.R
###############################################################################
#!/usr/bin/env Rscript
# ================================================================
# Title: Fledge Methylation Filtering Pipeline
# Author: Susan C. Anderson
# Date: 2026-02-11
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
#   all_fledge_filtered_for_model_2025_v1_with_xfoster_brood_careractual.csv.gz
# ================================================================

suppressPackageStartupMessages(library(data.table))

fledge_all <- fread("all_fledge_merged_long_20251019_v3_with_missing.csv.gz")
cat("Loaded dataset with", nrow(fledge_all), "rows and",
    uniqueN(fledge_all$sample), "individuals\n")

bad_ids <- c(79596, 79686, 79690, 62603, 62644, 79566, 79692, 79688)
fledge_all <- fledge_all[!sample %in% bad_ids]
cat("Removed bad individuals. Remaining:", uniqueN(fledge_all$sample), "samples\n")

fledge_10x <- fledge_all[!is.na(total_coverage) & total_coverage >= 10]
cat("After ≥10x filter:", nrow(fledge_10x), "rows\n")

fledge_99 <- fledge_10x[, {
  thresh <- quantile(total_coverage, probs = 0.999, na.rm = TRUE)
  .SD[total_coverage <= thresh]
}, by = sample]
cat("After per-sample 99.9% filter:", nrow(fledge_99), "rows\n")

site_group_counts <- fledge_99[, .(n_individuals = uniqueN(sample)),
                               by = .(chr, pos, carers)]

group_sizes <- fledge_99[, .(total_individuals = uniqueN(sample)), by = carers]
print(group_sizes)

min_thresholds <- group_sizes[, .(
  carers,
  min_individuals = ceiling(0.60 * total_individuals)
)]
print(min_thresholds)

site_group_counts <- merge(site_group_counts, group_sizes, by = "carers")
site_group_counts[, prop_present := n_individuals / total_individuals]
site_group_counts <- merge(site_group_counts, min_thresholds, by = "carers")

site_group_pass <- site_group_counts[
  prop_present >= 0.60 & n_individuals >= min_individuals
]

site_group_pass_counts <- site_group_pass[, .N, by = .(chr, pos)]
sites_to_keep <- site_group_pass_counts[N == uniqueN(fledge_99$carers)]

fledge_60 <- fledge_99[sites_to_keep, on = c("chr", "pos")]
cat("After 60% + min-individual filter:",
    uniqueN(fledge_60[, .(chr, pos)]), "unique sites retained\n")

site_var <- fledge_60[, .(var_meth = var(perc_meth, na.rm = TRUE)), by = .(chr, pos)]
threshold_var <- quantile(site_var$var_meth, probs = 0.10, na.rm = TRUE)
sites_high_var <- site_var[var_meth > threshold_var]

fledge_final <- fledge_60[sites_high_var, on = c("chr", "pos")]
cat("After variance filter:",
    uniqueN(fledge_final[, .(chr, pos)]), "unique sites retained\n")

out_file <- "all_fledge_filtered_for_model_2025_v1_with_xfoster_brood_careractual.csv.gz"
fwrite(fledge_final, out_file, compress = "gzip")

cat("\nFinal dataset saved:", out_file, "\n")
cat("Final unique sites:", uniqueN(fledge_final[, .(chr, pos)]), "\n")
cat("Final individuals:", uniqueN(fledge_final$sample), "\n")


###############################################################################
# 2) run_glmm_beta_carers_xfoster_age.slurm
###############################################################################
#!/bin/bash
#SBATCH --job-name=glmm_beta_car_xf_age
#SBATCH --output=logs/glmm_beta_car_xf_age_%j.out
#SBATCH --error=logs/glmm_beta_car_xf_age_%j.err
#SBATCH --partition=himem
#SBATCH --cpus-per-task=48
#SBATCH --mem=1400G
#SBATCH --time=48:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=susan.c.anderson@coyotes.usd.edu

module load R/4.4.2

export FLEDGE_INPUT="all_fledge_filtered_for_model_2025_v1_with_xfoster_brood_careractual.csv.gz"
export OUT_DIR="glmm_beta_out_carers_xfoster_age"
export OUT_ALL="glmm_results_beta_all_carers_xfoster_age.csv.gz"
export OUT_SIG="significant_glmm_sites_beta_global_FDR01_carers_xfoster_age.csv.gz"

# Match CPUs requested above
export N_CORES=48
export FDR_ALPHA=0.05
export SITE_LIMIT=0

mkdir -p logs

echo "Job started on $(hostname) at $(date)"
echo "Running GLMM Beta model using $N_CORES cores"

Rscript br_xf_carers_age.R

echo "Job finished at $(date)"


###############################################################################
# 3) br_xf_carers_age.R
###############################################################################
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(glmmTMB)
  library(parallel)
})

infile     <- Sys.getenv("FLEDGE_INPUT",
                         "all_fledge_filtered_for_model_2025_v1_with_xfoster_brood_careractual.csv.gz")
out_dir    <- Sys.getenv("OUT_DIR", "glmm_beta_out_carers_xfoster_age")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_all    <- file.path(out_dir, Sys.getenv("OUT_ALL", "glmm_results_beta_all_carers_xfoster_age.csv.gz"))
out_sig    <- file.path(out_dir, Sys.getenv("OUT_SIG", "significant_glmm_sites_beta_global_FDR01_carers_xfoster_age.csv.gz"))
parts_dir  <- file.path(out_dir, "parts")
dir.create(parts_dir, showWarnings = FALSE, recursive = TRUE)

n_cores    <- as.integer(Sys.getenv("N_CORES", "48"))
site_limit <- as.integer(Sys.getenv("SITE_LIMIT", "0"))

dt <- fread(infile)

dt[, carers := as.factor(carers)]
if ("2" %in% levels(dt$carers)) dt[, carers := relevel(carers, ref = "2")]

stopifnot("xfoster" %in% names(dt))
dt[, xfoster := as.factor(xfoster)]
if (all(c("n","y") %in% levels(dt$xfoster)))  dt[, xfoster := relevel(xfoster, ref = "n")]
if (all(c("no","yes") %in% levels(dt$xfoster))) dt[, xfoster := relevel(xfoster, ref = "no")]

dt[, natal_group := as.factor(natal_group)]
dt[, age_scaled  := as.numeric(scale(age))]

dt[, perc_meth_beta := fifelse(perc_meth <= 0,    0.0001,
                               fifelse(perc_meth >= 100, 0.9999, perc_meth/100))]

dt[, site_id := paste(chr, pos, sep = "_")]

sites <- unique(dt$site_id)
if (site_limit > 0) sites <- sites[seq_len(min(site_limit, length(sites)))]
message("Total sites to process: ", length(sites))

shards <- split(sites, cut(seq_along(sites), breaks = n_cores, labels = FALSE))

fit_one_site <- function(s) {
  sub <- dt[site_id == s]
  sub <- sub[complete.cases(perc_meth_beta, carers, xfoster, age_scaled, natal_group)]
  
  if (nrow(sub) < 10L ||
      length(unique(sub$carers)) < 2L ||
      length(unique(sub$xfoster)) < 2L ||
      var(sub$perc_meth_beta) == 0) return(NULL)
  
  fit <- try(
    glmmTMB(
      perc_meth_beta ~ carers + xfoster + age_scaled + (1 | natal_group),
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

if (.Platform$OS.type == "windows") {
  part_files <- lapply(seq_along(shards), run_shard)
} else {
  part_files <- mclapply(seq_along(shards), run_shard, mc.cores = n_cores)
}

parts <- list.files(parts_dir, pattern = "^part_\\d{3}\\.csv\\.gz$", full.names = TRUE)
if (!length(parts)) stop("No part files found.")

message("Combining ", length(parts), " shards ...")
all_res <- rbindlist(lapply(parts, fread), use.names = TRUE, fill = TRUE)
all_res <- all_res[!is.na(p) & !is.na(term) & !is.na(site_id)]

out_raw <- file.path(out_dir, "glmm_results_beta_raw_unadjusted_carers_xfoster_age.csv.gz")
fwrite(all_res, out_raw, compress = "gzip")
message("Saved raw model output: ", out_raw)

all_res[, p_adj_by_term := p.adjust(p, method = "BH"), by = term]
all_res[, p_adj_global  := p.adjust(p, method = "BH")]

fwrite(all_res, out_all, compress = "gzip")

sig_global <- all_res[p_adj_global <= 0.01 & term != "(Intercept)"]
fwrite(sig_global, out_sig, compress = "gzip")

out_sig_byterm <- file.path(out_dir, "significant_glmm_sites_beta_byterm_FDR01_carers_xfoster_age.csv.gz")
sig_byterm <- all_res[p_adj_by_term <= 0.01 & term != "(Intercept)"]
fwrite(sig_byterm, out_sig_byterm, compress = "gzip")

message("Done.\n  All results: ", out_all,
        "\n  Global FDR ≤ 0.01: ", out_sig,
        "\n  By-term FDR ≤ 0.01: ", out_sig_byterm)


###############################################################################
# 4) summarize_carers_byterm.R  (interactive post-filter)
###############################################################################
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
})

# NOTE: this is the model output produced above (gz)
infile <- "glmm_results_beta_all_carers_xfoster_age.csv.gz"
cat("Loading:", infile, "\n")
full <- fread(infile)

cat("Loaded", nrow(full), "rows across",
    uniqueN(full$site_id), "unique CpG sites\n\n")

full <- full[!is.na(site_id) & !is.na(term)]
before_n <- nrow(full)

# Deduplicate using BY-TERM FDR
setorder(full, site_id, term, p_adj_by_term)
full_min <- full[, .SD[1], by = .(site_id, term)]

after_n <- nrow(full_min)
cat(sprintf("Reduced from %d → %d unique site-term pairs\n\n",
            before_n, after_n))

wide <- full_min %>%
  filter(term %in% c("(Intercept)", "carers3+")) %>%
  select(site_id, term, estimate, p_adj_by_term) %>%
  pivot_wider(
    names_from = term,
    values_from = c(estimate, p_adj_by_term),
    names_sep = "_"
  ) %>%
  rename(
    estimate_intercept       = `estimate_(Intercept)`,
    estimate_carers3plus     = `estimate_carers3+`,
    p_adj_byterm_carers3plus = `p_adj_by_term_carers3+`
  )

inv_logit <- function(x) exp(x) / (1 + exp(x))

wide <- wide %>%
  mutate(
    pred_carers2     = inv_logit(estimate_intercept),
    pred_carers3plus = inv_logit(estimate_intercept + estimate_carers3plus),
    delta_beta       = pred_carers3plus - pred_carers2,
    delta_percent    = 100 * delta_beta
  )

summary_filtered <- wide %>%
  filter(
    p_adj_byterm_carers3plus <= 0.01,
    abs(delta_percent) >= 25
  ) %>%
  arrange(desc(abs(delta_percent)))

cat("Retained", nrow(summary_filtered),
    "sites with BY-TERM FDR ≤ 0.01 and |Δ%| ≥ 25%\n")
cat("Unique sites:", uniqueN(summary_filtered$site_id), "\n\n")

outfile <- "glmm_beta_significant_sites_CARERS_bytermFDR01_25pp.csv.gz"
fwrite(summary_filtered, outfile, compress = "gzip")
cat("Saved filtered results to:", outfile, "\n")


###############################################################################
# 5) annotate_and_dmrff_carers3plus.R
###############################################################################
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(rtracklayer)
  library(dmrff)
  library(readr)
})

# ---- Significant sites (must match summarize output name) ----
sig <- fread("glmm_beta_significant_sites_CARERS_bytermFDR01_25pp.csv.gz")
cat("Loaded significant CpGs:", uniqueN(sig$site_id), "\n")
print(head(sig))

sig[, c("chr","start") := tstrsplit(site_id, "_", fixed = TRUE)]
sig[, start := as.integer(start)]
sig[, end := start]

cpg_gr <- GRanges(seqnames = sig$chr,
                  ranges   = IRanges(sig$start, sig$end),
                  strand   = "*")
mcols(cpg_gr)$site_id <- sig$site_id

gff_file <- "ncbi_dataset/data/GCA_013400735.1/genomic.gff"
txdb <- makeTxDbFromGFF(gff_file, format = "gff")

gene_gr <- genes(txdb)
tss_pos <- ifelse(as.character(strand(gene_gr)) == "+",
                  start(gene_gr), end(gene_gr))

tss_gr <- GRanges(seqnames = seqnames(gene_gr),
                  ranges   = IRanges(tss_pos, tss_pos),
                  strand   = strand(gene_gr))
mcols(tss_gr)$gene_id <- names(gene_gr)

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
  site_id            = mcols(cpg_gr)$site_id[cpg_idx],
  gene_id            = mcols(tss_gr)$gene_id[tss_idx],
  distance_to_TSS    = as.integer(dist_signed),
  distance_direction = dir_label(dist_signed)
)

sig_annot_final <- nearest_dt[sig, on="site_id", allow.cartesian = TRUE]
fwrite(sig_annot_final, "significant_fledge_sites_TSS_annotated.csv", compress = "none")
cat("Saved: significant_fledge_sites_TSS_annotated.csv\n")

# Feature annotation
tx        <- transcripts(txdb)
utr5_gr   <- unlist(fiveUTRsByTranscript(txdb))
utr3_gr   <- unlist(threeUTRsByTranscript(txdb))
exon_gr   <- genes(txdb)

tss_window       <- promoters(tx, upstream = 300, downstream = 50)
promoter_window  <- promoters(tx, upstream = 2000, downstream = 200)
upstream_10k     <- promoters(gene_gr, upstream = 10000, downstream = 0)
downstream_10k   <- flank(gene_gr, width = 10000, start = FALSE)

cpg_gr2 <- GRanges(
  seqnames = sig_annot_final$chr,
  ranges   = IRanges(sig_annot_final$start, sig_annot_final$end),
  strand   = "*"
)
mcols(cpg_gr2)$site_id <- sig_annot_final$site_id

annot <- rep("intergenic", length(cpg_gr2))

annot_layer <- function(cpg, feat, label) {
  hits <- findOverlaps(cpg, feat)
  out <- rep(NA_character_, length(cpg))
  out[queryHits(hits)] <- label
  out
}

tss_hits <- annot_layer(cpg_gr2, tss_window, "TSS")
annot[!is.na(tss_hits)] <- "TSS"

prom_hits <- annot_layer(cpg_gr2, promoter_window, "promoter")
annot[annot == "intergenic" & !is.na(prom_hits)] <- "promoter"

utr5_hits <- annot_layer(cpg_gr2, utr5_gr, "5UTR")
annot[annot == "intergenic" & !is.na(utr5_hits)] <- "5UTR"

utr3_hits <- annot_layer(cpg_gr2, utr3_gr, "3UTR")
annot[annot == "intergenic" & !is.na(utr3_hits)] <- "3UTR"

gb_hits <- annot_layer(cpg_gr2, exon_gr, "gene_body")
annot[annot == "intergenic" & !is.na(gb_hits)] <- "gene_body"

up_hits <- annot_layer(cpg_gr2, upstream_10k, "upstream10k")
annot[annot == "intergenic" & !is.na(up_hits)] <- "upstream10k"

down_hits <- annot_layer(cpg_gr2, downstream_10k, "downstream10k")
annot[annot == "intergenic" & !is.na(down_hits)] <- "downstream10k"

sig_annot_final$annotation <- annot
fwrite(sig_annot_final, "significant_fledge_sites_annotated.csv", compress = "none")
cat("Saved: significant_fledge_sites_annotated.csv\n")

# ---------------- DMRFF (carers3+ ONLY) ----------------
stats <- fread("glmm_results_beta_all_carers_xfoster_age.csv.gz")

# carers3+ ONLY
stats <- stats[term == "carers3+"]
stopifnot(nrow(stats) > 0)

stats <- stats %>%
  mutate(site_id = as.character(site_id)) %>%
  separate(site_id, into = c("chr", "pos"), sep = "_", remove = FALSE) %>%
  mutate(
    pos = as.numeric(pos),
    chr = gsub("^chr", "", chr)
  )

anno <- fread("significant_fledge_sites_annotated.csv")
stats <- left_join(
  stats,
  dplyr::select(anno, site_id, chr, start, end, annotation),
  by = "site_id"
)

meth_raw <- fread("dms_methylation_matrix_FDR01_ES25_xf_carers.csv", header = FALSE)
header_row <- as.character(unlist(meth_raw[1, ]))
setnames(meth_raw, header_row)
meth <- meth_raw[-1, ]
setnames(meth, 1, "site_id")

num_cols <- names(meth)[-1]
meth[, (num_cols) := lapply(.SD, as.numeric), .SDcols = num_cols]

fwrite(meth, "dms_methylation_matrix_FDR01_ES25_xf_carers_clean2.csv.gz", compress = "gzip")
cat("Saved: dms_methylation_matrix_FDR01_ES25_xf_carers_clean2.csv.gz\n")

common_sites <- intersect(stats$site_id, meth$site_id)
stats <- stats %>% filter(site_id %in% common_sites)
meth  <- meth  %>% filter(site_id %in% common_sites)

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
  maxgap      = 2000,
  p.cutoff    = 0.05,
  minmem      = TRUE,
  verbose     = TRUE
)

dmr_filtered <- dmr_results %>%
  dplyr::filter(n >= 3, p.adjust < 0.05)

write.csv(dmr_filtered, "dmrff_fledge_final_DMRs2_carers3plus_only.csv", row.names = FALSE)
cat("Saved: dmrff_fledge_final_DMRs2_carers3plus_only.csv\n")

dmr <- fread("dmrff_fledge_final_DMRs2_carers3plus_only.csv")

promoter_gr <- promoters(transcripts(txdb), upstream = 2000, downstream = 200)

dmr_gr <- GRanges(
  seqnames = dmr$chr,
  ranges = IRanges(start = dmr$start, end = dmr$end),
  strand = "*"
)

hits_gene <- findOverlaps(dmr_gr, genes(txdb))
gene_annot <- data.frame(
  chr = as.character(seqnames(dmr_gr)[queryHits(hits_gene)]),
  start = start(dmr_gr)[queryHits(hits_gene)],
  end = end(dmr_gr)[queryHits(hits_gene)],
  gene_id = mcols(genes(txdb))$gene_id[subjectHits(hits_gene)],
  feature = "gene_body"
)

hits_prom <- findOverlaps(dmr_gr, promoter_gr)
prom_annot <- data.frame(
  chr = as.character(seqnames(dmr_gr)[queryHits(hits_prom)]),
  start = start(dmr_gr)[queryHits(hits_prom)],
  end = end(dmr_gr)[queryHits(hits_prom)],
  gene_id = mcols(promoter_gr)$tx_name[subjectHits(hits_prom)],
  feature = "promoter"
)

annot_combined <- bind_rows(gene_annot, prom_annot) %>% distinct()

nearest_hits <- nearest(dmr_gr, genes(txdb))
nearest_annot <- data.frame(
  chr = as.character(seqnames(dmr_gr)),
  start = start(dmr_gr),
  end = end(dmr_gr),
  nearest_gene = mcols(genes(txdb))$gene_id[nearest_hits]
)

dmr_annotated <- dmr %>%
  left_join(annot_combined, by = c("chr", "start", "end")) %>%
  left_join(nearest_annot, by = c("chr", "start", "end")) %>%
  mutate(feature = ifelse(is.na(feature), "intergenic", feature))

fwrite(dmr_annotated, "dmrff_fledge_final_DMRs_annotated2_carers3plus_only.csv.gz", compress = "gzip")
cat("Saved: dmrff_fledge_final_DMRs_annotated2_carers3plus_only.csv.gz\n")

cat("\n=== DONE ===\n")
