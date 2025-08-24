#!/usr/bin/env Rscript

# TwoSampleMR pipeline: robust, documented, data-free by default
# - Supports OpenGWAS IDs or local summary-statistics files
# - Performs clumping, harmonization, IVW / Egger / Weighted Median, heterogeneity,
#   pleiotropy test, Steiger directionality, and leave-one-out
# - Writes tidy outputs and example plots when requested

suppressPackageStartupMessages({
  if (!requireNamespace("TwoSampleMR", quietly = TRUE)) install.packages("TwoSampleMR")
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
  if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
  if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr")
})

library(TwoSampleMR)
library(dplyr)
library(ggplot2)
library(readr)

ensure_dir <- function(path) { if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE) }

# Helpers to read exposure/outcome from file if not using OpenGWAS
read_exposure_file <- function(path) {
  # Adjust column mappings to your file
  TwoSampleMR::read_exposure_data(
    filename = path,
    sep = ",",
    snp_col = "SNP", beta_col = "beta", se_col = "se",
    effect_allele_col = "ea", other_allele_col = "oa",
    eaf_col = "eaf", pval_col = "p",
    phenotype_col = "exposure"
  )
}

read_outcome_file <- function(path, snps) {
  TwoSampleMR::read_outcome_data(
    snps = snps,
    filename = path,
    sep = ",",
    snp_col = "SNP", beta_col = "beta", se_col = "se",
    effect_allele_col = "ea", other_allele_col = "oa",
    eaf_col = "eaf", pval_col = "p",
    phenotype_col = "outcome"
  )
}

run_pipeline <- function(
  exposure_id = NULL,
  outcome_id = NULL,
  exposure_file = NULL,
  outcome_file = NULL,
  pval_threshold = 5e-8,
  clump_r2 = 0.001,
  clump_kb = 10000,
  harmonise_action = 2,
  outdir = "mr_outputs",
  make_plots = TRUE
) {
  ensure_dir(outdir)

  # Instruments (exposure)
  if (!is.null(exposure_id)) {
    message("Extracting instruments from OpenGWAS: ", exposure_id)
    exposure <- TwoSampleMR::extract_instruments(exposure_id, p1 = pval_threshold)
  } else if (!is.null(exposure_file)) {
    message("Reading exposure instruments from file: ", exposure_file)
    exposure <- read_exposure_file(exposure_file)
  } else {
    stop("Provide exposure_id or exposure_file")
  }

  if (nrow(exposure) == 0) stop("No exposure instruments found at the specified threshold")

  # Clumping
  exposure <- TwoSampleMR::clump_data(exposure, clump_r2 = clump_r2, clump_kb = clump_kb)

  # Outcomes
  if (!is.null(outcome_id)) {
    message("Extracting outcome associations for clumped SNPs: ", outcome_id)
    outcome <- TwoSampleMR::extract_outcome_data(snps = exposure$SNP, outcomes = outcome_id)
  } else if (!is.null(outcome_file)) {
    outcome <- read_outcome_file(outcome_file, snps = exposure$SNP)
  } else {
    stop("Provide outcome_id or outcome_file")
  }

  # Harmonise
  harm <- TwoSampleMR::harmonise_data(exposure, outcome, action = harmonise_action)
  write_csv(harm, file.path(outdir, "harmonised_data.csv"))

  # Core MR methods
  methods <- c("mr_ivw", "mr_egger_regression", "mr_weighted_median")
  mr_results <- TwoSampleMR::mr(harm, method_list = methods)
  readr::write_csv(mr_results, file.path(outdir, "mr_results.csv"))

  # Heterogeneity and pleiotropy (Egger intercept)
  het <- TwoSampleMR::mr_heterogeneity(harm, method_list = methods)
  pleio <- TwoSampleMR::mr_pleiotropy_test(harm)
  readr::write_csv(het, file.path(outdir, "heterogeneity.csv"))
  readr::write_csv(pleio, file.path(outdir, "pleiotropy_egger_intercept.csv"))

  # Steiger directionality and LOO
  harm_steiger <- tryCatch(TwoSampleMR::steiger_filtering(harm), error = function(e) harm)
  readr::write_csv(harm_steiger, file.path(outdir, "harmonised_with_steiger.csv"))

  loo <- TwoSampleMR::mr_leaveoneout(harm)
  readr::write_csv(loo, file.path(outdir, "leave_one_out.csv"))

  if (make_plots) {
    p1 <- TwoSampleMR::mr_scatter_plot(mr_results, harm)[[1]]
    p2 <- TwoSampleMR::mr_leaveoneout_plot(loo)[[1]]
    ggplot2::ggsave(file.path(outdir, "mr_scatter.pdf"), p1, width = 6, height = 5)
    ggplot2::ggsave(file.path(outdir, "mr_leave_one_out.pdf"), p2, width = 6, height = 5)
  }

  list(mr = mr_results, heterogeneity = het, pleiotropy = pleio)
}

# CLI usage
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  # Example:
  # Rscript mr/TwoSampleMR_pipeline.R exposure_id=ieu-a-1189 outcome_id=ieu-a-966
  # or
  # Rscript mr/TwoSampleMR_pipeline.R exposure_file=exposure.csv outcome_file=outcome.csv

  kv <- strsplit(args, "=")
  opts <- setNames(lapply(kv, function(x) if (length(x) == 2) x[[2]] else NA_character_),
                   sapply(kv, function(x) x[[1]]))

  run_pipeline(
    exposure_id = opts[["exposure_id"]],
    outcome_id = opts[["outcome_id"]],
    exposure_file = opts[["exposure_file"]],
    outcome_file = opts[["outcome_file"]],
    pval_threshold = if (!is.null(opts[["pval_threshold"]])) as.numeric(opts[["pval_threshold"]]) else 5e-8,
    clump_r2 = if (!is.null(opts[["clump_r2"]])) as.numeric(opts[["clump_r2"]]) else 0.001,
    clump_kb = if (!is.null(opts[["clump_kb"]])) as.numeric(opts[["clump_kb"]]) else 10000,
    outdir = if (!is.null(opts[["outdir"]])) opts[["outdir"]] else "mr_outputs",
    make_plots = if (!is.null(opts[["make_plots"]])) as.logical(opts[["make_plots"]]) else TRUE
  )
}
