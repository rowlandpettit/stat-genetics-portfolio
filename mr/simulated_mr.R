#!/usr/bin/env Rscript

# Simulated two-sample MR with synthetic summary statistics
# - Generates M instruments with a true causal effect and optional pleiotropy
# - Demonstrates IVW, MR-Egger, and weighted median estimators

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

simulate_summary <- function(
  num_snps = 30,
  true_beta = 0.15,
  se_exposure = 0.02,
  se_outcome = 0.02,
  pleiotropy_sd = 0.00, # set >0 for directional or balanced pleiotropy
  seed = 42
) {
  set.seed(seed)
  snps <- paste0("rs", seq_len(num_snps))
  beta_x <- rnorm(num_snps, mean = 0.1, sd = 0.05)           # instrument strength
  beta_u <- rnorm(num_snps, mean = 0, sd = pleiotropy_sd)     # horizontal pleiotropy
  beta_y <- true_beta * beta_x + beta_u                       # causal + pleiotropy

  data.frame(
    SNP = snps,
    beta.exposure = beta_x + rnorm(num_snps, 0, se_exposure),
    se.exposure = rep(se_exposure, num_snps),
    beta.outcome = beta_y + rnorm(num_snps, 0, se_outcome),
    se.outcome = rep(se_outcome, num_snps),
    effect_allele.exposure = "A",
    other_allele.exposure = "G",
    effect_allele.outcome = "A",
    other_allele.outcome = "G",
    eaf.exposure = 0.4,
    eaf.outcome = 0.4,
    exposure = "SimExposure",
    outcome = "SimOutcome"
  )
}

run_simulated_mr <- function(
  num_snps = 30, true_beta = 0.15, pleiotropy_sd = 0.00,
  se_exposure = 0.02, se_outcome = 0.02, seed = 42,
  outdir = "mr_sim_outputs"
) {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  dat <- simulate_summary(num_snps, true_beta, se_exposure, se_outcome, pleiotropy_sd, seed)

  # TwoSampleMR expects a harmonised format; here alleles already match
  mr_res <- TwoSampleMR::mr(dat, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
  het <- TwoSampleMR::mr_heterogeneity(dat)
  pleio <- TwoSampleMR::mr_pleiotropy_test(dat)

  readr::write_csv(dat, file.path(outdir, "simulated_harmonised.csv"))
  readr::write_csv(mr_res, file.path(outdir, "simulated_mr_results.csv"))
  readr::write_csv(het, file.path(outdir, "simulated_heterogeneity.csv"))
  readr::write_csv(pleio, file.path(outdir, "simulated_pleiotropy_egger.csv"))

  p <- TwoSampleMR::mr_scatter_plot(mr_res, dat)[[1]]
  ggplot2::ggsave(file.path(outdir, "simulated_mr_scatter.pdf"), p, width = 6, height = 5)

  message("True beta = ", true_beta)
  print(mr_res)
}

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  kv <- strsplit(args, "=")
  opts <- setNames(lapply(kv, function(x) if (length(x) == 2) x[[2]] else NA_character_),
                   sapply(kv, function(x) x[[1]]))
  run_simulated_mr(
    num_snps = if (!is.null(opts[["num_snps"]])) as.integer(opts[["num_snps"]]) else 30,
    true_beta = if (!is.null(opts[["true_beta"]])) as.numeric(opts[["true_beta"]]) else 0.15,
    pleiotropy_sd = if (!is.null(opts[["pleiotropy_sd"]])) as.numeric(opts[["pleiotropy_sd"]]) else 0.00,
    seed = if (!is.null(opts[["seed"]])) as.integer(opts[["seed"]]) else 42,
    outdir = if (!is.null(opts[["outdir"]])) opts[["outdir"]] else "mr_sim_outputs"
  )
}
