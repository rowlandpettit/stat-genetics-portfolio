# Statistical Genetics Portfolio

Curated excerpts from PhD work: linkage disequilibrium score regression (LDSR) and Mendelian randomization (MR). Code is illustrative and data-free.

## Contents
- ldsr/ldsr_2021_reanalysis.Rmd — reproducible analysis outline for LDSR re-analysis (no data)
- mr/Conventional_MR.R — minimal MR workflow example using public packages

## Reproducibility
These examples are designed to run with public summary statistics. Replace file paths with public URLs or local files you control. No proprietary or individual-level data are included.

## Ethics & Privacy
- Excludes VCF/BGEN/PLINK, BAM/CRAM, and private sumstats.
- Only references publicly available resources.

## Background
This portfolio highlights two summary-statistics methods I used during my PhD work:

### Linkage Disequilibrium Score Regression (LDSC)
- **Goal**: Estimate SNP-heritability and genetic correlation using only GWAS summary statistics and LD patterns.
- **Idea**: SNPs with higher LD scores tag more variants and therefore aggregate more true signal; regressing GWAS test statistics on LD scores yields unbiased heritability under polygenicity.
- **Typical outputs**: SNP-heritability (h2), partitioned heritability by functional annotations, and cross-trait genetic correlation (rg).
- **What I did**: Re-analysis pipelines for complex traits (e.g., lung cancer subtypes; smoking-related traits), QC and munging of public summary stats, and replication across multiple trait panels.
- **References**: Bulik-Sullivan et al., 2015 (Nat Genet); LDSC software (`https://github.com/bulik/ldsc`).

### Mendelian Randomization (MR)
- **Goal**: Infer causal effects of an exposure on an outcome using genetic variants as instruments.
- **Assumptions**: Relevance, independence, and exclusion restriction (no horizontal pleiotropy affecting the outcome except via the exposure).
- **Common estimators**: IVW, MR-Egger, weighted median; sensitivity analyses include heterogeneity tests, leave-one-out, and MR-PRESSO (when available).
- **What I did**: Two-sample MR using public instruments for exposures (e.g., smoking traits) and disease outcomes (e.g., lung cancer subtypes), plus robustness checks and multiple-testing control.
- **References**: Davey Smith & Hemani, 2014 (Hum Mol Genet); MR-Base / TwoSampleMR (`https://github.com/MRCIEU/TwoSampleMR`).

## How to use these examples
- The code is purposefully data-free. To reproduce, substitute paths with public GWAS summary statistics and LD references (e.g., 1000 Genomes EUR LD scores from `https://alkesgroup.broadinstitute.org/ldsc/`).
- For MR, obtain instruments (exposure-associated SNPs) and outcome summary stats from public resources (e.g., GWAS Catalog) and load via TwoSampleMR where possible.
- All analyses should be run with appropriate QC, ancestry matching, and multiple-testing adjustments.

## References
- Bulik-Sullivan BK, et al. “LD score regression distinguishes confounding from polygenicity in genome-wide association studies.” Nat Genet (2015).
- Davey Smith G, Hemani G. “Mendelian randomization: genetic anchors for causal inference in epidemiological studies.” Hum Mol Genet (2014).
- LDSC software: `https://github.com/bulik/ldsc`
- TwoSampleMR: `https://github.com/MRCIEU/TwoSampleMR`

---

## Who / What / When / Where / Why
- **Who**: Statistical genetics work performed by Rowland Pettit during PhD research and collaborations.
- **What**: LDSC and MR pipelines for heritability, genetic correlation, and causal inference from GWAS summary statistics.
- **When**: PhD period and subsequent applied projects; methods remain current and widely used.
- **Where**: Analyses executed on Linux/macOS with R and Python; LDSC via Python; MR via R.
- **Why**: Quantify polygenic architecture and assess causal hypotheses while avoiding individual-level data requirements.

## Repository Guide
```
stat-genetics-portfolio/
├── ldsr/
│   ├── LDSC_workflow.Rmd              # End-to-end LDSC workflow (commands + theory)
│   └── ldsr_2021_reanalysis.Rmd       # Example re-analysis outline (no data)
├── mr/
│   ├── TwoSampleMR_pipeline.R         # Robust two-sample MR pipeline (IDs or files)
│   └── simulated_mr.R                 # Synthetic DGP to illustrate estimators
├── setup.R                            # Install R package dependencies
├── LICENSE
├── .gitignore
└── README.md
```

## Reproducible Setup
```r
# From repository root
source("setup.R")
```
- TwoSampleMR may use IEU OpenGWAS; set token if you have one: `Sys.setenv(IEU_OPENAPI_TOKEN = "<token>")`.
- For LDSC, install Python env and obtain LD reference (1KG EUR) from Broad.

## Notes on Data
- This repo includes no individual-level or proprietary data.
- Use public GWAS summary stats and LD references; cite sources appropriately.

