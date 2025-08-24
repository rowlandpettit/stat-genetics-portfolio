# Statistical Genetics Portfolio

Curated excerpts from PhD work: linkage disequilibrium score regression (LDSR) and Mendelian randomization (MR). Code is illustrative and data-free.

## Contents
- ldsr/LDSC_workflow.Rmd — end-to-end LDSC workflow (theory + commands + QC)
- ldsr/ldsr_2021_reanalysis_refactored.Rmd — refactored LDSR re-analysis notebook
- mr/ukbb_tricl_mr_refactored.Rmd — refactored UKB↔TRICL MR workflow (full diagnostics)
- mr/TwoSampleMR_pipeline.R — robust MR pipeline (IDs or local files)
- mr/simulated_mr.R — synthetic MR to illustrate estimators

## Reproducibility
These examples are designed to run with public summary statistics. Replace file paths with public URLs or local files you control. No proprietary or individual-level data are included.

## Ethics & Privacy
- Excludes VCF/BGEN/PLINK, BAM/CRAM, and private sumstats.
- Only references publicly available resources.

## Background
This portfolio condenses two core summary-statistics frameworks I applied during my PhD: LDSC to describe polygenic architecture and MR to probe causal effects. The goal here is clarity: explain what each method does, why it’s useful, and how I used it in practice.

### Linkage Disequilibrium Score Regression (LDSC)
- **What it is**: A regression of GWAS signal on local LD to separate true polygenic signal from confounding. It estimates SNP-heritability (how much variance is captured by common SNPs) and the genetic correlation between traits, using only summary statistics plus an external LD reference.
- **Why it works**: In polygenic traits, SNPs that sit in regions of high LD tend to “tag” more causal variation and therefore have larger expected test statistics. Regressing the observed statistics on LD quantifies how much signal scales with LD (heritability) versus a baseline that does not (confounding).
- **What you get**:
  - SNP-heritability with an intercept that helps diagnose inflation not explained by polygenicity.
  - Cross-trait genetic correlation (shared polygenic architecture) with uncertainty from block jackknife.
  - Partitioned heritability: enrichment of h2 within functional annotations (e.g., conserved elements, tissue marks).
- **How I used it**: Clean and standardize public GWAS, run h2/rg at scale, and compare partitioned heritability before/after excluding windows around smoking-associated loci to separate mediation from pleiotropy in lung cancer analyses. Emphasis on ancestry-matched LD references, HM3 SNP sets, and effective sample size for binary traits.

### Mendelian Randomization (MR)
- **What it is**: An instrumental-variables framework using genetic variants that influence an exposure to estimate its causal effect on an outcome from observational data.
- **Why it works**: Random allocation of alleles at conception makes genetic instruments largely independent of environmental confounders. Provided instruments are relevant and do not affect the outcome through other pathways, the genetic association with the outcome scales with the genetic association with the exposure.
- **What you get**:
  - Causal effect estimates via complementary estimators (IVW, MR-Egger, weighted median) that trade bias and robustness.
  - Diagnostics for instrument validity (F-statistics), heterogeneity (Cochran’s Q), directional pleiotropy (Egger intercept), directionality (Steiger), and sensitivity (leave-one-out, funnel plots).
  - Robust extensions (MR-PRESSO, RAPS, contamination mixture) to down-weight outliers or model overdispersion.
- **How I used it**: Two-sample MR linking UKB exposures (e.g., smoking intensity, education, income) to TRICL/OncoArray lung cancer subtypes, with full robustness checks and careful interpretation around sample overlap and multiple testing.

## How to use these examples
- The code is purposefully data-free. To reproduce, substitute paths with public GWAS summary statistics and LD references (e.g., 1000 Genomes EUR LD scores from `https://alkesgroup.broadinstitute.org/ldsc/`).
- For MR, obtain instruments (exposure-associated SNPs) and outcome summary stats from public resources (e.g., GWAS Catalog) and load via TwoSampleMR where possible.
- All analyses should be run with appropriate QC, ancestry matching, and multiple-testing adjustments.

## References
- Bulik-Sullivan BK, et al. “LD score regression distinguishes confounding from polygenicity in genome-wide association studies.” Nat Genet (2015).
- Finucane HK, et al. “Partitioning heritability by functional annotation using genome-wide association summary statistics.” Nat Genet (2015).
- Gazal S, et al. “Linkage disequilibrium–dependent architecture of human complex traits.” Nat Genet (2017).
- Davey Smith G, Hemani G. “Mendelian randomization: genetic anchors for causal inference in epidemiological studies.” Hum Mol Genet (2014).
- Bowden J, et al. “Mendelian randomization with invalid instruments: effect estimation and bias detection via Egger regression.” Int J Epidemiol (2015).
- Bowden J, et al. “Consistent estimation in Mendelian randomization with some invalid instruments using a weighted median estimator.” Genet Epidemiol (2016).
- Verbanck M, et al. “Detection of widespread horizontal pleiotropy in causal relationships inferred from MR.” Nat Genet (2018). (MR-PRESSO)
- Zhao Q, et al. “Statistical inference in two-sample summary-data MR using robust adjusted profile score.” Ann Stat (2020). (MR-RAPS)
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
│   ├── LDSC_workflow.Rmd                      # End-to-end LDSC workflow
│   └── ldsr_2021_reanalysis_refactored.Rmd    # Refactored re-analysis
├── mr/
│   ├── ukbb_tricl_mr_refactored.Rmd            # Full MR workflow (UKB↔TRICL)
│   ├── TwoSampleMR_pipeline.R                  # Robust two-sample MR pipeline
│   └── simulated_mr.R                           # Synthetic MR example
├── setup.R                            # Install R package dependencies
├── LICENSE
├── .gitignore
└── README.md
```

## Technical overview of thesis work (selected)
- Polygenic architecture characterization via LDSC (heritability, genetic correlation) at population scale; emphasis on careful QC and intercept interpretation.
- Partitioned heritability analyses to localize signal to functional annotations and test sensitivity to removal of smoking-associated windows in lung cancer.
- Two-sample MR pipelines connecting behavioral and socioeconomic exposures to cancer outcomes, with full robustness suite (IVW/Egger/median + PRESSO/RAPS) and directionality checks.
- Reproducibility practices: standardized munging, ancestry-matched LD references, effective N for binary traits, artifacted outputs (tables, plots) for auditability.

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

