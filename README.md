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
This portfolio condenses two core summary-statistics frameworks I applied during my PhD: LDSC for polygenic architecture and MR for causal inference. I’ve aimed for dense but precise explanations so a technical reader can assess both capability and judgment.

### Linkage Disequilibrium Score Regression (LDSC)
- **Problem addressed**: Differentiate true polygenic signal from confounding (e.g., cryptic relatedness, population structure) and summarize architecture without individual-level data.
- **LD score**: For SNP j, the LD score is \(l_j = \sum_k r_{jk}^2\) (within a window). Under polygenicity, SNPs with larger \(l_j\) tag more causal sites and have larger expected statistics.
- **Model**: For GWAS sample size N and M SNPs,
  \[ E[\chi_j^2] = a + \frac{N h^2}{M} l_j, \]
  where \(h^2\) is SNP-heritability and intercept \(a\) captures confounding and model misspecification. Weighted regression of \(\chi^2\) on \(l_j\) yields \(\hat h^2\) and \(\hat a\).
- **Genetic correlation (rg)**: For two traits with Z-scores \(z_{1j}, z_{2j}\),
  \[ E[z_{1j} z_{2j}] = c + \frac{\sqrt{N_1 N_2}\, r_g\, h_1 h_2}{M} l_j. \]
  The slope gives \(r_g\) (block jackknife for SEs) while \(c\) diagnoses cross-trait confounding.
- **Partitioned heritability**: Regress on annotation-specific LD scores to estimate enrichment in functional categories; enrichment = proportion of \(h^2\) divided by proportion of SNPs.
- **Key choices**: ancestry-matched reference LD (e.g., 1KG EUR), MAF filters, well-imputed SNP sets (HM3), effective N for case-control, and careful interpretation of intercepts.
- **What I built**: End-to-end pipelines to munge public sumstats, run h2/rg/partitioning, exclude ±500kb around smoking loci (sensitivity to mediation), and parse/visualize results reproducibly.

### Mendelian Randomization (MR)
- **Problem addressed**: Estimate causal effect of exposure X on outcome Y using genetic instruments G, minimizing confounding/reverse causation.
- **Assumptions**: (i) Relevance (G→X), (ii) Independence (G ⟂ confounders), (iii) Exclusion (G affects Y only via X). Violations (horizontal pleiotropy) are explicitly tested/handled.
- **Estimators**:
  - IVW (no-intercept weighted regression): consistent if horizontal pleiotropy is balanced and instruments are strong. A common parametrization is weighted regression of \(\beta_{Yj}\) on \(\beta_{Xj}\) with weights \(w_j = 1/\text{SE}(\beta_{Yj})^2\).
  - MR-Egger: includes intercept \(\alpha\); \(\alpha\)≈0 suggests no directional pleiotropy; slope gives causal effect but with lower power.
  - Weighted median: consistent if ≥50% of weight comes from valid instruments.
  - Robust extensions: MR-PRESSO (outlier removal), RAPS (overdispersion-robust), contamination-mixture.
- **Diagnostics**: clumping (LD), harmonization (allele alignment, palindromics), instrument strength (R²/F-stat), Cochran’s Q (heterogeneity), Egger intercept (pleiotropy), Steiger directionality (X→Y), leave-one-out and funnel plots.
- **Practical caveats**: winner’s curse in discovery GWAS, sample overlap bias in two-sample MR, ancestry mismatch, and multiple testing across traits.
- **What I built**: Scalable univariate and multivariable MR over UKB exposures and TRICL/OncoArray lung cancer outcomes, with full robustness suite and artifacted outputs (tables/plots) for auditability.

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

