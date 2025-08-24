.PHONY: help install simulate-mr mr-pipeline

help:
	@echo "Targets: install, simulate-mr, mr-pipeline"

install:
	Rscript setup.R

simulate-mr:
	Rscript mr/simulated_mr.R outdir=mr_sim_outputs

# Example with OpenGWAS IDs; replace with real IDs
mr-pipeline:
	Rscript mr/TwoSampleMR_pipeline.R exposure_id=ieu-a-1189 outcome_id=ieu-a-966 outdir=mr_outputs
