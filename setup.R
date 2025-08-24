# Setup script for statistical genetics examples

pkgs <- c(
  "TwoSampleMR",
  "dplyr",
  "ggplot2",
  "readr"
)

missing <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(missing)) install.packages(missing)

message("Optional: set IEU OpenGWAS token if available: \n  Sys.setenv(IEU_OPENAPI_TOKEN = '<token>')")
