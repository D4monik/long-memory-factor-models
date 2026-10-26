# Reproduces all results of the paper. Run from the repository root:
#   Rscript run_all.R
# The Monte Carlo study takes many hours; set RUN_SIMULATION <- FALSE to reproduce
# only the climate application. Step 2 needs the raw GHCN-M files
# (data/raw/README.md) and is skipped when the processed panel is present.

RUN_SIMULATION <- TRUE

if (RUN_SIMULATION) source("scripts/01_simulation.R")
if (!file.exists("data/processed/temperature_panel.csv") ||
    !file.exists("data/processed/station_metadata.csv")) {
  raw_zip <- "data/raw/ghcnm.tavg.v4.0.1.20260330.qcf.zip"
  
  if (file.exists(raw_zip)) {
    unzip(raw_zip, exdir = "data/raw")
  }
  
  source("scripts/02_prepare_data.R")
}
source("scripts/03_levels.R")
source("scripts/04_memory.R")
source("scripts/05_differences.R")
