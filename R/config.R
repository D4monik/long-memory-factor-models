# Settings for the Monte Carlo study (Section 4) and the climate application (Section 5).
# All paths are relative to the repository root.

FIG_DIR <- "output/figures"
TAB_DIR <- "output/tables"
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Monte Carlo study (Section 4) -------------------------------------------

SIM_N       <- c(50, 100, 200, 500, 1000)
SIM_T       <- c(50, 100, 200, 500, 1000)
SIM_D0      <- c(1.25, 0.75, 0)   # memory parameters of the three factors, eq. (28)
SIM_LAMBDA  <- c(0.5, 1.5)        # loadings drawn from U[0.5, 1.5]
SIM_S_MAX   <- 10L
SIM_REPS    <- 1000L
SIM_SIGMA2  <- c(0, 80)           # search interval for the idiosyncratic variance
SIM_TOL     <- 0.01               # stop when the selection frequency is within 0.5 +- SIM_TOL
SIM_SEED    <- 42L                # replication r uses the seed SIM_SEED + r
SIM_CORES   <- max(1L, parallel::detectCores() - 1L)
SIM_RESULTS <- "output/simulation_thresholds.rds"

# ---- Climate application (Section 5) ----------------------------------------

# Raw GHCN-M v4 files (see data/raw/README.md) and processed panel
GHCN_DAT      <- "data/raw/ghcnm.tavg.v4.0.1.20260330.qcf.dat"
GHCN_INV      <- "data/raw/ghcnm.tavg.v4.0.1.20260330.qcf.inv"
PANEL_FILE    <- "data/processed/temperature_panel.csv"
METADATA_FILE <- "data/processed/station_metadata.csv"
LEVELS_RESULTS <- "output/levels.rds"

# Sample (Section 5.1)
SAMPLE_START            <- as.Date("1960-01-01")
SAMPLE_END              <- as.Date("2025-12-01")
MAX_MISSING             <- 12L   # missing months per station
MAX_CONSECUTIVE_MISSING <- 2L
EUROPE_BBOX <- list(lat_min = 34, lat_max = 72, lon_min = -25, lon_max = 45)

# Warming: mean over RECENT minus mean over BASELINE (years, inclusive)
BASELINE <- c(1960, 1990)
RECENT   <- c(2021, 2025)

# Low-frequency transformation: keep cosine components with periods > 12.5 months
LOWPASS_CUTOFF <- 12.5

# Factor-number selection (Alessi, Barigozzi and Capasso, 2010)
ABC_S_MAX <- 20L
ABC_C     <- seq_len(30 * 500) / 500   # penalty multipliers c = 0.002, 0.004, ..., 30
ABC_REPS  <- 100L                      # random subpanels per subpanel size
ABC_SEED  <- 12345L

# Number of factors in levels and after first differencing
K_LEVELS      <- 2L
K_DIFFERENCES <- 1L

# Memory estimation (Mueller and Watson, 2008)
MW_Q      <- 13L
MW_D_GRID <- seq(0, 4, by = 0.01)
