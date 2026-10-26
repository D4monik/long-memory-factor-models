# Section 5.1: balanced panel of monthly mean temperatures at European GHCN-M stations.
# Input : GHCN-M v4 data and inventory files in data/raw/ (see data/raw/README.md)
# Output: data/processed/temperature_panel.csv, data/processed/station_metadata.csv

source("R/config.R")
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr); library(tidyr)
})

## Monthly records: station id, year, element, and 12 blocks (value, three flags)
month_cols <- sprintf("m%02d", 1:12)
widths <- fwf_widths(
  c(11, 4, 4, rep(c(5, 1, 1, 1), 12)),
  c("id", "year", "elem",
    as.vector(t(outer(month_cols, c("", "_mflag", "_qflag", "_sflag"), paste0))))
)
ghcn <- read_fwf(GHCN_DAT, widths, col_types = cols(.default = col_character())) |>
  filter(elem == "TAVG") |>
  mutate(id = trimws(id), year = as.integer(year))

# values in hundredths of a degree Celsius; -9999 denotes a missing value
for (m in month_cols) {
  v <- as.integer(ghcn[[m]])
  ghcn[[m]] <- ifelse(v == -9999L, NA_real_, v / 100)
  ghcn[[paste0(m, "_qflag")]] <- trimws(ghcn[[paste0(m, "_qflag")]])
}

## Station inventory: id, latitude, longitude, elevation, name
inv <- tibble(raw = read_lines(GHCN_INV))
fields <- str_match(
  substr(inv$raw, 12, nchar(inv$raw)),
  "^\\s*(-?\\d+(?:\\.\\d+)?)\\s+(-?\\d+(?:\\.\\d+)?)\\s+(-?\\d+(?:\\.\\d+)?)\\s*([A-Za-z0-9]{0,2})?\\s*(.*)$"
)
stations <- tibble(
  id = trimws(substr(inv$raw, 1, 11)),
  country_code = substr(id, 1, 2),
  lat = as.numeric(fields[, 2]),
  lon = as.numeric(fields[, 3]),
  elev = as.numeric(fields[, 4]),
  name = str_trim(fields[, 6])
)
europe <- stations |>
  filter(lat >= EUROPE_BBOX$lat_min, lat <= EUROPE_BBOX$lat_max,
         lon >= EUROPE_BBOX$lon_min, lon <= EUROPE_BBOX$lon_max)

## Station x month matrix; values with a quality flag other than blank or "A" are dropped
wide <- ghcn |>
  filter(id %in% europe$id) |>
  rename_with(~ paste0(.x, "_tavg"), all_of(month_cols)) |>
  pivot_longer(cols = c(paste0(month_cols, "_tavg"), paste0(month_cols, "_qflag")),
               names_to = c("mon", ".value"), names_pattern = "^(m\\d{2})_(tavg|qflag)$") |>
  mutate(date = as.Date(sprintf("%04d-%s-01", year, substr(mon, 2, 3))),
         tavg = ifelse(!is.na(qflag) & !(qflag %in% c("", "A")), NA_real_, as.numeric(tavg))) |>
  select(id, date, tavg) |>
  pivot_wider(id_cols = id, names_from = date, values_from = tavg)

dates <- sort(as.Date(colnames(wide)[-1]))
dates <- dates[dates >= SAMPLE_START & dates <= SAMPLE_END]
X <- as.matrix(wide[, as.character(dates)])
rownames(X) <- wide$id

## Stations observed at both endpoints with few (consecutive) missing months
max_run_missing <- function(x) {
  r <- rle(is.na(x))
  if (any(r$values)) max(r$lengths[r$values]) else 0L
}
keep <- rowSums(is.na(X)) <= MAX_MISSING &
  apply(X, 1, max_run_missing) <= MAX_CONSECUTIVE_MISSING &
  !is.na(X[, 1]) & !is.na(X[, ncol(X)])
X <- X[keep, , drop = FALSE]

# linear interpolation of the remaining internal gaps
X <- t(apply(X, 1, function(x) approx(which(!is.na(x)), x[!is.na(x)], seq_along(x))$y))
stopifnot(!anyNA(X))

## Export: months in rows, stations in columns; metadata in the same station order
dir.create(dirname(PANEL_FILE), recursive = TRUE, showWarnings = FALSE)
write.csv(data.frame(Date = dates, t(X), check.names = FALSE), PANEL_FILE, row.names = FALSE)
meta <- europe[match(rownames(X), europe$id), c("id", "country_code", "name", "lat", "lon", "elev")]
write.csv(meta, METADATA_FILE, row.names = FALSE)

message(sprintf("Balanced panel: N = %d stations, T = %d months", nrow(X), ncol(X)))
