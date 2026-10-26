# Section 5.2: low-frequency panel, number of factors, and the two-factor model in levels.
# Output: Figures 1-6, station and regional warming tables, output/levels.rds

source("R/config.R")
source("R/estimation.R")
source("R/factor_selection.R")
source("R/plots.R")

panel <- read.csv(PANEL_FILE, check.names = FALSE)
dates <- as.Date(panel$Date)
years <- as.integer(format(dates, "%Y"))
Y <- as.matrix(panel[, -1])
meta <- read.csv(METADATA_FILE, colClasses = c(id = "character"))
meta <- meta[match(colnames(Y), meta$id), ]
stopifnot(identical(meta$id, colnames(Y)))
meta$region <- region_of(meta$lat, meta$lon)

## Figure 1: monthly temperatures and low-frequency components
Y_low <- lowpass(Y, LOWPASS_CUTOFF)
save_pdf(plot_panel(Y, dates), "fig1a_temperatures.pdf", 8.5, 5)
save_pdf(plot_panel(Y_low, dates), "fig1b_lowfreq_components.pdf", 8.5, 5)

## Figure 2: number of factors
abc <- abc_selection(Y_low, ABC_S_MAX, ABC_C, ABC_REPS, ABC_SEED)
write.csv(abc$stable, file.path(TAB_DIR, "abc_stable_intervals_levels.csv"), row.names = FALSE)
save_pdf(plot_abc(abc$path, ABC_S_MAX), "fig2a_abc_levels.pdf", 7.5, 5)

fit <- fit_pca(Y_low, K_LEVELS)
save_pdf(plot_eigenvalues(fit$eigenvalues), "fig2b_eigenvalues_levels.pdf", 7.5, 5)

## Figure 3: estimated factors (separate panels)
Z <- standardize_factors(fit$factors, years)
save_pdf(plot_factors(Z[, 1, drop = FALSE], dates), "fig3a_factor1.pdf", 6, 4.5)
save_pdf(plot_factors(Z[, 2, drop = FALSE], dates), "fig3b_factor2.pdf", 6, 4.5)

## Figure 4: common component and factor contributions
xbar <- relative_mean(fit$common, years)
contrib <- sapply(seq_len(K_LEVELS), function(j)
  relative_mean(outer(fit$factors[, j], fit$loadings[, j]), years))
colnames(contrib) <- colnames(fit$factors)
save_pdf(plot_common_signal(xbar, dates), "fig4a_common_component.pdf", 8.2, 5)
save_pdf(plot_factor_contributions(contrib, dates), "fig4b_factor_contributions.pdf", 8.5, 5.5)

## Figures 5 and 6: loadings and station-specific warming
stations <- data.frame(meta[, c("id", "name", "country_code", "region", "lat", "lon")],
                       loading_1 = fit$loadings[, 1], loading_2 = fit$loadings[, 2],
                       warming = warming(fit$common, years), row.names = NULL)
regions <- aggregate(warming ~ region, stations, mean)
regions$n <- as.vector(table(stations$region)[regions$region])
write.csv(stations, file.path(TAB_DIR, "station_warming_levels.csv"), row.names = FALSE)
write.csv(regions, file.path(TAB_DIR, "regional_warming_levels.csv"), row.names = FALSE)

save_pdf(plot_loading_map(stations, "loading_1"), "fig5a_loadings_factor1.pdf", 6.2, 6.2)
save_pdf(plot_loading_map(stations, "loading_2"), "fig5b_loadings_factor2.pdf", 6.2, 6.2)
save_pdf(plot_warming_map(stations, "warming"), "fig6a_station_warming.pdf", 9, 6.2)
save_pdf(plot_regions(regions), "fig6b_regional_warming.pdf", 8, 5.2)

## Numbers reported in Section 5.2
decade <- 10 * (years %/% 10)
per_decade <- function(from, to) {
  i <- years >= from & years <= to
  120 * unname(coef(lm(xbar[i] ~ seq_len(sum(i))))[2])
}
cat("\nStable intervals of c (levels):\n"); print(abc$stable, row.names = FALSE)
cat(sprintf("Eigenvalue ratio 2nd/3rd: %.2f\n", fit$eigenvalues[2] / fit$eigenvalues[3]))
cat("\nDecade means of the standardized factors:\n")
print(round(aggregate(Z, list(decade = decade), mean), 2), row.names = FALSE)
cat("\nDecade means of the common component relative to the baseline:\n")
print(round(tapply(xbar, decade, mean), 2))
cat(sprintf("Trend per decade: %.2f (1960-1985), %.2f (1986-2025)\n",
            per_decade(1960, 1985), per_decade(1986, 2025)))
cat(sprintf("Average warming: %.2f; factor contributions: %s\n", mean(stations$warming),
            paste(sprintf("%.2f", colMeans(sapply(seq_len(K_LEVELS), function(j)
              warming(outer(fit$factors[, j], fit$loadings[, j]), years)))), collapse = ", ")))
cat(sprintf("Positive loadings on factor 1: %d of %d\n", sum(stations$loading_1 > 0), nrow(stations)))
cat("Station warming (min, 10%, median, 90%, max):\n")
print(round(quantile(stations$warming, c(0, 0.1, 0.5, 0.9, 1)), 2))
print(regions[order(-regions$warming), ], row.names = FALSE)

saveRDS(list(dates = dates, years = years, Y_low = Y_low, fit = fit, stations = stations),
        LEVELS_RESULTS)
