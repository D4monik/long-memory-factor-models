# Section 5.4: the same analysis after first differencing the low-frequency panel.
# Output: Figures 8 and 9, Table 2

source("R/config.R")
source("R/estimation.R")
source("R/factor_selection.R")
source("R/plots.R")

res <- readRDS(LEVELS_RESULTS)
years <- res$years
dY <- diff(res$Y_low)

## Figure 8(a): number of factors after differencing
abc <- abc_selection(dY, ABC_S_MAX, ABC_C, ABC_REPS, ABC_SEED)
write.csv(abc$stable, file.path(TAB_DIR, "abc_stable_intervals_differences.csv"), row.names = FALSE)
save_pdf(plot_abc(abc$path, ABC_S_MAX), "fig8a_abc_differences.pdf", 7.5, 5)

## One-factor model for the differenced panel; common component cumulated from zero
fit <- fit_pca(dY, K_DIFFERENCES, regression = FALSE)
C_lev <- res$fit$common
C_dif <- rbind(0, apply(fit$common, 2, cumsum))

## Figure 8(b): sample-average common component, levels vs. differences
save_pdf(plot_levels_vs_differences(relative_mean(C_lev, years), relative_mean(C_dif, years), res$dates),
         "fig8b_levels_vs_differences.pdf", 6, 4)

## Table 2: regional mean warming
S <- data.frame(res$stations[, c("id", "region", "lat", "lon")],
                observed = warming(res$Y_low, years),
                levels = warming(C_lev, years),
                differences = warming(C_dif, years), row.names = NULL)
reg <- aggregate(cbind(observed, levels, differences) ~ region, S, mean)
reg$n <- as.vector(table(S$region)[reg$region])
reg <- reg[order(-reg$observed), ]
total <- data.frame(region = "All stations", observed = mean(S$observed), levels = mean(S$levels),
                  differences = mean(S$differences), n = nrow(S))
write.csv(rbind(reg, total), file.path(TAB_DIR, "table2_regional_warming.csv"), row.names = FALSE)
row_tex <- function(d) sprintf("%s ($%d$) & $%.2f$ & $%.2f$ & $%.2f$ \\\\",
                               d$region, d$n, d$observed, d$levels, d$differences)
writeLines(c("\\begin{tabular}{@{}lccc@{}}", "\\toprule",
             "& Observed & Levels, two factors & Differences, one factor \\\\", "\\midrule",
             row_tex(reg), "\\midrule", row_tex(total), "\\bottomrule", "\\end{tabular}"),
           file.path(TAB_DIR, "table2_regional_warming.tex"))

## Figure 9: station warming after differencing minus warming in levels
D <- data.frame(S[, c("id", "lat", "lon")], diff = S$differences - S$levels)
save_pdf(plot_difference_map(D), "fig9_difference_map.pdf", 9, 6.2)

## Numbers reported in Section 5.4
ev <- fit$eigenvalues
cat("\nStable intervals of c (differences):\n"); print(abc$stable, row.names = FALSE)
cat(sprintf("Eigenvalue ratio 2nd/3rd after differencing: %.2f\n", ev[2] / ev[3]))
for (v in c("levels", "differences"))
  cat(sprintf("%-11s corr with observed %.2f, RMSE %.2f\n", v,
              cor(S[[v]], S$observed), sqrt(mean((S[[v]] - S$observed)^2))))
cat(sprintf("Differencing lowers warming at %d stations (by up to %.2f) and raises it at %d (by up to %.2f)\n",
            sum(D$diff < 0), -min(D$diff), sum(D$diff > 0), max(D$diff)))
print(rbind(reg, total), row.names = FALSE, digits = 3)
