# Section 4: largest idiosyncratic standard deviation at which IC_p2 selects the true
# number of factors in at least half of the replications (Table 1).
# The search can be interrupted: finished cells are stored in output/ and skipped
# when the script is run again.
# Output: output/tables/table1_thresholds.{tex,csv}

source("R/config.R")
source("R/factor_selection.R")
source("R/simulation.R")

designs <- data.frame(
  name = c("Standard signals, levels", "Standard signals, first differences",
           "Diminishing signals, levels", "Diminishing signals, first differences"),
  diminishing = c(FALSE, FALSE, TRUE, TRUE),
  differences = c(FALSE, TRUE, FALSE, TRUE)
)

# calibrated variances, T x N x design
thresholds <- if (file.exists(SIM_RESULTS)) readRDS(SIM_RESULTS) else
  array(NA_real_, dim = c(length(SIM_T), length(SIM_N), nrow(designs)),
        dimnames = list(paste0("T = ", SIM_T), paste0("N = ", SIM_N), designs$name))

cl <- parallel::makeCluster(SIM_CORES)
parallel::clusterExport(cl, c("simulate_selection", "frac_diff", "select_factors"))

for (k in seq_len(nrow(designs))) {
  for (i in seq_along(SIM_T)) {
    for (j in seq_along(SIM_N)) {
      if (!is.na(thresholds[i, j, k])) next

      # share of replications in which the true number of factors is selected
      selection_frequency <- function(sigma2) {
        s_hat <- parallel::parSapply(cl, seq_len(SIM_REPS), simulate_selection,
                                     d0 = SIM_D0, N = SIM_N[j], Tobs = SIM_T[i], sigma2 = sigma2,
                                     diminishing = designs$diminishing[k],
                                     differences = designs$differences[k],
                                     lambda_range = SIM_LAMBDA, s_max = SIM_S_MAX, seed = SIM_SEED)
        mean(s_hat == length(SIM_D0))
      }
      thresholds[i, j, k] <- calibrate_threshold(selection_frequency, SIM_SIGMA2[1], SIM_SIGMA2[2],
                                                 target = 0.5, tol = SIM_TOL)
      saveRDS(thresholds, SIM_RESULTS)
      message(sprintf("%s | T = %d, N = %d | sigma_u = %.3f",
                      designs$name[k], SIM_T[i], SIM_N[j], sqrt(thresholds[i, j, k])))
    }
  }
}
parallel::stopCluster(cl)

## Table 1: standard deviations, one panel per design
sd_u <- sqrt(thresholds)
table_rows <- function(M) {
  sprintf("$T = %d$ & %s \\\\", SIM_T,
          apply(M, 1, function(r) paste(sprintf("%.3f", r), collapse = " & ")))
}
tex <- unlist(lapply(seq_len(nrow(designs)), function(k) c(
  paste0("% ", designs$name[k]),
  sprintf("\\begin{tabular}{l%s}", strrep("c", length(SIM_N))),
  "\\toprule",
  paste0(" & ", paste(sprintf("$N = %d$", SIM_N), collapse = " & "), " \\\\"),
  "\\midrule",
  table_rows(sd_u[, , k]),
  "\\bottomrule", "\\end{tabular}", ""
)))
writeLines(tex, file.path(TAB_DIR, "table1_thresholds.tex"))

long <- as.data.frame.table(sd_u, responseName = "sigma_u")
names(long)[1:3] <- c("T", "N", "design")
write.csv(long, file.path(TAB_DIR, "table1_thresholds.csv"), row.names = FALSE)
