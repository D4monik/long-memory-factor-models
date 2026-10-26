# Section 5.3: memory of the estimated factors (low-frequency likelihood of
# Mueller and Watson, 2008, based on the first q = 13 cosine transforms).
# Output: Figure 7, table of estimates and 95% confidence intervals

source("R/config.R")
source("R/lowfreq_memory.R")
source("R/plots.R")

res <- readRDS(LEVELS_RESULTS)
factors <- res$fit$factors
prep <- lf_prepare(nrow(factors), MW_Q, MW_D_GRID)
est <- lapply(seq_len(ncol(factors)), function(j) lf_estimate(factors[, j], prep))

memory <- data.frame(factor = colnames(factors),
                     d_hat = sapply(est, `[[`, "d"),
                     ci_lower = sapply(est, function(e) e$ci[1]),
                     ci_upper = sapply(est, function(e) e$ci[2]))
write.csv(memory, file.path(TAB_DIR, "memory_estimates.csv"), row.names = FALSE)
print(memory, row.names = FALSE)

prof <- data.frame(factor = rep(colnames(factors), each = length(MW_D_GRID)),
                   d = MW_D_GRID, LR = unlist(lapply(est, `[[`, "lr")))
save_pdf(plot_memory(prof), "fig7_memory.pdf", 8, 3.8)
