# Factor-number selection, used in the Monte Carlo study and in the application.

# Number of factors selected by the logarithmic criterion
#   IC_c(s) = log V(s) + c s (n + T)/(nT) log(min(n, T)),   s = 0, ..., s_max,
# for every penalty multiplier in c_grid, given the eigenvalues of the sample
# covariance matrix of an (T x n) panel; V(s) is the average of the discarded
# eigenvalues, eq. (11). For c = 1 this is IC_p2 of Bai and Ng (2002).
select_factors <- function(eigenvalues, n, Tobs, s_max, c_grid = 1) {
  V <- rev(cumsum(rev(eigenvalues)))[seq_len(s_max + 1L)] / n   # V(0), ..., V(s_max)
  penalty <- (n + Tobs) / (n * Tobs) * log(min(n, Tobs))
  IC <- log(V) + outer(0:s_max, c_grid * penalty)
  max.col(-t(IC), ties.method = "first") - 1L
}

# Factor-number selection of Alessi, Barigozzi and Capasso (2010): selected number
# of factors in the full panel, and its standard deviation across the full panel
# and `reps` random subpanels of each size N - floor(N/2), ..., N - 1.
abc_selection <- function(X, s_max, c_grid, reps, seed) {
  X <- as.matrix(X)
  Tobs <- nrow(X)
  N <- ncol(X)
  Sigma <- crossprod(X) / Tobs
  select_subpanel <- function(idx) {
    ev <- eigen(Sigma[idx, idx, drop = FALSE], symmetric = TRUE, only.values = TRUE)$values
    select_factors(ev, length(idx), Tobs, s_max, c_grid)
  }

  full <- select_subpanel(seq_len(N))
  sizes <- rep(seq.int(N - floor(N / 2), N - 1L), each = reps)
  s1 <- as.numeric(full); s2 <- s1^2; lo <- hi <- full
  set.seed(seed)
  for (n in sizes) {
    s <- select_subpanel(sample.int(N, n))
    s1 <- s1 + s; s2 <- s2 + s^2
    lo <- pmin(lo, s); hi <- pmax(hi, s)
  }
  M <- length(sizes) + 1L
  sd_c <- sqrt(pmax(s2 - s1^2 / M, 0) / (M - 1))
  sd_c[lo == hi] <- 0
  path <- data.frame(c = c_grid, s_hat = full, sd = sd_c)

  # intervals of c on which all panels agree and s_hat < s_max
  stable <- ifelse(lo == hi & full < s_max, full, NA)
  runs <- rle(stable)
  end <- cumsum(runs$lengths)
  start <- end - runs$lengths + 1L
  keep <- which(!is.na(runs$values) & runs$lengths > 1L)
  list(path = path,
       stable = data.frame(c_from = c_grid[start[keep]], c_to = c_grid[end[keep]],
                           s_hat = runs$values[keep]))
}
