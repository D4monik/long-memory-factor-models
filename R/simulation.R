# Monte Carlo study (Section 4).

# Type-II fractional difference (1 - L)^d of each column of y, eq. (2); d may differ
# across columns, and a negative d gives fractional integration.
frac_diff <- function(y, d) {
  n <- nrow(y)
  vapply(seq_len(ncol(y)), function(j) {
    coef <- c(1, cumprod(((1:n) - 1 - d[j]) / (1:n)))
    stats::filter(c(rep(0, n), y[, j]), filter = coef, sides = 1)[-(1:n)]
  }, numeric(n))
}

# One replication of the data-generating process (28)-(29): the number of factors
# selected by IC_p2 in a panel with N series, Tobs periods and idiosyncratic
# variance sigma2. With diminishing = TRUE, the innovations of the nonstationary
# factors are scaled as in Assumption 6; with differences = TRUE, the criterion is
# applied to the first differences of the panel.
simulate_selection <- function(rep, d0, N, Tobs, sigma2, diminishing, differences,
                               lambda_range, s_max, seed) {
  set.seed(seed + rep)
  s0 <- length(d0)
  scale <- if (diminishing) 1 / Tobs^pmax(d0 - 1/2, 0) else rep(1, s0)
  eta <- sweep(matrix(rnorm(Tobs * s0), ncol = s0), 2, scale, "*")
  x <- frac_diff(eta, -d0)
  Lambda <- matrix(runif(N * s0, lambda_range[1], lambda_range[2]), ncol = N)
  y <- x %*% Lambda + matrix(rnorm(N * Tobs), nrow = Tobs) * sqrt(sigma2)
  if (differences) y <- diff(y)
  ev <- eigen(crossprod(y) / nrow(y), symmetric = TRUE, only.values = TRUE)$values
  select_factors(ev, N, nrow(y), s_max)
}

# Bisection for the idiosyncratic variance at which the selection frequency
# prob(sigma2), which decreases in sigma2, equals the target up to tol.
calibrate_threshold <- function(prob, lower, upper, target = 0.5, tol = 0.01, max_iter = 60L) {
  for (i in seq_len(max_iter)) {
    sigma2 <- (lower + upper) / 2
    p <- prob(sigma2)
    if (abs(p - target) <= tol) return(sigma2)
    if (p < target) upper <- sigma2 else lower <- sigma2
  }
  warning("bisection did not reach the tolerance; returning the last midpoint")
  sigma2
}
