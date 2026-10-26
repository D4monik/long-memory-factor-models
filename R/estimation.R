# Estimation functions for the climate application.

# Low-frequency component, eqs. (30)-(32): station mean plus the projection on the
# first J_L = floor(2T / cutoff) - 1 cosine functions (periods longer than `cutoff`).
lowpass <- function(Y, cutoff) {
  Y <- as.matrix(Y)
  n <- nrow(Y)
  q <- floor(2 * n / cutoff) - 1
  phi <- sqrt(2) * cos(pi * outer(seq_len(n) - 0.5, seq_len(q)) / n)
  Y_low <- sweep(phi %*% (crossprod(phi, Y) / n), 2, colMeans(Y), "+")
  dimnames(Y_low) <- dimnames(Y)
  Y_low
}

# Principal components of the uncentered panel, eqs. (6)-(8). Factor signs are
# chosen such that each factor has a positive linear trend. With regression = TRUE,
# the loadings are obtained by regressing each series on a constant and the factors;
# otherwise they are the PCA loadings.
fit_pca <- function(Y, k, regression = TRUE) {
  Y <- as.matrix(Y)
  pc <- prcomp(Y, center = FALSE, scale. = FALSE)
  Fx <- pc$x[, seq_len(k), drop = FALSE]
  sgn <- apply(Fx, 2, function(f) if (coef(lm(f ~ seq_along(f)))[2] >= 0) 1 else -1)
  Fx <- sweep(Fx, 2, sgn, "*")
  if (regression) {
    B <- lm.fit(cbind(1, Fx), Y)$coefficients
    intercept <- B[1, ]
    Lambda <- t(B[-1, , drop = FALSE])
  } else {
    intercept <- rep(0, ncol(Y))
    Lambda <- sweep(pc$rotation[, seq_len(k), drop = FALSE], 2, sgn, "*")
  }
  colnames(Fx) <- colnames(Lambda) <- paste("Factor", seq_len(k))
  rownames(Lambda) <- colnames(Y)
  list(factors = Fx, loadings = Lambda, intercept = intercept,
       common = sweep(Fx %*% t(Lambda), 2, intercept, "+"),
       eigenvalues = pc$sdev^2)
}

# Warming per column: mean over the recent period minus mean over the baseline.
warming <- function(C, years) {
  b <- years >= BASELINE[1] & years <= BASELINE[2]
  r <- years >= RECENT[1] & years <= RECENT[2]
  colMeans(C[r, , drop = FALSE]) - colMeans(C[b, , drop = FALSE])
}

# Cross-sectional mean relative to its baseline mean.
relative_mean <- function(C, years) {
  m <- rowMeans(C)
  m - mean(m[years >= BASELINE[1] & years <= BASELINE[2]])
}

# Factors relative to their baseline means, in standard deviations.
standardize_factors <- function(Fx, years) {
  b <- years >= BASELINE[1] & years <= BASELINE[2]
  sweep(sweep(Fx, 2, colMeans(Fx[b, , drop = FALSE])), 2, apply(Fx, 2, sd), "/")
}

# Regions of footnote 5, defined by station coordinates.
region_of <- function(lat, lon) {
  ifelse(lat >= 55, "Northern Europe",
  ifelse(lat < 45 & lon < 5, "Southwestern Europe",
  ifelse(lat < 45, "Southeastern Europe",
  ifelse(lon >= 20, "Eastern Europe", "Central/Western Europe"))))
}
