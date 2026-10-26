# Low-frequency likelihood for the memory parameter d (Mueller and Watson, 2008),
# Section 5.3. The cosine transforms X_j = T^{-1} sum_t phi_j((t - 1/2)/T) x_t,
# j = 1, ..., q, are approximately N(0, s^2 Omega(d)) with
# Omega(d) = T^{-2} Phi' A_d A_d' Phi for the type-II I(d) process, where A_d is the
# lower-triangular Toeplitz matrix of the coefficients of (1 - L)^{-d}.
# Inference uses the scale-invariant likelihood of X / ||X||, eq. (34).

lf_basis <- function(Tobs, q) {
  sqrt(2) * cos(pi * outer(seq_len(Tobs) - 0.5, seq_len(q)) / Tobs)
}

# A_d' Phi, computed by FFT; a = first column of A_d.
lf_AtPhi <- function(a, Phi) {
  n <- nrow(Phi); L <- 2L * n
  Pr <- Phi[n:1, , drop = FALSE]
  fa <- fft(c(a, rep(0, L - n)))
  R <- Re(mvfft(mvfft(rbind(Pr, matrix(0, L - n, ncol(Phi)))) * fa, inverse = TRUE))[1:n, , drop = FALSE] / L
  R[n:1, , drop = FALSE]
}

# Omega(d) on a grid of d (normalized to trace q; the likelihood is scale invariant).
lf_prepare <- function(Tobs, q, d_grid) {
  Phi <- lf_basis(Tobs, q)
  Omega <- lapply(d_grid, function(d) {
    a <- c(1, cumprod(((1:(Tobs - 1)) - 1 + d) / (1:(Tobs - 1))))
    W <- crossprod(lf_AtPhi(a, Phi))
    W * q / sum(diag(W))
  })
  list(Phi = Phi, d_grid = d_grid, Omega = Omega)
}

# log L(d) = -1/2 log|Omega| - q/2 log(X' Omega^{-1} X)
lf_loglik <- function(X, W) {
  R <- tryCatch(chol(W), error = function(e) NULL)
  if (is.null(R)) return(-Inf)
  z <- backsolve(R, X, transpose = TRUE)
  -sum(log(diag(R))) - length(X) / 2 * log(sum(z^2))
}

# Maximum likelihood estimate of d, likelihood-ratio statistics on the grid, and
# the confidence set obtained by inverting the likelihood-ratio test.
lf_estimate <- function(x, prep, level = 0.95) {
  X <- drop(crossprod(prep$Phi, x)) / length(x)
  ll <- vapply(prep$Omega, function(W) lf_loglik(X, W), numeric(1))
  lr <- 2 * (max(ll) - ll)
  list(d = prep$d_grid[which.max(ll)],
       ci = range(prep$d_grid[lr <= qchisq(level, 1)]),
       lr = lr)
}
