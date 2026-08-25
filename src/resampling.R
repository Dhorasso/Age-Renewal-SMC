# ============================================================
# resampling.R
#
# Resampling schemes for (bootstrap) particle filters:
#   - multinomial
#   - systematic   (default; low variance, O(N))
#   - stratified
# ============================================================

#' Resample particle indices given weights
#'
#' @param w numeric vector of (normalized or unnormalized) particle weights
#' @param method one of "systematic" (default), "multinomial", "stratified"
#' @return integer vector of length `length(w)` with resampled indices
#'
#' @details Non-finite or negative weights are treated as zero; if all
#' weights collapse to zero (or are non-finite) the function falls back
#' to uniform weights rather than failing, which keeps particle filters
#' robust to occasional numerical degeneracies.
Resample <- function(w, method = c("systematic", "multinomial", "stratified")) {

  method <- match.arg(method)
  N <- length(w)

  # --- clean weights ---
  w[!is.finite(w) | w < 0] <- 0
  total <- sum(w)

  if (total <= 0 || !is.finite(total)) {
    w <- rep(1 / N, N)
  } else {
    w <- w / total
  }
  w <- w / sum(w)  # guarantee exact normalization

  if (method == "multinomial") {
    return(sample.int(N, size = N, replace = TRUE, prob = w))
  }

  cs    <- cumsum(w)
  cs[N] <- 1.0

  positions <- switch(
    method,
    systematic = (runif(1) + (0:(N - 1))) / N,
    stratified = ((0:(N - 1)) + runif(N)) / N
  )

  .cdf_lookup(positions, cs, N)
}

#' Vectorised inverse-CDF lookup shared by systematic/stratified resampling
#'
#' @param positions sorted vector of N points in [0, 1)
#' @param cs cumulative weight vector (cs[N] == 1)
#' @param N number of particles
#' @return integer vector of resampled indices
.cdf_lookup <- function(positions, cs, N) {
  indexes <- integer(N)
  i <- 1L
  j <- 1L
  while (i <= N) {
    if (positions[i] <= cs[j]) {
      indexes[i] <- j
      i <- i + 1L
    } else {
      j <- j + 1L
      if (j > N) j <- N  # clamp instead of erroring on rounding edge cases
    }
  }
  indexes
}
