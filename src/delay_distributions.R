# ============================================================
# delay_distributions.R  
#
# Discretized delay-distribution PMFs used for reporting-delay /
# generation-time kernels.
# ============================================================


#' Discretized gamma or lognormal PMF, conditioned on delay >= 1
#'
#' Computes P(delay = x) for x = 1..max via differenced CDFs at
#' half-integer boundaries, then drops x = 0 and renormalizes. Used
#' for generation-time / reporting-delay kernels where a delay of
#' exactly zero days is excluded by construction.
#'
#' @param max maximum delay (in days) to retain in the PMF
#' @param mean mean of the underlying continuous distribution
#' @param sd standard deviation of the underlying continuous distribution
#' @param dist one of "gamma" (default) or "lognormal"
#' @return numeric vector of length `max`, the normalized PMF over 1:max
discretize_dist <- function(max, mean, sd, dist = c("gamma", "lognormal")) {
  dist <- match.arg(dist)
  x <- 0:max

  if (dist == "gamma") {
    shape <- (mean / sd)^2
    scale <- sd^2 / mean
    pmf <- pgamma(x + 0.5, shape = shape, scale = scale) -
      pgamma(pmax(x - 0.5, 0), shape = shape, scale = scale)

  } else {
    sigma2 <- log(1 + (sd^2 / mean^2))
    sigma  <- sqrt(sigma2)
    mu     <- log(mean) - sigma2 / 2
    pmf <- plnorm(x + 0.5, meanlog = mu, sdlog = sigma) -
      plnorm(pmax(x - 0.5, 0), meanlog = mu, sdlog = sigma)
  }

  pmf <- pmf[-1]         # drop x = 0
  pmf / sum(pmf)
}
