# ============================================================
# pmmh.R
#
# Particle Marginal Metropolis-Hastings (PMMH), using the bootstrap
# particle filter for an unbiased likelihood estimate and
# BayesianTools::runMCMC() (DEzs sampler) to explore theta-space.
#
# Reuses the repo's existing prior machinery (priors.R) rather than
# redefining prior density/sampling inline, and the existing bootstrap
# filter (particle_filter.R) for the likelihood estimate - this file
# only adds the BayesianTools wiring on top.
# ============================================================

source("src/utils.R")
source("src/particle_filter.R")

if (!requireNamespace("BayesianTools", quietly = TRUE)) {
  install.packages("BayesianTools")
}
library(BayesianTools)
library(compiler)

# ------------------------------------------------------------------
# PMMH
# ------------------------------------------------------------------
#' Particle Marginal Metropolis-Hastings for a state-space model
#'
#' @param SSM state-space model list (`InitState`, `StateProcess`, `ObsProcess`)
#' @param opts options list. In addition to the fields used by
#'   `BootstrapPF()`, relevant fields are:
#'   \describe{
#'     \item{paramPriors}{list of prior specs, see priors.R}
#'     \item{paramNames}{character/expression vector, parameter labels}
#'     \item{lower_bounds, upper_bounds}{box constraints on theta}
#'     \item{iterations}{total MCMC iterations across all chains (default 10000)}
#'     \item{nChains}{number of DEzs chains (default 2)}
#'     \item{burnin}{iterations discarded as burn-in (default 30% of `iterations`)}
#'     \item{message}{print sampler progress? (default TRUE)}
#'   }
#' @return a `bayesianOutput` object from `BayesianTools::runMCMC()`
PMMH_inner <- function(SSM, opts) {
  
  # --- unbiased log-likelihood estimate via the bootstrap PF ---
  # Not cached: the PF is stochastic (different value on repeat calls
  # at the same theta), MCMC explores a continuous space (exact repeat
  # evaluations are rare), and cache bookkeeping would cost more than
  # it saves here.
  likelihood <- function(theta) {
    if (!is.null(opts$lower_bounds) && !is.null(opts$upper_bounds)) {
      if (any(theta < opts$lower_bounds) || any(theta > opts$upper_bounds)) {
        return(-Inf)
      }
    }
    BootstrapPF(SSM, theta, opts)$Loglike
  }
  
  # --- prior density / sampler (delegates to priors.R) ---
  density_fn <- function(theta) logPriorDensity(theta, opts)
  sampler_fn <- function(n = 1) priorSampler(n, opts)
  
  prior <- createPrior(
    density = density_fn,
    sampler = sampler_fn,
    lower   = opts$lower_bounds %||% NULL,
    upper   = opts$upper_bounds %||% NULL
  )
  
  bayesianSetup <- createBayesianSetup(
    likelihood = likelihood,
    prior      = prior,
    names      = opts$paramNames
  )
  
  iterations <- opts$iterations %||% 10000
  settings <- list(
    iterations = iterations,
    nrChains   = opts$nChains %||% 2,
    burnin     = opts$burnin  %||% floor(0.3 * iterations),
    thin       = 3,
    message    = opts$message %||% TRUE
  )
  
  runMCMC(bayesianSetup = bayesianSetup, sampler = "DEzs", settings = settings)
}

PMMH <- cmpfun(PMMH_inner)
