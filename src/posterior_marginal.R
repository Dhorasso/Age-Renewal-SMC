# ============================================================
# posterior_marginal.R  (formerly PosteriorMarginal.R)
#
# Draws a posterior sample of latent trajectories (X) and predicted
# observations (X_obs), marginalising over both the theta-particles
# and their PF particle clouds. Always stores X and X_obs.
# ============================================================

library(parallel)

source("src/utils.R")
source("src/particle_filter.R")

#' Sample latent-state trajectories from the joint posterior
#'
#' @param theta matrix [N_theta, d] of posterior theta-particles
#'   (e.g. `SMC2()$theta[, T, ]`)
#' @param SSM state-space model list
#' @param opts options list. Relevant fields:
#'   \describe{
#'     \item{N}{PF particles per theta-draw (default 4000)}
#'     \item{Ns}{trajectories drawn per theta-draw (default 200)}
#'     \item{Nsub}{number of theta-particles to sub-sample (default 100)}
#'     \item{T, forecastingHorizon, state_dim, obs_dim}{as in `BootstrapPF()`}
#'   }
#' @return list with `X` [Nsub*Ns, T+h, state_dim] and
#'   `X_obs` [Nsub*Ns, T+h, obs_dim]
PosteriorMarginal <- function(theta, SSM, opts) {

  opts$N <- opts$N %||% 4000
  N      <- opts$N
  Ns     <- opts$Ns   %||% 200
  Nsub   <- opts$Nsub %||% 100
  opts$PredictedObservation <- TRUE

  T         <- opts$T
  h         <- opts$forecastingHorizon %||% 0
  state_dim <- opts$state_dim
  obs_dim   <- opts$obs_dim

  thetaSamples <- theta#[sample(nrow(theta), Nsub), , drop = FALSE]
  Ntheta_eff   <- nrow(thetaSamples)

  # pre-sample particle indices once, for reproducible sub-sampling
  selected_idx <- matrix(
    sample(1:N, Ns * Ntheta_eff, replace = TRUE),
    nrow = Ntheta_eff, ncol = Ns
  )

  X_all     <- matrix(NA_real_, nrow = Ntheta_eff * Ns, ncol = (T + h) * state_dim)
  X_obs_all <- matrix(NA_real_, nrow = Ntheta_eff * Ns, ncol = (T + h) * obs_dim)

  for (i in seq_len(Ntheta_eff)) {
    theta_i <- thetaSamples[i, , drop = TRUE]
    bf_res  <- BootstrapPF(SSM, theta_i, opts)

    idx   <- ((i - 1) * Ns + 1):(i * Ns)
    s_idx <- selected_idx[i, ]

    X_all[idx, ]     <- matrix(bf_res$X[s_idx, , ],     nrow = Ns)
    X_obs_all[idx, ] <- matrix(bf_res$X_obs[s_idx, , ], nrow = Ns)
  }

  list(
    X     = array(X_all,     dim = c(Ntheta_eff * Ns, T + h, state_dim)),
    X_obs = array(X_obs_all, dim = c(Ntheta_eff * Ns, T + h, obs_dim))
  )
}
