# ============================================================
# mcmc_kernels.R
#
# Hastings-within-Gibbs move-step kernels used in the SMC2()
# resample-move step, plus supporting helpers (kNN log-likelihood
# surrogate, adaptive block proposal covariance, adaptive Nx).
# Extracted from the original SMC2.R.
#
# HWG_kernel     - exact block-wise Metropolis-within-Gibbs kernel.
# DA_HWG_kernel  - delayed-acceptance version: a cheap kNN surrogate
#                  screens proposals before the (expensive) exact
#                  particle-filter likelihood is evaluated.
#
# Both call PF_SMC2() (defined in smc2.R) to evaluate the exact
# incremental likelihood; smc2.R in turn calls these kernels, so the
# two files must both be sourced before SMC2() is run (the mutual
# call is resolved at call time, not at source time, so load order
# between these two files does not matter).
# ============================================================

library(MASS)  # mvrnorm
library(FNN)   # get.knnx

source("src/utils.R")


# ------------------------------------------------------------------
# kNN surrogate log-likelihood
# ------------------------------------------------------------------
#' Cheap surrogate for the exact log-likelihood via inverse-distance
#' weighted k-nearest-neighbour interpolation over previously
#' evaluated theta-particles.
#'
#' @param theta_new parameter vector (or 1-row matrix) to query
#' @param theta_unique matrix of previously evaluated unique parameters
#' @param logLik_unique numeric vector, log-likelihoods at `theta_unique`
#' @param k number of neighbours (default 3)
#' @return scalar surrogate log-likelihood
knn_surrogate <- function(theta_new, theta_unique, logLik_unique, k = 3L) {
  if (is.null(theta_unique) || nrow(theta_unique) == 0L) return(-Inf)
  k <- min(k, nrow(theta_unique))
  if (is.vector(theta_new)) theta_new <- matrix(theta_new, nrow = 1L)

  res       <- FNN::get.knnx(theta_unique, theta_new, k = k)
  indices   <- res$nn.index[1L, ]
  distances <- res$nn.dist[1L, ]

  exact_hit <- which(distances < 1e-10)
  if (length(exact_hit) > 0L) return(logLik_unique[indices[exact_hit[1L]]])

  distances <- pmax(distances, 1e-12)
  weights   <- 1 / distances
  weights   <- weights / sum(weights)
  sum(weights * logLik_unique[indices])
}

# ------------------------------------------------------------------
# Adaptive (Haario-type) block proposal covariance
# ------------------------------------------------------------------
#' @param theta_cov full parameter covariance estimate, [d, d]
#' @param block_idx integer indices of the block being updated
#' @param d total number of parameters
#' @return proposal covariance matrix for this block, [b, b]
block_proposal_cov <- function(theta_cov, block_idx, d) {
  b       <- length(block_idx)
  Sigma_b <- theta_cov[block_idx, block_idx, drop = FALSE]
  if (b == 1L) Sigma_b <- matrix(Sigma_b, 1L, 1L)
  (2.38^2 / d) * Sigma_b + 1e-6 * diag(b)
}

#' Resolve the parameter blocking scheme for the move step
#'
#' @param opts options list; `opts$param_blocks` overrides automatic blocking
#' @param d total number of parameters
#' @param n_blocks desired number of (roughly equal-sized) blocks
#' @return list of integer index vectors, one per block
.resolve_blocks <- function(opts, d, n_blocks) {
  if (!is.null(opts$param_blocks)) return(opts$param_blocks)
  bs <- ceiling(d / (n_blocks %||% d))
  split(seq_len(d), ceiling(seq_len(d) / bs))
}

# ------------------------------------------------------------------
# HWG KERNEL  (exact, block-wise Metropolis-within-Gibbs)
# ------------------------------------------------------------------
#' @param theta current parameter vector
#' @param logLik current log-likelihood
#' @param logPrior current log-prior
#' @param X current particle history for this theta-particle
#' @param SSM state-space model list
#' @param opts options list (must contain `N_move`)
#' @param theta_cov current parameter covariance estimate
#' @param t current real-time index
#' @param n_blocks number of blocks (ignored if `opts$param_blocks` is set)
#' @return list with updated `theta`, `logLik`, `logPrior`, `X`,
#'   `n_accept`, `n_proposals`
HWG_kernel <- function(theta, logLik, logPrior, X,
                        SSM, opts, theta_cov, t, n_blocks = NULL) {
  d      <- length(theta)
  blocks <- .resolve_blocks(opts, d, n_blocks)

  n_accept    <- 0L
  n_proposals <- 0L

  for (m in seq_len(opts$N_move)) {
    for (block_idx in blocks) {

      theta_p            <- theta
      Sigma_b             <- block_proposal_cov(theta_cov, block_idx, d)
      theta_p[block_idx]  <- as.numeric(MASS::mvrnorm(1L, theta[block_idx], Sigma_b))

      lp_p        <- logPriorDensity(theta_p, opts)
      n_proposals <- n_proposals + 1L
      if (!is.finite(lp_p)) next

      pf        <- PF_SMC2(SSM, theta_p, opts, t)
      ll_p      <- pf$Loglike
      log_alpha <- lp_p + ll_p - logPrior - logLik

      if (is.finite(log_alpha) && log(runif(1L)) < log_alpha) {
        theta    <- theta_p
        logLik   <- ll_p
        logPrior <- lp_p
        X        <- pf$X
        n_accept <- n_accept + 1L
      }
    }
  }

  list(theta = theta, logLik = logLik, logPrior = logPrior, X = X,
       n_accept = n_accept, n_proposals = n_proposals)
}

# ------------------------------------------------------------------
# DA-HWG KERNEL  (delayed-acceptance, block-wise)
# ------------------------------------------------------------------
#' @inheritParams HWG_kernel
#' @param theta_unique matrix of unique previously-evaluated parameters
#' @param logLik_unique numeric vector, log-likelihoods at `theta_unique`
#' @param k_neighbors number of neighbours for the kNN surrogate
#' @return list with updated `theta`, `logLik`, `logPrior`, `X`,
#'   `n_accept`, `n_proposals`, `ll_sur_current`
DA_HWG_kernel <- function(theta, logLik, logPrior, X,
                           SSM, opts, theta_cov, t,
                           theta_unique, logLik_unique,
                           k_neighbors = 3L, n_blocks = NULL) {
  d <- length(theta)

  ll_sur_current <- knn_surrogate(theta, theta_unique, logLik_unique, k_neighbors)
  blocks         <- .resolve_blocks(opts, d, n_blocks)

  n_accept    <- 0L
  n_proposals <- 0L

  for (m in seq_len(opts$N_move)) {
    for (block_idx in blocks) {

      theta_p            <- theta
      Sigma_b            <- block_proposal_cov(theta_cov, block_idx, d)
      theta_p[block_idx] <- as.numeric(MASS::mvrnorm(1L, theta[block_idx], Sigma_b))

      lp_p        <- logPriorDensity(theta_p, opts)
      n_proposals <- n_proposals + 1L
      if (!is.finite(lp_p)) next

      # Stage 1: cheap surrogate-vs-surrogate screen
      ll_sur_p   <- knn_surrogate(theta_p, theta_unique, logLik_unique, k_neighbors)
      log_alpha1 <- (lp_p - logPrior) + (ll_sur_p - ll_sur_current)

      if (is.finite(log_alpha1) && log(runif(1L)) < log_alpha1) {

        # Stage 2: exact correction against the true particle-filter likelihood
        pf         <- PF_SMC2(SSM, theta_p, opts, t)
        ll_p       <- pf$Loglike
        log_alpha2 <- (ll_p - ll_sur_p) - (logLik - ll_sur_current)

        if (is.finite(log_alpha2) && log(runif(1L)) < log_alpha2) {
          theta          <- theta_p
          logLik         <- ll_p
          logPrior       <- lp_p
          X              <- pf$X
          ll_sur_current <- ll_sur_p
          n_accept       <- n_accept + 1L
        }
      }
    }
  }

  list(theta = theta, logLik = logLik, logPrior = logPrior, X = X,
       n_accept = n_accept, n_proposals = n_proposals,
       ll_sur_current = ll_sur_current)
}

# ------------------------------------------------------------------
# ADAPTIVE Nx (number of PF particles)
# ------------------------------------------------------------------
#' Increase the number of PF particles when log-likelihood variance
#' at the mean theta-particle exceeds a target threshold, and
#' re-run the bootstrap filter for every theta-particle at the new Nx.
#'
#' @param SSM state-space model list
#' @param theta_particles matrix [N_theta, d] of current theta-particles
#' @param X_list current particle-history array [N_theta, Nx, T+h, state_dim]
#' @param logLik_vec current log-likelihood per theta-particle
#' @param opts options list (must contain `N`, `sigma2_threshold`, `N_max`)
#' @param t current real-time index
#' @param T total number of observed time steps
#' @param h forecasting horizon
#' @param n_cores number of cores for parallel dispatch
#' @param use_parallel logical, use `parallel::mclapply`?
#' @return list with updated `opts`, `X_list`, `logLik`
adapt_Nx_replace <- function(SSM, theta_particles, X_list, logLik_vec,
                              opts, t, T, h, n_cores, use_parallel) {
  r_var     <- 5L
  threshold <- opts$sigma2_threshold %||% 1.5
  N_max     <- opts$N_max %||% 1000L
  max_iter  <- 3L
  N_theta   <- nrow(theta_particles)

  X_list_new <- X_list
  logLik_new <- logLik_vec

  iter <- 0L
  repeat {
    iter <- iter + 1L
    if (iter > max_iter || opts$N >= N_max) break

    theta_bar <- colMeans(theta_particles)
    ll_reps   <- replicate(r_var, PF_SMC2(SSM, theta_bar, opts, t)$Loglike)
    v         <- var(ll_reps)

    if (is.finite(v) && v <= threshold) break

    N_new <- min(N_max, ceiling(v * opts$N))
    cat(sprintf(
      "\n  up Nx increased -> %d (Var logLik = %.3f, iter = %d)", N_new, v, iter
    ))

    opts_new      <- opts
    opts_new$N    <- N_new
    opts_new$T    <- t
    opts_new$Data <- opts$Data[1:t, , drop = FALSE]

    pf_res_list <- par_lapply(seq_len(N_theta), function(i) {
      BootstrapPF(SSM, theta_particles[i, ], opts_new)
    }, n_cores = n_cores, use_parallel = use_parallel)

    X_list_new <- array(0, dim = c(N_theta, N_new, T + h, opts$state_dim))
    for (i in seq_len(N_theta)) {
      X_list_new[i, , 1:t, ] <- pf_res_list[[i]]$X[, 1:t, , drop = FALSE]
    }
    logLik_new <- sapply(pf_res_list, `[[`, "Loglike")

    opts$N <- N_new
  }

  list(opts = opts, X_list = X_list_new, logLik = logLik_new)
}
