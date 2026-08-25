# ============================================================
# smc2.R
#
# SMC^2 with Hastings-within-Gibbs (block updates), no delayed-acceptance
# until `time_DA`. No pre-epidemic padding: all time indexing is real
# time t = 1 .. T. X arrays are [N, T+h, state_dim] throughout.
#
# This file holds the algorithm's control flow only:
#   - IncrementalPF()  advance a particle cloud by one time step
#   - PF_SMC2()        run a full BPF up to a given time t
#   - SMC2()           the main outer parameter-particle loop
#
# Prior handling lives in priors.R, MCMC move kernels in
# mcmc_kernels.R, the bootstrap filter in particle_filter.R, and
# generic helpers in utils.R / resampling.R. Source all of these
# before calling SMC2() (see README.md for the load order).
# ============================================================

library(parallel)

source("src/utils.R")
source("src/resampling.R")
source("src/particle_filter.R")
source("src/mcmc_kernels.R")

# ------------------------------------------------------------------
# INCREMENTAL PF: advance the particle cloud from t-1 to t
# ------------------------------------------------------------------
#' @param X_prev particle history [N, t-1, state_dim] (real time)
#' @param Y full observation matrix, [T, ...]
#' @param t real-time step to advance to (>= 2)
#' @param SSM state-space model list
#' @param theta parameter vector
#' @param opts options list (see `BootstrapPF` for relevant fields)
#' @return list with `X` [N, t, state_dim] and scalar `logLik_incr`
IncrementalPF <- cmpfun(function(X_prev, Y, t, SSM, theta, opts) {

  N         <- opts$N
  state_dim <- opts$state_dim
  A         <- opts$A %||% 1L
  L         <- opts$L %||% 0L

  state_per_age   <- ceiling(state_dim / A)
  cols_state_list <- lapply(seq_len(A), function(a)
    ((a - 1L) * state_per_age + 1L):(a * state_per_age))

  obs_dim       <- opts$obs_dim
  obs_per_age   <- ceiling(obs_dim / A)
  cols_obs_list <- lapply(seq_len(A), function(a)
    ((a - 1L) * obs_per_age + 1L):(a * obs_per_age))

  # extend history by one column for time t
  X_full <- array(0, dim = c(N, t, state_dim))
  for (s in seq_len(t - 1L)) X_full[, s, ] <- X_prev[, s, ]

  X_full[, t, ] <- OneStepPrediction(
    X_full[, 1L:(t - 1L), , drop = FALSE],
    t, SSM, theta, opts, A, cols_state_list
  )

  upd <- OneStepUpdate(
    X_hist          = X_full,
    Y               = Y,
    t               = t,
    SSM             = SSM,
    theta           = theta,
    opts            = opts,
    A               = A,
    cols_state_list = cols_state_list,
    cols_obs_list   = cols_obs_list,
    X_obs_hist      = NULL,
    resample_method = opts$resample_method %||% "systematic",
    resample_scope  = opts$resample_scope  %||% "block",
    L               = L,
    store_X_obs     = FALSE
  )

  tmin <- max(1L, t - L)
  X_full[, tmin:t, ] <- upd$X_hist

  list(X = X_full, logLik_incr = upd$logLik_incr)
})

# ------------------------------------------------------------------
# PF WRAPPER — full BPF up to time t, with box-constraint short-circuit
# ------------------------------------------------------------------
#' @param SSM state-space model list
#' @param theta parameter vector
#' @param opts options list
#' @param t truncate `opts$Data` / `opts$T` to this time step
#' @return list from `BootstrapPF()`, or `list(Loglike = -Inf, X = NULL)`
#'   if `theta` violates `opts$lower_bounds` / `opts$upper_bounds`
PF_SMC2 <- function(SSM, theta, opts, t) {
  pf_opts      <- opts
  pf_opts$T    <- t
  pf_opts$Data <- opts$Data[1:t, , drop = FALSE]

  if (t > 1L && !is.null(opts$lower_bounds) && !is.null(opts$upper_bounds)) {
    if (any(theta < opts$lower_bounds) || any(theta > opts$upper_bounds))
      return(list(Loglike = -Inf, X = NULL))
  }

  BootstrapPF(SSM = SSM, theta = theta, opts = pf_opts)
}

# ------------------------------------------------------------------
# MAIN SMC^2 LOOP
# ------------------------------------------------------------------
#' Sequential Monte Carlo squared for the age-stratified renewal model
#'
#' @param SSM state-space model list (`InitState`, `StateProcess`, `ObsProcess`)
#' @param opts options list. In addition to the fields used by
#'   `BootstrapPF()` (including `resample_scope`), relevant fields are:
#'   \describe{
#'     \item{T}{number of observed time steps}
#'     \item{N_theta}{number of theta-particles (default 1000)}
#'     \item{N_x_init}{initial number of PF particles (default 200)}
#'     \item{ESS_threshold}{resample-move trigger, as a fraction of N_theta (default 0.5)}
#'     \item{paramPriors}{list of prior specs, see priors.R}
#'     \item{time_DA}{time step at which delayed-acceptance kicks in (default 200)}
#'     \item{k_neighbors}{neighbours for the kNN surrogate (default 3)}
#'     \item{parallel, n_cores}{parallel dispatch controls}
#'   }
#' @return list with `theta` [N_theta, T, d] history, `logEvidence`,
#'   `ESS`, `Nx`, `acceptance_rate`, and `final_particles`
SMC2 <- function(SSM, opts) {

  T            <- opts$T
  h            <- opts$forecastingHorizon %||% 0L
  N_theta      <- opts$N_theta      %||% 1000L
  ESS_thr      <- opts$ESS_threshold %||% 0.5
  d            <- length(opts$paramPriors)
  time_DA      <- opts$time_DA      %||% 200L
  k_neighbors  <- opts$k_neighbors  %||% 3L
  n_blocks     <- d
  use_parallel <- opts$parallel     %||% TRUE
  n_cores      <- opts$n_cores      %||% 4L

  # fixed-lag smoothing window
  opts$L <- max(
    length(opts$GenTime),
    length(opts$InfReportDelay),
    length(opts$InfDeathDelay %||% integer(0L)),
    50L
  ) + 1L
  L <- opts$L

  n_cores <- if (use_parallel) min(n_cores, parallel::detectCores() - 1L) else 1L

  cat(sprintf("Using %d theta-particles, %d PF-particles\n", N_theta, opts$N))
  cat(sprintf("Parallel : %s%s\n",
              if (use_parallel) "YES" else "NO",
              if (use_parallel) sprintf(" (%d cores)", n_cores) else " (sequential)"))
  cat(sprintf("Resample : scope = %s, method = %s\n",
              opts$resample_scope %||% "block", opts$resample_method %||% "systematic"))
  cat(sprintf("DA-HWG activates at t = %d\n", time_DA))

  # --- storage ---
  theta_hist       <- array(NA, dim = c(N_theta, T, d))
  logZ             <- numeric(T)
  ESS_hist         <- numeric(T)
  Nx_hist          <- numeric(T)
  accept_rate_hist <- numeric(T)

  opts$N   <- opts$N_x_init %||% 200L
  theta    <- priorSampler(N_theta, opts)
  w        <- rep(1 / N_theta, N_theta)
  logLik   <- numeric(N_theta)
  logPrior <- numeric(N_theta)

  X_list <- array(0, dim = c(N_theta, opts$N, T + h, opts$state_dim))

  # --- initialise: full BPF at t = 1 for each theta-particle ---
  pf_init <- par_lapply(seq_len(N_theta), function(i) {
    pf <- PF_SMC2(SSM, theta[i, ], opts, 1L)
    list(X        = pf$X,
         logLik   = pf$Loglike,
         logPrior = logPriorDensity(theta[i, ], opts))
  }, n_cores = n_cores, use_parallel = use_parallel)

  for (i in seq_len(N_theta)) {
    X_list[i, , 1L, ] <- pf_init[[i]]$X[, 1L, ]
    logLik[i]          <- pf_init[[i]]$logLik
    logPrior[i]        <- pf_init[[i]]$logPrior
  }

  logZ[1L]             <- mean(logLik)
  theta_hist[, 1L, ]   <- theta
  ESS_hist[1L]         <- N_theta
  Nx_hist[1L]          <- opts$N
  accept_rate_hist[1L] <- NA_real_

  start_time <- Sys.time()
  print_progress(1L, T, start_time, ESS = N_theta, Nx = opts$N, mode = "Init")

  # --- main loop t = 2 .. T ---
  for (t in 2L:T) {

    tmin <- max(1L, t - L)

    # incremental PF for every theta-particle
    pf_results <- par_lapply(seq_len(N_theta), function(i) {
      X_prev <- array(
        X_list[i, , 1L:(t - 1L), ],
        dim = c(opts$N, t - 1L, opts$state_dim)
      )
      IncrementalPF(X_prev, opts$Data, t, SSM, theta[i, ], opts)
    }, n_cores = n_cores, use_parallel = use_parallel)

    incr_logLik <- pmax(sapply(pf_results, `[[`, "logLik_incr"), -1e2)

    for (i in seq_len(N_theta)) {
      X_list[i, , tmin:t, ] <- pf_results[[i]]$X[, tmin:t, , drop = FALSE]
    }
    logLik <- logLik + incr_logLik

    # theta-particle weight update
    logZ[t] <- logZ[t - 1L] + log(sum(w * exp(incr_logLik)))

    log_w <- log(w) + incr_logLik
    if (all(!is.finite(log_w))) {
      warning(sprintf("All weights non-finite at t=%d, resetting", t))
      w <- rep(1 / N_theta, N_theta)
    } else {
      w                <- exp(log_w - max(log_w[is.finite(log_w)]))
      w[!is.finite(w)] <- 1e-6
      w                <- w / sum(w)
    }

    ESS             <- 1 / sum(w^2)
    use_DA          <- (t > time_DA)
    mode_str        <- if (use_DA) "DA-HWG" else "HWG"
    total_accepts   <- 0L
    total_proposals <- 0L

    # --- resample-move step ---
    if (ESS < ESS_thr * N_theta) {

      idx      <- Resample(w, "systematic")
      theta    <- theta[idx,   , drop = FALSE]
      logLik   <- logLik[idx]
      logPrior <- logPrior[idx]
      X_list   <- X_list[idx, , , , drop = FALSE]
      w        <- rep(1 / N_theta, N_theta)

      theta_unique  <- NULL
      logLik_unique <- NULL
      if (use_DA) {
        uniq          <- !duplicated(theta)
        theta_unique  <- theta[uniq, , drop = FALSE]
        logLik_unique <- logLik[uniq]
      }

      theta_cov   <- cov(theta)
      opts$N_move <- opts$N_move %||% 5L

      pf_move_results <- par_lapply(seq_len(N_theta), function(i) {
        X_particle <- X_list[i, , 1L:t, , drop = FALSE]

        if (use_DA) {
          DA_HWG_kernel(
            theta = theta[i, ], logLik = logLik[i], logPrior = logPrior[i],
            X = X_particle, SSM = SSM, opts = opts, theta_cov = theta_cov, t = t,
            theta_unique = theta_unique, logLik_unique = logLik_unique,
            k_neighbors = k_neighbors, n_blocks = n_blocks
          )
        } else {
          HWG_kernel(
            theta = theta[i, ], logLik = logLik[i], logPrior = logPrior[i],
            X = X_particle, SSM = SSM, opts = opts, theta_cov = theta_cov, t = t,
            n_blocks = n_blocks
          )
        }
      }, n_cores = n_cores, use_parallel = use_parallel)

      bad <- which(sapply(pf_move_results, Negate(is.list)))
      if (length(bad) > 0L) {
        cat("\nFailed workers:", bad, "\n")
        print(pf_move_results[bad])
      }

      for (i in seq_len(N_theta)) {
        res                  <- pf_move_results[[i]]
        theta[i, ]           <- res$theta
        logLik[i]            <- res$logLik
        logPrior[i]          <- res$logPrior
        X_list[i, , 1L:t, ]  <- res$X
        total_accepts        <- total_accepts   + res$n_accept
        total_proposals      <- total_proposals + res$n_proposals
      }

      # adaptive Nx
      if (opts$N < (opts$N_max %||% 1000L)) {
        adapt_res <- adapt_Nx_replace(
          SSM, theta, X_list, logLik, opts, t, T, h, n_cores, use_parallel
        )
        opts   <- adapt_res$opts
        X_list <- adapt_res$X_list
        logLik <- adapt_res$logLik
      }
    }

    idx                  <- Resample(w, "systematic")
    theta_hist[, t, ]    <- theta[idx, , drop = FALSE]
    ESS_hist[t]          <- ESS
    Nx_hist[t]           <- opts$N
    accept_rate_hist[t]  <- if (total_proposals > 0L)
      total_accepts / total_proposals
    else NA_real_

    print_progress(t, T, start_time, ESS = ESS, Nx = opts$N, mode = mode_str)
  }

  list(
    theta           = theta_hist,
    logEvidence     = logZ,
    ESS             = ESS_hist,
    Nx              = Nx_hist,
    acceptance_rate = accept_rate_hist,
    final_particles = X_list
  )
}
