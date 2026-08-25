# ============================================================
# particle_filter.R 
#
# Bootstrap particle filter (BPF) for the age-stratified renewal
# model. X arrays have dimensions [N, T+h, state_dim]; all indexing
# is in plain real time t = 1 .. T+h (no pre-epidemic padding).
#
# Age-structure is handled via `cols_state_list` / `cols_obs_list`,
# which map each age group a = 1..A onto its block of state / obs
# columns.
#
# Two resampling scopes are supported through `opts$resample_scope`:
#   "block"  (default) - each age group is resampled independently,
#                         using only its own observation weights.
#                         This improves particle efficiency by letting
#                         age groups replenish diversity separately.
#   "global" - a single set of indices, drawn from the *joint*
#                         (product-across-age) weight, is applied to
#                         the entire state vector. This preserves
#                         full posterior coherence across age groups
#                         at the cost of faster degeneracy.
#
# ============================================================

source("src/utils.R")
source("src/resampling.R")
library(compiler)

# ------------------------------------------------------------------
# ONE-STEP PREDICTION  (real time t)
# ------------------------------------------------------------------
#' Propagate all age blocks forward one time step
#'
#' @param X_hist particle history array, [N, t-1, state_dim]
#' @param t real-time index to predict into
#' @param SSM state-space model list; must provide `StateProcess()`
#' @param theta parameter vector
#' @param opts options list (must contain `N`, `state_dim`)
#' @param A number of age groups
#' @param cols_state_list list (length A) of state column indices per age
#' @return matrix [N, state_dim], the predicted state at time t
OneStepPrediction <- function(X_hist, t, SSM, theta, opts, A, cols_state_list) {
  N         <- opts$N
  state_dim <- opts$state_dim
  X_t       <- matrix(0, N, state_dim)

  for (a in seq_len(A)) {
    cols        <- cols_state_list[[a]]
    X_t[, cols] <- SSM$StateProcess(X_hist, t, theta, opts, cols, a)
  }
  X_t
}
OneStepPrediction <- cmpfun(OneStepPrediction)

# ------------------------------------------------------------------
# ONE-STEP UPDATE  (real time t)
# ------------------------------------------------------------------
#' Weight, resample and (optionally) window the particle history
#'
#' @param X_hist particle history array, [N, t, state_dim]
#' @param Y observation data
#' @param t real-time index to update at
#' @param SSM state-space model list; must provide `ObsProcess()`
#' @param theta parameter vector
#' @param opts options list
#' @param A number of age groups
#' @param cols_state_list list (length A) of state column indices per age
#' @param cols_obs_list list (length A) of observation column indices per age
#' @param X_obs_hist optional predicted-observation history array to update
#' @param resample_method "systematic" (default), "multinomial", or "stratified"
#' @param resample_scope "block" (default, per-age resampling) or "global"
#'   (single joint resampling step across all age groups)
#' @param L fixed-lag smoothing window (0 = no smoothing lag)
#' @param store_X_obs logical, whether to populate `X_obs_hist`
#' @return list with `X_hist`, `X_obs_hist`, `logLik_incr`, `logw_age`
OneStepUpdate <- function(X_hist, Y, t, SSM, theta, opts, A,
                           cols_state_list, cols_obs_list, X_obs_hist,
                           resample_method = "systematic",
                           resample_scope  = c("block", "global"),
                           L = 0, store_X_obs = FALSE) {

  resample_scope <- match.arg(resample_scope)

  N         <- opts$N
  state_dim <- opts$state_dim
  obs_dim   <- opts$obs_dim

  logw_age    <- vector("list", A)
  X_obs_t     <- if (store_X_obs) matrix(0, N, obs_dim) else NULL
  logLik_incr <- 0

  # --- per-age observation weights ---
  for (a in seq_len(A)) {
    cols_state <- cols_state_list[[a]]
    cols_obs   <- cols_obs_list[[a]]

    obs <- SSM$ObsProcess(
      X_hist[, 1:t, cols_state, drop = FALSE],
      Y, t, theta, opts, a
    )

    logw_age[[a]] <- obs$logw
    if (store_X_obs) X_obs_hist[, t, cols_obs] <- obs$X_obs
  }

  tmin <- max(1L, t - L)  # fixed-lag window

  if (t <= opts$T) {

    # Log-likelihood increment is identical regardless of resampling
    # scope: it is the sum, across age groups, of the log-mean-weight
    # for that group's observation.
    logLik_incr <- sum(vapply(logw_age, function(lw) {
      m <- max(lw)
      m + log(sum(exp(lw - m))) - log(N)
    }, numeric(1)))

    if (resample_scope == "global") {
      # --- single joint resampling step across all age groups ---
      logw_total <- Reduce(`+`, logw_age)
      m          <- max(logw_total)
      w_vec      <- exp(logw_total - m)
      w_vec      <- w_vec / sum(w_vec)
      indx       <- Resample(w_vec, method = resample_method)

      X_hist[, tmin:t, ] <- X_hist[indx, tmin:t, , drop = FALSE]
      if (store_X_obs) {
        X_obs_hist[, tmin:t, ] <- X_obs_hist[indx, tmin:t, , drop = FALSE]
      }

    } else {
      # --- independent resampling within each age block ---
      for (a in seq_len(A)) {
        lw    <- logw_age[[a]]
        m     <- max(lw)
        w_vec <- exp(lw - m)
        w_vec <- w_vec / sum(w_vec)

        indx       <- Resample(w_vec, method = resample_method)
        cols_state <- cols_state_list[[a]]
        X_hist[, tmin:t, cols_state] <-
          X_hist[indx, tmin:t, cols_state, drop = FALSE]

        if (store_X_obs) {
          cols_obs                      <- cols_obs_list[[a]]
          X_obs_hist[, tmin:t, cols_obs] <-
            X_obs_hist[indx, tmin:t, cols_obs, drop = FALSE]
        }
      }
    }
  }

  list(
    X_hist      = X_hist[, tmin:t, , drop = FALSE],
    X_obs_hist  = if (!is.null(X_obs_hist))
      X_obs_hist[, tmin:t, , drop = FALSE]
    else NULL,
    logLik_incr = logLik_incr,
    logw_age    = logw_age
  )
}
OneStepUpdate <- cmpfun(OneStepUpdate)

# ------------------------------------------------------------------
# FULL BOOTSTRAP PF
# ------------------------------------------------------------------
#' Run a full bootstrap particle filter over t = 1 .. T (+ forecast horizon)
#'
#' @param SSM state-space model list (`InitState`, `StateProcess`, `ObsProcess`)
#' @param theta parameter vector
#' @param opts options list. Relevant fields:
#'   \describe{
#'     \item{N}{number of particles}
#'     \item{T}{number of observed time steps}
#'     \item{forecastingHorizon}{extra steps to simulate beyond T (default 0)}
#'     \item{L}{fixed-lag smoothing window (default 0)}
#'     \item{A}{number of age groups (default 1)}
#'     \item{state_dim, obs_dim}{total state / observation dimension}
#'     \item{PredictedObservation}{store predicted observations? (default FALSE)}
#'     \item{resample_method}{"systematic" (default), "multinomial", "stratified"}
#'     \item{resample_scope}{"block" (default) or "global"; see file header}
#'   }
#' @param X_init optional initial particle matrix [N, state_dim]; otherwise
#'   sampled from `SSM$InitState()`
#' @param ... passed through to `SSM$InitState()`
#' @return list with `X` [N, T+h, state_dim], `Loglike`, and optionally `X_obs`
BootstrapPF <- function(SSM, theta, opts, X_init = NULL, ...) {

  # --- options ---
  Y         <- opts$Data
  N         <- opts$N
  T         <- opts$T
  h         <- opts$forecastingHorizon %||% 0
  L         <- opts$L %||% 0
  state_dim <- opts$state_dim
  obs_dim   <- opts$obs_dim
  A         <- opts$A %||% 1

  store_X_obs     <- opts$PredictedObservation %||% FALSE
  resample_method <- opts$resample_method %||% "systematic"
  resample_scope  <- opts$resample_scope  %||% "block"

  # --- column index helpers ---
  state_per_age <- ceiling(state_dim / A)
  obs_per_age   <- ceiling(obs_dim   / A)

  cols_state_list <- lapply(seq_len(A), function(a)
    ((a - 1) * state_per_age + 1):(a * state_per_age))
  cols_obs_list   <- lapply(seq_len(A), function(a)
    ((a - 1) * obs_per_age + 1):(a * obs_per_age))

  opts$I_cols <- seq(2, 2 * A, by = 2)

  # --- pre-allocate storage ---
  X       <- array(0, dim = c(N, T + h, state_dim))
  X_obs   <- if (store_X_obs) array(0, dim = c(N, T + h, obs_dim)) else NULL
  Loglike <- 0

  # -- t = 1: initialise + update --
  X_t1     <- if (is.null(X_init)) SSM$InitState(N, theta, opts, Y, ...) else X_init
  X[, 1, ] <- X_t1

  upd1 <- OneStepUpdate(
    X[, 1, , drop = FALSE], Y, 1L, SSM, theta, opts, A,
    cols_state_list, cols_obs_list, X_obs,
    resample_method, resample_scope, L, store_X_obs
  )

  X[, 1, ] <- upd1$X_hist[, 1, ]
  if (store_X_obs && !is.null(upd1$X_obs_hist)) X_obs[, 1, ] <- upd1$X_obs_hist[, 1, ]
  Loglike <- Loglike + upd1$logLik_incr

  # -- t = 2 .. T: main filtering loop --
  if (T > 1) {
    for (t in 2:T) {
      tmin <- max(1L, t - L)

      X[, t, ] <- OneStepPrediction(
        X[, 1:(t - 1), , drop = FALSE], t, SSM, theta, opts, A, cols_state_list
      )

      upd <- OneStepUpdate(
        X[, 1:t, , drop = FALSE], Y, t, SSM, theta, opts, A,
        cols_state_list, cols_obs_list, X_obs,
        resample_method, resample_scope, L, store_X_obs
      )

      X[, tmin:t, ] <- upd$X_hist
      if (store_X_obs && !is.null(upd$X_obs_hist)) X_obs[, tmin:t, ] <- upd$X_obs_hist
      Loglike <- Loglike + upd$logLik_incr
    }
  }

  # -- t = T+1 .. T+h: forecasting horizon (no further weighting) --
  if (h > 0) {
    for (t in (T + 1):(T + h)) {
      tmin <- max(1L, t - L)

      X[, t, ] <- OneStepPrediction(
        X[, 1:(t - 1), , drop = FALSE], t, SSM, theta, opts, A, cols_state_list
      )

      upd <- OneStepUpdate(
        X[, 1:t, , drop = FALSE], Y, t, SSM, theta, opts, A,
        cols_state_list, cols_obs_list, X_obs,
        resample_method, resample_scope, L, store_X_obs
      )

      X[, tmin:t, ] <- upd$X_hist
      if (store_X_obs && !is.null(upd$X_obs_hist)) X_obs[, tmin:t, ] <- upd$X_obs_hist
    }
  }

  res <- list(X = X, Loglike = Loglike)
  if (store_X_obs) res$X_obs <- X_obs
  res
}
BootstrapPF <- cmpfun(BootstrapPF)
