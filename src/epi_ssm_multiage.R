# ============================================================
# epi_ssm_multiage.R  
#
# Age-structured epidemic state-space model.
#
# State per particle per time step: [beta_1, I_1, beta_2, I_2, ..., beta_A, I_A]
#   beta_a  - log-scale transmission rate for age group a
#   I_a     - new infections in age group a
#
# Observation model: negative-binomial on reported cases, with optional
# day-of-week (DOW) multiplier and infection-to-report delay convolution.
#
# Key design decisions
# --------------------
#  1. Forecast-safe throughout: every read of Y_hist is gated on t <= opts$T;
#     DOW for t > opts$T is inferred from the last observed date.
#  2. Renewal force computed via matrix multiply (no nested R loops).
#  3. Cholesky base matrix precomputed once in opts (call precompute_chol_base).
#  4. All index arithmetic is centralised in small helpers to avoid drift.
# ============================================================

library(stats)
source("src/utils.R")

# ============================================================
# COLUMN / PARAMETER INDEX HELPERS
# ============================================================
# theta layout: [(DOW/contact extras)]
.dow_start   <- function(opts) 2L * opts$A + 2L        # first DOW weight, just past kappa

.extra_start <- function(opts) {                       # first "extra" (contact) param
  base <- 2L * opts$A + 1L                              # A sigmas + A rhos + kappa
  if (isTRUE(opts$week_effect) && !isTRUE(opts$precompute_dow)) base <- base + 6L
  base + 1L
}

.n_chol      <- function(A) as.integer(A * (A + 1L) / 2L)
.state_cols  <- function(a) c(2L * (a - 1L) + 1L, 2L * (a - 1L) + 2L)
.all_I_cols  <- function(A) seq(2L, 2L * A, by = 2L)

# ============================================================
# DAY-OF-WEEK HELPERS
# ============================================================

#' Extract day-of-week (1=Mon ... 7=Sun) from Y_hist row t.
#' Y_hist[, 1] may be Date or a numeric day index.
.dow_from_yhist <- function(Y_hist, t) {
  raw <- Y_hist[t, 1L]
  if (inherits(raw, "Date"))  return(as.integer(format(raw, "%u")))
  if (is.numeric(raw))        return(((as.integer(raw) - 1L) %% 7L) + 1L)
  stop("Y_hist[,1] must be Date or numeric.")
}

#' Infer DOW for a forecast step t > opts$T by advancing from the last
#' observed date. Never reads Y_hist beyond row opts$T.
.dow_forecast <- function(Y_hist, t, opts) {
  offset <- t - opts$T
  last   <- Y_hist[opts$T, 1L]
  if (inherits(last, "Date"))
    return(as.integer(format(last + offset, "%u")))
  if (is.numeric(last)) {
    last_dow <- ((as.integer(last) - 1L) %% 7L) + 1L
    return(((last_dow - 1L + offset) %% 7L) + 1L)
  }
  stop("Y_hist[,1] must be Date or numeric to infer forecast DOW.")
}

#' Unified DOW lookup: works for both observed and forecast steps.
.get_dow <- function(Y_hist, t, opts) {
  if (t <= opts$T) .dow_from_yhist(Y_hist, t) else .dow_forecast(Y_hist, t, opts)
}

#' Compute pre-averaged DOW weight matrix from recent observed data.
#'
#' @param Y_hist observation data; col 1 = date/day index, cols 2:(A+1) = counts
#' @param A number of age groups
#' @param max_weeks maximum number of trailing complete weeks to average over
#' @return `7 x A` matrix; rows = Mon-Sun, cols = age groups, each column sums to 7
compute_dow_weights <- function(Y_hist, A, max_weeks = 16L) {
  T_obs   <- nrow(Y_hist)
  dow_vec <- vapply(seq_len(T_obs), function(i) .dow_from_yhist(Y_hist, i), integer(1L))

  max_days <- (min(T_obs, max_weeks * 7L) %/% 7L) * 7L
  if (max_days == 0L) {
    warning("compute_dow_weights: fewer than 7 observations; returning uniform weights.")
    return(matrix(1, 7L, A))
  }

  idx     <- (T_obs - max_days + 1L):T_obs
  Y_win   <- Y_hist[idx, 2L:(A + 1L), drop = FALSE]
  dow_win <- dow_vec[idx]

  omega <- matrix(NA_real_, 7L, A)
  for (a in seq_len(A)) {
    x_a             <- as.numeric(Y_win[, a])
    x_a[is.na(x_a)] <- 0
    total           <- sum(x_a)
    if (total <= 0) { omega[, a] <- 1; next }
    for (d in 1:7) omega[d, a] <- 7 * sum(x_a[dow_win == d]) / total
  }
  omega
}

# ============================================================
# CONTACT MATRIX HELPERS  (Cholesky NCP reparameterisation)
# ============================================================

#' Precompute the Cholesky base matrix once. Call after building opts.
#' @param opts options list; must contain `N_pop`, `C_syth`
#' @return `opts` with `L_syth` added
precompute_chol_base <- function(opts) {
  opts$L_syth <- t(chol(diag(opts$N_pop) %*% opts$C_syth))
  opts
}

#' Build the effective contact matrix from Cholesky perturbation parameters.
#' Uses `opts$L_syth` if precomputed; falls back to computing it on the fly.
.chol_contact_matrix <- function(L_tilde_vec, opts) {
  A      <- opts$A
  L_syth <- opts$L_syth %||% t(chol(diag(opts$N_pop) %*% opts$C_syth))

  L_tilde <- matrix(0, A, A)
  L_tilde[lower.tri(L_tilde, diag = TRUE)] <- L_tilde_vec

  L <- L_syth + 0.05 * L_syth * L_tilde
  diag(1 / opts$N_pop) %*% (L %*% t(L))
}

#' Extract contact row a from theta given the chosen parameterisation.
.contact_row <- function(a, theta, opts) {
  off <- .extra_start(opts)
  A   <- opts$A

  if (isTRUE(opts$chol_contact)) {
    nc <- .n_chol(A)
    return(.chol_contact_matrix(theta[off:(off + nc - 1L)], opts)[a, ])
  }
  if (isTRUE(opts$scale_contact))
    return(theta[off + a - 1L] * opts$ContMatrix[a, ])
  if (isTRUE(opts$estimate_contact))
    return(theta[off + (a - 1L) * A + seq(0L, A - 1L)])

  opts$ContMatrix[a, ]  # fixed contact matrix
}

# ============================================================
# 1. InitState
# ============================================================
#' @param N number of particles
#' @param theta parameter vector (unused here, kept for interface consistency)
#' @param opts options list; must contain `A`, `pBeta0`, and, if `Y_hist`
#'   is NULL, `pI0`
#' @param Y_hist optional observation history used to centre each age
#'   group's I_0 prior on its early observed case counts
#' @return matrix `[N, 2*A]`
InitState <- function(N, theta, opts, Y_hist = NULL, ...) {
  A  <- opts$A
  X0 <- matrix(0, nrow = N, ncol = 2L * A)
  
  for (a in seq_len(A)) {
    cols     <- .state_cols(a)
    beta_col <- cols[1L]
    I_col    <- cols[2L]
    
    X0[, beta_col] <- opts$pBeta0[[a]](N)
    
    if (!is.null(opts$CAR)) {
      rho_a <- CAR[a]
    } else if (length(theta) > 2 * A) {
      rho_a <- theta[A + a]
    } else {
      rho_a <- 1
    }
    
    if (!is.null(Y_hist)) {
      Y_a             <- as.numeric(Y_hist[1:min(5L, nrow(Y_hist)), a + 1L])
      Y_a[is.na(Y_a)] <- 0
      X0[, I_col]     <- rpois(N, max(mean(Y_a) / rho_a, 1))
    } else {
      X0[, I_col] <- opts$pI0[[a]](N)
    }
  }
  X0
}
# ============================================================
# 2. StateProcess
# ============================================================
#' @param X_hist particle history array, `[N, t-1, state_dim]`
#' @param t current real-time step (>= 2)
#' @param theta parameter vector; `theta[a]` = RW2 volatility for age a
#' @param opts options list; must contain `A`, `GenTime`, `ContMatrix`, `T`
#' @param cols state columns for this age group (`(beta_col, I_col)`)
#' @param a current age-group index
#' @return matrix `[N, 2]` with columns `beta`, `I`
StateProcess <- function(X_hist, t, theta, opts, cols, a) {
  N        <- dim(X_hist)[1L]
  A        <- opts$A
  beta_col <- cols[1L]
  I_cols   <- .all_I_cols(A)

  # ---- Transmission rate (random-walk on log scale) ----
  sigma_a   <- theta[a]
  beta_prev <- X_hist[, t - 1L, beta_col]

  beta_log <- if (t <= 2L) {
    rnorm(N, log(beta_prev), sigma_a)
  } else if (t > opts$T) {
    log(beta_prev)  # forecast: dampen volatility (freeze at last level)
  } else {
    beta_prev2 <- X_hist[, t - 2L, beta_col]
    rnorm(N, 2 * log(beta_prev) - log(beta_prev2), sigma_a)
  }
  beta <- exp(beta_log)

  # ---- Renewal force via matrix multiply ----
  # Lambda_a[n] = sum_b C[b,a] * sum_s g[s] * I_b(t - s)
  g       <- opts$GenTime
  max_lag <- min(t - 1L, length(g))
  g_sub   <- g[seq_len(max_lag)]
  lag_idx <- seq(t - 1L, by = -1L, length.out = max_lag)

  Lambda_mat <- matrix(0, N, max_lag)
  for (b in seq_len(A)) {
    I_b_mat    <- matrix(X_hist[, lag_idx, I_cols[b]], N, max_lag)
    Lambda_mat <- Lambda_mat + opts$ContMatrix[b, a] * I_b_mat
  }
  Lambda <- as.vector(Lambda_mat %*% (g_sub / sum(g_sub)))

  # ---- New infections ----
  mu_I                    <- pmax(beta * Lambda, 0)
  mu_I[!is.finite(mu_I)] <- 0

  cbind(beta = beta, I = rpois(N, mu_I))
}

# ============================================================
# 3. ObsProcess  -  forecast-safe
# ============================================================
#' @param X_hist particle history array, `[N, t, state_dim]`
#' @param Y_hist observation data; col 1 = date/day index, cols 2:(A+1) = counts
#' @param t current real-time step
#' @param theta parameter vector; 
#' @param opts options list
#' @param a current age-group index
#' @return list with `logw` (log-weight, zero during forecast) and
#'   `X_obs` (predicted observations, or NULL if not requested)
ObsProcess <- function(X_hist, Y_hist, t, theta, opts, a, ...) {
  N      <- dim(X_hist)[1L]
  I_hist <- matrix(X_hist[, , 2L], nrow = N)  # [N x T_stored]
  in_obs <- (t >= 1L && t <= opts$T)
  A      <- opts$A
  
  # ---- Day-of-week multiplier ----
  omega_t <- if (isTRUE(opts$precompute_dow)) {
    if (is.null(opts$omega_dow))
      stop("opts$omega_dow must be set when precompute_dow = TRUE.")
    opts$omega_dow[.get_dow(Y_hist, t, opts), a]
    
  } else if (isTRUE(opts$week_effect)) {
    omega_raw <- theta[.dow_start(opts):(.dow_start(opts) + 5L)]
    S6        <- sum(omega_raw)
    omega_7   <- c(7 * omega_raw / (S6 + 1), 7 / (S6 + 1))
    omega_7[.get_dow(Y_hist, t, opts)]
    
  } else {
    1
  }
  
  # ---- Expected reported cases ----
  f_c       <- opts$InfReportDelay
  max_delay <- min(t - 1L, length(f_c))
  delay_idx <- pmax(seq(t - 1L, by = -1L, length.out = max(max_delay, 1L)), 1L)
  f_sub     <- f_c[seq_len(max(max_delay, 1L))]
  
  mu_infections <- as.vector(I_hist[, delay_idx, drop = FALSE] %*% (f_sub / sum(f_sub)))
  
  if (!is.null(opts$CAR)) {
    rho_a <- CAR[a]
    kappa   <- theta[A + 1]
  } else if (length(theta) > 2 * A) {
    rho_a <- theta[A + a]
    kappa   <- theta[2 * A + 1]
  } else {
    rho_a <- 1
    kappa   <- theta[A + 1]
  }
  
  mu_case       <- pmax(rho_a * mu_infections, 0)
  
  # ---- Log-weights (zero during forecast - no observations) ----
  logw <- rep(0, N)
  if (in_obs) {
    Y_case <- Y_hist[t, a + 1L]
    if (!is.na(Y_case))
      logw <- dnbinom(x = Y_case, size = kappa, mu = omega_t * mu_case, log = TRUE)
  }
  
  # ---- Predicted observations (optional) ----
  X_obs <- if (isTRUE(opts$PredictedObservation))
    rnbinom(N, size = kappa, mu = omega_t * mu_case)
  else NULL
  
  list(logw = logw, X_obs = X_obs)
}

# ============================================================
# Assemble SSM module
# ============================================================
EpiSSM <- list(
  InitState            = InitState,
  StateProcess         = StateProcess,
  ObsProcess           = ObsProcess,
  compute_dow_weights  = compute_dow_weights,
  precompute_chol_base = precompute_chol_base
)
