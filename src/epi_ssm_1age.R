# ============================================================
# epi_ssm_1age.R 
#
# Single-population (non-age-stratified) epidemic state-space model,
# compatible with BootstrapPF() / SMC2().
#
#   X_t ~ p(X_t | X_{1:t-1}, theta)      state: (R_t, I_t)
#   Y_t ~ p(Y_t | X_{1:t}, theta)        obs:   reported cases
#
# theta layout:
#   theta[1]      sigma  - RW2 volatility for log(R_t)
#   week_effect == TRUE:
#     theta[2:7]  omega_1..6 - raw day-of-week weights (7th derived)
#     theta[8]    kappa  - negative-binomial overdispersion
#   week_effect == FALSE:
#     theta[2]    kappa  - negative-binomial overdispersion
#
# All indexing is plain real time t = 1 .. T (no padding).
# ============================================================

library(stats)
source("src/utils.R")
source("src/delay_distributions.R")

# ------------------------------------------------------------
# 1. Initial state sampler
# ------------------------------------------------------------
#' @param N number of particles
#' @param theta parameter vector (unused here, kept for interface consistency)
#' @param opts options list; must contain `pR0` (list of 1 sampler function)
#'   and, if `Y_hist` is NULL, `pI0`
#' @param Y_hist optional observation history used to centre the I_0 prior
#'   on the first observed case count (mean-matched log-normal-Poisson mix)
#' @return matrix `[N, 2]` with columns `R`, `I`
InitState <- function(N, theta, opts, Y_hist = NULL, ...) {
  X0 <- matrix(0, nrow = N, ncol = 2L)

  if (!is.null(Y_hist)) {
    Y_1  <- Y_hist[1L, 2L]
    Y_1[is.na(Y_1)] <- 0
    mean_I <- max(mean(Y_1, na.rm = TRUE), 1)
    iota   <- rnorm(N, mean = log(mean_I), sd = 2)
    X0[, 2L] <- rpois(N, exp(iota))
  } else {
    X0[, 2L] <- opts$pI0[[1L]](N)
  }

  X0[, 1L] <- opts$pR0[[1L]](N)
  X0
}

# ------------------------------------------------------------
# 2. State transition process (RW2 on log R_t; renewal equation for I_t)
# ------------------------------------------------------------
#' @param X_hist particle history array, `[N, t-1, state_dim]`
#' @param t current real-time step (>= 2)
#' @param theta parameter vector (`theta[1]` = sigma)
#' @param opts options list; must contain `GenTime`
#' @return matrix `[N, 2]` with columns `R`, `I`
StateProcess <- function(X_hist, t, theta, opts, ...) {
  N <- dim(X_hist)[1L]

  R_prev <- X_hist[, t - 1L, 1L]

  # ---- R_t evolution: second-order random walk on log scale ----
  R <- if (t == 2L) {
    exp(rnorm(N, log(R_prev), theta[1L]))
  } else {
    R_prev2 <- X_hist[, t - 2L, 1L]
    exp(rnorm(N, 2 * log(R_prev) - log(R_prev2), theta[1L]))
  }

  # ---- Renewal equation ----
  g       <- opts$GenTime
  max_lag <- min(max(1L, t - 1L), length(g))
  lag_idx <- pmax(seq(from = max(1L, t - 1L), by = -1L, length.out = max_lag), 1L)
  g_sub   <- g[seq_len(max_lag)]

  I_hist <- X_hist[, , 2L, drop = TRUE]
  if (is.vector(I_hist)) I_hist <- matrix(I_hist, nrow = N)

  Lambda <- as.vector(I_hist[, lag_idx, drop = FALSE] %*% (g_sub / sum(g_sub)))

  cbind(R = R, I = rpois(N, R * Lambda))
}

# ------------------------------------------------------------
# 3. Observation process: negative-binomial reported cases,
#    with optional day-of-week reporting effect
# ------------------------------------------------------------
#' @param X_hist particle history array, `[N, t, state_dim]`
#' @param Y_hist observation data frame/matrix (col 1 = date/day index,
#'   col 2 = reported cases)
#' @param t current real-time step
#' @param theta parameter vector (see file header for layout)
#' @param opts options list; must contain `T`, `InfReportDelay`, and
#'   optionally `week_effect`
#' @return list with `logw` (log-weight, 0 outside the observed window)
#'   and `X_obs` (matrix `[N, 1]` of simulated reported cases)
ObsProcess <- function(X_hist, Y_hist, t, theta, opts, ...) {
  N           <- dim(X_hist)[1L]
  week_effect <- opts$week_effect %||% FALSE

  d         <- opts$InfReportDelay
  max_delay <- min(max(1L, t - 1L), length(d))
  lag_idx   <- pmax(seq(from = max(1L, t - 1L), by = -1L, length.out = max_delay), 1L)
  d_sub     <- d[seq_len(max_delay)]

  I_hist <- X_hist[, , 2L, drop = TRUE]
  if (is.vector(I_hist)) I_hist <- matrix(I_hist, nrow = N)

  # ---- Day-of-week multiplier ----
  if (week_effect) {
    omega_raw <- theta[2L:7L]
    S6        <- sum(omega_raw)
    omega     <- c(7 * omega_raw / (S6 + 1), 7 / (S6 + 1))

    raw_day <- Y_hist[t, 1L]
    day_of_week <- if (inherits(raw_day, "Date")) {
      as.integer(format(raw_day, "%u"))
    } else if (is.numeric(raw_day)) {
      ((as.integer(raw_day) - 1L) %% 7L) + 1L
    } else {
      stop("Date column in Y_hist must be Date or numeric")
    }

    omega_t <- omega[day_of_week]
    kappa   <- theta[8L]
  } else {
    omega_t <- 1
    kappa   <- theta[2L]
  }

  # ---- Expected reported cases ----
  mu <- if (t <= opts$T) {
    as.vector(I_hist[, lag_idx, drop = FALSE] %*% (d_sub / sum(d_sub)))
  } else {
    I_hist[, t, drop = FALSE]
  }
  mu <- pmax(mu, 0)

  C_obs <- rnbinom(N, size = kappa, mu = omega_t * mu)

  # ---- Log-weights (zero outside the observed window) ----
  logw <- 0
  if (t >= 1L && t <= opts$T) {
    Y_t  <- Y_hist[t, 2L]
    logw <- dnbinom(x = Y_t, size = kappa, mu = omega_t * mu, log = TRUE)
  }

  list(logw = logw, X_obs = matrix(C_obs, nrow = N, ncol = 1L))
}

# ============================================================
# Assemble SSM module
# ============================================================
EpiSSM_1age <- list(
  InitState    = InitState,
  StateProcess = StateProcess,
  ObsProcess   = ObsProcess
)
