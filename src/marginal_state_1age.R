# ============================================================
# marginal_state_1age.R
#
# Draws posterior state trajectories (and simulated observations)
# directly from an already-completed SMC2() run of the
# SINGLE-POPULATION (non-age-stratified) model in epi_ssm_1age.R.
#
# It is a drop-in alternative to PosteriorMarginal() for that model:
# same two outputs (X, X_obs), just computed a different way.
#
# This is the 1-age analogue of marginal_state.R. Differences from
# that file, driven by epi_ssm_1age.R:
#   - state_dim is fixed at 2 (columns R, I) - no age blocks, so
#     .state_cols()/A are gone entirely.
#   - StateProcess() is called with the same signature it actually
#     has in epi_ssm_1age.R: StateProcess(X_hist, t, theta, opts).
#   - theta layout has no per-age reporting rate rho_a. Expected
#     cases are the raw renewal-equation convolution, exactly as in
#     ObsProcess().
#   - kappa / the day-of-week weights come from theta[8]/theta[2:7]
#     when opts$week_effect is TRUE, or theta[2] when it's FALSE -
#     matching ObsProcess() rather than the multi-age theta layout.
# ============================================================

source("src/utils.R")
source("src/epi_ssm_1age.R")

#' Average reporting fraction per day-of-week, estimated from observed data
#'
#' Self-contained fallback used only when `opts$precompute_dow` is TRUE
#' and `opts$omega_dow` wasn't supplied: divides each day-of-week's mean
#' reported count by the overall mean, giving a multiplicative
#' day-of-week effect centred near 1.
#'
#' @param Y_hist data frame, col 1 = Date (or numeric day index),
#'   col 2 = reported counts
#' @param max_weeks cap on how many trailing weeks of data to use
#' @return numeric vector length 7, Mon..Sun
.precompute_dow_weights_1age <- function(Y_hist, max_weeks = 16L) {
  n_use     <- min(nrow(Y_hist), max_weeks * 7L)
  tail_rows <- utils::tail(Y_hist, n_use)
  
  raw_day <- tail_rows[[1]]
  dow <- if (inherits(raw_day, "Date")) {
    as.integer(format(raw_day, "%u"))
  } else {
    ((as.integer(raw_day) - 1L) %% 7L) + 1L
  }
  
  omega   <- rep(1, 7L)
  counts  <- tail_rows[[2]]
  day_means <- tapply(counts, dow, mean, na.rm = TRUE)
  overall   <- mean(counts, na.rm = TRUE)
  if (is.finite(overall) && overall > 0) {
    omega[as.integer(names(day_means))] <- day_means / overall
  }
  omega
}

############################################################
## STEP 1: MARGINAL STATE SAMPLES (with optional forecast
## extension baked in, driven by opts$forecastingHorizon)
##
## X_marginal: (N_s * Nx) x T_out x 2   (columns: R, I)
##   T_out = T                      if forecastingHorizon <= 0
##   T_out = T + forecastingHorizon otherwise
##
## The h-step forecast is pure forward simulation through
## StateProcess only (no reweighting/resampling), since there is no
## data beyond opts$T to condition on.
############################################################

#' Sample posterior state trajectories directly from an SMC2() result
#' (single-population model)
#'
#' @param result output of `SMC2()`; must contain `theta` and `final_particles`
#' @param opts options list (same as used for the `SMC2()` call, plus
#'   optionally an updated `forecastingHorizon`)
#' @param N_s number of theta-particles to sub-sample (default 100)
#' @return list with `X_marginal` ((N_s * Nx) x T_out x 2),
#'   `m_idx` (sampled theta-particle indices), and `theta_sub`
#'   (N_s x T x d, unchanged/un-forecast)
get_marginal_states_1age <- function(result, opts, N_s = 100) {
  
  theta_hist <- result$theta            # N_theta x T x d
  X_all      <- result$final_particles  # N_theta x Nx x T x 2 (T == opts$T)
  
  N_theta   <- dim(X_all)[1]
  Nx        <- dim(X_all)[2]
  T         <- dim(X_all)[3]
  state_dim <- dim(X_all)[4]  # == 2 (R, I)
  
  h     <- opts$forecastingHorizon %||% 0
  T_out <- T + max(h, 0)
  
  m_idx <- sample(seq_len(N_theta), N_s, replace = FALSE)
  
  X_marginal <- array(NA_real_, dim = c(N_s * Nx, T_out, state_dim))
  theta_sub  <- theta_hist[m_idx, , , drop = FALSE]  # N_s x T x d
  
  for (s in seq_len(N_s)) {
    m         <- m_idx[s]
    row_start <- (s - 1L) * Nx + 1L
    row_end   <- s * Nx
    
    # ---- fitted part: 1:T, copied straight from the SMC2 particle cloud ----
    X_marginal[row_start:row_end, 1:T, ] <- X_all[m, , , ]
    
    # ---- forecast part: T+1 .. T+h, pure forward simulation ----
    if (h > 0) {
      theta_t <- theta_sub[s, T, ]  # freeze theta at the last filtering time
      
      for (t in (T + 1L):T_out) {
        X_marginal[row_start:row_end, t, ] <- EpiSSM_1age$StateProcess(
          X_marginal[row_start:row_end, 1:(t - 1L), , drop = FALSE],
          t, theta_t, opts
        )
      }
    }
  }
  
  list(
    X_marginal = X_marginal,  # (N_s * Nx) x T_out x 2
    m_idx      = m_idx,
    theta_sub  = theta_sub    # N_s x T x d (theta itself is never forecast, only states are)
  )
}

############################################################
## STEP 2: SIMULATE OBSERVATIONS FROM MARGINAL STATE SAMPLES
##
## Runs over the full T_out (= dim(X_marginal)[2]), forecast-safe:
## day-of-week for t > opts$T is inferred by advancing from the last
## observed date, theta is frozen at theta_sub[s, opts$T, ] for any t
## beyond the fitted horizon. No reporting-rate scaling (rho == 1),
## matching ObsProcess() in epi_ssm_1age.R.
##
## Returns: C_sim, (N_s * Nx) x T_out x 1
############################################################

#' Simulate reported-case observations from sampled state trajectories
#' (single-population model)
#'
#' @param X_marginal output of `get_marginal_states_1age()$X_marginal`
#' @param theta_sub output of `get_marginal_states_1age()$theta_sub`
#' @param opts options list. Relevant fields: `Data`, `InfReportDelay`,
#'   `T`, `week_effect`, `precompute_dow`, `omega_dow` (optional, else
#'   estimated from `Data`)
#' @return list with `C_sim`, (N_s * Nx) x T_out x 1
simulate_observations_1age <- function(X_marginal, theta_sub, opts) {
  
  N_total <- dim(X_marginal)[1]
  T_out   <- dim(X_marginal)[2]
  N_s     <- dim(theta_sub)[1]
  Nx      <- N_total / N_s
  Y_hist  <- opts$Data
  f_c     <- opts$InfReportDelay
  Tfit    <- opts$T
  week_effect    <- opts$week_effect %||% FALSE
  precompute_dow <- isTRUE(opts$precompute_dow)
  
  if (precompute_dow && is.null(opts$omega_dow)) {
    opts$omega_dow <- .precompute_dow_weights_1age(Y_hist, max_weeks = 16L)
  }
  
  C_sim <- array(NA_real_, dim = c(N_total, T_out, 1L))
  
  I_idx <- 2L  # column layout is (R, I)
  
  for (s in seq_len(N_s)) {
    
    row_start <- (s - 1L) * Nx + 1L
    row_end   <- s * Nx
    theta_t   <- theta_sub[s, Tfit, ]  # frozen theta, matches the forecast branch in get_marginal_states_1age()
    
    if (week_effect) {
      omega_raw <- theta_t[2L:7L]
      S6        <- sum(omega_raw)
      omega_wk  <- c(7 * omega_raw / (S6 + 1), 7 / (S6 + 1))
      kappa     <- max(theta_t[8L], 1e-4)
    } else {
      kappa <- max(theta_t[2L], 1e-4)
    }
    
    I_real <- matrix(X_marginal[row_start:row_end, 1:T_out, I_idx], nrow = Nx, ncol = T_out)
    I_real <- pmax(I_real, 0)
    I_real[!is.finite(I_real)] <- 0
    
    for (t in seq_len(T_out)) {
      
      # -------------------------------------------------------
      # Day-of-week, forecast-safe
      # -------------------------------------------------------
      if (t <= Tfit) {
        raw_day <- Y_hist[t, 1]
        day_of_week <- if (inherits(raw_day, "Date")) {
          as.integer(format(raw_day, "%u"))
        } else if (is.numeric(raw_day)) {
          ((as.integer(raw_day) - 1L) %% 7L) + 1L
        } else {
          stop("Y_hist[,1] must be Date or numeric")
        }
      } else {
        offset <- t - Tfit
        last   <- Y_hist[Tfit, 1]
        day_of_week <- if (inherits(last, "Date")) {
          as.integer(format(last + offset, "%u"))
        } else if (is.numeric(last)) {
          last_dow <- ((as.integer(last) - 1L) %% 7L) + 1L
          ((last_dow - 1L + offset) %% 7L) + 1L
        } else {
          stop("Y_hist[,1] must be Date or numeric")
        }
      }
      
      omega_t <- if (precompute_dow) {
        opts$omega_dow[day_of_week]
      } else if (week_effect) {
        omega_wk[day_of_week]
      } else {
        1
      }
      
      # -------------------------------------------------------
      # Expected reported cases: same convolution as ObsProcess()
      # -- no padding; delay_idx is simply clamped to >= 1, which
      # reuses the earliest available I when history is short.
      # -------------------------------------------------------
      max_delay <- min(t - 1L, length(f_c))
      delay_idx <- pmax(seq(t - 1L, by = -1L, length.out = max(max_delay, 1L)), 1L)
      f_sub     <- f_c[seq_len(max(max_delay, 1L))]
      
      mu_infections <- as.vector(I_real[, delay_idx, drop = FALSE] %*% (f_sub / sum(f_sub)))
      mu_c <- pmax(omega_t * mu_infections, 0)
      mu_c[!is.finite(mu_c)] <- 0
      
      # -------------------------------------------------------
      # Sample negative binomial observations
      # -------------------------------------------------------
      c_draws <- rnbinom(Nx, size = kappa, mu = mu_c)
      c_draws[is.na(c_draws)] <- 0
      
      C_sim[row_start:row_end, t, 1L] <- c_draws
    }
  }
  
  list(C_sim = C_sim)
}
