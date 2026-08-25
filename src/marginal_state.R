# ============================================================
# marginal_state.R
#
# Draws posterior state trajectories (and simulated observations)
# directly from an already-completed SMC2()
#
# It is a drop-in alternative to PosteriorMarginal() for the
# age-stratified model: same two outputs (X, X_obs), just computed a
# different way. Typical call site (see run_ireland_analysis.R):
#
# ============================================================

source("src/utils.R")
source("src/epi_ssm_multiage.R")

#' Map age group `a` onto its block of state columns
#'
#' @param a age-group index (1-indexed)
#' @param A total number of age groups
#' @param state_dim total state dimension
#' @return integer vector of column indices for age group `a`
.state_cols <- function(a, A, state_dim) {
  per <- ceiling(state_dim / A)
  ((a - 1L) * per + 1L):(a * per)
}

#' Average reporting fraction per day-of-week, estimated from observed data
#'
#' Self-contained fallback used only when `opts$precompute_dow` is TRUE
#' and `opts$omega_dow` wasn't supplied: for each age column, divides
#' each day-of-week's mean reported count by the overall mean, giving a
#' multiplicative day-of-week effect centred near 1.
#'
#' @param Y_hist data frame, first column Date (or numeric day index),
#'   remaining columns reported counts per age group
#' @param A number of age groups
#' @param max_weeks cap on how many trailing weeks of data to use
#' @return numeric matrix [7 x A], rows Mon..Sun
.precompute_dow_weights <- function(Y_hist, A, max_weeks = 16L) {
  n_use <- min(nrow(Y_hist), max_weeks * 7L)
  tail_rows <- utils::tail(Y_hist, n_use)
  
  raw_day <- tail_rows[[1]]
  dow <- if (inherits(raw_day, "Date")) {
    as.integer(format(raw_day, "%u"))
  } else {
    ((as.integer(raw_day) - 1L) %% 7L) + 1L
  }
  
  omega <- matrix(1, 7L, A)
  for (a in seq_len(A)) {
    counts    <- tail_rows[[a + 1L]]
    day_means <- tapply(counts, dow, mean, na.rm = TRUE)
    overall   <- mean(counts, na.rm = TRUE)
    if (is.finite(overall) && overall > 0) {
      omega[as.integer(names(day_means)), a] <- day_means / overall
    }
  }
  omega
}

############################################################
## STEP 1: MARGINAL STATE SAMPLES (with optional forecast
## extension baked in, driven by opts$forecastingHorizon)
##
## X_marginal: (N_s * Nx) x T_out x state_dim
##   T_out = T                      if forecastingHorizon <= 0
##   T_out = T + forecastingHorizon otherwise
##
## The h-step forecast is pure forward simulation through
## StateProcess only (no reweighting/resampling), since there is no
## data beyond opts$T to condition on.
############################################################

#' Sample posterior state trajectories directly from an SMC2() result
#'
#' @param result output of `SMC2()`; must contain `theta` and `final_particles`
#' @param opts options list (same as used for the `SMC2()` call, plus
#'   optionally an updated `forecastingHorizon`)
#' @param N_s number of theta-particles to sub-sample (default 100)
#' @return list with `X_marginal` ((N_s * Nx) x T_out x state_dim),
#'   `m_idx` (sampled theta-particle indices), and `theta_sub`
#'   (N_s x T x d, unchanged/un-forecast)
get_marginal_states <- function(result, opts, N_s = 100) {
  
  theta_hist <- result$theta            # N_theta x T x d
  X_all      <- result$final_particles  # N_theta x Nx x T x state_dim (T == opts$T; see file header)
  
  N_theta   <- dim(X_all)[1]
  Nx        <- dim(X_all)[2]
  T         <- dim(X_all)[3]
  state_dim <- dim(X_all)[4]
  A         <- opts$A
  
  h     <- opts$forecastingHorizon %||% 0
  T_out <- T + max(h, 0)
  
  m_idx <- sample(seq_len(N_theta), N_s, replace = FALSE)
  
  X_marginal <- array(NA_real_, dim = c(N_s * Nx, T_out, state_dim))
  theta_sub  <- theta_hist[m_idx, , , drop = FALSE]  # N_s x T x d
  
  cols_state_list <- if (h > 0) lapply(seq_len(A), function(a) .state_cols(a, A, state_dim)) else NULL
  
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
        for (a in seq_len(A)) {
          cols <- cols_state_list[[a]]
          X_marginal[row_start:row_end, t, cols] <- EpiSSM$StateProcess(
            X_marginal[row_start:row_end, 1:(t - 1L), , drop = FALSE],
            t, theta_t, opts, cols, a
          )
        }
      }
    }
  }
  
  list(
    X_marginal = X_marginal,  # (N_s * Nx) x T_out x state_dim
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
## beyond the fitted horizon, and expected cases are scaled by the
## age-specific reporting rate rho_a (see file header for the
## assumed theta layout).
##
## Returns: C_sim, (N_s * Nx) x T_out x A
############################################################

#' Simulate reported-case observations from sampled state trajectories
#'
#' @param X_marginal output of `get_marginal_states()$X_marginal`
#' @param theta_sub output of `get_marginal_states()$theta_sub`
#' @param opts options list. Relevant fields: `A`, `Data`, `InfReportDelay`,
#'   `T`, `precompute_dow`, `omega_dow` (optional, else estimated
#'   from `Data`), `week_effect` (used only as a fallback if
#'   `precompute_dow` is FALSE)
#' @return list with `C_sim`, (N_s * Nx) x T_out x A
simulate_observations <- function(X_marginal, theta_sub, opts) {
  
  N_total <- dim(X_marginal)[1]
  T_out   <- dim(X_marginal)[2]
  A       <- opts$A
  N_s     <- dim(theta_sub)[1]
  Nx      <- N_total / N_s
  Y_hist  <- opts$Data
  f_c     <- opts$InfReportDelay
  Tfit    <- opts$T
  precompute_dow <- isTRUE(opts$precompute_dow)
  
 
  
  if (precompute_dow && is.null(opts$omega_dow)) {
    opts$omega_dow <- .precompute_dow_weights(Y_hist, A, max_weeks = 16L)
  }
  
  C_sim <- array(NA_real_, dim = c(N_total, T_out, A))
  
  for (s in seq_len(N_s)) {
    
    row_start <- (s - 1L) * Nx + 1L
    row_end   <- s * Nx
    theta_t   <- theta_sub[s, Tfit, ]  # frozen theta, matches the forecast branch in get_marginal_states()
    
    
    
    for (a in seq_len(A)) {
      
      if (!is.null(opts$CAR)) {
        rho_a <- CAR[a]
        kappa   <- theta_t[A + 1]
      } else if (length(theta_t) > 2 * A) {
        rho_a <- theta_t[A + a]
        kappa   <- theta_t[2 * A + 1]
      } else {
        rho_a <- 1
        kappa   <- theta_t[A + 1]
      }
      
      I_idx  <- 2L * (a - 1L) + 2L
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
          opts$omega_dow[day_of_week, a]
        } else if (isTRUE(opts$week_effect)) {
          omega_idx <- (A + 2L):(A + 7L)
          omega_raw <- theta_t[omega_idx]
          S6        <- sum(omega_raw)
          omega     <- c(7 * omega_raw / (S6 + 1), 7 / (S6 + 1))
          omega[day_of_week]
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
        mu_case <- pmax(rho_a * mu_infections, 0)
        mu_c    <- omega_t * mu_case
        mu_c[!is.finite(mu_c)] <- 0
        
        # -------------------------------------------------------
        # Sample negative binomial observations
        # -------------------------------------------------------
        c_draws <- rnbinom(Nx, size = kappa, mu = mu_c)
        c_draws[is.na(c_draws)] <- 0
        
        C_sim[row_start:row_end, t, a] <- c_draws
      }
    }
  }
  
  list(C_sim = C_sim)
}
