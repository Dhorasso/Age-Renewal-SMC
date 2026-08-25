# ============================================================
# sim_data_multi.R
#
# Synthetic age-stratified renewal dataset for the simulation study,
# ported from the original `generate_synthetic_data_age4_epinow_style_cases.R`.
#
# This is a faithful port, not a re-derivation: the beta_true
# trajectories (oscillating cosine/sine curves with hand-tuned
# per-age-group constants), the exact renewal-equation loop for
# I_true, and the exact negative-binomial observation model are all
# unchanged from the original. Only the following are adapted to fit
# this repo's structure:
#   - population/contact-matrix loading now goes through
#     `read_age_distribution()` / `build_contact_matrix()`
#     (contact_matrix_utils.R) instead of reading the raw PEA11 CSV
#     directly - both produce the same reciprocity-corrected,
#     age-aggregated matrix, since Ireland_country_level_age_distribution_85.csv
#     is already the same 85-band population vector the original
#     script derived from PEA11.
#   - `discretize_gamma()` -> `discretize_dist(..., dist = "gamma")`
#     (the consolidated equivalent in delay_distributions.R); values
#     are identical.
#   - wrapped in a function returning the same named list the rest of
#     the pipeline expects, so it can be cached via
#     `load_or_simulate()` in simulated_data_cache.R.
#   - the loop variable `t` was renamed `t_idx` to avoid shadowing R's
#     `t()` (transpose); no other change to the logic.
# ============================================================

source("src/delay_distributions.R")
source("src/contact_matrix_utils.R")

#' Forward-simulate the age-4 synthetic dataset (exact port of the
#' original generator's beta_true / I_true / Y logic)
#'
#' @param T number of days to simulate (original default: 100)
#' @param groups named list of age-group index vectors into the 85
#'   single-year-of-age contact matrix (default: four 20-year bands)
#' @param I0 length-A vector, new infections at t = 1 for each age group
#'   (original default: c(3, 8, 5, 2))
#' @param kappa_true negative-binomial overdispersion for reported cases
#'   (original default: 15)
#' @param rho reporting-fraction multiplier on the observation mean
#'   (original default: 1)
#' @param start_date first date of the simulated series (original: "2020-01-01")
#' @param contact_dir directory containing the Ireland contact-matrix CSVs
#' @param seed random seed (original default: 123)
#' @return list with `Y_df` (data frame: Date + reported cases per age
#'   group), `Y` (matrix version, no Date column), `gen_pmf`,
#'   `case_delay_pmf`, `beta_true` (`T x A`), `I_true` (`T x A`),
#'   `ContMatrix`, `N_pop`, `kappa_true`, `A`, `T`
simulate_epissm_multi <- function(T = 100L,
                                  groups = list(G1 = 1:25, G2 = 26:45,
                                                G3 = 46:65, G4 = 66:85),
                                  I0          = c(3, 8, 5, 2),
                                  kappa_true  = 15,
                                  rho         = 1,
                                  start_date  = "2020-01-01",
                                  contact_dir = "data/contact_matrices",
                                  seed        = 123L) {
  
  set.seed(seed)
  A <- length(groups)
  stopifnot(length(I0) == A)
  
  # ---- contact matrix + population (same reciprocity + aggregation math
  #      as the original script's steps 2-4) ----
  # Population by single year of age, 0–99
  N_pop_full <- read.csv("data/contact_matrices/PEA11.20260324T160347.csv")
  N_pop_full <- N_pop_full$VALUE
  # Collapse ages 85+ into a single open-ended group → length-85 vector
  N_pop <- c(N_pop_full[1:84], sum(N_pop_full[85:length(N_pop_full)]))
  
  # 85×85 country-level contact matrix (Mistry al. 2021)
  F_fine     <- as.matrix(read.csv(
    "data/contact_matrices/Ireland_country_level_M_overall_contact_matrix_85.csv",
    header = FALSE
  ))
  
  # Aggregation +  Enforce reciprocity
  cm       <- build_contact_matrix(F_fine, N_pop, groups)
  M_4x4    <- cm$M
  N_4group <- cm$N_group
  M_4x4 
  
  # ---- delay kernels (same distributional choices as the original) ----
  gen_pmf    <- discretize_dist(max = 14, mean = 3.3, sd = 2.1, dist = "gamma")
  incub_pmf  <- discretize_dist(max = 14, mean = 3.0, sd = 2.3, dist = "gamma")
  report_pmf <- discretize_dist(max = 10, mean = 2.0, sd = 1.0, dist = "gamma")
  
  # kept exactly as the original: no post-hoc clipping/renormalization here -
  # each day's window is renormalized locally in the observation loop below
  case_delay_pmf <- convolve(incub_pmf, rev(report_pmf), type = "open")
  
  # ---- true beta trajectories (exact form from the original script) ----
  t_idx      <- seq_len(T)
  beta_true  <- matrix(NA_real_, T, A)
  beta_const <- 0.25 * c(0.8, 0.75, 0.76, 1.1)
  pp         <- c(55, 60, 60, 45)
  
  for (a in seq_len(A)) {
    A_amp <- 0.5 + 0.15 * a
    P1    <- pp[a]
    P2    <- 180 - 4 * a
    t0    <- 3 + a
    if (a == 4) t0 <- 15
    tt <- pmax(t_idx - t0, 0)
    
    if (a == 2) {
      beta_true[, a] <- beta_const[a] * exp(
        A_amp * (sin(2 * pi / 1 * tt / P1) - 1.3) - tt / P2
      )
    } else {
      beta_true[, a] <- beta_const[a] * exp(
        A_amp * (cos(2 * pi / (max(1, 0.37 * a)) * tt / P1) -
                   max(1.05, 0.4 * a)) - tt / P2
      )
    }
  }
  
  # ---- latent infections: exact renewal equation from the original ----
  I_true <- matrix(0, T, A)
  I_true[1, ] <- I0
  
  # tilde M_ab = M_ab , via transpose (as in the original)
  M_tilde <- t(M_4x4)
  
  for (ti in 2:T) {
    infection_pressure <- rep(0, A)
    
    max_lag <- max(1, min(length(gen_pmf), ti - 1))
    g <- gen_pmf[seq_len(max_lag)]
    
    for (s in seq_len(max_lag)) {
      infectious_vector <- I_true[max(1, ti - s), ]
      contact_term      <- M_tilde %*% infectious_vector
      infection_pressure <- infection_pressure + (g[s] / sum(g)) * contact_term
    }
    
    I_true[ti, ] <- rpois(A, pmax(beta_true[ti, ] * infection_pressure, 1e-10))
  }
  
  # ---- reported cases: exact observation model from the original ----
  Y <- matrix(0, T, A)
  for (ti in seq_len(T)) {
    max_delay <- max(1, min(ti - 1, length(case_delay_pmf)))
    cdelay    <- case_delay_pmf[seq_len(max_delay)]
    
    if (max_delay > 0) {
      mu <- rep(0, A)
      for (d in seq_len(max_delay)) {
        mu <- mu + I_true[max(1, ti - d), ] * cdelay[d] / sum(cdelay)
      }
      Y[ti, ] <- rnbinom(A, size = kappa_true, mu = rho * mu)
    }
  }
  
  Y_df <- data.frame(Date = seq.Date(as.Date(start_date), by = "day", length.out = T), Y)
  
  list(
    Y_df           = Y_df,
    Y              = Y,
    gen_pmf        = gen_pmf,
    case_delay_pmf = case_delay_pmf,
    beta_true      = beta_true,
    I_true         = I_true,
    ContMatrix     = M_4x4,
    N_pop          = N_4group,
    kappa_true     = kappa_true,
    A              = A,
    T              = T
  )
}
