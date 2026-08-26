# ============================================================
# run_multi_simulation.R  
#
# Fits the age-stratified EpiSSM to SIMULATED data via SMC^2 and
# validates the fit against known ground-truth trajectories
# (beta_true, I_true).
#
# Simulated data is generated once via simulate_epissm_multi() and then
# cached to CSV under data/simulated/ (see src/simulated_data_cache.R):
# subsequent runs just reload the cached CSVs instead of re-simulating.
# Delete data/simulated/sim_multiage_* to force a fresh simulation.
#
# Run from the repository root:  Rscript scripts/run_multi_simulation.R
# ============================================================

library(ggplot2)
library(tidyverse)
library(patchwork)
library(cowplot)

source("src/smc2.R")
source("src/epi_ssm_multiage.R")
source("src/sim_data_multi.R")
source("src/simulated_data_cache.R")
source("src/diagnostics.R")
source("src/epi_diagnostics.R")
source("src/posterior_marginal.R")
source("src/marginal_state.R")

dir.create("figures/sim_data", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# 1. Simulated data (cached to CSV after the first run)
# ------------------------------------------------------------
sim <- load_or_simulate("sim_multiage", function() simulate_epissm_multi(T = 100))

Y_df       <- sim$Y_df
Y_df$Date  <- as.Date(Y_df$Date)
Y          <- as.matrix(Y_df[, -1])
A          <- sim$A
N_4group   <- sim$N_pop
M_4x4      <- sim$ContMatrix
gen_pmf        <- sim$gen_pmf
case_delay_pmf <- sim$case_delay_pmf
beta_true  <- sim$beta_true
I_true     <- sim$I_true

# ------------------------------------------------------------
# 2. SMC^2 options
# ------------------------------------------------------------
opts <- list(
  Data              = Y_df,
  T                 = nrow(Y_df),
  A                 = A,
  state_dim         = 2 * A,
  obs_dim           = A,
  N_pop             = N_4group,
  resample_method   = "systematic",
  resample_scope    = "block",
  L                 = 31,
  forecastingHorizon   = 0,
  PredictedObservation = FALSE,
  
  GenTime        = gen_pmf,
  InfReportDelay = case_delay_pmf,
  
  chol_contact     = FALSE,
  scale_contact    = FALSE,
  estimate_contact = FALSE,
  C_syth           = M_4x4,
  ContMatrix       = M_4x4,
  
  week_effect    = FALSE,
  precompute_dow = FALSE,
  
  pBeta0 = replicate(A, function(N) rbeta(N, shape1 = 1, shape2 = 10), simplify = FALSE),
  pI0    = replicate(A, function(N) sample(1:200, N, replace = TRUE), simplify = FALSE),
  
  # sigma[a] : RW volatility for log-transmissibility in age group a
  # kappa    : negative-binomial overdispersion for case counts
  paramNames = c(
    expression(sigma[1]), expression(sigma[2]),
    expression(sigma[3]), expression(sigma[4]),
    expression(kappa)
  ),
  paramPriors = list(
    list(dist = "lnorm", meanlog = log(0.02), sdlog = 0.5),
    list(dist = "lnorm", meanlog = log(0.02), sdlog = 0.5),
    list(dist = "lnorm", meanlog = log(0.02), sdlog = 0.5),
    list(dist = "lnorm", meanlog = log(0.02), sdlog = 0.5),
    list(dist = "exp",   rate = 0.1)
  ),
  lower_bounds = c(0, 0, 0, 0,   0),
  upper_bounds = c(1, 1, 1, 1, 100),
  param_blocks = list(c(1), c(2), c(3), c(4), c(5)),
  
  N_x_init         = 300,
  N_theta          = 800,
  ESS_threshold    = 0.4,
  k_neighbors      = 3,
  sigma2_threshold = 1.5,
  N_max            = 1000,
  n_cores          = 20,
  parallel         = TRUE
)

# ------------------------------------------------------------
# 3. Run SMC^2
# ------------------------------------------------------------
result <- SMC2(SSM = EpiSSM, opts = opts)
saveRDS(result, file = "data/simulated/smc2_multiage_sim_result.rds")
# result <-  readRDS("data/real/smc2_multiage_sim_result.rds") # load saved result
# ------------------------------------------------------------
# 4. Algorithm diagnostics
# ------------------------------------------------------------
fig_diagnostic <- plot_SMC2_diagnostics(result)
ggsave("figures/sim_data/fig_sim_diagnostic.pdf", fig_diagnostic, width = 24, height = 14)
print(fig_diagnostic)

# ------------------------------------------------------------
# 5. Parameter filtering estimates
# ------------------------------------------------------------
Nparam      <- length(opts$paramNames)
param_plots <- lapply(seq_len(Nparam), function(j) {
  plot_filter_one(theta_mat = result$theta[, , j], y_label = opts$paramNames[[j]])
})
fig_param <- .grid_shared_legend(param_plots, ncol = 2)
ggsave("figures/sim_data/fig_sim_param.pdf", fig_param,
       width = 24, height = 6 * ceiling(Nparam / 2))
print(fig_param)


# ------------------------------------------------------------
# 6. Posterior marginal state 
# ------------------------------------------------------------

#source("src/marginal_state.R")
opts$forecastingHorizon <- 0
marg       <- get_marginal_states(result, opts, N_s = 100)
X_marginal <- marg$X_marginal
theta_sub  <- marg$theta_sub
obs_sim    <- simulate_observations(X_marginal, theta_sub, opts)
C_sim      <- obs_sim$C_sim


#============================================
# we can also use final theta to run the BPF
#============================================

# Ns        <- 100
# theta_sub <- result$theta[sample(seq_len(opts$N_theta), Ns), opts$T, ]
# 
# opts$precompute_dow <- FALSE
# marg       <- PosteriorMarginal(theta_sub, EpiSSM, opts)
# X_marginal <- marg$X        # (Ns x Nx) x (T + h) x state_dim
# C_sim      <- marg$X_obs    # (Ns x Nx) x (T + h) x A

age_labels <- c("0-24", "25-44", "45-64", "65+")
age_colors <- c("0-24" = "#2166ac", "25-44" = "#4dac26",
                "45-64" = "#8B6914", "65+" = "#762a83")
n_col <- 2
w_age <- 12 * n_col
h_age <- 8 * ceiling(A / n_col)

# ------------------------------------------------------------
# 7. Age-specific transmissibility (beta) vs. ground truth
# ------------------------------------------------------------
beta_plots <- lapply(seq_len(A), function(a) {
  plot_filter_one_age(
    theta_mat       = X_marginal[, , 2 * a - 1],
    init_date       = as.character(min(Y_df$Date)),
    plot_title      = paste("Age group:", age_labels[a]),
    age_color       = age_colors[age_labels[a]],
    true_trajectory = beta_true[, a]
  )
})
fig_beta <- .assemble_age_fig(
  beta_plots, n_col,
  bquote(bold("Age-specific susceptibility parameter") ~ ~ bolditalic(beta)[a * "," ~ t])
)
ggsave("figures/sim_data/fig_sim_beta.pdf", fig_beta, width = w_age, height = h_age)
print(fig_beta)

# ------------------------------------------------------------
# 8. Effective reproduction number R_{a,t} / R_t vs. ground truth
# ------------------------------------------------------------
rt_res      <- compute_Rt_Rat(X_marginal[, , 2 * seq_len(A) - 1], M_4x4)
rt_res_true <- compute_Rt_Rat_single(beta_true, M_4x4)

rat_plots <- lapply(seq_len(A), function(a) {
  plot_filter_one_age(
    theta_mat       = rt_res$Rat[, , a],
    plot_title      = paste("Age group:", age_labels[a]),
    age_color       = age_colors[age_labels[a]],
    true_trajectory = rt_res_true$Rat[, a]
  ) + geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.7)
})
fig_rat <- .assemble_age_fig(
  rat_plots, n_col,
  expression(bold("Age-specific reproduction number") ~ ~ bolditalic(R)[a * "," ~ t])
)
ggsave("figures/sim_data/fig_sim_Rat.pdf", fig_rat, width = w_age, height = h_age)
print(fig_rat)

# ------------------------------------------------------------
# 9. New infections (I) by age group vs. ground truth
# ------------------------------------------------------------
I_plots <- lapply(seq_len(A), function(a) {
  plot_filter_one_age(
    theta_mat       = X_marginal[, 1:opts$T, 2 * (a - 1) + 2],
    plot_title      = paste("Age group:", age_labels[a]),
    age_color       = age_colors[age_labels[a]],
    true_trajectory = I_true[, a]
  )
})
fig_I <- .assemble_age_fig(I_plots, n_col, "New infections")
ggsave("figures/sim_data/fig_sim_I.pdf", fig_I, width = w_age, height = h_age)
print(fig_I)

# ------------------------------------------------------------
# 10. Reported cases by age group vs. simulated observations
# ------------------------------------------------------------
C_plots <- lapply(seq_len(A), function(a) {
  plot_filter_one_age(
    theta_mat   = C_sim[, , a],
    plot_title  = paste("Age group:", age_labels[a]),
    age_color   = age_colors[age_labels[a]],
    observed_points = Y[, a]
  )
})
fig_C <- .assemble_age_fig(C_plots, n_col, "Reported cases")
ggsave("figures/sim_data/fig_sim_C.pdf", fig_C, width = w_age, height = h_age)
print(fig_C)

# ------------------------------------------------------------
# 11. Population-level totals: Rt, total infections, total cases
# ------------------------------------------------------------
T_obs        <- opts$T
I_total_mat  <- sort_then_sum(lapply(seq_len(A), function(a) X_marginal[, 1:T_obs, 2 * (a - 1) + 2]))
I_total_true <- rowSums(I_true)
C_total_mat  <- sort_then_sum(lapply(seq_len(A), function(a) C_sim[, 1:T_obs, a]))
C_total_true <- rowSums(Y)

p_Rt <- plot_filter_one(theta_mat = rt_res$Rt, plot_title = "Effective reproduction number") +
  geom_hline(yintercept = 1, linetype = "dotted", linewidth = 1.0) +
  geom_line(data = data.frame(Time = seq_len(T_obs), value = rt_res_true$Rt),
            aes(x = Time, y = value), color = "red", linewidth = 1.3, inherit.aes = FALSE)

p_I <- plot_filter_one(theta_mat = I_total_mat, plot_title = "New infections") +
  geom_line(data = data.frame(Time = seq_len(T_obs), value = I_total_true),
            aes(x = Time, y = value), color = "red", linewidth = 1.3, inherit.aes = FALSE)

p_C <- plot_filter_one(theta_mat = C_total_mat, plot_title = "Reported cases") +
  geom_point(data = data.frame(Time = seq_len(T_obs), value = C_total_true),
             aes(x = Time, y = value), color = "red", size = 1.8, alpha = 0.8, inherit.aes = FALSE)

fig_total <- .grid_shared_legend(list(p_I, p_Rt, p_C), ncol = 1)
ggsave("figures/sim_data/fig_sim_total.pdf", fig_total, width = 12, height = 15)
print(fig_total)

