# ============================================================
# run_pmmh_ireland.R
#
# Example: fit the age-stratified EpiSSM to real Irish COVID-19 data
# via PMMH, as an alternative to the sequential SMC^2 fit in
# `run_ireland_analysis.R`.
#
# This reuses the *exact same* SSM, delay kernels, contact matrix, and
# most of the `opts` list as run_ireland_analysis.R — only the
# MCMC-specific opts fields (iterations, nChains, burnin, thin) are
# added, and SMC2(...) is swapped for PMMH(...). The same pattern
# applies to any other model/data combination already set up for
# SMC2() in this repo (e.g. EpiSSM_1age on the real 1-age series, or
# EpiSSM on the simulated multi-age dataset) — just build `opts` the
# same way you already do for that SMC2 script, add the MCMC fields
# below, and call PMMH() instead of SMC2().
#
# Run from the repository root:  Rscript scripts/run_pmmh_ireland.R
# ============================================================

library(tidyverse)
library(BayesianTools)

source("src/pmmh.R")
source("src/epi_ssm_multiage.R")
source("src/delay_distributions.R")
source("src/contact_matrix_utils.R")

dir.create("figures/real_data", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# 1. Load and window the real case data (identical to run_ireland_analysis.R)
# ------------------------------------------------------------
covid_data      <- read.csv("data/real/Covid_Data_Ireland.csv")
covid_data$Date <- as.Date(covid_data$Date)

start_date <- as.Date("2022-03-01")
end_date   <- as.Date("2022-07-17")
covid      <- covid_data %>% filter(Date >= start_date & Date <= end_date)

# ------------------------------------------------------------
# 2. Delay kernels + contact matrix (identical to run_ireland_analysis.R)
# ------------------------------------------------------------
gen_pmf        <- discretize_dist(max = 14, mean = 3.3, sd = 2.1, dist = "gamma")
incub_pmf      <- discretize_dist(max = 14, mean = 3.0, sd = 2.3, dist = "gamma")
report_pmf     <- discretize_dist(max = 10, mean = 2.0, sd = 1.0, dist = "gamma")
case_delay_pmf <- convolve(incub_pmf, rev(report_pmf), type = "open")
case_delay_pmf <- pmax(case_delay_pmf, 0)
case_delay_pmf <- case_delay_pmf / sum(case_delay_pmf)

groups <- list(G1 = 1:25, G2 = 26:45, G3 = 46:65, G4 = 66:85)
A      <- length(groups)

N_pop_full <- read_age_distribution("data/contact_matrices/Ireland_country_level_age_distribution_85.csv")
F_fine     <- as.matrix(read.csv(
  "data/contact_matrices/Ireland_country_level_M_overall_contact_matrix_85.csv",
  header = FALSE
))
cm       <- build_contact_matrix(F_fine, N_pop_full, groups)
M_4x4    <- cm$M
N_4group <- cm$N_group

omega_dow_mat <- compute_dow_weights(covid, A, max_weeks = 16)

# ------------------------------------------------------------
# 3. opts: identical to run_ireland_analysis.R, plus PMMH-only fields
# ------------------------------------------------------------
opts <- list(
  Data              = covid,
  T                 = nrow(covid),
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
  week_effect      = FALSE,
  precompute_dow   = TRUE,
  omega_dow        = omega_dow_mat,

  pBeta0 = replicate(A, function(N) rbeta(N, shape1 = 1, shape2 = 10), simplify = FALSE),
  pI0    = replicate(A, function(N) sample(400:1200, N, replace = TRUE), simplify = FALSE),

  paramNames = c("sigma[1]", "sigma[2]", "sigma[3]", "sigma[4]", "kappa"),
  paramPriors = list(
    list(dist = "lnorm", meanlog = log(0.02), sdlog = 0.5),
    list(dist = "lnorm", meanlog = log(0.02), sdlog = 0.5),
    list(dist = "lnorm", meanlog = log(0.02), sdlog = 0.5),
    list(dist = "lnorm", meanlog = log(0.02), sdlog = 0.5),
    list(dist = "exp",   rate = 0.1)
  ),
  lower_bounds = c(0, 0, 0, 0,   0),
  upper_bounds = c(1, 1, 1, 1, 200),

  # ── PF particle count (used by BootstrapPF - same field SMC2() reads) ──────
  N = 1000,

  # ── PMMH-only fields ────────────────────────────────────────────────────────
  iterations = 20000,  # total PMMH iterations
  nChains    = 3,       # DEzs requires >= 3 chains
  burnin     = 4500,    # discarded before convergence checks / posterior draws
  thin       = 3,
  message    = TRUE
)

# ------------------------------------------------------------
# 4. Sanity check: single BPF run at a fixed theta before committing to PMMH
# ------------------------------------------------------------
theta_fixed <- c(0.02, 0.02, 0.02, 0.02, 10)

t0  <- Sys.time()
res <- BootstrapPF(SSM = EpiSSM, theta = theta_fixed, opts = opts)
cat("BPF runtime:", format(Sys.time() - t0), "\n")
cat("Log-likelihood at theta_fixed:", res$Loglike, "\n")

# ------------------------------------------------------------
# 5. Run PMMH
# ------------------------------------------------------------
mcmcResults <- PMMH(EpiSSM, opts)
saveRDS(mcmcResults, file = "data/real/pmmh_ireland_result.rds")

# mcmcResults <- readRDS("data/real/pmmh_ireland_result.rds")

# ------------------------------------------------------------
# 6. Convergence diagnostics (BayesianTools' own plotting methods)
# ------------------------------------------------------------
pdf("figures/real_data/fig_pmmh_diagnostics.pdf", width = 10, height = 8)
plot(mcmcResults)
correlationPlot(mcmcResults)
marginalPlot(mcmcResults)
dev.off()

print(summary(mcmcResults))

# ------------------------------------------------------------
# 7. Posterior samples, after burn-in
# ------------------------------------------------------------
thetaSamples <- getSample(
  mcmcResults,
  parametersOnly = TRUE,
  thin           = 1,
  burnin         = opts$burnin,
  numSamples     = 100,
  coda           = FALSE
)



#============================================
# Marginal state estimate
#============================================

opts$forecastingHorizon <- 14
marg       <- PosteriorMarginal(thetaSamples, EpiSSM, opts)
X_marginal <- marg$X        # (Ns x Nx) x (T + h) x state_dim
C_sim      <- marg$X_obs    # (Ns x Nx) x (T + h) x A

age_labels <- c("0-24", "25-44", "45-64", "65+")
age_colors <- c("0-24" = "#2166ac", "25-44" = "#4dac26",
                "45-64" = "#8B6914", "65+" = "#762a83")
n_col <- 2
w_age <- 12 * n_col
h_age <- 8 * ceiling(A / n_col)

# ------------------------------------------------------------
# 9. Age-specific transmissibility (beta)
# ------------------------------------------------------------
beta_plots <- lapply(seq_len(A), function(a) {
  plot_filter_one_age(
    theta_mat            = X_marginal[, , 2 * a - 1],
    init_date            = as.character(start_date),
    plot_title           = paste("Age group:", age_labels[a]),
    age_color            = age_colors[age_labels[a]],
    forecast_start_date  = forecast_start,
    ymin = 0.02, ymax = 0.25
  )
})
fig_beta <- .assemble_age_fig(
  beta_plots, n_col,
  y_label = bquote(bold("Age-specific susceptibility parameter") ~ ~ bolditalic(beta)[a * "," ~ t]),
  x_label = "Date"
)
ggsave("figures/real_data/fig_age_beta.pdf", fig_beta, width = w_age, height = h_age)
print(fig_beta)

# ------------------------------------------------------------
# 10. New infections (I) by age group
# ------------------------------------------------------------
I_plots <- lapply(seq_len(A), function(a) {
  plot_filter_one_age(
    theta_mat           = X_marginal[, , 2 * (a - 1) + 2],
    init_date           = as.character(start_date),
    plot_title          = paste("Age group:", age_labels[a]),
    age_color           = age_colors[age_labels[a]],
    forecast_start_date = forecast_start
  )
})
fig_I <- .assemble_age_fig(I_plots, n_col, "New infections")
ggsave("figures/real_data/fig_age_I.pdf", fig_I, width = w_age, height = h_age)
print(fig_I)

# ------------------------------------------------------------
# 11. Reported cases (C) by age group: train fit + held-out forecast
# ------------------------------------------------------------
C_plots <- lapply(seq_len(A), function(a) {
  train_series <- covid_data %>% filter(Date >= start_date & Date <= cut_date) %>% pull(a + 1)
  test_series  <- covid_data %>% filter(Date > cut_date & Date <= cut_date + horizon) %>% pull(a + 1)
  
  plot_filter_one_age(
    theta_mat           = C_sim[, , a],
    init_date           = as.character(start_date),
    plot_title          = paste("Age group:", age_labels[a]),
    age_color           = age_colors[age_labels[a]],
    forecast_start_date = forecast_start,
    observed_points     = train_series,
    true_points         = test_series
  )
})
fig_C <- .assemble_age_fig(C_plots, n_col, "Reported cases")
ggsave("figures/real_data/fig_age_C.pdf", fig_C, width = w_age, height = h_age)
print(fig_C)

# ------------------------------------------------------------
# 12. CRPS at each forecast horizon, by age group
# ------------------------------------------------------------
N_total          <- dim(X_marginal)[1]
T_obs            <- dim(X_marginal)[2]
n_particles_crps <- 2000L
set.seed(42)
particle_idx <- sample(seq_len(N_total), min(n_particles_crps, N_total))

forecast_dates   <- covid_data %>% filter(Date > cut_date & Date <= cut_date + horizon) %>% pull(Date)
H                <- length(forecast_dates)
t_forecast_start <- opts$T + 1L  # first forecast column in X_marginal / C_sim

crps_mat <- coverage_mat <- matrix(NA_real_, nrow = H, ncol = A)

for (a in seq_len(A)) {
  obs_series <- covid_data %>% filter(Date > cut_date & Date <= cut_date + horizon) %>% pull(a + 1)
  for (h in seq_len(H)) {
    t_idx <- t_forecast_start + h - 1L
    if (t_idx > T_obs) next
    sim_log <- log1p(C_sim[particle_idx, t_idx, a])
    obs_log <- log1p(obs_series[h])
    crps_mat[h, a] <- scoringRules::crps_sample(y = obs_log, dat = sim_log)
    q <- quantile(sim_log, c(0.025, 0.975))
    coverage_mat[h, a] <- obs_log >= q[1] && obs_log <= q[2]
  }
}

crps_df <- as.data.frame(crps_mat) %>%
  setNames(age_labels) %>%
  mutate(date = forecast_dates, days = as.numeric(forecast_dates - cut_date)) %>%
  tidyr::pivot_longer(cols = all_of(age_labels), names_to = "age_group", values_to = "crps") %>%
  mutate(age_group = factor(age_group, levels = age_labels))

coverage95 <- colMeans(coverage_mat, na.rm = TRUE) * 100

age_labels_crps <- crps_df %>%
  group_by(age_group) %>%
  summarise(mean_crps = mean(crps, na.rm = TRUE), .groups = "drop") %>%
  mutate(label = sprintf("%s (CRPS = %.2f,  Cov95 = %.1f%%)",
                         age_group, mean_crps,
                         coverage95[match(age_group, age_labels)])) %>%
  { setNames(.$label, .$age_group) }

fig_crps <- ggplot(crps_df, aes(x = days, y = crps, colour = age_group, group = age_group)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = age_colors, name = "Age group", labels = age_labels_crps) +
  scale_x_continuous(breaks = unique(crps_df$days)) +
  labs(x = "Forecast horizon (days)", y = "CRPS (log scale)") +
  theme_bw(base_size = 11) +
  theme(
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background = element_rect(colour = "black", fill = "white"),
    legend.margin = margin(4, 6, 4, 6),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold")
  )

ggsave("figures/real_data/fig_crps.pdf", fig_crps, width = 8, height = 4)
print(fig_crps)
# ------------------------------------------------------------
# 13. Effective reproduction number (Rt / R_{a,t}) via NGM
# ------------------------------------------------------------
rt_res <- compute_Rt_Rat(X_marginal[, , 2 * seq_len(A) - 1], M_4x4)

rat_plots <- lapply(seq_len(A), function(a) {
  plot_filter_one_age(
    theta_mat           = rt_res$Rat[, , a],
    init_date           = as.character(start_date),
    plot_title          = paste("Age group:", age_labels[a]),
    age_color           = age_colors[age_labels[a]],
    forecast_start_date = forecast_start,
    ymin = 0.5, ymax = 1.5
  ) + geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.7)
})
fig_rat <- .assemble_age_fig(
  rat_plots, n_col,
  bquote(bold("Age-specific reproduction contribution") ~ ~ bolditalic(R)[a * "," ~ t])
)
ggsave("figures/real_data/fig_age_Rat.pdf", fig_rat, width = w_age, height = h_age)
print(fig_rat)

# ------------------------------------------------------------
# 14. Population-level totals: Rt, total infections, total cases
# ------------------------------------------------------------
I_total_mat <- sort_then_sum(lapply(seq_len(A), function(a) X_marginal[, 1:T_obs, 2 * (a - 1) + 2]))
C_total_mat <- sort_then_sum(lapply(seq_len(A), function(a) C_sim[, 1:T_obs, a]))

theme_strip_x <- theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
                       axis.ticks.x = element_blank())
theme_panel <- theme(axis.text = element_text(size = 18), axis.title = element_text(size = 18),
                     plot.margin = margin(12, 5, 12, 5))

fig_Rt <- plot_filter_one(theta_mat = rt_res$Rt, init_date = as.character(start_date),
                          y_label = "Eff. reproduction number", forecast_start_date = forecast_start) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 1.0) + theme_panel + theme_strip_x

fig_I_total <- plot_filter_one(theta_mat = I_total_mat, init_date = as.character(start_date),
                               y_label = "New infections", forecast_start_date = forecast_start) +
  theme_panel + theme_strip_x

fig_C_total <- plot_filter_one(theta_mat = C_total_mat, init_date = as.character(start_date),
                               y_label = "Reported cases", forecast_start_date = forecast_start,
                               observed_points = train_data$Total_Daily_Cases,
                               true_points     = test_data$Total_Daily_Cases) +
  coord_cartesian(ylim = c(0, 14000)) + theme_panel

fig_combined <- fig_Rt / fig_I_total / fig_C_total + plot_layout(heights = c(1, 1, 1))
ggsave("figures/real_data/fig_combined.pdf", fig_combined, width = 12, height = 13, limitsize = FALSE)
print(fig_combined)

