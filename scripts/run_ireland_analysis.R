# ============================================================
# run_ireland_analysis.R
#
# Fits the age-stratified EpiSSM to real Irish COVID-19 case data via
# SMC^2, produces a 14-day-ahead forecast, and scores it with CRPS.
#
# Run from the repository root:  scripts/run_ireland_analysis.R
# ============================================================

library(tidyverse)
library(lubridate)
library(scoringRules)

source("src/smc2.R")
source("src/epi_ssm_multiage.R")
source("src/delay_distributions.R")
source("src/contact_matrix_utils.R")
source("src/diagnostics.R")
source("src/epi_diagnostics.R")
source("src/marginal_state.R")
source("src/posterior_marginal.R")


dir.create("figures/real_data", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# 1. Load and window the real case data
# ------------------------------------------------------------
covid_data      <- read.csv("data/real/Covid_Data_Ireland.csv")
covid_data$Date <- as.Date(covid_data$Date)

start_date <- as.Date("2022-03-01")
end_date   <- as.Date("2022-07-17")
horizon    <- 14L
cut_date   <- end_date

covid <- covid_data %>% filter(Date >= start_date & Date <= end_date)

train_data <- covid_data %>% filter(Date >= start_date & Date <= cut_date)
test_data  <- covid_data %>% filter(Date > cut_date & Date <= cut_date + horizon)

forecast_start <- as.numeric(cut_date)  # numeric date for geom_vline xintercept

# ------------------------------------------------------------
# 2. Delay kernels
# ------------------------------------------------------------
gen_pmf        <- discretize_dist(max = 14, mean = 3.3, sd = 2.1, dist = "gamma")
incub_pmf      <- discretize_dist(max = 14, mean = 3.0, sd = 2.3, dist = "gamma")
report_pmf     <- discretize_dist(max = 10, mean = 2.0, sd = 1.0, dist = "gamma")
case_delay_pmf <- convolve(incub_pmf, rev(report_pmf), type = "open")
case_delay_pmf <- pmax(case_delay_pmf, 0)
case_delay_pmf <- case_delay_pmf / sum(case_delay_pmf)

# ------------------------------------------------------------
# 3. Age-group contact matrix (real Ireland contact data)
# ------------------------------------------------------------
groups <- list(G1 = 1:25, G2 = 26:45, G3 = 46:65, G4 = 66:85)  # 0-24, 25-44, 45-64, 65+
A      <- length(groups)

# Population by single year of age, 0–99
N_pop_full <- read.csv("data/contact_matrices/PEA11.20260324T160347.csv")
N_pop_full <-N_pop_full$VALUE
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
# ------------------------------------------------------------
# 4. Day-of-week weights (precomputed from observed data)
# ------------------------------------------------------------
omega_dow_mat <- compute_dow_weights(covid, A, max_weeks = 16)

# ------------------------------------------------------------
# 5. SMC^2 options
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
  
  # sigma[1..4] : RW volatility for log-transmissibility per age group
  # kappa       : negative-binomial overdispersion for case counts
  # theta layout: [sigma_1..sigma_4, rho_1..rho_4, kappa]
  
  ##  we can fix value of rho_a to the mean for fast infrence
  ##  by defining CAR <-c(rho_1, rho_2, rho_3, rho_4)
  paramNames = c(
    expression(sigma[1]), expression(sigma[2]),
    expression(sigma[3]), expression(sigma[4]),
    expression(rho[1]),   expression(rho[2]),
    expression(rho[3]),   expression(rho[4]),
    expression(kappa)
  ),

  paramPriors = list(
    list(dist = "lnorm", meanlog = log(0.02), sdlog = 0.5),
    list(dist = "lnorm", meanlog = log(0.02), sdlog = 0.5),
    list(dist = "lnorm", meanlog = log(0.02), sdlog = 0.5),
    list(dist = "lnorm", meanlog = log(0.02), sdlog = 0.5),
    list(dist = "beta",  shape1 = 42, shape2 = 193),  # rho_1 (0-24):   mean 0.18, sd 0.025
    list(dist = "beta",  shape1 = 60, shape2 = 213),  # rho_2 (25-44):  mean 0.22, sd 0.025
    list(dist = "beta",  shape1 = 77, shape2 = 164),  # rho_3 (45-64):  mean 0.32, sd 0.030
    list(dist = "beta",  shape1 = 83, shape2 = 115),  # rho_4 (65+):    mean 0.42, sd 0.035
    list(dist = "exp",   rate = 0.1)
  ),
  
  param_blocks <- list(c(1,5), c(2,6), c(3,7), c(4,8), c(9)),
  lower_bounds = c(0, 0, 0, 0, 0, 0, 0, 0,   0),
  upper_bounds = c(1, 1, 1, 1, 1, 1, 1, 1, 200),
  N_x_init         = 300,
  N_theta          = 800,
  ESS_threshold    = 0.4,
  k_neighbors      = 3,
  sigma2_threshold = 1.5,
  N_max            = 1000,
  n_cores          = 16,
  parallel         = TRUE
)

# ------------------------------------------------------------
# 6. Run SMC^2
# ------------------------------------------------------------
result <- SMC2(SSM = EpiSSM, opts = opts)
# saveRDS(result, file = "data/real/smc2_ireland_result.rds") # Save results
# result <-  readRDS("data/real/smc2_ireland_result_rho.rds") # load saved result
# ------------------------------------------------------------
# 7. Algorithm diagnostics + parameter filtering estimates
# ------------------------------------------------------------
fig_diagnostic <- plot_SMC2_diagnostics(result)
ggsave("figures/real_data/fig_diagnostic.pdf", fig_diagnostic, width = 24, height = 14)
print(fig_diagnostic)

param_plots <- lapply(seq_along(opts$paramNames), function(j) {
  plot_filter_one(theta_mat = result$theta[, , j], y_label = opts$paramNames[[j]])
})
fig_param <- .grid_shared_legend(param_plots, ncol = 2)
ggsave("figures/real_data/fig_param.pdf", fig_param,
       width = 24, height = 6 * ceiling(length(param_plots) / 2))
print(fig_param)

# ------------------------------------------------------------
# 8. Posterior marginal states, with 14-day forecast
# ------------------------------------------------------------

# ------------------------------------------------------------
Ns        <- 100
theta_sub <- result$theta[sample(seq_len(opts$N_theta), Ns), opts$T, ]

# # Overwrite rho_1..rho_A columns with fresh prior draws
# for (a in seq_len(opts$A)) {
#   idx    <- .rho_idx(a, opts)
#   prior  <- opts$paramPriors[[idx]]
#   stopifnot(prior$dist == "beta")  # guard against index drift if theta layout changes
#   theta_sub[, idx] <- rbeta(Ns, shape1 = prior$shape1, shape2 = prior$shape2)
# }



#source("src/marginal_state.R")
opts$forecastingHorizon <- 14
marg       <- get_marginal_states(result, opts, N_s = 100)
X_marginal <- marg$X_marginal
theta_sub  <- marg$theta_sub
obs_sim    <- simulate_observations(X_marginal, theta_sub, opts)
C_sim      <- obs_sim$C_sim


#============================================
# we can also use final theta to run the BPF
#============================================

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
fig_I <- .assemble_age_fig(I_plots, n_col, y_label = "New infections", x_label = "Date")
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
fig_C <- .assemble_age_fig(C_plots, n_col, y_label = "Reported cases", x_label = "Date")
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



















