# ============================================================
# run_1age_covid.R 
#
# Fits the single-population (non-age-stratified) EpiSSM to the
# total daily case series from Covid_Data_Ireland.csv via SMC^2, and
# plots filtering-distribution estimates only (no PMMH comparison -
# that has been dropped per project decision).
#

# Run from the repository root:  scripts/run_1age_covid.R
# ============================================================

library(dplyr)
library(ggplot2)
library(cowplot)
library(patchwork)

source("src/smc2.R")                # sources utils/resampling/particle_filter/priors/mcmc_kernels
source("src/epi_ssm_1age.R")
source("src/sim_data_1age.R")       # generation-time / report-delay kernels
source("src/diagnostics.R")
source("src/posterior_marginal.R")

dir.create("figures/real_data", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# 1. Load and window the real case data
# ------------------------------------------------------------
covid <- read.csv("data/real/Covid_Data_Ireland.csv")
covid$Date <- as.Date(covid$Date)

start_date <- as.Date("2022-03-01")
end_date   <- as.Date("2022-07-17")

covid <- covid %>% filter(Date >= start_date & Date <= end_date)
covid_sub <- covid %>% dplyr::select(Date, Total_Daily_Cases)


# ------------------------------------------------------------
#  Delay kernels
# ------------------------------------------------------------
gen_pmf        <- discretize_dist(max = 14, mean = 3.3, sd = 2.1, dist = "gamma")
incub_pmf      <- discretize_dist(max = 14, mean = 3.0, sd = 2.3, dist = "gamma")
report_pmf     <- discretize_dist(max = 10, mean = 2.0, sd = 1.0, dist = "gamma")
case_delay_pmf <- convolve(incub_pmf, rev(report_pmf), type = "open")
case_delay_pmf <- pmax(case_delay_pmf, 0)
case_delay_pmf <- case_delay_pmf / sum(case_delay_pmf)

# ------------------------------------------------------------
# 2. SMC^2 options
# ------------------------------------------------------------
opts <- list(
  Data = covid_sub,
  T    = nrow(covid_sub),
  
  N_x_init         = 300,
  L                = 30,
  resample_method  = "systematic",
  resample_scope   = "block",
  state_dim        = 2,
  obs_dim          = 1,
  forecastingHorizon   = 0,
  PredictedObservation = FALSE,
  
  GenTime        = gen_pmf,
  InfReportDelay = case_delay_pmf,
  
  pR0 = list(function(N) rnorm(N, 1.5, 0.2)),
  pI0 = list(function(N) sample(2500:5000, N, replace = TRUE)),
  
  A          = 1,
  ContMatrix = diag(1),
  
  week_effect = TRUE,
  
  N_theta       = 800,
  N_max         = 1000,
  ESS_threshold = 0.4,
  
  # theta[1]   : sigma  ~ Gamma(shape=1, scale=0.02)
  # theta[2:7] : omega_1..6 ~ LogNormal(0, 0.2), bounds [0, 2]
  # theta[8]   : kappa  ~ Exponential(rate=0.1)
  paramNames = c(
    expression(sigma),
    expression(omega[1]), expression(omega[2]), expression(omega[3]),
    expression(omega[4]), expression(omega[5]), expression(omega[6]),
    expression(kappa)
  ),
  paramPriors = list(
    list(dist = "gamma", shape = 1, scale = 0.02),
    list(dist = "lnorm", meanlog = 0, sdlog = 0.2),
    list(dist = "lnorm", meanlog = 0, sdlog = 0.2),
    list(dist = "lnorm", meanlog = 0, sdlog = 0.2),
    list(dist = "lnorm", meanlog = 0, sdlog = 0.2),
    list(dist = "lnorm", meanlog = 0, sdlog = 0.2),
    list(dist = "lnorm", meanlog = 0, sdlog = 0.2),
    list(dist = "exp",   rate = 0.1)
  ),
  lower_bounds = c(0,   0, 0, 0, 0, 0, 0,   0),
  upper_bounds = c(0.5, 2, 2, 2, 2, 2, 2, 100),
  
  param_blocks = list(c(1), c(2:7), c(8)),
  
  n_cores  = 20,
  parallel = TRUE
)

# ------------------------------------------------------------
# 3. Run SMC^2
# ------------------------------------------------------------
result <- SMC2(SSM = EpiSSM_1age, opts = opts)
saveRDS(result, file = "data/real/smc2_1age_covid_result.rds")
# result <-  readRDS("data/real/smc2_1age_covid_result.rds") # load saved result
# ------------------------------------------------------------
# 4. Diagnostics + parameter filtering estimates
# ------------------------------------------------------------
fig_diagnostic <- plot_SMC2_diagnostics(result)
ggsave("figures/real_data/fig_1age_diagnostic.pdf", fig_diagnostic, width = 24, height = 14)
print(fig_diagnostic)


T_last <- dim(result$theta)[2]

# DoW transformation theta[2:7]
dow <- result$theta[, , 2:7]
den <- 1 + apply(dow, c(1,2), sum)
omega <- abind::abind(sweep(7*dow, c(1,2), den, "/"), 7/den, along=3)

# ---- Parameter plots with transformed DoW effects ----
theta_plot <- result
theta_plot$theta[, , 2:7] <- omega[, , 1:6]

fig_params <- .grid_shared_legend(
  plot_filter_all(theta_plot, opts), ncol = 2
)

ggsave("figures/real_data/fig_1age_params.pdf", fig_params,
       width = 24, height = 6*ceiling(length(fig_params$grobs)/2))
print(fig_params)

# ------------------------------------------------------------
# 5. Day-of-week reporting effect (final-time posterior)
# ------------------------------------------------------------
omega_df <- as.data.frame(omega[, T_last, ]) |>
  setNames(c("Mon","Tue","Wed","Thu","Fri","Sat","Sun")) |>
  tidyr::pivot_longer(everything(), names_to="Day", values_to="omega") |>
  mutate(Day = factor(Day, levels=c("Mon","Tue","Wed","Thu","Fri","Sat","Sun")))

fig_week_effect <- ggplot(omega_df, aes(Day, omega)) +
  geom_boxplot(fill="orange3", outlier.shape=NA,
               staplewidth=.7, width=.4) +
  geom_hline(yintercept=1, linetype="dashed") +
  labs(x="", y="Reporting effect") +
  theme_bw()

ggsave("figures/real_data/fig_1age_week_effect.pdf",
       fig_week_effect, width=8, height=4)

print(fig_week_effect)


# ------------------------------------------------------------
# 6. Posterior marginal states + simulated observations
# ------------------------------------------------------------
source("src/marginal_state_1age.R")
opts$forecastingHorizon <- 14
marg       <- get_marginal_states_1age(result, opts, N_s = 100)
X_marginal <- marg$X_marginal
theta_sub  <- marg$theta_sub
obs_sim    <- simulate_observations_1age(X_marginal, theta_sub, opts)
C_sim      <- obs_sim$C_sim


#============================================
# we can also use final theta to run the BPF
#============================================
# 
# theta_sub <- result$theta[sample(seq_len(opts$N_theta), Ns), opts$T, ]
# 
# marg       <- PosteriorMarginal(theta_sub, EpiSSM_1age, opts)
# X_marginal <- marg$X        # (Ns * Nx) x (T + h) x state_dim
# C_sim      <- marg$X_obs    # (Ns * Nx) x (T + h) x 1

beta_idx <- 1L  # R_t (state column 1)
I_idx    <- 2L  # new infections (state column 2)

w <- 12; h <- 5

fig_Rt <- .grid_shared_ylabel(
  list(plot_filter_one(theta_mat = X_marginal[, 1:opts$T, beta_idx]) +
         geom_hline(yintercept = 1, linetype = "dashed")),
  ncol = 1, y_label = "Eff. reproduction number"
)
ggsave("figures/real_data/fig_1age_Rt.pdf", fig_Rt, width = w, height = h)
print(fig_Rt)

fig_I <- .grid_shared_ylabel(
  list(plot_filter_one(theta_mat = X_marginal[, 1:opts$T, I_idx])),
  ncol = 1, y_label = "New infections"
)
ggsave("figures/real_data/fig_1age_I.pdf", fig_I, width = w, height = h)
print(fig_I)

fig_C <- .grid_shared_ylabel(
  list(plot_filter_one(theta_mat = matrix(C_sim[, 1:opts$T, 1], nrow = dim(C_sim)[1]),
                       true_points = covid$Total_Daily_Cases)),
  ncol = 1, y_label = "Reported cases"
)
ggsave("figures/real_data/fig_1age_reported_cases.pdf", fig_C, width = w, height = h)
print(fig_C)
