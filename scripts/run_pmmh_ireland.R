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
