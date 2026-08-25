# ===========================================================================
# Rt Estimation using EpiNow2
# ===========================================================================
# This script estimates the time-varying reproduction number (Rt) from daily
# case counts using the EpiNow2 package.
#
# Install required packages (run once):
#   install.packages(c("readr", "dplyr", "lubridate"))
#   install.packages("EpiNow2")
#
# Data input:
#   A CSV file with at least two columns: "Date" and "Total_Daily_Cases"
# ===========================================================================

library(readr)
library(dplyr)
library(lubridate)
library(EpiNow2)

# ===========================================================================
# Load and prepare data
# ===========================================================================
start_date <- as.Date("2022-03-01")
end_date   <- as.Date("2022-07-17")

case_data <- read.csv("data/real/Covid_Data_Ireland.csv") %>%
  mutate(
    date    = as.Date(Date),
    confirm = as.integer(Total_Daily_Cases)
  ) %>%
  filter(date >= start_date & date <= end_date) %>%
  select(date, confirm) %>%
  arrange(date)

# View(case_data)

# ===========================================================================
# Helper: convert mean & sd to Gamma shape & rate
# ===========================================================================
gamma_shape_rate <- function(mean, sd) {
  stopifnot(mean > 0, sd > 0)
  shape <- (mean / sd)^2
  rate  <- mean / (sd^2)
  list(shape = shape, rate = rate)
}

# ===========================================================================
# Transform mean/sd into shape/rate for EpiNow2::Gamma
# ===========================================================================
reporting_delay_pars   <- gamma_shape_rate(mean = 2,   sd = 1)
generation_time_pars   <- gamma_shape_rate(mean = 3.3, sd = 2.1)
incubation_period_pars <- gamma_shape_rate(mean = 3,   sd = 2.3)

# ===========================================================================
# Build distributions
# ===========================================================================
reporting_delay <- Gamma(
  shape = reporting_delay_pars$shape,
  rate  = reporting_delay_pars$rate,
  max   = 10
)

generation_time <- Gamma(
  shape = generation_time_pars$shape,
  rate  = generation_time_pars$rate,
  max   = 14
)

incubation_period <- Gamma(
  shape = incubation_period_pars$shape,
  rate  = incubation_period_pars$rate,
  max   = 14
)

# ===========================================================================
# Run EpiNow2
# ===========================================================================
estimates <- epinow(
  data            = case_data,
  generation_time = gt_opts(generation_time),
  delays          = delay_opts(incubation_period + reporting_delay),
  rt              = rt_opts(prior = LogNormal(mean = 1.5, sd = 0.2)),
  stan            = stan_opts(cores = 4),
  verbose         = interactive()
)

plot(estimates)

# ===========================================================================
# Extract and save Rt trajectory
# ===========================================================================

# Extract summarised parameter estimates
rt_trajectory <- summary(estimates, type = "parameters")
# Keep only the Rt estimates
rt_trajectory <- subset(rt_trajectory, variable == "R")
# Keep the date and uncertainty estimates
rt_trajectory <- rt_trajectory[, c("date", "median", "lower_90", "upper_90")]

saveRDS(rt_trajectory, file = "data/real/rt_epinow2.rds")
