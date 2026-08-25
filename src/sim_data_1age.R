# ============================================================
# sim_data_1age.R
#
# NOTE: this replaces the missing `src/sim_data1age.R` referenced by
# the original `test_1age_smc2.R`. That script's active code path
# fits the *real* Covid_Data_Ireland.csv series (not simulated data) -
# it only ever used the missing file for two delay kernels: `w`
# (generation time) and `inf_report_pmf` (infection-to-report delay).
# Everything else the comment in the original file promised
# ("generates Y, Cmat, w, d") was dead weight: `Y`, `Cmat`, `d` were
# never referenced downstream.
#
# This file supplies just the two kernels that are actually used, via
# the discretized distributions already in delay_distributions.R.
# The specific mean/sd choices below are illustrative defaults from
# the equivalent block in Ireland_analysis.R - replace them if your
# 1-age study used different assumptions.
# ============================================================

source("src/delay_distributions.R")

#' Generation-time and infection-to-report delay kernels for the 1-age model
#'
#' @param gen_mean,gen_sd mean/sd (days) of the generation-time distribution
#' @param gen_max maximum generation-time lag to retain
#' @param incub_mean,incub_sd mean/sd (days) of the incubation period
#' @param report_mean,report_sd mean/sd (days) of onset-to-report delay
#' @param delay_max maximum lag to retain for the combined report delay
#' @return list with `w` (generation-time PMF) and `inf_report_pmf`
#'   (infection-to-report delay PMF, the convolution of incubation and
#'   onset-to-report delay)
build_1age_kernels <- function(gen_mean = 3.3, gen_sd = 2.1, gen_max = 14,
                                incub_mean = 3.0, incub_sd = 2.3,
                                report_mean = 2, report_sd = 1,
                                delay_max = 14) {

  w <- discretize_dist(max = gen_max, mean = gen_mean, sd = gen_sd, dist = "gamma")

  incub_pmf  <- discretize_dist(max = delay_max, mean = incub_mean,  sd = incub_sd,  dist = "gamma")
  report_pmf <- discretize_dist(max = delay_max, mean = report_mean, sd = report_sd, dist = "gamma")

  inf_report_pmf <- convolve(incub_pmf, rev(report_pmf), type = "open")
  inf_report_pmf <- pmax(inf_report_pmf, 0)
  inf_report_pmf <- inf_report_pmf / sum(inf_report_pmf)

  list(w = w, inf_report_pmf = inf_report_pmf)
}
