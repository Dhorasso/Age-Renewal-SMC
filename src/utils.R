# ============================================================
# utils.R
#
# Small, dependency-light helpers shared across the package:
#   - %||%           null-coalescing operator
#   - logSumExp      numerically stable log-sum-exp
#   - format_time    seconds -> "H:MM:SS" string
#   - print_progress console progress bar for SMC2()
#   - par_lapply     lapply / mclapply switch

# Prior specification, sampling, and log-density evaluation for the
# theta-particles in SMC2(). Extracted from the original SMC2.R.
#
# A prior for parameter j (`opts$paramPriors[[j]]`) can be either:
#   - a list `list(dist = "norm", mean = 0, sd = 1)` -> uses
#     `dnorm`/`rnorm` (or any `d<dist>`/`r<dist>` pair on the search
#     path, including `dhalfnorm`/`rhalfnorm` defined below), or
#   - a raw function `function(x) ...` returning a log-density
#     (used only by `computeLogPrior`, not by `priorSampler`).
#
# Sourced by every other file in src/. Keep this file dependency-free
# (base R only) so it can be loaded first with no side effects.
# ============================================================

#' Null-coalescing operator
#'
#' @param x value to test
#' @param y fallback used when `x` is NULL
#' @return `x` if not NULL, otherwise `y`
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Numerically stable log-sum-exp
#'
#' @param lx numeric vector of log-values
#' @return scalar, log(sum(exp(lx)))
logSumExp <- function(lx) {
  m <- max(lx)
  m + log(sum(exp(lx - m)))
}

#' Format a duration in seconds as H:MM:SS
#'
#' @param seconds numeric, duration in seconds
#' @return character string
format_time <- function(seconds) {
  if (!is.finite(seconds)) return("--:--:--")
  h <- floor(seconds / 3600)
  m <- floor((seconds %% 3600) / 60)
  s <- floor(seconds %% 60)
  sprintf("%d:%02d:%02d", h, m, s)
}

#' Print a single-line SMC2 progress bar
#'
#' @param t current time step
#' @param T total number of time steps
#' @param start_time `Sys.time()` value taken at the start of the run
#' @param ESS optional effective sample size to display
#' @param Nx optional number of PF particles to display
#' @param mode label for the current move-step kernel (e.g. "HWG", "DA-HWG")
print_progress <- function(t, T, start_time, ESS = NULL, Nx = NULL,
                            mode = "Standard") {
  pct       <- t / T
  elapsed   <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  bar_width <- 30L
  filled    <- floor(pct * bar_width)

  bar <- paste0(
    "\rSMC^2 progress: [",
    paste(rep("=", filled),             collapse = ""),
    ">",
    paste(rep(" ", bar_width - filled), collapse = ""),
    "] ",
    sprintf("%3.0f%%", pct * 100),
    " | Time: ", format_time(elapsed),
    " | Kernel: ", mode
  )
  if (!is.null(ESS)) bar <- paste0(bar, sprintf(" | ESS: %4.0f", ESS))
  if (!is.null(Nx))  bar <- paste0(bar, sprintf(" | Nx: %.1f",   Nx))

  cat(bar)
  if (t == T) cat("\n")
  flush.console()
}

#' Dispatch to lapply() or parallel::mclapply()
#'
#' @param X list/vector to iterate over
#' @param FUN function to apply
#' @param n_cores number of cores to use when `use_parallel = TRUE`
#' @param use_parallel logical switch
#' @return list of results, same as lapply()
par_lapply <- function(X, FUN, n_cores, use_parallel) {
  if (use_parallel) parallel::mclapply(X, FUN, mc.cores = n_cores)
  else               lapply(X, FUN)
}




## helper functions for prior sampling
# ------------------------------------------------------------------
# Half-normal distribution (scale sigma; mean = sigma * sqrt(2/pi))
# ------------------------------------------------------------------
dhalfnorm <- function(x, sigma = 1, log = FALSE) {
  lp <- ifelse(x < 0, -Inf, dnorm(x, mean = 0, sd = sigma, log = TRUE) + log(2))
  if (log) lp else exp(lp)
}

rhalfnorm <- function(n, sigma = 1) abs(rnorm(n, mean = 0, sd = sigma))

# ------------------------------------------------------------------
# Log-prior density, given a parameter vector and prior spec list
# ------------------------------------------------------------------
#' @param theta numeric vector of parameter values
#' @param priors list of per-parameter prior specs (see file header)
#' @return scalar log-prior density (sum across parameters)
computeLogPrior <- function(theta, priors) {
  nParams     <- length(theta)
  logPriorSum <- 0
  
  for (i in seq_len(nParams)) {
    x         <- theta[i]
    priorSpec <- priors[[i]]
    
    if (is.list(priorSpec) && !is.null(priorSpec$dist)) {
      distFun   <- get(paste0("d", priorSpec$dist))
      args      <- priorSpec
      args$x    <- x
      args$log  <- TRUE
      args$dist <- NULL
      logp      <- do.call(distFun, args)
      if (is.na(logp)) logp <- -Inf
      
    } else if (is.function(priorSpec)) {
      logp <- priorSpec(x)
      if (is.na(logp)) logp <- -Inf
      
    } else {
      stop("Invalid prior specification")
    }
    
    logPriorSum <- logPriorSum + logp
  }
  
  logPriorSum
}

# ------------------------------------------------------------------
# Sample N draws from the joint prior, respecting box constraints
# ------------------------------------------------------------------
#' @param n number of draws
#' @param opts options list; must contain `paramPriors`, may contain
#'   `lower_bounds` / `upper_bounds`
#' @return matrix [n, d] of prior draws
priorSampler <- function(n, opts) {
  d   <- length(opts$paramPriors)
  out <- matrix(NA, n, d)
  
  lb <- opts$lower_bounds %||% rep(-Inf, d)
  ub <- opts$upper_bounds %||% rep( Inf, d)
  
  for (j in seq_len(d)) {
    p <- opts$paramPriors[[j]]
    
    sampler <- if (is.list(p)) {
      f    <- get(paste0("r", p$dist))
      args <- p
      args$dist <- NULL
      function(n) do.call(f, c(list(n = n), args))
    } else {
      p  # assume p is already an r-function
    }
    
    samples  <- sampler(n)
    in_range <- samples >= lb[j] & samples <= ub[j]
    
    # rejection-resample until all n draws are in bounds
    max_iter <- 1000L
    iter     <- 0L
    while (any(!in_range) && iter < max_iter) {
      iter               <- iter + 1L
      n_bad              <- sum(!in_range)
      samples[!in_range] <- sampler(n_bad)
      in_range           <- samples >= lb[j] & samples <= ub[j]
    }
    
    if (any(!in_range)) {
      warning(sprintf(
        "Parameter %d: %d sample(s) still out of [%g, %g] after %d iterations.",
        j, sum(!in_range), lb[j], ub[j], max_iter
      ))
    }
    
    out[, j] <- samples
  }
  out
}

# ------------------------------------------------------------------
# Log-prior density with box-constraint short-circuit
# ------------------------------------------------------------------
#' @param theta numeric parameter vector
#' @param opts options list; must contain `paramPriors`, may contain
#'   `lower_bounds` / `upper_bounds`
#' @return scalar log-prior density, or -Inf if outside bounds
logPriorDensity <- function(theta, opts) {
  if (!is.null(opts$lower_bounds) && !is.null(opts$upper_bounds)) {
    if (any(theta < opts$lower_bounds) || any(theta > opts$upper_bounds))
      return(-Inf)
  }
  computeLogPrior(theta, opts$paramPriors)
}
