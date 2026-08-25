# ============================================================
# simulated_data_cache.R
#
# Generic "simulate once, then cache to CSV" wrapper. Re-running a
# stochastic data-generating process every time a script is sourced is
# slow and (without a fixed seed) non-reproducible across runs. This
# wrapper simulates once, writes the result to data/simulated/ as CSV,
# and on every subsequent call just reads the CSV back — no need to
# re-run the simulation, and the cached data is easy to inspect,
# version, or share independently of R.
#
# Works for any simulator that returns a *named list of data frames*
# (e.g. `list(Y_df = ..., beta_true = ..., I_true = ...)`); each
# element is written to / read from its own CSV file
# `data/simulated/<label>_<name>.csv`.
# ============================================================

#' Simulate (once) or load cached simulated data from CSV
#'
#' @param label short identifier for this dataset, used as a filename
#'   prefix (e.g. "sim1age", "sim_multiage")
#' @param simulate_fn zero-argument function that runs the simulation
#'   and returns a named list of data frames / matrices
#' @param dir directory to cache CSVs in (default "data/simulated")
#' @param overwrite if TRUE, re-run `simulate_fn()` and overwrite any
#'   existing cache (default FALSE)
#' @return the named list returned by `simulate_fn()` (either freshly
#'   simulated, or reconstructed from the CSV cache)
#'
#' @details Each list element is coerced with `as.data.frame()` before
#'   writing, and read back with `read.csv(..., stringsAsFactors = FALSE)`.
#'   A single-column element is returned as a plain vector; anything
#'   with >1 column, or with a "Date" column, is returned as a data
#'   frame (with "Date" parsed back to `Date`, if present) rather than
#'   silently coerced to a matrix, since coercion rules depend on the
#'   consumer downstream. If a given element started life as a plain
#'   numeric matrix, wrap it back with `as.matrix()` at the call site.
#'
#' @examples
#' \dontrun{
#' sim <- load_or_simulate("sim1age", function() {
#'   list(Y_df = simulate_one_age(...))
#' })
#' }
load_or_simulate <- function(label, simulate_fn,
                              dir = "data/simulated", overwrite = FALSE) {

  dir.create(dir, showWarnings = FALSE, recursive = TRUE)

  manifest_path <- file.path(dir, paste0(label, "_manifest.txt"))

  if (!overwrite && file.exists(manifest_path)) {
    names_cached <- readLines(manifest_path)
    out <- setNames(vector("list", length(names_cached)), names_cached)

    for (nm in names_cached) {
      path <- file.path(dir, paste0(label, "_", nm, ".csv"))
      df   <- read.csv(path, stringsAsFactors = FALSE)
      if ("Date" %in% names(df)) df$Date <- as.Date(df$Date)
      out[[nm]] <- if (ncol(df) == 1L) df[[1]] else df
    }

    message(sprintf(
      "[load_or_simulate] Loaded cached '%s' from %s (%d object(s)).",
      label, dir, length(names_cached)
    ))
    return(out)
  }

  message(sprintf("[load_or_simulate] Simulating '%s' ...", label))
  sim <- simulate_fn()

  for (nm in names(sim)) {
    path <- file.path(dir, paste0(label, "_", nm, ".csv"))
    write.csv(as.data.frame(sim[[nm]]), path, row.names = FALSE)
  }
  writeLines(names(sim), manifest_path)

  message(sprintf(
    "[load_or_simulate] Cached '%s' to %s (%d object(s)).",
    label, dir, length(sim)
  ))
  sim
}
