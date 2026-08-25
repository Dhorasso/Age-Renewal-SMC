# =============================================================================
# diagnostics.R 
# =============================================================================
# Purpose : Reusable ggplot2 helpers for SMC^2 / EpiSSM visualisation.
#           All functions return ggplot objects that can be composed or saved
#           by the caller. No side-effects (no ggsave calls here).
#
#
# Public API
# ----------
#   .smc2_theme()               - shared ggplot2 theme
#   plot_filter_one()           - credible-interval ribbon for a scalar state
#   plot_filter_one_age()       - same, coloured by age group
#   plot_filter_all()           - list of plots for every parameter
#   plot_prior_posterior_one()  - prior vs posterior density overlay
#   plot_prior_posterior_all()  - list of prior/posterior plots per parameter
#   plot_SMC2_diagnostics()     - 4-panel ESS / Nx / acceptance / log-evidence
#   .assemble_age_fig()         - lay out age-group panels with shared axis labels
#   .grid_shared_legend()       - grid of panels + one shared (enlarged) legend
#   .grid_shared_ylabel()       - grid of panels + shared y-axis label
#
# Internal helpers (prefixed with `.`)
#   .make_x_axis()               - build date or integer x-axis spec
#   .small_generic()             - reduce axis-title font size
# =============================================================================

library(ggplot2)
library(ggdist)
library(cowplot)
library(patchwork)

source("src/utils.R")

# =============================================================================
# SHARED THEME
# =============================================================================

#' Publication-ready ggplot2 theme for SMC^2 outputs.
.smc2_theme <- function() {
  theme_bw() +
    theme(
      axis.text    = element_text(size = 20),
      axis.title   = element_text(size = 18, face = "bold"),
      text         = element_text(size = 16),
      legend.text  = element_text(size = 16),
      legend.title = element_text(size = 16, face = "bold"),
      plot.title   = element_text(size = 18, face = "bold", hjust = 0.5)
    )
}

#' Reduce axis-title font size to plain 16 pt (used in diagnostic panels).
#' @param which "x", "y", or "both" (default)
.small_generic <- function(which = "both") {
  small <- element_text(size = 16, face = "plain")
  switch(which,
         x    = theme(axis.title.x = small),
         y    = theme(axis.title.y = small),
         both = theme(axis.title.x = small, axis.title.y = small))
}

# =============================================================================
# X-AXIS HELPER
# =============================================================================

#' Build x-axis specification (date or integer).
#'
#' @param T number of time points
#' @param init_date "YYYY-MM-DD" string, or NULL for integer axis
#' @return list with `x_vals`, `scale` (ggplot2 scale or NULL), `x_label`
.make_x_axis <- function(T, init_date = NULL) {
  if (!is.null(init_date)) {
    init_date <- as.Date(init_date)
    list(
      x_vals  = seq(init_date, by = "day", length.out = T),
      scale   = scale_x_date(date_labels = "%b", date_breaks = "1 month"),
      x_label = "Date"
    )
  } else {
    list(x_vals = seq_len(T), scale = NULL, x_label = "Time")
  }
}

# =============================================================================
# FILTERING DISTRIBUTION - SINGLE SCALAR STATE (population-level)
# =============================================================================

#' Plot credible-interval ribbons for one scalar latent state over time.
#'
#' @param theta_mat N x T matrix (rows = particles, cols = time)
#' @param y_label y-axis label (string or expression)
#' @param plot_title optional plot title
#' @param init_date "YYYY-MM-DD" string; if given, x-axis shows dates
#' @param forecast_start_date numeric date (or integer) for a vertical dashed line
#' @param observed_points training observations, length T, plotted as solid red points
#' @param true_points test/forecast observations, plotted as open circles at
#'   the tail of the time axis (last `length(true_points)` time points)
#' @return a ggplot object
plot_filter_one <- function(theta_mat,
                             y_label             = NULL,
                             plot_title          = NULL,
                             init_date           = NULL,
                             forecast_start_date = NULL,
                             observed_points     = NULL,
                             true_points         = NULL) {

  N  <- nrow(theta_mat)
  T  <- ncol(theta_mat)
  ax <- .make_x_axis(T, init_date)

  df <- data.frame(Time = rep(ax$x_vals, each = N), value = as.vector(theta_mat))

  p <- ggplot(df, aes(Time, value)) +
    stat_lineribbon(
      aes(fill = after_stat(level)),
      .width    = c(0.50, 0.75, 0.90, 0.95),
      color     = "black",
      linewidth = 1.2
    ) +
    scale_fill_brewer(palette = "Blues") +
    labs(x = ax$x_label, y = y_label, title = plot_title) +
    .smc2_theme()

  if (!is.null(ax$scale)) p <- p + ax$scale

  if (!is.null(forecast_start_date)) {
    p <- p + geom_vline(xintercept = forecast_start_date, linetype = "dashed", linewidth = 1.1)
  }

  if (!is.null(observed_points)) {
    ax_tr <- .make_x_axis(length(observed_points), init_date)
    p <- p + geom_point(
      data = data.frame(Time = ax_tr$x_vals, value = observed_points),
      aes(Time, value), color = "red", size = 1.6, alpha = 0.7, inherit.aes = FALSE
    )
  }

  if (!is.null(true_points)) {
    x_pts <- .anchor_points(true_points, forecast_start_date, init_date, ax)
    p <- p + geom_point(
      data = data.frame(Time = x_pts, value = true_points),
      aes(Time, value), shape = 21, fill = "white", color = "black",
      size = 2, stroke = 0.8, inherit.aes = FALSE
    )
  }

  p
}

# =============================================================================
# FILTERING DISTRIBUTION - SINGLE SCALAR STATE (age-stratified)
# =============================================================================

#' Anchor a set of forecast/test points on the time axis.
#'
#' If `forecast_start_date` is supplied, points are placed starting the
#' day after it (matching how `forecast_start_date` is itself used for
#' the vertical cutoff line: `cut_date` in the driver scripts). This
#' replaces a previous version that referenced an undeclared global
#' `cut_date` and only worked by accident when that variable happened
#' to exist in the calling environment. If `forecast_start_date` is
#' NULL, points fall back to the tail of the panel's own time axis
#' (the same convention `plot_filter_one()` uses).
#'
#' @param points the values being placed (only `length()` is used)
#' @param forecast_start_date numeric date (or integer time index) marking
#'   the last observed time point, or NULL
#' @param init_date "YYYY-MM-DD" string if the panel uses a date axis, else NULL
#' @param ax the axis spec returned by `.make_x_axis()` for this panel
#' @return a vector of x-positions, length `length(points)`
.anchor_points <- function(points, forecast_start_date, init_date, ax) {
  n_pts <- length(points)
  if (is.null(forecast_start_date)) return(tail(ax$x_vals, n_pts))

  if (!is.null(init_date)) {
    seq(as.Date(forecast_start_date, origin = "1970-01-01") + 1,
        by = "day", length.out = n_pts)
  } else {
    forecast_start_date + seq_len(n_pts)
  }
}

#' Age-stratified variant of `plot_filter_one()`.
#'
#' Draws quantile ribbons (25-75% inner band, 2.5-97.5% outer band) coloured
#' by the supplied `age_color`. Optionally overlays a dashed "true" trajectory
#' (for simulation studies), training/test observation points, and a vertical
#' forecast boundary.
#'
#' @param theta_mat N x T matrix for one age group
#' @param plot_title panel title (e.g. "Age group: 25-44")
#' @param age_color hex colour string for this age group
#' @param init_date "YYYY-MM-DD" string for date x-axis
#' @param forecast_start_date numeric date for vertical forecast boundary
#' @param observed_points training observations (solid black points), plotted
#'   at the start of the time axis
#' @param true_points test observations (open circles), plotted at the tail
#'   of the time axis (last `length(true_points)` time points)
#' @param true_trajectory known/simulated trajectory (dashed line), length T
#' @param true_color colour for `true_trajectory` / `true_points` (default "black")
#' @param ymin,ymax optional y-axis limits passed to `coord_cartesian()`
#' @return a ggplot object
plot_filter_one_age <- function(theta_mat,
                                 plot_title          = NULL,
                                 age_color            = "steelblue",
                                 init_date            = NULL,
                                 forecast_start_date  = NULL,
                                 observed_points      = NULL,
                                 true_points          = NULL,
                                 true_trajectory      = NULL,
                                 true_color           = "black",
                                 ymin                 = NULL,
                                 ymax                 = NULL) {

  T  <- ncol(theta_mat)
  ax <- .make_x_axis(T, init_date)

  quant_df <- do.call(rbind, lapply(seq_len(T), function(t) {
    v <- theta_mat[, t]
    data.frame(
      Time   = ax$x_vals[t],
      q025   = quantile(v, 0.025, na.rm = TRUE),
      q25    = quantile(v, 0.250, na.rm = TRUE),
      q75    = quantile(v, 0.750, na.rm = TRUE),
      q975   = quantile(v, 0.975, na.rm = TRUE),
      median = median(v, na.rm = TRUE)
    )
  }))

  p <- ggplot(quant_df, aes(Time)) +
    geom_ribbon(aes(ymin = q025, ymax = q975), fill = age_color, alpha = 0.22) +
    geom_ribbon(aes(ymin = q25,  ymax = q75),  fill = age_color, alpha = 0.35) +
    geom_line(aes(y = median), color = age_color, linewidth = 1.4) +
    labs(x = NULL, y = NULL, title = plot_title) +
    .smc2_theme() +
    theme(legend.position = "none",
          axis.title.x    = element_blank(),
          plot.title      = element_text(size = 24, hjust = 0.5))

  if (!is.null(ax$scale)) p <- p + ax$scale

  if (!is.null(forecast_start_date)) {
    p <- p + geom_vline(xintercept = forecast_start_date, linetype = "dashed", linewidth = 1.1)
  }

  if (!is.null(observed_points)) {
    ax_obs <- .make_x_axis(length(observed_points), init_date)
    p <- p + geom_point(
      data = data.frame(Time = ax_obs$x_vals, value = observed_points),
      aes(Time, value), color = "black", size = 1.6, alpha = 0.7, inherit.aes = FALSE
    )
  }

  if (!is.null(true_points)) {
    x_pts <- .anchor_points(true_points, forecast_start_date, init_date, ax)
    p <- p + geom_point(
      data = data.frame(Time = x_pts, value = true_points),
      aes(Time, value), shape = 21, fill = "white", color = true_color,
      size = 2.2, stroke = 0.8, inherit.aes = FALSE
    )
  }

  if (!is.null(true_trajectory)) {
    p <- p + geom_line(
      data = data.frame(Time = ax$x_vals, value = true_trajectory),
      aes(Time, value), color = true_color, linewidth = 1.3,
      linetype = "dashed", inherit.aes = FALSE
    )
  }

  if (!is.null(ymin) || !is.null(ymax)) p <- p + coord_cartesian(ylim = c(ymin, ymax))

  p
}

# =============================================================================
# FILTERING DISTRIBUTIONS - ALL PARAMETERS
# =============================================================================

#' Apply `plot_filter_one()` to every parameter dimension.
#'
#' @param smc2_result SMC^2 result list with `$theta` array (N x T x d)
#' @param opts options list; uses `opts$paramNames`
#' @param init_date "YYYY-MM-DD" string for date x-axis
#' @return named list of ggplot objects, one per parameter
plot_filter_all <- function(smc2_result, opts, init_date = NULL) {
  theta_hist  <- smc2_result$theta  # N x T x d
  d           <- dim(theta_hist)[3]
  param_names <- opts$paramNames %||% paste0("\u03b8_", seq_len(d))

  plots <- lapply(seq_len(d), function(j) {
    plot_filter_one(
      theta_mat  = theta_hist[, , j],
      y_label    = param_names[j],
      plot_title = if (is.character(param_names[j])) param_names[j] else NULL,
      init_date  = init_date
    )
  })

  names(plots) <- as.character(param_names)
  plots
}

# =============================================================================
# PRIOR vs POSTERIOR - SINGLE PARAMETER
# =============================================================================

#' Overlay prior density and posterior histogram for one parameter.
#'
#' @param samples numeric vector of posterior draws
#' @param weights optional importance weights (same length as `samples`)
#' @param prior_fn function `f(x)` returning prior density at `x`;
#'   if NULL only the posterior is shown
#' @param param_name label for the x-axis
#' @return a ggplot object
plot_prior_posterior_one <- function(samples, weights = NULL, prior_fn = NULL,
                                      param_name = "theta") {

  df        <- data.frame(value = samples)
  df$weight <- weights %||% 1

  p <- ggplot(df, aes(x = value)) +
    stat_density(
      aes(y = after_stat(density), fill = "Posterior", weight = weight),
      geom = "area", alpha = 0.5, color = "steelblue", linewidth = 1.1
    ) +
    labs(x = param_name, y = "Density") +
    .smc2_theme()

  if (!is.null(prior_fn)) {
    x_seq    <- seq(min(samples), max(samples), length.out = 200)
    prior_df <- data.frame(x = x_seq, y = prior_fn(x_seq))
    p <- p +
      geom_ribbon(data = prior_df, aes(x = x, ymin = 0, ymax = y, fill = "Prior"), alpha = 0.3) +
      geom_line(data = prior_df, aes(x = x, y = y), color = "firebrick",
                linewidth = 1.2, show.legend = FALSE)
  }

  p + scale_fill_manual(name = "", values = c("Posterior" = "steelblue", "Prior" = "firebrick"))
}

# =============================================================================
# PRIOR vs POSTERIOR - ALL PARAMETERS
# =============================================================================

#' Apply `plot_prior_posterior_one()` to every parameter dimension.
#'
#' Attempts to reconstruct prior density functions from `opts$paramPriors`
#' if `prior_fns` is not supplied directly.
#'
#' @param smc2_result SMC^2 result list with `$theta` array
#' @param opts options list; uses `opts$paramNames` and `opts$paramPriors`
#' @param prior_fns optional list of length d, each a density function; if
#'   NULL it is derived from `opts$paramPriors`
#' @return named list of ggplot objects, one per parameter
plot_prior_posterior_all <- function(smc2_result, opts, prior_fns = NULL) {
  theta_hist  <- smc2_result$theta
  d           <- dim(theta_hist)[3]
  param_names <- opts$paramNames %||% paste0("\u03b8_", seq_len(d))

  if (is.null(prior_fns) && !is.null(opts$paramPriors)) {
    prior_fns <- lapply(opts$paramPriors, function(p) {
      if (!is.list(p) || is.null(p$dist)) return(NULL)
      dfun <- tryCatch(get(paste0("d", p$dist)), error = function(e) NULL)
      if (is.null(dfun)) return(NULL)
      args      <- p
      args$dist <- NULL
      function(x) do.call(dfun, c(list(x = x), args))
    })
  }

  plots <- lapply(seq_len(d), function(j) {
    plot_prior_posterior_one(
      samples    = as.vector(theta_hist[, ncol(theta_hist), j]),
      prior_fn   = prior_fns[[j]],
      param_name = param_names[j]
    )
  })

  names(plots) <- as.character(param_names)
  plots
}

# =============================================================================
# SMC^2 ALGORITHM DIAGNOSTICS (4-panel summary figure)
# =============================================================================

#' Four-panel diagnostic figure for an SMC^2 run.
#'
#' Panels: effective sample size (ESS), number of particles (Nx),
#' MCMC acceptance rate (%), and log-marginal evidence.
#'
#' @param smc2_result list with `$ESS`, `$Nx`, `$acceptance_rate`, `$logEvidence`
#' @param init_date "YYYY-MM-DD" string for date x-axis; NULL for integer
#' @param ESS_thres ESS threshold (as a fraction of `opts$N_theta`) used to
#'   draw the reference line on the ESS panel
#' @return a single cowplot grid object containing all four panels
plot_SMC2_diagnostics <- function(smc2_result, init_date = NULL, ESS_thres = 0.4) {

  T  <- length(smc2_result$ESS)
  ax <- .make_x_axis(T, init_date)

  df <- data.frame(
    Time = ax$x_vals,
    ESS  = smc2_result$ESS,
    Nx   = smc2_result$Nx,
    acc  = smc2_result$acceptance_rate * 100,
    logZ = smc2_result$logEvidence
  )

  .panel <- function(y_col, color, ylab, title_str, extra = NULL) {
    p <- ggplot(df, aes(x = Time, y = .data[[y_col]])) +
      geom_line(color = color, linewidth = 1.2) +
      labs(title = title_str, x = ax$x_label, y = ylab) +
      .smc2_theme() + .small_generic("both")
    if (!is.null(ax$scale)) p <- p + ax$scale
    if (!is.null(extra))    p <- p + extra
    p
  }

  N_theta_max <- max(df$ESS, na.rm = TRUE)  # ESS starts at N_theta

  p1 <- .panel("ESS", "steelblue", "ESS", "Effective Sample Size") +
    geom_hline(yintercept = ESS_thres * N_theta_max, linetype = "dashed", color = "gray40")
  p2 <- .panel("Nx", "darkgreen", "Nx", "Number of Particles")

  p3 <- ggplot(df[!is.na(df$acc), ], aes(x = Time, y = acc)) +
    geom_line(color = "firebrick", linewidth = 1.2) +
    geom_hline(yintercept = 23.4, linetype = "dashed", color = "gray40") +
    ylim(0, 100) +
    labs(title = "MCMC Acceptance Rate", x = ax$x_label, y = "%") +
    .smc2_theme() + .small_generic("both")
  if (!is.null(ax$scale)) p3 <- p3 + ax$scale

  p4 <- .panel("logZ", "purple", "log Z", "Log Marginal Evidence")

  plot_grid(p1, p2, p3, p4, ncol = 2, align = "hv")
}

# =============================================================================
# GRID ASSEMBLY HELPERS
# =============================================================================

#' Arrange a list of plots in a grid with a single shared (enlarged) legend.
#'
#' @param plot_list list of ggplot objects
#' @param ncol number of columns in the grid
#' @return a cowplot grid object combining the panel grid and shared legend
.grid_shared_legend <- function(plot_list, ncol) {
  legend <- get_legend(
    plot_list[[1]] + theme(
      legend.position = "right",
      legend.title    = element_text(size = 24, face = "bold"),
      legend.text     = element_text(size = 24),
      legend.key.size = unit(2, "lines")
    )
  )

  plots_no_leg <- lapply(plot_list, function(p) p + theme(legend.position = "none"))
  panel_grid   <- plot_grid(plotlist = plots_no_leg, ncol = ncol, align = "hv")
  plot_grid(panel_grid, legend, ncol = 2, rel_widths = c(3, 0.15))
}

#' Arrange a list of plots in a grid with one shared y-axis label (and legend).
#'
#' Strips each panel's own axis titles and legend, lays them out in a grid,
#' and overlays a single rotated y-axis label and x-axis label plus one
#' shared legend (taken from the first plot).
#'
#' @param plot_list list of ggplot objects
#' @param ncol number of columns in the grid
#' @param y_label shared y-axis label
#' @param x_label shared x-axis label (default "Time")
#' @return a cowplot grid object
.grid_shared_ylabel <- function(plot_list, ncol, y_label, x_label = "Time") {
  no_lab <- lapply(plot_list, function(p) {
    p + theme(axis.title.y = element_blank(), axis.title.x = element_blank(),
              legend.position = "none")
  })

  legend <- get_legend(
    plot_list[[1]] + theme(
      legend.position = "right",
      legend.title    = element_text(size = 16, face = "bold"),
      legend.text     = element_text(size = 16),
      legend.key.size = unit(1.5, "lines")
    )
  )

  main_grid <- wrap_plots(no_lab, ncol = ncol)

  labeled <- ggdraw() +
    draw_plot(main_grid, x = 0.03, y = 0.06, width = 0.97, height = 0.94) +
    draw_label(y_label, x = 0.01, y = 0.53, angle = 90, fontface = "bold", size = 20) +
    draw_label(x_label, x = 0.52, y = 0.04, fontface = "bold", size = 20)

  plot_grid(labeled, legend, ncol = 2, rel_widths = c(1, 0.08))
}

#' Assemble age-stratified panels with shared axis labels (no legend).
#'
#' @param plot_list list of ggplot objects, one per age group
#' @param ncol number of columns in the grid
#' @param y_label shared y-axis label (string, expression, or bquote())
#' @param x_label shared x-axis label (default "Time"; use "Date" for
#'   date-indexed panels)
#' @return a cowplot grid object
.assemble_age_fig <- function(plot_list, ncol, y_label, x_label = "Time") {
  grid <- wrap_plots(plot_list, ncol = ncol) &
    theme(plot.margin = margin(20, 20, 20, 20))

  ggdraw() +
    draw_plot(grid, x = 0.02, y = 0.07, width = 0.96, height = 0.92) +
    draw_label(y_label, x = 0.015, y = 0.53, angle = 90, fontface = "bold", size = 30) +
    draw_label(x_label, x = 0.52,  y = 0.07, fontface = "bold", size = 28)
}
