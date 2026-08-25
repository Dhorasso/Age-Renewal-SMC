# ============================================================
# epi_diagnostics.R
#
# Derived epidemiological quantities computed from posterior particle
# trajectories: the next-generation matrix (NGM), the effective
# reproduction number R_t (its dominant eigenvalue), and the
# age-specific reproduction contributions R_{a,t}.

# ============================================================

#' Next-generation matrix at a single time point
#'
#' @param beta_t numeric vector of length A, age-specific transmissibility
#' @param C contact matrix, `C[b, a]` = contacts received by a from b
#' @return matrix K, `A x A`, with `K = diag(beta_t) %*% t(C)`
next_gen_matrix <- function(beta_t, C) {
  
  beta_t <- as.numeric(beta_t)
  C <- as.matrix(C)
  storage.mode(C) <- "numeric"
  
  sweep(t(C), 1, beta_t, "*")
}

#' Compute R_t and R_{a,t} across all particles and time points
#'
#' @param beta_array array `[N, T, A]` of age-specific transmissibility
#'   trajectories (e.g. `X_marginal[, , 2*seq_len(A) - 1]`)
#' @param C contact matrix, `A x A`, `C[b, a]` = contacts received by a from b
#' @return list with `Rt` (`N x T` matrix) and `Rat` (`N x T x A` array)
compute_Rt_Rat <- function(beta_array, C) {
  N <- dim(beta_array)[1]
  T <- dim(beta_array)[2]
  A <- dim(beta_array)[3]

  Rt  <- matrix(NA_real_, N, T)
  Rat <- array(NA_real_, dim = c(N, T, A))

  for (n in seq_len(N)) {
    for (t in seq_len(T)) {
      K          <- next_gen_matrix(beta_array[n, t, ], C)
      Rt[n, t]   <- Re(eigen(K, only.values = TRUE, symmetric = FALSE)$values[1])
      Rat[n, t, ] <- colSums(K)
    }
  }

  list(Rt = Rt, Rat = Rat)
}

#' Compute R_t and R_{a,t} for a single (e.g. "true") beta trajectory
#'
#' @param beta_mat matrix `T x A` of age-specific transmissibility
#' @param C contact matrix, `A x A`
#' @return list with `Rt` (numeric vector, length T) and `Rat`
#'   (matrix `T x A`)
compute_Rt_Rat_single <- function(beta_mat, C) {
  T <- nrow(beta_mat)
  A <- ncol(beta_mat)

  Rt  <- numeric(T)
  Rat <- matrix(NA_real_, T, A)

  for (t in seq_len(T)) {
    K       <- next_gen_matrix(beta_mat[t, ], C)
    Rt[t]   <- Re(eigen(K, only.values = TRUE, symmetric = FALSE)$values[1])
    Rat[t, ] <- colSums(K)
  }

  list(Rt = Rt, Rat = Rat)
}

#' Sort-then-sum aggregation of per-age particle matrices into a
#' population-level total
#'
#' Sorting each age-group slice along the particle dimension before
#' summing preserves rank-ordering within each time step, avoiding the
#' cross-particle pairing artefacts that a naive sum would introduce
#' when constructing credible intervals for the aggregate.
#'
#' @param mats list of `A` matrices, each `N x T` (one per age group)
#' @return matrix `N x T`, the population-level total
sort_then_sum <- function(mats) {
  N <- nrow(mats[[1]])
  T <- ncol(mats[[1]])
  A <- length(mats)

  sorted <- array(0, dim = c(N, T, A))
  for (a in seq_len(A)) {
    for (t in seq_len(T)) {
      sorted[, t, a] <- sort(mats[[a]][, t], decreasing = TRUE)
    }
  }
  apply(sorted, c(1, 2), sum)
}
