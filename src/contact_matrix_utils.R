# ============================================================
# contact_matrix_utils.R
#
# Builds an age-group-aggregated, reciprocity-corrected contact matrix
# from a single-year-of-age contact matrix and population vector.
# This logic was previously duplicated verbatim in the real-data and
# simulated-data analysis scripts; it now lives in one place.
# ============================================================

#' Enforce population reciprocity on a single-year-of-age contact matrix
#'
#' Symmetrizes contacts so that total contacts from age i to age j equal
#' total contacts from age j to age i: `N[i] * F[i,j] == N[j] * F[j,i]`.
#'
#' @param F square contact matrix (rows = contactor age, cols = contactee age)
#' @param N population vector, `length(N) == nrow(F)`
#' @return reciprocity-corrected contact matrix, same dimensions as `F`
enforce_reciprocity <- function(F, N) {
  n <- nrow(F)
  F_rec <- matrix(0, n, n)
  for (i in seq_len(n)) {
    F_rec[i, ] <- (N * F[, i] + N[i] * F[i, ]) / (2 * N[i])
  }
  F_rec
}

#' Aggregate a single-year-of-age contact matrix into broad age groups
#'
#' @param F_rec reciprocity-corrected single-year contact matrix
#' @param N population vector, `length(N) == nrow(F_rec)`
#' @param groups named list of integer index vectors, one per age group
#'   (e.g. `list(G1 = 1:25, G2 = 26:45, ...)`)
#' @return list with `M` (population-weighted aggregated contact matrix,
#'   `length(groups) x length(groups)`) and `N_group` (aggregated
#'   population sizes, named by `groups`)
aggregate_contact_matrix <- function(F_rec, N, groups) {
  A <- length(groups)
  M <- matrix(0, A, A)

  for (a in seq_len(A)) {
    rows <- groups[[a]]
    for (b in seq_len(A)) {
      cols      <- groups[[b]]
      M[a, b]   <- sum(F_rec[rows, cols] * N[cols]) / sum(N[rows])
    }
  }

  list(
    M       = M,
    N_group = vapply(groups, function(idx) sum(N[idx]), numeric(1))
  )
}

#' End-to-end: single-year contact matrix + population -> age-group matrix
#'
#' @param F single-year-of-age contact matrix (rows = contactor, cols = contactee)
#' @param N single-year-of-age population vector
#' @param groups named list of age-group index vectors
#' @return list with `M` (aggregated contact matrix) and `N_group`
#'   (aggregated population sizes)
build_contact_matrix <- function(F, N, groups) {
  F_rec <- enforce_reciprocity(F, N)
  aggregate_contact_matrix(F_rec, N, groups)
}
