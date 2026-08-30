# The rational (fractional) operator is the one model family whose sparsity
# pattern is not fixed over the course of an optimization. The pattern of
# K = Pl (and of Z = Pr) is built from floor(beta) and from whether beta is an
# integer, so it changes exactly when the smoothness alpha = 2 * beta crosses
# an integer -- everywhere else only the numbers move.
#
# That is what the Cholesky caches have to respect: the symbolic (analyze)
# phase must be redone on a pattern change, and may be skipped otherwise. These
# tests pin down the two facts the caching logic relies on.

library(ngme2)
library(testthat)

sparsity_pattern <- function(M) {
  M <- methods::as(methods::as(M, "generalMatrix"), "CsparseMatrix")
  list(dim = dim(M), i = M@i, p = M@p)
}

fractional_ops_fixture <- function(n_node = 41) {
  x <- seq(0, 1, length.out = n_node)
  fem <- rSPDE::rSPDE.fem1d(x)
  lump <- mass_lump_matrices(fem$C)
  list(C = lump$C, Ci = lump$Ci, G = fem$G,
       B_kappa = Matrix::Matrix(matrix(1, nrow(fem$C), 1), sparse = TRUE))
}

test_that("the rational operator keeps its sparsity pattern within an integer interval", {
  skip_on_cran()
  skip_if_not_installed("rSPDE")
  load_cpp_bindings()
  fx <- fractional_ops_fixture()
  ops <- function(beta, m = 1L) cpp_fractional_ops_bind(
    fx$C, fx$Ci, fx$G, beta, m, 1, log(5), fx$B_kappa)

  # Two smoothnesses inside the same integer interval: same pattern, different
  # numbers. This is the common case, and the reason the symbolic phase can be
  # computed once and reused for the whole optimization.
  a <- ops(1.2)
  b <- ops(1.8)
  expect_identical(sparsity_pattern(a$Pl), sparsity_pattern(b$Pl))
  expect_identical(sparsity_pattern(a$Pr), sparsity_pattern(b$Pr))
  expect_false(isTRUE(all.equal(as.numeric(a$Pl@x), as.numeric(b$Pl@x))))
})

test_that("the rational operator changes sparsity when the smoothness crosses an integer", {
  skip_on_cran()
  skip_if_not_installed("rSPDE")
  load_cpp_bindings()
  fx <- fractional_ops_fixture()
  ops <- function(beta, m = 1L) cpp_fractional_ops_bind(
    fx$C, fx$Ci, fx$G, beta, m, 1, log(5), fx$B_kappa)

  # An integer beta takes the non-rational branch: no root factors at all, and
  # Pr collapses to the diagonal Phi.
  integer_beta <- ops(1)
  fractional_beta <- ops(1.2)
  expect_false(identical(sparsity_pattern(integer_beta$Pl),
                         sparsity_pattern(fractional_beta$Pl)))
  expect_lt(Matrix::nnzero(integer_beta$Pr), Matrix::nnzero(fractional_beta$Pr))

  # floor(beta) enters as the power A^(floor(beta) - 1), so stepping over
  # beta = 2 (alpha = 4) adds a factor and widens Pl.
  expect_gt(Matrix::nnzero(ops(2.4)$Pl), Matrix::nnzero(ops(1.8)$Pl))
  expect_false(identical(sparsity_pattern(ops(1.8)$Pl),
                         sparsity_pattern(ops(2.4)$Pl)))

  # The rational order m also fixes how many factors Pl and Pr carry, so a
  # different m is a different pattern for the same smoothness.
  expect_false(identical(sparsity_pattern(ops(1.2, 1L)$Pl),
                         sparsity_pattern(ops(1.2, 2L)$Pl)))
})
