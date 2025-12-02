cpp_bindings_env <- new.env(parent = emptyenv())
load_cpp_bindings <- function() {
  testthat::skip_on_cran()
  if (!requireNamespace("Rcpp", quietly = TRUE)) skip("Rcpp not available")
  if (!requireNamespace("RcppEigen", quietly = TRUE)) skip("RcppEigen not available")

  cpp_file <- testthat::test_path("cpp", "fractional_ops_bindings.cpp")

  if (!file.exists(cpp_file)) {
    skip("C++ binding file not found")
  }

  # Set include path for Rcpp to find headers in src/
  # We need to point to the package root's src directory
  # testthat runs from tests/testthat, so src is at ../../src
  # But we also need to handle if it's run from package root

  src_dir <- normalizePath(file.path(testthat::test_path("..", "..", "src")), mustWork = FALSE)
  # If that didn't work (e.g. installed package?), try to find it relative to the cpp file
  if (!dir.exists(src_dir)) {
    src_dir <- normalizePath(file.path(dirname(cpp_file), "..", "..", "..", "src"), mustWork = FALSE)
  }

  rcpp_eigen_include <- system.file("include", package = "RcppEigen")
  Sys.setenv(PKG_CPPFLAGS = paste0(
    "-I\"", src_dir, "\" ",
    "-I\"", file.path(src_dir, "latents", "fractional"), "\" ",
    "-I\"", rcpp_eigen_include, "\""
  ))

  if (!isTRUE(cpp_bindings_env$loaded)) {
    Rcpp::sourceCpp(cpp_file, verbose = FALSE, rebuild = TRUE)
    cpp_bindings_env$loaded <- TRUE
  }
}

mass_lump_matrices <- function(C) {
  diag_vals <- as.numeric(Matrix::rowSums(C))
  C_lump <- methods::as(Matrix::Diagonal(x = diag_vals), "dgCMatrix")
  Ci_lump <- methods::as(Matrix::Diagonal(x = 1 / diag_vals), "dgCMatrix")
  list(C = C_lump, Ci = Ci_lump)
}


test_that("C++ Pl/Pr align with fractional.operators for nonstationary kappa", {
  load_cpp_bindings()
  if (!requireNamespace("rSPDE", quietly = TRUE)) skip("rSPDE not available")

  # Build a simple 1D FEM
  x <- seq(0, 1, length.out = 31)
  fem <- rSPDE::rSPDE.fem1d(x)
  C <- fem$C
  G <- fem$G
  n <- nrow(C)
  lump <- mass_lump_matrices(C)
  C_lump <- lump$C
  Ci_lump <- lump$Ci

  # Nonstationary kappa(s) = exp(B_kappa * theta_kappa)
  # Use intercept + linear term in s
  B_kappa <- Matrix::Matrix(cbind(1, x), sparse = TRUE)
  theta_kappa <- c(log(10), 0.2) # mild spatial variation

  # Compute kappa and L for the R reference
  kappa <- as.numeric(exp(B_kappa %*% theta_kappa))
  kp <- Matrix::Diagonal(n, kappa^2)
  L <- G + C_lump %*% kp

  # Model parameters
  beta <- 0.9
  m <- 2
  tau <- 1.3

  # Roots for the rational approximation from package internals
  roots <- rSPDE:::get.roots(m, beta)

  # R reference using fractional.operators
  c_sf <- min(kappa)^2
  op <- rSPDE::fractional.operators(L = L, beta = beta, C = C_lump, scale.factor = c_sf, m = m, tau = tau)

  # Compile and call C++ wrapper via sourceCpp
  res <- cpp_fractional_ops_with_roots_bind(
    C_lump, Ci_lump, G, beta, m, tau, theta_kappa, B_kappa,
    roots$rb, roots$rc, roots$factor
  )

  # Compare Pl and Pr
  expect_equal(as.matrix(op$Pl), as.matrix(res$Pl), tolerance = 1e-7)
  expect_equal(as.matrix(op$Pr), as.matrix(res$Pr), tolerance = 1e-7)

  betas <- c(0.9, 1.05, 2.8)
  ms <- c(1, 2, 3)
  tau <- 1.3 # Assuming tau is constant for these tests

  for (beta_val in betas) {
    for (m_val in ms) {
      # Roots for the rational approximation from package internals
      roots <- rSPDE:::get.roots(m_val, beta_val)

      # R reference using fractional.operators
      c_sf <- min(kappa)^2
      op <- rSPDE::fractional.operators(L = L, beta = beta_val, C = C_lump, scale.factor = c_sf, m = m_val, tau = tau)

      res_cpp <- cpp_fractional_ops_bind(
        C_lump, Ci_lump, G, beta_val, m_val, tau, theta_kappa, B_kappa
      )
      expect_equal(as.matrix(op$Pl), as.matrix(res_cpp$Pl),
        tolerance = 1e-7,
        info = paste0("Pl for beta=", beta_val, ", m=", m_val)
      )
      expect_equal(as.matrix(op$Pr), as.matrix(res_cpp$Pr),
        tolerance = 1e-7,
        info = paste0("Pr for beta=", beta_val, ", m=", m_val)
      )
    }
  }
})


test_that("test fractional operator", {
  # Build a simple 1D FEM
  x <- seq(0, 1, length.out = 31)
  fem <- rSPDE::rSPDE.fem1d(x)
  C <- fem$C
  G <- fem$G
  n <- nrow(C)
  lump <- mass_lump_matrices(C)
  C_lump <- lump$C
  Ci_lump <- lump$Ci

  # Nonstationary kappa(s) = exp(B_kappa * theta_kappa)
  # Use intercept + linear term in s
  B_kappa <- Matrix::Matrix(cbind(1, x), sparse = TRUE)
  theta_kappa <- c(log(10), 0.2) # mild spatial variation

  # Compute kappa and L for the R reference
  kappa <- as.numeric(exp(B_kappa %*% theta_kappa))
  kp <- Matrix::Diagonal(n, kappa^2)
  L <- G + C_lump %*% kp

  # Model parameters
  betas <- c(0.3, 1.2, 2.1, 3.2)
  ms <- c(1, 2, 3)
  tau <- 1.3

  for (beta in betas) {
    for (m in ms) {
      # Roots for the rational approximation from package internals
      roots <- rSPDE:::get.roots(m, beta)

      # R reference using fractional.operators
      c_sf <- min(kappa)^2
      op <- rSPDE::fractional.operators(L = L, beta = beta, C = C_lump, scale.factor = c_sf, m = m, tau = tau)

      K <- op$Pl
      Z <- op$Pr
      # compute Q = t(K) D K + t(Z) D Z
      d <- rexp(n)
      D <- diag(d)
      KtK <- Matrix::t(K) %*% D %*% K
      ZtZ <- Matrix::t(Z) %*% D %*% Z
      Q <- KtK + ZtZ

      expect_no_error(chol_Q <- Matrix::Cholesky(Q))
    }
  }
})
