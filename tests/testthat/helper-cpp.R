cpp_bindings_env <- new.env(parent = emptyenv())

load_cpp_bindings <- function() {
  testthat::skip_on_cran()
  if (!requireNamespace("Rcpp", quietly = TRUE)) skip("Rcpp not available")
  if (!requireNamespace("RcppEigen", quietly = TRUE)) skip("RcppEigen not available")

  cpp_file <- testthat::test_path("cpp", "fractional_ops_bindings.cpp")

  if (!file.exists(cpp_file)) {
    skip("C++ binding file not found")
  }

  src_dir <- normalizePath(
    file.path(testthat::test_path("..", "..", "src")),
    mustWork = FALSE
  )

  if (!dir.exists(src_dir)) {
    src_dir <- normalizePath(
      file.path(dirname(cpp_file), "..", "..", "..", "src"),
      mustWork = FALSE
    )
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
