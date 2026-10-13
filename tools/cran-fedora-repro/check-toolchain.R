cat("Checking Rcpp error conversion with", R.version.string, "\n")
cat("Rcpp version:", as.character(packageVersion("Rcpp")), "\n")
Rcpp::cppFunction('void ngme2_toolchain_exception(bool rcpp_error) {
  if (rcpp_error) Rcpp::stop("ngme2 toolchain probe");
  throw std::runtime_error("ngme2 toolchain probe");
}')
for (rcpp_error in c(FALSE, TRUE)) {
  message <- tryCatch({
    ngme2_toolchain_exception(rcpp_error)
    NA_character_
  }, error = conditionMessage)
  stopifnot(identical(message, "ngme2 toolchain probe"))
}
cat("Rcpp and standard C++ exceptions converted to R errors successfully\n")
