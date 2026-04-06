fixture_path <- function(...) {
  testthat::test_path("fixtures", ...)
}

has_fixture <- function(...) {
  file.exists(fixture_path(...))
}
