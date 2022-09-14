# Install released version from CRAN
vignette("pkgdown")

# Run once to configure package to use pkgdown
usethis::use_pkgdown()
# Run to build the website
pkgdown::build_site(
    examples = FALSE,
    devel = TRUE,
    lazy = TRUE
)
