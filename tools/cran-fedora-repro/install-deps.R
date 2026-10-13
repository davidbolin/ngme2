args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L, file.exists(args[[1L]]))
source_dir <- tempfile("ngme2-deps-")
dir.create(source_dir)
untar(args[[1L]], files = "ngme2/DESCRIPTION", exdir = source_dir)
description <- read.dcf(file.path(source_dir, "ngme2", "DESCRIPTION"))
fields <- intersect(c("Depends", "Imports", "LinkingTo", "Suggests"),
                    colnames(description))
packages <- unique(trimws(gsub("\\s*\\([^)]*\\)", "",
                              unlist(strsplit(description[1L, fields], ",")))))
# INLA is not distributed by CRAN. Install every other direct dependency,
# including Suggests, rather than silently skipping tests for absent packages.
packages <- setdiff(packages, c("R", "INLA", ""))
options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 600)
cat("Installing dependencies for ngme2", description[1L, "Version"], "\n")
print(packages)
# Build dependencies from CRAN source using Fedora's default GCC toolchain,
# not the R-hub binary repository. Docker switches to Clang only afterwards.
pak::pkg_install(packages, dependencies = NA)
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing dependencies: ", paste(missing, collapse = ", "))
