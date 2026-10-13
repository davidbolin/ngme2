# Fedora GCC and Clang checks

Build the current source package, including the short documentation vignette:

```sh
R CMD build .
./tools/cran-fedora-repro/run.sh ngme2_1.0.1.tar.gz /absolute/output/gcc check gcc
./tools/cran-fedora-repro/run.sh ngme2_1.0.1.tar.gz /absolute/output/clang check clang
```

Both checks use the same pinned Fedora 44 / R-devel R-hub GCC 16 image.
The Clang configuration installs Fedora's Clang, libc++ and libomp, and uses
them for ngme2 and its R dependencies. It builds a pinned GDAL 3.12.4 with
Clang/libc++ so that sf uses a compatible C++ interface; optional GDAL drivers
are disabled in this check environment. Mixing GCC-built Rcpp with a
Clang/libc++ package can crash when a C++ exception is converted to an R error,
so both images run a small Rcpp exception preflight. R itself remains GCC-built;
this is not an exact copy of CRAN's Clang machine. The R-hub `clang*` images are
Ubuntu images and are not used here.

All CRAN-hosted Imports, LinkingTo and Suggests are installed. INLA is omitted
because it comes from a separate repository; its absence can produce a NOTE.
The check uses `R CMD check --no-manual --as-cran`, rebuilds the documentation
vignette, preserves full test output and limits the container to two CPUs.
The PDF manual is not checked by this harness.

The `CRAN Fedora checks` workflow builds the checked-out source once, then
checks that exact tarball under both compilers on native x86_64 runners.
It runs on pushes to `cran` or by manual dispatch. Artifacts include the source
tarball and SHA-256, environment details, dependency/build logs, `00install.out`,
`00check.log`, and complete `testthat.Rout` output. Failed checks and one-hour
check timeouts fail the workflow and still upload the available logs.
Temporary installed-package and unpacked-source directories are excluded from
artifacts. This avoids case-sensitive HTML help filenames colliding when the
logs are downloaded on a case-insensitive macOS filesystem.

For a toolchain-only retry, select one `compiler` and provide the previous
`source_run_id` in the manual workflow inputs. This reuses the same tarball
without recompiling the package to build its vignette. Leave the run ID empty
when package source files have changed.

For individual test-file timings, use `profile` instead of `check`. Each file
runs in a separate R session in CRAN mode, with a ten-minute timeout (exit
status 124). Results are written to `test-profiles/summary.tsv` and one log per
file. Override `NGME_TEST_TIMEOUT_SECS` or `NGME_CHECK_TIMEOUT_SECS` as needed.

Running x86_64 containers on Apple Silicon uses emulation, so local elapsed
times cannot be compared directly with CRAN's timings.
The pinned image's Rscript also exited with status 132 under the tested Apple
Silicon emulation setup. Use the native x86_64 workflow if this happens locally.

Select `toolchain-probe` in the manual workflow inputs to check two regressions
without compiling ngme2: GCC-built Rcpp with Clang/libc++ must fail the exception
probe, while Clang/libstdc++ must succeed. The RNG probe records whether the
native libc++ diagnoses an ignored result, then requires explicit discard to
compile under `-Werror=unused-result` and preserve the generator state. libc++
22 does not annotate this operator as nodiscard; the original CRAN warning
depends on its library version. The deliberately inconsistent ABI case does
not indicate a failure in ngme2.

## Validation for ngme2 1.0.1

On 2026-10-13, the same source tarball passed both native x86_64 Fedora 44
checks with R-devel (2026-10-10 r90659) and the two-CPU limit:

| Compiler | R CMD check result | Test elapsed time | Logs |
| --- | --- | ---: | --- |
| GCC 16.2.1 | 0 errors, 0 warnings, 1 NOTE | 92 seconds | [GCC job](https://github.com/davidbolin/ngme2/actions/runs/38379082554/job/115234811803) |
| Clang 22.1.8 with libc++ | 0 errors, 0 warnings, 1 NOTE | 61 seconds | [Clang job](https://github.com/davidbolin/ngme2/actions/runs/38386109659/job/115258280541) |

Source SHA-256:
`f45a9dbb561d4ae41b8a01c2244b7dd085e074c38b70cbb12db9320e3ae14242`.
The incoming NOTE concerns resubmission after archival and the optional INLA
dependency. These counts are R CMD check results, not a claim that every
compiler diagnostic has been eliminated: legacy Eigen/SuiteSparse and unused
variable diagnostics remain in the installation logs. The original CRAN
90-minute timeout has not been reproduced on its exact hosts.

The [isolated toolchain probes](https://github.com/davidbolin/ngme2/actions/runs/38386107195)
passed: the inconsistent Rcpp ABI crashed with status 139, the consistent ABI
converted both exception types correctly, and explicit RNG discard compiled
under strict warnings and preserved generator state.
