# CRAN resubmission for ngme2 1.0.1

Please consider ngme2 1.0.1 for reinstatement on CRAN following its archival on
2026-10-11. This patch release updates the tests and fixes the reported Clang
warning without changing the public R API.

## Changes

* Limit the multichain estimation fits in the automated tests to two threads,
  preserving the number of chains and the existing test cases.
* Explicitly discard an unused random-number draw. The generator is still
  advanced, so this change preserves its random stream.
* Update the CHOLMOD integration for Eigen 5.0.1 and declare Matrix in LinkingTo.

## Check results

Both Fedora checks ran on native x86_64 with a two-CPU container limit, using
R-devel (2026-10-10 r90659) and `R CMD check --no-manual --as-cran`.

| Environment | Errors | Warnings | Notes | Test elapsed time |
| --- | ---: | ---: | ---: | ---: |
| Fedora 44 with GCC 16.2.1 | 0 | 0 | 1 | 92 seconds |
| Fedora 44 with Clang 22.1.8 and libc++ | 0 | 0 | 1 | 61 seconds |

Examples, compiled-code checks, tests and vignette rebuilding all passed.
The non-finite optimizer step regression test passed with both compilers.
All CRAN-hosted direct dependencies, including Suggests, were available; INLA
was not installed. The same source tarball was checked in both environments:

`f45a9dbb561d4ae41b8a01c2244b7dd085e074c38b70cbb12db9320e3ae14242`

The PDF manual was generated successfully on macOS arm64 with R 4.6.0 and
TeX Live 2026. Windows checks are not included in these results.

## Incoming note

The incoming NOTE identifies this as a new submission of an archived package
and lists INLA as a suggested package outside the mainstream repositories.
INLA is optional and its repository is declared in Additional_repositories.
Checks that depend on INLA are skipped when it is unavailable.

## Reproduction scope

The Fedora environments are pinned containers, not exact copies of CRAN's
hosts. R itself is GCC-built. The Clang environment builds Rcpp and the other
installed C++ dependencies with Clang/libc++, including a matching GDAL for sf.
Optional GDAL drivers are disabled in this check environment.

The original 90-minute timeout was not reproduced on the exact CRAN hosts.
The results above establish that the current tarball completes in the stated
environments; they do not establish the sole cause of the original timeout.

Check logs: [GCC](https://github.com/davidbolin/ngme2/actions/runs/38379082554/job/115234811803),
[Clang](https://github.com/davidbolin/ngme2/actions/runs/38386109659/job/115258280541).
