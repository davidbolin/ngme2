# Reproduce the CRAN Fedora GCC tests

This uses R-hub's pinned Fedora 44 x86_64 image with R-devel r90590 and GCC
16.2.1. It installs hard dependencies, `testthat`, and the
`ngme2_1.0.0.tar.gz` source package. Other optional `Suggests` packages,
including INLA, are not installed and tests requiring them may skip.
By default it runs `R CMD check --as-cran` (without manual or vignette builds)
and leaves the full `testthat.Rout` on the host.

Run on a native x86_64 Docker host:

```sh
./tools/cran-fedora-repro/run.sh
```

To locate a slow test, rerun against the same output directory in `profile`
mode:

```sh
./tools/cran-fedora-repro/run.sh /absolute/path/ngme2_1.0.0.tar.gz /absolute/output/directory profile
```

The profiler runs each test file in CRAN mode with a ten-minute per-file
limit. It writes `test-profiles/summary.tsv` and one log per test file. Exit
status 124 means that file reached the limit. The first argument specifies
the tarball, the second the output directory, and the third the mode.

The image is close to CRAN's Fedora GCC machine, but is not identical: the
R-hub image, installed dependency versions, available CPU cores, and test
execution order can differ. Running x86_64 containers on Apple Silicon uses
emulation, so elapsed times there cannot be compared with CRAN's timings.

The `CRAN Fedora test profiler` GitHub Actions workflow downloads the released
1.0.0 source tarball from CRAN and runs this image on a native x86_64 runner.
It runs automatically when these diagnostic files are pushed to `devel` and
can also be started manually in `profile` or `check` mode. Per-file logs and
the timing summary are uploaded as a workflow artifact even when a test fails.
