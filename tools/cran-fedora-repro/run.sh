#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
package_root="$(cd "$script_dir/../.." && pwd)"
version="$(sed -n 's/^Version: //p' "$package_root/DESCRIPTION")"
tarball="${1:-$package_root/ngme2_${version}.tar.gz}"
output_dir="${2:-$(mktemp -d "${TMPDIR:-/tmp}/ngme2-fedora.XXXXXX")}"
mode="${3:-check}"
compiler="${4:-gcc}"

if [[ ! -f "$tarball" ]]; then
  printf 'Package tarball not found: %s\nRun R CMD build . first.\n' "$tarball" >&2
  exit 2
fi
if [[ "$mode" != check && "$mode" != profile ]]; then
  printf 'Mode must be check or profile, got: %s\n' "$mode" >&2
  exit 2
fi
if [[ "$compiler" != gcc && "$compiler" != clang ]]; then
  printf 'Compiler must be gcc or clang, got: %s\n' "$compiler" >&2
  exit 2
fi

mkdir -p "$output_dir"
output_dir="$(cd "$output_dir" && pwd)"
# Use a version-independent filename only inside the Docker build context.
# R CMD check below sees the original PACKAGE_VERSION.tar.gz filename.
cp "$tarball" "$output_dir/ngme2-source.tar.gz"
tarball_name="${tarball##*/}"
if [[ "$tarball" != "$output_dir/$tarball_name" ]]; then
  cp "$tarball" "$output_dir/$tarball_name"
fi
cp "$script_dir/Dockerfile" "$script_dir/Makevars.clang" \
  "$script_dir/install-deps.R" "$script_dir/profile-tests.sh" \
  "$script_dir/build-gdal-clang.sh" "$script_dir/check-toolchain.R" "$output_dir/"

image="ngme2-cran-fedora-$compiler:$version"
printf 'Output directory: %s\nCompiler: %s\n' "$output_dir" "$compiler"
docker build --platform linux/amd64 --progress plain \
  --build-arg "NGME_COMPILER=$compiler" --tag "$image" \
  "$output_dir" 2>&1 | tee "$output_dir/docker-build.log"

docker run --rm --init --platform linux/amd64 "$image" \
  sh -c 'cat /etc/os-release; R --version; R CMD config CC; R CMD config CXX; R CMD config CXXFLAGS; Rscript -e "sessionInfo()"; rpm -q gcc clang libcxx libomp; gdal-config --version; gdal-config --prefix; true' \
  > "$output_dir/environment.txt" 2>&1
shasum -a 256 "$output_dir/$tarball_name" > "$output_dir/source-sha256.txt"

if [[ "$mode" == check ]]; then
  # GNU timeout must not be PID 1: in this Fedora image it exits 125 before
  # starting the child command. Docker's init also forwards termination signals.
  docker run --rm --init --platform linux/amd64 --cpus 2 \
    --env "_R_CHECK_FORCE_SUGGESTS_=false" \
    --mount "type=bind,source=$output_dir,target=/check" "$image" \
    timeout --signal=TERM --kill-after=30s "${NGME_CHECK_TIMEOUT_SECS:-3600}s" \
    R CMD check --no-manual --as-cran "/check/$tarball_name" \
    2>&1 | tee "$output_dir/check-console.log"
  printf 'Full test output: %s/ngme2.Rcheck/tests/testthat.Rout\n' "$output_dir"
else
  docker run --rm --init --platform linux/amd64 --cpus 2 \
    --env "NGME_TEST_TIMEOUT_SECS=${NGME_TEST_TIMEOUT_SECS:-600}" \
    --mount "type=bind,source=$output_dir,target=/check" "$image" \
    bash -c 'R CMD INSTALL /check/ngme2-source.tar.gz && /usr/local/bin/profile-ngme2-tests' \
    2>&1 | tee "$output_dir/profile-console.log"
  printf 'Per-file timing: %s/test-profiles/summary.tsv\n' "$output_dir"
fi
