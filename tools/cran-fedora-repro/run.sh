#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
package_root="$(cd "$script_dir/../.." && pwd)"
tarball="${1:-$package_root/ngme2_1.0.0.tar.gz}"
output_dir="${2:-$(mktemp -d "${TMPDIR:-/tmp}/ngme2-fedora-gcc16.XXXXXX")}"
mode="${3:-check}"

if [[ ! -f "$tarball" ]]; then
  printf 'Package tarball not found: %s\n' "$tarball" >&2
  exit 2
fi

if [[ "$mode" != check && "$mode" != profile ]]; then
  printf 'Mode must be check or profile, got: %s\n' "$mode" >&2
  exit 2
fi

mkdir -p "$output_dir"
output_dir="$(cd "$output_dir" && pwd)"
cp "$tarball" "$output_dir/ngme2_1.0.0.tar.gz"
cp "$script_dir/Dockerfile" "$output_dir/Dockerfile"
cp "$script_dir/profile-tests.sh" "$output_dir/profile-tests.sh"

printf 'Output directory: %s\n' "$output_dir"
docker build --platform linux/amd64 \
  --tag ngme2-cran-fedora-gcc16:1.0.0 \
  "$output_dir"

if [[ "$mode" == check ]]; then
  docker run --rm --platform linux/amd64 \
    --env "_R_CHECK_FORCE_SUGGESTS_=false" \
    --mount "type=bind,source=$output_dir,target=/check" \
    ngme2-cran-fedora-gcc16:1.0.0 \
    R CMD check --no-manual --no-build-vignettes --as-cran /check/ngme2_1.0.0.tar.gz
  printf 'Full test output: %s/ngme2.Rcheck/tests/testthat.Rout\n' "$output_dir"
else
  docker run --rm --platform linux/amd64 \
    --env "NGME_TEST_TIMEOUT_SECS=${NGME_TEST_TIMEOUT_SECS:-600}" \
    --mount "type=bind,source=$output_dir,target=/check" \
    ngme2-cran-fedora-gcc16:1.0.0
  printf 'Per-file timing: %s/test-profiles/summary.tsv\n' "$output_dir"
fi
