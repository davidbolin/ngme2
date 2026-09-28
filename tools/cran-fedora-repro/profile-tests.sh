#!/usr/bin/env bash
set -u

test_dir="${NGME_TEST_SOURCE:-/opt/ngme2-source/ngme2}/tests/testthat"
log_dir="${NGME_TEST_LOG_DIR:-/check/test-profiles}"
limit="${NGME_TEST_TIMEOUT_SECS:-600}"

mkdir -p "$log_dir"
printf 'file\tstatus\telapsed_seconds\n' > "$log_dir/summary.tsv"

failed=0
for test_file in "$test_dir"/test-*.R; do
  base="${test_file##*/}"
  stem="${base#test-}"
  stem="${stem%.R}"
  start="$(date +%s)"
  printf 'START %s %s\n' "$stem" "$(date -u +%FT%TZ)"

  if env NGME_TEST_DIR="$test_dir" NGME_TEST_FILTER="$stem" \
    NOT_CRAN=false TESTTHAT_IS_CHECKING=true \
    timeout --signal=TERM --kill-after=10s "${limit}s" \
    Rscript --vanilla -e 'testthat::test_dir(Sys.getenv("NGME_TEST_DIR"), filter = paste0("^", Sys.getenv("NGME_TEST_FILTER"), "$"), package = "ngme2", load_package = "installed", reporter = "summary")' \
    > "$log_dir/$stem.log" 2>&1; then
    status=0
  else
    status=$?
    failed=1
  fi

  elapsed="$(($(date +%s) - start))"
  printf '%s\t%s\t%s\n' "$stem" "$status" "$elapsed" >> "$log_dir/summary.tsv"
  printf 'END %s status=%s elapsed=%ss\n' "$stem" "$status" "$elapsed"
done

exit "$failed"
