#ifndef __ngme2_factor_counters__
#define __ngme2_factor_counters__

#include <atomic>

// Instrumentation for the factorization caches.
//
// These counters exist so that the test-suite can verify that the caching
// introduced for QQ (BlockModel) and K (Operator) is actually effective:
// a Gaussian-noise model must assemble and factorize QQ once per optimizer
// iteration rather than once per Gibbs draw, and the symbolic (analyze) phase
// must only rerun when the sparsity pattern really changes.
//
// They are process-wide and are incremented from OpenMP worker threads, hence
// the relaxed atomics.
namespace ngme_counters {

extern std::atomic<long long> QQ_builds;    // QQ assembled + numerically factorized
extern std::atomic<long long> QQ_analyzes;  // symbolic phase for QQ
extern std::atomic<long long> K_analyzes;   // symbolic phase for a latent operator K

inline void bump(std::atomic<long long> &c) {
  c.fetch_add(1, std::memory_order_relaxed);
}

void reset_all();

// Test hook: when the environment variable NGME2_DISABLE_FACTOR_CACHE is set
// to a non-empty value other than "0"/"false", every cache introduced for
// speed (QQ, its Cholesky factor, the AZ block matrix, and the symbolic phase
// of the operator factorizations) is bypassed and the matrices are rebuilt on
// every use, reproducing the uncached code path. The test-suite fits the same
// models both ways and requires identical results, which keeps the caching
// honest without freezing golden numbers into the tests.
bool cache_disabled();


} // namespace ngme_counters

#endif
