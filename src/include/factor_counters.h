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

// Work actually performed, as opposed to time taken. Wall clock is not a
// property of the fit: it moves with machine load, other processes and CPU
// throttling, so two runs of the same code are not comparable by it. These two
// counters, together with the factorization counts above, are: the cost of an
// iteration is a fixed part plus a part linear in each of them, with
// coefficients that are constant for a given model and data, so the counts are
// directly comparable between runs however the machine behaved.
extern std::atomic<long long> probe_solves; // probe columns solved against QQ
extern std::atomic<long long> gibbs_passes; // Gibbs sweeps in the gradient
extern std::atomic<long long> fisher_solves; // solves for the information estimate
extern std::atomic<long long> k_probe_solves;  // probe columns solved against K

// Which solver is drawing probes right now. The two are the same unit of work
// at very different cost, so they are counted apart; a thread-local marker set
// by the caller avoids giving the solver a member and shifting the layout of
// everything that holds one.
enum class probe_role { qq = 0, op = 1 };
extern thread_local probe_role current_probe_role;

// Marks a scope as drawing against the operator rather than the block
// precision. Saves and restores the previous value rather than resetting to a
// fixed one, so the guards nest: an inner scope ending must not hand the rest
// of an enclosing operator region back to the qq counter.
struct probe_role_scope {
  probe_role prev;
  explicit probe_role_scope(probe_role r) : prev(current_probe_role) {
    current_probe_role = r;
  }
  ~probe_role_scope() { current_probe_role = prev; }
  probe_role_scope(const probe_role_scope &) = delete;
  probe_role_scope &operator=(const probe_role_scope &) = delete;
};

inline void add(std::atomic<long long> &c, long long n) {
  c.fetch_add(n, std::memory_order_relaxed);
}

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
