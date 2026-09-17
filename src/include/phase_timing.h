#ifndef NGME_PHASE_TIMING_H
#define NGME_PHASE_TIMING_H

// Cumulative wall time in the Cholesky factorizations themselves.
//
// The per-phase timing already in compute_grad_and_hessian buckets the QQ
// factorization inside "sampleW" and says nothing about the operator's, so it
// cannot answer how much of a fit is the factorizations -- the part no change
// to the estimator can remove. These accumulators separate exactly that.
//
// Atomics because the factorizations run inside the OpenMP loops over chains
// and replicates; this is a diagnostic, so contention on it is preferable to
// per-thread state that has to be reduced somewhere.
//
// COMPILED OUT BY DEFAULT. Scope sits inside loops that run once per Gibbs
// draw and once per parameter, and every one of them would take two clock
// reads and an atomic increment on a cache line shared by every OpenMP thread
// -- a diagnostic has no business costing that in an ordinary fit. Without
// NGME_PHASE_TIMING the Scope below is an empty object the compiler removes
// entirely and the counters stay at zero, so factorization_timing() still
// exists and simply reports nothing.
//
// To profile, rebuild with the macro defined, e.g.
//     PKG_CPPFLAGS=-DNGME_PHASE_TIMING R CMD INSTALL .
#include <atomic>
#ifdef NGME_PHASE_TIMING
#include <chrono>
#endif

namespace ngme_timing {

inline std::atomic<long long> &qq_numeric_us() {
  static std::atomic<long long> v{0};
  return v;
}
inline std::atomic<long long> &qq_symbolic_us() {
  static std::atomic<long long> v{0};
  return v;
}
inline std::atomic<long long> &k_numeric_us() {
  static std::atomic<long long> v{0};
  return v;
}
inline std::atomic<long long> &k_symbolic_us() {
  static std::atomic<long long> v{0};
  return v;
}

// Whole-iteration phases, so one optimizer step can be accounted end to end
// rather than only inside compute_grad_and_hessian. The factorization counters
// above are NESTED inside these: qq_numeric within grad_sampleW, k_numeric
// within op_trace.
#define NGME_PHASE(name)                                                       \
  inline std::atomic<long long> &name() {                                      \
    static std::atomic<long long> v{0};                                        \
    return v;                                                                  \
  }
NGME_PHASE(op_build_us)   // build_KZ at base theta
NGME_PHASE(op_dK_us)      // dK / dZ, analytic or by differencing
NGME_PHASE(op_trace_us)   // tr(K^-1 dK): shortcut, selinv or probes
NGME_PHASE(grad_sampleV_us)
NGME_PHASE(grad_sampleW_us)
NGME_PHASE(grad_rbtrace_us)
NGME_PHASE(grad_total_us)
// Inside the RB traces: the probe block solve QU = QQ^-1 U (once per
// factorization, N_iter right-hand sides) against the sparse-times-dense
// products each trace then does with it. Splitting these says whether the cost
// is the solve or the per-parameter algebra, which have different remedies.
NGME_PHASE(rb_qu_solve_us)
NGME_PHASE(rb_product_us)
NGME_PHASE(rb_calls)
// Sections of compute_rb_trace, to locate the >90% that is neither the probe
// solve nor the estimator products.
NGME_PHASE(rb_sec_K_us)       // per-latent theta_K traces
NGME_PHASE(rb_sec_sigma_us)   // per-latent theta_sigma traces
NGME_PHASE(rb_sec_Z_us)       // the dZ / ADA block
NGME_PHASE(rb_sec_noise_us)   // the final theta_sigma loop (rebuilds AZ)
NGME_PHASE(qq_assemble_us)    // Q = K^T diag(1/SV) K, and QQ = Q + measure
NGME_PHASE(qq_prod_us)        // just the K^T D K triple product
NGME_PHASE(qq_add_us)         // just QQ = Q + measurement block
NGME_PHASE(qq_measure_us)     // rebuilding H^T H, nested inside qq_add
NGME_PHASE(rmvn_us)           // the rMVN draw and conditional mean solves
NGME_PHASE(set_param_us)      // Ngme::set_parameter_and_update, whole
NGME_PHASE(samplew_us)        // BlockModel::sampleW_VY, whole
NGME_PHASE(sw_ensureQQ_us)    // ensure_QQ(): assemble + factorize when stale
NGME_PHASE(sw_M_us)           // the right-hand side M
NGME_PHASE(sw_G_us)           // get_G(inv_SV)
NGME_PHASE(sw_H_us)           // get_sqrt_AtSVA()
// The remaining sections of BlockModel::grad(). The std::chrono locals beside
// them accumulate MILLISECONDS, so on a fast model every pass truncates to
// zero and those sections read as an unattributed residual; these do not.
NGME_PHASE(grad_V_us)         // sample_cond_V + sample_cond_noise_V
NGME_PHASE(grad_score_us)     // building the observation score s_full
NGME_PHASE(grad_assemble_us)  // per-latent gradient aggregation
NGME_PHASE(grad_prec_lat_us)  // preconditioner: latent block
NGME_PHASE(grad_prec_ZGN_us)  // preconditioner: Z / GN block
NGME_PHASE(grad_prec_merr_us) // preconditioner: measurement block
NGME_PHASE(opt_step_us)       // the SGD step, which ENCLOSES every grad phase
#undef NGME_PHASE

// Scoped accumulator: adds its lifetime to the given counter.
#ifdef NGME_PHASE_TIMING
class Scope {
public:
  explicit Scope(std::atomic<long long> &sink)
      : sink_(&sink), t0_(std::chrono::steady_clock::now()) {}
  Scope(const Scope &) = delete;
  Scope &operator=(const Scope &) = delete;
  // End the measurement early. For a section that cannot simply be wrapped in
  // braces because the variables it declares are used after it.
  void stop() {
    if (sink_ == nullptr)
      return;
    sink_->fetch_add(std::chrono::duration_cast<std::chrono::microseconds>(
                         std::chrono::steady_clock::now() - t0_)
                         .count(),
                     std::memory_order_relaxed);
    sink_ = nullptr;
  }
  ~Scope() { stop(); }

private:
  std::atomic<long long> *sink_;
  std::chrono::steady_clock::time_point t0_;
};
// Counter increments made outside a Scope go through this, so they vanish with
// everything else rather than leaving stray atomics in the hot paths.
inline void add(std::atomic<long long> &sink, long long v) {
  sink.fetch_add(v, std::memory_order_relaxed);
}
#else
class Scope {
public:
  explicit Scope(std::atomic<long long> &) {}
  void stop() {}
};
inline void add(std::atomic<long long> &, long long) {}
#endif

} // namespace ngme_timing

#endif // NGME_PHASE_TIMING_H
