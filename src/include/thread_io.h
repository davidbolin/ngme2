#ifndef NGME_THREAD_IO_H
#define NGME_THREAD_IO_H

#include <Rcpp.h>
#include <ostream>
#include <streambuf>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace ngme_io {

// True only on the thread R itself runs on.
//
// OpenMP numbers the thread that encounters a parallel construct 0 within the
// new team, so the main thread is the one whose ancestor is 0 at every nesting
// level. Checking omp_get_thread_num() alone is not enough: estimate.cpp calls
// omp_set_max_active_levels(2), so chains can be nested inside replicates and
// thread 0 of an inner team may itself be an outer worker.
inline bool on_main_thread() {
#ifdef _OPENMP
  if (!omp_in_parallel())
    return true;
  for (int level = omp_get_level(); level > 0; --level)
    if (omp_get_ancestor_thread_num(level) != 0)
      return false;
  return true;
#else
  return true;
#endif
}

namespace detail {
class null_buffer : public std::streambuf {
protected:
  int overflow(int c) override { return c; }
};

// One sink per thread, so the discarding path shares no mutable state.
inline std::ostream &null_stream() {
  static thread_local null_buffer buf;
  static thread_local std::ostream sink(&buf);
  return sink;
}
} // namespace detail

// Rcpp::Rcout and Rcpp::Rcerr write through Rprintf, which is an R API call.
// Used off the main thread that risks interleaved output, corrupted console
// state in front ends that assume single-threaded access, and -- if R signals
// during the write -- a longjmp across a thread boundary, which is undefined
// behaviour. The debug tracing in block.cpp / latent.cpp / ngme.cpp sits in
// code the parallel SGD loop reaches from its workers, so give those threads a
// sink that discards. On the main thread these are exactly Rcout / Rcerr, so
// single-threaded runs and everything outside a parallel region are unchanged;
// with several chains you see the output of the one running on the main thread
// rather than all of them interleaved.
inline std::ostream &out() {
  return on_main_thread() ? Rcpp::Rcout : detail::null_stream();
}

inline std::ostream &err() {
  return on_main_thread() ? Rcpp::Rcerr : detail::null_stream();
}

} // namespace ngme_io

#endif // NGME_THREAD_IO_H
