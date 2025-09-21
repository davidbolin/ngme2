#ifdef USEMKL
  #define EIGEN_USE_MKL_ALL
#endif

#define COMPLEX R_COMPLEX
#include <Rcpp.h>
#include <RcppEigen.h>
#undef COMPLEX

#include "block.h"

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
double compute_log_like_cpp(const Rcpp::List& R_ngme) {
    const unsigned long base_seed = 0;
    double total_log_like = 0.0;

    Rcpp::List replicates = R_ngme["replicates"];
    const int n_repl = replicates.size();

    for (int i = 0; i < n_repl; ++i) {
        Rcpp::List ngme_repl = replicates[i];
        BlockModel block(ngme_repl, base_seed + static_cast<unsigned long>(i));

        if (!block.is_all_gaussian()) {
            return 0.0;
        }

        total_log_like += block.log_likelihood();
    }

    return total_log_like;
}
