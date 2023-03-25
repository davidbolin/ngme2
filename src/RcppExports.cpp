// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// estimate_cpp
Rcpp::List estimate_cpp(const Rcpp::List& list_replicates, const Rcpp::List& control_opt);
RcppExport SEXP _ngme2_estimate_cpp(SEXP list_replicatesSEXP, SEXP control_optSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type list_replicates(list_replicatesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type control_opt(control_optSEXP);
    rcpp_result_gen = Rcpp::wrap(estimate_cpp(list_replicates, control_opt));
    return rcpp_result_gen;
END_RCPP
}
// sampling_cpp
Rcpp::List sampling_cpp(const Rcpp::List& ngme_replicate, int n, bool posterior, unsigned long seed);
RcppExport SEXP _ngme2_sampling_cpp(SEXP ngme_replicateSEXP, SEXP nSEXP, SEXP posteriorSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type ngme_replicate(ngme_replicateSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< bool >::type posterior(posteriorSEXP);
    Rcpp::traits::input_parameter< unsigned long >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(sampling_cpp(ngme_replicate, n, posterior, seed));
    return rcpp_result_gen;
END_RCPP
}
// rGIG_cpp
Eigen::VectorXd rGIG_cpp(Eigen::VectorXd p, Eigen::VectorXd a, Eigen::VectorXd b, unsigned long seed);
RcppExport SEXP _ngme2_rGIG_cpp(SEXP pSEXP, SEXP aSEXP, SEXP bSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type p(pSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type a(aSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type b(bSEXP);
    Rcpp::traits::input_parameter< unsigned long >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(rGIG_cpp(p, a, b, seed));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ngme2_estimate_cpp", (DL_FUNC) &_ngme2_estimate_cpp, 2},
    {"_ngme2_sampling_cpp", (DL_FUNC) &_ngme2_sampling_cpp, 4},
    {"_ngme2_rGIG_cpp", (DL_FUNC) &_ngme2_rGIG_cpp, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_ngme2(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
