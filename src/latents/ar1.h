#include <Eigen/SparseLU>
#include "../include/solver.h"
#include "../latent.h"
#include "../var.h"

class AR : public Latent {
    SparseLU<SparseMatrix<double> > solver;
    
    double alpha;
public:
    AR(Rcpp::List ar1_in) 
    : Latent(ar1_in)
    {
        // Init operator
        Rcpp::List ope_in = Rcpp::as<Rcpp::List> (ar1_in["operator_in"]); // containing C and G
        ope = new GC(ope_in);

        // solver.analyzePattern(getK());
    }

    void sample_V() {};
    void sample_cond_V() {};

    void setTheta(VectorXd& theta) { 
        ope->update(theta);
    } 

    // double _grad() {
    //     // 1. To compute tr(dK * K^-1)
    //     SparseMatrix<double> I (n_obs, n_obs); 
    //     I.setIdentity();
    //     double g = (dK * solver.solve(I)).eval().diagonal().sum();

    //     // 2. Compute the rest
    //     double rhs = W.transpose() * dK.transpose() * 
    //                 (VectorXd::Constant(n_obs, 1).cwiseQuotient(V).asDiagonal()) * (K * W + (h - V) * mu);

    //     return -(g - rhs);
    // }

};
// std::cout << "K=" << K; std::cout << "dK=" << dK;  // std::cout << "W=" << W;  std::cout << "V=" << V;