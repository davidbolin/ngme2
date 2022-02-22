#include <Eigen/SparseLU>
#include "../include/solver.h"
#include "../latent.h"
#include "../var.h"

class AR : public Latent {
    SparseLU<SparseMatrix<double> > solver_K;
    
public:
    AR(Rcpp::List ar1_in) 
    : Latent(ar1_in)
    {
        // Init operator
        Rcpp::List ope_in = Rcpp::as<Rcpp::List> (ar1_in["operator_in"]); // containing C and G
        ope = new GC(ope_in);

        solver_K.analyzePattern(getK());
    }

    void sample_V() {};
    void sample_cond_V() {};

    double _grad_kappa() {
        SparseMatrix<double> K = getK();
        SparseMatrix<double> dK = get_dK();
        VectorXd V = getV();
        solver_K.factorize(K);

        // 1. To compute tr(dK * K^-1)
        SparseMatrix<double> I (n_reg, n_reg); 
        I.setIdentity();
        double lhs = (dK * solver_K.solve(I)).eval().diagonal().sum();

        // 2. Compute the rest
        double rhs = W.transpose() * dK.transpose() * 
                    (VectorXd::Constant(n_reg, 1).cwiseQuotient(V).asDiagonal()) * (K * W + (h - V) * mu);

        return -(lhs - rhs) / n_reg;
    }

};
