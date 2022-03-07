#include <Eigen/SparseLU>
#include "../include/solver.h"
#include "../latent.h"
#include "../var.h"
#include <cmath>

using std::exp;
using std::log;
using std::pow;

class AR : public Latent {
    
    
public:
    AR(Rcpp::List ar1_in) 
    : Latent(ar1_in)
    {
        // Init operator
        Rcpp::List ope_in = Rcpp::as<Rcpp::List> (ar1_in["operator_in"]); // containing C and G
        ope = new GC(ope_in);

        solver_K.init(n_reg, 0,0,0);
        solver_K.analyze(getK());
    }

    void sample_cond_V() {
        var->sample_cond_V(getK(), W, h, mu, sigma);
    };
    

    double th2a(double th) const {
        return (-1 + 2*exp(th) / (1+exp(th)));
    }
    double a2th(double a) const {
        return (log((-1-a)/(-1+a)));
    }

    double get_theta_kappa() const {
        double a = getKappa();
        return a2th(a);
    }

    void set_theta_kappa(double v) {
        double a = th2a(v);
        setKappa(a);
    }

    // grad_kappa * dkappa/dtheta
    double grad_theta_kappa() {
        double a = getKappa();
        double th = a2th(a);
        double dth = 2 * exp(th) / pow(1+exp(th), 2);
        return _grad_kappa() * dth;
    }

};
