#ifndef NGME_LATANT_H
#define NGME_LATANT_H

#include <Eigen/Dense>
#include <string>
#include "operator.h"
#include "var.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class Latent {
friend class Var;
protected:
    double size;
    VectorXd theta_K, theta_m, theta_V;

    // VectorXd m = 0;
    MatrixXd A; // comes from data ()

    Operator K;
    Var var;
public:
    // initFromList(Rlist)
    Latent(MatrixXd A) {}
    
    
    virtual void sample_V();
    // { ngme2::rig(n_obs, nu, nu) for AR}

    // sample V given w and Y
    virtual void sample_cond_V();
};

class AR : public Latent {
friend class Var;

public:
    void sample_V() {
        if (var.var_type == "IG") {
            VectorXd v(size);
            // ngme2::rig(n_obs, nu, nu) for AR

            // v <- rGIG(size, nu, nu)
            var.V = v;
        }
        
        
    }
};


#endif