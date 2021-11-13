#ifndef NGME_VAR_H
#define NGME_VAR_H

#include <Eigen/Dense>
#include <string>

using Eigen::VectorXd;
using std::string;

class Var {
public:
    VectorXd V;
    string var_type;
    VectorXd theta_V;
    
    VectorXd grad;
    MatrixXd hess;

    virtual void sample_V();
    // { ngme2::rig(n_obs, nu, nu) for AR}

    // sample V given w and Y
    virtual void sample_cond_V();
};

// 
class indepent_IG : Var {

};

#endif
