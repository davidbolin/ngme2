#include <Eigen/SparseLU>
#include "../include/solver.h"
#include "../latent.h"
#include "../var.h"
#include <cmath>

using Eigen::MatrixXd;
using std::exp;
using std::log;
using std::pow;

class Matern_ns {
public:
    // parameter
    int alpha, n_mesh; // 2 or 4
    
    VectorXd paras;  // theta = (1 para1 para2 ...)
    MatrixXd Btau, Bkappa;
    SparseMatrix<double> C, G;
    SparseMatrix<double> Tau, Kappa;

    void update(VectorXd paras) { // update K, dK 
        VectorXd theta;
            theta.resize(paras.size() + 1);
            theta << 1, paras;
        if (alpha==2) {
            Tau = (Btau * theta).diagonal();
            // Kappa
            // K = ...
            
        } else if (alpha==4) {

        }
        
    }

    
    // const VectorXd Latent::getGrad() {}
    // void Latent::setTheta(const VectorXd& theta) {}
    

};
