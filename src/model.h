#ifndef NGME_MODEL_H
#define NGME_MODEL_H

#include <Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class Model {
public:
    virtual VectorXd& get_parameter()=0;
    virtual void      set_parameter(const VectorXd&)=0;

    virtual MatrixXd& precond()=0;
    virtual VectorXd& grad()=0;
};





// -----  For testing, f(x, y) = 3x^2 + 2y^2 + x + 3y + 5;
class SomeFun : public Model {
friend class Optimizer;
private:
    Eigen::VectorXd x;
    Eigen::VectorXd g;
    Eigen::MatrixXd H;
public:
    SomeFun()            
    : x(2), g(2), H(2, 2) { 
        x << 0, 0; 
        H << 6, 0,
             0, 4;
    }

    SomeFun(VectorXd x0) : x(x0) {}

    // gradient of f(x, y)
    VectorXd& grad() {
        g(0) = 6 * x(0) + 1;
        g(1) = 4 * x(1) + 3;
        return g;
    }

    // hessian of f(x, y)
    MatrixXd& precond() {
        return H;
    };

    void set_parameter(const VectorXd& x) {
        (*this).x = x;
    }

    VectorXd& get_parameter() {
        return x;
    }
};


#endif