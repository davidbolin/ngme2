#ifndef NGME_MODEL_H
#define NGME_MODEL_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
using Eigen::SparseMatrix;
using Eigen::VectorXd;
using Eigen::Triplet;

class Model {
public:
    virtual VectorXd             get_parameter() const=0;
    virtual VectorXd             get_stepsizes() const=0;
    virtual void                 set_parameter(const VectorXd&)=0;
    virtual VectorXd             grad()=0;
    virtual SparseMatrix<double> precond() const=0;
};


// -----  For testing, f(x, y) = 3x^2 + 2y^2 + x + 3y + 5;

typedef Eigen::Triplet<double> T;

class SomeFun : public Model {
friend class Optimizer;
private:
    Eigen::VectorXd x;
    Eigen::SparseMatrix<double> H;
public:
    SomeFun()            
    : x(2), H(2, 2) { 
        x << 0, 0; 

        // set H
        std::vector<T> tripletList;
        tripletList.reserve(2);
        tripletList.push_back(T(0, 0, 6));
        tripletList.push_back(T(1, 1, 4));
        H.setFromTriplets(tripletList.begin(), tripletList.end());
        // H << 6, 0,
        //      0, 4;
    }

    SomeFun(VectorXd x0) : x(x0) {}

    // gradient of f(x, y)
    VectorXd grad() {
        VectorXd g (2);
        g(0) = 6 * x(0) + 1;
        g(1) = 4 * x(1) + 3;
        return g;
    }

    // hessian of f(x, y)
    SparseMatrix<double> precond() const {
        return H;
    };

    void set_parameter(const VectorXd& x) {
        (*this).x = x;
    }

    VectorXd get_parameter() const {
        return x;
    }

    VectorXd get_stepsizes() const {
        return Eigen::VectorXd::Constant(2, 0);
    }
};


#endif