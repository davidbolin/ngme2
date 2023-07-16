#ifndef NGME_MODEL_H
#define NGME_MODEL_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Triplet;

class Model {
public:
    virtual VectorXd             get_parameter()=0;
    virtual VectorXd             get_stepsizes() {
        return VectorXd::Ones(get_parameter().size());
    };
    virtual void                 set_parameter(const VectorXd&)=0;
    virtual VectorXd             grad()=0;
    virtual MatrixXd             precond(int strategy=0, double eps=1e-5)=0;

    virtual int                  get_n_params() const = 0;
    virtual ~Model() = default;
};


// -----  For testing, f(x, y) = 3x^2 + 2y^2 + x + 3y + 5;

typedef Eigen::Triplet<double> T;

class SomeFun : public Model {
friend class Optimizer;
private:
    Eigen::VectorXd x;
    Eigen::MatrixXd H;
public:
    SomeFun()
    : x(2), H(2, 2) {
        x << 0, 0;

        // set H
        std::vector<T> tripletList;
        // tripletList.reserve(2);
        // tripletList.push_back(T(0, 0, 6));
        // tripletList.push_back(T(1, 1, 4));
        // H.setFromTriplets(tripletList.begin(), tripletList.end());
        H << 6, 0,
             0, 4;
    }
    SomeFun(VectorXd x0) : x(x0) {}

    int get_n_params() const override {
        return 2;
    }

    // gradient of f(x, y)
    VectorXd grad() override {
        VectorXd g (2);
        g(0) = 6 * x(0) + 1;
        g(1) = 4 * x(1) + 3;
        return g;
    }

    // hessian of f(x, y)
    MatrixXd precond(int strategy, double eps) override {
        return H;
    };

    void set_parameter(const VectorXd& x) {
        (*this).x = x;
    }

    VectorXd get_parameter() override {
        return x;
    }

    VectorXd get_stepsizes() override {
        return Eigen::VectorXd::Constant(2, 0);
    }
};


#endif