#ifndef NGME_MODEL_H
#define NGME_MODEL_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Triplet;

class Model {
protected:
    bool computed_{false};
public:
    virtual VectorXd             get_parameter()=0;
    virtual VectorXd             get_stepsizes() {
        return VectorXd::Ones(get_parameter().size());
    };
    virtual void                 set_parameter_and_update(const VectorXd&, bool with_precond)=0;
    // Unified compute entry: perform all heavy computations (Gibbs/RB) here.
    // with_precond=true computes preconditioner caches in addition to gradients.
    virtual void                 compute(bool with_precond=false, double eps=1e-5)=0;
    // Accessors: return the last computed results
    virtual VectorXd             grad()=0;
    // Preconditioner getter.
    virtual MatrixXd             precond()=0;

    virtual int                  get_n_params() const = 0;
    virtual ~Model() = default;
    // Cache state helpers
    bool is_computed() const { return computed_; }
protected:
    void mark_computed() { computed_ = true; }
    void reset_computed() { computed_ = false; }
};


// -----  For testing, f(x, y) = 3x^2 + 2y^2 + x + 3y + 5;

typedef Eigen::Triplet<double> T;

class SomeFun : public Model {
friend class Optimizer;
private:
    Eigen::VectorXd x;
    Eigen::MatrixXd H;
    Eigen::VectorXd last_grad;
    bool grad_valid{false};
    Eigen::MatrixXd last_prec;
    bool prec_valid{false};
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

    void compute(bool with_precond=false, double eps=1e-5) override {
        (void)eps;
        last_grad.resize(2);
        last_grad(0) = 6 * x(0) + 1;
        last_grad(1) = 4 * x(1) + 3;
        grad_valid = true;
        if (with_precond) {
            last_prec = H;
            prec_valid = true;
        }
        mark_computed();
    }

    // gradient accessor
    VectorXd grad() override {
        if (!grad_valid) compute(false);
        return last_grad;
    }

    // hessian accessor
    MatrixXd precond() override {
        if (!prec_valid) compute(true);
        return last_prec;
    }

    void set_parameter_and_update(const VectorXd& x, bool with_precond) override {
        (*this).x = x;
        grad_valid = false;
        prec_valid = false;
        reset_computed();
    }

    VectorXd get_parameter() override {
        return x;
    }

    VectorXd get_stepsizes() override {
        return Eigen::VectorXd::Constant(2, 0);
    }
};


#endif
