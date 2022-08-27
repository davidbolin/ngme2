#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <cmath>

using namespace std;
using std::pow;
using namespace Eigen;
void help() {}

const double a = 0.5;

// timer
template <
    class result_t   = std::chrono::milliseconds,
    class clock_t    = std::chrono::steady_clock,
    class duration_t = std::chrono::milliseconds
>
auto since(std::chrono::time_point<clock_t, duration_t> const& start)
{
    return std::chrono::duration_cast<result_t>(clock_t::now() - start);
}

enum class Opt_flag {ope, mu, sigma, var};
enum class Fix_flag {ope, mu, sigma, var, V, W};

class Base {
protected:
    Vector2d vec {1,2};
public:
    Vector2d& getVec() {return vec;}
};

int main(int argc, char const *argv[])
{
    // Vector3d v (1,2,3);
    // VectorXd v2; v2.resize(v.size()+1);
    // v2 << 0, v;

    Base b;
// auto start = std::chrono::steady_clock::now();

// //     Matrix2f v (4,1,2,3);
//     Vector2d k (3,2);
// // v.segment(2, 2) = k;
//     double sigma_eps = 2.0;
//     unsigned n_obs = 100;
// //  g = -1.0 * n_obs / sigma_eps + pow(sigma_eps, -3) * tmp.dot(tmp);
//     double tmp = 100;
    // v = v.array().pow(3);
    // cout<< b.getVec()(1) << endl ;
    // Vector4d v (1,2,3,4);
    // Vector4d vv (1,2,3,4);
    // // v = v.array().log();
    // MatrixXd spm = v.asDiagonal();
    // spm.sparseView();
    // const int i = v.size();
    // // Eigen::DiagonalMatrix<double, i> D (1,2,3,4);
    // // cout << spm;
    // vv(0) += 10;


    VectorXd v1 (1);
    v1.resize(3);
    v1 << 1,2,3;

    // MatrixXd m1 (3, 3);
    SparseMatrix<double> m1 = {3,3};
    m1 = v1.asDiagonal();

    VectorXd v2 (2);
    v2 << 1,2;
    
    // v1.resize(2);
    v1 = v2;

    Matrix3f m;
m << 1, 2, 3,
     4, 5, 6,
     7, 8, 9;

    // cout << m.col(0); // m(Eigen::all, 1) * 5;
    VectorXd m9 (0);
    v1.segment(0, 0) = m9;
    cout << v1;

    return 0;
}


// Eigen::VectorXd cholesky_solver::rMVN(Eigen::VectorXd &mu, Eigen::VectorXd &z)
// {

//   // dest = R.permutationP() * mu;
//   Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> RP;
//   RP.resize(mu.size());
//   RP = R.permutationP();
//   Eigen::VectorXd dest = R.permutationP() * mu;
//   dest = R.matrixL().solve(dest);
//   dest = R.matrixU().solve(dest + z);
//   dest = R.permutationPinv() * dest;
//   return dest;
// }



