#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <cmath>

using std::pow;
using namespace Eigen;
using namespace std;
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

class Base {
protected:
    Vector2d vec {1,2};
public:
    Vector2d& getVec() {return vec;}
};

int main(int argc, char const *argv[])
{
    Vector3d v (1,2,3);

// auto start = std::chrono::steady_clock::now();

// //     Matrix2f v (4,1,2,3);
//     Vector2d k (3,2);
// // v.segment(2, 2) = k;
//     double sigma_eps = 2.0;
//     unsigned n_obs = 100;
// //  g = -1.0 * n_obs / sigma_eps + pow(sigma_eps, -3) * tmp.dot(tmp);
//     double tmp = 100;
    v = v.array().pow(3);
    cout<< v << endl ;

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