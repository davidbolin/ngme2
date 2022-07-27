#include <Eigen/Sparse>
// #include <bench/BenchTimer.h>
#include <iostream>
#include <vector>
using namespace Eigen;
using namespace std;
// Use RowMajor to make use of multi-threading
typedef SparseMatrix<double> SpMat;
typedef Triplet<double> T;
// Assemble sparse matrix from
// https://eigen.tuxfamily.org/dox/TutorialSparse_example_details.html
void insertCoefficient(int id, int i, int j, double w, vector<T>& coeffs,
                       VectorXd& b, const VectorXd& boundary)
{
  int n = int(boundary.size());
  int id1 = i+j*n;
        if(i==-1 || i==n) b(id) -= w * boundary(j); // constrained coefficient
  else  if(j==-1 || j==n) b(id) -= w * boundary(i); // constrained coefficient
  else  coeffs.push_back(T(id,id1,w));              // unknown coefficient
}
void buildProblem(vector<T>& coefficients, VectorXd& b, int n)
{
  b.setZero();
  ArrayXd boundary = ArrayXd::LinSpaced(n, 0,M_PI).sin().pow(2);
  for(int j=0; j<n; ++j)
  {
    for(int i=0; i<n; ++i)
    {
      int id = i+j*n;
      insertCoefficient(id, i-1,j, -1, coefficients, b, boundary);
      insertCoefficient(id, i+1,j, -1, coefficients, b, boundary);
      insertCoefficient(id, i,j-1, -1, coefficients, b, boundary);
      insertCoefficient(id, i,j+1, -1, coefficients, b, boundary);
      insertCoefficient(id, i,j,    4, coefficients, b, boundary);
    }
  }
}
int main()
{
  int n = 3000;  // size of the image
  int m = n*n;  // number of unknowns (=number of pixels)
  // Assembly:
  vector<T> coefficients;          // list of non-zeros coefficients
  VectorXd b(m);                   // the right hand side-vector resulting from the constraints
  buildProblem(coefficients, b, n);
  SpMat A(m,m);
  A.setFromTriplets(coefficients.begin(), coefficients.end());
  // Solving:
  // Use ConjugateGradient with Lower|Upper as the UpLo template parameter to make use of multi-threading
//   BenchTimer t;
//   t.reset(); t.start();
  ConjugateGradient<SpMat, Lower|Upper> solver(A);
  VectorXd x = solver.solve(b);         // use the factorization to solve for the given right hand side
//   t.stop();
//   cout << "Real time: " << t.value(1) << endl; // 0=CPU_TIMER, 1=REAL_TIMER
  return 0;
}





