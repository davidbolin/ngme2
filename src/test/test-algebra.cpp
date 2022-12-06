/*
 * This file uses the Catch unit testing library, alongside
 * testthat's simple bindings, to test a C++ function.
 *
 * For your own packages, ensure that your test files are
 * placed within the `src/` folder, and that you include
 * `LinkingTo: testthat` within your DESCRIPTION file.
 */

// All test files should include the <testthat.h>
// header file.
#include <testthat.h>

#include <Eigen/Sparse>
#include "../include/solver.h"

using namespace Eigen;

context("Test generating MVN") {
  cholesky_solver solver;

  test_that("two plus two equals four") {
    SparseMatrix<double> m1(4, 4);
    m1.coeffRef(0,0) = 1;
    m1.coeffRef(1,1) = 3.;
    m1.coeffRef(2,2) = 4.;
    // m1.coeffRef(2,3) = 5.;
    // m1.coeffRef(3,2) = 6.;
    m1.coeffRef(3,3) = 7.;
    m1.makeCompressed();

    VectorXd mu (4);
    mu << 1,1,1,1;

    std::cout << m1;
    solver.analyze(m1);
    std::cout << solver.rMVN(mu, mu);
  }
}


