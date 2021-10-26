// simple test case

#include <testthat.h>
#include <iostream>
#include "model.h"
#include "optimizer.h"
#include "optimizer.cpp"


context("tests for optimizer") {

    test_that("convergence for gd") {
        model f;
        optimizer opt;
        Eigen::Vector2d x0 (10, 10);

        VectorXd out = opt.sgd(f, x0, 0.05, 0.01, false);

        // optimum should be (-1/6, -3/4)
        expect_true(out(0) + 1.0/6 < 0.001);
        expect_true(out(1) + 3.0/4 < 0.001);
    }


    // Use Hessian as preconditoner
    test_that("convergence using hessian") {
        model f;
        optimizer opt;
        Eigen::Vector2d x0(10, 10);

        VectorXd out = opt.sgd(f, x0, 0.05, 0.01, true);

        expect_true(out(0) + 1.0/6 < 0.001);
        expect_true(out(1) + 3.0/4 < 0.001);
    }

}
