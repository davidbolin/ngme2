// simple test case

#include <testthat.h>
#include <iostream>
#include "model.h"
#include "optimizer.h"
#include "optimizer.cpp"


context("tests for optimizer") {

    test_that("simple test for gd") {
        model f;
        optimizer opt;
        Eigen::Vector2d x0 (5,5);
        VectorXd out = opt.gd(f, x0, 0.05, 0.01, false);

        // close to 0
        expect_true(out.norm() < 1);
    }

}
