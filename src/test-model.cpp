// simple test case

#include <testthat.h>
#include <iostream>
#include "model.h"
#include "optimizer.h"
#include "optimizer.cpp"

context("tests for models") {
    test_that("case 1") {
        expect_true(1 < 2);
    }

}
