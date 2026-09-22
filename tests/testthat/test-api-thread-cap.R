test_that("the OpenMP thread budget caps workers without dropping chains", {
  ctl <- control_opt(n_parallel_chain = 4, max_num_threads = 2)
  expect_equal(ctl$n_parallel_chain, 4)
  expect_equal(ctl$num_threads, c(2, 1))
  ctl <- control_opt(n_parallel_chain = 2, max_num_threads = 5)
  expect_equal(ctl$num_threads, c(2, 2))
})
