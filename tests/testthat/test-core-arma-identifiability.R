# A Gaussian latent ARMA(p, q) observed under Gaussian measurement noise is
# second-order equivalent to an ARMA(p, max(p, q)), so its likelihood has a
# one-dimensional flat ridge whenever q >= p. ngme() warns in that case.

sim_data <- function(n = 200) {
  set.seed(7)
  data.frame(x = 1:n, y = as.numeric(arima.sim(list(ar = 0.5), n)) + rnorm(n))
}

fit_quick <- function(model, family = "normal") {
  ngme(y ~ 0 + f(x, model = model, noise = noise_normal()),
       family = family, data = sim_data(),
       control_opt = control_opt(seed = 1, burnin = 5, iterations = 10,
                                 n_parallel_chain = 1,
                                 warn_no_convergence = FALSE))
}

test_that("non-identifiable Gaussian ARMA(p, q) with q >= p warns", {
  for (ord in list(c(1, 1), c(1, 2), c(2, 2))) {
    expect_warning(
      fit_quick(arma(ar_order = ord[1], ma_order = ord[2])),
      "not identifiable",
      info = paste0("ARMA(", ord[1], ",", ord[2], ")")
    )
  }
})

test_that("identifiable ARMA configurations do not warn", {
  # p > q is identified with free measurement noise
  expect_no_warning(fit_quick(arma(ar_order = 2, ma_order = 1)))
  expect_no_warning(fit_quick(arma(ar_order = 1, ma_order = 0)))

  # fixing the measurement scale removes the ridge
  expect_no_warning(fit_quick(
    arma(ar_order = 1, ma_order = 1),
    family = noise_normal(theta_sigma = 0, fix_theta_sigma = TRUE)
  ))

  # fixing the MA coefficients removes it too
  expect_no_warning(fit_quick(arma(ar_order = 1, ma_order = 1, fix_ma = TRUE)))
})

test_that("non-Gaussian noise identifies the model through higher moments", {
  # non-Gaussian latent noise
  expect_no_warning(
    ngme(y ~ 0 + f(x, model = arma(ar_order = 1, ma_order = 1),
                   noise = noise_nig()),
         data = sim_data(),
         control_opt = control_opt(seed = 1, burnin = 5, iterations = 10,
                                   n_parallel_chain = 1,
                                   warn_no_convergence = FALSE))
  )
  # non-Gaussian measurement noise
  expect_no_warning(fit_quick(arma(ar_order = 1, ma_order = 1), family = "nig"))
})
