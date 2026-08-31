test_that("prior constructors validate and carry target", {
  p <- prior_normal(mean = 1, sd = 2, target = "field")
  expect_s3_class(p, "ngme_prior_spec")
  expect_equal(p$dist, "normal")
  expect_equal(unname(p$hyper), c(1, 2))
  expect_equal(p$target, "field")

  expect_error(prior_normal(sd = -1))
  expect_error(prior_pc_sd(u = 0, alpha = 0.1))
  expect_error(prior_pc_sd(u = 1, alpha = 1.2))
  expect_error(prior_half_cauchy(scale = 0))
  expect_error(prior_inv_exponential(lambda = 0))
  expect_error(prior_inv_exponential(lower = -1))
})

test_that("prior_inv_exp is an alias of prior_inv_exponential", {
  p1 <- prior_inv_exponential(lambda = 1.2, lower = 0.3, target = "coef")
  p2 <- prior_inv_exp(lambda = 1.2, lower = 0.3, target = "coef")
  expect_equal(p2, p1)
})

test_that("ngme_noise compiles named priors to internal format", {
  p <- priors(
    mu = prior_normal(0, 2, target = "coef"),
    sigma = prior_pc_sd(1, 0.1, target = "field"),
    nu = prior_inv_exponential(2, lower = 0.5, target = "coef")
  )

  noise <- ngme_noise(
    noise_type = "nig",
    theta_mu = 0,
    theta_sigma = 0,
    theta_nu = 0,
    B_mu = matrix(1, 1, 1),
    B_sigma = matrix(1, 1, 1),
    B_nu = matrix(1, 1, 1),
    prior = p
  )

  expect_equal(noise$prior_mu$type, "normal")
  expect_equal(noise$prior_mu$target, "coef")
  expect_equal(noise$prior_mu$param[2], 1 / 4)

  expect_equal(noise$prior_sigma$type, "pc.sd")
  expect_equal(noise$prior_sigma$target, "field")
  expect_equal(noise$prior_sigma$param[1], -log(0.1) / 1)

  expect_equal(noise$prior_nu$type, "inv.exponential")
  expect_equal(noise$prior_nu$target, "coef")
  expect_equal(noise$prior_nu$param[1], 2)
  expect_equal(noise$prior_nu$param[2], 0.5)
})

test_that("f compiles operator priors from new API", {
  fit <- f(
    map = 1:10,
    model = ar1(),
    noise = noise_normal(),
    prior = priors(
      theta = prior_normal(0, 2),
      rho = prior_half_cauchy(1)
    )
  )

  expect_true(inherits(fit, "ngme_model"))
  expect_equal(length(fit$prior_theta_K), length(fit$operator$theta_K))
  expect_equal(fit$prior_theta_K[[1]]$type, "half.cauchy")
  expect_equal(fit$prior_theta_K[[1]]$target, "coef")

  expect_error(
    f(
      map = 1:10,
      model = ar1(),
      noise = noise_normal(),
      prior = priors(bad_name = prior_normal(0, 1))
    ),
    "Unknown operator prior names"
  )
})

test_that("prior_pc_nu calibrates lambda from P(1/nu > U) = alpha", {
  # Cabral, Bolin & Rue (2023), Sec 3.3: lambda = -log(alpha) / U.
  expect_equal(prior_pc_nu(U = 1, alpha = 0.05)$hyper[["lambda"]], -log(0.05))
  expect_equal(prior_pc_nu(U = 2, alpha = 0.01)$hyper[["lambda"]], -log(0.01) / 2)
  expect_equal(prior_pc_nu()$dist, "inv.exponential")
  expect_error(prior_pc_nu(U = 0))
  expect_error(prior_pc_nu(alpha = 0))
  expect_error(prior_pc_nu(alpha = 1))
})

test_that("prior_pc_nu defaults are settable through options", {
  withr::with_options(
    list(ngme2.pc_nu_U = 2, ngme2.pc_nu_alpha = 0.01),
    expect_equal(prior_pc_nu()$hyper[["lambda"]], -log(0.01) / 2)
  )
  # and revert once the option is gone
  expect_equal(prior_pc_nu()$hyper[["lambda"]], -log(0.01) / 2.5)
})

test_that("f sets the default nu prior from the PC calibration, not from h", {
  fit <- f(
    map = 1:20,
    model = ar1(),
    noise = noise_nig()
  )
  expect_equal(fit$noise$prior_nu$type, "inv.exponential")
  expect_equal(fit$noise$prior_nu$target, "coef")
  expect_equal(fit$noise$prior_nu$param[1], -log(0.01) / 2.5, tolerance = 1e-12)
  expect_equal(fit$noise$prior_nu$param[2], 0)
})

test_that("the default nu prior does not depend on the size of the domain", {
  # nu is a property of the process: observing more of it must not move the
  # prior, so lambda has to be free of any domain-extent term.
  lams <- vapply(c(20, 40, 80, 160), function(n) {
    f(map = 1:n, model = ar1(), noise = noise_nig())$noise$prior_nu$param[1]
  }, numeric(1))
  expect_equal(max(lams) / min(lams), 1, tolerance = 1e-12)
})

test_that("the default nu prior is invariant to mesh refinement", {
  skip_on_cran()
  skip_if_not_installed("fmesher")
  # The same domain at three resolutions must give (nearly) the same lambda.
  loc <- withr::with_seed(3, matrix(runif(600 * 2), 600, 2))
  d <- data.frame(c1 = loc[, 1], c2 = loc[, 2])
  lams <- vapply(c(0.13, 0.075, 0.045), function(ct) {
    msh <- fmesher::fm_mesh_2d(loc = loc, cutoff = ct,
                               max.edge = c(ct * 2.5, ct * 8))
    fit <- f(~ c1 + c2, model = matern(mesh = msh), noise = noise_nig(), data = d)
    fit$noise$prior_nu$param[1]
  }, numeric(1))
  # node counts differ several-fold, so this is a real test of the invariance
  expect_lt(max(lams) / min(lams), 1.1)
})

test_that("f keeps explicit nu prior and does not override it", {
  fit <- f(
    map = 1:20,
    model = ar1(),
    noise = noise_nig(
      prior = priors(
        nu = prior_normal(0, 2)
      )
    )
  )

  expect_equal(fit$noise$prior_nu$type, "normal")
  expect_equal(fit$noise$prior_nu$param[2], 1 / 4)
})

test_that("f keeps normal default prior for non-stationary nu", {
  B_nu <- cbind(1, seq(0, 1, length.out = 20))
  fit <- f(
    map = 1:20,
    model = ar1(),
    noise = noise_nig(
      B_nu = B_nu,
      theta_nu = c(0, 0)
    )
  )

  expect_equal(fit$noise$prior_nu$type, "normal")
  expect_equal(fit$noise$prior_nu$param[1], 0)
  expect_equal(fit$noise$prior_nu$param[2], 0.01)
})

test_that("beta prior compilation supports global and named forms", {
  b <- compile_beta_priors(
    priors(
      beta = prior_normal(0, 2),
      x = prior_half_cauchy(1),
      feff_3 = prior_normal(1, 0.5)
    ),
    c("(Intercept)", "x", "z")
  )

  expect_equal(length(b), 3)
  expect_equal(b[[1]]$type, "normal")
  expect_equal(b[[1]]$param[2], 1 / 4)
  expect_equal(b[[2]]$type, "half.cauchy")
  expect_equal(b[[3]]$type, "normal")
  expect_equal(b[[3]]$param[1], 1)

  expect_error(
    compile_beta_priors(priors(unknown = prior_normal(0, 1)), c("a", "b")),
    "Unknown beta prior names"
  )
  expect_error(
    compile_beta_priors(priors(beta = prior_normal(0, 1, target = "field")), c("a")),
    "target='coef' only"
  )
})

test_that("default beta prior leaves intercept unconstrained", {
  b <- compile_beta_priors(NULL, c("(Intercept)", "x", "(Intercept)_u_wind", "z"))

  expect_equal(length(b), 4)
  expect_equal(b[[1]]$type, "none")
  expect_equal(b[[1]]$target, "coef")
  expect_equal(length(b[[1]]$param), 0)
  expect_equal(b[[2]]$type, "normal")
  expect_equal(b[[3]]$type, "none")
  expect_equal(b[[4]]$type, "normal")
})

test_that("ngme accepts prior_beta and stores compiled priors per replicate", {
  dat <- data.frame(y = rnorm(20), x = rnorm(20), id = 1:20)

  fit <- ngme(
    y ~ x + f(id, model = ar1(), noise = noise_normal()),
    data = dat,
    family = noise_normal(),
    prior_beta = priors(
      beta = prior_normal(0, 2),
      x = prior_half_cauchy(1)
    ),
    control_opt = control_opt(estimation = FALSE, standardize_fixed = FALSE)
  )

  pb <- fit$replicates[[1]]$prior_beta
  expect_equal(length(pb), ncol(fit$replicates[[1]]$X))
  expect_equal(pb[[2]]$type, "half.cauchy")
  expect_equal(pb[[1]]$type, "normal")
})
