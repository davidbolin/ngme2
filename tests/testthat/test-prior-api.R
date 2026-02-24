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

test_that("f sets default nu prior to inv.exponential using median(h)", {
  fit <- f(
    map = 1:20,
    model = ar1(),
    noise = noise_nig()
  )

  expect_equal(fit$noise$prior_nu$type, "inv.exponential")
  expect_equal(fit$noise$prior_nu$target, "coef")
  expect_equal(fit$noise$prior_nu$param[1], log(2), tolerance = 1e-12)
  expect_equal(fit$noise$prior_nu$param[2], 0)
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
    control_opt = control_opt(estimation = FALSE)
  )

  pb <- fit$replicates[[1]]$prior_beta
  expect_equal(length(pb), ncol(fit$replicates[[1]]$X))
  expect_equal(pb[[2]]$type, "half.cauchy")
  expect_equal(pb[[1]]$type, "normal")
})
