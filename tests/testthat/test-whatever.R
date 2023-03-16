
test_that("test names", {
  expect_equal(f(1:10, name="ar")$name, "ar")

  bm <- ngme(
    Y ~ f(x=1:10, name="f1", model="ar1") +
      f(x=1:10, name="f2", model="ar1"),
    data = data.frame(Y = rnorm(10)),
    control = control_opt(estimation = FALSE)
  )

  expect_equal(names(bm$latents), c("f1", "f2"))
})

# run manually
test_that("predict traceplot", {
  n_obs <<- 500
  x1 <- rexp(n_obs); x2 <- rnorm(n_obs)
  beta <- c(-2, 4, 1)
  alpha <- 0.75
  mu <- 1.5; sigma <- 2.3; nu <- 2; sigma_eps <- 0.8

  ar1 <- model_ar1(1:n_obs, alpha=alpha, noise=noise_nig(mu=mu, sigma=sigma, nu=nu))
  W <- simulate(ar1)
  Y <- beta[1] + x1 * beta[2] + x2 * beta[3] + W + rnorm(n_obs, sd = sigma_eps)

  # make 1/10 observation NA
  idx_NA <- sample(1:n_obs, size=n_obs/10)
  Y2 <- Y; Y2[idx_NA] <- NA

  # first we test estimation
  out <- ngme(
    Y2 ~ x1 + x2
      + f(1:n_obs, model="ar1", noise = noise_nig())
      + f(1:n_obs, model="ar1", noise = noise_normal(), name="f2"),
    family = noise_normal(),
    control = control_opt(
      n_parallel_chain = 4,
      iteration = 20,
      print_check_info = FALSE
   ),
   seed = 100,
   data = data.frame(Y2 = Y2, x1 = x1, x2 = x2)
  )
  out

  expect_no_error(traceplot(out, name = "field1"))
  expect_no_error(traceplot(out, name = "f2"))
})


# tests on posterior sampling
test_that("test posterior sampling and model_validation()", {
# load_all()
  n_obs <<- 20
  alpha <- 0.3; mu = -3; sigma=2; nu=2; sigma_eps <- 0.5
  my_ar <<- model_ar1(1:n_obs, alpha=alpha, noise=noise_nig(mu=mu, sigma=sigma, nu=nu))
  W <- simulate(my_ar)
  Y <- W + rnorm(n=length(W), sd=sigma_eps)

  out <- ngme(Y ~ 0 + f(model=my_ar),
    data=list(Y=Y),
    control=control_opt(print_check_info = FALSE)
  )
  # traceplot(out, "field1")

  expect_no_error(samples <- sampling_cpp(out, 10, posterior = TRUE, seed=seed))
  # samples[["AW"]]
})


test_that("modify_ngme_with_idx_NA", {
  # load_all()
  X2 = c(3,4,5,9,2)
  spde1 <<- model_matern(mesh = 1:10, loc=X2)
  mm <- ngme(
    YY ~ 1 + X1 + f(X2, model="rw1") + f(model=spde1),
    data = list(
      YY=c(0,2,3,4,5),
      X1=c(3,4,5,6,7),
      X2=X2
    ),
    control=control_opt(iterations = 10, print_check_info = FALSE)
  )
  mm
  new_model <- modify_ngme_with_idx_NA(mm, idx_NA = 3)
  str(new_model$noise)
  # new_model$latents[[2]]$A_pred

  # new_model$latents[[1]]$A_pred
  expect_true(length(new_model$Y) == 4)
  # new_model$latents[[1]]
})

test_that("test lpo CV", {
  # load_all()
  n_obs <<- 100
  ar_mod <- f(1:n_obs, model="ar1", noise=noise_nig())
  yy <- simulate(ar_mod) + rnorm(n_obs, sd=0.5)
  ng_100 <- ngme(
    yy~0+f(1:n_obs, model="ar1", noise=noise_nig()),
    data = list(yy=yy),
    control = control_opt(iterations = 100)
  )
  ng_1000 <- ngme(
    yy~0+f(1:n_obs, model="ar1", noise=noise_nig()),
    data = list(yy=yy),
    control = control_opt(iterations = 1000)
  )
  cross_validation(ng_100, type="lpo", times=10)
  cross_validation(ng_1000, type="lpo", times=10)

  expect_true(TRUE)
})


test_that("merge_reps works", {
  repls <- c("1 1 2 2 2 3 4 5 5",
             "1 2 3 1 5 6 4 1 2",
             "1 1 1 1 1 2 2 3 3",
             "1 1 1 1 1 2 2 1 1")
  repls <- lapply(repls, function(x) as.numeric(strsplit(x, " ")[[1]]))
  repls

  # equal to 1 1 1 1 1 2 2 1 1
  expect_true(all(merge_repls(repls)
    == as.numeric(strsplit("1 1 1 1 1 2 2 1 1", " ")[[1]])))
})
