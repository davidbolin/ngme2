
test_that("R interface of replicate", {
  load_all()
  library(INLA)

  f(model=ar1(1:3))
  ar <- ar1(1:5)
  sub_fmodel(ar, 1:3)$map

  n_obs <<- 6; Y <- rnorm(n_obs)

  matern1d <- model_matern(map = sample(1:10, size=6), mesh = inla.mesh.1d(loc=1:10))

  arr <- model_ar1(1:n_obs)

  repl1 <- c(1,1,2,2,2,2); repl2 <- c(1,1,2,2,2,2)
  formula <- Y ~ f(model=arr, replicate = repl1) +
    f(model=matern1d, replicate = repl2)

  load_all()
  m1 <- ngme(
    Y ~ f(model=arr, replicate = repl1) + f(model=matern1d, replicate = repl2),
    data = list(Y=Y),
    control_opt = control_opt(
      estimation = F
    )
  )
  m1

  # latent:
  matern1d$A
  A1s <- split(ar1$A, repl)
  A2s <- split(matern1d$A, repl)

  # block: divide Y, X,
  Ys <- split(Y, repl); Ys
  # Xs <- split(X, repl), Xs
  # W_szs, V_szs

  m_noise <- m_noises
})

test_that("test create ngme block", {
  { # compute m1
  library(INLA)
  load_all()
    n_obs <<- 6; Y <- rnorm(n_obs)
    matern1d <- model_matern(loc = sample(1:10, size=6), mesh = inla.mesh.1d(loc=1:10))
    arr <- model_ar1(1:n_obs)

    repl1 <- c(1,1,2,2,2,2);
    repl2 <- c(1,1,2,2,2,2)
    formula <- Y ~ f(model=arr, replicate = repl1) +
      f(model=matern1d, replicate = repl2)

    repl <- merge_repls(list(repl1, repl2))

    m1 <- ngme(
      Y ~ x + f(model=arr) + f(model=matern1d),
      data = list(Y=Y, x=1:n_obs),
      control_opt = control_opt(estimation = F)
    )
  }

  load_all()
  out <- ngme(
    Y ~ x + f(model=arr, replicate = repl1) + f(model=matern1d, replicate = repl2),
    data = list(Y=Y, x=1:n_obs),
    control_opt = control_opt(estimation = F)
  )
  out$latents[[1]]$replicate

})

# test_that("test create ngme replicate")
y <- rnorm(100)

inla(
  y~ 1 + x,
  data = list(y=y, x=1:10)
)


# Y     1    2   3
# Year  201 202 203

test_that("basic ar1 case with different length", {
  n_obs1  <- 500
  n_obs2 <- 100
  alpha <- 0.75
  mu <- -3; sigma <- 2.3; nu <- 2; sigma_eps <- 0.8
  ar_1 <- model_ar1(1:n_obs1, alpha=alpha, noise=noise_nig(mu=mu, sigma=sigma, nu=nu))
  ar_2 <- model_ar1(1:n_obs2, alpha=alpha, noise=noise_nig(mu=mu, sigma=sigma, nu=nu))

  W1 <- simulate(ar_1)
  W2 <- simulate(ar_2)
  Y <- c(W1, W2) + rnorm(n_obs1 + n_obs2, sd = sigma_eps)
  beta <- c(3,1,-1)
  x1 <- rnorm(n_obs1 + n_obs2)
  x2 <- rexp(n_obs1 + n_obs2)
  Y <- beta[1] + beta[2] * x1 + beta[3] * x2 + Y
  # Y <- beta[1] + Y

  load_all()
  out <- ngme(
    Y ~ x1 + x2 + f(time,
      model="ar1",
      name="ar",
      alpha = -0.5,
      replicate = c(rep(1, n_obs1), rep(2, n_obs2)),
      noise=noise_nig(),
      control=control_f(
        numer_grad = T
      ),
      debug = FALSE
    ),
    data = list(Y = Y, time=c(1:n_obs1, 1:n_obs2), x1=x1, x2=x2),
    control_opt = control_opt(
      estimation = T,
      iterations = 50,
      n_parallel_chain = 2,
      print_check_info = FALSE,
      verbose = F
    ),
    debug = FALSE
  )
  out
  load_all()
cross_validation(out)
  traceplot(out)
  traceplot(out, "ar")
  plot(simulate(out[[1]]$latents[["ar"]]), type="l")
  plot(Y, type="l")
  prds <- predict(out, loc=list(ar=501:600))$mean

  traceplot(out, "ar")
  traceplot(out, "mn")
  plot(attr(W1, "noise"), out[[1]]$latents[[1]]$noise)

  # test on predict function
  predict(out,
    loc = list(c(1,2,3,4)),
    data = list(cbind(1, rnorm(4), rexp(4)))
  )

  cross_validation(out)

})
