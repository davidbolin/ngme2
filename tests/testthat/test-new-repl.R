
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