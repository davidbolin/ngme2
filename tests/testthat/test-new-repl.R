
test_that("R interface of replicate", {
  n_obs <<- 6; Y <- rnorm(n_obs)

  matern1d <- model_matern(loc = sample(1:10, size=6), mesh = inla.mesh.1d(loc=1:10))

  arr <- model_ar1(1:n_obs)

  repl1 <- c(1,1,2,2,2,2); repl2 <- c(1,1,2,2,2,2)
  formula <- Y ~ f(model=arr, replicate = repl1) +
    f(model=matern1d, replicate = repl2)

  repl <- merge_repls(list(repl1, repl2))

load_all()
  m1 <- ngme(
    Y ~ f(model=arr) + f(model=matern1d),
    data = list(Y=Y),
    control = control_opt(
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



test_that("test split block", {
  { # compute m1
    n_obs <<- 6; Y <- rnorm(n_obs)
    matern1d <- model_matern(loc = sample(1:10, size=6), mesh = inla.mesh.1d(loc=1:10))
    arr <- model_ar1(1:n_obs)

    repl1 <- c(1,1,2,2,2,2);
    repl2 <- c(1,1,2,2,2,2)
    formula <- Y ~ f(model=arr, replicate = repl1) +
      f(model=matern1d, replicate = repl2)

    repl <- merge_repls(list(repl1, repl2))

  load_all()
    m1 <- ngme(
      Y ~ x + f(model=arr) + f(model=matern1d),
      data = list(Y=Y, x=1:n_obs),
      control = control_opt(estimation = F)
    )
  }
  f(mesh = list(mesh))
  6 -> 11 2222
  m1$W_sizes 10+6
  m1$latents

  str(m1)
  m1$X
  split(m1$Y, repl)
  split_matrix(m1$X, repl)

  mat1 <- matrix(1:9, nrow = 3)
  mat1
  split_matrix(mat1, c(1,1,2))

  As <- lapply(m1$latents, function(latent) {
    split_matrix(latent$A, repl)
  })

  load_all()
  res <- split_block(m1, repl)

})

