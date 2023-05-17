# Y <- c(1:5, 1:5)
# group <- c(rep(1, 5), rep(2, 5))
# Y ~ Fixeff(~1+x, group=1) + matern_nd(..., group=group)

# test on random effects
test_that("the R interface", {
  load_all()
  f(1:10, model="ar1")
  m = f(~1+x, data=data.frame(x=c(1,2,3)), theta_K = c(2,1,1),
     effect_type="normal", eval=T)
  m
  # 1 random effect + 1 latent model
  fm0 <- Y ~ f(~1+x, effect_type="normal") +
          f(x, model="ar1")

  # 2 random effects + 1 latent model
  fm1 <- Y ~ f(~1, effect_type="normal") + f(~0+x, effect_type="normal") +
          f(x, model="ar1", group=group)

  # 2 random effects + 1 latent models (replicate = 2)
  fm2 <- Y ~ f(~1+x, effect_type="normal", replicate = repl) +
          # f(~0+x, effect_type="normal", replicate = repl) +
          f(x, model="ar1", group=group, replicate = repl)

load_all()
  out <- ngme(
    fm0,
    data = data.frame(Y=1:6, x=c(1,2,3,1,2,3), repl=c(1,1,1,2,2,2)),
    control_opt = control_opt(
      estimation = T,
      iterations = 5
    )
  )

out$replicates[[1]]$randeff[[1]]$B_reff


out$replicates[[1]]$randeff
out$replicates[[1]]$randeff[[1]]$n_params
out$n_params # beta (1) + ar1 (2) + reff(2) + reff(2) + mnoise(1)
out$replicates[[1]]$n_la_params

  out$n_params
  out
})

test_that("test on 1d random effects", {
  load_all()
  re(matrix(c(1:30), ncol=3), theta_K=c(1,2,3,1,2,3))

  # 1d rand eff
  group <- 50
  u <- rnorm(group, 0, 5)
  each_obs <- 10
  repl <- rep(1:group, each=each_obs)
  Z <- rnorm(group*each_obs, 0, 1)
  # simulate Y
  Y <- double(group*each_obs)
  for (i in 1:group) {
    Y[((i-1)*each_obs+1):(i*each_obs)] <-
      1 + Z[((i-1)*each_obs+1):(i*each_obs)]*u[i] +
        rnorm(each_obs)
  }
  library(lme4)
  lmer(Y ~ 0 + (0+Z|repl), data=data.frame(Y=Y, Z=Z, repl=repl))

  load_all()
  out <- ngme(
    Y ~ 0 + f(~0+Z, effect_type="normal", replicate=repl),
    data = data.frame(Y=Y, Z=Z, repl=repl),
    control_opt = control_opt(
      estimation = T,
      iterations = 500,
      n_parallel_chain = 1
    )
  )
  out
  traceplot(out, "effect1")

})

test_that("test estimation vs lmer", {

  # pure random effects model
  library(lme4)
  sleepstudy
  ggplot(sleepstudy, aes(x=Days, y=Reaction)) +
    geom_point() +
    facet_wrap(~Subject, scales="free_y") +
    geom_smooth(method="lm", se=FALSE)

  lmer(Reaction ~ Days + (Days|Subject), data=sleepstudy)

  load_all()
  ngme(
    Reaction ~ f(~Days, effect_type="normal"),
      # f(Days, model="ar1", group=Subject),
    data=sleepstudy,
    control_opt = control_opt(
      estimation = T,
      iterations = 10
    )
  )
})

test_that("test simulation", {
  Sigma <- matrix(c(20, 5, 5, 10), 2, 2)
  # simulate U ~ N(0, Sigma)
  U <- MASS::mvrnorm(100, mu=c(0,0), Sigma=Sigma)
  var(U[, 1])
  var(U[, 2])
  cov(U[, 1], U[, 2])

  # simulate Y = U + e
  replicate = 1:100
  x1 <- rnorm(100)
  x2 <- rexp(100)
  Y <- x1 * U[, 1] + x2 * U[, 2] + rnorm(100)


  # test
  Sigma
  category <- 100;
  each = 20

  U <- MASS::mvrnorm(category, mu=c(0,0), Sigma=Sigma)

  Z <- cbind(rep(1, category * each), rexp(category * each))
  Y <- double(category * each)
  for (i in 1:category) {
    Y[((i-1)*each+1) : (i*each)] <-
      Z[((i-1)*each+1) : (i*each), 1] * U[i, 1] +
      Z[((i-1)*each+1) : (i*each), 2] * U[i, 2] +
      rnorm(each)
  }
  replicate = rep(1:category, each=each)
  library(lme4)
  lmer(Y ~ 0 + (0+Z|replicate), data=data.frame(Y=Y, Z=Z, replicate=replicate))

  ngme(Y ~ f(
    ~0 + z1 + z2, effect_type="normal",
    replicate = replicate, W = U, fix_W = T),
       data=data.frame(Y=Y, z1=Z[, 1], z2=Z[, 2], replicate=replicate),
       control_opt = control_opt(
         estimation = T,
         iterations = 10,
         sampling_strategy = "all"
       )
  )
})
