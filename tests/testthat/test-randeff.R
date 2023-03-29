# Y <- c(1:5, 1:5)
# group <- c(rep(1, 5), rep(2, 5))
# Y ~ Fixeff(~1+x, group=1) + matern_nd(..., group=group)

# test on random effects
test_that("the R interface", {
  f(~1+x, data=data.frame(x=c(1,2,3)), effect_type="normal", eval=T)

  # 1 random effect + 1 latent model
  fm0 <- Y ~ f(~1+x, effect_type="normal") +
          f(x, model="ar1")

  # 2 random effects + 1 latent model
  fm1 <- Y ~ f(~1, effect_type="normal") + f(~0+x, effect_type="normal") +
          f(x, model="ar1", group=group)

  # 2 random effects + 1 latent models (replicate = 2)
  fm2 <- Y ~ f(~1, effect_type="normal", replicate = repl) +
          f(~0+x, effect_type="normal", replicate = repl) +
          f(x, model="ar1", group=group, replicate = repl)

load_all()
  out <- ngme(
    fm2,
    data = data.frame(Y=1:6, x=c(1,2,3,1,2,3), repl=c(1,1,1,2,2,2)),
    control_opt = control_opt(
      estimation = F,
      iterations = 10
    )
  )

out$replicates[[1]]$B_reffs

out$replicates[[1]]$randeff
out$replicates[[1]]$randeff[[1]]$n_params
out$n_params # beta (1) + ar1 (2) + reff(2) + reff(2) + mnoise(1)
out$replicates[[1]]$n_la_params

  out$n_params
  out
})

cbind(randeffs[[1]]$B_reff, randeffs[[2]]$B_reff)
