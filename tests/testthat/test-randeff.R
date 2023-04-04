# Y <- c(1:5, 1:5)
# group <- c(rep(1, 5), rep(2, 5))
# Y ~ Fixeff(~1+x, group=1) + matern_nd(..., group=group)

# test on random effects
test_that("the R interface", {
  load_all()
  f(1:10, model="ar1")
  m = f(~1+x, data=data.frame(x=c(1,2,3)), effect_type="normal", eval=T)
  str(m)
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
      estimation = F,
      iterations = 10
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

test_that("test estimation", {
  # U ~ N(0, 2)
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