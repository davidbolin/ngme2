test_that("type-G1 works", {
  ############################## Type-G1
  n <- 500
  load_all()
  bv_type_G1 <- f(
    zeta = 1, eta = 1,
    c(1:n, 1:n),
    model="bv",
    sub_models=list(
      f1 = list(model="rw1"),
      f2 = list(model="ar1", rho=0.5)
    ),
    group = c(rep("f1", n), rep("f2", n)),
    noise = list(
      f1 = noise_nig(single_V=T, nu=2, mu=-2),
      f2 = noise_nig(single_V=T, nu=2, mu=2),
      share_V = T
    ),
    eval = T
  )

  # simulate and estimate
  W <- simulate.ngme_model(bv_type_G1)
  Y <- W + rnorm(2*n)
  out <- ngme(
    Y ~ 0 + f(
      c(1:n, 1:n),
      model="bv",
      sub_models=list(
        f1 = list(model="rw1"),
        f2 = list(model="ar1")
      ),
      debug=T,
      group = c(rep("f1", n), rep("f2", n)),
      noise = list(
        f1 = noise_nig(single_V = T),
        f2 = noise_nig(single_V = T),
        share_V = T
      )
    ),
    data = data.frame(Y),
    family=noise_normal(),
    control_opt = control_opt(
      estimation = T,
      n_parallel_chain = 1,
      iterations = 10
    )
  )
  out
  traceplot(out, "field1")
})

test_that("type-G2 works", {
  ############################## Type-G2
  n <- 500
  bv_type_G2 <- f(
    c(1:n, 1:n),
    model="bv",
    sub_models=list(
      f1 = list(model="rw1"),
      f2 = list(model="ar1", rho=0.5)
    ),
    group = c(rep("f1", n), rep("f2", n)),
    noise = list(
      f1 = noise_nig(single_V=T),
      f2 = noise_nig(single_V=T),
      share_V = F
    ),
    eval = T
  )
  simulate(bv_type_G2)
})

test_that("type-G3 works", {
  ############################## Type-G3
  n <- 500
  bv_type_G2 <- f(
    c(1:n, 1:n),
    model="bv",
    sub_models=list(
      f1 = list(model="rw1"),
      f2 = list(model="ar1", rho=0.5)
    ),
    group = c(rep("f1", n), rep("f2", n)),
    noise = list(
      f1 = noise_nig(single_V=T),
      f2 = noise_nig(single_V=T),
      share_V = F
    ),
    eval = T
  )
  simulate(bv_type_G2)
})

test_that("type-G3 works", {
  load_all()
  n <- 500
  bv_type_G3 <- f(
    c(1:n, 1:n),
    model="bv",
    zeta = 1, eta=1,
    sub_models=list(
      f1 = list(model="rw1"),
      f2 = list(model="ar1", rho=0.5)
    ),
    group = c(rep("f1", n), rep("f2", n)),
    noise = list(
      f1 = noise_nig(mu=-2, nu=2),
      f2 = noise_nig(mu=2, nu=2),
      share_V = T
    ),
    eval = T
  )
  bv_type_G3

  # simulate and estimate
  W3 <- simulate.ngme_model(bv_type_G3)
  # share V
  all(head(attr(W3, "V"), n) == tail(attr(W3, "V"), n))
  Y <- W3 + rnorm(2*n)

  out3 <- ngme(
    Y ~ 0 + f(
      c(1:n, 1:n),
      model="bv",
      sub_models=list(
        f1 = list(model="rw1"),
        f2 = list(model="ar1")
      ),
      debug=T,
      group = c(rep("f1", n), rep("f2", n)),
      noise = list(
        f1 = noise_nig(single_V=F),
        f2 = noise_nig(single_V=F),
        share_V = T
        # ,fix_V = T, V = attr(W3, "V")
      )
    ),
    data = data.frame(Y),
    family=noise_normal(),
    control_opt = control_opt(
      estimation = T,
      n_parallel_chain = 4,
      iterations = 1000
    )
  )

  # 4 chain, 1000 iter 225s
  out3
  traceplot(out3, "field1")
})

test_that("type-G4 works", {
  load_all()
  n <- 500
  bv_type_G4 <- f(
    c(1:n, 1:n),
    model="bv",
    zeta = 1, eta=1,
    sub_models=list(
      f1 = list(model="rw1"),
      f2 = list(model="ar1", rho=0.5)
    ),
    group = c(rep("f1", n), rep("f2", n)),
    noise = list(
      f1 = noise_nig(mu=2, nu=2),
      f2 = noise_nig(mu=-2, nu=0.5)
    ),
    eval = T
  )

  # simulate and estimate
  W4 <- simulate.ngme_model(bv_type_G4)
  Y <- W4 + rnorm(2*n)
  out4 <- ngme(
    Y ~ 0 + f(
      c(1:n, 1:n),
      model="bv",
      sub_models=list(
        f1 = list(model="rw1"),
        f2 = list(model="ar1")
      ),
      debug=T,
      group = c(rep("f1", n), rep("f2", n)),
      noise = list(
        f1 = noise_nig(single_V=F),
        f2 = noise_nig(single_V=F),
        share_V = F
        # ,fix_V = T, V = attr(W4, "V")
      )
    ),
    data = data.frame(Y),
    family=noise_normal(),
    control_opt = control_opt(
      estimation = T,
      n_parallel_chain = 1,
      iterations = 300
    )
  )
  # 81s, good
  out4
  traceplot(out4, "field1")
})