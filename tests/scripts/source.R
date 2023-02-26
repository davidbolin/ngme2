  load_all()
  n <<- 10
  z1 = arima.sim(n, model = list(ar = 0.5), sd = 0.5)
  z2 = arima.sim(n, model = list(ar = 0.5), sd = 0.5)

  out2 <- ngme(
    z ~ f(c(1:n ,1:n), model="rw1", replicate=rep(1:2,each=n), noise=noise_normal()),
    data = list(z = c(z1, z2)),
    control = ngme_control(
      estimation = T,
      iterations = 20,
      n_parallel_chain = 1
    ),
    debug = TRUE
  )


check()
