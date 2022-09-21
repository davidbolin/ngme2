library(devtools); load_all()

y <- rnorm(10)
ngme.rw1(y)$C + ngme.rw1(y)$G

rw1 <- ngme.rw1(y)

str(f(y, model = rw1,
    noise = ngme.noise.nig(),
    debug = TRUE
))

out <- ngme(
  y ~ 0 + f(y, model = rw1,
    noise = ngme.noise.nig(),
    debug = TRUE
  ),
  data = list(y = y),
  noise = ngme.noise.normal(),
  debug = TRUE,
  control = ngme.control(
    iterations = 2,
    estimation = TRUE
  )
)
out
