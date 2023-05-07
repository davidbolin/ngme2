# ngme2 demo

# doing replicate using INLA
library(INLA)
n = 500
z1 = arima.sim(n, model = list(ar = 0.5), sd = 0.5) # independent replication
z2 = arima.sim(n, model = list(ar = 0.5), sd = 0.5) # from AR(1) process
y = c(z1, z2)

idx <- c(1:n, 1:n)
repl <- rep(1:2, each=n)

# use INLA
formula <- y ~ -1 + f(idx, model="ar1", replicate = repl)
res_inla <- inla(
  formula,
  family="gaussian",
  data = list(y=c(z1, z2), idx=idx, repl=repl)
)
summary(res_inla)

# use ngme2
formula <- y ~ -1 + f(idx, model="ar1", replicate = repl)
res_ngme <- ngme(
  formula,
  family="gaussian",
  data = list(y=c(z1, z2), idx=idx, repl=repl)
)
res_ngme

# Doing simulation using ngme2
myar <- f(1:n, model="ar1", alpha = 0.75,
  noise=noise_nig(
  mu = -3,
  sigma = 2,
  nu = 2
))

W1 <- simulate(myar)
W2 <- simulate(myar)
Y = c(W1, W2) + rnorm(2*n, sd=0.9)

res_ngme2 <- ngme(
  y ~ -1 + f(idx, model="ar1", replicate = repl, noise=noise_nig()),
  family="gaussian",
  data = list(y=Y, idx=idx, repl=repl),
  control = control_opt(
    iterations = 1000
  )
)
res_ngme2

plot(myar$noise,
  res_ngme2$models[["field1"]]$noise)
