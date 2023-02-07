# ngme2 demo

# doing replicate using INLA
library(INLA)
n = 300
z1 = arima.sim(n, model = list(ar = 0.5), sd = 0.5) # independent replication
z2 = arima.sim(n, model = list(ar = 0.5), sd = 0.5) # from AR(1) process
y = c(z1, z2)

idx <- c(1:n, 1:n)
repl <- rep(1:2, each=n)

formula <- y ~ -1 + f(idx, model="ar1", replicate = repl)
res_inla <- inla(
  formula,
  family="gaussian",
  data = list(y=c(z1, z2), idx=idx, repl=repl)
)
summary(res_inla)
sqrt(1 / 3.22)

load_all()
formula <- y ~ -1 + f(idx, model="ar1", replicate = repl)
res_ngme <- ngme(
  formula,
  family="gaussian",
  data = list(y=c(z1, z2), idx=idx, repl=repl)
)
res_ngme
