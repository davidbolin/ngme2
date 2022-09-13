######## test NA ########
library(devtools); load_all()

n_obs <- 5
sim <- ngme.simulate(
  f(1:n_obs, model = "ar1", noise = ngme.noise.nig())
)
Y <- sim$realization + rnorm(n_obs)

index_NA <- sample(1:n_obs, 0.2 * n_obs)
Y[index_NA] <- NA
index_data <- which(!is.na(Y))

# finally
ngme_out <- ngme(
  formula = Y ~ 0 + f(
    model = "ar1",
    theta_K = 0.4,
    noise = ngme.noise.nig(),
    debug = TRUE
  ),
  data = data.frame(Y = Y),
  noise = ngme.noise.normal(),
  control = ngme.control(
    iteration = 2
  ),
  debug = ngme.debug(not_run = F)
)

str(ngme_out)

index <- c(1,2,3)
index_NA <- 2

Filter(index, index %in% index_NA)

# model.frame(Y~0, na.action=NULL)

load_all()
ngme.simulate(f(model = "ar1"))
