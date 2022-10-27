# Simple scripts for test ngme function
{
  library(devtools);load_all()

  seed <- 10
  set.seed(seed)
}

{ ############  1. simulate AR with nig noise
  n_obs <- 1000
  ar_mu <- 4
  ar_sigma <- 1
  ar_eta <- 1

  ar1_process <- simulate(
    f(1:n_obs,
      model = "ar1",
      theta_K = 0.9,
      noise = noise_nig(
        mu = ar_mu,
        sigma = ar_sigma,
        nu = ar_eta
      )
    ),
    seed = 1
  )

  noise_theta_mu    <- -6
  noise_theta_sigma <- 2
  noise_theta_V     <- 1.7

  nig_noise <- simulate(
    noise_nig(
      theta_mu = noise_theta_mu,
      theta_sigma = noise_theta_sigma,
      theta_V = noise_theta_V,
    ),
    nsim = n_obs
  )

  # Y <- ar1_process + nig_noise

  # use normal noise
    Y <- ar1_process + rnorm(n_obs)
}
# f(1:n_obs,
#   model = "ar1",
#   theta_K = 0.9,
#   noise = noise_nig(
#     theta_mu = NULL,
#     B_sigma = matrix(c(1,2), ncol=2),
#     theta_sigma = c(1,2),
#     nu = ar_eta
#   ))$par_string

# { #simulate using old fashion
# n_obs <- 500
# alpha1 <- 0.9
# mu1 = 4; delta = -mu1
# sigma1 = 3
# nu1 = 1

# trueV1 <- ngme2::rig(n_obs, nu1, nu1)
# noise1 <- delta + mu1*trueV1 + sigma1 * sqrt(trueV1) * rnorm(n_obs)

# trueW1 <- Reduce(function(x,y){y + alpha1*x}, noise1, accumulate = T)
# Y = trueW1 + rnorm(n_obs, mean=0, sd = 0.5)
# }
load_all()
ngme_out <- ngme(
  Y ~ 0 +
  f(1:n_obs,
    model = "ar1",
    theta_K = 0.2,
    # fix_theta_K = TRUE,
    # W = as.numeric(ar1_process),
    # fix_W = TRUE,
    noise = noise_nig(
      theta_mu = 1.1,
      theta_sigma = 0.1,
      # theta_V = 2,
      # V = attr(ar1_process, "noise")$V,
      # fix_V = TRUE,
      fix_theta_mu      = F,
      fix_theta_sigma   = F,
      fix_theta_V       = F
    ),
    control = ngme_control_f(
      numer_grad       = F,
      use_precond      = T
    ),
    debug = T
  ),
  data = data.frame(Y = Y),
  control = ngme_control(
    estimation = T,
    exchange_VW = FALSE,
    n_parallel_chain = 4,
    stop_points = 50,
    burnin = 200,
    iterations = 1,
    gibbs_sample = 5,
    stepsize = 1,
    kill_var = FALSE,
    threshold = 1e-4,

    std_lim = 0.1,
    trend_lim = 0.1
  ),
  family = noise_normal(),
  # noise = attr(nig_noise, "noise"),
  seed = 10,
  # , last_fit = ngme_out
  debug = T
)
ngme_out

# noise
traceplot(ngme_out, parameter = "theta_mu", f_index = 0)
traceplot(ngme_out, parameter = "theta_sigma", f_index = 0)
traceplot(ngme_out, parameter = "theta_V", f_index = 0)

# ar1 model
traceplot(ngme_out, parameter = "theta_K",     f_index = 1)
traceplot(ngme_out, parameter = "theta_mu",    f_index = 1)
traceplot(ngme_out, parameter = "theta_sigma", f_index = 1)
traceplot(ngme_out, parameter = "theta_V",     f_index = 1)

# compare ar model
# plot(ngme_out$latents[[1]]$noise, col = "red")
plot(noise_nig(
      theta_mu = ar_mu,
      theta_sigma = ar_sigma,
      theta_V = ar_eta
    ),
    ngme_out$latents[[1]]$noise
)

load_all()

model_ar1(1:n_obs,
  theta_K = 0.9,
  noise = noise_nig(
    mu = ar_mu,
    theta_sigma = ar_sigma,
    nu = ar_eta
  ))

within(ll, rm(x))

fun <- function(a=0, b=0) {
  list(a=a,b=b)
}

ll <- list(a=1, b=2)
fun(unlist(ll))

?as.list

ll <- list(a=1,b=2,w=NULL)

ll
class(ll) <- "sth"
within(ll, {
  if (is.null(w)) k = 54
})


load_all()

rw1 <- model_rw1(1:3, noise=noise_normal())

rw1$C + rw1$G

?within
