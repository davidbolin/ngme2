a2th <- function(k) {log((-1-k)/(-1+k))}
th2a <- function(th) {-1 + (2*exp(th)) / (1+exp(th))}

library(devtools); load_all()
# load_all(reset = FALSE, recompile = FALSE)

############  0. generating fix effects and control
seed <- 8
set.seed(seed)

control = ngme.control(burnin = 100,
                       iterations = 5,
                       gibbs_sample = 5,
                       stepsize = 1,
                       kill_var = FALSE,
                       threshold = 1e-4,
                       opt_fix_effect = T)

# fix effects
# X <- (model.matrix(Y1 ~ x1 + x2))  # design matrix
# Y1 = as.numeric(Y1 + X %*% beta)

############  1. test AR with nig noise

n_obs <- 1000
alpha <- 0.5
mu = 2; delta = -mu
nu = 1
sigma = 3

# non-stationary mesurement noise
## sigma
B_noise_sigma <- cbind(1, (1:n_obs) / n_obs);
# theta_noise_sigma <- c(0, 0)
theta_noise_sigma <- c(1, -2)
(noise_sigma <- drop(exp(B_noise_sigma %*% theta_noise_sigma)))

# B_noise_mu <- cbind(1, (1:n_obs) / n_obs);
B_noise_mu <- (1:n_obs) / n_obs
theta_noise_mu <- c(4)
# theta_noise_mu <- c(0, 0)
(noise_mu <- (B_noise_mu * theta_noise_mu))
B_noise_mu <- matrix(B_noise_mu, ncol=1)

noise_V <- rep(1, n_obs)
noise_V <- ngme2::rig(n_obs, nu, nu, seed=seed)
m_noise <- noise_mu * (-1 + noise_V) + noise_sigma * sqrt(noise_V) * rnorm(n_obs)

n_obs1 <- n_obs
trueV1 <- ngme2::rig(n_obs1, nu, nu, seed=seed)
noise1 <- delta + mu*trueV1 + sigma * sqrt(trueV1) * rnorm(n_obs1)
trueW1 <- Reduce(function(x,y) {y + alpha*x}, noise1, accumulate = T)

Y1 = trueW1 + m_noise

# B_noise_mu
# ngme.noise(
#   type = "nig",
#   theta_sigma = c(0, 0),
#   B_sigma = B_noise_sigma
#   , theta_mu = c(1),
#   B_mu = B_noise_mu
# )

################################################################
# index <- c(1:n_obs1, 1:n_obs)
replicates = c(rep(1, n_obs1), rep(2, n_obs))
print(Y1)

ngme_out = ngme(
  Y1 ~ 0 +
  f(1:n_obs1,
    replicates = NULL,
    model = "ar1",
    theta_K = 0.5,
    noise = ngme.noise(
      type        = "nig",
      theta_V     = 1,
      theta_mu    = mu + 1,
      theta_sigma = log(sigma) + 1
    ),
    control = ngme.control.f(
      numer_grad       = FALSE,
      use_precond      = TRUE,
      fix_operator     = TRUE,
      fix_mu           = TRUE,
      fix_sigma        = TRUE,
      fix_noise        = TRUE,
      fix_V            = FALSE,
      fix_W            = FALSE,

      init_V           = NULL,
      init_W           = NULL
    ),
    debug = FALSE
  ),
  data=data.frame(),
  control=control,
  noise = ngme.noise(
    type = "nig",
    # type = "normal",
    theta_sigma = c(0, 0),
    B_sigma = B_noise_sigma
    , theta_mu = c(0),
    B_mu = B_noise_mu
  ),
  debug = ngme.debug(
    debug = TRUE,
    fix_merr = FALSE
  ),
  seed = 1
)
ngme_out$result

# result
# res0 <- c(alpha, mu, sigma, nu); names(res0) <- c("alpha", "mu", "log(sigma)", "nu"); res0
# res1 <- c(alpha, mu, log(sigma), nu); names(res1) <- c("alpha", "mu", "log(sigma)", "nu"); res1
# res2 <- c(beta, sigma_eps); res2
# str(ngme_out$output)


# plot alpha
#   plot_out(ngme_out$trajectory, start=1, n=1, transform = th2a)
# # plot mu
#   plot_out(ngme_out$trajectory, start=2, n=1)
# # plot sigma
#   plot_out(ngme_out$trajectory, start=3, n=1, transform = exp)
# # plot var
#   plot_out(ngme_out$trajectory, start=4, n=1, transform = exp)
# # plot fix effects
#   # plot_out(ngme_out$trajectory, start=5, n=3)
# # plot m err
#   plot_out(ngme_out$trajectory, start=5, n=1, transform = exp)
# # plot alpha
# plot_out(ngme_out$trajectory, start=1, n=1, transform = th2a)
# plot_out(ngme_out$trajectory, start=2, n=1)

############  1.2 construct AR with normal noise
# parameter for ar1
# alpha2 <- 0.7
# sigma2 = 2
#
# noise2 = sigma2 * rnorm(n_obs)
#
# trueW2 <- Reduce(function(x,y){y + alpha2*x}, noise2, accumulate = T)
#
# Y2 = trueW2 + rnorm(n_obs, mean=0, sd=sigma_eps)
#
# # fix effects
# X <- (model.matrix(Y2 ~ x1 + x2))  # design matrix
# Y2 = as.numeric(Y2 + X %*% beta)

########### 2. generate fixed effect and measurement noise

# Y <- trueW1 + trueW2 + rnorm(n_obs, mean=0, sd=sigma_eps)
#
# X <- (model.matrix(Y ~ x1 + x2))  # design matrix
# beta <- c(-3, -1, 2)
#
# Y = as.numeric(Y + X %*% beta)

########### 3. run ngme

# specify the control for ngme
# args(control.ngme)
# args(control.f)
# args(control.ngme)

##### ngme for 1 ar
# nig
# th2k(0.5)


# normal
# ngme_out = ngme(Y2 ~ x1 + x2 +
#                   f(Y2, model="ar1", var="normal",
#                     control=control.f(numer_grad = FALSE)),
#                 data=data.frame(Y2=(as.numeric(Y2)), x1=x1, x2=x2),
#                 control=control)


##### ngme for 2 ar
# ngme_out = ngme(Y1 ~ 0 +
#                   # f(1:length(Y2), model="ar1", var="normal"),
#                   f(1:length(Y1), model="ar1", var="nig", debug = TRUE),
#                 data     = data.frame(Y1=Y1),
#                 control  = control,
#                 debug    = debug)

########### 4. results

# xx <- seq(-10, 10, length=1000);
# plot(xx, th2k(xx))

# nig ar1: a=0.5, mu=2, sigma=2, nu=1
# normal ar1: a=0.7 sigma=3

# ngme_out
  # ngme_out$output
  # ngme_out$trajectory

# plot(ngme_out, param = "fe", type = "traj")
# plot(ngme_out, param = "me", type = "traj")
# plot(ngme_out, param = "la", type = "traj", which=1)
# plot(ngme_out, param = "la", type = "traj", which=2)
# summary(ngme_out)
# ngme_out$result
