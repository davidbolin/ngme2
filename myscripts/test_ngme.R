library(devtools)
# load_all(reset = FALSE, recompile = FALSE)
load_all()

a2th <- function(k) {log((-1-k)/(-1+k))}
th2a <- function(th) {-1 + (2*exp(th)) / (1+exp(th))}

############  0. generating fix effects and control

n_obs <- 1000
sigma_eps = 0.5
x1 = runif(n_obs)
x2 = rexp(n_obs)
beta <- c(-3, -1, 2)

control = ngme.control(burnin=30,
                       iterations = 1,
                       gibbs_sample = 5,
                       stepsize = 1,
                       kill_var = FALSE,
                       threshold = 1e-4,
                       opt_fix_effect = T)

############  1. test AR with nig noise
# parameter for ar1

alpha1 <- 0.5
mu1 = 2; delta = -mu1
nu1 = 1
sigma1 = 3

trueV1 <- ngme2::rig(n_obs, nu1, nu1)
noise1 <- delta + mu1*trueV1 + sigma1 * sqrt(trueV1) * rnorm(n_obs)

trueW1 <- Reduce(function(x,y){y + alpha1*x}, noise1, accumulate = T)
Y1 = trueW1 + rnorm(n_obs, mean=0, sd=sigma_eps)

# fix effects
X <- (model.matrix(Y1 ~ x1 + x2))  # design matrix
Y1 = as.numeric(Y1 + X %*% beta)

c(alpha1, mu1, sigma1, nu1); c(beta, sigma_eps)

ngme_out = ngme(Y1 ~ x1 + x2 +
                  f(1:length(Y1),
                    model="ar1",
                    noise=ngme.noise(type="nig", theta.noise=1),
                    control=ngme.control.f(
                      numer_grad       = FALSE,
                      theta.K          = 0.9,
                      opt_operator     = T,
                      opt_mu           = T,
                      opt_sigma        = T,
                      opt_var          = F,
                      use_precond      = TRUE
                    ),
                    theta.mu = mu1,
                    theta.sigma = log(sigma1),
                    debug=TRUE),
                family="normal",
                data=data.frame(Y1=(as.numeric(Y1)), x1=x1, x2=x2),
                control=control,
                debug=ngme.debug(
                  debug = TRUE,
                  fixSigEps = FALSE,
                  sigEps = 0.5
                ))
# plot alpha
  plot_out(ngme_out$trajectory, start=1, n=1, transform = th2a)
# plot mu
  plot_out(ngme_out$trajectory, start=2, n=1)
# plot sigma
  plot_out(ngme_out$trajectory, start=3, n=1, transform = exp)
# plot var
  plot_out(ngme_out$trajectory, start=4, n=1, transform = exp)
# plot fix effects
  plot_out(ngme_out$trajectory, start=5, n=3)
# plot m err
  plot_out(ngme_out$trajectory, start=8, n=1, transform = exp)
# plot alpha

plot_out(ngme_out$trajectory, start=1, n=1, transform = th2a)
  ngme_out$estimates

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
  ngme_out$estimates
  ngme_out$trajectory

# plot(ngme_out, param = "fe", type = "traj")
# plot(ngme_out, param = "me", type = "traj")
# plot(ngme_out, param = "la", type = "traj", which=1)
# plot(ngme_out, param = "la", type = "traj", which=2)
# summary(ngme_out)


