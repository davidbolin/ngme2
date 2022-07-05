library(devtools); load_all()
# load_all(reset = FALSE, recompile = FALSE)

a2th <- function(k) {log((-1-k)/(-1+k))}
th2a <- function(th) {-1 + (2*exp(th)) / (1+exp(th))}

############  0. generating fix effects and control
set.seed(7)



control = ngme.control(burnin=100,
                       iterations = 500,
                       gibbs_sample = 5,
                       stepsize = 1,
                       kill_var = FALSE,
                       threshold = 1e-4,
                       opt_fix_effect = T)

# fix effects
# X <- (model.matrix(Y1 ~ x1 + x2))  # design matrix
# Y1 = as.numeric(Y1 + X %*% beta)

############  1. test AR with nig noise
# parameter for ar1
# x1 = runif(n_obs)
# x2 = rexp(n_obs)
# beta <- c(-3, -1, 2)

n_obs <- 1000
sigma_eps = 0.5
alpha <- 0.5
mu = 2; delta = -mu
nu = 1
sigma = 3

n_obs1 <- 2*n_obs
trueV1 <- ngme2::rig(n_obs1, nu, nu)
noise1 <- delta + mu*trueV1 + sigma * sqrt(trueV1) * rnorm(n_obs1)
trueW1 <- Reduce(function(x,y){y + alpha*x}, noise1, accumulate = T)
Y1 = trueW1 + rnorm(n_obs1, mean=0, sd=sigma_eps)

trueV2 <- ngme2::rig(n_obs, nu, nu)
noise2 <- delta + mu*trueV2 + sigma * sqrt(trueV1) * rnorm(n_obs)
trueW2 <- Reduce(function(x,y){y + alpha*x}, noise2, accumulate = T)
Y2 = trueW2 + rnorm(n_obs, mean=0, sd=sigma_eps)

Y <- c(Y1, Y2)

################################################################
index <- c(1:n_obs1, 1:n_obs)
replicates = c(rep(1, n_obs1), rep(2, n_obs))

ngme_out = ngme(Y ~ 0 +
                  f(index,
                    replicates = replicates,
                    model=ngme.ar1(
                      index,
                      replicates = replicates,
                      alpha=0.9,
                      use_num_dK = FALSE
                    ),
                    noise=ngme.noise(
                      type="nig",
                      theta.noise=1
                    ),
                    control=ngme.control.f(
                      numer_grad       = FALSE,
                      use_precond      = TRUE,

                      fix_operator     = FALSE,
                      fix_mu           = FALSE,
                      fix_sigma        = FALSE,
                      fix_noise        = FALSE
                    ),
                    theta.mu = mu+2,
                    theta.sigma = log(sigma)+2,
                    theta.noise = 1.01,
                    debug=TRUE
                  ),
                family="normal",
                data=data.frame(
                  index=index
                  # ,
                  # Y1=Y1,
                  # x1=x1,
                  # x2=x2
                ),
                control=control,
                start=ngme.start(
                  # W = c(trueW1, trueW2)
                ),
                debug=ngme.debug(
                  debug = TRUE,
                  fix_merr = FALSE
                ))

ngme_out$result

# result
res0 <- c(alpha, mu, sigma, nu); names(res0) <- c("alpha", "mu", "log(sigma)", "nu"); res0
res1 <- c(alpha, mu, log(sigma), nu); names(res1) <- c("alpha", "mu", "log(sigma)", "nu"); res1
res2 <- c(beta, sigma_eps); res2
str(ngme_out$output)

# plot alpha
  plot_out(ngme_out$trajectory, start=1, n=1, transform = th2a)
# plot mu
  plot_out(ngme_out$trajectory, start=2, n=1)
# plot sigma
  plot_out(ngme_out$trajectory, start=3, n=1, transform = exp)
# plot var
  plot_out(ngme_out$trajectory, start=4, n=1, transform = exp)
# plot fix effects
  # plot_out(ngme_out$trajectory, start=5, n=3)
# plot m err
  plot_out(ngme_out$trajectory, start=5, n=1, transform = exp)
# plot alpha
plot_out(ngme_out$trajectory, start=1, n=1, transform = th2a)
plot_out(ngme_out$trajectory, start=2, n=1)

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
ngme_out$result
