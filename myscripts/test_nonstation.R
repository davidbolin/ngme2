library(devtools)
load_all()
k2th <- function(k) {log((-1-k)/(-1+k))}
th2k <- function(th) {-1 + (2*exp(th)) / (1+exp(th))}

######## 1. test AR with non-stationary mu and simga
n_obs <- 1000
sigma_eps = 1

alpha <- 0.5
nu = 1

# mu <- 2;
# B.mu <- matrix(1, nrow=n_obs, ncol=1)
B.mu <- cbind(1, rexp(n_obs))
mu <- drop(B.mu %*% c(3, 5))
delta = -mu

# sigma = 2;
# B.sigma <- cbind(1, runif(n_obs))
B.sigma <- cbind(1, (1:n_obs)/n_obs)
sigma <- drop(exp(B.sigma %*% c(1.5, -1.2)))

# B.sigma <- B.sigma[, 1]
# sigma <- drop(exp(B.sigma * -0.4))
# range(sigma)
# sigma <- 1.1

trueV <- ngme2::rig(n_obs, nu, nu)
noise <- delta + mu*trueV + sigma * sqrt(trueV) * rnorm(n_obs)

trueW <- Reduce(function(x,y){y + alpha*x}, noise, accumulate = T)
Y = trueW + rnorm(n_obs, mean=0, sd=sigma_eps)
range(Y)

# fitting
control = control.ngme(burnin=100, iterations = 1000,
                       gibbs_sample = 5, stepsize = 1,
                       kill_var = FALSE, threshold = 1e-4,
                       opt_fix_effect = T)

# debug = debug.ngme(fixW = TRUE, trueW=trueW)

ngme_out = ngme(Y ~ 0 +
                  f(1:length(Y),
                    model = "ar1",
                    noise = ngme.noise(type="nig",
                                       theta.noise=1),
                    B.mu=B.mu,
                    theta.mu=c(0, 0),
                    B.sigma = B.sigma,
                    theta.sigma = log(c(1, 1)),
                    control = control.f(numer_grad     = FALSE,
                                      init_operator    = 0.5,
                                      opt_operator     = TRUE,
                                      opt_mu           = TRUE,
                                      opt_sigma        = TRUE,
                                      opt_var          = TRUE),
                    debug = TRUE),
                family = "normal",
                data = data.frame(Y=(as.numeric(Y))),
                control = control #,
                # debug=debug
)

###### show results
ngme_out$trajectory
ngme_out$estimates

# # plot alpha
plot_out(ngme_out$trajectory, start=1, n=1, transform = th2k)
# # plot mu
plot_out(ngme_out$trajectory, start=2, n=2)
# # plot sigma
#plot_out(ngme_out$trajectory, start=4, n=2, type="grad")
plot_out(ngme_out$trajectory, start=4, n=2)
# # plot var
plot_out(ngme_out$trajectory, start=6, n=1, transform = exp)
# # plot m err
plot_out(ngme_out$trajectory, start=7, n=1, transform = exp)



