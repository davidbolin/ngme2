# test nig measurement nosie

n_obs <- 500
sigma_eps <- 0.5
alpha <- 0.5
mu <- 2; delta <- -mu
nu <- 1
sigma <- 3

n_obs <- n_obs
trueV <- ngme2::rig(n_obs, nu, nu)
noise <- delta + mu*trueV + sigma * sqrt(trueV) * rnorm(n_obs)
trueW <- Reduce(function(x, y) {y + alpha*x}, noise, accumulate = T)

noise_nu <- 2.5
noise_mu <- 1.5
noise_sigma <- 1.5
noise_V <- ngme2::rig(n_obs, noise_nu, noise_nu)
nig_merr <- -noise_mu + noise_mu*noise_V + 
    noise_sigma*sqrt(noise_V)*rnorm(n_obs)

Y <- trueW + nig_merr



