# The `ngme2` Package  <img src="man/figures/logo.png" align="right" width="80" height="80" />

## Description
`ngme2` (https://davidbolin.github.io/ngme2/) is a unified, efficient, and flexible framework for fitting latent **non-Gaussian** models in R. It extends the SPDE-based Gaussian modeling toolkit to handle skewness, heavy tails, and non-smooth behavior while keeping familiar workflows for estimation, prediction, and model assessment.

## Features
- Temporal processes: AR(1), Ornstein–Uhlenbeck, random walks, and ARMA models.
- Spatial processes: Matérn fields (including graph-based meshes) and separable/non-separable space–time variants.
- Non-Gaussian driven noise: NIG, generalized asymmetric Laplace, skew-\(t\); non-Gaussian measurement noise is also supported.
- Joint and custom structures: bivariate type-G fields, longitudinal random-effects models, and user-defined operators via `generic` / `generic_ns`, flexible combinations of different models and driven noise.
- Practical tools: kriging-style prediction, cross-validation helpers, and diagnostics for convergence/fit quality.

## Quick start
```r
library(ngme2)
time_index <- seq(1, 1000, by = 1)
n <- length(time_index)

# Define the AR(1) model with NIG driven noise
ar1_nig <- f(time_index,
  model = ar1(rho = 0.7),
  noise = noise_nig(mu = 3, sigma = 2, nu = 0.5)
)

# Simulate the AR(1) process with NIG driven noise
nig_field <- simulate(ar1_nig, seed = 123, nsim = 1)[[1]]
Y <- nig_field + rnorm(n, mean = 0, sd = 1)
plot(time_index, nig_field, type = "l")

# Fit the model
fit <- ngme(
  formula = Y ~ 0 + f(time_index, model = ar1(), noise = noise_nig()),
  data    = data.frame(Y = Y, time_index = time_index),
  family  = "normal" # likelihood family
)

summary(fit)

```

Use `ngme_optimizers()` to see available optimizers and configure stochastic gradient settings via `control_opt`.

## Installation
The stable version can be installed with:
```r
install.packages("ngme2", repos = "https://davidbolin.github.io/ngme2/")
```
See the [Installation and Configuration](https://davidbolin.github.io/ngme2/articles/Installation.html) vignette if compilation tools are needed.

## Learn more
- [Get Started](https://davidbolin.github.io/ngme2/articles/ngme2.html): for a guided tour of the modeling framework.
- [Online vignettes](https://davidbolin.github.io/ngme2/articles/)
- [Online documentations](https://davidbolin.github.io/ngme2/reference/index.html)
