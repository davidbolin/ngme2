# ngme

The development version can be installed using the command

```r
remotes::install_github("davidbolin/ngme2", ref = "devel")
```

## done

* matern model
* non-stationary matern
* non-stationery nig noise (e.g. mu)
* replicates of the process (in the f function)

## to-do

* (To be done jointly): homepage, vignette + documentation

* non-gaussian measurement noise + analytical non-stationary matern
    family = "nig", parameter to be estimate : noise_mu, noise_sigma
    1. hardcode every noise in block model
    2. make sampleW -> sampleW(noise_type)
* analyze the Argo dataset
* handle temporal dependence
* random effects
* parallel chains
* hessian from all other chains
* convergence checks
* multivariate matern model
* multivariate measurement noise
* prior for the parameters
* confidence intervals for the parameters
* sharing parameters between different models
* rational approximation
* write a book
