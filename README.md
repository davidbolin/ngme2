# ngme

The development version can be installed using the command

```r
remotes::install_github("davidbolin/ngme2", ref = "devel")
```

## to-do

* (To be done jointly): homepage, vignette + documentation

* non-gaussian measurement noise + analytical non-stationary matern
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

## doing now

**Non-gaussian measurement noise**
family = "nig", parameter to be estimate : noise_mu, noise_sigma, noise_nu
    1. hardcode every noise in block model or split into classes


## bug found

1. In printing grad. for the non-stationary mu.
2. In ngme.start, change the name of the output to theta_merr.

## Meeting notes

## done

* Matern model
* non-stationary Matern
* non-stationery nig noise (e.g. mu)
* replicates of the process (in the f function)