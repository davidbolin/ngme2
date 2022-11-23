# The `ngme2` Package

`ngme2` is an R package used for fitting non-gaussian mixed effects models. These models are fitted using maximum likelihood estimation and preconditioned stochastic gradient descent.

Basic statistical operations such as likelihood evaluations and kriging predictions are also implemented.

## Introduction

Several popular Gaussian random field models can be represented as solutions to stochastic partial differential equations (SPDEs) of the form
$$
L^\beta (\tau u) = \mathcal{W}.
$$
Here $\mathcal{W}$ is a Gaussian white noise, $L$ is a second-order differential operator, the fractional power $\beta > 0$ determines the smoothness of u. See [An explicit link between Gaussian fields and Gaussian Markov random fields: the stochastic partial differential equation approach](https://rss.onlinelibrary.wiley.com/doi/full/10.1111/j.1467-9868.2011.00777.x) for further details.

This package aims to address the non-Gaussian extension to the SPDE approach model by replacing the driven noise $\mathcal{W}$ to be a non-Gaussian noise $\mathcal{M}$.
More specificly, a type-G Lévy process.

The increment of a type-G Lévy process can be represented as
$$
\gamma + \mu V + \sigma \sqrt{V}Z,
$$
where $\gamma, \mu$ are parameters, $Z\sim N(0,1)$ and is independent of $V$, and $V$ is a positive infinitely divisible random variable.

One example is the normal inverse Gaussian (NIG) noise. (See `vignette("SPDE-approach", package = "ngme2")` for more details)

## Installation instructions #
The development version can be installed using the command
```r
remotes::install_github("davidbolin/ngme2", ref = "devel")
```

See also the [Installation and Configuration][ref] vignette.

[ref]: https://davidbolin.github.io/ngme2/articles/Installation_and_configuration.html "Installation and Configuration"