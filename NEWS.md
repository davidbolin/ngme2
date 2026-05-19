# ngme2 0.9.8 (2026-05-19)

* Add `pm25_quarterly_2022`, a processed example dataset of quarterly
  PM2.5 monitor averages from EPA AQS/AirData for the contiguous United States
  in 2022. The dataset includes monitor location, quarter, quarterly average
  PM2.5, contributing day counts, and a replicate index for quarterly models.

# ngme2 0.9.7 (2026-04-29)

* Add `noise_moments()` for stationary `ngme_noise` objects, with direct
  helpers `noise_nig_moments()` and `noise_gal_moments()` returning variance,
  standard deviation, skewness, kurtosis, excess kurtosis, and central moments.
* Improve prediction for composite latent noises (`normal_nig`, `normal_gal`,
  and `nig_gal`) by expanding prediction projection matrices when required and
  using the output-sample prediction path for non-mean estimators.
* Make inverse-exponential prior calibration less intrusive to the global RNG
  state by using local seeding for generated Gaussian samples.
* Update pkgdown configuration for the current reference site.

# ngme2 0.9.6 (2026-04-17)

* Improve estimation and sampling console output: `ngme()` now respects
  `control_opt(verbose = ...)` more consistently and uses messages instead of
  unconditional `cat()` output.
* Improve C++ estimation robustness for serial/OpenMP builds, including
  single-chain variance handling and guarded OpenMP-only calls.
* Make print and plotting helpers more S3-friendly: print methods now return
  their input invisibly, and `traceplot()` no longer prints last estimates as a
  side effect.
* Make `simulate.ngme()` compatible with the S3 `simulate()` interface by
  accepting `posterior` and `m_noise` through `...`.
* Export `get_data_from_formula()` and improve diagnostics that depend on
  captured model/design output.
* Refresh package metadata, examples, and roxygen documentation, including
  SuiteSparse attribution and a fuller package description.

# ngme2 0.9.5 (2026-03-17)

## New feature: VAR(1) bivariate latent model

* Add `var1()` — a Vector Autoregressive order-1 latent field for bivariate
  time series modelling.
* Stationarity is **guaranteed by construction** via the Cayley reparameterization:
  four unconstrained parameters $(p_1, p_2, p_3, p_4)$ are mapped to a $2\times2$
  VAR coefficient matrix $A$ with spectral radius $\rho(A) < 1$ at every iteration.
* The C++ backend implements the Cayley transform natively (class `RCallback` in
  `src/latents/rcallback.cpp`).
* Supports all `ngme2` noise distributions (NIG, GAL, normal) as a single
  shared innovation noise.
* `print()` method displays the recovered $A$ matrix, its spectral radius, and
  the raw $(p_1, p_2, p_3, p_4)$ values.
* Add vignette *"VAR(1) Bivariate Model in ngme2"* (`vignettes/var1-model.Rmd`)
  covering: model specification, Cayley reparameterization, simulation study
  with parameter recovery, convergence trace plots, and NIG vs Gaussian
  model comparison.

# ngme2 0.9.4 (2026-03-12)
* Add posterior distribution plotting for SGLD samples via
  `posterior_plot()`.
* Add `plot()` support for `ngme_sgld_ci` objects, reusing stored SGLD
  samples to visualize marginal posterior distributions.

# ngme2 0.9.3 (2026-02-26)
* Refine fixed-effect standardization: SVD now applies only to non-intercept
  columns; intercept columns remain on their original parameterization.
* Preserve no-intercept model semantics (`~ 0 + ...`): skip fixed-effect
  centering when no intercept is present.
* Improve `fe()` centering with structural zeros: grouped `fe()` columns are
  centered using in-group rows only, so out-of-group structural zeros remain
  zero.
* Fix fixed-effect scale restoration for multi-replicate models: when needed,
  reconstruction now uses replicate row indices (`data_idx`) instead of always
  taking the first `n` rows.
* Improve warm-start robustness (`start = previous_fit`) across
  standardization settings by remapping fixed effects through the current model
  parameterization.
* Update default fixed-effect priors: intercept-like columns
  (`"(Intercept)"*`) now default to `prior_none()`, while non-intercept
  columns keep the default `N(0,10)` prior.
* For `standardize_fixed = TRUE`, add prior compatibility handling:
  isotropic normal priors on standardized columns are transformed to the SVD
  basis; incompatible custom `prior_beta` specifications now automatically
  disable fixed-effect standardization with a warning.

# ngme2 0.9.2 (2026-02-24)
* Add `prior_inv_exponential(lambda, lower)` for `nu`, implementing
  `kappa = 1 / nu ~ Exp(lambda)` as a first-class prior option.
* Add shorthand alias `prior_inv_exp(...)` for the same prior.
* Add calibration helpers `calibrate_inv_exp_lambda_driven_nig()` and
  `calibrate_inv_exp_lambda()` for choosing `lambda` from a driven-noise
  tail-inflation target.
* Improve calibration robustness for non-monotone `R_c(nu)` curves: the
  helper now scans for crossings and reports observed `R_c` range when the
  requested target is unattainable.
* Update default `nu` prior in `f()` for NIG-driven noise: when `nu` prior is
  not explicitly set and `nu` is stationary, use
  `prior_inv_exp(lambda = log(2)/median(h), lower = nu_lower_bound)`.
  For non-stationary `nu`, keep the legacy `N(0,10)` default prior.

# ngme2 0.9.1 (2026-02-19)
* Harden error handling in `ngme()` estimation/sampling path: C++ exceptions
  are now propagated as R errors (including OpenMP parallel regions) instead of
  potentially terminating the R session.
* Fix `nu` initialization in noise helper constructors to respect
  `nu_lower_bound`, using `theta_nu = log(nu - nu_lower_bound)` and validating
  `nu > nu_lower_bound`.
* Align `normal_nig` conversion, printing, and plotting with effective
  parameterization `nu = nu_lower_bound + exp(theta_nu)`.

# ngme2 0.9.0 (2026-02-19)
* Refactor prior API (breaking change):
  use `prior_normal()`, `prior_pc_sd()`, `prior_half_cauchy()`, `prior_none()`, and `priors(...)`.
* Update `f()` and `ngme_noise()` to accept unified `prior = ...` inputs
  (remove `prior_theta_K` and `prior_mu/prior_sigma/prior_nu` arguments).
* Remove legacy `ngme_prior()` interface and its documentation entry.
* Add prior target support (`coef`/`field`) for noise parameter priors and
  per-parameter operator prior compilation.
* Add fixed-effect prior support via `ngme(..., prior_beta = ...)`, using the
  same `prior_*()`/`priors(...)` API.
* Add user-facing vignette: `Prior Templates for Stationary and Non-Stationary Models`.

# ngme2 0.8.5 (2026-02-11)
* Add grad.norm plateau-based step decay via `control_opt(stepsize_decay = "grad_norm_plateau")` (epoch-level, synchronized across chains)
* Add `stepsize_decay()` helper for configuring decay options
* Update verbose output to include stepsize decay scale and effective stepsize
* Improve `cross_validation(data = ...)` model rebuild for refit-on-new-data workflows:
  it now resolves external formula symbols (for example `mesh`, `B`, `n_basis`) from the fitted object when needed, and falls back to rebuild-without-start plus hyperparameter transplant if `start` state dimensions differ.
* Add chain-aware prediction/CV aggregation via `chain_combine = "predictive_average"` in `predict()` and `cross_validation()`, which averages predictions across optimization chains instead of averaging parameters first.

# ngme2 0.8.4 (2026-02-01)
* Fix iid model using argument mesh instead of map

# ngme2 0.8.3 (2025-12-10)
* Update nu parameterization to be relative to nu_lower_bound (nu = nu_lower_bound + exp(B_nu * theta_nu))
* Allow NULL specification for noise parameter in non-stationary case
* Update fixed effects initialization
* Update start=fit logic

# ngme2 0.8.2 (2025-12-02)
* Add support for different type of ldlt solvers
* Use solver_type and solver_backend in `control_opt` to specify the solver
* Give print when hit nu lower bound

# ngme2 0.8.1 (2025-12-01)
* Add ldlt solver (solver_type = "ldlt")
* Minor fix on matern model on handling integer alpha case

# ngme2 0.8.0 (2025-11-30)
* Simplify solver structure
* Update to use model() instead of string for model definition
* Big update on the R interface, use model() instead of string for model definition
* Refactor on the codebase
* Add more vignettes

# ngme2 0.7.1
* Documentation updates
* Add print log likelihood for Gaussian models
* Improvements to prediction with condition on specific data points
* Enhanced preconditioner with Gibbs samples
* Updates to AR1 vignette
* Various optimization workflow improvements

# ngme2 0.6.0
* First version of the package

# ngme2 0.3.0
* Add replicate feature
* Add OU process
* Add tensor-product model
