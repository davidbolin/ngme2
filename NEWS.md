# ngme2 (development version)

* Rao-Blackwellisation is now on by default (`control_opt(rao_blackwellization = )`).
* Fixed the Gibbs sweep order that made Rao-Blackwellisation unusable with a non-Gaussian
  `family`. The measurement `V` was redrawn after `W`, which left `QQ` and `cond_W`
  conditioned on the previous `V` while the gradient used the new one; the fit diverged until
  `QQ` lost positive-definiteness. Both `V` blocks are now drawn before `W`. Gaussian
  families were unaffected, since there `V` is fixed.
* The Rao-Blackwell traces are taken from an exact selected (Takahashi) inverse when the
  Cholesky factor of `QQ` is cheap enough, and from Hutchinson probes otherwise, controlled
  by `control_opt(selinv_max_fill = )`. Triangular operators (`ar1`, `ou`, `arma`) and 1-d
  meshes qualify; a 2-d mesh does not. For those models the traces are now exact and the fit
  no longer depends on `n_trace_iter` at all.
* The Hutchinson probe count adapts during the run (`control_opt(trace_adapt = )`) so that
  the trace estimator contributes a target share of the gradient variance instead of being
  fixed a priori. 
* The convergence trend test is now scale-free: it compares the fitted slope to its own
  standard error rather than to an absolute bound, so `trend_lim` is a t-value (default 2)
  and no longer depends on the parameter's units or on `n_batch`. The old absolute test never
  passed on any model tried; `trend_use_tstat = FALSE` restores it. The relative-standard-
  deviation part is off by default (`use_std_check`), as it conflates chain disagreement with
  parameter imprecision.
* A run that exhausts its iteration budget without converging now warns and names the
  parameters still failing, instead of returning silently. Set
  `control_opt(warn_no_convergence = FALSE)` for short runs where convergence is not expected.
* New post-convergence polish phase (`control_opt(polish_iterations = )`, default 50, capped
  at the number of iterations the fit ran). Once the stopping rule fires the step size is
  scaled by `polish_stepsize_factor` and the reported estimate becomes the average of the
  iterates over that window (Polyak-Ruppert), which cut the estimator's standard deviation by
  about a fifth in testing. Set to 0 for the previous last-batch average.
* Removed the Pflug convergence diagnostic (`pflug_conv_check`, `pflug_alpha`). It could
  never fire: the rule looks for sign changes in the running sum of `<g_t, g_{t-1}>`, but the
  gradient is evaluated on a Gibbs sample of the latent field, so consecutive gradients stay
  positively correlated however close the optimum is. Its intent -- detect that systematic
  progress has stopped -- is served by the slope t-statistic, which asks the same question of
  the parameter trajectory.
* `src/Makevars` now tracks header dependencies, so editing a header rebuilds every object
  that includes it. Previously only each source file was compared against its own object,
  which could leave objects compiled against two different layouts of the same class.
* Fixed `single_V` gradients for the NIG and t noise types, which did not match the score of
  the corresponding density and had the wrong sign over part of the range.
* The `"ws"` replicate sampling strategy now applies the importance weight it was missing, so
  it no longer differs in scale from `"all"` by roughly the number of replicates.
* Update estimation code to only reassemble and refactorize the latent precision
  matrix `QQ` when it can actually have changed. The matrix is also split into 
  `K' diag(1/SV) K` and its measurement block `(AZ)' D (AZ)`, and these are 
  caches separately. Results are unchanged by the caching and setting the 
  environment variable `NGME2_DISABLE_FACTOR_CACHE` restores the uncached behaviour.
* Fixed a bug for non-symmetric operators, which made the preconditioner Hessian
  incorrect and the first-order trace far noisier than necessary. These are now 
  handled by factorizing `K` directly with a sparse LU. The new 
  `control_opt(nonsym_solver = )` selects it: `"lu"` (default) or `"normal_equations"` 
  to reproduce results from 0.9.8 and earlier. 
* Fixed the default prior on the NIG/GAL `nu` parameter, which was strongly
  biasing every non-Gaussian latent field on a fine mesh. The prior is now 
  mesh dependent so that it is equally uninformative independently of mesh width.
* Increased the burnin for the Gibbs chain for the posterior samples of `W` and `V` 
  that `ngme()` returns from 1 to 100, and before a posterior `simulate()` draw. 
* Default random seeds are now drawn from the ambient R random number stream
  instead of the system clock, so `set.seed()` makes `simulate()`,
  `control_opt()`, `ngme()`, `predict()`, `cross_validation()`,
  `compute_ngme_ci()`, `compute_ngme_sgld_samples()` and `test_ngme()`
  reproducible in the usual R way. 
* Restore AR1 model to use stationary initial conditions and fix all models with 
  triangular operators to use exact trace calculations. 

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
