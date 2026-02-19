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
