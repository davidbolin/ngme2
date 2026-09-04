# ngme2 (development version)

* Remove `moving_window` from `ngme()`. It was declared and documented ("return
  the average estimation of last .. iterations") but never read anywhere in the
  function, so it silently did nothing. 
* `spacetime()` gains `cc_variance_free`, DEFAULT TRUE. The operator now puts
  the whole `cc` factor on the temporal term (`K = cc*BtCs + Ls`) instead of
  splitting it `sqrt(cc)` between the temporal and spatial terms. This makes the
  variance less dependent on `cc` and typically improves mixing. 
* `cross_validation()` accepts a NAMED list of `transform`s and scores the same
  posterior draws on each scale. `transform` is applied after the draws exist
  and never affects sampling, so calling once per scale re-drew identical
  samples and doubled the cost of a cross-validation. Verified to reproduce the
  two-call scores to 10 significant figures. An unnamed list keeps its previous
  meaning of one entry per model.
* Fix `n_burnin` being silently dropped in `cross_validation()`. 
* Live diagnostics (`[iter ...]`, `[stationarity]`, `[polish]`) now go to
  stderr via `REprintf` rather than stdout. 
* The polish phase's geometric decay now starts from the step size the search
  was using. `polish_stepsize_factor` defaults to 1.0.
* Trajectory continuation now requires the same model.
* Parallelise the gradient computation over replicates. 
* Fix the advection term of the non-separable `spacetime()` model, which was
  quadratic in the advection field and simply wrong. This changes fitted results 
  for every non-separable model with free advection.
* Fix `spacetime()` assembling a different operator in C++ than in R. 
* `spacetime()` now computes its operator traces from the diagonal blocks,
  needing neither Hutchinson probes nor a factorization of the full
  `(nt*ns) x (nt*ns)` operator. 
* The `spacetime()` operator now starts from its approximatestationary distribution 
  instead of in zero. `spacetime(stationary_init = FALSE)` restores the old first block for
  reproducing earlier fits.
* Convergence is now decided in two phases, which separates two questions the
  old rule ran together. `control_opt(conv_criterion = )` selects between them
  and defaults to the new `"stationarity"`.
  - The SEARCH phase stops as soon as the chains have settles: a Geweke
    split-half comparison of the running mean, per parameter, computed on the
    second half of the history so the burn-in transient leaves the comparison
    window as the run proceeds.
  - The POLISH phase then reduces the step size geometrically and averages the 
    iterates (Polyak-Ruppert, which was already what the reported estimate is) 
    and runs until that average is precise. it stops when the Monte Carlo error
    of the estimate falls below `mc_se_lim * se_stat`, i.e. when the Monte Carlo
    error is a fraction of the statistical error the data can support.
* Very short runs now report non-convergence rather than claiming it. The
  stationarity test needs enough sub-batch history to split into two halves, so
  a fit given only a few dozen iterations warns instead of declaring success. 
* The stationarity test's comparison window is derived from the measured
  correlation length of the chains rather than being a fixed count
  (`control_opt(stationarity_window = 0, stationarity_eff = )`). tau is
  estimated as the ratio of the batch-means variance to the naive variance of
  the same mean and the window is `2 * stationarity_eff * tau`. 
* The optional trend-triggered step-size decay
  (`stepsize_decay(method = "on_trend")`) used the same |t| <= `trend_lim` test
  to decide the chains had stopped moving, and so would decay soonest when the
  window was noisiest. It now uses the drift rate as well.
* `ngme(start = previous_fit)` is now a true continuation. `start_sd` is no
  longer applied on a warm start at all. On top of that, each chain resumes
  from its own final parameters instead of from the chain-averaged fit. 
  Each chain's endpoint is now stored on the fit (`attr(fit, "chain_params")`) 
  and handed back to the optimizer, and the storedoptimisation trajectory is 
  concatenated across the restarts rather than replaced, so `traceplot()` and any 
  convergence diagnostic see the whole path.
  Each chain's own latent state (`W`, and the mixing variables `V`) is carried
  over too: `start` holds the chain average of those.
  `control_opt(continue_chains = FALSE)` restores the previous behaviour.
* Convergence gained a third clause: the estimate must be precise, not merely
  stationary. R-hat and the trend test say the chains reached the same stationary
  distribution and stopped drifting; neither says the answer is pinned down, because
  constant-step SGD settles into a distribution whose width is set by the step size
  rather than by the data. `control_opt(mc_se_conv_check = , mc_se_lim = )` requires the
  Monte Carlo error of the reported estimate to be at most `mc_se_lim` (default 0.5) of
  `se_stat = sqrt(diag(P^-1))`, the standard error the data supports, read from the
  preconditioner.
* Measurement-noise parameter names are prefixed `meas_`. A latent noise and the
  measurement noise both have `mu`/`sigma`/`nu`, so `par_names` could contain `sigma_1`
  twice and anything keyed by name silently merged two unrelated parameters.
* R-hat is now computed over the same window as the trend test (sub-batch means
  spanning the second half of the run) instead of over the current batch alone, and
  `n_conv_batch` defaults to 2 rather than 1. 
* The convergence trend test gained a magnitude floor, `control_opt(trend_rel_lim = )`
  (default 0.01). A parameter now fails the test only when its drift is both
  statistically detectable (`|slope|/se > trend_lim`) and materially fast, moving more
  than `trend_rel_lim` of its own scale (`|mean| + sd`) per 100 iterations. It is a rate,
  not a movement per regression window.
* The trend window is now built from sub-batch means spread evenly over the second half
  of the run so far, instead of the last `n_slope_check` whole checkpoints. It is
  therefore available at the first checkpoint rather than the `n_slope_check`-th. 
  `n_min_batch` now defaults to 1 rather than `min(n_batch, 3)`, since it is no 
  longer what holds a converged fit back.
* With `verbose = TRUE` and `print_check_info = FALSE`, each convergence checkpoint now
  prints a progress summary: how many parameters have converged, and the worst
  R-hat, t-statistic and drift rate with the parameter responsible. 
* `predict()` gained a `replicate` argument. Every prediction previously came
  from the first replicate whatever the caller asked for, since the latent field
  `W` is per-replicate. The argument takes a replicate level (the name used in
  the fit) or a position, and an unmatched value is an error listing what the
  fit actually has.
* A namespace-qualified `ngme2::f(...)` in a model formula is now recognised as
  a latent term. `terms.formula(specials = )` matches specials by name, so the
  qualified call fell through to the fixed-effect path and the latent field was
  dropped from the model silently. The qualification is stripped from the
  function position only, and the bare name is bound for evaluation, so this
  keeps working when `ngme2` is not attached.
* Fixed `f(model = ngme2::<constructor>(...), noise = noise_normal())`, which
  failed with "the condition has length > 1". `f()` read the constructor name as
  `as.character(model_expr[[1]])`, which is `c("::", "ngme2", "tp")` for a
  qualified call, so the `%in%` test that follows was not a single logical.
* Removed the dead preconditioner plumbing from `control_opt()`
  (`precond_strategy`, `precond_by_diff_chain`, `compute_precond_each_iter`).
  `precond_sgd()` has not carried a `preconditioner` field since the
  preconditioner became the analytic second-order Hessian, so `control_opt()`
  overwrote its own default with `NULL` and passed `precond_strategy =
  numeric(0)` to C++ -- which never read it. `numerical_eps` is unaffected.
* The default `nu` prior for stationary NIG latent noise is now the penalised-complexity
  prior of Cabral, Bolin and Rue (2023) calibrated as in their Section 3.3, from
  `P(1/nu > U) = alpha` with `U = 2.5` and `alpha = 0.01`, i.e. `P(nu < 0.4) = 0.01`.
* A failed Cholesky of the measurement precision `QQ` no longer aborts the fit outright.
  Such a failure is nearly always a transient SGD excursion so the diagonal is now
  nudged by an escalating amount (from 1e-10 of the diagonal scale, six attempts) and the
  factorization retried. Only if that fails does estimation stop, and the error now reports
  the iteration, the smallest diagonal entry and the asymmetry instead of printing them to
  the console from inside the parallel SGD region.
* `sf` and `inlabru` are declared in `Suggests`; `spacetime()` now checks for `sf`
  up front rather than failing inside `sf::st_coordinates()`.
* The `Parana` section of the introductory vignette is skipped when `INLA` cannot be loaded,
  instead of failing the whole vignette.
* Rao-Blackwellisation is now on by default (`control_opt(rao_blackwellization = )`).
* Fixed the Gibbs sweep order that made Rao-Blackwellisation unusable with a non-Gaussian
  `family`. The measurement `V` was redrawn after `W`, which left `QQ` and `cond_W`
  conditioned on the previous `V` while the gradient used the new one; the fit diverged until
  `QQ` lost positive-definiteness. Both `V` blocks are now drawn before `W`. Gaussian
  families were unaffected, since there `V` is fixed.
* The Rao-Blackwell traces are taken from an exact selected (Takahashi) inverse when the
  Cholesky factor of `QQ` is cheap enough, and from Hutchinson probes otherwise, controlled
  by `control_opt(selinv_max_fill = )`. Triangular operators (`ar1`, `ou`, `arma`) and 1d
  meshes qualify; a 2d mesh does not. For those models the traces are now exact and the fit
  no longer depends on `n_trace_iter` at all.
* The Hutchinson probe count adapts during the run (`control_opt(trace_adapt = )`) so that
  the trace estimator contributes a target share of the gradient variance instead of being
  fixed a priori. 
* The convergence trend test is now scale-free: it compares the fitted slope to its own
  standard error rather than to an absolute bound, so `trend_lim` is a t-value (default 2)
  and no longer depends on the parameter's units or on `n_batch`. 
* A run that exhausts its iteration budget without converging now warns and names the
  parameters still failing, instead of returning silently. Set
  `control_opt(warn_no_convergence = FALSE)` for short runs where convergence is not expected.
* Removed the Pflug convergence diagnostic (`pflug_conv_check`, `pflug_alpha`). It could
  never fire: the rule looks for sign changes in the running sum of `<g_t, g_{t-1}>`, but the
  gradient is evaluated on a Gibbs sample of the latent field, so consecutive gradients stay
  positively correlated however close the optimum is. 
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
* Fixed the trace estimators for non-symmetric operators, which made the
  preconditioner Hessian incorrect and the first-order trace far noisier than
  necessary. `control_opt(nonsym_solver = )` chooses the factorization for this case.
  `"normal_equations"` (default) for a Cholesky of `t(K) K`, or `"lu"` for a sparse LU 
  of `K`. The Cholesky is cheaper for every non-symmetric operator in the package. 
  Use `"lu"` if `K` is so ill-conditioned that the Cholesky of `t(K) K` fails.
* A failed factorization of `K` for a non-symmetric operator is now reported for
  both routes, not only for the LU. 
* Speed-ups for space-time models. `spacetime()` builds its spatial block once per
  assembly of `K` rather than once per time node: with the default `B_gamma_x` /
  `B_gamma_y` the advection field is the same at every time node, so all `nt - 1`
  blocks are the same matrix. The Rao-Blackwell traces no longer assemble the
  matrix they trace either: `tr(QQ^-1 A' diag(d) B)` is now taken as
  `A' (d * (B QU))` on the probe block instead of forming the triple product,
  which is as dense as `QQ` and was built once per parameter per iteration.
  `NGME2_NO_FACTORED_TRACE` restores the assembled trace for comparison.
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
