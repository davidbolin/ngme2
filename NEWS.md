# ngme2 0.8.3 (2025-12-10)
* Update nu parameterization to be relative to nu_lower_bound (nu = nu_lower_bound + exp(B_nu * theta_nu))
* Allow NULL specification for noise parameter in non-stationary case
* Update fixed effects initialization
* Update start=fit logic
* Add Eigen SparseQR backend option (solver_backend = "eigen", solver_type = "qr")

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
