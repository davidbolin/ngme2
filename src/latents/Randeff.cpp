// // -------------- random effects implementation ----------------
// Randeff::Randeff(const Rcpp::List& R_randeff, unsigned long seed) :
//   randeff_rng      (seed),
//   effect_type      (Rcpp::as<string>        (R_randeff["effect_type"])),
//   B_reff           (Rcpp::as<MatrixXd>      (R_randeff["B_reff"])),
//   Sigma            (Rcpp::as<MatrixXd>      (R_randeff["Sigma"])),
//   invSigma         (Sigma.inverse()),
//   n_reff           (B_reff.cols()),
//   n_cov_params     (n_reff * (n_reff + 1) / 2),
//   n_params         (n_cov_params + n_reff + 1), // mu, Sigma_vech, nu
//   U                (VectorXd::Zero(n_reff)),
//   mu               (VectorXd::Zero(n_reff)),
//   Dd               (duplicatematrix(n_reff)),
//   var              (Rcpp::as<Rcpp::List> (R_randeff["mix_var"]), randeff_rng())
// {}

// VectorXd Randeff::precond_grad() {
//   VectorXd g = VectorXd::Zero(n_params);

//   // gradient of mu (of length n_reff)
//   VectorXd V = var.getV();
//   VectorXd g_mu = ((-1 + V(0)) / V(0)) * (invSigma * U);

//   // gradient of Sigma_vech (of length n_cov_params)
//   MatrixXd iSkroniS = kroneckerProduct(invSigma, invSigma);
//   VectorXd UUt = vec(U * U.transpose());
//   VectorXd dSigma_vech = 0.5 * Dd.transpose() * iSkroniS  * (UUt - vec(Sigma));

//   g << g_mu, dSigma_vech, var.grad_log_nu();
//   return g;
// }

// VectorXd Randeff::get_parameter() {
//   VectorXd p = VectorXd::Zero(n_params);
//   p << mu, vech(Sigma), var.get_nu();
//   return p;
// }

// void Randeff::set_parameter(const VectorXd& p) {
//   mu = p.head(n_reff);
//   VectorXd vech = p.segment(n_reff, n_cov_params);
//   VectorXd v = Dd * vech;
//   Sigma = veci(v, n_reff, n_reff);
//   invSigma = Sigma.inverse();
//   var.set_log_nu(p.tail(1)(0));

//   // sample U|V ~ N(0, 1/V * Sigma^-1)
// }
