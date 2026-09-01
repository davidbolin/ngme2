#include "latent.h"
#include "include/thread_io.h"
#include "prior.h"
#include <algorithm>
#include <chrono>
#include <stdexcept>

namespace {
std::string parse_prior_target(const Rcpp::List &prior_list,
                               const std::string &default_target = "coef") {
  std::string target = prior_list.containsElementNamed("target")
                           ? Rcpp::as<std::string>(prior_list["target"])
                           : default_target;
  if (target != "coef" && target != "field") {
    throw std::invalid_argument("prior target must be 'coef' or 'field'");
  }
  return target;
}

VectorXd prior_score_vec(const std::string &type, const VectorXd &param,
                         const VectorXd &x) {
  VectorXd out(x.size());
  for (int i = 0; i < x.size(); ++i) {
    out(i) = PriorUtil::d_log_dens(type, param, x(i));
  }
  return out;
}

void parse_prior_spec(const Rcpp::List &prior_list, std::string &type,
                      VectorXd &param, std::string &target) {
  type = Rcpp::as<std::string>(prior_list["type"]);
  param = Rcpp::as<VectorXd>(prior_list["param"]);
  target = parse_prior_target(prior_list);
}
} // namespace

// get theta[unfixed]
VectorXd get_parameter_unfixed(const VectorXd &theta,
                               const vector<bool> &fix_theta) {
  // Count number of unfixed parameters
  int n_unfixed = 0;
  for (size_t i = 0; i < fix_theta.size(); i++) {
    if (!fix_theta[i])
      n_unfixed++;
  }

  // Create vector of unfixed parameters
  VectorXd unfixed = VectorXd::Zero(n_unfixed);
  int j = 0;
  for (size_t i = 0; i < fix_theta.size(); i++) {
    if (!fix_theta[i]) {
      unfixed(j++) = theta(i);
    }
  }
  return unfixed;
}

// set theta[unfixed] = new_theta
void set_parameter_unfixed(VectorXd &theta, const vector<bool> &fix_theta,
                           const VectorXd &new_theta) {
  int j = 0;
  for (size_t i = 0; i < fix_theta.size(); i++) {
    if (!fix_theta[i]) {
      theta(i) = new_theta(j++);
    }
  }
}

// K is V_size * W_size matrix
Latent::Latent(const Rcpp::List &model_list, unsigned long seed)
    : latent_rng(seed), model_type(Rcpp::as<string>(model_list["model"])),
      noise_type(Rcpp::as<vector<string>>(model_list["noise_type"])),
      debug(Rcpp::as<bool>(model_list["debug"])), n_noise(noise_type.size()),
      W_size(Rcpp::as<int>(model_list["W_size"])),
      V_size(Rcpp::as<int>(model_list["V_size"])),
      n_params(Rcpp::as<int>(model_list["n_params"])),

      // operator
      ope(OperatorFactory::create(
          Rcpp::as<Rcpp::List>(model_list["operator"]))),
      h(ope->get_h()), theta_K(Rcpp::as<VectorXd>(model_list["theta_K"])),
      n_theta_K(theta_K.size()), symmetricK(ope->is_symmetric()),
      zero_trace(ope->is_zero_trace()),

      // K             (V_size, W_size),
      // dK            (n_theta_K),
      trace(n_theta_K, 0.0), eps(0.001),

      W(W_size), prevW(W_size), cond_W(W_size), V(h), prevV(h),
      A(Rcpp::as<SparseMatrix<double, 0, int>>(model_list["A"])),

      p_vec(V_size), a_vec(V_size), b_vec(V_size) {
  if (debug)
    ngme_io::out() << "begin constructor of latent" << std::endl;
  assert(W_size == V_size);

  // construct from ngme_noise
  fix_flag[latent_fix_theta_K] = Rcpp::as<bool>(model_list["fix_theta_K"]);
  Rcpp::List noise_in = Rcpp::as<Rcpp::List>(model_list["noise"]);
  fix_flag[latent_fix_theta_mu] = Rcpp::as<bool>(noise_in["fix_theta_mu"]);
  // Handle vector-based fix_theta_sigma
  Rcpp::LogicalVector fix_theta_sigma_r =
      Rcpp::as<Rcpp::LogicalVector>(noise_in["fix_theta_sigma"]);
  fix_theta_sigma_vec =
      std::vector<bool>(fix_theta_sigma_r.begin(), fix_theta_sigma_r.end());
  fix_flag[latent_fix_theta_nu] = noise_in.containsElementNamed("fix_theta_nu")
                                      ? Rcpp::as<bool>(noise_in["fix_theta_nu"])
                                      : false;

  single_V = Rcpp::as<bool>(noise_in["single_V"]);

  B_mu = Rcpp::as<MatrixXd>(noise_in["B_mu"]);
  B_sigma = Rcpp::as<MatrixXd>(noise_in["B_sigma"]);
  B_nu = Rcpp::as<MatrixXd>(noise_in["B_nu"]);

  n_theta_mu = Rcpp::as<int>(noise_in["n_theta_mu"]);
  n_theta_sigma = Rcpp::as<int>(noise_in["n_theta_sigma"]);
  n_theta_nu = Rcpp::as<int>(noise_in["n_theta_nu"]);
  nu_lower_bound = noise_in.containsElementNamed("nu_lower_bound")
                       ? Rcpp::as<double>(noise_in["nu_lower_bound"])
                       : 0.0;

  theta_mu = Rcpp::as<VectorXd>(noise_in["theta_mu"]);
  theta_sigma = Rcpp::as<VectorXd>(noise_in["theta_sigma"]);
  theta_nu = Rcpp::as<VectorXd>(noise_in["theta_nu"]);

  share_V = Rcpp::as<bool>(noise_in["share_V"]);
  single_V = Rcpp::as<bool>(noise_in["single_V"]);

  // init priors for parameter of K and noise
  Rcpp::List prior_k_list = Rcpp::as<Rcpp::List>(model_list["prior_theta_K"]);
  prior_K_type.resize(n_theta_K, "none");
  prior_K_param.resize(n_theta_K, VectorXd::Zero(0));
  prior_K_target.resize(n_theta_K, "coef");
  if (n_theta_K > 0 && prior_k_list.size() > 0) {
    if (prior_k_list.containsElementNamed("type")) {
      std::string t;
      VectorXd p;
      std::string target;
      parse_prior_spec(prior_k_list, t, p, target);
      for (int i = 0; i < n_theta_K; ++i) {
        prior_K_type[i] = t;
        prior_K_param[i] = p;
        prior_K_target[i] = target;
      }
    } else {
      if (prior_k_list.size() != n_theta_K) {
        throw std::invalid_argument(
            "prior_theta_K must have length equal to n_theta_K");
      }
      for (int i = 0; i < n_theta_K; ++i) {
        Rcpp::List one = Rcpp::as<Rcpp::List>(prior_k_list[i]);
        parse_prior_spec(one, prior_K_type[i], prior_K_param[i],
                         prior_K_target[i]);
      }
    }
  }

  for (int i = 0; i < n_theta_K; ++i) {
    if (prior_K_target[i] != "coef") {
      throw std::invalid_argument(
          "theta_K prior target='field' is not supported; use target='coef'");
    }
  }

  Rcpp::List prior_list = Rcpp::as<Rcpp::List>(noise_in["prior_mu"]);
  parse_prior_spec(prior_list, prior_mu_type, prior_mu_param, prior_mu_target);
  prior_list = Rcpp::as<Rcpp::List>(noise_in["prior_sigma"]);
  parse_prior_spec(prior_list, prior_sigma_type, prior_sigma_param,
                   prior_sigma_target);
  prior_list = Rcpp::as<Rcpp::List>(noise_in["prior_nu"]);
  parse_prior_spec(prior_list, prior_nu_type, prior_nu_param, prior_nu_target);

  // init W
  if (model_list["W"] != R_NilValue) {
    W = Rcpp::as<VectorXd>(model_list["W"]);
  } else {
    W = VectorXd::Zero(W_size);
  }
  prevW = W;
  fix_flag[latent_fix_W] = Rcpp::as<bool>(model_list["fix_W"]); // fixW

  // init K, Q, dK
  int n_trace_iter = Rcpp::as<int>(model_list["n_trace_iter"]);
  int solver_type = Rcpp::as<int>(model_list["solver_type"]);
  n_trace_iter_ = n_trace_iter;
  solver_type_ = solver_type;
  nonsym_solver_ = model_list.containsElementNamed("nonsym_solver")
                       ? Rcpp::as<int>(model_list["nonsym_solver"])
                       : 0;
  robust_ = model_list.containsElementNamed("robust")
                ? Rcpp::as<bool>(model_list["robust"])
                : false;

  // build mu, sigma, compute trace, ...
  update_each_iter(true);

  // Initialize V
  if (!Rf_isNull(noise_in["V"])) {
    V = Rcpp::as<VectorXd>(noise_in["V"]);
    prevV = V;
  } else {
    // V=h at init.
    // to-fix (for normal noise case)
    // sample_uncond_V(); // sample_uncond_V();
    for (int i = 0; i < noise_type.size(); i++) {
      int n = V_size / n_noise;
      if (noise_type[i] == "normal")
        V.segment(i * n, n) = h.segment(i * n, n);
    }
  }

  if (debug)
    ngme_io::out() << "End constructor of latent" << std::endl;
  last_gradient_ = VectorXd::Zero(n_params);
  last_precond_ = MatrixXd::Zero(n_params, n_params);
  invalidate_derivatives();
}

void Latent::invalidate_derivatives() {
  deriv_cache.initialized = false;
  deriv_cache.grad_K_ready = deriv_cache.grad_mu_ready = false;
  deriv_cache.grad_sigma_ready = deriv_cache.grad_nu_ready = false;
  hess_cache.ready = false;
  grad_cache_valid_ = false;
  precond_cache_valid_ = false;
}


Rcpp::List Latent::output() const {
  return Rcpp::List::create(
      Rcpp::Named("model") = model_type, Rcpp::Named("noise_type") = noise_type,
      Rcpp::Named("theta_K") = theta_K, // same parameterization as input
      Rcpp::Named("theta_mu") = theta_mu,
      Rcpp::Named("theta_sigma") = theta_sigma,
      Rcpp::Named("theta_nu") = theta_nu, Rcpp::Named("V") = V,
      Rcpp::Named("W") = W);
}

const VectorXd Latent::get_parameter() {
  VectorXd parameter = VectorXd::Zero(n_params);

  if (!fix_flag[latent_fix_theta_K])
    parameter.segment(0, n_theta_K) = theta_K;
  if (!fix_flag[latent_fix_theta_mu])
    parameter.segment(n_theta_K, n_theta_mu) = theta_mu;
  // Handle vector-based fixing for theta_sigma
  parameter.segment(n_theta_K + n_theta_mu, n_theta_sigma) =
      get_parameter_unfixed(theta_sigma, fix_theta_sigma_vec);
  if (!fix_flag[latent_fix_theta_nu])
    parameter.segment(n_theta_K + n_theta_mu + n_theta_sigma, n_theta_nu) =
        theta_nu;

  if (nig_std_active())
    parameter.segment(n_theta_K, 3) = nig_std_from_native();

  if (debug) {
    if (std::isnan(parameter(0)) || std::isnan(-parameter(0)))
      throw std::runtime_error("isnan");
    ngme_io::out() << "parameter= " << parameter << std::endl;
  }

  return parameter;
}

const VectorXd Latent::get_grad() {
  if (!grad_cache_valid_) {
    throw std::runtime_error(
        "Latent::get_grad() called before compute_grad_and_hessian()");
  }
  return last_gradient_;
}

void Latent::compute_grad_and_hessian(bool rao_blackwell, bool with_precond) {
  if (!state_ready_) {
    throw std::runtime_error(
        "Latent::compute_grad_and_hessian() called before update_each_iter()");
  }
  if (with_precond && !state_has_precond_terms_) {
    throw std::runtime_error("Latent preconditioner requested but second-order "
                             "terms not prepared; call update_each_iter(true)");
  }
  bool need_grad =
      (!grad_cache_valid_) || (grad_cache_rb_mode_ != rao_blackwell);
  bool need_precond = with_precond && !precond_cache_valid_;
  if (!need_grad && !need_precond)
    return;

  if (need_grad) {
    if (debug)
      ngme_io::out() << "Start latent gradient compute" << std::endl;
    VectorXd grad = VectorXd::Zero(n_params);

    bool need_K = !fix_flag[latent_fix_theta_K];
    bool need_mu = (n_theta_mu > 0) && !fix_flag[latent_fix_theta_mu];
    bool need_sigma = n_theta_sigma > 0;
    bool need_nu = !fix_flag[latent_fix_theta_nu];

    update_derivatives(need_K, need_mu, need_sigma, need_nu, rao_blackwell);

    if (need_K)
      grad.segment(0, n_theta_K) = deriv_cache.grad_theta_K;
    if (need_mu)
      grad.segment(n_theta_K, n_theta_mu) = deriv_cache.grad_theta_mu;
    if (n_theta_sigma > 0)
      grad.segment(n_theta_K + n_theta_mu, n_theta_sigma) =
          deriv_cache.grad_theta_sigma;
    if (need_nu)
      grad.segment(n_theta_K + n_theta_mu + n_theta_sigma, n_theta_nu) =
          deriv_cache.grad_theta_nu;

    last_gradient_ = grad;
    grad_cache_valid_ = true;
    grad_cache_rb_mode_ = rao_blackwell;
    if (debug)
      ngme_io::out() << "finish latent gradient" << std::endl;
  }

  if (need_precond) {
    auto t_start = std::chrono::steady_clock::now();

    // Build analytic Hessian blocks (state already contains necessary
    // derivatives)
    compute_hessian_blocks(false);

    MatrixXd precond_full = MatrixXd::Zero(n_params, n_params);
    const int off_K = 0;
    const int off_mu = n_theta_K;
    const int off_sigma = n_theta_K + n_theta_mu;
    const int off_nu = n_theta_K + n_theta_mu + n_theta_sigma;

    if (n_theta_K > 0)
      precond_full.block(off_K, off_K, n_theta_K, n_theta_K) = hess_cache.H_K;
    if (n_theta_K > 0 && n_theta_mu > 0)
      precond_full.block(off_K, off_mu, n_theta_K,
                         n_theta_mu) = hess_cache.H_K_mu,
                         precond_full.block(off_mu, off_K, n_theta_mu,
                                            n_theta_K) =
                             hess_cache.H_K_mu.transpose();
    if (n_theta_K > 0 && n_theta_sigma > 0)
      precond_full.block(off_K, off_sigma, n_theta_K,
                         n_theta_sigma) = hess_cache.H_K_sigma,
                         precond_full.block(off_sigma, off_K, n_theta_sigma,
                                            n_theta_K) =
                             hess_cache.H_K_sigma.transpose();
    if (n_theta_mu > 0)
      precond_full.block(off_mu, off_mu, n_theta_mu, n_theta_mu) =
          hess_cache.H_mu;
    if (n_theta_sigma > 0)
      precond_full.block(off_sigma, off_sigma, n_theta_sigma, n_theta_sigma) =
          hess_cache.H_sigma;
    if (n_theta_mu > 0 && n_theta_sigma > 0)
      precond_full.block(off_mu, off_sigma, n_theta_mu,
                         n_theta_sigma) = hess_cache.H_mu_sigma,
                         precond_full.block(off_sigma, off_mu, n_theta_sigma,
                                            n_theta_mu) =
                             hess_cache.H_mu_sigma.transpose();
    if (n_theta_nu > 0)
      precond_full.block(off_nu, off_nu, n_theta_nu, n_theta_nu) =
          hess_cache.H_nu;
    if (nig_std_active()) {
      MatrixXd J = nig_std_jacobian();
      // noise-noise block
      precond_full.block(off_mu, off_mu, 3, 3) =
          J.transpose() * precond_full.block(off_mu, off_mu, 3, 3).eval() * J;
      // cross blocks with theta_K: J on the noise side only
      if (n_theta_K > 0) {
        precond_full.block(off_K, off_mu, n_theta_K, 3) =
            precond_full.block(off_K, off_mu, n_theta_K, 3).eval() * J;
        precond_full.block(off_mu, off_K, 3, n_theta_K) =
            precond_full.block(off_K, off_mu, n_theta_K, 3).transpose();
      }
    }

    // Jitter to ensure PD in placeholder form.
    // These blocks hold the Hessian, and Ngme::precond() returns its negative
    // as the information matrix that the LLT solve factorises. Ridging the
    // information therefore means SUBTRACTING here.
    precond_full.diagonal().array() -= 1e-5;

    last_precond_ = precond_full;
    precond_cache_valid_ = true;

    if (debug) {
      auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                    std::chrono::steady_clock::now() - t_start)
                    .count();
      ngme_io::out() << "[latent] compute_precond_matrix (analytic skeleton) "
                     "timing (ms): total="
                  << ms << std::endl;
    }
  }
}

void Latent::set_parameter_and_update(const VectorXd &theta,
                                      bool with_precond) {
  // The optimiser may be working in the standardised coordinates; everything
  // below (and every derivative) stays native, so convert on the way in.
  if (nig_std_active()) {
    double tm, ts, tn;
    nig_std_to_native(theta.segment(n_theta_K, 3), tm, ts, tn);
    if (!fix_flag[latent_fix_theta_K])
      theta_K = theta.segment(0, n_theta_K);
    theta_mu(0) = tm;
    theta_sigma(0) = ts;
    theta_nu(0) = tn;
    state_ready_ = false;
    state_has_precond_terms_ = false;
    invalidate_derivatives();
    update_each_iter(with_precond);
    return;
  }
  // nig, gal and normal+nig
  if (!fix_flag[latent_fix_theta_K])
    theta_K = theta.segment(0, n_theta_K);
  if (!fix_flag[latent_fix_theta_mu])
    theta_mu = theta.segment(n_theta_K, n_theta_mu);
  // Handle vector-based fixing for theta_sigma setting
  set_parameter_unfixed(theta_sigma, fix_theta_sigma_vec,
                        theta.segment(n_theta_K + n_theta_mu, n_theta_sigma));
  if (!fix_flag[latent_fix_theta_nu])
    theta_nu =
        theta.segment(n_theta_K + n_theta_mu + n_theta_sigma, n_theta_nu);

  state_ready_ = false;
  state_has_precond_terms_ = false;
  invalidate_derivatives();

  update_each_iter(with_precond);
}

void Latent::sample_cond_V() {
  if (fix_flag[latent_fix_V])
    return;
  // update b_inc (p,a_inc already built)
  b_inc = (getK() * W + mu.cwiseProduct(h)).cwiseQuotient(sigma).array().pow(2);

  int n = V_size / n_noise; // n is equal to h_size (of each mesh)
  // sample conditional V
  if (single_V) {
    if (share_V) {
      // type-G1 model (n_nu == 2, but keep same nu)
      double v1 = rGIG_cpp(-0.5 - n, nu[0] + a_inc.dot(h),
                           nu[0] + b_inc.dot(h.cwiseInverse()), latent_rng());
      V = v1 * h;
    } else {
      // type-G2 (n_nu==2) and also univariate single noise (n_nu==1)
      for (int i = 0; i < n_noise; i++) {
        if (noise_type[i] == "normal")
          continue;

        double v1 =
            rGIG_cpp(-0.5 - n * 1.0 / n_noise,
                     nu[i] + a_inc.segment(i * n, n).dot(h.segment(i * n, n)),
                     nu[i] + b_inc.segment(i * n, n).dot(
                                 h.segment(i * n, n).cwiseInverse()),
                     latent_rng());
        V.segment(i * n, n) = v1 * h.segment(i * n, n);
      }
    }
  } else {
    if (share_V) {
      // type-G3 (n_nu == 2, but keep same nu)
      // a_vec = rep(nu, n), b_vec = nu * h
      V.head(n) =
          rGIG_cpp(VectorXd::Constant(n, -1.5),
                   a_vec.head(n) + a_inc.head(n) + a_inc.tail(n),
                   b_vec.head(n) + b_inc.head(n) + b_inc.tail(n), latent_rng());
      V.tail(n) = V.head(n);
    } else {
      // type-G4 (n_nu==2) and also univariate case general noise (n_nu==1)
      for (int i = 0; i < n_noise; i++) {
        if (noise_type[i] == "normal")
          continue;

        // sample conditional V
        V.segment(i * n, n) =
            rGIG_cpp((p_vec + p_inc).segment(i * n, n),
                     (a_vec + a_inc).segment(i * n, n),
                     (b_vec + b_inc).segment(i * n, n), latent_rng());
      }

      // consider the case of normal + nig/gal
      if (noise_type[0] == "normal_nig" && n_noise == 1) {
        // Set the 2nd half of V to be 1
        V.segment(V_size / 2, V_size / 2) = VectorXd::Ones(V_size / 2);
      }
    }
  }

  invalidate_derivatives();
}

void Latent::sample_uncond_V() {
  if (fix_flag[latent_fix_V])
    return;
  prevV = V;
  int n = V_size / n_noise;

  // same logic as in simulation.R
  if (single_V) {
    for (int i = 0; i < n_noise; i++) {
      double v;
      if (noise_type[i] == "nig" || noise_type[i] == "normal_nig")
        v = rGIG_cpp(-0.5, nu[i], nu[i], latent_rng());
      else if (noise_type[i] == "gal")
        v = rGIG_cpp(h[i] * nu[i], 2 * nu[i], 1e-14, latent_rng());
      else if (noise_type[i] == "t")
        v = rGIG_cpp(-nu[i] / 2, 1e-14, nu[i], latent_rng());
      V.segment(i * n, n) = v * h.segment(i * n, n);
    }
  } else {
    for (int i = 0; i < n_noise; i++) {
      // sample unconditional V
      V.segment(i * n, n) =
          rGIG_cpp(p_vec.segment(i * n, n), a_vec.segment(i * n, n),
                   b_vec.segment(i * n, n), latent_rng());
    }
  }

  // keep both V same for bivariate noise with share_V case
  if (n_theta_nu == 2 && share_V) {
    V.segment(n, n) = V.segment(0, n);
  }
  invalidate_derivatives();
}

void Latent::update_derivatives(bool need_grad_theta_K, bool need_grad_theta_mu,
                                bool need_grad_theta_sigma,
                                bool need_grad_theta_nu, bool rao_blackwell) {
  bool reset = (!deriv_cache.initialized) ||
               (deriv_cache.rao_blackwell != rao_blackwell);
  if (reset) {
    deriv_cache.initialized = true;
    deriv_cache.rao_blackwell = rao_blackwell;
    deriv_cache.grad_K_ready = deriv_cache.grad_mu_ready = false;
    deriv_cache.grad_sigma_ready = deriv_cache.grad_nu_ready = false;
  }

  if (need_grad_theta_K)
    compute_theta_K(true, rao_blackwell);
  if (need_grad_theta_mu)
    compute_theta_mu(true, rao_blackwell);
  if (need_grad_theta_sigma)
    compute_theta_sigma(true, rao_blackwell);
  if (need_grad_theta_nu)
    compute_theta_nu(true);
}

void Latent::compute_theta_K(bool need_grad, bool rao_blackwell) {
  if (!need_grad || deriv_cache.grad_K_ready)
    return;
  auto t_total_start = std::chrono::steady_clock::now();
  long long t_chain_ms = 0;
  VectorXd grad = VectorXd::Zero(n_theta_K);
  if (!fix_flag[latent_fix_theta_K]) {
    VectorXd WW = (rao_blackwell) ? cond_W : W;
    VectorXd SV = sigma.array().pow(2).matrix().cwiseProduct(V);
    VectorXd tmp = (getK() * WW) - mu.cwiseProduct(V - h);
    for (int j = 0; j < n_theta_K; j++) {
      grad(j) =
          trace[j] - (get_dK(j) * WW).cwiseProduct(SV.cwiseInverse()).dot(tmp);
    }

    for (int l = 0; l < n_theta_K; l++) {
      grad(l) +=
          PriorUtil::d_log_dens(prior_K_type[l], prior_K_param[l], theta_K(l));
    }

    // apply per-parameter fixing mask from operator (if provided)
    std::vector<bool> kfix = ope->get_fix_mask_K();
    if ((int)kfix.size() == n_theta_K) {
      for (int i = 0; i < n_theta_K; ++i)
        if (kfix[i])
          grad(i) = 0.0;
    }
  }
  if (debug) {
    auto t_total_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                          std::chrono::steady_clock::now() - t_total_start)
                          .count();
    ngme_io::out() << "[latent] compute_theta_K timing (ms): total=" << t_total_ms
                << ", Z-chain=" << t_chain_ms << std::endl;
  }
  deriv_cache.grad_theta_K = grad;
  deriv_cache.grad_K_ready = true;
}

void Latent::compute_theta_mu(bool need_grad, bool rao_blackwell) {
  if (!need_grad || deriv_cache.grad_mu_ready)
    return;

  VectorXd grad = VectorXd::Zero(n_theta_mu);
  if (n_theta_mu > 0 && !fix_flag[latent_fix_theta_mu]) {
    bool purely_normal = (n_noise == 1 && noise_type[0] == "normal") ||
                         (n_noise == 2 && noise_type[0] == "normal" &&
                          noise_type[1] == "normal");
    if (!purely_normal) {
      VectorXd WW = (rao_blackwell) ? cond_W : W;
      // We do not Rao-Blackwellised over V, deliberately. With u = K WW this splits
      // as sigma^-2 [ (u + 2 mu h) - mu V - mu h^2/V - h u/V ]: the middle
      // terms are pure in V and could be replaced by E[.|W], but the cross
      // term h u / V cannot, since cond_W = QQ^-1 M is itself a function of V
      // through a solve. Averaging only part of a sum is not a
      // Rao-Blackwellisation and carries no variance guarantee: here the
      // pure-V terms are negatively correlated with the cross term, so
      // replacing them removes a cancellation and the optimiser wobble on
      // log sigma grows about threefold.
      VectorXd SV = sigma.array().pow(2).matrix().cwiseProduct(V);
      for (int l = 0; l < n_theta_mu; l++) {
        grad(l) = (V - h)
                      .cwiseProduct(B_mu.col(l).cwiseQuotient(SV))
                      .dot(getK() * WW - mu.cwiseProduct(V - h));
      }
    }
    if (prior_mu_target == "coef") {
      for (int l = 0; l < n_theta_mu; l++) {
        grad(l) += PriorUtil::d_log_dens(prior_mu_type, prior_mu_param,
                                         theta_mu(l));
      }
    } else {
      grad += B_mu.transpose() * prior_score_vec(prior_mu_type, prior_mu_param, mu);
    }
  }

  deriv_cache.grad_theta_mu = grad;
  deriv_cache.grad_mu_ready = true;
}

void Latent::compute_theta_sigma(bool need_grad, bool rao_blackwell) {
  if (!need_grad || deriv_cache.grad_sigma_ready)
    return;

  VectorXd grad = VectorXd::Zero(n_theta_sigma);
  bool all_fixed =
      std::all_of(fix_theta_sigma_vec.begin(), fix_theta_sigma_vec.end(),
                  [](bool b) { return b; });
  if (!all_fixed) {
    VectorXd WW = (rao_blackwell) ? cond_W : W;
    // Sampled V here too -- see the note in compute_theta_mu.
    VectorXd SV = sigma.array().pow(2).matrix().cwiseProduct(V);
    VectorXd V_minus_h = V - h;
    VectorXd tmp = (getK() * WW - mu.cwiseProduct(V_minus_h))
                       .array()
                       .pow(2)
                       .matrix()
                       .cwiseQuotient(SV);
    VectorXd all_grad =
        B_sigma.transpose() * (tmp - VectorXd::Ones(tmp.size()));
    int j = 0;
    for (int i = 0; i < all_grad.size(); i++) {
      if (!fix_theta_sigma_vec[i]) {
        grad(j++) = all_grad(i);
      }
    }
    if (prior_sigma_target == "coef") {
      j = 0;
      for (int l = 0; l < theta_sigma.size(); l++) {
        if (!fix_theta_sigma_vec[l]) {
          grad(j++) += PriorUtil::d_log_dens(prior_sigma_type, prior_sigma_param,
                                             theta_sigma(l));
        }
      }
    } else {
      VectorXd sigma_lp = B_sigma * theta_sigma;
      VectorXd prior_full = B_sigma.transpose() *
                            prior_score_vec(prior_sigma_type, prior_sigma_param,
                                            sigma_lp);
      j = 0;
      for (int l = 0; l < prior_full.size(); l++) {
        if (!fix_theta_sigma_vec[l]) {
          grad(j++) += prior_full(l);
        }
      }
    }
  }

  deriv_cache.grad_theta_sigma = grad;
  deriv_cache.grad_sigma_ready = true;
}

void Latent::compute_theta_nu(bool need_grad) {
  if (!need_grad || deriv_cache.grad_nu_ready)
    return;

  VectorXd grad = VectorXd::Zero(n_theta_nu);
  if (n_theta_nu > 0 && !fix_flag[latent_fix_theta_nu]) {
    if (noise_type[0] != "normal") {
      if (n_noise == 1) {
        if (noise_type[0] == "normal_nig") {
          MatrixXd B_nu_half = B_nu.block(0, 0, B_nu.rows() / 2, B_nu.cols());
          VectorXd nu_half = nu.segment(0, nu.size() / 2);
          VectorXd V_half = V.segment(0, V.size() / 2);
          VectorXd prevV_half = prevV.segment(0, prevV.size() / 2);
          VectorXd h_half = h.segment(0, h.size() / 2);
          grad = NoiseUtil::grad_theta_nu(noise_type[0], B_nu_half, nu_half,
                                          V_half, prevV_half, h_half,
                                          nu_lower_bound, single_V);
        } else {
          grad = NoiseUtil::grad_theta_nu(noise_type[0], B_nu, nu, V, prevV, h,
                                          nu_lower_bound, single_V);
        }
      } else {
        grad = NoiseUtil::grad_theta_nu(noise_type[0], B_nu, nu, V, prevV, h,
                                        nu_lower_bound, single_V);
      }
    }
  }

  VectorXd grad_out = -grad;
  if (n_theta_nu > 0 && !fix_flag[latent_fix_theta_nu]) {
    if (prior_nu_target == "coef") {
      for (int l = 0; l < n_theta_nu; ++l) {
        grad_out(l) +=
            PriorUtil::d_log_dens(prior_nu_type, prior_nu_param, theta_nu(l));
      }
    } else {
      VectorXd nu_lp = B_nu * theta_nu;
      grad_out +=
          B_nu.transpose() * prior_score_vec(prior_nu_type, prior_nu_param, nu_lp);
    }
  }

  deriv_cache.grad_theta_nu = grad_out;
  deriv_cache.grad_nu_ready = true;
}

// at init, and after each set parameter. need_precond toggles whether to build
// costly second-derivative information for preconditioning.
void Latent::update_each_iter(bool need_precond) {
  bool need_full_update = !state_ready_;
  bool need_upgrade = need_precond && !state_has_precond_terms_;
  if (!need_full_update && !need_upgrade) {
    return;
  }
  auto t_total_start = std::chrono::steady_clock::now();
  long long t_noise_ms = 0, t_trace_ms = 0;

  UpdateOptions uopts;
  uopts.compute_K = true;
  uopts.compute_Z = true;
  uopts.compute_dK = true; // prefer analytic dK when available
  uopts.compute_dZ = true;
  uopts.compute_d2K = need_precond;
  uopts.compute_d2Z = need_precond;
  uopts.compute_HK_trace = need_precond && !zero_trace;
  // The operator trace is read below only under these same conditions; saying
  // so lets update_all skip factorizing K when nothing will use it.
  uopts.compute_trace = !fix_flag[latent_fix_theta_K] && !zero_trace;
  uopts.robust_reanalyze = robust_;
  uopts.n_trace_iter = n_trace_iter_;
  uopts.solver_type = solver_type_;
  uopts.nonsym_solver = nonsym_solver_;
  uopts.fix_mask_thetaK = ope->get_fix_mask_K();
  // Fresh probes each iteration.
  uopts.trace_seed = static_cast<unsigned int>(latent_rng());
  ope->update_all(theta_K, uopts);

  mu = B_mu * theta_mu;
  sigma = (B_sigma * theta_sigma).array().exp();
  nu = nu_lower_bound + (B_nu * theta_nu).array().exp();

  // update p,a,b, depend on nu, h
  auto t_noise_start = std::chrono::steady_clock::now();
  for (int i = 0; i < n_noise; i++) {
    if (noise_type[i] == "normal")
      continue;

    // update p_vec, a_vec, b_vec
    NoiseUtil::update_gig(noise_type[i], nu, p_vec, a_vec, b_vec, h, single_V);

    // update p_inc, a_inc
    p_inc = VectorXd::Constant(V_size, -0.5 * dim);
    a_inc = mu.cwiseQuotient(sigma).array().pow(2);
  }
  t_noise_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                   std::chrono::steady_clock::now() - t_noise_start)
                   .count();

  // Update traces (pulled from operator)
  if (!fix_flag[latent_fix_theta_K] && !zero_trace) {
    auto t_trace_start = std::chrono::steady_clock::now();
    VectorXd tv = ope->get_trace_trK();
    if (tv.size() == n_theta_K) {
      for (int i = 0; i < n_theta_K; ++i)
        trace[i] = tv(i);
    }
    t_trace_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                     std::chrono::steady_clock::now() - t_trace_start)
                     .count();
  }

  prevV = V;
  prevW = W;
  state_ready_ = true;
  state_has_precond_terms_ = need_precond;
  invalidate_derivatives();
  if (debug) {
    auto t_total_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                          std::chrono::steady_clock::now() - t_total_start)
                          .count();
    ngme_io::out() << "[latent] update_each_iter timing (ms): total=" << t_total_ms
                << ", noise_update=" << t_noise_ms << ", trace=" << t_trace_ms
                << std::endl;
  }
}

// main function for computing analytic Hessian blocks
// The reparameterisation mixes mu, sigma and nu, so it only makes sense when
// all three are free, stationary, and nu is unshifted.
bool Latent::nig_std_active() const {
  if (nig_param_mode_ == 0)
    return false;
  if (n_theta_mu != 1 || n_theta_sigma != 1 || n_theta_nu != 1)
    return false;
  if (noise_type.empty() ||
      (noise_type[0] != "nig" && noise_type[0] != "normal_nig"))
    return false;
  if (nu_lower_bound != 0.0)
    return false;
  if (fix_flag[latent_fix_theta_mu] || fix_flag[latent_fix_theta_nu])
    return false;
  for (bool fixed : fix_theta_sigma_vec)
    if (fixed)
      return false;
  return true;
}

// xi from Cabral sec 2.2. Because zeta* = zeta sqrt(eta) gives zeta^2 eta =
// zeta*^2, xi depends on zeta* alone, and lies in [1/2, 1].
static inline double nig_xi(double zstar) {
  const double z2 = zstar * zstar;
  return 1.0 + z2 - std::fabs(zstar) * std::sqrt(1.0 + z2);
}

// native (mu, log sigma, log nu) -> t
VectorXd Latent::nig_std_from_native() const {
  const double mu_v = theta_mu(0);
  const double sigma_v = std::exp(theta_sigma(0));
  const double nu_v = std::exp(theta_nu(0));
  const double eta = 1.0 / nu_v;
  const double zeta = mu_v / sigma_v;
  VectorXd t(3);
  t(0) = std::log(std::sqrt(sigma_v * sigma_v + mu_v * mu_v * eta));
  if (nig_param_mode_ == 1) {
    t(1) = zeta;
    t(2) = std::log(eta);
  } else {
    const double zstar = zeta * std::sqrt(eta);
    const double xi = nig_xi(zstar);
    t(1) = (nig_param_mode_ == 3) ? std::asinh(zstar) : zstar;
    t(2) = std::log(eta / (xi * xi));
  }
  return t;
}

// t -> (theta_mu, theta_sigma, theta_nu)
VectorXd Latent::nig_native_from_t(const VectorXd &t) const {
  const double sm = std::exp(t(0));
  double zeta, eta;
  if (nig_param_mode_ == 1) {
    zeta = t(1);
    eta = std::exp(t(2));
  } else {
    // mode 3 puts zeta* on an asinh scale. xi -> 1/2 as |zeta*| grows, so the
    // excess kurtosis saturates at ~36 and the likelihood is flat out there;
    // asinh is linear near 0 and logarithmic in the tail, which matches that
    // sensitivity. Note zeta* = sinh(w) is exactly psi = tanh(w) under the
    // identity xi = 1/(1 + |psi|), psi = zeta*/sqrt(1 + zeta*^2).
    const double zstar = (nig_param_mode_ == 3) ? std::sinh(t(1)) : t(1);
    const double xi = nig_xi(zstar);
    eta = std::exp(t(2)) * xi * xi;   // eta = eta* xi^2
    zeta = zstar / std::sqrt(eta);
  }
  const double D2 = 1.0 + zeta * zeta * eta;
  const double sigma_v = sm / std::sqrt(D2);
  VectorXd out(3);
  out(0) = zeta * sigma_v;        // theta_mu = mu
  out(1) = std::log(sigma_v);     // theta_sigma = log sigma
  out(2) = -std::log(eta);        // theta_nu = log nu = -log eta
  return out;
}

void Latent::nig_std_to_native(const VectorXd &t, double &theta_mu_out,
                               double &theta_sigma_out,
                               double &theta_nu_out) const {
  VectorXd v = nig_native_from_t(t);
  theta_mu_out = v(0);
  theta_sigma_out = v(1);
  theta_nu_out = v(2);
}

MatrixXd Latent::nig_std_jacobian() const {
  VectorXd t = nig_std_from_native();
  MatrixXd J(3, 3);
  for (int j = 0; j < 3; ++j) {
    const double st = 1e-6 * std::max(1.0, std::fabs(t(j)));
    VectorXd tp = t, tm = t;
    tp(j) += st;
    tm(j) -= st;
    J.col(j) = (nig_native_from_t(tp) - nig_native_from_t(tm)) / (2.0 * st);
  }
  return J;
}

void Latent::compute_hessian_blocks(bool rao_blackwell) {
  (void)rao_blackwell;
  // Initialize blocks
  if (n_theta_K > 0) {
    hess_cache.H_K = MatrixXd::Zero(n_theta_K, n_theta_K);
  }
  if (n_theta_mu > 0) {
    hess_cache.H_mu = MatrixXd::Zero(n_theta_mu, n_theta_mu);
  }
  if (n_theta_sigma > 0) {
    hess_cache.H_sigma = MatrixXd::Zero(n_theta_sigma, n_theta_sigma);
  }
  if (n_theta_nu > 0) {
    hess_cache.H_nu = MatrixXd::Zero(n_theta_nu, n_theta_nu);
  }
  if (n_theta_K > 0 && n_theta_mu > 0) {
    hess_cache.H_K_mu = MatrixXd::Zero(n_theta_K, n_theta_mu);
  }
  if (n_theta_K > 0 && n_theta_sigma > 0) {
    hess_cache.H_K_sigma = MatrixXd::Zero(n_theta_K, n_theta_sigma);
  }
  if (n_theta_mu > 0 && n_theta_sigma > 0) {
    hess_cache.H_mu_sigma = MatrixXd::Zero(n_theta_mu, n_theta_sigma);
  }

  // Common terms
  VectorXd Dinv =
      (sigma.array().square().matrix().cwiseProduct(V)).cwiseInverse();
  VectorXd m = mu.cwiseProduct(V - h);
  VectorXd r = getK() * W - m;

  // H_K: W|V contribution exact; add operator-side trace terms later; Y|W
  // handled in Block
  if (n_theta_K > 0) {
    std::vector<VectorXd> KjW(n_theta_K);
    for (int j = 0; j < n_theta_K; ++j)
      KjW[j] = get_dK(j) * W;
    for (int j = 0; j < n_theta_K; ++j) {
      for (int k = j; k < n_theta_K; ++k) {
        double term_vec = -KjW[k].cwiseProduct(Dinv).dot(KjW[j]);
        const auto &Kjk = ope->get_d2K(j, k);
        if (Kjk.rows() > 0) {
          VectorXd KjkW = Kjk * W;
          term_vec += -r.cwiseProduct(Dinv).dot(KjkW);
        }
        hess_cache.H_K(j, k) = term_vec;
        if (j != k)
          hess_cache.H_K(k, j) = term_vec;
      }
    }
    // Add operator-side trace terms if available
    const MatrixXd &HKtr = ope->get_HK_trace();
    if (HKtr.rows() == n_theta_K && HKtr.cols() == n_theta_K) {
      hess_cache.H_K += HKtr;
    }
    // Same convention as above: subtract to ridge the information.
    hess_cache.H_K.diagonal().array() -= 1e-9;
  }

  // H_K_mu cross block from W|V: B^T diag(V-h) D K_j W (per j)
  if (n_theta_K > 0 && n_theta_mu > 0) {
    hess_cache.H_K_mu.setZero(n_theta_K, n_theta_mu);
    VectorXd Dinv =
        (sigma.array().square().matrix().cwiseProduct(V)).cwiseInverse();
    VectorXd diagVH =
        (V - h).cwiseProduct(Dinv); // elementwise (V-h)/(sigma^2 V)
    for (int j = 0; j < n_theta_K; ++j) {
      VectorXd KjW = get_dK(j) * W;          // n x 1
      VectorXd z = diagVH.cwiseProduct(KjW); // n x 1
      VectorXd row = B_mu.transpose() * z;   // n_theta_mu x 1
      hess_cache.H_K_mu.row(j) = row.transpose();
    }
  }

  // H_K_sigma cross block from W|V: 2 B_sigma^T ( r ∘ (D K_j W) )
  if (n_theta_K > 0 && n_theta_sigma > 0) {
    hess_cache.H_K_sigma.setZero(n_theta_K, n_theta_sigma);
    for (int j = 0; j < n_theta_K; ++j) {
      VectorXd KjW = get_dK(j) * W;                        // n x 1
      VectorXd z = r.cwiseProduct(Dinv).cwiseProduct(KjW); // n x 1
      VectorXd row = 2.0 * B_sigma.transpose() * z;        // n_theta_sigma x 1
      hess_cache.H_K_sigma.row(j) = row.transpose();
    }
  }

  // H_mu: - B_mu^T diag( (V-h)^2 / (sigma^2 V) ) B_mu
  if (n_theta_mu > 0) {
    VectorXd z = (V - h);
    VectorXd wmu = z.cwiseProduct(z).cwiseProduct(Dinv);
    hess_cache.H_mu = -(B_mu.transpose() * wmu.asDiagonal() * B_mu);
  }

  // H_sigma: - B_sigma^T diag( 2 (r.^2) / (sigma^2 V) ) B_sigma
  if (n_theta_sigma > 0) {
    VectorXd wsig = 2.0 * r.array().square().matrix().cwiseProduct(Dinv);
    hess_cache.H_sigma = -(B_sigma.transpose() * wsig.asDiagonal() * B_sigma);
  }

  // H_mu_sigma cross: -2 B_sigma^T diag( (Kz ⊙ (V-h)) / (σ^2 ⊙ V) ) B_mu
  if (n_theta_mu > 0 && n_theta_sigma > 0) {
    VectorXd wms = 2.0 * r.cwiseProduct(V - h).cwiseProduct(Dinv);
    // The provided form is d^2/d(theta_sigma d theta_mu^T) = -2 B_sigma^T
    // diag(...) B_mu Our cache stores mu x sigma; so take transpose
    MatrixXd H_sigma_mu = -(B_sigma.transpose() * wms.asDiagonal() * B_mu);
    hess_cache.H_mu_sigma = H_sigma_mu.transpose(); // (mu x sigma)
  }
  // H_nu: NIG prior contribution (no cross terms)
  if (n_theta_nu > 0 && !fix_flag[latent_fix_theta_nu]) {
    if (noise_type[0] != "normal") {
      if (n_noise == 1) {
        if (noise_type[0] == "normal_nig") {
          // normal+nig packaged in a single latent: use first half for NIG part
          MatrixXd B_nu_half = B_nu.block(0, 0, B_nu.rows() / 2, B_nu.cols());
          VectorXd nu_half = nu.segment(0, nu.size() / 2);
          VectorXd V_half = V.segment(0, V.size() / 2);
          VectorXd h_half = h.segment(0, h.size() / 2);
          hess_cache.H_nu = NoiseUtil::hess_theta_nu(
              noise_type[0], B_nu_half, nu_half, V_half, h_half, nu_lower_bound,
              single_V);
        } else {
          hess_cache.H_nu =
              NoiseUtil::hess_theta_nu(noise_type[0], B_nu, nu, V, h,
                                       nu_lower_bound, single_V);
        }
      } else {
        // Multi-noise: follow gradient path (use primary noise_type[0])
        hess_cache.H_nu =
            NoiseUtil::hess_theta_nu(noise_type[0], B_nu, nu, V, h,
                                     nu_lower_bound, single_V);
      }
      hess_cache.H_nu = -hess_cache.H_nu;
    }
  }
}

MatrixXd Latent::preconditioner() const {
  if (!precond_cache_valid_) {
    throw std::runtime_error(
        "Latent::preconditioner() called before compute_grad_and_hessian()");
  }
  return last_precond_;
}
