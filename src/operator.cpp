#include "operator.h"

// for initialize Latent models
Operator::~Operator() = default;

// Unified updater: update K, Z, (optionally) dK, dZ, factorization, and traces
void Operator::update_all(const VectorXd& theta, const UpdateOptions& opts) {
    // 1) Build K and Z once at base theta
    build_KZ(theta);

    // Ensure storage ready
    if ((int)dK.size() != n_theta_K) dK.assign(n_theta_K, SparseMatrix<double>(h.size(), h.size()));
    if ((int)dZ.size() != n_theta_K) dZ.assign(n_theta_K, SparseMatrix<double>(h.size(), h.size()));
    if ((int)d2K.size() != n_theta_K) d2K.assign(n_theta_K, std::vector<SparseMatrix<double>>(n_theta_K, SparseMatrix<double>(h.size(), h.size())));
    if ((int)d2Z.size() != n_theta_K) d2Z.assign(n_theta_K, std::vector<SparseMatrix<double,0,int>>(n_theta_K, SparseMatrix<double,0,int>(h.size(), h.size())));

    // 2) Derivatives: analytic or numeric
    bool have_analytic = false;
    if ((opts.compute_dK || opts.compute_dZ) && opts.prefer_analytic_dK) {
        have_analytic = update_dKdZ(theta);
    }
    if ((opts.compute_dK || opts.compute_dZ) && !have_analytic) {
        // Numeric differencing around theta
        // Backup base K/Z
        SparseMatrix<double> K_base = K;
        SparseMatrix<double> Z_base = Z;
        auto is_fixed = [&](int j){ return (!opts.fix_mask_thetaK.empty() && opts.fix_mask_thetaK[j]); };
        for (int j=0; j<n_theta_K; ++j) {
            if (is_fixed(j)) { dK[j].setZero(); dZ[j].setZero(); continue; }
            double eps = opts.eps_dK;
            if (opts.diff_dK_mode == DiffMode::Forward) {
                VectorXd th_f = theta; th_f(j) += eps;
                build_KZ(th_f);
                if (opts.compute_dK) dK[j] = (K - K_base) * (1.0/eps);
                if (opts.compute_dZ) dZ[j] = (Z - Z_base) * (1.0/eps);
            } else {
                VectorXd th_p = theta; th_p(j) += eps;
                build_KZ(th_p); SparseMatrix<double> Kp = K; SparseMatrix<double,0,int> Zp = Z;
                VectorXd th_m = theta; th_m(j) -= eps;
                build_KZ(th_m); SparseMatrix<double> Km = K; SparseMatrix<double,0,int> Zm = Z;
                if (opts.compute_dK) dK[j] = (Kp - Km) * (1.0/(2.0*eps));
                if (opts.compute_dZ) dZ[j] = (Zp - Zm) * (1.0/(2.0*eps));
            }
            // Restore base
            K = K_base; Z = Z_base;
        }
    }
    // 2b) Second derivatives (numeric unless analytic provided)
    bool have_analytic2 = false;
    if ((opts.compute_d2K || opts.compute_d2Z) && opts.prefer_analytic_d2K) {
        have_analytic2 = update_d2Kd2Z(theta);
    }
    if ((opts.compute_d2K || opts.compute_d2Z) && !have_analytic2) {
        // Mixed partials central difference for j!=k, and pure second for j==k
        SparseMatrix<double> K0 = K;
        SparseMatrix<double,0,int> Z0 = Z;
        auto step_for = [&](double th){ double s = std::max(1.0, std::abs(th)); return opts.eps_dK * s; };
        for (int j=0; j<n_theta_K; ++j) {
            double hj = step_for(theta(j));
            // Pure second K_{jj}
            if (opts.compute_d2K) {
                VectorXd thp = theta; thp(j) += hj; build_KZ(thp); SparseMatrix<double> Kp = K;
                VectorXd thm = theta; thm(j) -= hj; build_KZ(thm); SparseMatrix<double> Km = K;
                d2K[j][j] = (Kp - 2.0*K0 + Km) * (1.0/(hj*hj));
            }
            if (opts.compute_d2Z) {
                VectorXd thp = theta; thp(j) += hj; build_KZ(thp); auto Zp = Z;
                VectorXd thm = theta; thm(j) -= hj; build_KZ(thm); auto Zm = Z;
                d2Z[j][j] = (Zp - 2.0*Z0 + Zm) * (1.0/(hj*hj));
            }
            // Mixed K_{jk}, j<k
            for (int k=j+1; k<n_theta_K; ++k) {
                double hk = step_for(theta(k));
                VectorXd th_pp = theta; th_pp(j) += hj; th_pp(k) += hk; build_KZ(th_pp); auto Kpp = K; auto Zpp = Z;
                VectorXd th_pm = theta; th_pm(j) += hj; th_pm(k) -= hk; build_KZ(th_pm); auto Kpm = K; auto Zpm = Z;
                VectorXd th_mp = theta; th_mp(j) -= hj; th_mp(k) += hk; build_KZ(th_mp); auto Kmp = K; auto Zmp = Z;
                VectorXd th_mm = theta; th_mm(j) -= hj; th_mm(k) -= hk; build_KZ(th_mm); auto Kmm = K; auto Zmm = Z;
                if (opts.compute_d2K) {
                    d2K[j][k] = (Kpp - Kpm - Kmp + Kmm) * (1.0/(4.0*hj*hk));
                    d2K[k][j] = d2K[j][k];
                }
                if (opts.compute_d2Z) {
                    auto Zmix = (Zpp - Zpm - Zmp + Zmm) * (1.0/(4.0*hj*hk));
                    d2Z[j][k] = Zmix; d2Z[k][j] = Zmix;
                }
            }
            // Restore base for next j
            K = K0; Z = Z0;
        }
    }

    // 3) Factorization (init once with symmetric flag and control params)
    if (!llt_inited) {
        // If opts.solver_type not set by caller, stays default 0
        cholK_solver.init(K.rows(), opts.n_trace_iter, symmetric, opts.solver_type);
        llt_inited = true;
    }
    if (!analyzed_cholK) { cholK_solver.analyze(K); analyzed_cholK = true; }
    cholK_solver.compute(K);

    // 4) Traces (only if dK available)
    trace_ready = false;
    if (!dK.empty() && n_theta_K > 0) {
        if (trace_vals.size() != n_theta_K) trace_vals = VectorXd::Zero(n_theta_K);
        // tr(K^{-1} dK)
        for (int j=0; j<n_theta_K; ++j) {
            if (!opts.fix_mask_thetaK.empty() && opts.fix_mask_thetaK[j]) { trace_vals(j) = 0.0; continue; }
            trace_vals(j) = cholK_solver.trace(dK[j]);
        }
    }
    // global trace mode used; no per-call state needed
    trace_ready = true;

    // 5) H_K trace block: -tr(K^{-1}K_k K^{-1}K_j) + tr(K^{-1}K_{jk})
    HK_trace_ready = false;
    if (opts.compute_HK_trace && n_theta_K > 0) {
        if (HK_trace.rows() != n_theta_K || HK_trace.cols() != n_theta_K)
            HK_trace = MatrixXd::Zero(n_theta_K, n_theta_K);
        else
            HK_trace.setZero();
        for (int j=0; j<n_theta_K; ++j) {
            for (int k=j; k<n_theta_K; ++k) {
                double t1 = 0.0, t2 = 0.0;
                // -tr(K^{-1} K_k K^{-1} K_j)
                if (dK[k].rows()>0 && dK[j].rows()>0) {
                    t1 = - cholK_solver.trace2(dK[k], dK[j]);
                }
                // + tr(K^{-1} K_{jk})
                if (!d2K.empty() && d2K[j][k].rows()>0) {
                    t2 = cholK_solver.trace(d2K[j][k]);
                }
                double val = t1 + t2;
                HK_trace(j,k) = val;
                if (j!=k) HK_trace(k,j) = val;
            }
        }
        HK_trace_ready = true;
    }
}

std::shared_ptr<Operator> OperatorFactory::create(const Rcpp::List& operator_in) {
  string model_type = Rcpp::as<string> (operator_in["model"]);
  string generic_type = Rcpp::as<string> (operator_in["generic_type"]);

  if (model_type == "generic") {
    return std::make_shared<Generic>(operator_in);
  } else if (model_type == "generic_ns") {
    return std::make_shared<generic_ns>(operator_in);
  } else if (generic_type == "generic") {
    return std::make_shared<Generic>(operator_in);
  } else if (generic_type == "generic_ns") {
    return std::make_shared<generic_ns>(operator_in);
  } else if (model_type == "tp") {
    return std::make_shared<Tensor_prod>(operator_in);
  } else if (model_type == "spacetime") {
    return std::make_shared<Spacetime>(operator_in);
  } else if (model_type == "matern") {
    // Unify: always use the new Matern operator.
    // It handles stationary/non-stationary, integer/fractional alpha, and produces Z.
    return std::make_shared<Matern>(operator_in);
  } else if (model_type == "arma") {
    return std::make_shared<Arma>(operator_in);
  } else if (model_type == "re") {
    return std::make_shared<Randeff>(operator_in);
  } else if (model_type == "bv") {
    return std::shared_ptr<Operator>(new Bivar(operator_in));
  } else if (model_type == "bv_normal") {
    return std::shared_ptr<Operator>(new Bivar_normal_ope(operator_in));
  } else if (model_type == "bv_matern_normal") {
    return std::shared_ptr<Operator>(new bv_matern_normal(operator_in));
  } else if (model_type == "bv_matern_nig") {
    return std::shared_ptr<Operator>(new bv_matern_nig(operator_in));
  } else {
    throw std::runtime_error("Unknown model.");
  }
};
