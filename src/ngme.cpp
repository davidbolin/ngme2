// implement the Ngme class and rand effect class
#include "ngme.h"

#include "include/thread_io.h"
#include <atomic>

#ifdef _OPENMP
#include <omp.h>
#pragma omp declare reduction(vec_plus : Eigen::VectorXd : omp_out += omp_in)  \
    initializer(omp_priv = Eigen::VectorXd::Zero(omp_orig.size()))
#pragma omp declare reduction(mat_plus : Eigen::MatrixXd : omp_out += omp_in)  \
    initializer(omp_priv =                                                     \
                    Eigen::MatrixXd::Zero(omp_orig.rows(), omp_orig.cols()))
#endif

// --------------- Ngme class ----------------
Ngme::Ngme(const Rcpp::List &R_ngme, unsigned long seed, int sampling_strategy,
           int num_threads_repl, double start_sd)
    : n_params(Rcpp::as<int>(R_ngme["n_params"])),
      sampling_strategy(sampling_strategy), num_threads_repl(num_threads_repl),
      gen(seed), debug(false) {
  Rcpp::List control_ngme = R_ngme["control_ngme"];
  debug = Rcpp::as<bool>(control_ngme["debug"]);
  stepsize =
      VectorXd::Constant(n_params, Rcpp::as<double>(control_ngme["stepsize"]));

  Rcpp::List list_replicates = Rcpp::as<Rcpp::List>(R_ngme["replicates"]);
  n_repl = list_replicates.size();

  // build for each replicates
  for (int i = 0; i < n_repl; i++) {
    Rcpp::List ngme_repl = Rcpp::as<Rcpp::List>(list_replicates[i]);
    ngme_repls.push_back(std::make_shared<BlockModel>(ngme_repl, seed + i));
    num_each_repl.push_back(ngme_repls[i]->get_n_obs());
  }
  sum_num_each_repl =
      std::accumulate(num_each_repl.begin(), num_each_repl.end(), 0.0);

  // Init the random number generator
  weighted_sampler = std::discrete_distribution<int>(num_each_repl.begin(),
                                                     num_each_repl.end());

  VectorXd p = ngme_repls[0]->get_parameter();
  std::normal_distribution<double> distribution(0.0, start_sd);
  for (int i = 0; i < p.size(); ++i) {
    p[i] += distribution(gen);
  }
  ngme_repls[0]->set_parameter_and_update(p, true);
  current_param_ = p;
  repl_dirty_.assign(n_repl, static_cast<unsigned char>(0));
  repl_precond_ready_.assign(n_repl, static_cast<unsigned char>(0));
  if (sampling_strategy == Strategy::ws && n_repl > 0) {
    std::fill(repl_dirty_.begin(), repl_dirty_.end(),
              static_cast<unsigned char>(1));
    repl_dirty_[0] = static_cast<unsigned char>(0);
    repl_precond_ready_[0] = static_cast<unsigned char>(1);
  }
}

void Ngme::compute(bool with_precond, double eps) {
  // Aggregate over replicates
  last_grad_ = VectorXd::Zero(n_params);
  // If with_precond=true, also compute and cache preconditioner for current
  // strategy
  if (with_precond) {
    last_precond_ = MatrixXd::Zero(n_params, n_params);
    precond_valid_ = false;
  }
  if (sampling_strategy == Strategy::all) {
    VectorXd grad_acc = VectorXd::Zero(n_params);
    MatrixXd precond_acc = with_precond
                               ? MatrixXd::Zero(n_params, n_params)
                               : MatrixXd::Zero(0, 0);
    std::atomic<bool> compute_failed(false);
    std::string compute_error;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(num_threads_repl)        \
    reduction(vec_plus : grad_acc) reduction(mat_plus : precond_acc)
#endif
    for (int i = 0; i < n_repl; i++) {
      if (compute_failed.load(std::memory_order_relaxed))
        continue;
      try {
        ngme_repls[i]->compute_grad_and_hessian(with_precond, eps);
        grad_acc += ngme_repls[i]->get_gradient();
        if (with_precond)
          precond_acc += ngme_repls[i]->get_preconditioner();
      } catch (const std::exception &e) {
#pragma omp critical(ngme_parallel_exception)
        {
          if (!compute_failed.load(std::memory_order_relaxed))
            compute_error = e.what();
          compute_failed.store(true, std::memory_order_relaxed);
        }
      }
    }
    if (compute_failed.load(std::memory_order_relaxed))
      throw std::runtime_error(compute_error);
    last_grad_ = grad_acc;
    if (with_precond)
      last_precond_ = precond_acc;
  } else { // ws
    int idx = weighted_sampler(gen);
    sync_repl_if_needed(idx, with_precond);
    ngme_repls[idx]->compute_grad_and_hessian(with_precond, eps);
    // Replicate idx is drawn with probability p_idx = n_idx / N, so the raw
    // single-replicate gradient estimates sum_i p_i g_i, not sum_i g_i as the
    // "all" branch does. Without the 1/p_idx importance weight the two
    // strategies differ in scale by roughly n_repl, which silently rescales the
    // effective step size when the strategy is switched.
    double w = 1.0;
    if (num_each_repl[idx] > 0 && sum_num_each_repl > 0)
      w = sum_num_each_repl / num_each_repl[idx];
    last_grad_ = w * ngme_repls[idx]->get_gradient();
    if (with_precond)
      last_precond_ = w * ngme_repls[idx]->get_preconditioner();
  }
  grad_valid_ = true;
  if (with_precond) {
    precond_valid_ = true;
  }
  mark_computed();
}


MatrixXd Ngme::precond() {
  if (!grad_valid_)
    throw std::runtime_error("precond() called before compute()");
  return -last_precond_;
}

VectorXd Ngme::grad() {
  if (!grad_valid_)
    throw std::runtime_error("grad() called before compute()");
  return -last_grad_;
}

void Ngme::burn_in(int iterations) {
  sync_all_repls_if_needed(false);
  std::atomic<bool> burnin_failed(false);
  std::string burnin_error;
#pragma omp parallel for schedule(static) num_threads(num_threads_repl)
  for (int i = 0; i < n_repl; i++) {
    if (burnin_failed.load(std::memory_order_relaxed)) {
      continue;
    }
    try {
      ngme_repls[i]->burn_in(iterations);
    } catch (const std::exception &e) {
#pragma omp critical(ngme_parallel_exception)
      {
        if (!burnin_failed.load(std::memory_order_relaxed)) {
          burnin_failed.store(true, std::memory_order_relaxed);
          burnin_error = e.what();
        }
      }
    } catch (...) {
#pragma omp critical(ngme_parallel_exception)
      {
        if (!burnin_failed.load(std::memory_order_relaxed)) {
          burnin_failed.store(true, std::memory_order_relaxed);
          burnin_error = "unknown exception";
        }
      }
    }
  }
  if (burnin_failed.load(std::memory_order_relaxed)) {
    Rcpp::stop("C++ exception in replicate burn-in: %s", burnin_error.c_str());
  }
}

VectorXd Ngme::get_parameter() {
  VectorXd p = ngme_repls[0]->get_parameter();
  return p;
}

void Ngme::set_parameter_and_update(const VectorXd &p, bool with_precond) {
  current_param_ = p;
  if (sampling_strategy == Strategy::ws) {
    std::fill(repl_dirty_.begin(), repl_dirty_.end(),
              static_cast<unsigned char>(1));
    std::fill(repl_precond_ready_.begin(), repl_precond_ready_.end(),
              static_cast<unsigned char>(0));
  } else {
    // set the same parameter for all latent
    // std::chrono::steady_clock::time_point begin =
    // std::chrono::steady_clock::now();
    std::atomic<bool> update_failed(false);
    std::string update_error;
#pragma omp parallel for schedule(static) num_threads(num_threads_repl)
    for (int i = 0; i < n_repl; i++) {
      if (update_failed.load(std::memory_order_relaxed)) {
        continue;
      }
      try {
        ngme_repls[i]->set_parameter_and_update(p, with_precond);
        repl_dirty_[i] = static_cast<unsigned char>(0);
        repl_precond_ready_[i] =
            static_cast<unsigned char>(with_precond ? 1 : 0);
      } catch (const std::exception &e) {
#pragma omp critical(ngme_parallel_exception)
        {
          if (!update_failed.load(std::memory_order_relaxed)) {
            update_failed.store(true, std::memory_order_relaxed);
            update_error = e.what();
          }
        }
      } catch (...) {
#pragma omp critical(ngme_parallel_exception)
        {
          if (!update_failed.load(std::memory_order_relaxed)) {
            update_failed.store(true, std::memory_order_relaxed);
            update_error = "unknown exception";
          }
        }
      }
    }
    if (update_failed.load(std::memory_order_relaxed)) {
      Rcpp::stop("C++ exception while setting parameters: %s",
                 update_error.c_str());
    }
  }

  // set the different parameter for each random effect
  if (debug)
    ngme_io::out() << "set_parameter() in ngme class" << std::endl;
  // Invalidate caches
  grad_valid_ = false;
  precond_valid_ = false;
  reset_computed();
}

void Ngme::sync_repl_if_needed(int idx, bool with_precond) const {
  if (sampling_strategy != Strategy::ws || idx < 0 || idx >= n_repl) {
    return;
  }
  const bool need_sync = (repl_dirty_[idx] != 0) ||
                         (with_precond && repl_precond_ready_[idx] == 0);
  if (!need_sync) {
    return;
  }
  ngme_repls[idx]->set_parameter_and_update(current_param_, with_precond);
  repl_dirty_[idx] = static_cast<unsigned char>(0);
  repl_precond_ready_[idx] = static_cast<unsigned char>(with_precond ? 1 : 0);
}

void Ngme::sync_all_repls_if_needed(bool with_precond) const {
  if (sampling_strategy != Strategy::ws) {
    return;
  }
  std::atomic<bool> sync_failed(false);
  std::string sync_error;
#pragma omp parallel for schedule(static) num_threads(num_threads_repl)
  for (int i = 0; i < n_repl; i++) {
    if (sync_failed.load(std::memory_order_relaxed)) {
      continue;
    }
    const bool need_sync = (repl_dirty_[i] != 0) ||
                           (with_precond && repl_precond_ready_[i] == 0);
    if (!need_sync) {
      continue;
    }
    try {
      ngme_repls[i]->set_parameter_and_update(current_param_, with_precond);
      repl_dirty_[i] = static_cast<unsigned char>(0);
      repl_precond_ready_[i] =
          static_cast<unsigned char>(with_precond ? 1 : 0);
    } catch (const std::exception &e) {
#pragma omp critical(ngme_parallel_exception)
      {
        if (!sync_failed.load(std::memory_order_relaxed)) {
          sync_failed.store(true, std::memory_order_relaxed);
          sync_error = e.what();
        }
      }
    } catch (...) {
#pragma omp critical(ngme_parallel_exception)
      {
        if (!sync_failed.load(std::memory_order_relaxed)) {
          sync_failed.store(true, std::memory_order_relaxed);
          sync_error = "unknown exception";
        }
      }
    }
  }
  if (sync_failed.load(std::memory_order_relaxed)) {
    Rcpp::stop("C++ exception while syncing replicates: %s",
               sync_error.c_str());
  }
}
