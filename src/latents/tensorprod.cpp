// getK from matern class
// getK from ar class
// newK = K1 (time) %x% k2 (space)
// dKtx_t = dKt_t %x% Kx
// dKtx_x = kt %x% dKx_x

#include "../operator.h"
#include "thread_io.h"
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include "MatrixAlgebra.h"
#include <unsupported/Eigen/KroneckerProduct>

// tensor product for the C G class
Tensor_prod::Tensor_prod(const Rcpp::List &operator_list)
    : Operator(operator_list),
      first(OperatorFactory::create(operator_list["first"])),
      second(OperatorFactory::create(operator_list["second"])),
      n_theta_1(first->get_n_theta_K()), n_theta_2(second->get_n_theta_K()) {}

void Tensor_prod::build_KZ(const VectorXd &theta_K) {
  // build_KZ runs at the top of every update_all, so this is where the factors
  // stop being current for the theta the rest of the call will use.
  factors_current_ = false;
  // report the time for this function
  double time = 0;
  auto timer_computeg = std::chrono::steady_clock::now();
  // std::cout << "update K now" << std::endl;
  first->build_KZ(theta_K.segment(0, n_theta_1));
  // std::cout << "update K now1" << std::endl;
  second->build_KZ(theta_K.segment(n_theta_1, n_theta_2));
  // std::cout << "new K size = " << first->getK().rows() *
  // second->getK().rows() << " " << first->getK().cols() *
  // second->getK().cols() << std::endl;

  // use Eigen kronecker product
  KroneckerProductSparse<SparseMatrix<double>, SparseMatrix<double>>
      kroneckerEigen(first->getK(), second->getK());

  kroneckerEigen.evalTo(K);
  // std::cout << "update K now3" << std::endl;

  time = std::chrono::duration_cast<std::chrono::milliseconds>(
             std::chrono::steady_clock::now() - timer_computeg)
             .count();
  // std::cout << "size and time for kronecker product is " << K.rows() << " "
  // << K.cols() << " " << time << std::endl; Z remains identity (set once in
  // Operator base)
}


// Analytic dK for a tensor product.
//
// K = K_1 (x) K_2 is linear in each factor, so
//     dK/dtheta_1j = (dK_1/dtheta_1j) (x) K_2
//     dK/dtheta_2j = K_1 (x) (dK_2/dtheta_2j),
// which is the same algebra compute_traces_structured() already relies on for
// the traces. Falling back to the base class instead meant rebuilding the full
// Kronecker product once per parameter per iteration by finite differences and
// carrying that difference's O(eps) error into the gradient.
//
// Returns false if either factor lacks analytic derivatives, leaving the base
// class to difference as before.
// Bring both factors up to date at theta, with the options the parent's
// derivatives and structured traces between them need. Runs at most once per
// update_all: update_dKdZ gets here first, and update_d2Kd2Z and
// compute_traces_structured then reuse what it left behind.
bool Tensor_prod::update_factors(const VectorXd &theta,
                                 const UpdateOptions &opts) {
  if (factors_current_)
    return true;
  if (n_theta_1 + n_theta_2 != n_theta_K)
    return false;
  const int n1 = static_cast<int>(first->getK().rows());
  const int n2 = static_cast<int>(second->getK().rows());
  if (n1 <= 0 || n2 <= 0)
    return false;

  // Probe budget per factor: a trace of the product needs n_other times fewer
  // probes on the factor for the same accuracy, capped at the factor's size.
  auto budget = [&](int n_self, int n_other) {
    const long long want =
        static_cast<long long>(std::max(1, opts.n_trace_iter)) * n_other;
    return static_cast<int>(std::min<long long>(want, n_self));
  };

  UpdateOptions sub = opts;
  sub.compute_K = true;
  sub.compute_Z = true;
  sub.compute_dK = true;
  sub.compute_dZ = false;
  // A superset of what either consumer needs, so one update serves both.
  sub.compute_d2K = opts.compute_HK_trace || opts.compute_d2K;
  sub.compute_d2Z = false;
  // The structured traces need each factor's own tr(F^-1 dF) whenever EITHER
  // trace family is wanted -- theta_K can be fixed (no first-order trace) while
  // a preconditioner still asks for H_K.
  sub.compute_trace = opts.compute_trace || opts.compute_HK_trace;

  const bool have_mask =
      static_cast<int>(opts.fix_mask_thetaK.size()) == n_theta_K;

  UpdateOptions o1 = sub;
  o1.n_trace_iter = budget(n1, n2);
  if (have_mask)
    o1.fix_mask_thetaK.assign(opts.fix_mask_thetaK.begin(),
                              opts.fix_mask_thetaK.begin() + n_theta_1);
  else
    o1.fix_mask_thetaK.clear();

  UpdateOptions o2 = sub;
  o2.n_trace_iter = budget(n2, n1);
  if (have_mask)
    o2.fix_mask_thetaK.assign(opts.fix_mask_thetaK.begin() + n_theta_1,
                              opts.fix_mask_thetaK.end());
  else
    o2.fix_mask_thetaK.clear();

  // This also restores each factor's K to base theta, which matters because a
  // factor that differences numerically leaves itself at the last perturbed
  // value.
  first->update_all(theta.segment(0, n_theta_1), o1);
  second->update_all(theta.segment(n_theta_1, n_theta_2), o2);
  factors_current_ = true;
  return true;
}

bool Tensor_prod::update_dKdZ(const VectorXd &theta_K) {
  // NGME_TP_NUMERIC_DK=1 falls back to the base class's full-size differencing,
  // so the two can be compared directly.
  static const bool disabled = [] {
    const char *e = std::getenv("NGME_TP_NUMERIC_DK");
    return e && *e && std::string(e) != "0";
  }();
  if (disabled)
    return false;
  if (n_theta_1 + n_theta_2 != n_theta_K)
    return false;
  const VectorXd th1 = theta_K.segment(0, n_theta_1);
  const VectorXd th2 = theta_K.segment(n_theta_1, n_theta_2);

  // build_KZ() has already left both factors at this theta.
  const SparseMatrix<double> K1 = first->getK();
  const SparseMatrix<double> K2 = second->getK();
  if ((int)dK.size() != n_theta_K)
    dK.assign(n_theta_K, SparseMatrix<double>(h.size(), h.size()));

  const double eps = last_eps_dK_;
  auto kron_into = [](const SparseMatrix<double> &A,
                      const SparseMatrix<double> &B,
                      SparseMatrix<double> &out) {
    KroneckerProductSparse<SparseMatrix<double>, SparseMatrix<double>> kp(A, B);
    kp.evalTo(out);
  };

  // Perturbing a parameter of one factor leaves the other factor untouched, so
  //     dK/dtheta_1j = (dK_1/dtheta_1j) (x) K_2
  //     dK/dtheta_2j = K_1 (x) (dK_2/dtheta_2j)
  // and the derivative is needed only on the SMALL factor, then lifted once.
  // The base class instead rebuilds the whole Kronecker product per parameter
  // and subtracts two full-size sparse matrices.
  //
  // Preferred route: ask each factor for its OWN dK. A factor with a closed
  // form (matern, generic) then hands back an exact derivative, and the
  // product inherits it. Falls back to differencing the factor when no options
  // are in scope to drive the update.
  bool from_factors = false;
  if (cur_opts_ != nullptr && update_factors(theta_K, *cur_opts_)) {
    from_factors = true;
    for (int j = 0; j < n_theta_1; ++j)
      kron_into(first->get_dK(j), K2, dK[j]);
    for (int j = 0; j < n_theta_2; ++j)
      kron_into(K1, second->get_dK(j), dK[n_theta_1 + j]);
  }

  if (!from_factors) {
    for (int j = 0; j < n_theta_1; ++j) {
      VectorXd t = th1;
      t(j) += eps;
      first->build_KZ(t);
      SparseMatrix<double> d1 = (first->getK() - K1) * (1.0 / eps);
      kron_into(d1, K2, dK[j]);
    }
    if (n_theta_1 > 0)
      first->build_KZ(th1); // restore

    for (int j = 0; j < n_theta_2; ++j) {
      VectorXd t = th2;
      t(j) += eps;
      second->build_KZ(t);
      SparseMatrix<double> d2 = (second->getK() - K2) * (1.0 / eps);
      kron_into(K1, d2, dK[n_theta_1 + j]);
    }
    if (n_theta_2 > 0)
      second->build_KZ(th2); // restore
  }

  // Z is the identity for a tensor product and carries no parameters, so every
  // dZ is zero. Setting them here also spares the base class differencing an
  // n x n identity once per parameter.
  if ((int)dZ.size() != n_theta_K)
    dZ.assign(n_theta_K, SparseMatrix<double, 0, int>(h.size(), h.size()));
  for (int j = 0; j < n_theta_K; ++j)
    dZ[j].setZero();
  return true;
}

// Second derivatives, from the same Kronecker algebra as the first.
//
//     d2K/dtheta_1j dtheta_1k = (d2K_1/dtheta_1j dtheta_1k) (x) K_2
//     d2K/dtheta_1j dtheta_2k = (dK_1/dtheta_1j) (x) (dK_2/dtheta_2k)
//     d2K/dtheta_2j dtheta_2k = K_1 (x) (d2K_2/dtheta_2j dtheta_2k)
//
// each a single Kronecker lift of something the factors already hold, against
// the base class's four full-size rebuilds per pair.
bool Tensor_prod::update_d2Kd2Z(const VectorXd &theta_K) {
  static const bool disabled = [] {
    const char *e = std::getenv("NGME_TP_NUMERIC_DK");
    return e && *e && std::string(e) != "0";
  }();
  if (disabled)
    return false;
  if (n_theta_1 + n_theta_2 != n_theta_K)
    return false;
  if (cur_opts_ == nullptr || !update_factors(theta_K, *cur_opts_))
    return false;

  const SparseMatrix<double> &K1 = first->getK();
  const SparseMatrix<double> &K2 = second->getK();
  const long long n = (long long)K1.rows() * K2.rows();

  // A factor that could not supply second derivatives hands back an empty
  // matrix; without them the same-factor blocks are unavailable, so give the
  // whole job back to the base class rather than fill in zeros.
  for (int j = 0; j < n_theta_1; ++j)
    if (first->get_d2K(j, j).rows() != K1.rows())
      return false;
  for (int j = 0; j < n_theta_2; ++j)
    if (second->get_d2K(j, j).rows() != K2.rows())
      return false;

  if ((int)d2K.size() != n_theta_K)
    d2K.assign(n_theta_K, std::vector<SparseMatrix<double>>(
                              n_theta_K, SparseMatrix<double>(n, n)));
  if ((int)d2Z.size() != n_theta_K)
    d2Z.assign(n_theta_K, std::vector<SparseMatrix<double, 0, int>>(
                              n_theta_K, SparseMatrix<double, 0, int>(n, n)));

  auto kron_into = [](const SparseMatrix<double> &A,
                      const SparseMatrix<double> &B,
                      SparseMatrix<double> &out) {
    KroneckerProductSparse<SparseMatrix<double>, SparseMatrix<double>> kp(A, B);
    kp.evalTo(out);
  };
  auto put = [&](int a, int b, const SparseMatrix<double> &A,
                 const SparseMatrix<double> &B) {
    kron_into(A, B, d2K[a][b]);
    d2Z[a][b].setZero();
    if (a != b) {
      d2K[b][a] = d2K[a][b];
      d2Z[b][a].setZero();
    }
  };

  for (int j = 0; j < n_theta_1; ++j)
    for (int k = j; k < n_theta_1; ++k)
      put(j, k, first->get_d2K(j, k), K2);

  for (int j = 0; j < n_theta_1; ++j)
    for (int k = 0; k < n_theta_2; ++k)
      put(j, n_theta_1 + k, first->get_dK(j), second->get_dK(k));

  for (int j = 0; j < n_theta_2; ++j)
    for (int k = j; k < n_theta_2; ++k)
      put(n_theta_1 + j, n_theta_1 + k, K1, second->get_d2K(j, k));

  return true;
}

// Traces of a tensor product operator, without ever factorizing K.
//
// With K = K_1 (x) K_2 the derivatives keep the Kronecker structure,
//     dK/dtheta_1j = (dK_1/dtheta_1j) (x) K_2,
//     dK/dtheta_2j = K_1 (x) (dK_2/dtheta_2j),
// and (A (x) B)^-1 = A^-1 (x) B^-1, (A (x) B)(C (x) D) = AC (x) BD give
//     tr(K^-1 (dK_1 (x) K_2)) = tr(K_1^-1 dK_1) tr(I_n2) = n_2 tr(K_1^-1 dK_1)
//     tr(K^-1 (K_1 (x) dK_2)) = n_1 tr(K_2^-1 dK_2).
// The same algebra closes the H_K block. For j, k in the same factor,
//     tr(K^-1 dK_k K^-1 dK_j) = n_other tr(F^-1 dF_k F^-1 dF_j)
//     tr(K^-1 d2K_jk)         = n_other tr(F^-1 d2F_jk),
// so that block is just n_other times the factor's own H_K block. For j in one
// factor and k in the other the two terms are equal and opposite,
//     tr(K^-1 dK_k K^-1 dK_j) = tr(K_1^-1 dK_1j) tr(K_2^-1 dK_2k)
//     d2K/dtheta_1j dtheta_2k = dK_1j (x) dK_2k, same trace,
// so the cross block vanishes identically.
bool Tensor_prod::compute_traces_structured(const VectorXd &theta,
                                            const UpdateOptions &opts) {
  const int n1 = static_cast<int>(first->getK().rows());
  const int n2 = static_cast<int>(second->getK().rows());
  if (n1 <= 0 || n2 <= 0 || n_theta_1 + n_theta_2 != n_theta_K)
    return false;

  // The factors were brought up to date by update_dKdZ earlier in this same
  // update_all; update_factors() is a no-op then, and does the work only when
  // the derivative hooks did not run (no options in scope, or derivatives not
  // requested).
  if (!update_factors(theta, opts))
    return false;

  if (!first->traces_ready() || !second->traces_ready())
    return false;
  const VectorXd &t1 = first->get_trace_trK();
  const VectorXd &t2 = second->get_trace_trK();
  if (t1.size() != n_theta_1 || t2.size() != n_theta_2)
    return false;

  if (trace_vals.size() != n_theta_K)
    trace_vals = VectorXd::Zero(n_theta_K);
  for (int j = 0; j < n_theta_1; ++j)
    trace_vals(j) = n2 * t1(j);
  for (int j = 0; j < n_theta_2; ++j)
    trace_vals(n_theta_1 + j) = n1 * t2(j);

  if (opts.compute_HK_trace) {
    const MatrixXd &H1 = first->get_HK_trace();
    const MatrixXd &H2 = second->get_HK_trace();
    if (H1.rows() != n_theta_1 || H2.rows() != n_theta_2)
      return false;
    if (HK_trace.rows() != n_theta_K || HK_trace.cols() != n_theta_K)
      HK_trace = MatrixXd::Zero(n_theta_K, n_theta_K);
    else
      HK_trace.setZero(); // cross blocks stay zero
    HK_trace.topLeftCorner(n_theta_1, n_theta_1) = n2 * H1;
    HK_trace.bottomRightCorner(n_theta_2, n_theta_2) = n1 * H2;
  }
  return true;
}

// Non-separable Space-time model

Spacetime::Spacetime(const Rcpp::List &operator_list)
    : Operator(operator_list),
      Ct_diag(Rcpp::as<VectorXd>(operator_list["Ct_diag"])),
      Cs_diag(Rcpp::as<VectorXd>(operator_list["Cs_diag"])),
      BtCs(Rcpp::as<SparseMatrix<double, 0, int>>(operator_list["BtCs"])),
      Gs(Rcpp::as<SparseMatrix<double, 0, int>>(operator_list["Gs"])),
      Ct(Rcpp::as<SparseMatrix<double, 0, int>>(operator_list["Ct"])),
      Cs(Rcpp::as<SparseMatrix<double, 0, int>>(operator_list["Cs"])),
      Bx(Rcpp::as<SparseMatrix<double, 0, int>>(operator_list["Bx"])),
      By(Rcpp::as<SparseMatrix<double, 0, int>>(operator_list["By"])),
      S(Rcpp::as<SparseMatrix<double, 0, int>>(operator_list["S"])),
      Bs(Rcpp::as<SparseMatrix<double, 0, int>>(operator_list["Bs"])),
      Hxx(Rcpp::as<SparseMatrix<double, 0, int>>(operator_list["Hxx"])),
      Hyy(Rcpp::as<SparseMatrix<double, 0, int>>(operator_list["Hyy"])),
      Hxy(Rcpp::as<SparseMatrix<double, 0, int>>(operator_list["Hxy"])),
      Hyx(Rcpp::as<SparseMatrix<double, 0, int>>(operator_list["Hyx"])),
      // B_gamma_x (Rcpp::as<MatrixXd> (operator_list["B_gamma_x"])),
      // B_gamma_y (Rcpp::as<MatrixXd> (operator_list["B_gamma_y"])),
      B_gamma_x_list_input(
          Rcpp::as<Rcpp::List>(operator_list["B_gamma_x_list"])),
      B_gamma_y_list_input(
          Rcpp::as<Rcpp::List>(operator_list["B_gamma_y_list"])),
      theta_gamma_x(Rcpp::as<VectorXd>(operator_list["theta_gamma_x"])),
      theta_gamma_y(Rcpp::as<VectorXd>(operator_list["theta_gamma_y"])),
      n_theta_gamma_x(Rcpp::as<int>(operator_list["n_theta_gamma_x"])),
      n_theta_gamma_y(Rcpp::as<int>(operator_list["n_theta_gamma_y"])),
      lambda(Rcpp::as<double>(operator_list["lambda"])),
      alpha(Rcpp::as<double>(operator_list["alpha"])),
      method(Rcpp::as<string>(operator_list["method"])),
      stabilization(Rcpp::as<bool>(operator_list["stabilization"])),
      fix_gamma(Rcpp::as<bool>(operator_list["fix_gamma"])),
      shared_theta_gamma(Rcpp::as<bool>(operator_list["shared_theta_gamma"])),
      nt(Rcpp::as<int>(operator_list["nt"])), B_gamma_x_list(nt - 1),
      B_gamma_y_list(nt - 1) {
  cc_variance_free = operator_list.containsElementNamed("cc_variance_free")
                        ? Rcpp::as<bool>(operator_list["cc_variance_free"])
                        : false;
  stationary_init = operator_list.containsElementNamed("stationary_init")
                        ? Rcpp::as<bool>(operator_list["stationary_init"])
                        : true;
  ns_ = (int)Cs.rows();
  // BtCs without its first block row. The rw1 operator's first row is the
  // trapezoid (0.5, 1, ..., 1, 0.5), an aggregate over EVERY time slice, and it
  // is the only reason K is not block lower-triangular. It also leaves slice 1
  // with no operator block of its own, so that slice is pinned only through the
  // aggregate. Dropping the row and giving slice 1 its stationary block fixes
  // both, and makes tr(K^-1 dK) exact from the diagonal blocks alone.
  {
    std::vector<Eigen::Triplet<double>> trips;
    trips.reserve(BtCs.nonZeros());
    for (int k = 0; k < BtCs.outerSize(); ++k)
      for (SparseMatrix<double, 0, int>::InnerIterator it(BtCs, k); it; ++it)
        if (it.row() >= ns_)
          trips.emplace_back((int)it.row(), (int)it.col(), it.value());
    BtCs_st.resize(BtCs.rows(), BtCs.cols());
    BtCs_st.setFromTriplets(trips.begin(), trips.end());
  }
  // turn B_gamma_x_list_input into B_gamma_x_list
  if (!fix_gamma) {
    for (int i = 0; i < nt - 1; i++) {
      B_gamma_x_list[i] = Rcpp::as<MatrixXd>(B_gamma_x_list_input[i]);
      B_gamma_y_list[i] = Rcpp::as<MatrixXd>(B_gamma_y_list_input[i]);
    }
    // Does the advection design actually vary over time? With the default
    // B_gamma_x / B_gamma_y it does not -- the same matrix is repeated at every
    // time node -- and then so does every spatial block of K.
    gamma_time_invariant = true;
    for (int i = 1; i < nt - 1 && gamma_time_invariant; i++)
      gamma_time_invariant =
          B_gamma_x_list[i].rows() == B_gamma_x_list[0].rows() &&
          B_gamma_x_list[i].cols() == B_gamma_x_list[0].cols() &&
          B_gamma_y_list[i].rows() == B_gamma_y_list[0].rows() &&
          B_gamma_y_list[i].cols() == B_gamma_y_list[0].cols() &&
          B_gamma_x_list[i] == B_gamma_x_list[0] &&
          B_gamma_y_list[i] == B_gamma_y_list[0];
  }
}

void Spacetime::build_KZ(const VectorXd &theta_K) {
  K.setZero();
  double c = exp(theta_K[0]);
  double kappa = exp(theta_K[1]);
  // spatial operator of the FIRST interior block, needed for the stationary
  // initial condition below
  SparseMatrix<double> Ls_first;

  if (!fix_gamma) {
    if (shared_theta_gamma) {
      theta_gamma_x = theta_K.segment(2, n_theta_gamma_x);
      theta_gamma_y = theta_gamma_x; // Use same parameter for both x and y
    } else {
      theta_gamma_x = theta_K.segment(2, n_theta_gamma_x);
      theta_gamma_y = theta_K.segment(2 + n_theta_gamma_x, n_theta_gamma_y);
    }

    // One spatial block of K for the advection field at time node i. When the
    // advection design does not vary over time this is called once and the
    // result reused for every block, rather than rebuilt nt - 1 times: each
    // call is several diagonal-scaled sparse products, and build_KZ() itself
    // runs once per parameter per iteration under numeric differencing, so the
    // repeat was the dominant cost of assembling K.
    auto build_Ls = [&](int i) {
      const VectorXd gamma_x = B_gamma_x_list[i] * theta_gamma_x;
      const VectorXd gamma_y = B_gamma_y_list[i] * theta_gamma_y;

      SparseMatrix<double, 0, int> Ls =
          kappa * kappa * Cs + lambda * Gs +
          (gamma_x.asDiagonal() * Bx + gamma_y.asDiagonal() * By);

      if (alpha == 4)
        Ls = Ls * Cs_diag.cwiseInverse().asDiagonal() * Ls.transpose();

      if (stabilization && !(gamma_x.norm() < 1e-8 && gamma_y.norm() < 1e-8)) {
        VectorXd gamma_xx = gamma_x.array().square();
        VectorXd gamma_yy = gamma_y.array().square();
        VectorXd gamma_xy = gamma_x.array() * gamma_y.array();

        SparseMatrix<double, 0, int> Si =
            gamma_xx.asDiagonal() * Hxx * gamma_xx.asDiagonal() +
            gamma_yy.asDiagonal() * Hyy * gamma_yy.asDiagonal() +
            gamma_xy.asDiagonal() * (Hxy + Hyx) * gamma_xy.asDiagonal();

        const double gamma_norm =
            sqrt((gamma_x.array().square() + gamma_y.array().square()).sum());
        Ls = Ls + Cs_diag.asDiagonal() * Si / gamma_norm;
      }
      return Ls;
    };

    if (gamma_time_invariant) {
      const SparseMatrix<double, 0, int> Ls = build_Ls(0);
      for (int i = 1; i < nt; ++i)
        setSparseBlock(&K, i * Ls.rows(), i * Ls.cols(), Ls);
      Ls_first = Ls;
    } else {
      for (int i = 1; i < nt; ++i) {
        const SparseMatrix<double, 0, int> Ls = build_Ls(i - 1);
        setSparseBlock(&K, i * Ls.rows(), i * Ls.cols(), Ls);
        if (i == 1) Ls_first = Ls;
      }
    }
  } else { // fix_gamma = TRUE
    SparseMatrix<double> Ls = (kappa * kappa * Cs + lambda * Gs + Bs);
    // alpha=4, L = L %*% solve(Ct %x% Cs, L)
    if (alpha == 4)
      Ls = Ls * Cs_diag.cwiseInverse().asDiagonal() * Ls.transpose();
    for (int i = 1; i < nt; ++i) {
      setSparseBlock(&K, i * Ls.rows(), i * Ls.cols(), Ls);
    }
    Ls_first = Ls;
  }
  K = (stationary_init ? BtCs_st : BtCs) + K / c;

  // sqrt(c): the symmetric split, spatial term ~ 1/sqrt(c), variance ~ c.
  // c:       the whole factor on the temporal term, spatial term c-free.
  K = (cc_variance_free ? c : std::sqrt(c)) * K;

  // STATIONARY INITIAL CONDITION.
  // Interior block rows read  M W_t - N W_{t-1},  with
  //     M = sqrt(c) Cs + Ls / sqrt(c),      N = sqrt(c) Cs,
  // so the one-step map is A = M^-1 N. AR(1) starts from stationarity by
  // putting sqrt(1 - rho^2) in K[1,1]; the matrix statement of the same thing
  // is K_1' K_1 = M'(I - A'A)M = M'M - N'N, which reduces to exactly 1 - rho^2
  // in the scalar case. It is sparse, and it depends on cc as rho does.
  if (stationary_init && nt > 1 && Ls_first.rows() == ns_) {
    // K_1 = gamma * M, with M the interior diagonal block.
    //
    // The exact stationary block is a factor of P1 = M'M - N'N, but a factor is
    // the wrong object for this model. Under non-Gaussian noise each ROW of K
    // carries its own mixing variable, V_i ~ IG(nu h_i, nu h_i^2) with h_i the
    // mesh weight of that row (see latent.cpp, sample_V). A Cholesky factor's
    // rows are arbitrary combinations of nodes, so the noise would be attached
    // to combinations rather than to increments, and the fitted model would
    // depend on the FACTORIZATION ORDERING -- a numerical choice, not a
    // modelling one. For Gaussian noise that is invisible, since only K_1' K_1
    // enters; for NIG it changes the likelihood.
    //
    // Taking K_1 proportional to M instead keeps the first block's sparsity and
    // row structure identical to every interior block, so its rows are the same
    // kind of object and h_i keeps its meaning. The scale follows the AR(1)
    // rule by matching traces,
    //     gamma^2 = 1 - tr(N'N)/tr(M'M),
    // which is exactly 1 - rho^2 when M and N are scalars.
    // Must mirror the interior blocks under whichever scaling is in force,
    // or slice 1 is stationary for a different operator than the one it feeds.
    //   split form : M = sqrt(c) Cs + Ls/sqrt(c),  N = sqrt(c) Cs
    //   variance-free: M = c Cs + Ls,              N = c Cs
    // gamma^2 = 1 - ||N||^2/||M||^2 therefore changes with the scaling too.
    const double sq = cc_variance_free ? c : std::sqrt(c);
    SparseMatrix<double> Msp =
        cc_variance_free ? (SparseMatrix<double>)(sq * Cs + Ls_first)
                         : (SparseMatrix<double>)(sq * Cs + Ls_first / sq);
    SparseMatrix<double> Nsp = sq * Cs;
    const double mm = Msp.squaredNorm();
    double g2 = (mm > 0) ? (1.0 - Nsp.squaredNorm() / mm) : 1.0;
    if (!(g2 > 0.0))
      g2 = 1e-12; // N dominates M: the step is nearly a pure copy
    SparseMatrix<double> K1 = std::sqrt(g2) * Msp;
    {
      // Add, never insert: K_1 carries fill beyond K's pattern here too, and
      // setSparseBlock() writes through coeffRef, which is O(nnz(K)) per new
      // entry and made the assembly scale with nt.
      std::vector<Eigen::Triplet<double>> tp;
      tp.reserve(K1.nonZeros());
      for (int cc2 = 0; cc2 < K1.outerSize(); ++cc2)
        for (SparseMatrix<double>::InnerIterator it(K1, cc2); it; ++it)
          tp.emplace_back((int)it.row(), (int)it.col(), it.value());
      SparseMatrix<double> K1pad(K.rows(), K.cols());
      K1pad.setFromTriplets(tp.begin(), tp.end());
      K = K + K1pad;
    }
    // If the factorization fails, K keeps the zero first block: degenerate, but
    // no worse than before this change, and it does not abort the fit.
  }

  // SparseMatrix<double> Ls = (kappa*kappa * Cs + lambda * Gs + Bs);

  // alpha=4, L = L %*% solve(Ct %x% Cs, L)
  // if (alpha == 4)
  //   Ls = Ls * Cs_diag.cwiseInverse().asDiagonal() * Ls.transpose();
  // Ct is diagonal

  // if (method == "galerkin") {
  //   KroneckerProductSparse<SparseMatrix<double>, SparseMatrix<double> >
  //   kroneckerEigen(Ct, Ls);

  //   // update K
  //   K = kroneckerEigen.eval();
  // } else if (method == "euler") {
  //   // if (stabilization) Ls = Ls + S;
  //   // Build K = bdiag(0, Ls, ..., Ls)
  //   // 1st approach: Using Kronecker product
  //   // create diag(0, 1, 1, ..., 1) (nt-1) 1
  //   // Eigen::SparseMatrix<double> I_sparse(Ct.rows(), Ct.rows());
  //   // for (int i = 1; i < Ct.rows(); ++i) {
  //   //   I_sparse.insert(i, i) = 1;
  //   // }

  //   // KroneckerProductSparse<SparseMatrix<double>, SparseMatrix<double> >
  //   kroneckerEigen(I_sparse, Ls);
  //   // K = kroneckerEigen.eval();

  //   // 2nd approach: build Bdiag(0, Ls, ..., Ls) directly
  //   K.setZero();
  //   for (int i = 0; i < nt-1; ++i) {
  //     setSparseBlock(&K, i * Ls.rows(), i * Ls.cols(), Ls_list[i]);
  //   }
  // }
}

void Spacetime::update_dK(const VectorXd &theta_K) {}

// ---------------------------------------------------------------------------
// Closed-form derivatives of the space-time operator.
//
//   K(theta) = s * Base + (s/c) * blockdiag(0, Ls, ..., Ls) + sqrt(g2) * pad(M)
//
// with c = exp(theta_0), kappa = exp(theta_1), s = c^a (a = 1 when the cc
// factor rides entirely on the temporal term, 1/2 for the symmetric split),
// Base a constant matrix, and
//
//   Ls = kappa^2 Cs + lambda Gs + diag(gamma_x) Bx + diag(gamma_y) By,
//
// where gamma_x = B_gamma_x theta_gamma_x is LINEAR in its parameters. So Ls is
// linear in every gamma parameter, depends on theta_1 only through
// kappa^2 = e^{2 theta_1}, and is free of theta_0; the one awkward piece is the
// stationary scalar g2 = 1 - ||N||^2/||M||^2, which is smooth and handled by
// the chain rule below.
//
// The point of doing this at all: every derivative has the SAME SHAPE as K --
// a multiple of Base, plus a block-diagonal spatial part, plus a first-block
// correction -- so one assembly routine serves all of them, and all the
// theta-dependent work is ns x ns rather than (nt*ns) x (nt*ns). Numeric
// differencing instead calls build_KZ() once per parameter for dK and four
// times per parameter PAIR for d2K, each one assembling the full operator.
//
// Not covered, and handed back to numeric differencing: alpha == 4 (Ls becomes
// quadratic in itself) and stabilization (nonlinear in gamma, with a norm).
// ---------------------------------------------------------------------------

namespace {
// Frobenius inner product of two sparse matrices.
inline double fro_dot(const SparseMatrix<double> &A,
                      const SparseMatrix<double> &B) {
  return A.cwiseProduct(B).sum();
}
} // namespace

void Spacetime::assemble_shaped(SparseMatrix<double> &out, double base_coef,
                                const std::vector<SparseMatrix<double>> &blk,
                                bool uniform, const SparseMatrix<double> &first,
                                bool has_first) const {
  const SparseMatrix<double, 0, int> &Base = stationary_init ? BtCs_st : BtCs;
  std::vector<Eigen::Triplet<double>> tp;
  size_t cap = has_first ? (size_t)first.nonZeros() : 0;
  for (int i = 1; i < nt; ++i)
    cap += (size_t)blk[uniform ? 0 : i - 1].nonZeros();
  if (base_coef != 0.0)
    cap += (size_t)Base.nonZeros();
  tp.reserve(cap);
  if (base_coef != 0.0)
    for (int c = 0; c < Base.outerSize(); ++c)
      for (SparseMatrix<double, 0, int>::InnerIterator it(Base, c); it; ++it)
        tp.emplace_back((int)it.row(), c, base_coef * it.value());
  for (int i = 1; i < nt; ++i) {
    const SparseMatrix<double> &B = blk[uniform ? 0 : i - 1];
    const int off = i * ns_;
    for (int c = 0; c < B.outerSize(); ++c)
      for (SparseMatrix<double>::InnerIterator it(B, c); it; ++it)
        tp.emplace_back((int)it.row() + off, c + off, it.value());
  }
  if (has_first)
    for (int c = 0; c < first.outerSize(); ++c)
      for (SparseMatrix<double>::InnerIterator it(first, c); it; ++it)
        tp.emplace_back((int)it.row(), c, it.value());
  out.resize(K.rows(), K.cols());
  out.setFromTriplets(tp.begin(), tp.end());
  out.makeCompressed();
}

// dLs/dtheta_m for the interior block at time node i, complete (the kappa
// chain rule included). False means identically zero.
bool Spacetime::dLs_block(int m, int i, double k2,
                          SparseMatrix<double> &out) const {
  if (m == 1) { // d(kappa^2)/d theta_1 = 2 kappa^2
    out = (2.0 * k2) * Cs;
    return true;
  }
  if (fix_gamma || m < 2)
    return false;
  const int g = m - 2;
  if (shared_theta_gamma) {
    if (g >= n_theta_gamma_x)
      return false;
    // One parameter drives both components, so both terms contribute.
    out = SparseMatrix<double>(
              (SparseMatrix<double>)(B_gamma_x_list[i].col(g).asDiagonal() * Bx)) +
          SparseMatrix<double>(
              (SparseMatrix<double>)(B_gamma_y_list[i].col(g).asDiagonal() * By));
    return true;
  }
  if (g < n_theta_gamma_x) {
    out = B_gamma_x_list[i].col(g).asDiagonal() * Bx;
    return true;
  }
  const int gy = g - n_theta_gamma_x;
  if (gy >= n_theta_gamma_y)
    return false;
  out = B_gamma_y_list[i].col(gy).asDiagonal() * By;
  return true;
}

bool Spacetime::d2Ls_block(int m, int n, int i, double k2,
                           SparseMatrix<double> &out) const {
  (void)i;
  // Ls is linear in every gamma parameter and free of theta_0, so the only
  // second derivative that survives is d2(kappa^2)/d theta_1^2 = 4 kappa^2.
  if (m == 1 && n == 1) {
    out = (4.0 * k2) * Cs;
    return true;
  }
  return false;
}

// Gather what both derivative routines need from theta. Returns false for a
// configuration the closed form does not cover, leaving numeric differencing
// in place.
// NGME_SPACETIME_NO_ANALYTIC=1 forces numeric differencing, as an escape hatch
// and so the two routes can be compared in one binary.
static bool spacetime_analytic_disabled() {
  static const bool v = [] {
    const char *e = std::getenv("NGME_SPACETIME_NO_ANALYTIC");
    return e && *e && std::string(e) != "0";
  }();
  return v;
}

bool Spacetime::analytic_state(const VectorXd &theta_K, double &c_out,
                               double &k2_out, double &a_out, double &s_out,
                               double &q_out, bool &uniform_out,
                               std::vector<SparseMatrix<double>> &Ls,
                               bool &has_stat, SparseMatrix<double> &M,
                               double &mm, double &nn, double &u) const {
  if (spacetime_analytic_disabled())
    return false;
  if (alpha != 2 || stabilization || nt <= 1 || ns_ <= 0)
    return false;
  if ((int)theta_K.size() != n_theta_K || n_theta_K < 2)
    return false;

  const double c = std::exp(theta_K[0]);
  const double kappa = std::exp(theta_K[1]);
  const double k2 = kappa * kappa;
  const double a = cc_variance_free ? 1.0 : 0.5;
  const double s = std::pow(c, a);
  const bool uniform = fix_gamma || gamma_time_invariant;
  const int nblk = uniform ? 1 : nt - 1;

  // The interior spatial blocks, formed exactly as build_KZ() forms them.
  Ls.assign(nblk, SparseMatrix<double>());
  for (int b = 0; b < nblk; ++b) {
    if (fix_gamma) {
      Ls[b] = k2 * Cs + lambda * Gs + Bs;
    } else {
      VectorXd tgx = theta_K.segment(2, n_theta_gamma_x);
      VectorXd tgy = shared_theta_gamma
                         ? tgx
                         : (VectorXd)theta_K.segment(2 + n_theta_gamma_x,
                                                     n_theta_gamma_y);
      const VectorXd gx = B_gamma_x_list[b] * tgx;
      const VectorXd gy = B_gamma_y_list[b] * tgy;
      Ls[b] =
          k2 * Cs + lambda * Gs + (gx.asDiagonal() * Bx + gy.asDiagonal() * By);
    }
  }

  has_stat = stationary_init;
  if (has_stat) {
    M = cc_variance_free ? (SparseMatrix<double>)(s * Cs + Ls[0])
                         : (SparseMatrix<double>)(s * Cs + Ls[0] / s);
    SparseMatrix<double> N = s * Cs;
    mm = M.squaredNorm();
    nn = N.squaredNorm();
    if (!(mm > 0))
      return false;
    const double g2 = 1.0 - nn / mm;
    // build_KZ() clamps a non-positive g2 to a constant, which is a kink;
    // hand that case to numeric differencing rather than differentiate it.
    if (!(g2 > 0.0))
      return false;
    u = std::sqrt(g2);
  }
  c_out = c; k2_out = k2; a_out = a; s_out = s; q_out = s / c;
  uniform_out = uniform;
  return true;
}

bool Spacetime::update_dKdZ(const VectorXd &theta_K) {
  double c, k2, a, s, q, mm = 0, nn = 0, u = 1;
  bool uniform = false, has_stat = false;
  std::vector<SparseMatrix<double>> Ls;
  SparseMatrix<double> M;
  if (!analytic_state(theta_K, c, k2, a, s, q, uniform, Ls, has_stat, M, mm, nn,
                      u))
    return false;
  const int nblk = (int)Ls.size();

  if ((int)dK.size() != n_theta_K)
    dK.assign(n_theta_K, SparseMatrix<double>(K.rows(), K.cols()));
  if ((int)dZ.size() != n_theta_K)
    dZ.assign(n_theta_K, SparseMatrix<double>(K.rows(), K.cols()));

  std::vector<SparseMatrix<double>> blk(nblk);
  SparseMatrix<double> first, dM;
  for (int m = 0; m < n_theta_K; ++m) {
    // build_KZ() never writes Z, so Z does not depend on theta at all.
    dZ[m].setZero();

    double base_coef = 0.0;
    if (m == 0) {
      // K = s Base + q P, with s = c^a and q = c^(a-1).
      base_coef = a * s;
      for (int b = 0; b < nblk; ++b)
        blk[b] = ((a - 1.0) * q) * Ls[b];
      if (has_stat)
        dM = cc_variance_free
                 ? (SparseMatrix<double>)((a * s) * Cs)
                 : (SparseMatrix<double>)((a * s) * Cs - (a / s) * Ls[0]);
    } else {
      bool any = false;
      for (int b = 0; b < nblk; ++b) {
        SparseMatrix<double> d;
        if (dLs_block(m, b, k2, d)) {
          blk[b] = q * d;
          any = true;
        } else {
          blk[b] = SparseMatrix<double>(ns_, ns_);
        }
      }
      if (!any) { // this parameter does not enter K
        dK[m].setZero();
        continue;
      }
      if (has_stat) {
        SparseMatrix<double> d0;
        if (dLs_block(m, 0, k2, d0))
          dM = cc_variance_free ? d0 : (SparseMatrix<double>)(d0 / s);
        else
          dM = SparseMatrix<double>(ns_, ns_);
      }
    }

    if (has_stat) {
      // g2 = 1 - nn/mm, K_1 = sqrt(g2) M. N = s Cs, so nn = s^2 ||Cs||^2
      // depends on theta_0 alone and d nn / d theta_0 = 2 a nn.
      const double dnn = (m == 0) ? 2.0 * a * nn : 0.0;
      const double dmm = 2.0 * fro_dot(M, dM);
      const double dr = dnn / mm - nn * dmm / (mm * mm);
      const double du = -dr / (2.0 * u); // d sqrt(1-r) = -dr / (2 sqrt(1-r))
      first = du * M + u * dM;
    }
    assemble_shaped(dK[m], base_coef, blk, uniform, first, has_stat);
  }
  return true;
}

bool Spacetime::update_d2Kd2Z(const VectorXd &theta_K) {
  double c, k2, a, s, q, mm = 0, nn = 0, u = 1;
  bool uniform = false, has_stat = false;
  std::vector<SparseMatrix<double>> Ls;
  SparseMatrix<double> M;
  if (!analytic_state(theta_K, c, k2, a, s, q, uniform, Ls, has_stat, M, mm, nn,
                      u))
    return false;
  const int nblk = (int)Ls.size();

  // First derivatives are needed by every pair, so they are formed once here
  // rather than inside the double loop.
  std::vector<SparseMatrix<double>> dM(n_theta_K);
  std::vector<double> dmm(n_theta_K, 0.0), dnn(n_theta_K, 0.0),
      dr(n_theta_K, 0.0);
  std::vector<std::vector<SparseMatrix<double>>> dblk(
      n_theta_K, std::vector<SparseMatrix<double>>(nblk));
  for (int m = 0; m < n_theta_K; ++m) {
    if (m == 0) {
      for (int b = 0; b < nblk; ++b)
        dblk[m][b] = ((a - 1.0) * q) * Ls[b];
      if (has_stat)
        dM[m] = cc_variance_free
                    ? (SparseMatrix<double>)((a * s) * Cs)
                    : (SparseMatrix<double>)((a * s) * Cs - (a / s) * Ls[0]);
      dnn[m] = 2.0 * a * nn;
    } else {
      for (int b = 0; b < nblk; ++b) {
        SparseMatrix<double> d;
        if (dLs_block(m, b, k2, d))
          dblk[m][b] = q * d;
        else
          dblk[m][b] = SparseMatrix<double>(ns_, ns_);
      }
      if (has_stat) {
        SparseMatrix<double> d0;
        if (dLs_block(m, 0, k2, d0))
          dM[m] = cc_variance_free ? d0 : (SparseMatrix<double>)(d0 / s);
        else
          dM[m] = SparseMatrix<double>(ns_, ns_);
      }
    }
    if (has_stat) {
      dmm[m] = 2.0 * fro_dot(M, dM[m]);
      dr[m] = dnn[m] / mm - nn * dmm[m] / (mm * mm);
    }
  }

  if ((int)d2K.size() != n_theta_K)
    d2K.assign(n_theta_K,
               std::vector<SparseMatrix<double>>(
                   n_theta_K, SparseMatrix<double>(K.rows(), K.cols())));
  if ((int)d2Z.size() != n_theta_K)
    d2Z.assign(n_theta_K, std::vector<SparseMatrix<double, 0, int>>(
                              n_theta_K,
                              SparseMatrix<double, 0, int>(K.rows(), K.cols())));

  std::vector<SparseMatrix<double>> blk(nblk);
  SparseMatrix<double> first, d2M;
  for (int m = 0; m < n_theta_K; ++m) {
    for (int n2 = m; n2 < n_theta_K; ++n2) {
      d2Z[m][n2].setZero();
      if (n2 != m)
        d2Z[n2][m].setZero();

      // ---- block-diagonal part ----
      double base_coef = 0.0;
      if (m == 0 && n2 == 0) {
        base_coef = a * a * s;
        for (int b = 0; b < nblk; ++b)
          blk[b] = ((a - 1.0) * (a - 1.0) * q) * Ls[b];
      } else if (m == 0) {
        // d/dtheta_0 of (q dLs_n) = (a-1) q dLs_n
        for (int b = 0; b < nblk; ++b)
          blk[b] = (a - 1.0) * dblk[n2][b];
      } else {
        for (int b = 0; b < nblk; ++b) {
          SparseMatrix<double> d2;
          if (d2Ls_block(m, n2, b, k2, d2))
            blk[b] = q * d2;
          else
            blk[b] = SparseMatrix<double>(ns_, ns_);
        }
      }

      // ---- stationary first block ----
      if (has_stat) {
        if (m == 0 && n2 == 0) {
          d2M = cc_variance_free
                    ? (SparseMatrix<double>)((a * a * s) * Cs)
                    : (SparseMatrix<double>)((a * a * s) * Cs +
                                             (a * a / s) * Ls[0]);
        } else if (m == 0) {
          // M_n is dLs_n (cc-free) or dLs_n / s; only the latter sees theta_0.
          d2M = cc_variance_free ? SparseMatrix<double>(ns_, ns_)
                                 : (SparseMatrix<double>)((-a) * dM[n2]);
        } else {
          SparseMatrix<double> d2;
          if (d2Ls_block(m, n2, 0, k2, d2))
            d2M = cc_variance_free ? d2 : (SparseMatrix<double>)(d2 / s);
          else
            d2M = SparseMatrix<double>(ns_, ns_);
        }
        const double d2nn = (m == 0 && n2 == 0) ? 4.0 * a * a * nn : 0.0;
        const double d2mm =
            2.0 * fro_dot(dM[m], dM[n2]) + 2.0 * fro_dot(M, d2M);
        const double d2r = d2nn / mm -
                           (dnn[m] * dmm[n2] + dnn[n2] * dmm[m]) / (mm * mm) -
                           nn * d2mm / (mm * mm) +
                           2.0 * nn * dmm[m] * dmm[n2] / (mm * mm * mm);
        // u = sqrt(1 - r):  u_m = -r_m/(2u),  u_mn = -r_mn/(2u) - r_m r_n/(4u^3)
        const double du_m = -dr[m] / (2.0 * u);
        const double du_n = -dr[n2] / (2.0 * u);
        const double d2u =
            -d2r / (2.0 * u) - dr[m] * dr[n2] / (4.0 * u * u * u);
        first = d2u * M + du_m * dM[n2] + du_n * dM[m] + u * d2M;
      }

      assemble_shaped(d2K[m][n2], base_coef, blk, uniform, first, has_stat);
      if (n2 != m)
        d2K[n2][m] = d2K[m][n2];
    }
  }
  return true;
}

// Exact traces from the diagonal blocks alone.
//
// With the stationary initial condition K is block lower-BIDIAGONAL: block row
// 0 holds only K_1, and row t holds -N at (t, t-1) and M_t at (t, t). K^-1 is
// then block lower-triangular, so in
//     tr(K^-1 dK) = sum_{s,t} tr[ (K^-1)_{t,s} (dK)_{s,t} ]
// the subdiagonal term would need (K^-1)_{t,t+1}, which is zero above the
// diagonal. Only the diagonal blocks survive:
//     tr(K^-1 dK) = sum_t tr( K_tt^-1 (dK)_tt ).
// This needs no probes and, more importantly, no factorization of the full
// (nt*ns) x (nt*ns) operator -- Operator::update_all builds that solver only
// when no structural shortcut supplies the traces.
//
// Under fix_gamma, or free gamma with a time-invariant advection design, every
// interior block is the SAME matrix, so one factorization and one block trace
// serve all nt - 1 of them however long the series.
// tr(M_b^-1 B) from the selected inverse of diagonal block b.
//
// selinv_trace() computes tr(Q^-1 B) for the Q this solver factorized, which is
// exactly the block, and returns false when B touches an entry the factor's
// pattern does not cover, so the failure mode is a fallback, not a silently
// wrong trace.
bool Spacetime::block_trace_selinv(int b,
                                   const SparseMatrix<double, 0, int> &B,
                                   double &out) {
  if (b < 0 || b >= (int)blk_solver_.size() || !blk_solver_[b])
    return false;
  return blk_solver_[b]->selinv_trace(B, out);
}

bool Spacetime::compute_traces_structured(const VectorXd &theta,
                                          const UpdateOptions &opts) {
  if (!stationary_init || !opts.compute_trace)
    return false;
  if (n_theta_K <= 0 || (int)dK.size() < n_theta_K || ns_ <= 0 || nt <= 0)
    return false;
  if (K.rows() != (long long)ns_ * nt || K.cols() != K.rows())
    return false;
  // NGME_SPACETIME_NO_HK=1 restores the pre-Hessian behaviour  so the
  // two can be A/B'd in one binary instead of across two builds.
  static const bool no_hk = [] {
    const char *e = std::getenv("NGME_SPACETIME_NO_HK");
    return e && *e && std::string(e) != "0";
  }();
  if (no_hk && opts.compute_HK_trace)
    return false;

  const int n = ns_;
  auto get_block = [&](const SparseMatrix<double> &M, int r0, int c0) {
    std::vector<Eigen::Triplet<double>> tp;
    for (int c = c0; c < c0 + n; ++c)
      for (SparseMatrix<double>::InnerIterator it(M, c); it; ++it)
        if (it.row() >= r0 && it.row() < r0 + n)
          tp.emplace_back((int)it.row() - r0, c - c0, it.value());
    SparseMatrix<double> B(n, n);
    B.setFromTriplets(tp.begin(), tp.end());
    B.makeCompressed();
    return B;
  };

  // Interior blocks coincide unless the advection design varies over time.
  const bool uniform_interior = fix_gamma || gamma_time_invariant;
  const int n_distinct = uniform_interior ? std::min(2, nt) : nt;

  std::vector<int> rep(n_distinct);
  std::vector<SparseMatrix<double>> Kblk(n_distinct);
  for (int b = 0; b < n_distinct; ++b) {
    rep[b] = (b == 0) ? 0 : b; // block index this factorization represents
    Kblk[b] = get_block(K, rep[b] * n, rep[b] * n);
  }

  // The derivative blocks, hoisted out of the parameter loop below: each was
  // being rebuilt from triplets on every use.
  auto dK_block = [&](int j, int b) { return get_block(dK[j], rep[b] * n, rep[b] * n); };

  // tr(M^-1 B) = sum_{(i,j) in nnz(B)} (M^-1)_{ji} B_{ij}. Given a dense M_b^-1
  // this is O(nnz(B)) and needs no solve, which is the whole point of building
  // the inverse once per block below instead of solving once per quantity.
  auto trace_with_inv = [](const MatrixXd &Zi, const SparseMatrix<double> &B) {
    double s = 0.0;
    for (int c = 0; c < B.outerSize(); ++c)
      for (SparseMatrix<double>::InnerIterator it(B, c); it; ++it)
        s += Zi(c, (int)it.row()) * it.value();
    return s;
  };

  // Dense M_b^-1, built only when the Hessian is wanted (see below). Empty
  // otherwise: the first-order traces alone do not justify an ns x ns inverse.
  std::vector<MatrixXd> Zinv;

  // ---------------------------------------------------------------------
  // Preferred route: selected inverse.
  //
  //     tr(M^-1 B) = sum_{(i,j) in nnz(B)} (M^-1)_{ji} B_{ij}
  //
  // so M^-1 is needed only on the pattern of B, never as a full inverse. The
  // previous implementation instead densified B (ns x ns) and solved against
  // all ns right-hand sides to keep only the ns diagonal entries of the
  // result -- once PER PARAMETER, per iteration. Measured against a single
  // factorization of the same block that is ~40x at ns=1126 and ~144x at
  // ns=2277, and it grows like ns^2 while the factorization does not.
  //
  // A selected inverse costs a small multiple of one factorization and serves
  // every parameter, so the work becomes O(n_distinct) factorizations instead
  // of O(n_distinct * n_theta_K) dense solves.
  //
  // It needs a symmetric positive-definite block: that holds without advection
  // (M = sqrt(cc) Cs + Ls/sqrt(cc), both symmetric), and not with a free
  // advection field, which is why the fallback below is kept rather than
  // removed. sparse_llt_solver::selinv_trace also reports failure when a
  // derivative entry falls outside the factor pattern, so a wrong answer is
  // not among the outcomes.
  // ---------------------------------------------------------------------
  // NGME_SPACETIME_NO_SELINV=1 forces the old dense-solve route. Kept as an
  // escape hatch, and so the two can be compared in one binary.
  static const bool selinv_disabled = [] {
    const char *e = std::getenv("NGME_SPACETIME_NO_SELINV");
    return e && *e && std::string(e) != "0";
  }();
  bool use_selinv = !selinv_disabled && (blk_selinv_state_ != 0);
  if (use_selinv) {
    if ((int)blk_solver_.size() != n_distinct) {
      blk_solver_.clear();
      blk_solver_.resize(n_distinct);
      blk_nnz_.assign(n_distinct, -1);
    }
    for (int b = 0; b < n_distinct && use_selinv; ++b) {
      const SparseMatrix<double> &A = Kblk[b];
      // Symmetry is a property of the assembled block, so it is checked on the
      // numbers rather than inferred from fix_gamma.
      SparseMatrix<double> D = SparseMatrix<double>(A.transpose()) - A;
      if (D.nonZeros() > 0 &&
          D.coeffs().abs().maxCoeff() >
              1e-12 * std::max(1.0, A.coeffs().abs().maxCoeff())) {
        selinv_fail_ = 1; // asymmetric block
        use_selinv = false;
        break;
      }
      if (!blk_solver_[b]) {
        blk_solver_[b].reset(new sparse_llt_solver());
        blk_nnz_[b] = -1;
      }
      if (blk_nnz_[b] != (long long)A.nonZeros()) {
        // stype 0 (SimplicialLLT) so the selected inverse reuses this factor
        // instead of building a second one; see ensure_selinv_factor().
        blk_solver_[b]->init(n, /*Ntrace*/ 1, /*symmetric*/ true, /*stype*/ 0);
        blk_solver_[b]->analyze(A);
        blk_nnz_[b] = (long long)A.nonZeros();
      }
      blk_solver_[b]->compute(A);
      if (!blk_solver_[b]->factorization_success()) {
        selinv_fail_ = 2; // Cholesky rejected the block
        use_selinv = false;
      }
    }
  }

  if (trace_vals.size() != n_theta_K)
    trace_vals = VectorXd::Zero(n_theta_K);

  if (use_selinv) {
    VectorXd tv = VectorXd::Zero(n_theta_K);
    for (int j = 0; j < n_theta_K && use_selinv; ++j) {
      if (!opts.fix_mask_thetaK.empty() && opts.fix_mask_thetaK[j])
        continue;
      double tr = 0.0, t0 = 0.0;
      if (!block_trace_selinv(0, dK_block(j, 0), t0)) { selinv_fail_ = 3; use_selinv = false; break; }
      tr += t0;
      if (uniform_interior) {
        if (nt > 1) {
          const int b1 = std::min(1, n_distinct - 1);
          double t1 = 0.0;
          if (!block_trace_selinv(b1, dK_block(j, b1), t1)) { selinv_fail_ = 3; use_selinv = false; break; }
          tr += (double)(nt - 1) * t1;
        }
      } else {
        for (int t = 1; t < nt; ++t) {
          double tt = 0.0;
          if (!block_trace_selinv(t, dK_block(j, t), tt)) { selinv_fail_ = 3; use_selinv = false; break; }
          tr += tt;
        }
      }
      tv(j) = tr;
    }
    if (use_selinv) {
      // NGME_SPACETIME_TRACE_CHECK=1 recomputes every trace the old way and
      // compares. This checks the quantity itself on every iteration, which is
      // a far sharper test than comparing converged estimates.
      static const bool check = [] {
        const char *e = std::getenv("NGME_SPACETIME_TRACE_CHECK");
        return e && *e && std::string(e) != "0";
      }();
      if (check) {
        std::vector<Eigen::SparseLU<SparseMatrix<double>>> lu(n_distinct);
        bool lu_ok = true;
        for (int b = 0; b < n_distinct && lu_ok; ++b) {
          lu[b].compute(Kblk[b]);
          lu_ok = (lu[b].info() == Eigen::Success);
        }
        if (lu_ok) {
          double worst = 0.0;
          int worst_j = -1;
          for (int j = 0; j < n_theta_K; ++j) {
            if (!opts.fix_mask_thetaK.empty() && opts.fix_mask_thetaK[j])
              continue;
            double ref = 0.0;
            if (uniform_interior) {
              ref += MatrixXd(lu[0].solve(MatrixXd(dK_block(j, 0)))).trace();
              if (nt > 1) {
                const int b1 = std::min(1, n_distinct - 1);
                ref += (double)(nt - 1) *
                       MatrixXd(lu[b1].solve(MatrixXd(dK_block(j, b1)))).trace();
              }
            } else {
              for (int t = 0; t < nt; ++t)
                ref += MatrixXd(lu[t].solve(MatrixXd(dK_block(j, t)))).trace();
            }
            const double rel =
                std::abs(tv(j) - ref) / std::max(1.0, std::abs(ref));
            if (rel > worst) { worst = rel; worst_j = j; }
          }
          // fprintf, not ngme_io::out(): that discards anything written off
          // the main thread, which would hide exactly the calls most worth
          // seeing. This is a debug path behind an env var, not normal output.
          std::fprintf(stderr,
                       "[spacetime trace check] call %d  max rel diff %.3e "
                       "(theta_K[%d])\n",
                       ++selinv_calls_, worst, worst_j);
        }
      }
      trace_vals = tv;
      blk_selinv_state_ = 1;
    }
  }

  // Set when this call has already factorized every diagonal block into
  // hk_lu_, so the Hessian section below does not redo the work.
  bool hk_lu_ready = false;

  if (!use_selinv) {
    if (std::getenv("NGME_SPACETIME_TRACE_CHECK"))
      std::fprintf(stderr,
                   "[spacetime trace check] FALLBACK to dense LU path "
                   "(reason %d: 1=asymmetric block, 2=chol failed, "
                   "3=derivative outside factor pattern)\n", selinv_fail_);
    // Fallback: the original dense-solve route. Reached for a non-symmetric
    // block (free advection), a block the Cholesky rejects, or a derivative
    // reaching outside the factor pattern.
    blk_selinv_state_ = 0;
    blk_solver_.clear();
    blk_nnz_.clear();
    // Factorize into hk_lu_ rather than a local vector: the Hessian section
    // below needs exactly these factorizations, and on this path (asymmetric
    // block -- free advection) there is no selinv solver for it to use.
    if ((int)hk_lu_.size() != n_distinct) {
      hk_lu_.clear();
      for (int b = 0; b < n_distinct; ++b)
        hk_lu_.emplace_back(new Eigen::SparseLU<SparseMatrix<double>>());
    }
    for (int b = 0; b < n_distinct; ++b) {
      hk_lu_[b]->compute(Kblk[b]);
      if (hk_lu_[b]->info() != Eigen::Success)
        return false; // fall back to the generic path rather than guess
    }
    hk_lu_ready = true;

    // When the Hessian is also wanted, one dense inverse per block replaces
    // every dense solve on this path. The n_theta_K first-order solves here
    // and the n_theta_K + n_theta_K*(n_theta_K+1)/2 solves the Hessian section
    // would otherwise do. At n_theta_K = 4 that is 18 solves down to 1.
    if (opts.compute_HK_trace) {
      Zinv.resize(n_distinct);
      MatrixXd I = MatrixXd::Identity(n, n);
      for (int b = 0; b < n_distinct; ++b)
        Zinv[b] = hk_lu_[b]->solve(I);
    }

    trace_vals.setZero();
    for (int j = 0; j < n_theta_K; ++j) {
      if (!opts.fix_mask_thetaK.empty() && opts.fix_mask_thetaK[j])
        continue;
      double tr = 0.0;
      if (!Zinv.empty()) {
        if (uniform_interior) {
          tr += trace_with_inv(Zinv[0], dK_block(j, 0));
          if (nt > 1) {
            const int b1 = std::min(1, n_distinct - 1);
            tr += (double)(nt - 1) * trace_with_inv(Zinv[b1], dK_block(j, b1));
          }
        } else {
          for (int t = 0; t < nt; ++t)
            tr += trace_with_inv(Zinv[t], dK_block(j, t));
        }
        // NGME_SPACETIME_TRACE_CHECK=1 also covers THIS route: the inverse
        // reduction is a different computation from the dense solves below, so
        // it gets compared against them rather than assumed equivalent.
        if (std::getenv("NGME_SPACETIME_TRACE_CHECK")) {
          double ref = 0.0;
          if (uniform_interior) {
            ref += MatrixXd(hk_lu_[0]->solve(MatrixXd(dK_block(j, 0)))).trace();
            if (nt > 1) {
              const int b1 = std::min(1, n_distinct - 1);
              ref += (double)(nt - 1) *
                     MatrixXd(hk_lu_[b1]->solve(MatrixXd(dK_block(j, b1)))).trace();
            }
          } else {
            for (int t = 0; t < nt; ++t)
              ref += MatrixXd(hk_lu_[t]->solve(MatrixXd(dK_block(j, t)))).trace();
          }
          std::fprintf(stderr,
                       "[spacetime inv trace check] theta_K[%d] rel diff %.3e\n",
                       j, std::abs(tr - ref) / std::max(1.0, std::abs(ref)));
        }
        trace_vals(j) = tr;
        continue;
      }
      if (uniform_interior) {
        // block 0 once, then one interior block scaled by how many there are
        MatrixXd X0 = hk_lu_[0]->solve(MatrixXd(dK_block(j, 0)));
        tr += X0.trace();
        if (nt > 1) {
          const int b1 = std::min(1, n_distinct - 1);
          MatrixXd X1 = hk_lu_[b1]->solve(MatrixXd(dK_block(j, b1)));
          tr += (double)(nt - 1) * X1.trace();
        }
      } else {
        for (int t = 0; t < nt; ++t) {
          MatrixXd Xt = hk_lu_[t]->solve(MatrixXd(dK_block(j, t)));
          tr += Xt.trace();
        }
      }
      trace_vals(j) = tr;
    }
  }

  // ---------------------------------------------------------------------
  // Hessian traces, from the same block structure.
  //
  // K is block lower-bidiagonal, so K^-1 is block lower-triangular and both
  // A = K^-1 dK_k and B = K^-1 dK_j are block lower-triangular. In
  //     tr(AB) = sum_{s,t} tr(A_ts B_st)
  // a term with t > s needs B_st, which is above the diagonal and therefore
  // zero, so only the diagonal blocks survive. And
  //     A_tt = sum_u (K^-1)_{tu} (dK_k)_{ut},
  // where K^-1 lower forces u <= t while dK_k block lower-bidiagonal allows
  // only u in {t, t+1}; the two leave u = t. Hence
  //     tr(K^-1 dK_k K^-1 dK_j) = sum_t tr( M_t^-1 (dM_k)_t M_t^-1 (dM_j)_t )
  //     tr(K^-1 d2K_jk)         = sum_t tr( M_t^-1 (d2M_jk)_t ),
  // the same reduction the first-order traces above already use, and resting
  // on the same block-lower-triangular property checked below.
  //
  // X_j = M^-1 (dM_j) is formed once per parameter rather than per pair: the
  // pair then costs only tr(X_k X_j), a dense reduction.
  // ---------------------------------------------------------------------
  if (opts.compute_HK_trace) {
    if (HK_trace.rows() != n_theta_K || HK_trace.cols() != n_theta_K)
      HK_trace = MatrixXd::Zero(n_theta_K, n_theta_K);
    else
      HK_trace.setZero();

    // Everything below is read off ONE dense M_b^-1 per block: X_j is a
    // dense-times-sparse product and the d2 trace is an O(nnz) reduction, so
    // the whole Hessian costs n_distinct solves rather than one per quantity.
    // The selinv path takes its inverse from the Cholesky factor it already
    // built; only the dense fallback needs an LU, and it built one above.
    if (Zinv.empty()) {
      Zinv.resize(n_distinct);
      MatrixXd I = MatrixXd::Identity(n, n);
      if (use_selinv) {
        for (int b = 0; b < n_distinct; ++b) {
          MatrixXd rhs = I;
          Zinv[b] = blk_solver_[b]->solve(rhs);
        }
      } else {
        if (!hk_lu_ready) {
          if ((int)hk_lu_.size() != n_distinct) {
            hk_lu_.clear();
            for (int b = 0; b < n_distinct; ++b)
              hk_lu_.emplace_back(new Eigen::SparseLU<SparseMatrix<double>>());
          }
          for (int b = 0; b < n_distinct; ++b) {
            hk_lu_[b]->compute(Kblk[b]);
            if (hk_lu_[b]->info() != Eigen::Success)
              return false;
          }
          hk_lu_ready = true;
        }
        for (int b = 0; b < n_distinct; ++b)
          Zinv[b] = hk_lu_[b]->solve(I);
      }
    }

    // X[b][j] = M_b^-1 (dM_j)_b
    std::vector<std::vector<MatrixXd>> X(
        n_distinct, std::vector<MatrixXd>(n_theta_K));
    for (int b = 0; b < n_distinct; ++b)
      for (int j = 0; j < n_theta_K; ++j)
        X[b][j] = Zinv[b] * dK_block(j, b);

    // tr(X_k X_j) = sum_{ab} (X_k)_{ab} (X_j)_{ba}
    auto tr_prod = [](const MatrixXd &Xk, const MatrixXd &Xj) {
      return (Xk.array() * Xj.transpose().array()).sum();
    };

    for (int j = 0; j < n_theta_K; ++j) {
      for (int k = j; k < n_theta_K; ++k) {
        double t1 = 0.0, t2 = 0.0;
        auto add_block = [&](int b, double w) {
          t1 += w * tr_prod(X[b][k], X[b][j]);
          const SparseMatrix<double, 0, int> &d2 = get_d2K(j, k);
          if (d2.rows() == (long long)ns_ * nt)
            t2 += w * trace_with_inv(Zinv[b],
                                     get_block(d2, rep[b] * n, rep[b] * n));
        };
        if (uniform_interior) {
          add_block(0, 1.0);
          if (nt > 1)
            add_block(std::min(1, n_distinct - 1), (double)(nt - 1));
        } else {
          for (int t = 0; t < nt; ++t)
            add_block(t, 1.0);
        }
        const double val = -t1 + t2;
        HK_trace(j, k) = val;
        if (j != k)
          HK_trace(k, j) = val;
      }
    }
  }

  // NGME_SPACETIME_HK_CHECK=1 recomputes the whole H_K block the generic way --
  // Hutchinson probes against a factorization of the FULL nt*ns operator -- and
  // reports the worst disagreement. Those are stochastic, so the tolerance is
  // the probe noise, not machine precision; what this catches is a structurally
  // wrong identity, which would be off by a factor, not by a few percent.
  if (opts.compute_HK_trace && std::getenv("NGME_SPACETIME_HK_CHECK")) {
    Eigen::SparseLU<SparseMatrix<double>> fullK;
    fullK.compute(K);
    if (fullK.info() == Eigen::Success) {
      const int m = (int)K.rows();
      double worst = 0.0;
      int wj = -1, wk = -1;
      for (int j = 0; j < n_theta_K; ++j) {
        for (int k = j; k < n_theta_K; ++k) {
          // tr(K^-1 dK_k K^-1 dK_j) densely, as the reference.
          MatrixXd Xk = fullK.solve(MatrixXd(dK[k]));
          MatrixXd Xj = fullK.solve(MatrixXd(dK[j]));
          double ref = -(Xk.array() * Xj.transpose().array()).sum();
          const SparseMatrix<double, 0, int> &d2 = get_d2K(j, k);
          if (d2.rows() == m)
            ref += MatrixXd(fullK.solve(MatrixXd(d2))).trace();
          const double rel =
              std::abs(HK_trace(j, k) - ref) / std::max(1.0, std::abs(ref));
          if (rel > worst) { worst = rel; wj = j; wk = k; }
        }
      }
      std::fprintf(stderr,
                   "[spacetime HK check] max rel diff %.3e at (%d,%d)\n",
                   worst, wj, wk);
    }
  }

  // Verify once, by checking the STRUCTURE the identity rests on rather than
  // recomputing the identity itself. tr(K^-1 dK) = sum_t tr(K_tt^-1 (dK)_tt)
  // holds because K is block lower-triangular; and the shortcut above, which
  // scales one interior block by nt-1, additionally needs the interior diagonal
  // blocks to be equal. Both are O(nnz) to check.
  if (!block_trace_checked_) {
    block_trace_checked_ = true;
    // (a) nothing above the block diagonal
    double upper = 0.0;
    for (int c2 = 0; c2 < K.outerSize(); ++c2)
      for (SparseMatrix<double>::InnerIterator it(K, c2); it; ++it) {
        const int bi = (int)it.row() / n, bj = (int)it.col() / n;
        if (bj > bi)
          upper = std::max(upper, std::abs(it.value()));
      }
    if (upper > 1e-10 * std::max(1.0, K.coeffs().abs().maxCoeff()))
      throw std::runtime_error(
          "spacetime: K is not block lower-triangular (largest entry above the "
          "block diagonal is " + std::to_string(upper) +
          "), so the diagonal-block trace identity does not hold.");
    // (b) if one interior block is being scaled by nt-1, they must be equal
    if (uniform_interior && nt > 2) {
      const SparseMatrix<double> B1 = get_block(K, n, n);
      for (int t = 2; t < nt; ++t) {
        SparseMatrix<double> Bt = get_block(K, t * n, t * n);
        SparseMatrix<double> D = Bt - B1;
        if (D.norm() > 1e-8 * std::max(1.0, B1.norm()))
          throw std::runtime_error(
              "spacetime: interior diagonal blocks differ (block " +
              std::to_string(t) +
              "), so scaling one of them by nt-1 is not valid.");
      }
    }
  }

  return true;
}
