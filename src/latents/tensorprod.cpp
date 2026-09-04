// getK from matern class
// getK from ar class
// newK = K1 (time) %x% k2 (space)
// dKtx_t = dKt_t %x% Kx
// dKtx_x = kt %x% dKx_x

#include "../operator.h"
#include <algorithm>
#include "MatrixAlgebra.h"
#include <unsupported/Eigen/KroneckerProduct>

// tensor product for the C G class
Tensor_prod::Tensor_prod(const Rcpp::List &operator_list)
    : Operator(operator_list),
      first(OperatorFactory::create(operator_list["first"])),
      second(OperatorFactory::create(operator_list["second"])),
      n_theta_1(first->get_n_theta_K()), n_theta_2(second->get_n_theta_K()) {}

void Tensor_prod::build_KZ(const VectorXd &theta_K) {
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
  sub.compute_d2K = opts.compute_HK_trace;
  sub.compute_d2Z = false;
  sub.compute_trace = true;

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

  // Also restores each factor's K to base theta: the numeric differencing in
  // Operator::update_all rebuilds the composite from perturbed factors and
  // restores only the composite, leaving the factors at the last perturbed
  // value.
  first->update_all(theta.segment(0, n_theta_1), o1);
  second->update_all(theta.segment(n_theta_1, n_theta_2), o2);

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
bool Spacetime::compute_traces_structured(const VectorXd &theta,
                                          const UpdateOptions &opts) {
  if (!stationary_init || !opts.compute_trace)
    return false;
  // The Hessian traces are not implemented here; let the generic path have them
  // rather than returning a partial answer.
  if (opts.compute_HK_trace)
    return false;
  if (n_theta_K <= 0 || (int)dK.size() < n_theta_K || ns_ <= 0 || nt <= 0)
    return false;
  if (K.rows() != (long long)ns_ * nt || K.cols() != K.rows())
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

  std::vector<Eigen::SparseLU<SparseMatrix<double>>> lu(n_distinct);
  std::vector<int> rep(n_distinct);
  for (int b = 0; b < n_distinct; ++b) {
    rep[b] = (b == 0) ? 0 : b; // block index this factorization represents
    SparseMatrix<double> A = get_block(K, rep[b] * n, rep[b] * n);
    lu[b].compute(A);
    if (lu[b].info() != Eigen::Success)
      return false; // fall back to the generic path rather than guess
  }

  if (trace_vals.size() != n_theta_K)
    trace_vals = VectorXd::Zero(n_theta_K);
  trace_vals.setZero();

  for (int j = 0; j < n_theta_K; ++j) {
    if (!opts.fix_mask_thetaK.empty() && opts.fix_mask_thetaK[j])
      continue;
    double tr = 0.0;
    if (uniform_interior) {
      // block 0 once, then one interior block scaled by how many there are
      MatrixXd X0 = lu[0].solve(MatrixXd(get_block(dK[j], 0, 0)));
      tr += X0.trace();
      if (nt > 1) {
        MatrixXd X1 = lu[std::min(1, n_distinct - 1)].solve(
            MatrixXd(get_block(dK[j], n, n)));
        tr += (double)(nt - 1) * X1.trace();
      }
    } else {
      for (int t = 0; t < nt; ++t) {
        MatrixXd Xt = lu[t].solve(MatrixXd(get_block(dK[j], t * n, t * n)));
        tr += Xt.trace();
      }
    }
    trace_vals(j) = tr;
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
