#include "../operator.h"

OU::OU(const Rcpp::List &operator_list) : Operator(operator_list) {
  // Extract mesh and compute time differences
  if (operator_list.containsElementNamed("mesh") &&
      !Rf_isNull(operator_list["mesh"])) {
    Rcpp::List mesh = operator_list["mesh"];
    if (mesh.containsElementNamed("loc") && !Rf_isNull(mesh["loc"])) {
      Rcpp::NumericVector loc = mesh["loc"];
      n = loc.size();

      // Compute time differences dt = diff(loc)
      dt.resize(n - 1);
      for (int i = 0; i < n - 1; ++i) {
        dt(i) = loc[i + 1] - loc[i];
        if (dt(i) <= 0) {
          Rcpp::stop("OU model requires strictly increasing mesh locations");
        }
      }
    } else {
      Rcpp::stop("OU model requires mesh with 'loc' field");
    }
  } else {
    Rcpp::stop("OU model requires a mesh");
  }

  // Initialize K and Z (Z remains identity by default in base class)
  K.resize(n, n);
  Z.resize(n, n);
  Z.setIdentity();

  // Build initial K from theta_K
  Eigen::VectorXd theta = Rcpp::as<Eigen::VectorXd>(operator_list["theta_K"]);
  build_KZ(theta);
}

void OU::build_KZ(const VectorXd &theta_K) {
  // theta_K[0] is log(theta), where theta is the mean-reversion rate
  // rho[i] = exp(-theta * dt[i])
  double theta = std::exp(theta_K(0));

  // Compute rho values for each time interval
  VectorXd rho(n - 1);
  for (int i = 0; i < n - 1; ++i) {
    rho(i) = std::exp(-theta * dt(i));
    // Validate rho is in (0, 1)
    if (rho(i) <= 0 || rho(i) >= 1) {
      Rcpp::warning("OU: rho value %f is outside (0,1) for dt=%f, theta=%f",
                    rho(i), dt(i), theta);
    }
  }

  // Build band sparse matrix K using triplet form
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(n + (n - 1)); // main diagonal + sub-diagonal

  // Main diagonal
  // First element: sqrt(max(0, 1 - rho[0]^2))
  double first_diag = std::sqrt(std::max(0.0, 1.0 - rho(0) * rho(0)));
  triplets.emplace_back(0, 0, first_diag);

  // Remaining diagonal elements: 1
  for (int i = 1; i < n; ++i) {
    triplets.emplace_back(i, i, 1.0);
  }

  // Sub-diagonal: -rho[i] at position (i+1, i)
  for (int i = 0; i < n - 1; ++i) {
    triplets.emplace_back(i + 1, i, -rho(i));
  }

  // Build sparse matrix
  K.setFromTriplets(triplets.begin(), triplets.end());

  // Z remains identity (already set in constructor)
}
