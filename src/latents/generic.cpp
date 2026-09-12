#include "../operator.h"
#include "../include/transform_utils.h"

Generic::Generic(const Rcpp::List& operator_list):
    Operator(operator_list),
    matrices(Rcpp::as<std::vector<SparseMatrix<double, 0, int>>>(operator_list["matrices"]))
{
    if (operator_list.containsElementNamed("trans") && !Rf_isNull(operator_list["trans"])) {
        Rcpp::List trans_list = operator_list["trans"];
        
        // Initialize the trans_map structure
        if (Rf_isNewList(trans_list)) {
            Rcpp::CharacterVector list_names = trans_list.names();
            
            // Store parameter names for later use
            for (int i = 0; i < list_names.size(); i++) {
                std::string param_name = Rcpp::as<std::string>(list_names[i]);
                param_names.push_back(param_name);
                
                // Get the transformations for this parameter
                Rcpp::CharacterVector param_trans = trans_list[param_name];
                std::vector<std::string> trans_vec;
                
                // Convert CharacterVector to std::vector<std::string>
                for (int j = 0; j < param_trans.size(); j++) {
                    trans_vec.push_back(Rcpp::as<std::string>(param_trans[j]));
                }
                
                // Store in trans_map
                trans_map[param_name] = trans_vec;
            }
        }
    }
}

void Generic::build_KZ(const VectorXd& theta_K) {
    // Default coefficients: all ones
    VectorXd coef = VectorXd::Ones(matrices.size());
    
    // If there are parameters, apply their transformations
    if (theta_K.size() > 0) {
        // For each parameter in theta_K
        for (int p = 0; p < param_names.size(); p++) {
            // Get parameter name and value
            std::string param_name = param_names[p];
            double param_value = theta_K[p];
            // Skip if parameter not in trans_map
            if (trans_map.find(param_name) == trans_map.end()) {
                continue;
            }
            const std::vector<std::string>& param_trans = trans_map[param_name];
            
            // Apply transformations for this parameter to each matrix
            for (size_t i = 0; i < param_trans.size() && i < matrices.size(); i++) {
                const std::string& trans_type = param_trans[i];
                if (trans_type != "null") {
                    double transformed_value = ngme::transforms::apply_transform(param_value, trans_type);
                    coef[i] *= transformed_value;
                }
            }
        }
    }
    // Build K as a linear combination of matrices with coefficients
    K.setZero();
    for (int i = 0; i < matrices.size(); i++) {
        K += coef[i] * matrices[i];
    }
}

double Generic::apply_transform(double value, const std::string& trans_type) const {
    return ngme::transforms::apply_transform(value, trans_type);
}

// ---------------------------------------------------------------------------
// Closed-form derivatives of the generic operator.
//
//   K = sum_i coef_i M_i,      coef_i = prod_p T(theta_p, trans[p][i])
//
// The matrices M_i are fixed data, so the only theta-dependence is in the
// scalar coefficients and
//
//   d coef_i / d theta_p        = T'(theta_p) prod_{q != p} T(theta_q)
//   d2 coef_i / dtheta_p dtheta_r = T'(theta_p) T'(theta_r) prod_{q != p,r} T
//   d2 coef_i / dtheta_p^2      = T''(theta_p) prod_{q != p} T(theta_q).
//
// The products are formed by explicit loops rather than by dividing coef_i
// through by one factor: a transform can legitimately return zero (sech does,
// at large |theta|), and dividing would turn that into a NaN.
// ---------------------------------------------------------------------------

void Generic::coef_factors(const VectorXd &theta_K, std::vector<VectorXd> &f,
                           std::vector<VectorXd> &f1,
                           std::vector<VectorXd> &f2) const {
  const int nm = (int)matrices.size();
  const int np = (int)param_names.size();
  f.assign(np, VectorXd::Ones(nm));
  f1.assign(np, VectorXd::Zero(nm));
  f2.assign(np, VectorXd::Zero(nm));
  for (int p = 0; p < np; ++p) {
    auto it = trans_map.find(param_names[p]);
    if (it == trans_map.end())
      continue; // absent: contributes a constant 1, derivative 0
    const std::vector<std::string> &tr = it->second;
    const double v = theta_K[p];
    for (int i = 0; i < nm && i < (int)tr.size(); ++i) {
      if (tr[i] == "null")
        continue; // constant 1
      f[p][i] = ngme::transforms::apply_transform(v, tr[i]);
      f1[p][i] = ngme::transforms::transform_derivative(v, tr[i], 1);
      f2[p][i] = ngme::transforms::transform_derivative(v, tr[i], 2);
    }
  }
}

bool Generic::update_dKdZ(const VectorXd &theta_K) {
  const int nm = (int)matrices.size();
  const int np = (int)param_names.size();
  if (np != n_theta_K || nm == 0 || (int)theta_K.size() != np)
    return false;

  std::vector<VectorXd> f, f1, f2;
  coef_factors(theta_K, f, f1, f2);

  if ((int)dK.size() != n_theta_K)
    dK.assign(n_theta_K, SparseMatrix<double>(K.rows(), K.cols()));
  if ((int)dZ.size() != n_theta_K)
    dZ.assign(n_theta_K, SparseMatrix<double>(K.rows(), K.cols()));

  for (int p = 0; p < np; ++p) {
    dK[p].setZero();
    for (int i = 0; i < nm; ++i) {
      double co = f1[p][i];
      for (int q = 0; q < np && co != 0.0; ++q)
        if (q != p)
          co *= f[q][i];
      if (co != 0.0)
        dK[p] += co * matrices[i];
    }
    // Z carries no parameters in this operator.
    dZ[p].setZero();
  }
  return true;
}

bool Generic::update_d2Kd2Z(const VectorXd &theta_K) {
  const int nm = (int)matrices.size();
  const int np = (int)param_names.size();
  if (np != n_theta_K || nm == 0 || (int)theta_K.size() != np)
    return false;

  std::vector<VectorXd> f, f1, f2;
  coef_factors(theta_K, f, f1, f2);

  if ((int)d2K.size() != n_theta_K)
    d2K.assign(n_theta_K,
               std::vector<SparseMatrix<double>>(
                   n_theta_K, SparseMatrix<double>(K.rows(), K.cols())));
  if ((int)d2Z.size() != n_theta_K)
    d2Z.assign(n_theta_K, std::vector<SparseMatrix<double, 0, int>>(
                              n_theta_K,
                              SparseMatrix<double, 0, int>(K.rows(), K.cols())));

  for (int p = 0; p < np; ++p)
    for (int r = p; r < np; ++r) {
      SparseMatrix<double> out(K.rows(), K.cols());
      for (int i = 0; i < nm; ++i) {
        double co = (p == r) ? f2[p][i] : f1[p][i] * f1[r][i];
        for (int q = 0; q < np && co != 0.0; ++q)
          if (q != p && q != r)
            co *= f[q][i];
        if (co != 0.0)
          out += co * matrices[i];
      }
      d2K[p][r] = out;
      d2Z[p][r].setZero();
      if (r != p) {
        d2K[r][p] = out;
        d2Z[r][p].setZero();
      }
    }
  return true;
}
