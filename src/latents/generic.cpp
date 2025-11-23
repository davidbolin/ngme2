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
