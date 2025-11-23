#include "../operator.h"
#include "../include/transform_utils.h"
#include <Rcpp.h>
#include <cmath>
#include <string>
#include <vector>
#include <map>

generic_ns::generic_ns(const Rcpp::List& operator_list) :
    Operator(operator_list)
{
    // Extract matrices
    if (operator_list.containsElementNamed("matrices")) {
        matrices = Rcpp::as<std::vector<SparseMatrix<double, 0, int>>>(operator_list["matrices"]);
        
        // Check if matrices list is empty
        if (matrices.empty()) {
            Rcpp::stop("Matrices list cannot be empty");
        }
    } else {
        Rcpp::stop("Matrices list is required");
    }
    
    // Extract parameter map
    if (operator_list.containsElementNamed("param_map")) {
        Rcpp::List param_map_list = operator_list["param_map"];
        Rcpp::CharacterVector param_names = param_map_list.names();
        
        for (int i = 0; i < param_names.size(); i++) {
            std::string param_name = Rcpp::as<std::string>(param_names[i]);
            std::vector<int> indices = Rcpp::as<std::vector<int>>(param_map_list[param_name]);
            
            // Adjust indices to be 0-based
            for (int& idx : indices) {
                idx -= 1;
            }
            
            param_map[param_name] = indices;
        }
    }
    
    // Extract basis matrices and create param matrices
    // If B_theta_K is not supplied, initialize with matrices of 1s
    if (operator_list.containsElementNamed("B_theta_K") && !Rf_isNull(operator_list["B_theta_K"])) {
        Rcpp::List B_theta_K_list = operator_list["B_theta_K"];
        Rcpp::CharacterVector B_names = B_theta_K_list.names();
        
        // Store basis matrices
        for (int i = 0; i < B_names.size(); i++) {
            std::string param_name = Rcpp::as<std::string>(B_names[i]);
            MatrixXd B_mat = Rcpp::as<MatrixXd>(B_theta_K_list[param_name]);
            B_theta_K[param_name] = B_mat;
        }
    }
    
    // Extract position combinations - this is required
    if (operator_list.containsElementNamed("position") && !Rf_isNull(operator_list["position"])) {
        Rcpp::List position_list = operator_list["position"];
        for (int i = 0; i < position_list.size(); i++) {
            std::vector<int> pos = Rcpp::as<std::vector<int>>(position_list[i]);
            
            // Adjust indices to be 0-based
            for (int& idx : pos) {
                idx -= 1;
            }
            
            position.push_back(pos);
        }
    } else {
        Rcpp::stop("Position parameter is required for generic_ns");
    }
    
    // Extract parameter transformations
    if (operator_list.containsElementNamed("trans") && !Rf_isNull(operator_list["trans"])) {
        Rcpp::List trans_list = operator_list["trans"];
        Rcpp::CharacterVector trans_names = trans_list.names();
        
        for (int i = 0; i < trans_names.size(); i++) {
            std::string param_name = Rcpp::as<std::string>(trans_names[i]);
            std::string trans_type = Rcpp::as<std::string>(trans_list[param_name]);
            trans_map[param_name] = trans_type;
        }
    }
}

void generic_ns::build_KZ(const VectorXd& theta_K) {
    int n = matrices[0].rows();
    
    // Create parameter matrices using basis expansions
    std::map<std::string, SparseMatrix<double, 0, int>> param_matrices;
    std::vector<SparseMatrix<double, 0, int>> all_matrices;
    
    // First, create parameter matrices
    for (const auto& param_entry : param_map) {
        std::string param_name = param_entry.first;
        std::vector<int> indices = param_entry.second;
        
        // Get basis matrix
        MatrixXd B;
        if (B_theta_K.find(param_name) != B_theta_K.end()) {
            B = B_theta_K[param_name];
        } else {
            // Default to matrix of ones (1 column)
            B = MatrixXd::Ones(n, 1);
        }
        
        // Extract parameter values from theta_K
        VectorXd param_values(indices.size());
        for (size_t i = 0; i < indices.size(); i++) {
            int idx = indices[i];
            if (idx >= 0 && idx < theta_K.size()) {
                param_values(i) = theta_K(idx);
            }
        }
        
        // Apply basis expansion
        VectorXd expanded_values = B * param_values;
        
        // Apply transformation
        std::string trans_type = "identity";  // Default
        if (trans_map.find(param_name) != trans_map.end()) {
            trans_type = trans_map[param_name];
        }
        
        // Transform values
        for (int i = 0; i < expanded_values.size(); i++) {
            expanded_values(i) = ngme::transforms::apply_transform(expanded_values(i), trans_type);
        }
        
        // Create diagonal matrix from expanded values
        SparseMatrix<double, 0, int> diag_mat(n, n);
        diag_mat.reserve(n);
        
        for (int i = 0; i < n; i++) {
            diag_mat.insert(i, i) = expanded_values(i);
        }
        
        // Store parameter matrix
        param_matrices[param_name] = diag_mat;
    }
    
    // Create combined list of matrices (parameters first, then fixed matrices)
    for (const auto& param_pair : param_matrices) {
        all_matrices.push_back(param_pair.second);
    }
    
    for (const auto& mat : matrices) {
        all_matrices.push_back(mat);
    }
    
    // Initialize K as zero matrix
    K.resize(n, n);
    K.setZero();
    
    // Apply matrix combinations according to position
    for (const auto& pos : position) {
        // Start with identity matrix
        SparseMatrix<double, 0, int> term(n, n);
        term.setIdentity();
        
        // Multiply matrices in the combination
        for (int idx : pos) {
            if (idx >= 0 && idx < all_matrices.size()) {
                term = term * all_matrices[idx];
            } else {
                Rcpp::stop("Position index out of bounds: " + std::to_string(idx+1) + 
                           " > " + std::to_string(all_matrices.size()));
            }
        }
        
        // Add to K
        K += term;
    }
}

double generic_ns::apply_transform(double value, const std::string& trans_type) const {
    return ngme::transforms::apply_transform(value, trans_type);
} 
