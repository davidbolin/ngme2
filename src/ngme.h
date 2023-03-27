#ifndef NGME_H
#define NGME_H

#include <Rcpp.h>
#include <vector>
#include <string>
#include <iostream>
#include <numeric>

#include "include/MatrixAlgebra.h"
#include "model.h"
#include "block.h"

using std::vector;

enum Strategy {all, ws};

// Ngme = a list of model over replicates
class Ngme : public Model {
private:
  // for each replicate
  vector<std::unique_ptr<BlockModel>> ngme_repls;
  int n_repl;
  int n_params;
  vector<int> num_each_repl;
  double sum_num_each_repl;

  int sampling_strategy;
  std::mt19937 gen;
  std::discrete_distribution<int> weighted_sampler;

public:
  Ngme(const Rcpp::List& list_ngmes, unsigned long seed, int sampling_strategy);

  VectorXd get_parameter() const {
    return ngme_repls[0]->get_parameter();
  }

  VectorXd get_stepsizes() const {
    return ngme_repls[0]->get_stepsizes();
  }

  void set_parameter(const VectorXd& x) {
    for (int i=0; i < n_repl; i++) {
      ngme_repls[i]->set_parameter(x);
    }
  }

  MatrixXd precond() const {
    MatrixXd precond = MatrixXd::Zero(n_params, n_params);
    for (int i=0; i < n_repl; i++) {
      precond += ngme_repls[i]->precond();
    }
    return precond;
  }

  VectorXd grad() {
    VectorXd g = VectorXd::Zero(n_params);
    for (int i=0; i < n_repl; i++) {
      g += ngme_repls[i]->grad();
    }
    return g;
  }

  // add subsampling for some ngme_repls (SGD-IS sample only 1 replicate each time, according to the weights)
  VectorXd precond_grad();

  std::string get_par_string() const {
    return ngme_repls[0]->get_par_string();
  }

  int get_n_params() const {
    return n_params;
  }

  vector<Rcpp::List> output() const {
    vector<Rcpp::List> output;
    for (int i=0; i < n_repl; i++) {
      output.push_back(ngme_repls[i]->output());
    }
    return output;
  }

  vector<vector<VectorXd>> get_VW() const {
    vector<vector<VectorXd>> output;
    for (int i=0; i < n_repl; i++) {
      output.push_back(ngme_repls[i]->get_VW());
    }
    return output;
  }

  void set_prev_VW(const vector<vector<VectorXd>>& prev_VW) {
    for (int i=0; i < n_repl; i++) {
      ngme_repls[i]->set_prev_VW(prev_VW[i]);
    }
  }
};

#endif