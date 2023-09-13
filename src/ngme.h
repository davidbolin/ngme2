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
// #include <chrono>

using std::vector;

enum Strategy {all, ws};

// Ngme = a list of repliacte models. Each replicate model is a BlockModel
// parameters(Ngme) = parameters(Ngme_replicate) (share the same parameter over replicate)
//    + parameters(Ngme_replicate_reff) (for each replicate)
// n_params = replicate.n_params + n_re_params * (n_repl - 1)
class Ngme : public Model {
private:
  // for each replicate
  vector<std::shared_ptr<BlockModel>> ngme_repls;
  int n_repl, n_params;

  vector<int> num_each_repl;
  double sum_num_each_repl;

  int sampling_strategy;
  std::mt19937 gen;
  std::discrete_distribution<int> weighted_sampler;

  bool debug;

  // for time profiling
  // std::chrono::milliseconds ngme_precond_time, ngme_grad_time, ngme_set_time, ngme_get_time;
public:
  Ngme(const Rcpp::List& list_ngmes, unsigned long seed, int sampling_strategy);

  VectorXd get_parameter() override;
  void set_parameter(const VectorXd& p) override;
  VectorXd get_stepsizes() override {
    return VectorXd::Ones(n_params);
  }

  MatrixXd precond(int strategy=0, double eps=1e-5) override;
  VectorXd grad() override;
  // add subsampling for some ngme_repls (SGD-IS sample only 1 replicate each time, according to the weights)

  std::string get_par_string() const {
    return ngme_repls[0]->get_par_string();
  }

  int get_n_params() const override {
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

  void burn_in(int iterations);
};

#endif