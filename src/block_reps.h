#ifndef BLOCK_REPS_H
#define BLOCK_REPS_H

#include <Rcpp.h>
#include <vector>
#include <string>
#include <iostream>
#include <numeric>

#include "model.h"
#include "block.h"

using std::vector;

// replicates of block models
class Block_reps : public Model {
private:
  vector<std::unique_ptr<BlockModel>> blocks;
  int n_blocks;
  int n_params;
  vector<int> n_obs;
  double sum_n_obs;

  int sampling_strategy;
  std::mt19937 gen;
  std::discrete_distribution<int> importance_sampler;
public:
  Block_reps(const Rcpp::List& list_ngmes, unsigned long seed, int sampling_strategy) :
    n_blocks(list_ngmes.size()),
    sampling_strategy(sampling_strategy),
    gen(seed)
  {
    for (int i=0; i < n_blocks; i++) {
      Rcpp::List block_model = Rcpp::as<Rcpp::List> (list_ngmes[i]);
      blocks.push_back(std::make_unique<BlockModel>(block_model, seed));
      n_obs.push_back(blocks[i]->get_n_obs());
    }
    sum_n_obs = std::accumulate(n_obs.begin(), n_obs.end(), 0.0);
    n_params = blocks[0]->get_n_params();

    // Init the random number generator
    importance_sampler = std::discrete_distribution<int>(n_obs.begin(), n_obs.end());
  }

  VectorXd get_parameter() const {
    return blocks[0]->get_parameter();
  }

  VectorXd get_stepsizes() const {
    return blocks[0]->get_stepsizes();
  }

  void set_parameter(const VectorXd& x) {
    for (int i=0; i < n_blocks; i++) {
      blocks[i]->set_parameter(x);
    }
  }

  MatrixXd precond() const {
    MatrixXd precond = MatrixXd::Zero(n_params, n_params);
    for (int i=0; i < n_blocks; i++) {
      precond += blocks[i]->precond();
    }
    return precond;
  }

  VectorXd grad() {
    VectorXd g = VectorXd::Zero(n_params);
    for (int i=0; i < n_blocks; i++) {
      g += blocks[i]->grad();
    }
    return g;
  }

  // add subsampling for some blocks (SGD-IS sample only 1 replicate each time, according to the weights)
  VectorXd precond_grad();

  std::string get_par_string() const {
    return blocks[0]->get_par_string();
  }

  int get_n_params() const {
    return n_params;
  }

  vector<Rcpp::List> output() const {
    vector<Rcpp::List> output;
    for (int i=0; i < n_blocks; i++) {
      output.push_back(blocks[i]->output());
    }
    return output;
  }

  vector<vector<VectorXd>> get_VW() const {
    vector<vector<VectorXd>> output;
    for (int i=0; i < n_blocks; i++) {
      output.push_back(blocks[i]->get_VW());
    }
    return output;
  }

  void set_prev_VW(const vector<vector<VectorXd>>& prev_VW) {
    for (int i=0; i < n_blocks; i++) {
      blocks[i]->set_prev_VW(prev_VW[i]);
    }
  }
};

#endif