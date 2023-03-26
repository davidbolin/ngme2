#ifndef NGME_RANDEFF_H
#define NGME_RANDEFF_H

#include <Rcpp.h>
#include <vector>
#include <string>
#include <iostream>
#include <numeric>

#include "model.h"
#include "block.h"

using std::vector;
enum Rand_eff_type {normal, nig};
enum Strategy {all, ws};

// Structure for random effects
class Randeff {
  private:
    Rand_eff_type family;
    VectorXd parameter;
    int n_params;
  public:
    Randeff(const Rcpp::List& list_ngmes, unsigned long seed);
    const VectorXd& get_parameter() const {return parameter;}
    void set_parameter(VectorXd& parameter) {this->parameter = parameter;}
    VectorXd grad();
};

// Ngme = a list of model over replicates
class Ngme : public Model {
private:
  // for each replicate
  vector<std::unique_ptr<BlockModel>> blocks;
  int n_repl;
  int n_params;
  vector<int> num_each_repl;
  double sum_num_each_repl;

  int sampling_strategy;
  std::mt19937 gen;
  std::discrete_distribution<int> weighted_sampler;

  // random effects
  vector<int> sample_indices;
public:
  Ngme(const Rcpp::List& list_ngmes, unsigned long seed, int sampling_strategy) :
    n_repl(list_ngmes.size()),
    sampling_strategy(sampling_strategy),
    gen(seed)
  {
    for (int i=0; i < n_repl; i++) {
      Rcpp::List block_model = Rcpp::as<Rcpp::List> (list_ngmes[i]);
      blocks.push_back(std::make_unique<BlockModel>(block_model, seed));
      num_each_repl.push_back(blocks[i]->get_n_obs());
    }
    sum_num_each_repl = std::accumulate(num_each_repl.begin(), num_each_repl.end(), 0.0);
    n_params = blocks[0]->get_n_params();

    // Init the random number generator
    weighted_sampler = std::discrete_distribution<int>(num_each_repl.begin(), num_each_repl.end());
  }

  VectorXd get_parameter() const {
    return blocks[0]->get_parameter();
  }

  VectorXd get_stepsizes() const {
    return blocks[0]->get_stepsizes();
  }

  void set_parameter(const VectorXd& x) {
    for (int i=0; i < n_repl; i++) {
      blocks[i]->set_parameter(x);
    }
  }

  MatrixXd precond() const {
    MatrixXd precond = MatrixXd::Zero(n_params, n_params);
    for (int i=0; i < n_repl; i++) {
      precond += blocks[i]->precond();
    }
    return precond;
  }

  VectorXd grad() {
    VectorXd g = VectorXd::Zero(n_params);
    for (int i=0; i < n_repl; i++) {
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
    for (int i=0; i < n_repl; i++) {
      output.push_back(blocks[i]->output());
    }
    return output;
  }

  vector<vector<VectorXd>> get_VW() const {
    vector<vector<VectorXd>> output;
    for (int i=0; i < n_repl; i++) {
      output.push_back(blocks[i]->get_VW());
    }
    return output;
  }

  void set_prev_VW(const vector<vector<VectorXd>>& prev_VW) {
    for (int i=0; i < n_repl; i++) {
      blocks[i]->set_prev_VW(prev_VW[i]);
    }
  }
};

#endif