#ifndef NGME_PROBING_H
#define NGME_PROBING_H

#include <Eigen/SparseCore>
#include <vector>

// Structured probe vectors for the Hutchinson trace estimators.
//
// The plain estimator draws dense Rademacher probes and averages u^T A u. Its
// error comes entirely from the off-diagonal terms A_ik, i != k, which survive
// because u_i u_k is +-1 rather than zero. Nothing about the DISTRIBUTION of
// the probes removes those terms: Gaussian, Rademacher and orthogonalised
// probes all leave the same sum of squared off-diagonals, differing only in a
// constant that is 1 + O(m/n).
//
// What does remove them is the SUPPORT. Partition the indices into classes and
// let each probe carry signs on one class only: the surviving terms are then
// just the pairs lying in the SAME class. Colour the graph of the matrix so
// that two indices share a colour only when they are more than d edges apart,
// and every surviving pair is one whose A_ik has already decayed -- which for
// the inverse of a local operator means geometrically in d. The estimator stays
// exactly unbiased for any A, since the classes partition the index set and the
// signs are independent; only its variance changes.
//
// The number of colours does not grow with n -- it is set by the local
// connectivity of the graph, so a distance-2 colouring of a 2-d mesh needs the
// same ~19 colours at n = 4e3 as at n = 2e4. That is what makes this affordable:
// the probe count is unchanged, and the accuracy improves geometrically in the
// colouring distance instead of as 1/sqrt(probe count).
namespace ngme_probing {

// Undirected structural graph of a square matrix: the edges are the nonzero
// positions of A or of A^T. Only the pattern matters, so no values are kept.
struct Adjacency {
  int n{0};
  std::vector<int> outer; // size n + 1
  std::vector<int> inner; // size nnz
  bool empty() const { return n == 0; }
  void clear() {
    n = 0;
    outer.clear();
    inner.clear();
    outer.shrink_to_fit();
    inner.shrink_to_fit();
  }
};

// Build the graph above. Returns an empty Adjacency for a non-square matrix,
// which is the caller's signal that probing does not apply.
Adjacency build_adjacency(const Eigen::SparseMatrix<double, 0, int> &A);

// Greedy distance-d colouring: two indices carry the same colour only if the
// shortest path between them is longer than d edges.
//
// Returns the per-index colours in [0, n_colours) and writes the count to
// n_colours. Returns an empty vector (and n_colours = 0) if the colouring would
// need more than max_colours -- the caller has a probe budget, and a colouring
// larger than it is of no use, so the search is abandoned rather than run to
// completion on a graph that is too dense for this to pay.
std::vector<int> distance_colouring(const Adjacency &adj, int d,
                                    int max_colours, int &n_colours);

} // namespace ngme_probing

#endif
