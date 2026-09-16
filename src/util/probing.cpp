#include "../include/probing.h"

#include <algorithm>
#include <numeric>

namespace ngme_probing {

Adjacency build_adjacency(const Eigen::SparseMatrix<double, 0, int> &A) {
  Adjacency adj;
  const int n = A.rows();
  if (n <= 0 || A.cols() != n)
    return adj; // a rectangular operator has no index graph to colour

  // Degree of each vertex first, so the edge arrays can be filled in one pass
  // with no per-vertex vector. Both orientations are counted: the pattern may
  // be stored one-sided (a symmetric matrix held as a triangle) or genuinely
  // non-symmetric, and the graph wants the edge either way. Duplicates from a
  // matrix that stores both halves are removed when the rows are sorted below.
  std::vector<int> deg(n, 0);
  for (int c = 0; c < A.outerSize(); ++c)
    for (Eigen::SparseMatrix<double, 0, int>::InnerIterator it(A, c); it; ++it) {
      const int r = it.row();
      if (r == c)
        continue; // self-loops carry no information for a colouring
      ++deg[c];
      ++deg[r];
    }

  adj.n = n;
  adj.outer.assign(n + 1, 0);
  for (int i = 0; i < n; ++i)
    adj.outer[i + 1] = adj.outer[i] + deg[i];
  adj.inner.assign(adj.outer[n], 0);

  std::vector<int> fill(adj.outer.begin(), adj.outer.begin() + n);
  for (int c = 0; c < A.outerSize(); ++c)
    for (Eigen::SparseMatrix<double, 0, int>::InnerIterator it(A, c); it; ++it) {
      const int r = it.row();
      if (r == c)
        continue;
      adj.inner[fill[c]++] = r;
      adj.inner[fill[r]++] = c;
    }

  // Sort and unique each adjacency list. The BFS below stamps what it visits,
  // so duplicates would not corrupt the result, but they would be walked once
  // per copy at every level of every expansion.
  int out = 0;
  std::vector<int> new_outer(n + 1, 0);
  for (int i = 0; i < n; ++i) {
    const int b = adj.outer[i], e = adj.outer[i + 1];
    std::sort(adj.inner.begin() + b, adj.inner.begin() + e);
    const int keep =
        (int)(std::unique(adj.inner.begin() + b, adj.inner.begin() + e) -
              (adj.inner.begin() + b));
    for (int q = 0; q < keep; ++q)
      adj.inner[out + q] = adj.inner[b + q];
    out += keep;
    new_outer[i + 1] = out;
  }
  adj.inner.resize(out);
  adj.inner.shrink_to_fit();
  adj.outer.swap(new_outer);
  return adj;
}

std::vector<int> distance_colouring(const Adjacency &adj, int d,
                                    int max_colours, int &n_colours) {
  n_colours = 0;
  const int n = adj.n;
  if (n <= 0 || d < 1 || max_colours < 1)
    return {};

  std::vector<int> colour(n, -1);
  // stamp[w] holds the source vertex whose expansion last reached w, which
  // makes "already visited in THIS expansion" an integer compare instead of a
  // clear over the neighbourhood. Each vertex is a source exactly once, so the
  // stamps cannot collide.
  std::vector<int> stamp(n, -1);
  std::vector<int> frontier, next_frontier, seen;
  std::vector<char> used;

  // Colour the most constrained vertices first: greedy colouring is sensitive
  // to the order, and descending degree is what keeps the count near the
  // minimum for the mesh graphs this is used on.
  std::vector<int> order(n);
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(), order.end(), [&adj](int a, int b) {
    return (adj.outer[a + 1] - adj.outer[a]) > (adj.outer[b + 1] - adj.outer[b]);
  });

  for (int oi = 0; oi < n; ++oi) {
    const int v = order[oi];
    seen.clear();
    frontier.clear();
    stamp[v] = v;
    frontier.push_back(v);
    seen.push_back(v);
    for (int lev = 0; lev < d && !frontier.empty(); ++lev) {
      next_frontier.clear();
      for (size_t t = 0; t < frontier.size(); ++t) {
        const int u = frontier[t];
        for (int q = adj.outer[u]; q < adj.outer[u + 1]; ++q) {
          const int w = adj.inner[q];
          if (stamp[w] == v)
            continue;
          stamp[w] = v;
          next_frontier.push_back(w);
          seen.push_back(w);
        }
      }
      frontier.swap(next_frontier);
    }

    used.assign((size_t)n_colours + 1, 0);
    for (size_t t = 0; t < seen.size(); ++t) {
      const int c = colour[seen[t]];
      if (c >= 0 && c < n_colours)
        used[c] = 1;
    }
    int c0 = 0;
    while (c0 < n_colours && used[c0])
      ++c0;
    colour[v] = c0;
    if (c0 >= n_colours)
      n_colours = c0 + 1;
    if (n_colours > max_colours) {
      // Past the budget there is nothing to gain by finishing: the caller
      // cannot afford one probe per colour.
      n_colours = 0;
      return {};
    }
  }
  return colour;
}

} // namespace ngme_probing
