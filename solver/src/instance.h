// Instance data: two-stage facility location with endogenous (affine, decision-dependent) scenario probabilities.
//
// Recourse LP for scenario s (flows x_ik, unmet demand u_k):
//   Q_s(y) = min  sum_ik C_ik x_ik + pen * sum_k u_k   s.t.  sum_i x_ik + u_k = d_sk,  sum_k x_ik <= cap_i y_i,  x,u >= 0
// Dual:   max alpha^T d_s - sum_i beta_i cap_i y_i   s.t.  alpha_k - beta_i <= C_ik,  alpha_k <= pen,  0 <= beta <= pen.
// Probabilities p_s(y) = b_s + a_s^T y with sum_s b_s = 1, sum_s a_s = 0, and p_s(y) >= 0 at every binary y
// (bpos_s := b_s + sum_j min(a_sj, 0) >= 0).
#pragma once
#include <string>
#include <vector>

struct Instance {
  int n = 0, m = 0, S = 0;
  double pen = 10.0;
  std::vector<double> cap, f;     // (n)
  std::vector<double> C;          // (n*m) row-major: C[i*m + k]
  std::vector<double> d;          // (S*m): d[s*m + k]
  std::vector<double> b;          // (S)
  std::vector<double> A;          // (S*n): A[s*n + j]
  std::vector<double> bpos;       // (S), derived
  std::string name, endo = "strong";

  void finalize();                                   // derived fields + validity checks
  double p(int s, const double* y) const;            // p_s(y)
  const double* ds(int s) const { return &d[(size_t)s * m]; }
  double a(int s, int j) const { return A[(size_t)s * n + j]; }
  std::string describe() const;
};

Instance load_instance(const std::string& path);
void save_instance(const Instance& inst, const std::string& path);

// Instance generator. Costs C ~ U(1,5), demands U(1.5,4.5) scaled along a linear ramp over the scenarios, capacities
// ~1.8 x mean total demand / n, fixed costs calibrated from the recourse saving of each site times U(0.6,1.6).
// strength: probability mass one site can move; density: fraction of sites that influence demand; fmult: multiplier
// on the calibrated fixed costs; endo: "weak" | "strong". bpos0 (strong only): the high-demand scenarios get base
// mass b_s = 0, so bpos_s = 0 for them and the linking rows of those scenarios are vacuous.
// rng: "pcg64" is numpy.random.default_rng(seed) re-implemented (nprng.h) and names the instance gen_*: the paper's table
// instances n{25..50}_s{1..5} are make_instance(n, 15, 20, s) with this stream; "mt" is a std::mt19937_64 stream,
// names cgen_*, used for the sensitivity, sparse, bpos0 and customer-count instances.
Instance make_instance(int n, int m, int S, int seed, double strength = 0.9, double density = 1.0, double fmult = 1.0,
                       const std::string& endo = "strong", bool bpos0 = false, const std::string& rng = "pcg64");
