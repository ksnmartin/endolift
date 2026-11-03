// Closed-form scenario blocks of the linked lifted-cut master (Section 4 of the paper):
//   L0 = max(0, max_k ad_k - g_k.y),  Lp_j = max(0, max_k ad_k y_j - g_k.Y_j),  Lm_j = max(0, max_k ad_k (1-y_j) - g_k.(y - Y_j))
//   e0* = max(L0, max_j Lp_j + Lm_j),  phi_s = bpos_s e0* + sum_{A_sj>0} A_sj Lp_j + sum_{A_sj<0} |A_sj| Lm_j.
#pragma once
#include <vector>
#include "instance.h"
#include "recourse.h"

class ScenarioBlocks {
 public:
  // infl: the sites whose factors (ep, em, linking rows) are carried; empty = all sites. With influence
  // sparsity (a_sj = 0 for every s) a site's factors have zero objective weight and are dropped: this is a valid
  // relaxation (the e0 cuts alone give exactness at binary points) that needs pair columns (j in infl) x (any l) only.
  ScenarioBlocks(const Instance& inst, const CutPool& pool, std::vector<int> infl = {});
  std::vector<int> infl;
  bool unlinked = false;   // drop the linking rows: e0* = L0 (the unlinked per-factor model C2R, an ablation)
  int add_cuts(const CutPool& pool);                 // append pool cuts not yet in the blocks; returns the number added
  // Block levels at (y, Y): e0*, Lp (n), Lm (n). Y is n x n row-major, symmetric with diag y.
  void levels(int s, const double* y, const double* Y, double& e0, std::vector<double>& Lp, std::vector<double>& Lm) const;
  // phi_s(y, Y) and a subgradient: gy (n), gY (n x n, w.r.t. the full symmetric lifted matrix). Both are overwritten.
  double value_grad(int s, const double* y, const double* Y, std::vector<double>& gy, std::vector<double>& gY) const;
  int count(int s) const { return count_[s]; }

 private:
  const Instance& inst_;
  std::vector<std::vector<double>> ad_;   // per s: (K)
  std::vector<std::vector<double>> G_;    // per s: (K x n) row-major, G_k = beta_k * cap
  std::vector<int> count_;
  void append(int s, const DualVertex& v);
  // scratch (per call; sized n or K*n) -- value_grad/levels are const and allocate their own scratch
};
