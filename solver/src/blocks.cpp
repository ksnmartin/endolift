#include "blocks.h"
#include <algorithm>
#include <cmath>

ScenarioBlocks::ScenarioBlocks(const Instance& inst, const CutPool& pool, std::vector<int> infl_) : infl(std::move(infl_)), inst_(inst), ad_(inst.S), G_(inst.S), count_(inst.S, 0) {
  if (infl.empty()) for (int j = 0; j < inst.n; ++j) infl.push_back(j);
  add_cuts(pool);
}

void ScenarioBlocks::append(int s, const DualVertex& v) {
  double ad = 0.0;
  const double* ds = inst_.ds(s);
  for (int k = 0; k < inst_.m; ++k) ad += v.alpha[k] * ds[k];
  ad_[s].push_back(ad);
  for (int i = 0; i < inst_.n; ++i) G_[s].push_back(v.beta[i] * inst_.cap[i]);
}

int ScenarioBlocks::add_cuts(const CutPool& pool) {
  int added = 0;
  for (int s = 0; s < inst_.S; ++s) {
    const auto& cs = pool.cuts(s);
    for (size_t k = count_[s]; k < cs.size(); ++k) { append(s, cs[k]); ++added; }
    count_[s] = (int)cs.size();
  }
  return added;
}

// GY[k*n + j] = sum_l G[k,l] Y[j,l]
static void gy_gY(const std::vector<double>& G, int K, int n, const std::vector<int>& infl, const double* y, const double* Y, std::vector<double>& Gy, std::vector<double>& GY) {
  Gy.assign(K, 0.0); GY.assign((size_t)K * n, 0.0);
  for (int k = 0; k < K; ++k) {
    const double* g = &G[(size_t)k * n];
    double acc = 0.0;
    for (int l = 0; l < n; ++l) acc += g[l] * y[l];
    Gy[k] = acc;
    double* out = &GY[(size_t)k * n];
    for (int j : infl) {
      const double* yr = Y + (size_t)j * n;
      double a = 0.0;
      for (int l = 0; l < n; ++l) a += g[l] * yr[l];
      out[j] = a;
    }
  }
}

void ScenarioBlocks::levels(int s, const double* y, const double* Y, double& e0, std::vector<double>& Lp, std::vector<double>& Lm) const {
  const int n = inst_.n, K = (int)ad_[s].size();
  Lp.assign(n, 0.0); Lm.assign(n, 0.0); e0 = 0.0;
  if (K == 0) return;
  const auto& ad = ad_[s];
  std::vector<double> Gy, GY;
  gy_gY(G_[s], K, n, infl, y, Y, Gy, GY);
  double L0 = 0.0;
  for (int k = 0; k < K; ++k) L0 = std::max(L0, ad[k] - Gy[k]);
  for (int j : infl) {
    double lp = -INFINITY, lm = -INFINITY;
    for (int k = 0; k < K; ++k) {
      double gyj = GY[(size_t)k * n + j];
      lp = std::max(lp, ad[k] * y[j] - gyj);
      lm = std::max(lm, ad[k] * (1.0 - y[j]) - Gy[k] + gyj);
    }
    Lp[j] = std::max(0.0, lp); Lm[j] = std::max(0.0, lm);
  }
  e0 = L0;
  if (!unlinked) for (int j : infl) e0 = std::max(e0, Lp[j] + Lm[j]);
}

double ScenarioBlocks::value_grad(int s, const double* y, const double* Y, std::vector<double>& gy, std::vector<double>& gY) const {
  const int n = inst_.n, K = (int)ad_[s].size();
  gy.assign(n, 0.0); gY.assign((size_t)n * n, 0.0);
  if (K == 0) return 0.0;
  const auto& ad = ad_[s]; const auto& G = G_[s];
  std::vector<double> Gy, GY;
  gy_gY(G, K, n, infl, y, Y, Gy, GY);
  int k0 = 0; double c0best = -INFINITY;
  for (int k = 0; k < K; ++k) { double c = ad[k] - Gy[k]; if (c > c0best) { c0best = c; k0 = k; } }
  const double L0 = std::max(0.0, c0best);
  std::vector<int> kp(n, 0), km(n, 0);
  std::vector<double> Lp(n), Lm(n);
  for (int j : infl) {
    double bp = -INFINITY, bm = -INFINITY;
    for (int k = 0; k < K; ++k) {
      double gyj = GY[(size_t)k * n + j];
      double cp = ad[k] * y[j] - gyj, cm = ad[k] * (1.0 - y[j]) - Gy[k] + gyj;
      if (cp > bp) { bp = cp; kp[j] = k; }
      if (cm > bm) { bm = cm; km[j] = k; }
    }
    Lp[j] = std::max(0.0, bp); Lm[j] = std::max(0.0, bm);
  }
  auto grad_p = [&](int j, double w) {
    if (Lp[j] > 0.0) {
      int k = kp[j]; gy[j] += w * ad[k];
      const double* g = &G[(size_t)k * n]; double* row = &gY[(size_t)j * n];
      for (int l = 0; l < n; ++l) row[l] -= w * g[l];
    }
  };
  auto grad_m = [&](int j, double w) {
    if (Lm[j] > 0.0) {
      int k = km[j]; gy[j] -= w * ad[k];
      const double* g = &G[(size_t)k * n]; double* row = &gY[(size_t)j * n];
      for (int l = 0; l < n; ++l) { gy[l] -= w * g[l]; row[l] += w * g[l]; }
    }
  };
  int jbest = infl[0]; double sbest = -INFINITY;
  for (int j : infl) { double sm = Lp[j] + Lm[j]; if (sm > sbest) { sbest = sm; jbest = j; } }
  const double bpos = inst_.bpos[s];
  double e0;
  if (!unlinked && sbest > L0) {
    e0 = sbest; grad_p(jbest, bpos); grad_m(jbest, bpos);
  } else {
    e0 = L0;
    if (L0 > 0.0) { const double* g = &G[(size_t)k0 * n]; for (int l = 0; l < n; ++l) gy[l] -= bpos * g[l]; }
  }
  double val = bpos * e0;
  for (int j : infl) {
    double a = inst_.a(s, j);
    if (a > 0.0) { val += a * Lp[j]; grad_p(j, a); }
    else if (a < 0.0) { val += -a * Lm[j]; grad_m(j, -a); }
  }
  return val;
}
