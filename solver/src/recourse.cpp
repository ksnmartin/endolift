#include "recourse.h"
#include <cmath>
#include <cstring>
#include <stdexcept>

RecourseDual::RecourseDual(const Instance& inst) : inst_(inst) {
  const int n = inst.n, m = inst.m;
  h_.setOptionValue("output_flag", false);
  h_.setOptionValue("presolve", "off");
  // columns: alpha_k (free, <= pen), beta_i in [0, pen]; objective min -r^T alpha + c^T beta
  std::vector<double> lo(m + n), hi(m + n);
  for (int k = 0; k < m; ++k) { lo[k] = -kHighsInf; hi[k] = inst.pen; }
  for (int i = 0; i < n; ++i) { lo[m + i] = 0.0; hi[m + i] = inst.pen; }
  h_.addVars(m + n, lo.data(), hi.data());
  // rows: alpha_k - beta_i <= C_ik
  std::vector<double> rlo((size_t)n * m, -kHighsInf), rhi((size_t)n * m);
  std::vector<HighsInt> starts((size_t)n * m), idx(2 * (size_t)n * m);
  std::vector<double> val(2 * (size_t)n * m);
  HighsInt r = 0;
  for (int i = 0; i < n; ++i)
    for (int k = 0; k < m; ++k) {
      rhi[r] = inst.C[(size_t)i * m + k]; starts[r] = 2 * r;
      idx[2 * r] = k; val[2 * r] = 1.0; idx[2 * r + 1] = m + i; val[2 * r + 1] = -1.0; ++r;
    }
  h_.addRows(r, rlo.data(), rhi.data(), 2 * r, starts.data(), idx.data(), val.data());
  cost_.assign(m + n, 0.0); c_.resize(n);
  idx_.resize(m + n);
  for (int j = 0; j < m + n; ++j) idx_[j] = j;
}

double RecourseDual::solve(const double* r, const double* c, std::vector<double>& alpha, std::vector<double>& beta) {
  const int n = inst_.n, m = inst_.m;
  for (int k = 0; k < m; ++k) cost_[k] = -r[k];
  for (int i = 0; i < n; ++i) cost_[m + i] = c[i];
  h_.changeColsCost(m + n, idx_.data(), cost_.data());
  h_.run(); ++solves;
  if (h_.getModelStatus() != HighsModelStatus::kOptimal) {
    h_.clearSolver(); h_.run();
    if (h_.getModelStatus() != HighsModelStatus::kOptimal)
      throw std::runtime_error("recourse dual failed: " + h_.modelStatusToString(h_.getModelStatus()));
  }
  const auto& x = h_.getSolution().col_value;
  alpha.assign(x.begin(), x.begin() + m); beta.assign(x.begin() + m, x.begin() + m + n);
  return -h_.getInfo().objective_function_value;
}

double RecourseDual::Q(int s, const double* y, std::vector<double>& alpha, std::vector<double>& beta) {
  for (int i = 0; i < inst_.n; ++i) c_[i] = inst_.cap[i] * y[i];
  return solve(inst_.ds(s), c_.data(), alpha, beta);
}

DualSet::DualSet(const Instance& inst, int workers) {
  for (int t = 0; t < std::max(1, workers); ++t) duals.emplace_back(new RecourseDual(inst));
}
long DualSet::solves() const { long k = 0; for (auto& d : duals) k += d->solves; return k; }

bool CutPool::add(int s, const std::vector<double>& alpha, const std::vector<double>& beta) {
  std::string key((alpha.size() + beta.size()) * sizeof(int64_t), '\0');
  char* p = &key[0];
  for (double v : alpha) { int64_t q = (int64_t)std::llround(v / tol_); std::memcpy(p, &q, sizeof q); p += sizeof q; }
  for (double v : beta) { int64_t q = (int64_t)std::llround(v / tol_); std::memcpy(p, &q, sizeof q); p += sizeof q; }
  if (!keys_[s].insert(std::move(key)).second) return false;
  pools_[s].push_back({alpha, beta});
  return true;
}
size_t CutPool::total() const { size_t k = 0; for (auto& p : pools_) k += p.size(); return k; }

double F_exact(const Instance& inst, DualSet& ds, const double* y) {
  std::vector<double> Qv(inst.S);
  for_each_scenario(inst.S, ds, [&](int s, RecourseDual& dual) {
    std::vector<double> al, be; Qv[s] = dual.Q(s, y, al, be);
  });
  double v = 0.0;
  for (int i = 0; i < inst.n; ++i) v += inst.f[i] * y[i];
  for (int s = 0; s < inst.S; ++s) v += inst.p(s, y) * Qv[s];
  return v;
}
