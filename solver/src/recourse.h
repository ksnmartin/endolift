// Recourse dual LP (one warm-started HiGHS object, costs change per solve) and the per-scenario cut pool.
#pragma once
#include <cstdint>
#include <memory>
#include <string>
#include <unordered_set>
#include <vector>
#include "Highs.h"
#include "instance.h"

class RecourseDual {
 public:
  explicit RecourseDual(const Instance& inst);
  // max alpha^T r - beta^T c over the dual polyhedron; returns the value, fills alpha (m) and beta (n) with a vertex.
  double solve(const double* r, const double* c, std::vector<double>& alpha, std::vector<double>& beta);
  // Q_s(y): r = d_s, c = cap * y.
  double Q(int s, const double* y, std::vector<double>& alpha, std::vector<double>& beta);
  long solves = 0;

 private:
  const Instance& inst_;
  Highs h_;
  std::vector<double> cost_, c_;
  std::vector<HighsInt> idx_;
};

// A set of duals, one per worker thread (HiGHS objects are not thread-safe).
struct DualSet {
  DualSet(const Instance& inst, int workers);
  std::vector<std::unique_ptr<RecourseDual>> duals;
  int size() const { return (int)duals.size(); }
  long solves() const;
};

// Run fn(s, dual) for s = 0..S-1 across the duals in parallel (fn must only touch per-s state).
template <class Fn>
void for_each_scenario(int S, DualSet& ds, Fn fn);

struct DualVertex { std::vector<double> alpha, beta; };

class CutPool {
 public:
  explicit CutPool(int S, double tol = 1e-9) : pools_(S), keys_(S), tol_(tol) {}
  bool add(int s, const std::vector<double>& alpha, const std::vector<double>& beta);   // false if a duplicate
  const std::vector<DualVertex>& cuts(int s) const { return pools_[s]; }
  size_t size(int s) const { return pools_[s].size(); }
  size_t total() const;
  int S() const { return (int)pools_.size(); }

 private:
  std::vector<std::vector<DualVertex>> pools_;
  std::vector<std::unordered_set<std::string>> keys_;
  double tol_;
};

// f^T y + sum_s p_s(y) Q_s(y) by direct evaluation (S recourse LPs, parallel over the duals).
double F_exact(const Instance& inst, DualSet& ds, const double* y);

#include <atomic>
#include <thread>
template <class Fn>
void for_each_scenario(int S, DualSet& ds, Fn fn) {
  int T = std::min(ds.size(), S);
  if (T <= 1) { for (int s = 0; s < S; ++s) fn(s, *ds.duals[0]); return; }
  std::atomic<int> next(0);
  std::vector<std::thread> th;
  for (int t = 0; t < T; ++t)
    th.emplace_back([&, t] { for (int s = next++; s < S; s = next++) fn(s, *ds.duals[t]); });
  for (auto& x : th) x.join();
}
