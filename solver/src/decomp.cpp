#include "decomp.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <string>

static double now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }

static std::vector<int> influencing(const Instance& inst, bool sparse) {
  std::vector<int> I;
  for (int j = 0; j < inst.n; ++j) {
    bool any = !sparse;
    for (int s = 0; s < inst.S && !any; ++s) any = inst.a(s, j) != 0.0;
    if (any) I.push_back(j);
  }
  return I;
}

DecomposedBound::DecomposedBound(const Instance& inst, const CutPool& pool, double lam, int max_age, int edge_weight, bool reset_each, bool sparse, bool lazy_mc)
    : blocks(inst, pool, influencing(inst, sparse)), n(inst.n), S(inst.S), inst_(inst), lam_(lam), max_age_(max_age), reset_each_(reset_each), lazy_mc_(lazy_mc) {
  std::vector<char> infl(n, 0); for (int j : blocks.infl) infl[j] = 1;
  pidx.assign((size_t)n * n, -1);
  for (int j = 0; j < n; ++j) for (int l = j + 1; l < n; ++l)
    if (infl[j] || infl[l]) { pidx[(size_t)j * n + l] = (int)pj.size(); pj.push_back(j); pl.push_back(l); }
  P = (int)pj.size(); nvar = n + P + S;
  h_.setOptionValue("output_flag", false);
  h_.setOptionValue("simplex_dual_edge_weight_strategy", edge_weight);   // 1 devex (~20% faster than steepest edge here)
  std::vector<double> lo(nvar, 0.0), hi(nvar, kHighsInf);
  for (int j = 0; j < n + P; ++j) hi[j] = 1.0;
  h_.addVars(nvar, lo.data(), hi.data());
  cost_.assign(nvar, 0.0);
  for (int j = 0; j < n; ++j) cost_[j] = inst.f[j];
  for (int s = 0; s < S; ++s) cost_[n + P + s] = 1.0;
  std::vector<HighsInt> all(nvar); for (int j = 0; j < nvar; ++j) all[j] = j;
  h_.changeColsCost(nvar, all.data(), cost_.data());
  // McCormick rows per pair k: Yp_k - y_j <= 0, Yp_k - y_l <= 0, Yp_k - y_j - y_l >= -1 (all up front, or lazily)
  mc_have_.assign((size_t)3 * P, lazy_mc_ ? 0 : 1);
  if (!lazy_mc_) {
    std::vector<double> L, U, V; std::vector<HighsInt> st, I;
    for (int k = 0; k < P; ++k) {
      HighsInt c = n + k;
      st.push_back((HighsInt)I.size()); L.push_back(-kHighsInf); U.push_back(0.0); I.insert(I.end(), {c, pj[k]}); V.insert(V.end(), {1.0, -1.0});
      st.push_back((HighsInt)I.size()); L.push_back(-kHighsInf); U.push_back(0.0); I.insert(I.end(), {c, pl[k]}); V.insert(V.end(), {1.0, -1.0});
      st.push_back((HighsInt)I.size()); L.push_back(-1.0); U.push_back(kHighsInf); I.insert(I.end(), {c, pj[k], pl[k]}); V.insert(V.end(), {1.0, -1.0, -1.0});
    }
    add_rows(L, U, st, I, V, false);
  }
  lb_.assign(n, 0.0); ub_.assign(n, 1.0);
}

// lazy_mc: add the McCormick rows violated at the master point (y, Yp); returns the number added.
int DecomposedBound::add_violated_mc(const double* y, const double* Yp) {
  std::vector<double> L, U, V; std::vector<HighsInt> st, I;
  const double tol = 1e-7;
  for (int k = 0; k < P; ++k) {
    HighsInt c = n + k; int j = pj[k], l = pl[k];
    if (!mc_have_[3 * k] && Yp[k] > y[j] + tol) { mc_have_[3 * k] = 1; st.push_back((HighsInt)I.size()); L.push_back(-kHighsInf); U.push_back(0.0); I.insert(I.end(), {c, (HighsInt)j}); V.insert(V.end(), {1.0, -1.0}); }
    if (!mc_have_[3 * k + 1] && Yp[k] > y[l] + tol) { mc_have_[3 * k + 1] = 1; st.push_back((HighsInt)I.size()); L.push_back(-kHighsInf); U.push_back(0.0); I.insert(I.end(), {c, (HighsInt)l}); V.insert(V.end(), {1.0, -1.0}); }
    if (!mc_have_[3 * k + 2] && Yp[k] < y[j] + y[l] - 1.0 - tol) { mc_have_[3 * k + 2] = 1; st.push_back((HighsInt)I.size()); L.push_back(-1.0); U.push_back(kHighsInf); I.insert(I.end(), {c, (HighsInt)j, (HighsInt)l}); V.insert(V.end(), {1.0, -1.0, -1.0}); }
  }
  if (!L.empty()) add_rows(L, U, st, I, V, false);
  mc_rows_added += (long)L.size();
  return (int)L.size();
}

void DecomposedBound::add_rows(const std::vector<double>& L, const std::vector<double>& U, const std::vector<HighsInt>& starts,
                               const std::vector<HighsInt>& idx, const std::vector<double>& val, bool is_cut) {
  HighsStatus st = h_.addRows((HighsInt)L.size(), L.data(), U.data(), (HighsInt)idx.size(), starts.data(), idx.data(), val.data());
  if (st == HighsStatus::kError) throw std::runtime_error("addRows failed");
  row_lo_.insert(row_lo_.end(), L.begin(), L.end());
  row_cut_.insert(row_cut_.end(), L.size(), is_cut ? 1 : 0);
  row_age_.insert(row_age_.end(), L.size(), 0);
  if (keep_cuts) row_mc_.resize(row_lo_.size());
}

int DecomposedBound::age_rows() {
  if (max_age_ < 0) return 0;
  const auto& st = h_.getBasis().row_status;
  const auto& rv = h_.getSolution().row_value;
  const size_t R = row_lo_.size();
  std::vector<HighsInt> kill;
  for (size_t r = 0; r < R; ++r) {
    bool basic = st[r] == HighsBasisStatus::kBasic;
    bool slack = (rv[r] - row_lo_[r]) > 1e-7 * std::max(1.0, std::fabs(row_lo_[r]));
    row_age_[r] = (basic && slack && row_cut_[r]) ? row_age_[r] + 1 : 0;
    if (row_age_[r] > max_age_) kill.push_back((HighsInt)r);
  }
  if (kill.empty()) return 0;
  h_.deleteRows((HighsInt)kill.size(), kill.data());
  size_t w = 0, ki = 0;
  for (size_t r = 0; r < R; ++r) {
    if (ki < kill.size() && (size_t)kill[ki] == r) { ++ki; continue; }
    row_lo_[w] = row_lo_[r]; row_cut_[w] = row_cut_[r]; row_age_[w] = row_age_[r];
    if (keep_cuts && r < row_mc_.size()) row_mc_[w] = std::move(row_mc_[r]);
    ++w;
  }
  row_lo_.resize(w); row_cut_.resize(w); row_age_.resize(w);
  if (keep_cuts) row_mc_.resize(w);
  deleted += (long)kill.size();
  return (int)kill.size();
}

void DecomposedBound::set_y_bounds(const double* lb, const double* ub) {
  lb_.assign(lb, lb + n); ub_.assign(ub, ub + n);
  std::vector<HighsInt> idx(n); for (int j = 0; j < n; ++j) idx[j] = j;
  h_.changeColsBounds(n, idx.data(), lb_.data(), ub_.data());
}

void DecomposedBound::project(std::vector<double>& y, std::vector<double>& Yp) const {
  for (int j = 0; j < n; ++j) y[j] = std::min(std::max(y[j], lb_[j]), ub_[j]);
  for (int k = 0; k < P; ++k) {
    double lo = std::max(0.0, y[pj[k]] + y[pl[k]] - 1.0), hi = std::min(y[pj[k]], y[pl[k]]);
    Yp[k] = std::min(std::max(Yp[k], lo), hi);
  }
}

void DecomposedBound::unfold(const double* y, const double* Yp, std::vector<double>& Y) const {
  Y.assign((size_t)n * n, 0.0);
  for (int j = 0; j < n; ++j) Y[(size_t)j * n + j] = y[j];
  for (int k = 0; k < P; ++k) { Y[(size_t)pj[k] * n + pl[k]] = Yp[k]; Y[(size_t)pl[k] * n + pj[k]] = Yp[k]; }
  // pairs absent from the sparse lifting are never read by the blocks (their rows j are not carried)
}

double DecomposedBound::evaluate(const double* y, const double* Yp, std::vector<MasterCut>& cuts) {
  unfold(y, Yp, Y_);
  cuts.resize(S);
  double tot = 0.0;
  for (int s = 0; s < S; ++s) {
    double v = blocks.value_grad(s, y, Y_.data(), gy_, gY_);
    MasterCut& c = cuts[s]; c.s = s; c.gy.resize(n); c.gp.resize(P);
    double dot = 0.0;
    for (int j = 0; j < n; ++j) { c.gy[j] = gy_[j] + gY_[(size_t)j * n + j]; dot += c.gy[j] * y[j]; }
    for (int k = 0; k < P; ++k) { c.gp[k] = gY_[(size_t)pj[k] * n + pl[k]] + gY_[(size_t)pl[k] * n + pj[k]]; dot += c.gp[k] * Yp[k]; }
    c.c = v - dot; tot += v;
  }
  return tot;
}

long DecomposedBound::add_cuts(const std::vector<MasterCut>& cuts) {
  if (cuts.empty()) return 0;
  std::vector<double> L, U, V; std::vector<HighsInt> st, I;
  for (const auto& c : cuts) {
    st.push_back((HighsInt)I.size()); L.push_back(c.c); U.push_back(kHighsInf);
    I.push_back(n + P + c.s); V.push_back(1.0);
    for (int j = 0; j < n; ++j) if (std::fabs(c.gy[j]) > 1e-9) { I.push_back(j); V.push_back(-c.gy[j]); }
    for (int k = 0; k < P; ++k) if (std::fabs(c.gp[k]) > 1e-9) { I.push_back(n + k); V.push_back(-c.gp[k]); }
  }
  add_rows(L, U, st, I, V, true);
  if (keep_cuts) { size_t R0 = row_lo_.size() - cuts.size(); for (size_t i = 0; i < cuts.size(); ++i) row_mc_[R0 + i] = cuts[i]; }
  ncuts += (long)cuts.size();
  return (long)cuts.size();
}

void DecomposedBound::export_state(std::vector<MasterCut>& cuts, std::vector<double>& cy, std::vector<double>& cYp, bool& has_center) const {
  cuts.clear();
  for (size_t r = 0; r < row_lo_.size() && r < row_mc_.size(); ++r) if (row_cut_[r] && !row_mc_[r].gy.empty()) cuts.push_back(row_mc_[r]);
  cy = cy_; cYp = cYp_; has_center = has_center_;
}

void DecomposedBound::import_state(const std::vector<MasterCut>& cuts, const std::vector<double>& cy, const std::vector<double>& cYp, bool has_center) {
  add_cuts(cuts);
  if (has_center) { cy_ = cy; cYp_ = cYp; has_center_ = true; }
}

void DecomposedBound::dump_stalled(const char* tag) {
  std::string base = std::string("stall_") + tag + "_" + std::to_string(h_.getNumRow());
  h_.writeModel(base + ".mps"); h_.writeBasis(base + ".bas");
  std::fprintf(stderr, "  [lp] stalled LP written to %s.mps\n", base.c_str());
}

// 0: warm dual simplex capped at first_cap; 1: same basis, reset solver state; 2: interior point + crossover from cold.
HighsModelStatus DecomposedBound::run_lp_escalating(const char* tag) {
  const double first_cap = 60.0; double used = 0.0; HighsModelStatus st = HighsModelStatus::kNotset;
  for (int attempt = 0; attempt < 3; ++attempt) {
    if (attempt == 1) { HighsBasis b = h_.getBasis(); h_.clearSolver(); h_.setBasis(b); }
    else if (attempt == 2) { h_.clearSolver(); h_.setOptionValue("solver", "ipm"); h_.setOptionValue("run_crossover", "on"); }
    double rem = std::max(1.0, lp_time_limit - used);
    double cap = attempt == 0 ? std::min(rem, first_cap) : (attempt == 1 ? first_cap : std::max(rem, 2 * first_cap));
    h_.setOptionValue("time_limit", h_.getRunTime() + cap);     // HiGHS checks time_limit against its cumulative clock
    double t = now(); h_.run(); double dt = now() - t; used += dt;
    st = h_.getModelStatus();
    if (attempt == 2) h_.setOptionValue("solver", "simplex");
    if (st == HighsModelStatus::kOptimal || st == HighsModelStatus::kInfeasible) {
      if (attempt > 0) std::fprintf(stderr, "  [lp] %s: recovered at attempt %d (%s, rows %d, %.1fs)\n", tag, attempt, h_.modelStatusToString(st).c_str(), (int)h_.getNumRow(), dt);
      return st;
    }
    std::fprintf(stderr, "  [lp] %s: attempt %d status %s rows %d simplex it %d %.0fs\n", tag, attempt, h_.modelStatusToString(st).c_str(), (int)h_.getNumRow(), (int)h_.getInfo().simplex_iteration_count, dt);
  }
  return st;
}

BoundResult DecomposedBound::run(const double* lb, const double* ub, double tol, double cutoff, int max_iter, bool verbose) {
  const double t0 = now();
  BoundResult R;
  if (lb || ub) set_y_bounds(lb ? lb : lb_.data(), ub ? ub : ub_.data());
  double LB = -INFINITY, UB = INFINITY;
  std::vector<MasterCut> cuts, viol;
  std::vector<double> y(n), Yp(P), th(S), ze_y(n), ze_Yp(P);
  bool have_center = false; std::vector<double> cy, cYp;
  auto fy = [&](const double* v) { double a = 0.0; for (int j = 0; j < n; ++j) a += inst_.f[j] * v[j]; return a; };
  auto violated = [&](const std::vector<MasterCut>& cs, std::vector<MasterCut>& out) {
    out.clear();
    for (const auto& c : cs) {
      double v = c.c;
      for (int j = 0; j < n; ++j) v += c.gy[j] * y[j];
      for (int k = 0; k < P; ++k) v += c.gp[k] * Yp[k];
      if (v > th[c.s] + 1e-9) out.push_back(c);
    }
  };
  if (has_center_) {
    cy = cy_; cYp = cYp_; project(cy, cYp);
    double t1 = now(); double tot = evaluate(cy.data(), cYp.data(), cuts); R.t_eval += now() - t1;
    UB = fy(cy.data()) + tot; have_center = true; add_cuts(cuts);
  }
  int it = 0;
  for (it = 1; it <= max_iter; ++it) {
    if (reset_each_) { HighsBasis b = h_.getBasis(); h_.clearSolver(); h_.setBasis(b); }
    double t1 = now(); HighsModelStatus st = run_lp_escalating("decomp-master"); R.t_master += now() - t1;
    if (st == HighsModelStatus::kInfeasible) { R.lb = R.ub = INFINITY; R.iters = it; R.infeasible = true; R.cuts = ncuts; return R; }
    if (st != HighsModelStatus::kOptimal) { dump_stalled("decomp-master"); throw LPTimeout("decomp master: " + h_.modelStatusToString(st)); }
    const auto& x = h_.getSolution().col_value;
    LB = h_.getInfo().objective_function_value;
    age_rows();
    std::copy(x.begin(), x.begin() + n, y.begin()); std::copy(x.begin() + n, x.begin() + n + P, Yp.begin()); std::copy(x.begin() + n + P, x.begin() + nvar, th.begin());
    if (lazy_mc_ && add_violated_mc(y.data(), Yp.data()) > 0) continue;     // master point outside the McCormick set: not a valid bound yet
    if (std::isfinite(UB) && UB - LB <= tol * std::max(1.0, std::fabs(UB))) break;
    if (std::isfinite(cutoff) && LB >= cutoff - tol * std::max(1.0, std::fabs(cutoff))) { R.pruned = true; break; }
    bool stabilized = have_center;
    if (stabilized) {
      for (int j = 0; j < n; ++j) ze_y[j] = lam_ * y[j] + (1 - lam_) * cy[j];
      for (int k = 0; k < P; ++k) ze_Yp[k] = lam_ * Yp[k] + (1 - lam_) * cYp[k];
    } else { ze_y = y; ze_Yp = Yp; }
    t1 = now(); double tot = evaluate(ze_y.data(), ze_Yp.data(), cuts); R.t_eval += now() - t1;
    double val = fy(ze_y.data()) + tot;
    if (val < UB) { UB = val; cy = ze_y; cYp = ze_Yp; have_center = true; }
    violated(cuts, viol);
    if (viol.empty() && stabilized) {       // stabilized point did not cut the master point: Kelley step
      t1 = now(); tot = evaluate(y.data(), Yp.data(), cuts); R.t_eval += now() - t1;
      val = fy(y.data()) + tot;
      if (val < UB) { UB = val; cy = y; cYp = Yp; }
      violated(cuts, viol);
    }
    if (verbose && (it % 10 == 0 || it == 1))
      std::printf("    it %4d LB %.6f UB %.6f gap %.4f%% rows %d (%.1fs)\n", it, LB, UB, (UB - LB) / std::max(1.0, std::fabs(UB)) * 100, (int)h_.getNumRow(), now() - t0);
    if ((std::isfinite(UB) && UB - LB <= tol * std::max(1.0, std::fabs(UB))) || viol.empty()) break;
    add_cuts(viol);
  }
  if (have_center) { cy_ = cy; cYp_ = cYp; has_center_ = true; }
  R.lb = LB; R.ub = UB; R.iters = std::min(it, max_iter); R.cuts = ncuts; R.rows = h_.getNumRow(); R.deleted = deleted; R.y = y; R.Yp = Yp;
  return R;
}
