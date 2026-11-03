#include "mcb.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <queue>
#include <set>

static double now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }

namespace {

// The MC-B LP held incrementally in one HiGHS object: columns y (n) | theta (S) | z (one per (s, j) with a_sj != 0).
class McbMaster {
 public:
  McbMaster(const Instance& inst, const std::vector<double>& Qmin, const std::vector<double>& Qmax)
      : n(inst.n), S(inst.S), inst_(inst), Qmin_(Qmin), Qmax_(Qmax) {
    h_.setOptionValue("output_flag", false);
    zidx_.assign((size_t)S * n, -1);
    int nz = 0;
    for (int s = 0; s < S; ++s) for (int j = 0; j < n; ++j) if (inst.a(s, j) != 0.0) zidx_[(size_t)s * n + j] = n + S + nz++;
    nvar = n + S + nz;
    std::vector<double> lo(nvar, 0.0), hi(nvar), cost(nvar, 0.0);
    for (int j = 0; j < n; ++j) { hi[j] = 1.0; cost[j] = inst.f[j]; }
    for (int s = 0; s < S; ++s) { lo[n + s] = Qmin[s]; hi[n + s] = Qmax[s]; cost[n + s] = inst.b[s]; }
    for (int s = 0; s < S; ++s) for (int j = 0; j < n; ++j) { int c = zidx_[(size_t)s * n + j]; if (c >= 0) { hi[c] = Qmax[s]; cost[c] = inst.a(s, j); } }
    h_.addVars(nvar, lo.data(), hi.data());
    std::vector<HighsInt> all(nvar); for (int j = 0; j < nvar; ++j) all[j] = j;
    h_.changeColsCost(nvar, all.data(), cost.data());
    // four envelope rows per z column: z - Qmin y >= 0, z - Qmax y <= 0, z - th - Qmax y >= -Qmax, z - th - Qmin y <= -Qmin
    std::vector<double> L, U, V; std::vector<HighsInt> st, I;
    mcrow_.assign((size_t)S * n, -1);
    HighsInt r = 0;
    for (int s = 0; s < S; ++s) for (int j = 0; j < n; ++j) {
      int c = zidx_[(size_t)s * n + j]; if (c < 0) continue;
      mcrow_[(size_t)s * n + j] = r;
      st.push_back((HighsInt)I.size()); L.push_back(0.0); U.push_back(kHighsInf); I.insert(I.end(), {(HighsInt)c, (HighsInt)j}); V.insert(V.end(), {1.0, -Qmin[s]}); ++r;
      st.push_back((HighsInt)I.size()); L.push_back(-kHighsInf); U.push_back(0.0); I.insert(I.end(), {(HighsInt)c, (HighsInt)j}); V.insert(V.end(), {1.0, -Qmax[s]}); ++r;
      st.push_back((HighsInt)I.size()); L.push_back(-Qmax[s]); U.push_back(kHighsInf); I.insert(I.end(), {(HighsInt)c, (HighsInt)(n + s), (HighsInt)j}); V.insert(V.end(), {1.0, -1.0, -Qmax[s]}); ++r;
      st.push_back((HighsInt)I.size()); L.push_back(-kHighsInf); U.push_back(-Qmin[s]); I.insert(I.end(), {(HighsInt)c, (HighsInt)(n + s), (HighsInt)j}); V.insert(V.end(), {1.0, -1.0, -Qmin[s]}); ++r;
    }
    h_.addRows(r, L.data(), U.data(), (HighsInt)I.size(), st.data(), I.data(), V.data());
    count_.assign(S, 0);
  }

  void set_y_bounds(const double* lb, const double* ub) {
    std::vector<HighsInt> idx(n); for (int j = 0; j < n; ++j) idx[j] = j;
    h_.changeColsBounds(n, idx.data(), lb, ub);
  }

  // Node-local recourse bounds: theta_s in [Qmin_s, Qmax_s], z_sj <= Qmax_s, envelope rows rewritten.
  void set_theta_bounds(const std::vector<double>& Qmin, const std::vector<double>& Qmax) {
    std::vector<HighsInt> idx(S); for (int s = 0; s < S; ++s) idx[s] = n + s;
    h_.changeColsBounds(S, idx.data(), Qmin.data(), Qmax.data());
    for (int s = 0; s < S; ++s) for (int j = 0; j < n; ++j) {
      int c = zidx_[(size_t)s * n + j]; if (c < 0) continue;
      HighsInt r = mcrow_[(size_t)s * n + j];
      h_.changeColBounds(c, 0.0, Qmax[s]);
      h_.changeCoeff(r, j, -Qmin[s]);
      h_.changeCoeff(r + 1, j, -Qmax[s]);
      h_.changeCoeff(r + 2, j, -Qmax[s]); h_.changeRowBounds(r + 2, -Qmax[s], kHighsInf);
      h_.changeCoeff(r + 3, j, -Qmin[s]); h_.changeRowBounds(r + 3, -kHighsInf, -Qmin[s]);
    }
    Qmin_ = Qmin; Qmax_ = Qmax;
  }

  // Benders rows theta_s + (beta o cap)^T y >= alpha^T d_s for every pool vertex not yet in the master.
  long add_cuts(const CutPool& pool) {
    std::vector<double> L, U, V; std::vector<HighsInt> st, I; long added = 0;
    for (int s = 0; s < S; ++s) {
      const auto& cs = pool.cuts(s);
      for (size_t k = count_[s]; k < cs.size(); ++k) {
        double ad = 0.0; for (int q = 0; q < inst_.m; ++q) ad += cs[k].alpha[q] * inst_.ds(s)[q];
        st.push_back((HighsInt)I.size()); L.push_back(ad); U.push_back(kHighsInf);
        I.push_back(n + s); V.push_back(1.0);
        for (int j = 0; j < n; ++j) { double g = cs[k].beta[j] * inst_.cap[j]; if (std::fabs(g) > 1e-12) { I.push_back(j); V.push_back(g); } }
        ++added;
      }
      count_[s] = (long)cs.size();
    }
    if (added) h_.addRows((HighsInt)L.size(), L.data(), U.data(), (HighsInt)I.size(), st.data(), I.data(), V.data());
    return added;
  }

  // Warm dual simplex; on a stall, the basis in a fresh solver state, then IPM.
  HighsModelStatus solve(double& obj, std::vector<double>& x, double time_limit) {
    h_.setOptionValue("time_limit", time_limit > 0 ? h_.getRunTime() + time_limit : kHighsInf);   // HiGHS checks time_limit against its cumulative clock
    h_.run(); HighsModelStatus st = h_.getModelStatus();
    if (st != HighsModelStatus::kOptimal && st != HighsModelStatus::kInfeasible) {
      HighsBasis b = h_.getBasis(); h_.clearSolver(); h_.setBasis(b); h_.run(); st = h_.getModelStatus();
    }
    if (st != HighsModelStatus::kOptimal && st != HighsModelStatus::kInfeasible) {
      h_.setOptionValue("solver", "ipm"); h_.run(); st = h_.getModelStatus(); h_.setOptionValue("solver", "simplex");
    }
    if (st == HighsModelStatus::kOptimal) { obj = h_.getInfo().objective_function_value; x = h_.getSolution().col_value; }
    return st;
  }

  long num_rows() const { return h_.getNumRow(); }
  int n, S, nvar;

  // ---- MIP mode: y integer, HiGHS branch-and-bound on the whole master
  void make_integer(int threads) {
    std::vector<HighsInt> idx(n); std::vector<HighsVarType> vt(n, HighsVarType::kInteger);
    for (int j = 0; j < n; ++j) idx[j] = j;
    h_.changeColsIntegrality(n, idx.data(), vt.data());
    h_.setOptionValue("mip_rel_gap", 1e-6); h_.setOptionValue("mip_abs_gap", 1e-6);
    (void)threads;   // HiGHS' global scheduler is already initialised by the recourse LP objects; the option cannot be changed
    if (std::getenv("ENDOG_MIP_LOG")) h_.setOptionValue("output_flag", true);
  }
  // Returns the model status; obj = objective of the best integer solution (if any), bound = MIP dual bound.
  HighsModelStatus solve_mip(double& obj, double& bound, std::vector<double>& x, double time_limit) {
    h_.clearSolver();   // every round is a fresh MIP solve; also drops any stale solution
    h_.setOptionValue("time_limit", time_limit > 0 ? h_.getRunTime() + time_limit : kHighsInf);
    h_.run(); HighsModelStatus st = h_.getModelStatus();
    const HighsInfo& info = h_.getInfo();
    bound = info.mip_dual_bound;
    bool feas = info.primal_solution_status == kSolutionStatusFeasible;
    if (feas) { obj = info.objective_function_value; x = h_.getSolution().col_value; } else { obj = INFINITY; x.clear(); }
    return st;
  }

 private:
  const Instance& inst_;
  Highs h_;
  std::vector<double> Qmin_, Qmax_;
  std::vector<int> zidx_; std::vector<HighsInt> mcrow_;
  std::vector<long> count_;
};

struct McbNode { double bound; long id; std::vector<double> lb, ub, qlo, qhi; int side; };   // side: -1 root, 0 ub changed, 1 lb changed
struct McbCmp { bool operator()(const McbNode& a, const McbNode& b) const { return a.bound > b.bound || (a.bound == b.bound && a.id > b.id); } };

}  // namespace

BncStats branch_and_cut_mcb(const Instance& inst, const BncOptions& opt, bool node_local, const CutPool* init) {
  const double t0 = now();
  const int n = inst.n, S = inst.S, m = inst.m; const double tol = opt.tol;
  BncStats st; st.name = inst.name; st.n = n; st.m = m; st.S = S; st.endo = inst.endo; st.workers = opt.workers; st.opt = opt;
  DualSet ds(inst, opt.workers);
  CutPool pool(S);
  std::vector<double> ones(n, 1.0), zeros(n, 0.0), Qmin(S), Qmax(S);
  {
    std::vector<double> al, be;
    for (int s = 0; s < S; ++s) { Qmin[s] = ds.duals[0]->Q(s, ones.data(), al, be); pool.add(s, al, be); }
    for (int s = 0; s < S; ++s) { ds.duals[0]->Q(s, zeros.data(), al, be); pool.add(s, al, be); Qmax[s] = 0.0; for (int k = 0; k < m; ++k) Qmax[s] += inst.pen * inst.ds(s)[k]; }
  }
  if (init) for (int s = 0; s < S; ++s) for (const auto& v : init->cuts(s)) pool.add(s, v.alpha, v.beta);
  McbMaster M(inst, Qmin, Qmax);
  M.add_cuts(pool);
  auto wall_left = [&] { return opt.wall > 0 ? opt.wall - (now() - t0) : -1.0; };
  auto wall_hit = [&] { return opt.wall > 0 && now() - t0 > opt.wall; };
  auto ub_cut = [&](double UB) { return UB - tol * std::max(1.0, std::fabs(UB)); };
  // separate Benders cuts at y against theta: S recourse LPs; fills Qv; returns the number of new vertices
  auto separate = [&](const std::vector<double>& x, std::vector<double>& Qv) {
    std::vector<std::vector<DualVertex>> found(S); Qv.assign(S, 0.0);
    for_each_scenario(S, ds, [&](int s, RecourseDual& dual) {
      std::vector<double> al, be; Qv[s] = dual.Q(s, x.data(), al, be);
      if (Qv[s] > x[n + s] + 1e-7) found[s].push_back({al, be});
    });
    long added = 0;
    for (int s = 0; s < S; ++s) for (auto& v : found[s]) added += pool.add(s, v.alpha, v.beta);
    st.sep_lps += S;
    return added;
  };
  double obj = 0.0; std::vector<double> x;
  bool proven = true;
  // ---- root: converged LP relaxation with Benders cuts separated at the fractional point
  int root_iters = 0;
  {
    std::vector<double> Qv;
    for (root_iters = 1; root_iters <= opt.root_max_rounds; ++root_iters) {
      double t1 = now(); HighsModelStatus hs = M.solve(obj, x, wall_left()); st.lp_time += now() - t1; ++st.lps;
      if (hs != HighsModelStatus::kOptimal) throw LPTimeout("mcb root LP: " + std::string(hs == HighsModelStatus::kTimeLimit ? "time limit" : "not optimal"));
      t1 = now(); long added = separate(x, Qv); st.sep_time += now() - t1;
      if (added == 0) break;
      M.add_cuts(pool);
      if (wall_hit()) break;
    }
  }
  st.root_value = obj; st.root_iters = root_iters; st.root_time = now() - t0; st.root_cuts = (long)pool.total();
  if (opt.verbose) std::printf("  [%s] root %.6f after %d rounds, %zu cuts, %.0fs\n", node_local ? "mcb-dyn" : "mcb", obj, root_iters, pool.total(), now() - t0);
  // ---- tree
  double UB = INFINITY; std::vector<double> inc;
  std::priority_queue<McbNode, std::vector<McbNode>, McbCmp> heap;
  long counter = 0;
  heap.push({obj, counter++, zeros, ones, Qmin, Qmax, -1});
  std::set<std::vector<char>> tried;
  double global_lb = obj;
  auto round_once = [&](const std::vector<double>& y, std::vector<double>& yb, double& F) {
    yb.resize(n); std::vector<char> key(n);
    for (int j = 0; j < n; ++j) { yb[j] = std::nearbyint(std::min(std::max(y[j], 0.0), 1.0)); key[j] = (char)yb[j]; }
    if (!tried.insert(key).second) return false;
    double t1 = now();
    std::vector<double> Qv(S); std::vector<DualVertex> vert(S);
    for_each_scenario(S, ds, [&](int s, RecourseDual& dual) { Qv[s] = dual.Q(s, yb.data(), vert[s].alpha, vert[s].beta); });
    F = 0.0; for (int j = 0; j < n; ++j) F += inst.f[j] * yb[j];
    for (int s = 0; s < S; ++s) F += inst.p(s, yb.data()) * Qv[s];
    long added = 0; for (int s = 0; s < S; ++s) added += pool.add(s, vert[s].alpha, vert[s].beta);   // free Benders cuts
    if (added) M.add_cuts(pool);
    st.sep_time += now() - t1; st.sep_lps += S;
    return true;
  };
  std::vector<double> Qv;
  while (!heap.empty() && proven && !opt.root_only) {
    if (wall_hit()) { proven = false; break; }
    McbNode nd = heap.top(); heap.pop();
    global_lb = nd.bound;
    if (nd.bound >= ub_cut(UB)) continue;
    ++st.nodes;
    M.set_y_bounds(nd.lb.data(), nd.ub.data());
    if (node_local) {
      double t1 = now();
      std::vector<double> al, be;
      if (nd.side == 0) { for_each_scenario(S, ds, [&](int s, RecourseDual& dual) { std::vector<double> a2, b2; nd.qlo[s] = dual.Q(s, nd.ub.data(), a2, b2); }); st.sep_lps += S; }
      else if (nd.side == 1) { for_each_scenario(S, ds, [&](int s, RecourseDual& dual) { std::vector<double> a2, b2; nd.qhi[s] = dual.Q(s, nd.lb.data(), a2, b2); }); st.sep_lps += S; }
      M.set_theta_bounds(nd.qlo, nd.qhi);
      st.sep_time += now() - t1;
    }
    int rounds = 0;
    while (true) {
      if (wall_hit()) { proven = false; break; }
      double t1 = now(); HighsModelStatus hs = M.solve(obj, x, wall_left()); st.lp_time += now() - t1; ++st.lps;
      st.max_rows = std::max(st.max_rows, M.num_rows());
      if (hs == HighsModelStatus::kInfeasible) break;
      if (hs != HighsModelStatus::kOptimal) { proven = false; if (opt.verbose) std::printf("  [mcb] node %ld: LP not optimal, stopping\n", st.nodes); break; }
      if (obj >= ub_cut(UB)) break;
      double fmax = 0.0; int jmax = 0;
      for (int j = 0; j < n; ++j) { double fr = std::fabs(x[j] - std::nearbyint(x[j])); if (fr > fmax) { fmax = fr; jmax = j; } }
      if (fmax < 1e-6) {
        // integer point: S recourse LPs against theta; new cuts re-solve the node, otherwise it is an incumbent
        t1 = now(); long added = separate(x, Qv); st.sep_time += now() - t1;
        if (added) { M.add_cuts(pool); continue; }
        double F = 0.0; std::vector<double> yb(n);
        for (int j = 0; j < n; ++j) { yb[j] = std::nearbyint(x[j]); F += inst.f[j] * yb[j]; }
        for (int s = 0; s < S; ++s) F += inst.p(s, yb.data()) * Qv[s];
        if (F < UB) { UB = F; inc = yb; if (opt.verbose) std::printf("  [mcb] node %ld: incumbent %.6f (lb %.6f) %.0fs\n", st.nodes, UB, global_lb, now() - t0); }
        break;
      }
      if (opt.heur) {
        std::vector<double> yb; double F;
        if (round_once(x, yb, F) && F < UB) {
          UB = F; inc = yb; ++st.heur_hits;
          if (opt.verbose) std::printf("  [mcb] node %ld: rounding incumbent %.6f (lb %.6f) %.0fs\n", st.nodes, UB, global_lb, now() - t0);
          if (obj >= ub_cut(UB)) break;
        }
      }
      if (rounds < 1) {   // one fractional separation round per node
        t1 = now(); long added = separate(x, Qv); st.sep_time += now() - t1; ++rounds;
        if (added) { M.add_cuts(pool); continue; }
      }
      McbNode c0{obj, counter++, nd.lb, nd.ub, nd.qlo, nd.qhi, 0}; c0.ub[jmax] = 0.0;
      McbNode c1{obj, counter++, nd.lb, nd.ub, nd.qlo, nd.qhi, 1}; c1.lb[jmax] = 1.0;
      heap.push(std::move(c0)); heap.push(std::move(c1));
      break;
    }
    if (opt.verbose && st.nodes % 50 == 0)
      std::printf("  [mcb] nodes %ld open %zu lb %.6f ub %.6f rows %ld %.0fs\n", st.nodes, heap.size(), global_lb, UB, M.num_rows(), now() - t0);
  }
  if (!heap.empty()) { proven = false; global_lb = std::min(global_lb, heap.top().bound); } else if (proven) global_lb = UB;
  st.OPT = UB; st.lb = global_lb; st.y = inc; st.proven = proven;
  st.n_open = inc.empty() ? -1 : (int)std::count(inc.begin(), inc.end(), 1.0);
  st.total_time = now() - t0; st.cuts = (long)pool.total(); st.rows = M.num_rows(); st.vars = M.nvar;
  return st;
}

BncStats mcb_mip_loop(const Instance& inst, const BncOptions& opt) {
  const double t0 = now();
  const int n = inst.n, S = inst.S, m = inst.m; const double tol = opt.tol;
  BncStats st; st.name = inst.name; st.n = n; st.m = m; st.S = S; st.endo = inst.endo; st.workers = opt.workers; st.opt = opt;
  DualSet ds(inst, opt.workers);
  CutPool pool(S);
  std::vector<double> ones(n, 1.0), zeros(n, 0.0), Qmin(S), Qmax(S);
  {
    std::vector<double> al, be;
    for (int s = 0; s < S; ++s) { Qmin[s] = ds.duals[0]->Q(s, ones.data(), al, be); pool.add(s, al, be); }
    for (int s = 0; s < S; ++s) { ds.duals[0]->Q(s, zeros.data(), al, be); pool.add(s, al, be); Qmax[s] = 0.0; for (int k = 0; k < m; ++k) Qmax[s] += inst.pen * inst.ds(s)[k]; }
  }
  McbMaster M(inst, Qmin, Qmax);
  M.add_cuts(pool);
  auto wall_left = [&] { return opt.wall > 0 ? opt.wall - (now() - t0) : -1.0; };
  auto wall_hit = [&] { return opt.wall > 0 && now() - t0 > opt.wall; };
  // S recourse LPs at y against theta; fills Qv; returns the number of new pool vertices
  auto separate = [&](const std::vector<double>& x, std::vector<double>& Qv) {
    std::vector<std::vector<DualVertex>> found(S); Qv.assign(S, 0.0);
    for_each_scenario(S, ds, [&](int s, RecourseDual& dual) {
      std::vector<double> al, be; Qv[s] = dual.Q(s, x.data(), al, be);
      if (Qv[s] > x[n + s] + 1e-7) found[s].push_back({al, be});
    });
    long added = 0;
    for (int s = 0; s < S; ++s) for (auto& v : found[s]) added += pool.add(s, v.alpha, v.beta);
    st.sep_lps += S;
    return added;
  };
  double obj = 0.0; std::vector<double> x, Qv;
  // ---- root: the same converged LP relaxation as MC-B (Table 1), its pool seeds the MIP
  int root_iters = 0;
  for (root_iters = 1; root_iters <= opt.root_max_rounds; ++root_iters) {
    double t1 = now(); HighsModelStatus hs = M.solve(obj, x, wall_left()); st.lp_time += now() - t1; ++st.lps;
    if (hs != HighsModelStatus::kOptimal) throw LPTimeout("mcb-mip root LP: not optimal");
    t1 = now(); long added = separate(x, Qv); st.sep_time += now() - t1;
    if (added == 0) break;
    M.add_cuts(pool);
    if (wall_hit()) break;
  }
  st.root_value = obj; st.root_iters = root_iters; st.root_time = now() - t0; st.root_cuts = (long)pool.total();
  if (opt.verbose) std::printf("  [mcb-mip] root LP %.6f after %d rounds, %zu cuts, %.0fs\n", obj, root_iters, pool.total(), now() - t0);
  // ---- MIP rounds
  M.make_integer(opt.workers);
  double UB = INFINITY, LB = obj; std::vector<double> inc; bool proven = false;
  for (;;) {
    if (wall_hit()) break;
    double limit = opt.mip_time; if (opt.wall > 0) limit = std::min(limit, wall_left());
    double mobj, bound; double t1 = now();
    HighsModelStatus hs = M.solve_mip(mobj, bound, x, limit);
    st.lp_time += now() - t1; ++st.lps; ++st.nodes;
    st.max_rows = std::max(st.max_rows, M.num_rows());
    if (hs == HighsModelStatus::kInfeasible) { if (opt.verbose) std::printf("  [mcb-mip] MIP infeasible\n"); break; }
    bool optimal = hs == HighsModelStatus::kOptimal;
    if (std::isfinite(bound)) LB = std::max(LB, bound);   // the master is a relaxation: its dual bound is a global lower bound
    if (x.empty()) { if (opt.verbose) std::printf("  [mcb-mip] round %ld: no integer solution within the MIP time limit (bound %.6f) %.0fs\n", st.nodes, bound, now() - t0); if (wall_hit()) break; else continue; }
    std::vector<double> yb(n); for (int j = 0; j < n; ++j) { yb[j] = std::nearbyint(x[j]); x[j] = yb[j]; }
    t1 = now(); long added = separate(x, Qv); st.sep_time += now() - t1;
    double F = 0.0; for (int j = 0; j < n; ++j) F += inst.f[j] * yb[j];
    for (int s = 0; s < S; ++s) F += inst.p(s, yb.data()) * Qv[s];
    if (F < UB) { UB = F; inc = yb; }
    if (opt.verbose) std::printf("  [mcb-mip] round %ld: mip %.6f (%s, bound %.6f) F %.6f best %.6f lb %.6f cuts +%ld rows %ld %.0fs\n",
                                 st.nodes, mobj, optimal ? "optimal" : "time limit", bound, F, UB, LB, added, M.num_rows(), now() - t0);
    if (UB - LB <= tol * std::max(1.0, std::fabs(UB))) { proven = true; LB = UB; break; }
    if (added == 0 && optimal) { proven = true; LB = UB; break; }
    if (added) M.add_cuts(pool);
  }
  st.OPT = UB; st.lb = LB; st.y = inc; st.proven = proven;
  st.n_open = inc.empty() ? -1 : (int)std::count(inc.begin(), inc.end(), 1.0);
  st.total_time = now() - t0; st.cuts = (long)pool.total(); st.rows = M.num_rows(); st.vars = M.nvar;
  return st;
}
