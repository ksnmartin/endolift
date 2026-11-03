#include "bnc.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <queue>
#include <sstream>

static double now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }

// Lifted-cut separation at the decomposed point: S (1 + 2n) recourse LPs, parallel over scenarios. Returns cuts added.
static long separate(const Instance& inst, DualSet& ds, CutPool& pool, DecomposedBound& D, std::vector<double> y, std::vector<double> Yp, long& lps, double tol = 1e-7) {
  const int n = inst.n, S = inst.S;
  D.project(y, Yp);          // clean LP tolerances: Y_jl < 0 or > y_l by 1e-9 makes a factor recourse LP infeasible
  std::vector<double> Y; D.unfold(y.data(), Yp.data(), Y);
  std::vector<std::vector<DualVertex>> found(S);
  std::vector<long> cnt(S, 0);
  for_each_scenario(S, ds, [&](int s, RecourseDual& dual) {
    double e0; std::vector<double> Lp, Lm, al, be, r(inst.m), c(n);
    D.blocks.levels(s, y.data(), Y.data(), e0, Lp, Lm);
    const double* d = inst.ds(s);
    if (dual.Q(s, y.data(), al, be) > e0 + tol) found[s].push_back({al, be});
    cnt[s] = 1;
    for (int j : D.blocks.infl) {          // sparse lifting: only the factors the master carries are separated
      const double* Yj = &Y[(size_t)j * n];
      for (int k = 0; k < inst.m; ++k) r[k] = y[j] * d[k];
      for (int i = 0; i < n; ++i) c[i] = inst.cap[i] * Yj[i];
      if (dual.solve(r.data(), c.data(), al, be) > Lp[j] + tol) found[s].push_back({al, be});
      for (int k = 0; k < inst.m; ++k) r[k] = (1.0 - y[j]) * d[k];
      for (int i = 0; i < n; ++i) c[i] = inst.cap[i] * (y[i] - Yj[i]);
      if (dual.solve(r.data(), c.data(), al, be) > Lm[j] + tol) found[s].push_back({al, be});
      cnt[s] += 2;
    }
  });
  long added = 0;
  for (int s = 0; s < S; ++s) { lps += cnt[s]; for (auto& v : found[s]) added += pool.add(s, v.alpha, v.beta); }
  return added;
}

static double ub_cut(double UB, double tol) { return UB - tol * std::max(1.0, std::fabs(UB)); }

NodeWorker::NodeWorker(const Instance& inst_, const BncOptions& opt_) : inst(inst_), opt(opt_), ds(inst_, opt_.workers), pool(inst_.S) {
  const int n = inst.n, S = inst.S;
  std::vector<double> ones(n, 1.0), zeros(n, 0.0), al, be;
  for (const double* y0 : {ones.data(), zeros.data()})
    for (int s = 0; s < S; ++s) { ds.duals[0]->Q(s, y0, al, be); pool.add(s, al, be); }
  D.reset(new DecomposedBound(inst, pool, opt.lam, opt.max_age, 1, true, opt.sparse, opt.lazy_mc));
  D->blocks.unlinked = opt.unlinked;
}

BoundResult NodeWorker::run_root(BncStats& st, double t0) {
  const double tol = opt.tol;
  int root_rounds = 0; double rtol = opt.root_loose > 0 ? opt.root_loose : tol;
  BoundResult r;
  while (true) {
    double t1 = now(); r = D->run(nullptr, nullptr, rtol); st.lp_time += now() - t1; st.lps += r.iters;
    ++root_rounds;
    t1 = now(); long added = separate(inst, ds, pool, *D, r.y, r.Yp, st.sep_lps); st.sep_time += now() - t1;
    if (opt.verbose && (root_rounds % 10 == 0 || added == 0))
      std::printf("  [dec] root round %d: bound %.6f +%ld cuts (pool %zu, master rows %ld, tol %g) %.0fs\n", root_rounds, r.lb, added, pool.total(), r.rows, rtol, now() - t0);
    if (added == 0 && rtol > tol) { rtol = tol; continue; }       // no cut at the loose bound: tighten and separate once more
    if (added == 0 || root_rounds >= opt.root_max_rounds) break;
    D->blocks.add_cuts(pool);
    if (opt.wall > 0 && now() - t0 > opt.wall) break;
  }
  st.root_value = r.lb; st.root_iters = root_rounds; st.root_time = now() - t0; st.root_cuts = (long)pool.total();
  st.max_rows = std::max(st.max_rows, r.rows); st.vars = D->nvar;
  return r;
}

bool NodeWorker::try_round(const std::vector<double>& y, std::vector<double>& yb, double& F, NodeResult* acc) {
  const int n = inst.n;
  yb.resize(n); std::vector<char> key(n);
  for (int j = 0; j < n; ++j) { yb[j] = std::nearbyint(std::min(std::max(y[j], 0.0), 1.0)); key[j] = (char)yb[j]; }
  if (!tried_.insert(key).second) return false;
  // Exact evaluation with S recourse LPs; the optimal dual vertices are free and go into the pool (referee item), so
  // the pool holds an optimal vertex for every integer point the heuristic has visited.
  const int S = inst.S;
  double t1 = now();
  std::vector<double> Qv(S); std::vector<DualVertex> vert(S);
  for_each_scenario(S, ds, [&](int s, RecourseDual& dual) { Qv[s] = dual.Q(s, yb.data(), vert[s].alpha, vert[s].beta); });
  F = 0.0;
  for (int i = 0; i < n; ++i) F += inst.f[i] * yb[i];
  for (int s = 0; s < S; ++s) F += inst.p(s, yb.data()) * Qv[s];
  long added = 0;
  for (int s = 0; s < S; ++s) if (pool.add(s, vert[s].alpha, vert[s].beta)) { ++added; if (acc && (int)acc->new_cuts.size() == S) acc->new_cuts[s].push_back(vert[s]); }
  if (added) D->blocks.add_cuts(pool);
  if (acc) { acc->sep_time += now() - t1; acc->sep_lps += S; }
  return true;
}

long NodeWorker::absorb_cuts(const std::vector<std::vector<DualVertex>>& cuts) {
  long added = 0;
  for (int s = 0; s < inst.S && s < (int)cuts.size(); ++s) for (const auto& v : cuts[s]) added += pool.add(s, v.alpha, v.beta);
  if (added) D->blocks.add_cuts(pool);
  return added;
}

NodeResult NodeWorker::process(const NodeTask& nd, double UB, long node_index, double t0) {
  const int n = inst.n, S = inst.S; const double tol = opt.tol;
  NodeResult R; R.new_cuts.assign(S, {});
  std::vector<double> yb; double F;
  while (true) {
    double t1 = now(); BoundResult r = D->run(nd.lb.data(), nd.ub.data(), opt.node_tol, std::isfinite(UB) ? UB : INFINITY);
    R.lp_time += now() - t1; R.iters += r.iters; R.rows = std::max(R.rows, r.rows); R.lb = r.lb;
    if (r.pruned) R.pruned = true;
    if (r.infeasible) { R.infeasible = true; break; }
    if (r.lb >= ub_cut(UB, tol)) break;
    const std::vector<double>& y = r.y;
    double fmax = 0.0; int jmax = 0;
    for (int j = 0; j < n; ++j) { double fr = std::fabs(y[j] - std::nearbyint(y[j])); if (fr > fmax) { fmax = fr; jmax = j; } }
    if (opt.heur && fmax >= 1e-6 && try_round(y, yb, F, &R) && F < UB) {
      UB = F; R.has_inc = true; R.F = F; R.y = yb; ++R.heur_hits;
      if (opt.verbose) std::printf("  [dec] node %ld: rounding incumbent %.6f (node lb %.6f) %.0fs\n", node_index, UB, r.lb, now() - t0);
      if (r.lb >= ub_cut(UB, tol)) break;
    }
    if (fmax < 1e-6) {
      // integer leaf: S recourse LPs vs the block epigraph values
      std::vector<double> ybin(n); for (int j = 0; j < n; ++j) ybin[j] = std::nearbyint(y[j]);
      std::vector<double> Y; D->unfold(ybin.data(), r.Yp.data(), Y);
      std::vector<std::vector<DualVertex>> found(S);
      t1 = now();
      for_each_scenario(S, ds, [&](int s, RecourseDual& dual) {
        double e0; std::vector<double> Lp, Lm, al, be;
        D->blocks.levels(s, ybin.data(), Y.data(), e0, Lp, Lm);
        if (dual.Q(s, ybin.data(), al, be) > e0 + 1e-7) found[s].push_back({al, be});
      });
      long added = 0;
      for (int s = 0; s < S; ++s) for (auto& v : found[s]) if (pool.add(s, v.alpha, v.beta)) { ++added; R.new_cuts[s].push_back(v); }
      R.sep_time += now() - t1; R.sep_lps += S;
      if (added) { D->blocks.add_cuts(pool); continue; }      // same node, re-bound from the warm center
      double Fy = F_exact(inst, ds, ybin.data());
      if (Fy < UB) {
        UB = Fy; R.has_inc = true; R.F = Fy; R.y = ybin;
        if (opt.verbose) std::printf("  [dec] node %ld: incumbent %.6f (node lb %.6f) %.0fs\n", node_index, UB, r.lb, now() - t0);
      }
      break;
    }
    R.branch_var = jmax;
    break;
  }
  return R;
}

struct NodeCmp { bool operator()(const NodeTask& a, const NodeTask& b) const { return a.bound > b.bound || (a.bound == b.bound && a.id > b.id); } };

BncStats branch_and_cut_dec(const Instance& inst, const BncOptions& opt) {
  const double t0 = now(); const int n = inst.n, S = inst.S; const double tol = opt.tol;
  BncStats st; st.name = inst.name; st.n = n; st.m = inst.m; st.S = S; st.endo = inst.endo; st.workers = opt.workers; st.opt = opt;
  NodeWorker W(inst, opt);
  double UB = INFINITY; std::vector<double> inc; double global_lb = -INFINITY; bool proven = true;
  auto wall_hit = [&] { return opt.wall > 0 && now() - t0 > opt.wall; };
  try {
    BoundResult r = W.run_root(st, t0);
    if (opt.heur) {
      std::vector<double> yb; double F; NodeResult acc;
      if (W.try_round(r.y, yb, F, &acc)) {
        inc = yb; UB = F; ++st.heur_hits; st.sep_time += acc.sep_time; st.sep_lps += acc.sep_lps;
        if (opt.verbose) std::printf("  [dec] root rounding incumbent %.6f (lb %.6f)\n", UB, r.lb);
      }
    }
    std::priority_queue<NodeTask, std::vector<NodeTask>, NodeCmp> heap;
    long counter = 0;
    heap.push({r.lb, counter++, std::vector<double>(n, 0.0), std::vector<double>(n, 1.0)});
    global_lb = r.lb;
    while (!heap.empty() && proven && !opt.root_only) {
      if (wall_hit()) { proven = false; break; }
      NodeTask nd = heap.top(); heap.pop();
      global_lb = nd.bound;
      if (nd.bound >= ub_cut(UB, tol)) continue;
      ++st.nodes;
      NodeResult R = W.process(nd, UB, st.nodes, t0);
      st.lp_time += R.lp_time; st.sep_time += R.sep_time; st.lps += R.iters; st.sep_lps += R.sep_lps;
      st.heur_hits += R.heur_hits; st.pruned_early += R.pruned; st.max_rows = std::max(st.max_rows, R.rows);
      if (R.has_inc && R.F < UB) { UB = R.F; inc = R.y; }
      if (R.branch_var >= 0) {
        NodeTask c0{R.lb, counter++, nd.lb, nd.ub}; c0.ub[R.branch_var] = 0.0;
        NodeTask c1{R.lb, counter++, nd.lb, nd.ub}; c1.lb[R.branch_var] = 1.0;
        heap.push(std::move(c0)); heap.push(std::move(c1));
      }
      if (opt.verbose && st.nodes % 50 == 0)
        std::printf("  [dec] nodes %ld open %zu lb %.6f ub %.6f rows %ld %.0fs\n", st.nodes, heap.size(), global_lb, UB, R.rows, now() - t0);
    }
    if (!heap.empty()) global_lb = std::min(global_lb, heap.top().bound);
    else global_lb = UB;
  } catch (const LPTimeout& e) {
    std::printf("  [dec] LP timeout, run left unproven: %s\n", e.what()); proven = false;
  }
  st.OPT = UB; st.lb = global_lb; st.y = inc; st.proven = proven;
  st.n_open = inc.empty() ? -1 : (int)std::count(inc.begin(), inc.end(), 1.0);
  st.total_time = now() - t0; st.cuts = (long)W.pool.total(); st.rows = W.D->num_rows(); st.rows_deleted = W.D->deleted;
  return st;
}

std::string BncStats::csv_header() const {
  return "name,n,m,S,endo,master,ranks,workers,OPT,lb,proven,n_open,root_value,root_gap,root_iters,root_time,root_cuts,nodes,lps,sep_lps,lp_time,sep_time,total_time,cuts,rows,vars,max_rows,rows_deleted,node_tol,heur,heur_hits,pruned_early,root_loose,max_age,lam,y";
}
std::string BncStats::csv_row(const char* master) const {
  std::ostringstream o; o.precision(15);
  double rg = std::isfinite(OPT) ? (OPT - root_value) / std::fabs(OPT) * 100.0 : NAN;
  o << name << ',' << n << ',' << m << ',' << S << ',' << endo << ',' << master << ',' << ranks << ',' << workers << ',' << OPT << ',' << lb << ',' << (proven ? 1 : 0) << ',' << n_open << ','
    << root_value << ',' << rg << ',' << root_iters << ',' << root_time << ',' << root_cuts << ',' << nodes << ',' << lps << ',' << sep_lps << ','
    << lp_time << ',' << sep_time << ',' << total_time << ',' << cuts << ',' << rows << ',' << vars << ',' << max_rows << ',' << rows_deleted << ','
    << opt.node_tol << ',' << (opt.heur ? 1 : 0) << ',' << heur_hits << ',' << pruned_early << ',' << opt.root_loose << ',' << opt.max_age << ',' << opt.lam << ',';
  for (double v : y) o << (v > 0.5 ? '1' : '0');
  return o.str();
}
std::string BncStats::summary(const char* master) const {
  char buf[512];
  double rg = std::isfinite(OPT) ? (OPT - root_value) / std::fabs(OPT) * 100.0 : NAN;
  std::snprintf(buf, sizeof buf, "%s: OPT %.6f lb %.6f %s  open %d/%d  root gap %.2f%% (%d it, %.0fs)  nodes %ld  lps %ld  sep LPs %ld  lp %.0fs sep %.0fs  total %.0fs  rows %ld (max %ld, deleted %ld)",
                master, OPT, lb, proven ? "proven" : "NOT PROVEN", n_open, n, rg, root_iters, root_time, nodes, lps, sep_lps, lp_time, sep_time, total_time, rows, max_rows, rows_deleted);
  return buf;
}
