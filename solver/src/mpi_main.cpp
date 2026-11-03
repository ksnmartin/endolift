// endog_mpi: node-parallel branch-and-cut. Rank 0 coordinates (heap, incumbent, global cut pool); ranks 1..R-1 bound
// nodes with their own master (NodeWorker). The root runs once on rank 0 (with --root-workers separation threads);
// its dual-vertex pool, active master cut rows and stabilization center are broadcast so every worker starts from the
// converged root state. Cuts found at integer leaves travel to the coordinator with the node result and reach the
// other workers with their next node assignment.
#include <mpi.h>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <queue>
#include <string>
#include <vector>
#include "bnc.h"
#include "instance.h"

static double now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }
enum Tag { TAG_REQ = 1, TAG_NODE = 2, TAG_STOP = 3, TAG_RESULT = 4 };

// Flat double buffer with cursor (everything is shipped as doubles: counts included).
struct Buf {
  std::vector<double> d; size_t pos = 0;
  void put(double v) { d.push_back(v); }
  void put(const std::vector<double>& v) { put((double)v.size()); d.insert(d.end(), v.begin(), v.end()); }
  double get() { return d.at(pos++); }
  std::vector<double> getv() { size_t k = (size_t)get(); std::vector<double> v(d.begin() + pos, d.begin() + pos + k); pos += k; return v; }
  void put_cuts(const std::vector<std::vector<DualVertex>>& cs) {
    put((double)cs.size());
    for (auto& per : cs) { put((double)per.size()); for (auto& v : per) { put(v.alpha); put(v.beta); } }
  }
  std::vector<std::vector<DualVertex>> get_cuts() {
    std::vector<std::vector<DualVertex>> cs((size_t)get());
    for (auto& per : cs) { per.resize((size_t)get()); for (auto& v : per) { v.alpha = getv(); v.beta = getv(); } }
    return cs;
  }
  void send(int to, int tag) { MPI_Send(d.data(), (int)d.size(), MPI_DOUBLE, to, tag, MPI_COMM_WORLD); }
  int recv(int from, int tag, MPI_Status* status) {   // returns the actual tag
    MPI_Status st; MPI_Probe(from, tag, MPI_COMM_WORLD, &st);
    int count; MPI_Get_count(&st, MPI_DOUBLE, &count);
    d.resize(count); pos = 0;
    MPI_Recv(d.data(), count, MPI_DOUBLE, st.MPI_SOURCE, st.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    if (status) *status = st;
    return st.MPI_TAG;
  }
};

static std::vector<std::vector<DualVertex>> pool_slice(const CutPool& pool, std::vector<long>& from) {
  std::vector<std::vector<DualVertex>> out(pool.S());
  for (int s = 0; s < pool.S(); ++s) {
    const auto& cs = pool.cuts(s);
    for (size_t k = from[s]; k < cs.size(); ++k) out[s].push_back(cs[k]);
    from[s] = (long)cs.size();
  }
  return out;
}

struct NodeCmp { bool operator()(const NodeTask& a, const NodeTask& b) const { return a.bound > b.bound || (a.bound == b.bound && a.id > b.id); } };

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  setvbuf(stdout, nullptr, _IOLBF, 0);   // progress lines must reach a redirected log as they happen
  int rank, R; MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_Comm_size(MPI_COMM_WORLD, &R);
  std::string inst_path, out, endo = "strong";
  int n = 15, m = 15, S = 20, seed = 1; double strength = 0.9, density = 1.0, fmult = 1.0;
  BncOptions opt; int root_workers = 0;   // 0: same as --workers
  auto need = [&](int& i) { if (++i >= argc) { if (rank == 0) std::fprintf(stderr, "missing value for %s\n", argv[i - 1]); MPI_Abort(MPI_COMM_WORLD, 2); } return argv[i]; };
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    if (a == "--inst") inst_path = need(i);
    else if (a == "--n") n = std::atoi(need(i));
    else if (a == "--m") m = std::atoi(need(i));
    else if (a == "--S") S = std::atoi(need(i));
    else if (a == "--seed") seed = std::atoi(need(i));
    else if (a == "--endo") endo = need(i);
    else if (a == "--strength") strength = std::atof(need(i));
    else if (a == "--density") density = std::atof(need(i));
    else if (a == "--fmult") fmult = std::atof(need(i));
    else if (a == "--workers") opt.workers = std::atoi(need(i));
    else if (a == "--root-workers") root_workers = std::atoi(need(i));
    else if (a == "--wall") opt.wall = std::atof(need(i));
    else if (a == "--tol") opt.tol = std::atof(need(i));
    else if (a == "--node-tol") opt.node_tol = std::atof(need(i));
    else if (a == "--root-loose") opt.root_loose = std::atof(need(i));
    else if (a == "--max-age") opt.max_age = std::atoi(need(i));
    else if (a == "--lam") opt.lam = std::atof(need(i));
    else if (a == "--no-heur") opt.heur = false;
    else if (a == "--sparse") opt.sparse = true;
    else if (a == "--lazy-mc") opt.lazy_mc = true;
    else if (a == "--out") out = need(i);
    else if (a == "-v") opt.verbose = true;
    else { if (rank == 0) std::fprintf(stderr, "unknown argument %s\n", argv[i]); MPI_Abort(MPI_COMM_WORLD, 2); }
  }
  if (R < 2) { if (rank == 0) std::fprintf(stderr, "endog_mpi needs at least 2 ranks (1 coordinator + workers)\n"); MPI_Abort(MPI_COMM_WORLD, 2); }
  const double t0 = now();
  Instance inst = inst_path.empty() ? make_instance(n, m, S, seed, strength, density, fmult, endo) : load_instance(inst_path);
  n = inst.n; S = inst.S;
  const double tol = opt.tol;
  BncOptions wopt = opt; wopt.verbose = opt.verbose && rank == 0;
  if (rank == 0 && root_workers > 0) wopt.workers = root_workers;
  if (rank == 0) { std::printf("%s\n", inst.describe().c_str()); std::fflush(stdout); }

  // ---- root once, on rank 0; pool + active master cuts + center go to the workers
  BncStats st; st.name = inst.name; st.n = n; st.m = inst.m; st.S = S; st.endo = inst.endo; st.workers = opt.workers; st.ranks = R; st.opt = opt;
  NodeWorker W(inst, wopt);
  BoundResult root;
  int root_ok = 1;
  if (rank == 0) {
    W.D->keep_cuts = true;
    try { root = W.run_root(st, t0); } catch (const LPTimeout& e) { root_ok = 0; std::printf("  [dec] root LP timeout: %s\n", e.what()); }
  }
  MPI_Bcast(&root_ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (!root_ok) { MPI_Finalize(); return 1; }
  {
    Buf b;
    if (rank == 0) {
      std::vector<long> z(S, 0); b.put_cuts(pool_slice(W.pool, z));
      std::vector<MasterCut> cuts; std::vector<double> cy, cYp; bool hc;
      W.D->export_state(cuts, cy, cYp, hc); W.D->keep_cuts = false;
      b.put((double)cuts.size());
      for (auto& c : cuts) { b.put((double)c.s); b.put(c.c); b.put(c.gy); b.put(c.gp); }
      b.put(hc ? 1.0 : 0.0); b.put(cy); b.put(cYp); b.put(root.lb); b.put(root.y);
      if (opt.verbose) std::printf("  [dec] broadcasting root state: %zu vertices, %zu master cuts, %.1f MB\n", W.pool.total(), cuts.size(), b.d.size() * 8.0 / 1e6);
    }
    int cnt = (int)b.d.size(); MPI_Bcast(&cnt, 1, MPI_INT, 0, MPI_COMM_WORLD);
    b.d.resize(cnt); MPI_Bcast(b.d.data(), cnt, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (rank != 0) {
      b.pos = 0; W.absorb_cuts(b.get_cuts());
      std::vector<MasterCut> cuts((size_t)b.get());
      for (auto& c : cuts) { c.s = (int)b.get(); c.c = b.get(); c.gy = b.getv(); c.gp = b.getv(); }
      bool hc = b.get() > 0.5; auto cy = b.getv(); auto cYp = b.getv();
      W.D->import_state(cuts, cy, cYp, hc);
      root.lb = b.get(); root.y = b.getv();
    }
  }

  if (rank == 0) {
    // ------------------------------------------------------------------ coordinator
    double UB = INFINITY; std::vector<double> inc; bool proven = true;
    if (opt.heur) {
      std::vector<double> yb; double F; NodeResult acc;
      if (W.try_round(root.y, yb, F, &acc)) { inc = yb; UB = F; ++st.heur_hits; st.sep_time += acc.sep_time; st.sep_lps += acc.sep_lps;
        if (opt.verbose) std::printf("  [dec] root rounding incumbent %.6f (lb %.6f)\n", UB, root.lb); }
    }
    auto ub_cut = [&] { return UB - tol * std::max(1.0, std::fabs(UB)); };
    auto wall_hit = [&] { return opt.wall > 0 && now() - t0 > opt.wall; };
    std::priority_queue<NodeTask, std::vector<NodeTask>, NodeCmp> heap;
    long counter = 0;
    heap.push({root.lb, counter++, std::vector<double>(n, 0.0), std::vector<double>(n, 1.0)});
    std::map<int, NodeTask> inflight;                 // worker rank -> node it is bounding
    std::vector<std::vector<long>> sent(R, std::vector<long>(S, 0));   // per worker: pool cuts already shipped
    for (int w = 1; w < R; ++w) for (int s = 0; s < S; ++s) sent[w][s] = (long)W.pool.size(s);   // everyone has the root pool
    std::queue<int> idle;
    { Buf b; for (int w = 1; w < R; ++w) { MPI_Status s; b.recv(MPI_ANY_SOURCE, TAG_REQ, &s); idle.push(s.MPI_SOURCE); } }
    double global_lb = root.lb; double lp_sum = 0.0; long busy_max = 0;
    while (true) {
      while (!idle.empty() && !heap.empty() && !wall_hit()) {
        NodeTask nd = heap.top(); heap.pop();
        if (nd.bound >= ub_cut()) continue;
        int w = idle.front(); idle.pop();
        ++st.nodes;
        Buf b; b.put(UB); b.put(nd.bound); b.put((double)st.nodes); b.put(nd.lb); b.put(nd.ub); b.put_cuts(pool_slice(W.pool, sent[w]));
        b.send(w, TAG_NODE);
        inflight[w] = std::move(nd);
      }
      busy_max = std::max(busy_max, (long)inflight.size());
      if (inflight.empty()) break;                      // heap empty (done) or wall hit with nothing in flight
      if (wall_hit()) proven = false;
      Buf b; MPI_Status s; b.recv(MPI_ANY_SOURCE, TAG_RESULT, &s); int w = s.MPI_SOURCE;
      NodeTask nd = std::move(inflight[w]); inflight.erase(w); idle.push(w);
      // parse result
      bool infeasible = b.get() > 0.5, pruned = b.get() > 0.5; double lb = b.get(); int bv = (int)b.get();
      bool has_inc = b.get() > 0.5; double F = b.get(); std::vector<double> y = b.getv();
      auto cuts = b.get_cuts();
      int iters = (int)b.get(), hh = (int)b.get(); long sep_lps = (long)b.get(), rows = (long)b.get(); double lp_t = b.get(), sep_t = b.get();
      st.lps += iters; st.heur_hits += hh; st.sep_lps += sep_lps; st.max_rows = std::max(st.max_rows, rows); lp_sum += lp_t; st.sep_time += sep_t;
      st.pruned_early += pruned; (void)infeasible;
      for (int s2 = 0; s2 < S; ++s2) for (auto& v : cuts[s2]) W.pool.add(s2, v.alpha, v.beta);
      if (has_inc && F < UB) { UB = F; inc = y; if (opt.verbose) std::printf("  [dec] incumbent %.6f from rank %d (node %ld) %.0fs\n", UB, w, st.nodes, now() - t0); }
      if (bv >= 0) {
        NodeTask c0{lb, counter++, nd.lb, nd.ub}; c0.ub[bv] = 0.0;
        NodeTask c1{lb, counter++, nd.lb, nd.ub}; c1.lb[bv] = 1.0;
        heap.push(std::move(c0)); heap.push(std::move(c1));
      }
      global_lb = heap.empty() ? UB : heap.top().bound;
      for (auto& kv : inflight) global_lb = std::min(global_lb, kv.second.bound);
      if (opt.verbose && st.nodes % 50 == 0)
        std::printf("  [dec] nodes %ld open %zu inflight %zu lb %.6f ub %.6f %.0fs\n", st.nodes, heap.size(), inflight.size(), global_lb, UB, now() - t0);
    }
    for (int w = 1; w < R; ++w) { Buf b; b.put(0.0); b.send(w, TAG_STOP); }
    if (!heap.empty()) { proven = false; global_lb = std::min(global_lb, heap.top().bound); } else if (proven) global_lb = UB;
    st.OPT = UB; st.lb = global_lb; st.y = inc; st.proven = proven;
    st.n_open = inc.empty() ? -1 : (int)std::count(inc.begin(), inc.end(), 1.0);
    st.lp_time += lp_sum; st.total_time = now() - t0; st.cuts = (long)W.pool.total(); st.rows = W.D->num_rows(); st.rows_deleted = W.D->deleted;
    std::printf("%s  [ranks %d, max busy %ld]\n", st.summary("dec-mpi").c_str(), R, busy_max);
    if (!out.empty()) {
      std::filesystem::path op(out);
      if (op.has_parent_path()) std::filesystem::create_directories(op.parent_path());
      bool fresh = !std::ifstream(out).good();
      std::ofstream f(out, std::ios::app);
      if (fresh) f << st.csv_header() << '\n';
      f << st.csv_row("dec-mpi") << '\n';
    }
  } else {
    // ------------------------------------------------------------------ worker
    { Buf b; b.put(0.0); b.send(0, TAG_REQ); }
    while (true) {
      Buf b; int tag = b.recv(0, MPI_ANY_TAG, nullptr);
      if (tag == TAG_STOP) break;
      double UB = b.get(); NodeTask nd; nd.bound = b.get(); long idx = (long)b.get(); nd.lb = b.getv(); nd.ub = b.getv();
      W.absorb_cuts(b.get_cuts());
      NodeResult Rz;
      try { Rz = W.process(nd, UB, idx, t0); }
      catch (const LPTimeout& e) { std::fprintf(stderr, "  [rank %d] LP timeout at node %ld: %s\n", rank, idx, e.what()); Rz.branch_var = -1; Rz.lb = nd.bound; }
      Buf r; r.put(Rz.infeasible); r.put(Rz.pruned); r.put(Rz.lb); r.put((double)Rz.branch_var);
      r.put(Rz.has_inc); r.put(Rz.F); r.put(Rz.y); r.put_cuts(Rz.new_cuts);
      r.put((double)Rz.iters); r.put((double)Rz.heur_hits); r.put((double)Rz.sep_lps); r.put((double)Rz.rows); r.put(Rz.lp_time); r.put(Rz.sep_time);
      r.send(0, TAG_RESULT);
    }
  }
  MPI_Finalize();
  return 0;
}
