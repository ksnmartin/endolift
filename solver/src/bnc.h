// Branch-and-cut with the decomposed node bound (Sections 4-5 of the paper).
// The root loop and the per-node work are factored out so the single-rank tree (branch_and_cut_dec) and the MPI
// coordinator/worker tree (mpi_main.cpp) run exactly the same node code.
#pragma once
#include <memory>
#include <set>
#include <string>
#include <vector>
#include "decomp.h"
#include "instance.h"
#include "recourse.h"

struct BncOptions {
  int workers = 1;
  double wall = -1.0;            // seconds, <= 0: none
  double tol = 1e-6, lam = 0.2, node_tol = 1e-5, root_loose = 1e-4;   // root_loose <= 0 disables
  int max_age = 5, root_max_rounds = 200;
  bool heur = true, root_only = false, verbose = false, sparse = false, lazy_mc = false, unlinked = false;
  double mip_time = 1800.0;      // mcb-mip: HiGHS time limit per MIP solve
};

struct BncStats {
  std::string name; int n = 0, m = 0, S = 0; std::string endo; int workers = 1, ranks = 1;
  double OPT = INFINITY, lb = -INFINITY; bool proven = false; int n_open = -1;
  double root_value = 0, root_time = 0; int root_iters = 0; long root_cuts = 0;
  long nodes = 0, lps = 0, sep_lps = 0, cuts = 0, rows = 0, vars = 0, max_rows = 0, rows_deleted = 0;
  double lp_time = 0, sep_time = 0, total_time = 0;
  int heur_hits = 0, pruned_early = 0;
  std::vector<double> y;
  BncOptions opt;
  std::string csv_header() const;
  std::string csv_row(const char* master = "dec-cpp") const;
  std::string summary(const char* master = "dec-cpp") const;
};

struct NodeTask { double bound; long id; std::vector<double> lb, ub; };

struct NodeResult {
  bool infeasible = false, pruned = false;
  double lb = -INFINITY;                 // node bound after the last bound run
  int branch_var = -1;                   // >= 0: two children (ub[j] = 0, lb[j] = 1) with bound lb
  bool has_inc = false; double F = INFINITY; std::vector<double> y;   // incumbent found at this node (F < UB at the time)
  std::vector<std::vector<DualVertex>> new_cuts;   // per scenario, dual vertices added to the pool at this node
  int iters = 0, heur_hits = 0; long sep_lps = 0, rows = 0; double lp_time = 0, sep_time = 0;
};

// Everything one rank needs to bound nodes: the master, the recourse duals, a replica of the cut pool.
class NodeWorker {
 public:
  NodeWorker(const Instance& inst, const BncOptions& opt);
  // Root: decomposed bound <-> lifted-cut separation until no cut is violated. Fills stats.root_*; returns the last bound.
  BoundResult run_root(BncStats& stats, double t0);
  // Round the given point once (per distinct rounded point); returns true with F if it was new.
  bool try_round(const std::vector<double>& y, std::vector<double>& yb, double& F, NodeResult* acc = nullptr);
  // Process one node against the current incumbent value UB (cutoff + pruning). Pool cuts found here are added
  // locally and reported in the result.
  NodeResult process(const NodeTask& task, double UB, long node_index, double t0);
  // Merge cuts received from elsewhere (dedup) into the pool and the blocks; returns the number new to this rank.
  long absorb_cuts(const std::vector<std::vector<DualVertex>>& cuts);

  const Instance& inst; BncOptions opt;
  DualSet ds; CutPool pool; std::unique_ptr<DecomposedBound> D;
 private:
  std::set<std::vector<char>> tried_;
};

BncStats branch_and_cut_dec(const Instance& inst, const BncOptions& opt);
