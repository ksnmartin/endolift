// Scenario-decomposed (Benders) node bound for the linked lifted-cut master: a multicut, in-out stabilized
// cutting-plane master over (y, Yp pairs, theta_s) with the McCormick (BQP) rows, scenario blocks in closed form.
// Scenario cuts do not depend on node bounds, so they stay valid everywhere in the tree; the stabilization center
// persists across nodes (projected onto the new box) and slack cut rows are aged out.
#pragma once
#include <stdexcept>
#include <vector>
#include "Highs.h"
#include "blocks.h"
#include "instance.h"
#include "recourse.h"

struct LPTimeout : std::runtime_error { using std::runtime_error::runtime_error; };

struct MasterCut { int s; double c; std::vector<double> gy, gp; };   // theta_s >= c + gy.y + gp.Yp

struct BoundResult {
  double lb = -INFINITY, ub = INFINITY;
  int iters = 0;
  long cuts = 0, rows = 0, deleted = 0;
  double t_master = 0.0, t_eval = 0.0;
  std::vector<double> y, Yp;
  bool pruned = false, infeasible = false;
};

class DecomposedBound {
 public:
  // sparse: carry only the factors of influencing sites (a_sj != 0 for some s) and the pair columns (j in infl) x (any l).
  // lazy_mc: start without McCormick rows and add violated ones after each master solve.
  DecomposedBound(const Instance& inst, const CutPool& pool, double lam = 0.2, int max_age = 5, int edge_weight = 1, bool reset_each = true,
                  bool sparse = false, bool lazy_mc = false);
  BoundResult run(const double* lb, const double* ub, double tol, double cutoff = INFINITY, int max_iter = 500, bool verbose = false);
  void set_y_bounds(const double* lb, const double* ub);
  void project(std::vector<double>& y, std::vector<double>& Yp) const;    // onto the node box and the McCormick set
  void unfold(const double* y, const double* Yp, std::vector<double>& Y) const;
  double evaluate(const double* y, const double* Yp, std::vector<MasterCut>& cuts);   // sum_s phi_s and one cut per scenario
  long add_cuts(const std::vector<MasterCut>& cuts);

  ScenarioBlocks blocks;
  int n, S, P, nvar;
  std::vector<int> pj, pl;          // pair k = (pj[k], pl[k]), pj < pl
  std::vector<int> pidx;            // (j*n + l) -> pair index or -1
  int pair(int j, int l) const { return pidx[(size_t)std::min(j, l) * n + std::max(j, l)]; }
  long mc_rows_added = 0;
  long ncuts = 0, deleted = 0;
  long num_rows() const { return h_.getNumRow(); }
  double lp_time_limit = 600.0;
  // MPI: the root is computed once on the coordinator and its master state shipped to the workers.
  bool keep_cuts = false;                       // record the active cut rows so export_state can return them
  void export_state(std::vector<MasterCut>& cuts, std::vector<double>& cy, std::vector<double>& cYp, bool& has_center) const;
  void import_state(const std::vector<MasterCut>& cuts, const std::vector<double>& cy, const std::vector<double>& cYp, bool has_center);

 private:
  const Instance& inst_;
  Highs h_;
  double lam_; int max_age_; bool reset_each_, lazy_mc_;
  std::vector<char> mc_have_;        // lazy_mc: which of the 3P McCormick rows are in the master
  int add_violated_mc(const double* y, const double* Yp);
  std::vector<double> cost_, lb_, ub_;
  std::vector<double> row_lo_; std::vector<char> row_cut_; std::vector<int> row_age_;
  std::vector<MasterCut> row_mc_;    // parallel to row_lo_ while keep_cuts (empty entries for non-cut rows)
  bool has_center_ = false; std::vector<double> cy_, cYp_;
  std::vector<double> gy_, gY_, Y_;   // scratch
  void add_rows(const std::vector<double>& L, const std::vector<double>& U, const std::vector<HighsInt>& starts,
                const std::vector<HighsInt>& idx, const std::vector<double>& val, bool is_cut);
  int age_rows();
  HighsModelStatus run_lp_escalating(const char* tag);
  void dump_stalled(const char* tag);
};
