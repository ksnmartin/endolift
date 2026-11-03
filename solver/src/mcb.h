// McCormick master (MC-B) branch-and-cut: the baseline of the paper, in the same code base as the decomposed solver.
// theta_s >= Q_s(y) by Benders cuts, products y_j theta_s replaced by z_sj in the McCormick envelope over
// [Qmin_s, Qmax_s]. node_local = false: static bounds Qmin = Q_s(1), Qmax = pen * sum_k d_sk;
// node_local = true: at every node theta_s in [Q_s(u), Q_s(l)] for the node box l <= y <= u and the envelope rows are
// rewritten with the node constants (one side inherited from the parent, so S recourse LPs per node).
// Tree: best-first on the parent bound, one fractional separation round per node, rounding heuristic, incumbent
// cutoff, branching on the most fractional coordinate (the same tree as bnc.cpp).
#pragma once
#include "bnc.h"

// init: an initial pool of dual vertices (e.g. the converged root pool of the linked master, --seed-c2l); may be null.
BncStats branch_and_cut_mcb(const Instance& inst, const BncOptions& opt, bool node_local, const CutPool* init = nullptr);

// MC-B inside a full MIP solver: the static McCormick master, seeded with the cut pool of its converged root LP, is
// handed to the HiGHS MIP solver with y integer; the integer point is separated with S recourse LPs, the new Benders
// rows are appended and the MIP re-solved (no lazy-constraint callback in HiGHS: every round is a fresh MIP solve,
// opt.mip_time seconds each). Proven when the MIP is optimal and its integer point has no violated cut.
BncStats mcb_mip_loop(const Instance& inst, const BncOptions& opt);
