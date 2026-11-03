# Experiment record

Every run behind a number in the paper, in order, with the CSV it wrote under `experiments/solver/`. All runs use the
C++ solver in `solver/` (HiGHS 1.15.1 built from source); the instance files are in `solver/instances/` and are
regenerated from their seeds by the solver's generator (`experiments/scripts/reproduce.sh instances`). Earlier
development history (an interpreted prototype, 2026-10-27 to 10-30) is not part of the record; its instance stream
is re-implemented in `solver/src/nprng.h`.

## Solver (solver/, CMake, HiGHS v1.15.1 built from source), 2026-10-30

`solver/src/{instance,recourse,blocks,decomp,bnc}.cpp`, binary `solver/build/endog`. Defaults: devex pricing, cut age
5, node tol 1e-5, in-out lam 0.2, incumbent cutoff, rounding heuristic, loose-then-tight root, solver reset before
every master solve. Recourse duals are solved warm in one HiGHS object per worker thread. Validation
(`tools/validate.sh 3 25 30 35 40`, then the 50-site instances): all 15 optima identical to
`instances/reference_optima.csv` to machine precision (rel diff <= 7e-16), all proven; node counts vary by a few
percent between runs because warm recourse solves return different degenerate dual vertices. Separation is 2-11 s
per run; master LP time is the remainder.

### MPI node-parallel tree (solver/src/mpi_main.cpp, Open MPI 5.0.11 from source), 2026-10-31

Rank 0 coordinates (heap, incumbent, global pool), ranks 1..R-1 bound nodes; root redundant on all ranks. First
runs with 1+4 ranks on a machine already loaded with three 100-site runs (6 threads busy of 8 cores):

| instance | single-rank C++ | 1+4 ranks | root (both) | tree speedup |
|---|---|---|---|---|
| n25 s1 | 21 s (91 nodes) | 15 s (89 nodes) | 9 s | 2.0x |
| n35 s2 | 299 s (527 nodes) | 148 s (523 nodes) | 35 s | 2.3x |
| n40 s3 | 168 s (97 nodes) | 125 s (125 nodes) | 38 s | 1.5x |

Same optima, all proven. Speedups are on an oversubscribed machine and understate the method; a clean measurement
needs an idle box. Node counts differ slightly because the best-first order interleaves across workers.

Clean measurement on the idle laptop (results/solver/mpi_4050.log, 2026-10-31 22:08): single rank with 8 separation
threads vs 1 coordinator + 7 workers.

| instance | single: nodes / root / total | 1+7: nodes / root / total | tree speedup | total speedup |
|---|---|---|---|---|
| n40 s1 | 21 / 19 s / 53 s | 19 / 40 s / 62 s | 1.5x | 0.85x |
| n40 s2 | 85 / 17 s / 89 s | 87 / 37 s / 82 s | 1.6x | 1.1x |
| n40 s3 | 131 / 12 s / 189 s | 127 / 32 s / 80 s | 3.7x | 2.4x |
| n50 s1 | 149 / 28 s / 456 s | 161 / 70 s / 190 s | 3.6x | 2.4x |
| n50 s2 | 61 / 23 s / 214 s | 71 / 45 s / 106 s | 3.1x | 2.0x |
| n50 s3 | 409 / 42 s / 970 s | 405 / 92 s / 362 s | 3.4x | 2.7x |

The redundant root is ~2x slower under MPI (8 copies of the master LP share memory bandwidth); tree speedup on 7
workers is 1.5-3.7x, poor on 19-87-node trees (best-first rarely has 7 comparable nodes open), 3.1-3.7x on
127-409-node trees. Net 0.85-2.7x. Single-rank idle-machine times are also the clean C++ numbers for the paper
(40: 53/89/189 s; 50: 456/214/970 s).

## Referee-report experiments (2026-10-31)

### Lazy McCormick rows and sparse lifting (C++ `--lazy-mc`, `--sparse`; results/solver/lazy_sparse.log)
n50 s2: default 277 s (87 nodes) vs --lazy-mc 226 s (65 nodes), same optimum. cgen n40 density 0.2: dense lifting
2 s / 2,581 rows, --sparse 5 s / 986 rows, root gap 0.22% vs 0.27% (sparse drops the factors of non-influencing
sites: a valid relaxation of the linked model). Same optimum 53.648998 in all variants.

### Sparse influence lifting at 100-500 sites (C++ `--sparse`; results/solver/sparse.log), 8 threads, wall 7200 s
Instances cgen_n{100,200,500}: m = 15, S = 20, influence density 0.1 / 0.05 / 0.02 (about 10 influencing sites per
scenario); the m = 100 variant has 100 customers. Separation is restricted to the influencing set I (a bug fix: with
non-influencing factors in the separation loop the root loop never converged).

| n | density | seed | OPT | root gap (time) | nodes | total | master rows |
|---|---|---|---|---|---|---|---|
| 100 | 0.10 | 1 | 47.829147 | 0.08% (1 s) | 255 | 8 s | 2,164 |
| 100 | 0.10 | 2 | 46.275180 | 0.04% (1 s) | 121 | 4 s | 2,184 |
| 100 | 0.10 | 3 | 50.156085 | 0.04% (1 s) | 43 | 6 s | 3,610 |
| 200 | 0.05 | 1 | 46.184683 | 0.01% (3 s) | 35 | 8 s | 6,112 |
| 200 | 0.05 | 2 | 46.601454 | 0.01% (3 s) | 2,155 | 114 s | 4,901 |
| 200 | 0.05 | 3 | 47.438137 | 0.01% (3 s) | 667 | 59 s | 4,448 |
| 500 | 0.02 | 1 | 47.766355 | 0.00% (12 s) | 57 | 26 s | 12,207 |
| 500 | 0.02 | 2 | 46.519291 | 0.00% (31 s) | 385 | 222 s | 20,922 |
| 500 | 0.02 | 3 | 45.413226 | 0.00% (17 s) | 167 | 82 s | 19,723 |
| 500, m=100 | 0.02 | 1 | 262.941205 | 0.00% (108 s) | 2,865 | 7201 s, gap 0.003% | 19,788 |

Dense lifting on the n100 seed 1 instance: root 28 s, 14,983 rows, gap 0.07% (vs 1 s, 2,164 rows, 0.08% sparse).
All proven optima are the same whether the lifting is sparse or dense (checked on n100 s1 and cgen n40). With
~10 influencing sites per scenario the master has O(S * |I| * n) rather than O(S * n^2) lifted columns, so 500 sites is
a 26-222 s job on a laptop. The m = 100 variant is not harder for the bound (root gap 0.00%) but each node costs
~1.5 s of recourse LPs and the tree is wide (2,865 nodes in 2 h, lb/ub 262.9340/262.9412).

### Sensitivity in S and endogeneity (C++ decomposed, cgen n=30, results/solver/sensitivity.csv)
| instance | OPT | root gap | nodes | time |
|---|---|---|---|---|
| S=50 s1 | 81.125149 | 0.26% | 51 | 65 s |
| S=50 s2 | 77.403381 | 0.26% | 83 | 56 s |
| S=50 s3 | 78.232978 | 0.24% | 35 | 41 s |
| S=100 s1 | 78.092443 | 0.25% | 159 | 196 s |
| S=100 s2 | 77.837251 | 0.18% | 11 | 35 s |
| S=100 s3 | 76.412461 | 0.30% | 37 | 65 s |
| weak S=20 s1 | 77.402408 | 0.18% | 141 | 22 s |
| weak S=20 s2 | 80.962734 | 0.00% | 0 | 0 s |
| weak S=20 s3 | 68.534553 | 0.19% | 123 | 43 s |
Root gap flat in S (0.18-0.30%); time grows about linearly in S (S=20 ~50 s, S=50 ~55 s, S=100 ~100 s mean).

Node-local McCormick on the same instances: see the C++ sens_mcb runs below.

## Referee report 2 (2026-10-31/11-01)

Code changes: `--master mcb | mcb-dyn` (src/mcb.cpp, McCormick masters in the same tree code); rounding heuristic
adds its dual vertices to the pool (bnc.cpp try_round); MPI root once on rank 0 (`--root-workers`,
DecomposedBound::export_state/import_state); node-local bounds Q_s(l)/Q_s(u) cached per node (S LPs per node).

### Decomposed C++ on new seeds 4-5 (results/solver/extra_seeds.log, idle machine, 8 threads)
| n | seed | OPT | root gap | nodes | time |
|---|---|---|---|---|---|
| 25 | 4 | 83.363274 | 0.37% | 21 | 12 s |
| 25 | 5 | 72.512428 | 0.55% | 129 | 26 s |
| 30 | 4 | 82.762817 | 0.25% | 125 | 82 s |
| 30 | 5 | 70.902920 | 0.34% | 109 | 39 s |
| 35 | 4 | 84.090611 | 0.44% | 311 | 133 s |
| 35 | 5 | 82.510571 | 0.45% | 67 | 80 s |
| 40 | 4 | 74.792125 | 0.24% | 117 | 163 s |
| 40 | 5 | 70.538514 | 0.19% | 317 | 452 s |
| 50 | 4 | 75.970369 | 0.12% | 113 | 517 s |
| 50 | 5 | 70.837130 | 0.08% | 205 | 448 s |
All proven; same ranges as seeds 1-3 (root gap 0.08-0.55%, 12-517 s).

### MPI with the root once on rank 0 (results/solver/mpi_4050b.log, 1+7 ranks, --root-workers 8)
| instance | single: root / total | 1+7: root / total | net |
|---|---|---|---|
| n40 s1 | 17 / 53 s | 28 / 50 s | 1.06x |
| n40 s2 | 18 / 94 s | 24 / 79 s | 1.2x |
| n40 s3 | 12 / 155 s | 19 / 64 s | 2.4x |
| n50 s1 | 35 / 394 s | 43 / 173 s | 2.3x |
| n50 s2 | 22 / 193 s | 34 / 96 s | 2.0x |
| n50 s3 | 38 / 706 s | 59 / 286 s | 2.5x |
The root is still ~1.5x slower than single-rank: the seven idle ranks busy-wait in MPI_Bcast (Open MPI default
progress). Rerun with `--mca mpi_yield_when_idle 1` in mpi_4050c.log.

### MPI with idle ranks yielding (results/solver/mpi_4050c.log, `--mca mpi_yield_when_idle 1`)
| instance | single root / tree / total | 1+7 root / tree / total | tree speedup | net |
|---|---|---|---|---|
| n40 s1 | 17 / 36 / 53 s | 16 / 17 / 33 s | 2.1x | 1.6x |
| n40 s2 | 18 / 76 / 94 s | 20 / 36 / 56 s | 2.1x | 1.7x |
| n40 s3 | 12 / 143 / 155 s | 12 / 37 / 49 s | 3.9x | 3.2x |
| n50 s1 | 35 / 359 / 394 s | 29 / 92 / 121 s | 3.9x | 3.3x |
| n50 s2 | 22 / 171 / 193 s | 24 / 54 / 78 s | 3.2x | 2.5x |
| n50 s3 | 38 / 668 / 706 s | 43 / 169 / 212 s | 4.0x | 3.3x |
Root now equal to single-rank; this is the version in paper Table 5.

### C++ McCormick masters (results/solver/mcb_cpp.log/.csv; `--master mcb-dyn`, `--master mcb`), idle machine
Node counts at 25-40 sites: 25: 1,211/485/425; 30: 2,191/663/1,531; 35: 2,479/3,471/1,827; 40: 5,675/-/8,303.
n40 s2 unproven in 7200 s (0.41%). Static MC-B has far fewer nodes than a tree without the heuristic cuts
(n25 s1 4,859) because the heuristic adds its dual vertices as Benders cuts for every master; both masters get the
same treatment. Table 2 of the paper is generated from these CSVs by solver/build/make_tables. Static seeds 1-3
(phase B), node-local seeds 4-5 (phase C) and the 50-site McCormick runs are queued in the same chain.

**Bug found 2026-11-01 19:50:** McbMaster::solve set HiGHS `time_limit` to the remaining wall, but HiGHS checks it
against its cumulative run clock, so every unproven C++ McCormick run stopped when LP time reached half the wall
(totals 3,640-3,660 s for wall 7,200; 7,306-7,320 s for wall 14,400). The decomposed solver already used
getRunTime() + cap. Fixed (same idiom); all unproven McCormick runs rerun in two parallel chains
(mcb_cpp_rerun_dyn.log: node-local n40 s2, s4, n50 s1-5; mcb_cpp_rerun_static.log: static n30 s1, n35 s1-2, n40 s1-3,
n50 s1-3). Two runs share the machine, so proven times from these reruns carry ~10% noise. The proven runs before the
fix are unaffected (they never hit a limit).


### 50-site McCormick in C++ (mcb_cpp_rerun_dyn.log, mcb_cpp_rerun_static.log; wall 14,400 s, after the time-limit fix)
| seed | static nodes / result | node-local nodes / result | decomposed nodes / time |
|---|---|---|---|
| 1 | 16,048 / unproven 3.39%, best sol. not optimal | 23,233 / proven 13,863 s | 125 / 394 s |
| 2 | 21,859 / unproven 2.50% | 27,386 / unproven 0.73% (clean rerun) | 63 / 193 s |
| 3 | 19,222 / unproven 3.23%, best sol. not optimal | 16,852 / unproven 1.32% | 281 / 706 s |
| 4 | not run | 19,995 / unproven 0.84% | 113 / 517 s |
| 5 | not run | 23,902 / unproven 0.40% | 205 / 448 s |
Node-local always holds the optimal solution; static misses it on 40 s1 and 50 s1/s3. Laptop hit zero battery
2026-11-02 ~04:51-07:31 (hibernate/resume): static n40 s3 overran to 9,602 s (unproven either way); the n50 s2
node-local run that spanned the gap was rerun clean (above). Table 2 of the paper is fully regenerated from the CSVs.

### Sensitivity McCormick columns in C++ (sens_mcb.log/.csv; `--master mcb-dyn --wall 3600`, idle machine, 2026-11-03)
| instance | root gap | nodes | time |
|---|---|---|---|
| S50 s1 | 7.89% | 1,239 | 164 s |
| S50 s2 | 8.77% | 2,535 | 608 s |
| S50 s3 | 8.04% | 1,599 | 239 s |
| S100 s1 | 8.41% | 4,270 | unproven 0.20% |
| S100 s2 | 9.22% | 3,189 | 2,235 s |
| S100 s3 | 7.49% | 2,111 | 1,351 s |
| weak s1 | 0.79% | 425 | 5 s |
| weak s2 | 0.00% | 1 | 0 s |
| weak s3 | 0.71% | 1,317 | 37 s |
These are the McCormick columns of paper Table 3. Nothing running as of 2026-11-03 morning.

## Unlinked roots and the HiGHS-MIP cut loop (2026-11-03)

Unlinked roots (`--master dec --unlinked --root-only`, tag `c2r-cpp`, experiments/solver/unlinked_cpp.csv): the 15 root
LP values at 25-50 sites, seeds 1-3; Table 1 is generated by solver/build/make_tables from mcb_cpp.csv (MC-B),
unlinked_cpp.csv and the dec runs.

HiGHS-MIP cut loop (`--master mcb-mip --wall 7200 --mip-time 1800`, tag `mcbmip-cpp`, experiments/solver/mcb_mip_cpp.csv):
static MC-B master with the converged root pool, y integer, HiGHS 1.15.1 MIP each round, integer point separated with
S recourse LPs, repeat. The loop proves all six:

| n | seed | rounds | MIP loop (s) | node-local MC-B tree (s) | decomposed (s) | MIP / decomposed |
|---|---|---|---|---|---|---|
| 30 | 1 | 7 | 180 | 4,499 | 24 | 7.5 |
| 30 | 2 | 3 | 36 | 326 | 22 | 1.6 |
| 30 | 3 | 3 | 30 | 337 | 10 | 3.0 |
| 40 | 1 | 3 | 234 | 763 | 53 | 4.4 |
| 40 | 2 | 6 | 1,264 | unproven (0.15%) | 94 | 13.4 |
| 40 | 3 | 5 | 584 | 1,185 | 155 | 3.8 |

Summed 2,328 s against 358 s. Reading: the MIP solver (its own cuts, presolve, heuristics) closes the McCormick root
gap 2-25x faster than our plain node-local MC-B tree, but is still 1.6-13x slower than the decomposed solver on every
instance; the per-round MIP solves at 40 sites take 75-310 s each (n40 s2: 75, 181, 193, 234, 271, 310 s), each a tree
from scratch, which a lazy-constraint callback would avoid. The paper's MIP paragraph, abstract and conclusion were
rewritten accordingly (an earlier "does not change this" claim is withdrawn). The gap
the MIP has to close is the unchanged 6.7-9.7% root gap of Table 1.

## Referee report 4 experiments (2026-11-03, scratchpad referee4.sh; CSVs bpos0.csv, seedc2l.csv, msens.csv)

**b_s = 0 variant** (`--bpos0`: the high-demand half of the scenarios gets base mass 0, so bpos_s = 0 and their linking
rows are vacuous; instances solver/instances/cgen_n{25..40}_bpos0_s{1,2,3}.txt; roots of mcb-dyn / dec --unlinked /
dec, dec trees with wall 3600 at 25-40, mcb-dyn trees at 25-30). Paper Table 4. Root gaps: MC-B 17-32%, unlinked
0.3-3.2%, linked within 0.1 pt of unlinked (as Remark rem:bpos predicts). Decomposed trees are much slower than on
the standard instances: 77-771 nodes, 21-2,834 s, n35 s3 unproven at 1 h (0.23%). Node-local MC-B at 25-30: faster on
4 of 6 (16-156 s vs 21-595 s), unproven on n30 s3 (9,738 nodes, 0.74%; dec 122 s). Honest reading in the paper: the
advantage is the root; at a 1-3% root gap the per-node master LP does not pay at small n.

**Seeded McCormick** (`--master mcb-dyn --seed-c2l`, the converged linked root pool merged into the MC-B pool, 25-35
sites seeds 1-3, wall 3600): root value identical to the unseeded MC-B root to the printed digits (the converged MC-B
root satisfies every Benders cut at its own point), nodes within 3% of Table 2's node-local column
(25: 1,231/489/421 vs 1,211/485/425; 30: 2,157/667/1,579 vs 2,191/663/1,531; 35: 2,409/3,469/1,893 vs 2,479/3,471/1,827),
times within seconds. One sentence in the MIP paragraph.

**Customers** (n = 30, m in {15, 50, 100, 200}, fresh C++-generated instances cgen_n30_m*_s*.txt, dec, 8 threads,
wall 3600). Paper Table 5.

| m | root gap | root (s) | pool at root | added in tree | nodes | time |
|---|---|---|---|---|---|---|
| 15 | 0.30 / 0.00 / 0.31% | 10 / 3 / 7 | 8.0k / 2.2k / 3.8k | 56 / 0 / 329 | 111 / 3 / 233 | 61 / 3 / 92 s |
| 50 | 0.50 / 0.57 / 0.78% | 24 / 58 / 58 | 16.0k / 19.6k / 16.8k | 1.0k / 0.5k / 2.1k | 985 / 131 / 397 | 325 / 142 / 282 s |
| 100 | 1.29 / 0.44 / 0.83% | 167 / 77 / 58 | 31.6k / 29.4k / 20.7k | 23.3k / 0.2k / 4.3k | 3,485 / 63 / 781 | 1,950 / 115 / 382 s |
| 200 | 1.24† / 1.09 / 1.52†% | 203 / 285 / 272 | 42.6k / 51.5k / 48.9k | 53.9k / 9.4k / 45.3k | 3,858 / 1,431 / 5,302 | unproven 0.18% / 1,113 s / unproven 0.02% |

(† against the incumbent.) Separation < 6% of the time at every m; on the unproven m = 200 runs 3,446 / 3,397 of
3,600 s are master LP solves (20-30 per node). The customer count reaches the solver through the root gap and the pool
size, not through the recourse LP size.

Paper now has 16 pages, 8 tables; every number from the C++ solver. Nothing running.
