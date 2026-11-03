# How the method was found

A record of the path from the question to the solver in `solver/` and the paper in `paper/main.tex`, with the
dead ends kept in. Dates are from the commit log (`git log`) and the experiment record (`experiments/RESULTS.md`).
Section numbers refer to the paper.

## 0. The question (2026-10-27)

Two-stage stochastic facility location in which the scenario probabilities depend on the opening decision:

    min_y  f'y + sum_s p_s(y) Q_s(y),   y in {0,1}^n,   p_s(y) = b_s + a_s'y,

Q_s(y) a transportation LP with capacities cap_i y_i and a penalty on unmet demand. The product p_s(y) Q_s(y) of
a decision-dependent weight and a decision-dependent LP value is bilinear, and the published exact method for it
is to introduce theta_s >= Q_s(y) and linearise y_j theta_s with McCormick envelopes over global bounds on Q_s.
The starting position was that this is weak because the bounds are global, and that other fields that multiply a
0/1 quantity by a value function might have a tighter device.

## 1. The cross-field survey (2026-10-27 to 10-28): `reframing-bilinearity-report.md` and `reframing-bilinearity/`

Seven parallel surveys (operations research, change of measure, variational physics, MDP/RL, global optimisation
and submodularity, economics and bilevel, experimental design) catalogued every reformulation device each field
uses for a decision-dependent expectation. The report's conclusion is in its first table: every published exact
method for this objective is one of five ways of making a 0/1 variable multiply a bounded continuous one, and
McCormick is already the convex hull of that single product. Fields that reach convexity do so by making the payoff
constant (choice models, persuasion), by making the distribution itself the decision (occupation measures, rational
inattention), or by making the energy equal the cost (annealing). None applies when the payoff is an LP value in y.

Four candidates came out of the survey as derived, not published:

1. a Gibbs/KL exponential-cone program (variational physics);
2. Benders cuts of Q_s multiplied by the probability and lifted to monomials (global optimisation);
3. Charnes-Cooper applied to the whole two-stage program (economics);
4. a Hammer-Rubin convexification (performative prediction / RL).

## 2. Verification on a toy (2026-10-28): `reframing-bilinearity/verification_of_derived_reformulations.md`

All four were proved correct and checked numerically on a 3-site, 4-customer instance with every cut loop run to
convergence, so that each number is the bound of the formulation and not of a truncated loop. Three lost badly to
McCormick (Candidate 1: 33-42% root gap, Candidate 3: 7%, Candidate 4: 40-50%). Candidate 2 in its "whole p times
cut" form neither dominated nor was dominated by McCormick. Two refinements of it decided the project:

- **2b, per-factor cuts.** Write p_s(y) as a posiform, bpos_s + sum_j a_sj y_j + sum_j |a_sj| (1 - y_j), and multiply
  each Benders cut of Q_s by each factor separately. Because the LP value is positively homogeneous in its
  right-hand side, y_j Q_s(y) and (1 - y_j) Q_s(y) are again recourse LPs at scaled right-hand sides, and the same
  dual vertex gives a valid cut for each of them. In the lifted space Y_jl = y_j y_l the cuts are linear, need no
  bound on Q_s, and dominate McCormick for every valid bound. Toy gaps 1.0-2.0% against 1.8-2.0%.
- **2c, the explicit level-1 RLT LP** (keep S(1+n) copies of the recourse) was integral on the toy but has S(1+n)
  recourse copies, which is not a route to large n.

These are Section 4 of the paper: the factor cuts (eq. cut0/cutp/cutm) and Remark 3 (the factor cuts are
Benders cuts at the conditional configurations).

## 3. The linking rows (2026-10-28)

The per-factor master has three epigraph variables per (s, j): e0_s for Q_s, e+_sj for y_j Q_s and e-_sj for
(1 - y_j) Q_s. They are not independent: y_j Q_s(y) + (1 - y_j) Q_s(y) = Q_s(y), so the row e0_s = e+_sj + e-_sj is
valid for every j. These "linking rows" are the level-1 RLT identity in epigraph form. On 5-10 sites they closed
about two thirds of the McCormick gap (average gaps: McCormick 2.7%, per-factor 2.2%, linked 0.9%, explicit RLT LP
0.26%). The linked lifted-cut master, C2L, was from then on the candidate contribution. What the rows buy, and
when they buy nothing, is Remark 4 of the paper: the tightening of scenario s is exactly bpos_s (e0*_s - L_0), so a
scenario whose whole base mass can be moved away gains nothing (this is the b_s = 0 experiment of Section 7, Table 4).

## 4. Kill criterion and the strong generator (2026-10-29)

The toy numbers were for n = 3, where good relaxations are nearly integral, so the question was whether the
advantage survives n. The criterion set in advance: if C2L kept less than half of the McCormick gap closed at
15-20 sites, stop. It kept 64-98% (gaps 0.01-1.23% against 0.73-3.44%). The generator was then changed from
"weak" endogeneity (uniform base probabilities, small shifts) to "strong" (the low-demand half of the scenarios
carries the mass at y = 0 and opened sites move up to 90% of it to the high-demand scenarios: facilities attract
demand), because weak instances are easy for every master (Table 3). On strong instances the McCormick root gap is
6-10% and the linked gap 0.1-0.4%; those are the instances of the paper.

## 5. The monolithic master hits a wall (2026-10-29)

The first branch-and-cut kept the whole linked LP in one solver object. Three devices made it reach 40 sites:
lazy lifted rows (only the e0 row of each cut is added, the e+/e- rows when violated), row aging (rows basic and
strictly slack for five solves are deleted), and an activation cap per scenario per pass. Two pitfalls are worth
recording. Aging that deletes rows merely because they are basic cycles forever on degenerate tight rows
(deleted, violated, re-activated, basic again); the fix is to age only rows that are basic and strictly slack, and
to make rows re-activated three times permanent. And a warm dual simplex occasionally cycles in-process on a
5,800-row LP that solves in 0.4 s from the same basis in a fresh solver state; the fix is to reset the solver and
set the basis before every master solve, which costs nothing. Even so, the node LP was 30-200 ms at 30 sites and up
to 9 s at 40 (33,000 rows), and McCormick's tree of 6 ms nodes won on two of six instances at 35-40 sites. The
tree sizes were right (17-617 nodes against 38,000-688,000); the per-node cost was not.

## 6. The decisive step: scenario decomposition with closed-form blocks (2026-10-29, commit a75ca06)

The linked LP is dual block-angular: for fixed (y, Y) the scenario blocks separate, and each block's value is a
maximum of finitely many affine functions of (y, Y) with a closed form, no LP. For scenario s the level of the e0
epigraph is L_0 = max_k cut_k(y), the factor levels L_j^+ and L_j^- are the same maxima over the lifted cuts, and the
linking rows force e0*_s = max(L_0, max_j (L_j^+ + L_j^-)); a subgradient in (y, Y) comes from the arg-max vertices
(Section 5, eq. levels/grad). The node bound is therefore a small multicut cutting-plane master over (y, Y, theta)
with S closed-form oracles. Stabilised in-out steps (lambda 0.2) with centres warm-started across nodes, and cut
aging in the master, made it 2.3x faster than the warm monolithic LP per node at 30 sites and 3.3x at 40 (root 8x),
with bounds equal to 4e-5. Because the cuts are valid at every node, the pool is shared by the whole tree.

With the tuned decomposed node bound (devex pricing, cut age 5, node tolerance 1e-5, incumbent cutoff in the
bound loop, rounding heuristic whose dual vertices join the pool) the branch-and-cut beat McCormick on 11 of 12
instances at 25-40 sites and proved all three 50-site instances in 6-37 min where McCormick was unproven after 4 h
(Section 7, Table 2). Profiling showed 97% of the time in the master LP solve and separation under 4%, which
settled the implementation plan: the master, not a compiled scenario loop, is the cost.

## 7. The solver (2026-10-30 to 10-31, commit 156f240)

The algorithm was then written in C++ (`solver/`): HiGHS 1.15.1 built from source, one warm-started recourse dual
per worker thread, the closed-form blocks in `blocks.cpp`, the decomposed bound in `decomp.cpp`, the tree in
`bnc.cpp`. Same optima to 1e-16 on every instance, 2.1x faster single-rank, separation negligible. Two extensions
followed: an MPI node-parallel tree (`mpi_main.cpp`, one coordinator, workers, root computed once on rank 0 and
broadcast; Section 8, Table 7) and a sparse-influence lifting for instances in which only a fraction of the sites
influence demand (`--sparse`; the lifting is restricted to (influencing) x (any) pairs, and the separation must be
restricted to the influencing set or the root loop never converges). The sparse variant proves 500-site instances
in minutes (Table 8); the dense master is quadratic in n and stops at 100 sites on a laptop (Table 6).

## 8. Four referee rounds (2026-10-31 to 11-03, commits c247866 to 02eff36)

The paper was reviewed four times and each round added an experiment rather than an argument:

- **A fair McCormick.** McCormick with node-local recourse bounds [Q_s(u), Q_s(l)] on the node box, recomputed
  with S LPs per node, has 10-40x fewer nodes than the static master and is time-competitive with ours up to 35
  sites. It became the baseline of Table 2; the static master is reported beside it. The growth rate is what it
  does not change: its root gap stays at 6-10%, and it is 6-14x slower at 40 sites and unproven at 50.
- **McCormick masters in the same C++ tree** (`mcb.cpp`), so that no comparison crosses implementations.
- **A full MIP solver.** The static McCormick master, seeded with its converged root pool, handed to the HiGHS MIP
  solver in a cut loop. The first version of this experiment, in the early code, left four of six instances
  unproven and the paper said "a full MIP solver does not change this". The C++ loop proves all six, 2-25x faster
  than our own McCormick tree and still 1.6-13x slower than the decomposed solver. The claim was withdrawn and the
  paragraph rewritten with the new numbers.
- **The identity behind the linking rows** was pushed to its limit (b_s = 0 on the high-demand scenarios, Table 4):
  linked and unlinked roots coincide, the factor cuts alone hold 0.3-3.2% against 17-32%, and on that gap the cheap
  McCormick node wins at 25-30 sites. The paper states this as the edge of the method.
- **Customers** (Table 5): the pool, not the recourse LP size, grows with m; at 200 customers two of three
  instances are unproven in an hour, with a root gap of 1.1-1.5%.
- Smaller items: the zero vertex in the dominance proposition, the Motzkin condition in the scope remark, the
  conditional-configuration reading of the factor cuts, a seeded-McCormick check (no measurable change).

## 9. What the project taught

- The lever was never a tighter envelope for y_j theta_s; it was removing the bound from the product by using
  positive homogeneity, so that the cut pool of one scenario serves every factor.
- The linking rows are the RLT level-1 identity in epigraph form, and their value is exactly bpos_s (e0* - L_0).
  Where that weight is zero the method degrades gracefully to the unlinked lifted cuts, which are still far tighter
  than McCormick, but the tree then has to work.
- A tight root is worth a slow node only up to a point: our node costs 0.07-0.20 s, McCormick's 6 ms. Above a
  root gap of about one percent the comparison at small n turns.
- Decomposition was the step that made the method practical, and it was found by looking at the dual structure of
  the linked LP rather than at the solver.
- Every claim that rested on a weaker implementation of the baseline eventually had to be withdrawn. The paper
  now compares three masters in one tree code and one MIP solver, on instances the generator reproduces from the seed.

## 10. Withdrawn or narrowed along the way

| claim | status |
|---|---|
| "RLT-1 explicit LP is the right master" (toy integrality) | narrowed: S(1+n) recourse copies; replaced by the linked lifted cuts |
| "handing McCormick to a full MIP solver does not change this" | withdrawn (Section 7): the MIP loop proves all six, 1.6-13x slower than ours |
| "a compiled scenario loop would help" | withdrawn: separation is under 5% of the time in the C++ solver |
| "the linking rows close two thirds of the gap" as a general statement | narrowed: exactly bpos_s (e0* - L_0), zero when bpos_s = 0 |
| McCormick with static bounds as the baseline | replaced by node-local bounds; both reported |
