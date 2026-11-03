# endolift

*endo*genous probabilities, *lift*ed cuts.

Two-stage stochastic facility location with **decision-dependent (endogenous) scenario probabilities**, solved to
proven optimality by a branch-and-cut whose relaxation is a *linked lifted-cut* master, evaluated by scenario
decomposition. The code base contains the C++ solver (generator, three masters, MPI tree, table builder), every
experiment of the paper, and the paper itself (`paper/main.tex`, author Karyamsetty Martin).

Headline: on 25 instances with 25–50 sites and 20 scenarios the solver proves optimality in 4 s to 12 min with a
root gap of 0.07–0.55 %, where the standard McCormick branch-and-cut (same tree code) has a 6–10 % root gap, leaves
13 of 40 runs unproven after 2–4 h, and needs 3–400× more nodes. With sparse influence the same solver proves
500-site instances in 26–222 s.

## 1. The problem

Sites `i = 1..n` with capacity `cap_i` and fixed cost `f_i` serve customers `k = 1..m`. Scenario `s = 1..S` has
demands `d_s` and transport costs `C`; unmet demand costs `pen` per unit. First stage `y ∈ {0,1}^n`; the recourse is
the transportation LP

    Q_s(y) = min { Σ C_ik x_ik + pen Σ u_k :  Σ_i x_ik + u_k ≥ d_sk,  Σ_k x_ik ≤ cap_i y_i,  x,u ≥ 0 }

with dual `Q_s(y) = max_{(α,β)∈D} αᵀd_s − Σ β_i cap_i y_i` over the **y-independent polytope**
`D = {α_k − β_i ≤ C_ik, 0 ≤ α ≤ pen, 0 ≤ β ≤ pen}` (the cap `β ≤ pen` is a valid truncation, paper §2).

The scenario probabilities are affine in the decision, `p_s(y) = b_s + a_sᵀy`, a probability vector at every binary
point. The problem is

    min_y  fᵀy + Σ_s p_s(y) Q_s(y)

whose objective is bilinear in `(p_s, Q_s)` and therefore non-convex: `y_j · Q_s(y)` is a binary variable times the
value of an LP. The textbook treatment (Basciftci, Ahmed & Shen 2021) replaces `y_j θ_s` by a McCormick envelope
over a static box `[Q_s(1), pen Σ d_s]`, which is why that relaxation is loose.

## 2. The algorithm

### 2.1 Posiform and factor recourse LPs

Write `p_s(y)` as a posiform: `p_s(y) = bpos_s + Σ_{a_sj>0} a_sj y_j + Σ_{a_sj<0} |a_sj| (1 − y_j)` with
`bpos_s = b_s + Σ_j min(a_sj, 0) ≥ 0`. Then `p_s(y) Q_s(y)` is a nonnegative combination of three kinds of
**factors**: `Q_s(y)`, `y_j Q_s(y)` and `(1 − y_j) Q_s(y)`.

Two facts carry everything. The dual polytope `D` does not depend on `y`, so each dual vertex gives a Benders cut
valid everywhere; and `Q` is positively homogeneous in its right-hand side. Hence, for binary `y`,

    y_j Q_s(y)       = Q(y_j d_s,        cap ∘ (y_j y))       (demand and capacities scaled by y_j)
    (1 − y_j) Q_s(y) = Q((1 − y_j) d_s,  cap ∘ (y − y_j y))

are again values of the *same* recourse LP, at right-hand sides that are affine in `(y, Y)` with `Y = yyᵀ` lifted
(`Y_jj = y_j`). Every vertex `(α, β)` of `D`, with `g = β ∘ cap`, therefore yields three cuts that are **linear in
the lifted variables**, one per factor epigraph `e⁰_s, e⁺_sj, e⁻_sj`:

    e⁰_s  ≥ αᵀd_s        − gᵀy
    e⁺_sj ≥ αᵀd_s · y_j  − gᵀY_j
    e⁻_sj ≥ αᵀd_s (1−y_j) − gᵀ(y − Y_j)

These are the Adams–Sherali (RLT level-1) multiplications of the Benders cut by `y_j` and `1 − y_j`, applied on the
dual side so that one dual vertex is reused by all `2n + 1` factors of a scenario and across scenarios.

### 2.2 Linking and the linked master (C2L)

For every `j`, `y_j Q_s(y) + (1 − y_j) Q_s(y) = Q_s(y)`. The master carries this as the **linking rows**
`e⁺_sj + e⁻_sj ≥ e⁰_s`, which tie the factor epigraphs of a scenario together. The linked master minimises
`fᵀy + Σ_s (bpos_s e⁰_s + Σ a_sj e⁺_sj + Σ |a_sj| e⁻_sj)` over the cuts, the linking rows and the McCormick rows of
`Y`. Paper Prop. 1: the relaxation is exact at every binary point once the pool holds the optimal vertices there;
Prop. 2 (dominance): `v_MC ≤ v_C2R ≤ v_C2L ≤ F` pointwise, where C2R is the unlinked per-factor model. On the test
instances the root gap drops from 6–10 % (McCormick) to 0.26–1.14 % (unlinked) to 0.07–0.43 % (linked). A remark in
the paper quantifies what linking buys: exactly `bpos_s · (e⁰*_s − L_0)` per scenario.

### 2.3 Scenario decomposition with closed-form blocks

With `(y, Y)` fixed the linked LP separates by scenario, and each scenario block has a closed form (Prop. 3). Define
the **levels**, maxima of affine functions of `(y, Y)` over the vertex pool:

    L_0   = max_v ( αᵀd_s − gᵀy )
    L⁺_j  = max_v ( αᵀd_s y_j − gᵀY_j )
    L⁻_j  = max_v ( αᵀd_s (1−y_j) − gᵀ(y − Y_j) )

Then the block value is `φ_s(y,Y) = bpos_s · e⁰*_s + Σ_{a>0} a_sj L⁺_j + Σ_{a<0} |a_sj| L⁻_j` with
`e⁰*_s = max(L_0, max_j (L⁺_j + L⁻_j))`, and a subgradient `γ_s` is read off the active pieces (paper eq. (16)–(17);
any element of the subdifferential is valid at ties). Evaluating all `S` blocks costs one matrix product
`G_s Yᵀ` per scenario and **no LP**.

The node bound is therefore computed by a **multicut cutting-plane master** over `(y, Y, θ)` with `θ_s ≥ φ_s`
cuts `θ_s ≥ φ_s(ŷ,Ŷ) + γ_sᵀ((y,Y) − (ŷ,Ŷ))`, stabilised in-out (`--lam 0.2`), with cut aging (`--max-age 5`),
a loose root tolerance (`--root-loose 1e-4`) and a node tolerance (`--node-tol 1e-5`). Its cuts are valid at every
node of the tree because they are valid for every `(y, Y)`, so the pool and the master rows are inherited by the
children; the node bound recorded is the master LP value.

### 2.4 Branch-and-cut

    root:  repeat { solve master; separate: S recourse LPs at ŷ (and at the factor points for the influencing set);
                    add new dual vertices to the pool } until no violated cut   → root gap 0.07–0.55 %
    tree:  best-first on the parent bound; at a node: inherit cuts, re-solve the master, one separation round,
           rounding heuristic (its dual vertices also enter the pool), incumbent cutoff,
           branch on the most fractional y_j.  Proven when the heap is empty.

Separation solves one warm-started dual recourse LP per scenario in a thread pool (`--workers`). The same tree code
runs three masters (`--master`):

| master | relaxation | per-node cost |
|---|---|---|
| `dec` | linked lifted-cut master via scenario decomposition (the paper's method) | master LP + `S` recourse LPs |
| `mcb` | McCormick envelope with static bounds `[Q_s(1), pen Σ d_s]` + Benders cuts (Basciftci et al.) | master LP + `S` recourse LPs |
| `mcb-dyn` | McCormick with **node-local** bounds `[Q_s(u), Q_s(l)]` for the node box `l ≤ y ≤ u`, one side inherited from the parent | master LP + `2S` recourse LPs |
| `mcb-mip` | static McCormick master handed to the HiGHS MIP solver each round, integer point separated, cuts appended (the "full MIP solver" baseline) | one MIP solve + `S` recourse LPs per round |

`--unlinked` with `--master dec` drops the linking rows (the C2R ablation of Table 1).

### 2.5 Extensions in the code

- **Node-parallel MPI tree** (`endog_mpi`): rank 0 solves the root once (`--root-workers`), broadcasts the pool, the
  active master cuts and the stabilisation centre, then coordinates heap, incumbent and global pool while ranks
  `1..R−1` bound nodes. Tree speedup 2.1–4.0× on 7 workers, net 1.6–3.3×.
- **Sparse influence lifting** (`--sparse`): when only a set `I` of sites influences the probabilities
  (`a_sj ≠ 0`), the lifting keeps pairs `(j ∈ I) × (any)` only, so the master has `O(S·|I|·n)` instead of
  `O(S·n²)` lifted columns; separation must be restricted to `I`. This is what reaches 500 sites.
- **Lazy McCormick rows** (`--lazy-mc`): the McCormick rows of `Y` added only when violated.
- **Instance generator** with weak/strong endogeneity, influence density, scenario count (`--endo --strength
  --density --S`).

## 3. Repository layout

    chandelure/
    ├── README.md                      this file
    ├── paper/
    │   └── main.tex, main.pdf         the paper (tectonic); Tables 1, 2, 4 and 5 are written between % TABk markers by make_tables
    ├── solver/                        C++ solver (CMake, C++17)
    │   ├── CMakeLists.txt             fetches and builds HiGHS v1.15.1 from source (FetchContent)
    │   ├── README.md                  build and CLI details
    │   ├── src/
    │   │   ├── instance.*             data, text format, generator (PCG64 stream, nprng.h)
    │   │   ├── nprng.h                SeedSequence + PCG64 (XSL-RR) stream: seeds reproduce the table instances bit for bit
    │   │   ├── recourse.*             dual recourse LP per worker thread, cut pool, exact F(y)
    │   │   ├── blocks.*               closed-form scenario blocks: levels, value, subgradient (§2.3)
    │   │   ├── decomp.*               decomposed node bound: in-out stabilised multicut master, aging, export/import
    │   │   ├── bnc.*                  root loop, best-first tree, threaded separation, heuristic, CSV stats (dec)
    │   │   ├── mcb.*                  McCormick masters, static and node-local, same tree (mcb, mcb-dyn)
    │   │   ├── main.cpp               CLI  → build/endog
    │   │   └── mpi_main.cpp           MPI coordinator/worker → build/endog_mpi
    │   ├── instances/                 n{25..50}_s{1..5}.txt table instances (`--n N --seed s`), cgen_*.txt sensitivity,
    │   │                              sparse, bpos0 and customer instances (`--rng mt`), reference_optima.csv
    │   ├── tools/                     make_tables.cpp (→ build/make_tables), validate.sh, aws_bootstrap.sh
    │   └── third_party/get_openmpi.sh builds Open MPI 5.0.11 from its tarball into third_party/openmpi
    ├── experiments/
    │   ├── RESULTS.md                 the running record of every run, with the numbers quoted in the paper
    │   ├── scripts/reproduce.sh       reruns every C++ experiment stage by stage (§6)
    │   └── solver/*.csv               results, one row per run (logs are not tracked)
    └── notes/                         how the method was found, the cross-field survey and its verification, literature sizes

## 4. Build

Everything is built from source into the repository tree. The only things installed on the machine are a C++
compiler, CMake and git; HiGHS and Open MPI are downloaded and compiled by the build itself, so no package-manager
binaries enter the solver stack.

### 4.1 What you need

| requirement | why | tested with |
|---|---|---|
| C++17 compiler with 128-bit integers (`unsigned __int128`: clang or gcc) | the solver; the PCG64 stream in `nprng.h` | Apple clang 21 (macOS 26) |
| CMake ≥ 3.24 | the build; `FetchContent` for HiGHS | CMake 4.4 |
| git and network access on the first configure | fetches HiGHS v1.15.1 (~50 MB) | — |
| POSIX shell, `curl`, `shasum`, `tar`, `make` (optional) | `get_openmpi.sh` downloads, checks and builds Open MPI 5.0.11 | macOS built-ins |
| `tectonic` (optional) | compiles `paper/main.tex`; it fetches the LaTeX packages it needs on first run | Tectonic 0.17 |
| 8 cores, 16 GB | the paper's timings; the solver runs on fewer cores, `--workers` sets the thread count | Apple silicon laptop |

Disk: about 90 MB for the solver build tree and 120 MB for Open MPI. No Python, no system LP solver, no license.

### 4.2 Solver

    cd solver
    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release     # first run: downloads HiGHS v1.15.1 and builds it statically, ~3 min
    cmake --build build -j

This produces `build/endog` (the solver and generator) and `build/make_tables` (the paper's table builder). Check it:

    ./build/endog --inst instances/n25_s2.txt --workers 4           # proves OPT 79.111870 in a few seconds
    tools/validate.sh 4 25                                            # the three 25-site instances against reference_optima.csv
    ./build/endog --n 25 --seed 2 --gen-only --save-inst /tmp/i.txt   # regenerates instances/n25_s2.txt from the seed

`-DHIGHS_SOURCE_DIR=/path/to/HiGHS` builds from a local HiGHS checkout (offline machines); `-DHIGHS_GIT_TAG` pins
another tag, but the paper's numbers are with v1.15.1 and a different HiGHS can change node counts.

### 4.3 MPI tree (optional)

The node-parallel tree needs an MPI. The script builds Open MPI 5.0.11 from its tarball (SHA-256 checked) into
`solver/third_party/openmpi`, without Fortran or OpenSHMEM, with its internal PMIx, hwloc and libevent:

    third_party/get_openmpi.sh                           # ~10 min; JOBS=8 to use more cores
    cmake -S . -B build && cmake --build build -j        # CMake now finds mpicc there and also builds build/endog_mpi

Any MPI that `find_package(MPI)` can see works instead (set `MPI_HOME`); the paper used this Open MPI. The install
prefix is compiled into the Open MPI binaries, so if the repository is moved the script has to be run again. Run the
MPI tree with `--mca mpi_yield_when_idle 1` as in §5, otherwise idle ranks slow the shared root computation.

### 4.4 Paper

    cd paper && tectonic main.tex                        # main.pdf; first run downloads the LaTeX packages

The tables between the `% TABk-BEGIN/END` markers are written by `solver/build/make_tables <repo root>` from the
CSVs in `experiments/solver/` (`reproduce.sh tables` does both steps).

### 4.5 Linux notes

The code is portable C++17 and the scripts are POSIX shell or zsh (`validate.sh` and `reproduce.sh` use zsh).
On Linux install `g++` or `clang`, `cmake`, `git`, `curl`, `bzip2`, `make` and, for the two scripts, `zsh`;
`solver/tools/aws_bootstrap.sh` installs the toolchain on a fresh Ubuntu 24.04 host (add `zsh` by hand) and then
runs the two builds above. `caffeinate` in the reproduction notes is
macOS-only; use `systemd-inhibit` or nothing.

### 4.6 If something fails

- *"Could not find HiGHS" or a download error on the first configure*: the fetch needs network access or
  `-DHIGHS_SOURCE_DIR`; delete `build/` and configure again after fixing it.
- *CMake complains the cache was created in another directory*: the tree was moved or copied; delete `build/`
  and configure again (and rerun `get_openmpi.sh` if the MPI tree is used).
- *`mpirun` fails with "Library not loaded"*: same cause, the Open MPI prefix changed; rerun `get_openmpi.sh`.
- *HiGHS "threads" or "time_limit" warnings*: HiGHS prints them when a run is restarted in-process; they are harmless
  and the solver sets its own limits.

## 5. Running the solver

All commands run from `solver/`. An instance comes either from a file (`--inst`) or from the generator (`--n` and
friends); the master is chosen with `--master`; the run ends with one summary line and, with `--out`, one CSV row.

    ./build/endog --inst instances/n25_s1.txt --workers 8                        # decomposed linked master (default), 8 threads
    ./build/endog --inst instances/n40_s2.txt --master mcb-dyn --workers 8 --wall 7200 --out ../experiments/solver/x.csv
    ./build/endog --inst instances/n30_s1.txt --master mcb-mip --workers 8 --wall 7200 --mip-time 1800
    ./build/endog --n 100 --density 0.1 --seed 1 --rng mt --sparse --workers 8 -v   # generate a sparse instance and solve it
    ./build/endog --n 30 --seed 4 --gen-only --save-inst instances/mine.txt          # generate only (this is instances/n30_s4.txt)
    third_party/openmpi/bin/mpirun --mca mpi_yield_when_idle 1 -np 8 ./build/endog_mpi \
        --inst instances/n50_s3.txt --workers 1 --root-workers 8                    # 1 coordinator + 7 workers, root once on rank 0

| flag | meaning | default |
|---|---|---|
| `--inst file` | load an instance file; overrides the generator flags | — |
| `--n N --m M --S S --seed K` | generate: sites, customers, scenarios, seed | 15, 15, 20, 1 |
| `--endo strong\|weak`, `--strength x`, `--density x`, `--fmult x` | generator: endogeneity model, mass a site can move, fraction of influencing sites, fixed-cost multiplier | strong, 0.9, 1.0, 1.0 |
| `--rng pcg64\|mt` | generator stream: `pcg64` reproduces the table instances `n*_s*.txt` (names `gen_*`), `mt` the `cgen_*` files | pcg64 |
| `--bpos0` | generator, strong endogeneity only: the high-demand scenarios get base mass 0, so their linking rows are vacuous (Table 4) | off |
| `--gen-only` | write the instance with `--save-inst` and exit without solving | off |
| `--master dec\|mcb\|mcb-dyn\|mcb-mip` | decomposed linked master; static McCormick; node-local McCormick; static McCormick inside the HiGHS MIP solver (cut loop) | dec |
| `--unlinked` | `dec` without the linking rows (the C2R model of Table 1) | off |
| `--seed-c2l` | `mcb`/`mcb-dyn` only: start the McCormick master from the converged linked root pool | off |
| `--mip-time sec` | `mcb-mip` only: HiGHS time limit per MIP solve | 1800 |
| `--workers T` | threads for the recourse LPs (separation, node-local bounds, heuristic) | 1 |
| `--wall sec` | wall-clock limit; the run reports `NOT PROVEN` with the bound reached | none |
| `--tol x` | optimality tolerance on the absolute gap | 1e-6 |
| `--node-tol x`, `--root-loose x` | stopping tolerance of the node cutting-plane loop; loose tolerance of the first root phase (0 disables) | 1e-5, 1e-4 |
| `--max-age k`, `--lam x` | delete master rows slack for k solves; in-out stabilisation weight | 5, 0.2 |
| `--sparse` | `dec`: lift only pairs with an influencing site; needed for density < 1 at 100+ sites | off |
| `--lazy-mc` | `dec`: add the McCormick rows on the lifted variables only when violated | off |
| `--no-heur`, `--root-only` | no rounding heuristic; stop after the root cutting-plane loop | off |
| `--out file.csv`, `--save-inst file`, `-v` | append one result row (header written if the file is new); save the instance; print every round and node | — |
| `--root-workers T` | `endog_mpi` only: threads rank 0 uses for the root before broadcasting the pool | 1 |

`endog_mpi` runs the decomposed master only and takes the instance, generator, tolerance and `--sparse`/`--lazy-mc`
flags above; it needs at least two ranks (one coordinator). The `--mca mpi_yield_when_idle 1` option matters on a
laptop: without it the idle worker ranks slow the root computation by half.

The last output line summarises the run and starts with the master's tag (`dec-cpp`, `c2r-cpp` for `--unlinked`,
`mcb-cpp`, `mcbdyn-cpp`, `mcbmip-cpp`):

    dec-cpp: OPT 82.708614 lb 82.708614 proven  open 42/50  root gap 0.10% (12 it, 35s)  nodes 125  lps 981  sep LPs 26400  lp 341s sep 48s  total 394s  rows 3961 (max 4031, deleted 20714)

`OPT` is the incumbent, `lb` the proven bound, `open` the opened sites, `root gap` is (OPT − root LP value)/OPT with
the root loop's iterations and time, `lps` the master LP solves, `sep LPs` the recourse LPs, `lp`/`sep`/`total` the
times. The CSV row (`--out`) has the columns

    name, n, m, S, endo, master, ranks, workers, OPT, lb, proven, n_open, root_value, root_gap, root_iters, root_time,
    root_cuts, nodes, lps, sep_lps, lp_time, sep_time, total_time, cuts, rows, vars, max_rows, rows_deleted,
    node_tol, heur, heur_hits, pruned_early, root_loose, max_age, lam, y

where `root_cuts`/`cuts` are the pool sizes after the root and at the end, and `y` the incumbent as a 0/1 string.
`tools/validate.sh 8 25 30 35` solves the table instances of those sizes and checks the optima against
`instances/reference_optima.csv`.

Instance file format (`endog-instance 1`): line 1 the format tag, then the name, the endogeneity model, `n m S pen`,
then `cap` (one row of n), `f` (one row of n), `C` (n rows of m), `d` (S rows of m), `b` (one row of S), `A` (S rows
of n). The generator writes it with `--save-inst`; `--inst` reads it.

## 6. Reproducing the paper's experiments

All runs were made on an idle 8-core, 16 GB laptop, 8 threads, walls 7,200 s at 25–40 sites and 14,400 s at 50.
`experiments/scripts/reproduce.sh` runs each stage with exactly these settings and writes fresh CSVs under
`experiments/solver/`; `experiments/RESULTS.md` records the runs behind the committed CSVs, with dates.

| paper table | stage | what it runs | time |
|---|---|---|---|
| instances | `reproduce.sh instances` | regenerates `n{25..50}_s{1..5}.txt` with `--n N --seed s --gen-only` and checks them against the committed files | 1 min |
| Table 2, decomposed columns | `reproduce.sh decomposed` | `--master dec`, 25 instances | ~1 h |
| Table 2, McCormick columns | `reproduce.sh mccormick` | `--master mcb-dyn` on seeds 1–5, `--master mcb` on seeds 1–3 | 3–4 days (most runs hit the wall) |
| Table 7, MPI | `reproduce.sh mpi` | single rank vs `mpirun -np 8 … --root-workers 8` at 40 and 50 sites | ~1 h |
| Table 3, sensitivity | `reproduce.sh sensitivity` | generates `cgen_n30_{S50,S100,weak}_s{1,2,3}`, runs `dec` and `mcb-dyn` with wall 3,600 s | ~4 h |
| Table 8, sparse influence | `reproduce.sh sparse` | generates `cgen_n{100,200,500}` at density 0.1/0.05/0.02, runs `dec --sparse` | ~3 h |
| Table 1, unlinked column | `reproduce.sh unlinked` | `--master dec --unlinked --root-only` on seeds 1–3 | 10 min |
| "McCormick inside a full MIP solver" | `reproduce.sh mip` | `--master mcb-mip --wall 7200 --mip-time 1800` on the 30- and 40-site seeds 1–3 | ~40 min |
| Table 4, b_s = 0 variant | `reproduce.sh bpos0` | generates `cgen_n{25..40}_bpos0_s{1,2,3}`, roots of all three masters, `dec` trees, `mcb-dyn` trees at 25–30 | ~3 h |
| seeded McCormick sentence | `reproduce.sh seedc2l` | `--master mcb-dyn --seed-c2l` on the 25–35 seeds 1–3 | ~2 h |
| Table 5, customers | `reproduce.sh msens` | generates `cgen_n30_m{15,50,100,200}_s{1,2,3}`, `dec` with wall 3,600 s | ~2 h |
| Tables 1, 2, 4, 5 → `main.tex` | `reproduce.sh tables` | `solver/build/make_tables` then `tectonic` | 1 min |

Run `caffeinate -i -s experiments/scripts/reproduce.sh all` on macOS to keep the machine awake; the stages are
independent, so run them one at a time if the machine is shared. `make_tables` reads the committed CSV names
(`mcb_cpp.csv`, `extra_seeds.csv`, …); to build Table 2 from a fresh rerun, point its `load(...)` lines at the new
files (`decomposed.csv`, `mccormick.csv`). Expect ~10 % timing noise between runs and small node-count differences
when the machine is not idle; optima and root gaps are deterministic.

Every number in the paper comes from this code, and every instance is reproduced by it from the seed: the table
instances by the default PCG64 stream of `nprng.h`, the `cgen_*` instances by `--rng mt`.
The dense 100-site runs (`experiments/solver/n100_s*.log`, laptop limit, not in a table) are
`./build/endog --inst instances/n100_s1.txt --workers 8 --wall 14400 -v` and the MPI variant with `-np 8`.

## 7. Results at a glance (from the paper)

| | root gap | 25–50 sites, 25 instances | 500 sites, sparse |
|---|---|---|---|
| decomposed linked master | 0.07–0.55 % | all proven, 4 s – 12 min, 9–431 nodes | proven, 26–222 s |
| McCormick, node-local bounds | 6–10 % | 19 of 25 proven; 700 → 22,000 nodes from 25 to 50 sites | — |
| McCormick, static bounds | 6–10 % | 8 of 15 proven; best solution not optimal on 3 of the 6 unproven 40/50-site runs | — |

See `experiments/RESULTS.md` for every run and `paper/main.pdf` for the full account.

## 8. AI assistance

This project was carried out with Claude (Anthropic) as an AI assistant: the literature survey in `notes/`, the C++
solver, the experiment scripts, the experiment record and the paper were written with it, in sessions directed by
the author. The research question, the method, the experimental standards and the final decisions on every claim
are the author's, who reviewed and edited all generated content, verified the proofs and the reported numbers
against the solver's output, and takes full responsibility for the code and the paper. The same disclosure appears in
the paper (section "Declaration on the use of AI assistance").

## 9. Citation and license

If you use the code or the results, cite the paper and the software (`CITATION.cff` has both; GitHub-style
"Cite this repository" tools read it):

    @unpublished{martin2026endolift,
      author = {Martin, Karyamsetty},
      title  = {Facility Location with Decision-Dependent Scenario Probabilities:
                Linked Lifted Cuts and a Scenario-Decomposed Branch-and-Cut},
      year   = {2026},
      note   = {Manuscript. Code, instances and results: the endolift repository}
    }

    @software{martin2026endolift_code,
      author  = {Martin, Karyamsetty},
      title   = {endolift: a scenario-decomposed branch-and-cut with linked lifted cuts
                 for facility location with decision-dependent scenario probabilities},
      year    = {2026},
      version = {1.0.0}
    }

The code, scripts, instances, results and notes are released under the MIT License (`LICENSE`). The manuscript in
`paper/` is copyright the author, all rights reserved until its publication terms are settled. HiGHS (MIT) and
Open MPI (BSD-3-Clause) are built from source by the build and keep their own licenses.
