# endog (C++)

Branch-and-cut for two-stage stochastic facility location with endogenous (decision-dependent, affine) scenario
probabilities, using the scenario-decomposed linked lifted-cut bound. Generator, masters, tree, MPI tree and the paper's table
builder are all here; nothing else is needed.

## Build

HiGHS is fetched from source at the pinned tag in `CMakeLists.txt` (v1.15.1, the version of every result in the paper) and
built statically into the tree. First configure downloads it; later builds are cached.

    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build -j

Use `-DHIGHS_SOURCE_DIR=/path/to/HiGHS` for a local checkout (offline builds, patches).

## Run

    ./build/endog --n 50 --S 20 --seed 1 --workers 4 --wall 14400 --out ../experiments/solver/n50.csv -v     # generate and solve
    ./build/endog --inst instances/n25_s1.txt --workers 3                                       # solve an instance file

Flags: `--master dec|mcb|mcb-dyn|mcb-mip` (decomposed linked master, static McCormick, node-local McCormick, all
sharing the tree code; static McCormick inside the HiGHS MIP solver), `--unlinked` (dec without linking rows),
`--mip-time sec` (per MIP solve, 1800), `--n --m --S --seed --endo strong|weak --strength --density --fmult --rng pcg64|mt`
(generator), `--gen-only` (write `--save-inst` and exit), `--workers T` (threads
for the recourse LPs), `--wall sec`, `--tol`, `--node-tol` (1e-5), `--root-loose` (1e-4), `--max-age` (5), `--lam`
(0.2), `--sparse` (sparse influence lifting), `--lazy-mc`, `--no-heur`, `--root-only`, `--out file.csv` (appends one
row), `--save-inst file`, `-v`.

The generator's default stream (`--rng pcg64`, `src/nprng.h`: SeedSequence + PCG64 XSL-RR, the stream of
numpy's `default_rng`) is such that `--n 25 --seed 1 --gen-only --save-inst f` reproduces `instances/n25_s1.txt` exactly;
`--rng mt` (std::mt19937_64) is the stream of the `cgen_*` instances. `tools/validate.sh [workers] [sizes...]`
solves the table instances and checks the optima against `instances/reference_optima.csv`; `tools/make_tables.cpp`
(`build/make_tables <repo root>`) writes Tables 1, 2, 4 and 5 of the paper from the result CSVs.

## Layout

- `src/instance.*`  instance data, text file format, generator; `src/nprng.h` the numpy-compatible random stream
- `src/recourse.*`  recourse dual LP (one warm-started HiGHS object per worker thread), cut pool, `F_exact`
- `src/blocks.*`    closed-form scenario blocks: levels and value/subgradient
- `src/decomp.*`    decomposed node bound: in-out stabilized multicut master over (y, Yp, theta) with cut aging
- `src/bnc.*`       root cut loop, best-first tree, threaded lifted-cut separation, rounding heuristic, CSV stats
- `src/mcb.*`       McCormick masters (static and node-local bounds) on the same tree
- `src/main.cpp`    CLI

## MPI node-parallel tree

Open MPI is built from its release tarball into `third_party/openmpi` by `third_party/get_openmpi.sh` (downloads the release tarball, checks its SHA-256, configures with internal pmix/hwloc/libevent, installs into `third_party/openmpi`).
CMake picks it up automatically and builds `endog_mpi`. Rank 0 solves the root once (`--root-workers T` threads),
broadcasts the cut pool, the active master cuts and the stabilisation centre, then coordinates (heap, incumbent,
global cut pool) while ranks 1..R-1 bound nodes. Same flags as `endog` (no `--root-only`, `--save-inst`). Pass
`--mca mpi_yield_when_idle 1` to `mpirun` so idle ranks do not spin during the root.

    third_party/openmpi/bin/mpirun --mca mpi_yield_when_idle 1 -np 8 ./build/endog_mpi --inst instances/n50_s3.txt --workers 1 --root-workers 8 --wall 14400 --out ../experiments/solver/mpi.csv -v

Use `--oversubscribe` when ranks exceed cores. `--workers` is threads per rank for the recourse LPs; with one rank per
core keep it at 1.
