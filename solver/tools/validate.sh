#!/bin/zsh
# Run the solver on the table instances and compare the optima with instances/reference_optima.csv.
# usage: tools/validate.sh [workers] [sizes...]   e.g. tools/validate.sh 3 25 30
cd "$(dirname "$0")/.." || exit 1
W=${1:-3}; shift; SIZES=${@:-25 30}
mkdir -p ../experiments/solver
for n in ${=SIZES}; do for s in 1 2 3; do
  echo "=== n$n s$s"
  ./build/endog --inst instances/n${n}_s${s}.txt --workers $W --out ../experiments/solver/validate.csv | tail -1
  ref=$(grep "^n${n}_s${s}," instances/reference_optima.csv)
  tail -1 ../experiments/solver/validate.csv | awk -F, -v ref="$ref" '{
    split(ref, r, ","); opt = $9 + 0; pr = $11; d = (opt - r[2]) / r[2]; if (d < 0) d = -d
    printf "    reference OPT %.6f nodes %s %ss | this run OPT %.6f proven %s | rel diff %.2e %s\n", r[2], r[3], r[4], opt, pr, d, (d < 1e-6 && pr == "1") ? "OK" : "MISMATCH" }'
done; done
