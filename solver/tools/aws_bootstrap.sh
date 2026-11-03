#!/bin/sh
# One-shot setup on a fresh Ubuntu 24.04 EC2 instance: toolchain, Open MPI from source, solver build.
# usage (on the server, after the repo is in ~/chandelure):   sh ~/chandelure/solver/tools/aws_bootstrap.sh
# Then run, inside tmux, with one rank per PHYSICAL core (half the vCPUs on c7i/c7a):
#   cd ~/chandelure/solver && third_party/openmpi/bin/mpirun -np $(nproc --all | awk '{print int($1/2)}') --bind-to core \
#     ./build/endog_mpi --inst instances/n100_s1.txt --workers 1 --wall 43200 -v --out ../experiments/solver/n100_aws.csv
set -e
sudo apt-get update -q && sudo apt-get install -y -q build-essential cmake git curl bzip2 tmux htop
cd "$(dirname "$0")/.."
JOBS=$(nproc) third_party/get_openmpi.sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j"$(nproc)"
echo "built: $(pwd)/build/endog and endog_mpi; physical cores: $(lscpu | awk '/^Core\(s\) per socket/{c=$4} /^Socket\(s\)/{s=$2} END{print c*s}')"
