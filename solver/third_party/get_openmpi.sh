#!/bin/sh
# Download, verify and build Open MPI into solver/third_party/openmpi (no system packages).
set -e
VER=5.0.11
SHA=e668a3c4acd50c41dc204c8a6dd98a611e0f26af89cf677577fa9be8a2698003
HERE="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$HERE/src" && cd "$HERE/src"
[ -f openmpi-$VER.tar.bz2 ] || curl -sL -O https://download.open-mpi.org/release/open-mpi/v5.0/openmpi-$VER.tar.bz2
echo "$SHA  openmpi-$VER.tar.bz2" | shasum -a 256 -c -
rm -rf openmpi-$VER && tar xjf openmpi-$VER.tar.bz2   # always a clean tree: the install prefix is baked into the build
cd openmpi-$VER
./configure --prefix="$HERE/openmpi" --disable-mpi-fortran --disable-oshmem --with-pmix=internal --with-hwloc=internal --with-libevent=internal
make -j"${JOBS:-4}"
make install
echo "Open MPI $VER installed in $HERE/openmpi"
