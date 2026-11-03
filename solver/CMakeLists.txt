cmake_minimum_required(VERSION 3.24)
project(endog LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

# ---------------------------------------------------------------- HiGHS: built from source, pinned to v1.15.1 (the
# version of every result in the paper). Set HIGHS_SOURCE_DIR to use a local checkout instead.
set(HIGHS_GIT_TAG "v1.15.1" CACHE STRING "HiGHS git tag to fetch")
set(HIGHS_SOURCE_DIR "" CACHE PATH "Local HiGHS source checkout (skips the download)")
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
set(FAST_BUILD ON CACHE BOOL "" FORCE)        # HiGHS: library + executable only, no tests
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(PYTHON_BUILD_SETUP OFF CACHE BOOL "" FORCE)
include(FetchContent)
if(HIGHS_SOURCE_DIR)
  FetchContent_Declare(highs SOURCE_DIR "${HIGHS_SOURCE_DIR}")
else()
  FetchContent_Declare(highs
    GIT_REPOSITORY https://github.com/ERGO-Code/HiGHS.git
    GIT_TAG        ${HIGHS_GIT_TAG}
    GIT_SHALLOW    TRUE)
endif()
FetchContent_MakeAvailable(highs)

# ---------------------------------------------------------------- solver
find_package(Threads REQUIRED)
add_library(endog_core
  src/instance.cpp
  src/recourse.cpp
  src/blocks.cpp
  src/decomp.cpp
  src/bnc.cpp
  src/mcb.cpp)
target_include_directories(endog_core PUBLIC src)
target_link_libraries(endog_core PUBLIC highs::highs Threads::Threads)
target_compile_options(endog_core PRIVATE -Wall -Wextra -O3)

add_executable(endog src/main.cpp)
target_link_libraries(endog PRIVATE endog_core)
target_compile_options(endog PRIVATE -Wall -Wextra -O3)

# paper tables from the result CSVs (no solver dependency): solver/build/make_tables <repo root>
add_executable(make_tables tools/make_tables.cpp)
target_compile_options(make_tables PRIVATE -Wall -Wextra -O2)

# ---------------------------------------------------------------- MPI node-parallel tree (optional)
# Open MPI built from source into third_party/openmpi (see README); any MPI found by find_package works.
if(EXISTS "${CMAKE_SOURCE_DIR}/third_party/openmpi/bin/mpicc" AND NOT DEFINED MPI_HOME)
  set(MPI_HOME "${CMAKE_SOURCE_DIR}/third_party/openmpi")
  set(MPI_C_COMPILER "${MPI_HOME}/bin/mpicc")
  set(MPI_CXX_COMPILER "${MPI_HOME}/bin/mpicxx")
endif()
find_package(MPI COMPONENTS C)
if(MPI_C_FOUND)
  add_executable(endog_mpi src/mpi_main.cpp)
  target_link_libraries(endog_mpi PRIVATE endog_core MPI::MPI_C)
  target_compile_options(endog_mpi PRIVATE -Wall -Wextra -O3)
  message(STATUS "MPI found: endog_mpi will be built (${MPI_C_COMPILER})")
else()
  message(STATUS "MPI not found: endog_mpi skipped")
endif()
