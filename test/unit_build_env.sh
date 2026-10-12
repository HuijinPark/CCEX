#!/usr/bin/env bash
# How to compile a unit test against the current build, for the drivers that compile one outside the build
# (check_rotation, check_pcce/unit_build.py; D26). Prints KEY='value' lines:
#   CCEX_UNIT_CXX       the compiler (the MPI wrapper for do_compile.sh, CMake's own compiler for CMake)
#   CCEX_UNIT_ENV       NAME=value pairs that compiler needs in its environment (I_MPI_CXX=icpx on idm)
#   CCEX_UNIT_CXXFLAGS  the flags of bin/main.out
#   CCEX_UNIT_INCLUDES  the -I flags (with MPI's for CMake)
#   CCEX_UNIT_LIBS      the libraries
#   CCEX_UNIT_OBJS      the objects of main.out but main's, with the refinement policies' and SQLite's
# The build writes them itself: the Makefile of do_compile.sh through its unit-env target, a CMake build
# (CCEX_BUILD_DIR=<its build folder>, any generator) as unit-build.env, or unit-build-<config>.env for a
# multi-config generator (CCEX_BUILD_CONFIG=<config>). The build is brought up to date first
# (make bin/main.out, cmake --build --target ccex), so the objects are those of the current sources; until
# 2026-10-02 each driver compared objects with sources and headers itself and refused stale ones.
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
if [ -n "${CCEX_BUILD_DIR:-}" ]; then
    # a multi-config generator (Ninja Multi-Config, Xcode, Visual Studio) writes one file per configuration:
    # CCEX_BUILD_CONFIG=<config> picks it and the configuration that is built
    if [ -n "${CCEX_BUILD_CONFIG:-}" ]; then
        config=(--config "$CCEX_BUILD_CONFIG"); env_file="$CCEX_BUILD_DIR/unit-build-$CCEX_BUILD_CONFIG.env"
    else
        config=(); env_file="$CCEX_BUILD_DIR/unit-build.env"
    fi
    cmake --build "$CCEX_BUILD_DIR" "${config[@]}" --target ccex >&2 \
        || { echo "ERROR: cmake --build $CCEX_BUILD_DIR ${config[*]} --target ccex failed" >&2; exit 1; }
    if [ ! -f "$env_file" ]; then
        echo "ERROR: no $env_file: configure again with this CMakeLists.txt (cmake -B ...)" >&2
        ls "$CCEX_BUILD_DIR"/unit-build-*.env > /dev/null 2>&1 && echo "       a multi-config build: set CCEX_BUILD_CONFIG to one of: $(cd "$CCEX_BUILD_DIR" && ls unit-build-*.env | sed 's/^unit-build-//; s/\.env$//' | tr '\n' ' ')" >&2
        exit 1
    fi
    cat "$env_file"
else
    [ -f "$REPO/Makefile" ] || { echo "ERROR: no Makefile: run \"bash do_compile.sh <preset>\" first, or set CCEX_BUILD_DIR to a CMake build" >&2; exit 1; }
    out="$(make -s -C "$REPO" --no-print-directory unit-env)" \
        || { echo "ERROR: make unit-env failed (a Makefile from before 2026-10-02 has no unit-env: run do_compile.sh again)" >&2; exit 1; }
    grep "^CCEX_UNIT_" <<< "$out"
fi
