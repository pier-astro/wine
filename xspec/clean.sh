#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

# Keep lmodel.dat, hand-written C++ sources, headers, scripts and documentation.
rm -f Makefile Makefile_wine pkgIndex.tcl lpack_wine.cxx wine.cxx wineFunctionMap.cxx wineFunctionMap.h
rm -f -- ./*.o ./*.pco ./*.a ./*.so ./*.dylib
echo "Removed local XSPEC build products from $PWD."
