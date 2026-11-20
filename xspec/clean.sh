#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

# Keep lmodel.dat, hand-written C++ sources, headers, scripts and documentation.
rm -f Makefile Makefile_redwine pkgIndex.tcl lpack_redwine.cxx redwine.cxx redwineFunctionMap.cxx redwineFunctionMap.h
# Generated wrappers from the former package name also enter hmake sources.
rm -f lpack_wine.cxx wineFunctionMap.cxx wineFunctionMap.h
rm -f -- ./*.o ./*.pco ./*.a ./*.so ./*.dylib
echo "Removed local XSPEC build products from $PWD."
