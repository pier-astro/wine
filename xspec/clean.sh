#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

# Keep lmodel.dat, hand-written C++ sources, headers, scripts and documentation.
rm -f Makefile Makefile_redwine pkgIndex.tcl lpack_redwine.cxx redwine.cxx redwineFunctionMap.cxx redwineFunctionMap.h
rm -f -- ./*.o ./*.pco ./*.a ./*.so ./*.dylib
echo "Removed local XSPEC build products from $PWD."
