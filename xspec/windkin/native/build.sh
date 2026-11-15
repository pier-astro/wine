#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

if [[ -z "${HEADAS:-}" ]] || ! command -v initpackage >/dev/null || ! command -v hmake >/dev/null; then
    echo "Initialize HEASoft/XSPEC before building windkin." >&2
    exit 1
fi

# Select the config for the Python that can import NumPy and SciPy.
python_config="${WINDKIN_PYTHON_CONFIG:-python3-config}"
export WINDKIN_PY_CPPFLAGS="$("$python_config" --includes)"
export WINDKIN_PY_LDFLAGS="$("$python_config" --embed --ldflags)"

mkdir -p ../build/native
cp windkin.cxx lmodel.dat ../build/native/
cd ../build/native
rm -f lpack_windkin.cxx lpack_windkin.o windkinFunctionMap.cxx \
    windkinFunctionMap.h windkinFunctionMap.o windkin.o Makefile pkgIndex.tcl \
    libwindkin.dylib libwindkin.so
initpackage windkin lmodel.dat .
# XSPEC generates its Makefile. Extend its existing compiler and linker flags.
awk '/^include / {
    print "HD_CXXFLAGS += $(WINDKIN_PY_CPPFLAGS)"
    print "HD_SHLIB_LIBS += $(WINDKIN_PY_LDFLAGS)"
}
{print}' Makefile > Makefile.windkin
mv Makefile.windkin Makefile
hmake
echo "Load windkin with: lmod windkin $PWD"
