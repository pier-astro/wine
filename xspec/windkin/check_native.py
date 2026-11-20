"""Check the native uniform photon kernel against Python and analytic limits.

Requires a C++14 compiler and FFTW (the initialized HEASoft installation suffices).
Builds a small temporary C interface to the same header used by XSPEC. No XSPEC
interpreter is needed. Run: python xspec/windkin/check_native.py
"""
import ctypes
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile

import numpy as np
from scipy.integrate import quad

root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(root / "src"))
from redwine.experimental.kernels import Wind, convolve


def main():
    native = root / "xspec/windkin/native"
    headas = Path(os.environ["HEADAS"]) if "HEADAS" in os.environ else None
    wrapper = r'''
#include "uniform.h"
#include <string>
static thread_local std::string message;
extern "C" int run(int n,const double* edges,const double* flux,const double* params,double* out) {
    try {
        const std::vector<double> e(edges,edges+n+1),f(flux,flux+n);
        std::array<double,6> p; std::copy(params,params+6,p.begin());
        const auto result=windkin_uniform::convolve(e,f,p);
        std::copy(result.begin(),result.end(),out); return 0;
    } catch (const std::exception& error) {message=error.what(); return 1;}
}
extern "C" const char* error() {return message.c_str();}
'''
    with tempfile.TemporaryDirectory(prefix="windkin-native-") as directory:
        directory = Path(directory)
        source, library = directory / "check.cxx", directory / "check.so"
        source.write_text(wrapper)
        command = shlex.split(os.environ.get("CXX", "c++"))
        command += ["-std=c++14", "-O2", "-shared", "-fPIC", str(source), "-I" + str(native)]
        if headas:
            command += ["-I" + str(headas / "include"), "-L" + str(headas / "lib"),
                        "-Wl,-rpath," + str(headas / "lib")]
        command += ["-lfftw3", "-o", str(library)]
        subprocess.run(command, check=True)
        backend = ctypes.CDLL(str(library))
        pointer = np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")
        backend.run.argtypes = [ctypes.c_int, pointer, pointer, pointer, pointer]
        backend.run.restype = ctypes.c_int
        backend.error.restype = ctypes.c_char_p

        def evaluate(edges, flux, params):
            edges, flux, params = [np.ascontiguousarray(v, dtype=float) for v in (edges, flux, params)]
            output = np.zeros_like(flux)
            if backend.run(len(flux), edges, flux, params, output):
                raise ValueError(backend.error().decode())
            return output

        rng = np.random.default_rng(742)
        cases = [(.28, .17, 15, 75, 45, 0), (.05, .10, 0, 180, 0, 0),
                 (.25, 0, 0, 180, 90, 0), (.25, 0, 20, 50, 45, 0),
                 (.25, .2, 0, 90, 90, 0), (.25, .2, 0, 90, 0, 0),
                 (.99, .019999, 0, 180, 90, 0), (.9999, 0, 20, 50, 45, 0),
                 (.001, .002, 0, 90, 90, 0), (0, 0, 20, 60, 40, 0),
                 (.25, 1e-15, 20, 60, 40, 0), (.28, .17, 20, 70, 45, 1),
                 (-.28, .17, 20, 70, 45, 0), (.25, .1, 40, 40, 50, 0),
                 (.25, .1, 30, 60, 90, 0), (.25, .1, 90, 180, 180, 0)]
        for j in range(30):
            lo, hi = np.sort(rng.uniform(0, .9999, 2))
            inner, outer = np.sort(rng.uniform(0, 90, 2))
            cases.append(((lo + hi)/2, hi-lo, inner, outer, rng.uniform(0, 180), j % 2))
        grids = [np.geomspace(.001, 60, 201), np.geomspace(.001, 60, 40001),
                 np.array([0., .01, 2., 5., 6.7, 6.700001, 7., 8., 10., 100., 1e5]),
                 np.geomspace(6., 6.2, 4097)]
        maximum = 0.
        for edges in grids:
            flux = np.diff(edges) / np.maximum(edges[:-1], .001)**1.8
            flux /= flux.sum()
            flux[len(flux)//2] += 1
            for params in cases:
                wind = Wind("uniform", abs(params[0]), params[1], *params[2:5],
                            inflow=params[0] < 0, bicone=bool(params[5]))
                expected = convolve(edges, flux, wind)
                actual = evaluate(edges, flux, params)
                error = abs(actual - expected).sum()/max(abs(expected).sum(), 1e-100)
                maximum = max(maximum, error)
                assert error < 1e-5, (len(flux), params, error)
                np.testing.assert_allclose(actual, expected, rtol=1e-5, atol=1e-7)
                # A changed input spectrum must not return a cached output.
                np.testing.assert_allclose(evaluate(edges, -flux, params), -actual,
                                           rtol=1e-10, atol=1e-13)
            print(f"Native/Python: {len(flux)} bins, {len(cases)} geometries passed", flush=True)

        # A full spherical shell preserves photons when the entire Doppler
        # support lies on the grid. Independently integrate its analytic line.
        edges = np.geomspace(.001, 1e5, 40001)
        flux = np.zeros(len(edges)-1)
        source = np.searchsorted(edges, 6.7)-1
        flux[source] = 1.
        for center, width in ((.25, 0), (.05, .1), (.99, .019999)):
            actual = evaluate(edges, flux, (center, width, 0, 180, 90, 0))
            np.testing.assert_allclose(actual.sum(), 1., atol=2e-7, rtol=0)
        # Mono-speed full-sphere CDF, integrated over a finite source bin.
        edges = grids[2]
        flux = np.zeros(len(edges)-1)
        source = 4
        flux[source] = 1.
        beta = .25
        s = np.sqrt(1-beta*beta)
        bounds = (s/(1+beta), s/(1-beta))
        def exact_cdf(y):
            a, b = edges[source:source+2]
            return quad(lambda rest: s/(2*beta)*np.clip(y/rest-bounds[0], 0, bounds[1]-bounds[0]),
                        a, b, epsabs=1e-12)[0]/(b-a)
        expected = np.diff([exact_cdf(y) for y in edges])
        np.testing.assert_allclose(evaluate(edges, flux, (beta, 0, 0, 180, 90, 0)),
                                   expected, atol=1e-8, rtol=1e-7)
        # Invalid inputs are rejected before FFT/kernel construction.
        for bad in ((.1, .3, 0, 90, 45, 0), (.9, .3, 0, 90, 45, 0),
                    (.2, .1, 60, 30, 45, 0), (.2, .1, 0, 180, 45, 1)):
            try:
                evaluate(edges, flux, bad)
            except ValueError:
                pass
            else:
                raise AssertionError(f"Invalid parameters accepted: {bad}")
        print(f"Native uniform validation passed; maximum relative L1={maximum:.3g}.")


if __name__ == "__main__":
    main()
