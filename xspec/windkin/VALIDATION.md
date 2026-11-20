# Validation

## Photon convention

The shared Wind and Sector kernels use RED-WINE's integrated photon weight
delta squared, following equation (3) of the XRISM Collaboration Nature 2025
prescription. The speed PDF has unit integral. A complete spherical shell
preserves the intrinsic photon integral. A cone retains its direction-dependent
beaming and geometric weight.

The independent three-dimensional Sobol reference integrates delta squared in
wind-axis coordinates. Exact full-sphere finite-bin integrals, uniform and
truncated-Gaussian velocity PDFs, nearly mono-speed limits, signed Gaussian
populations, inflow/bicone symmetries, input-bin overlaps, irregular grids and
logarithmic FFT grids are checked by check.py. Finite-grid clipping transitions
near |beta|=1 are explicit breakpoints in the scalar reference quadrature.

On 2026-11-20, the physics checks passed. The six conical line comparisons with
direct 3D integration had relative L1 errors below 1.9e-4. The direct convolution
comparisons had maximum absolute errors below 8e-8 for unit input photon flux.
The 29 sector and geometry unit tests passed. PyXSPEC uniform, Gaussian and
sector adapter checks passed with the real Resolve response. Adapter comparisons
allow 2e-8 absolute error per bin, within the declared 1e-7 tabulation budget;
round trips through XSPEC can select different adaptive integration panels.

## XSPEC runtime and grid convergence

Native `wuconv` evaluates the entire uniform-speed photon convolution in
C++/FFTW (`native/uniform.h`); other models retain their Python bridge.
The Python reference remains available through `windkin.xspec.register()`.
Use an explicit logarithmic calculation grid for FFT evaluation. A dense
irregular grid retains direct finite-bin integration, which is much slower.

On 2026-11-20, `check_native.py` passed 46 geometries on each of four grids:
200 and 40,000 logarithmic bins, 10 irregular bins including a zero energy
edge, and 4,096 very narrow logarithmic bins. Cases include zero width, a
PDF reaching beta=0, speeds near c, inflow, bicone, equatorial tangencies,
empty cones, and random geometries. Maximum native/reference relative L1
was 1.1e-7 (including clipped narrow-band cases). Full-sphere photon
conservation, an independent finite-source-bin analytic sphere calculation,
signed input spectra, cache invalidation, and invalid parameter rejection
also passed. The compiled XSPEC package passed the real-response adapter
checks, including the unchanged Gaussian and sector models; the 29 geometry
and sector unit tests passed. Cubic tabulation uses the same 1e-7 budget and velocity
quadrature uses the same 2e-9 criterion as the Python reference.

The reproducible `benchmark.py` changes speed, width, geometry, ionization
and turbulence between calls. The Python kernel cache is cleared, and native
kernel parameters change on every trial. On the local Apple Silicon host,
32 trials per configuration with the real Resolve response and hem_total
table gave these median latencies:

| Calculation bins, 0.001–60 keV | wemlos_Cf | wuconv C++ | wuconv Python |
|---|---:|---:|---:|
| 20,000 | 8.94 ms | 8.73 ms | 19.40 ms |
| 40,000 | 17.86 ms | 11.30 ms | 23.53 ms |
| 80,000 | 30.03 ms | 16.72 ms | 34.30 ms |

These include parameter updates, table interpolation, convolution and response
folding for one emission component. At every trial the native and Python folded
spectra were compared; maximum relative L1 was 1.3e-12. The mono-speed model is
a different physical calculation and is included only as a runtime reference.

An additional check loaded the PDS456 6abs1em broadband session with Resolve,
Xtend and NuSTAR, copied evaluated parameter values, and replaced only the
emission convolution by C++ or Python `wuconv`. Across 32 perturbations of its
independent thawed parameters on a 40,000-bin grid, median complete parameter
update/likelihood/folding time was 210.1 ms (C++) versus 236.8 ms (Python), an
11.3% reduction. Maximum folded relative L1 was 5.2e-14 and maximum statistic
difference was 2.8e-10. This validates evaluation performance and equivalence;
it is not a refit or a final BXA configuration. Sampling overhead and worker
hardware will affect the end-to-end gain; server timings remain unmeasured.

At the reference geometry (Cf=0.7, skew=0.5, i=45 degrees, vmean=0.28c,
width=0.17c), wuconv at 20,000 and 40,000 bins differed from 80,000 bins by
3.9e-6 and 1.0e-6 in folded relative L1. The zero-width additive profile agreed
with compiled wlinelos_Cf to relative L1 below 2.9e-4 at i=15,45,75 degrees.
The residual difference includes RED-WINE's center-bin integration versus
windkin's finite-bin integral. These reference checks support a logarithmic
grid for inference; check representative parameter extremes for a new dataset.

Run from the checkout with HEASoft initialized:

    python xspec/windkin/check_native.py
    PYTHONPATH=src:xspec/windkin:$PYTHONPATH python xspec/windkin/check.py
    PYTHONPATH=src python -m pytest tests/test_sector.py tests/test_geometry.py -q
    python xspec/windkin/benchmark.py --spectrum PATH_TO_PHA --table PATH_TO_AMOD --output timing.json

For compiled adapter and real-response checks, use --xspec --native, --redwine and --spectrum
options described in README.md. Start the real-response check from the PHA
directory so XSPEC resolves the response files referenced by its headers.

The native `wuconvlos_Cf` adapter uses RED-WINE's shared LOS geometry with
`(vmean, width, Cf, skew, incl, bicone)`. Three response-folded comparisons
against physical-angle `wuconv` passed. `tests/check_geometry_procedures.tcl`
checks the extra spread-parameter offset alongside the existing mono-speed
layout. xspeclab's export/covariance/profiling remain model-independent.
