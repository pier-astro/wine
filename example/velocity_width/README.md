# Experimental uniform-beta line profile

This prototype evaluates the RED-WINE conical line density at positive energy
centers (keV). The input normalization is the intrinsic line normalization for
a full sphere. The cone is a single hollow cap with half-opening angles
`inner` and `outer` and inclination `inclination`, in radians. Uniform
intrinsic emissivity per unit beta is assumed, independent of position in the
cap. The resulting profile includes the existing relativistic weights.

## Conventions and inversion

The derivation was checked against `01_Relativistic_Emission.md` and
`02_Emission_Profiles_and_Geometry.md` in the xflow-draft notes, the Python
`redwine.kernels.lineshape`, the XSPEC `redwine_lineshape.h`, and the
geometry normalization check. Write μ = cos θ, s = sqrt(1−β²),
g = E/E₀ = s/(1−βμ), and f = wφ/(2π). Then
dμ/dE = s/(β E₀ g²) for positive β. The photon density is
`norm * f * s/(2 beta E0)`. The energy density has an extra factor g.
Upon integrating over E, these are angular integrals weighted by g² and g³,
respectively, with measure dΩ/(4π). Statements about photon number
conservation in the legacy docstrings therefore do not imply an unboosted
observer-frame integrated photon flux for an arbitrary cone.

The mixture averages the original density with dβ/velocity_width. In the
zero-speed limit its integrated flux approaches
`norm * (cos(inner) - cos(outer)) / 2`. It is never renormalized by its
observed integral. A positive-speed zero-width request calls the original
Python kernel directly. A zero-speed, zero-width line is a delta measure
and the point evaluator returns infinity at E₀ for a nonempty line.
A nonzero width reaching β=0 can also give an integrable logarithmic
divergence at E₀ if a ring of nonzero length intersects the wind there.
Use bin integrals before fitting such cases. The experiment covers outward
speeds, with 0 ≤ β_min ≤ β_max < 1. Invalid intervals are rejected without
clipping or changing their weights. XSPEC supports signed fixed speeds;
signed mixtures are outside the experiment.

Set q = 1 − E₀/E and u = tanh(atanh(β)/2). The original Doppler equation
gives (2−q)u² − 2μu + q = 0. For each energy, both physical roots of each
tangency are retained. The tangencies are μ = cos(i±inner), cos(i±outer)
and ±1 for the kinematic limits. Between them the geometry is either
constant or partial. Repeated roots and roots outside the interval
produce no additional panels. Inverse cosine clipping is implemented with
an atan2 expression that retains the small angular span at the equatorial limb.

## Methods and validation

- **Reference** uses scalar adaptive QUADPACK on all split panels with
  absolute and relative tolerances of 2e−12 on the beta integral before
  the final normalization. It shares the analytic breakpoint solver with
  the other methods. Independent checks use the old fixed-beta kernel on
  smooth intervals and forward integration in wind-axis coordinates.
- **Direct** uses vectorized, adaptive 8/16-node Gauss–Legendre pairs in
  raw β on the same breakpoints. Sharing the splits prevents missing narrow
  support and isolates the effect of the semi-analytic integration.
- **Semi** integrates constant arcs with
  `F(beta) = sqrt(1-beta²) + log(beta/(1+sqrt(1-beta²)))`.
  For partial arcs it uses the same quadrature pair after the map
  β = a + (b−a) sin²(πt/2), which softens the square-root behavior at
  tangencies. Successive-rule tolerances are local error indicators.
  A global Legendre expansion was not benchmarked. The piecewise kernel
  is nonsmooth and positive-weight local quadrature already gives
  nonnegative profiles without truncation oscillations.

Run in the experiment directory with RED-WINE installed (or set
`PYTHONPATH` to its `src` directory).

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  python check_and_benchmark.py --output benchmark.csv
```

Omit `--output` to run only the focused checks. The benchmark uses
0.3–25 keV logarithmic grids and a 6.97 keV line. It measures median
geometry preparation over 101 repetitions and evaluation over 3
(reference) or 9 (other methods), after one warmup. Evaluation includes
energy-dependent breakpoint construction, integration, and the output
array. The CSV reports errors against the reference, with relative errors
restricted to values above the listed floor (10⁻⁶ of peak or 10⁻¹⁴,
whichever is greater) and an absolute error for the remaining samples.
The timings describe this Python experiment on an Apple M1 Pro with
Python 3.11.9, NumPy 1.26.4 and SciPy 1.14.1. They do not measure an
XSPEC fit or a compiled implementation.
