# Repository structure and development

The repository has two independent user interfaces:

- `src/wine/`: importable Python package with `geometry.py` (physical geometry and coordinate conversion), `kernels.py` (analytical profiles and convolution), `models.py` (observer-frame slab functions), `grid.py` (engine-independent spectral interpolation), and `ogip.py` (OGIP FITS reader).
- `xspec/`: C++ local XSPEC functions registered in `lmodel.dat`. Python and XSPEC can read the same precomputed spectra, but their default interpolation of spectral values differs; see the [Python guide](python.md#interpolation).

`tests/` checks Python grids, slab evaluation, and OGIP table reading. `xspec/test_geometry_models.py` checks that the constrained geometry wrappers agree with the base XSPEC models. Run the Python tests from the repository root:

```bash
python -m unittest discover -s tests -q
```

`python tests/check_geometry_normalization.py` isolates the geometric
azimuthal weight, without Doppler factors, convolution, or XSPEC. It checks
`integral w_phi(mu) dmu / (4*pi) = (cos(tin)-cos(tout))/2`, independently of
inclination, against the analytic solid angle. A full sphere gives 1; a single
hemisphere gives 1/2. For the public hemispheric covering fraction this is
`Cf/2`. These are nonrelativistic flux fractions relative to a full-sphere
input spectrum; the check makes no claim about finite-velocity boosting.

After an XSPEC build, run `python xspec/test_geometry_models.py` with a Python interpreter that provides PyXSPEC. The local build and cleanup commands are documented in the [XSPEC guide](xspec.md).

`python xspec/test_energy_coverage.py` checks shifts and emission broadening on
matched log-energy grids. It verifies that extending input coverage stabilizes
the common output band and that a constant transmission matches its analytic
reference. It uses built-in spectra, so no external photoionization tables or
observational data are required.

`python xspec/test_geometry_fit.py` tests known-spectrum recovery from a
constructor-like start, a native `Cf` profile, and a small rectangular
`(Cf,incl)` scan including physical boundary points. It uses expected counts
without fluctuations. It is not a global-convergence or significance test;
the [validation report](geometry-validation.md) records the distant-start
diagnostic separately.

`python xspec/test_fixed_angles.py` is the focused native-link prototype for
a LoS fit at fixed `tin=20`, `tout=60` degrees, with inclination free inside
that interval. It checks six link evaluations including endpoints, then fits
expected counts from a displaced inclination and continuum normalization.
It does not yet test singular geometries, linked components, or scan export.

`python xspec/test_angle_scan.py` extends that prototype to exact angular
limits and two identity-linked emitters. It checks a 3x3 fixed-opening-angle
profile grid in both traversal orders, including zero/full coverage and an
excluded cell. This remains a standalone constraint check, not the public
Tcl procedure or a session/export test.

`python xspec/test_geometry_procedures.py` exercises the public Tcl geometry
report and angular scan, with a file-backed unit-area diagonal RMF. It verifies
identity links, named models/two groups, caller-state preservation, simulated
fit failures versus geometric exclusions, and native/structured FITS export.
This integration check requires xspeclab's checkout beside WINE for the export
test; the WINE procedures themselves do not depend on xspeclab.
The same check now tests derived LoS opening-angle profiles, native inclination
errors, fixed-Cf constraints, monotonic fixed-geometry transformations, and
off-wind inclination with fixed geometry. Local covariance uncertainties agree
with xspeclab's independent numerical `coverror` propagation, including an
angular expression with cross terms. `tclsh procedures/test_wind_geometry.tcl`
checks analytic geometry derivatives against finite differences for generic,
LoS, and both off-wind sectors, plus covariance cross terms and profile-boundary
callbacks. The old `test_wine_profile.tcl` remains a regression check for search
tolerance, diagnostics, and layout.

`xspec/wine_geometry.h` shares the geometry conversion between all six emission
entry points. The Python geometry tests use an independent bracketed root as
a reference for the analytic quadratic. Accepted XSPEC coordinates use bounded
metadata; `wine_xspec.h` provides a warned one-sided continuation for finite
out-of-bounds derivative probes. This continuation is separate from Python's
strict constructor validation and does not make out-of-range links physical.

## Spectral boundary

The [response notebook](../example/README.md) builds an illustrative
50 eV Gaussian detector from a text effective-area curve, without XSPEC. Its
code cells have been checked sequentially with Python, including sparse row
probabilities, finite-band losses, and FITS integrity. Live XSPEC line folding
checks the centroid, binned width and area-times-flux normalization.
The [two-zone generator](../example/README.md) is now a short Bash recipe using
native XSPEC `fakeit none` and `ftgrouppha optmin` with ten minimum counts.
It requires the external photoionization tables, not PyXSPEC. Only the Poisson
example is retained. Analysis commands are documented in the example README;
`example/results/geometry.ipynb` reads the saved native and angular profiles.
Earlier expected-count recovery and timing checks remain historical validation,
not extra user-facing scripts. These examples do not certify global convergence
or scientific uncertainty calibration.

WINE receives rest-frame spectra as arrays or OGIP tables and does not import or execute CLOUDY or XSTAR. Spectral production workflows belong in a separate tool. This separation lets a user supply custom spectra while keeping wind kinematics and geometry in WINE. Check units, energy bins, parameter definitions, and any radiative-transfer assumptions before combining files from different producers.
