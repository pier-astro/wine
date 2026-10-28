# Geometry implementation: first validation pass

2026-10-27. Francesco and the author agreed on the parameterization. This
pass implements and checks the geometry/model interface; the consolidated
physical-angle procedures and scan export are the next stage. The project
rename has not been performed, and changes remain uncommitted.

## Implemented interface

| XSPEC family | Geometry coordinates after `vout` |
| --- | --- |
| `windem`, `windline` | `Cf, hole, incl` |
| `windemlos`, `windlinelos` | `Cf, skew, incl` |
| `windemoff`, `windlineoff` | `Cf, hole, view, region` |

Line models then take `LineE` and their native additive normalization. The
shared physical profile/convolution has not been replaced or optimized. The
normalization remains `Cf/2` relative to a full sphere, with hemispheric
`Cf=cos(tin)-cos(tout)`. The line-energy registration also now has consistent
hard/soft bounds (the former definition produced a bounds warning).

Python exports `Geometry.angles`, `.free`, `.los`, and `.off`. Their public
angles are degrees, and constructors broadcast for chain/scan conversion.
`emission_slab` takes a scalar `geometry` object and converts to radians at
the kernel boundary. The low-level kernel convention remains radians.

This is a breaking development API. In particular, saved parameter counts
can remain the same while meanings and bounds change. The two author-provided
XCMs were converted only in temporary diagnostic copies; the originals were
not edited. Generic xspeclab parameter/covariance exports do not hard-code the
retired `uincl`/`uview` coordinate order. The older WINE `wineincl*` procedures
are explicitly documented as inapplicable to the new coordinates; they are
retained as the tested search implementation for the upcoming consolidation.

## Numerical handling

`wine_geometry.h` is used by all six XSPEC emission entry points. It computes
the same physical angles as Python. The LoS quadratic uses a positive-term
discriminant and a stable root expression that also covers `skew=0.5`.
Tests compare it to a separate bracketed root, including endpoints and near
zero/full coverage. `Cf=0` clears the emitted flux rather than leaving stale
values from the previous evaluation.

XSPEC's finite derivative probes can leave individual hard ranges. The
adapter evaluates the nearest physical boundary for finite out-of-range
probes and issues a notice once per parameter/spectrum. This is a documented
one-sided numerical continuation, not permission for out-of-range accepted
parameters or links. Non-finite inputs remain errors. Python constructors
instead reject out-of-range coordinates; inverted explicit angles warn and
represent zero emission. Degenerate boundary coordinates do not imply an
angle measurement.

## Checks and results

- `python -m unittest discover -s tests -q`: 11 tests passed, including an
  independent bracketed reference for the LoS geometry and existing grid/table
  model checks.
- `python xspec/test_geometry_models.py`: passed generic versus LoS/off
  spectra at equivalent physical angles, boundary cases, registration,
  zero-coverage reuse, and the zero-velocity hemispheric normalization.
- Thirty spectra recorded from the pre-change explicit-angle implementation
  (convolution and line families, three velocities and five geometries) match
  the new generic coordinates. Maximum peak-scaled difference was about
  `1.9e-16`; the spectral physics is preserved at these reference points.
- `python xspec/test_energy_coverage.py`: passed; the extended-grid differences
  remain below `7e-13` for the five shift/broadening checks.
- `python xspec/test_geometry_fit.py`: known-spectrum recovery from a
  constructor-like, displaced start reached C approximately `1.1e-10` for
  expected counts generated at `Cf=0.4, skew=0.65, incl=35 deg`. Native
  `error 1.0` on `Cf` returned a bracketing interval with clear diagnostic
  flags, and a native profiled 3x3 `(Cf,incl)` scan completed with nine finite
  statistics, including coverage 0/1 and inclination 0/90 deg.
- The converted `err.xcm` reproducer started with the same statistic as the
  previous coordinates (`1167.31664`). One native Levenberg iteration completed
  without the formerly observed geometry exception, reaching `974.1843284`.
  Serial and five-process derivative runs gave the same result. This is an
  integration check, not a converged new PDS 456 fit.

Expected-count simulations have no Poisson fluctuations. Their fit/profile
checks demonstrate operation of the interface and constraints, not a detection
significance or scientific uncertainty calibration. Off-model fits, full
real-data convergence, and the future physical-angle profiling/scanning
procedures have not been certified by this pass.

## Shallow likelihood and distant starts

One more distant synthetic start (`Cf=0.6, skew=0.25, incl=60 deg`) stopped at
a secondary solution with C approximately 0.66 rather than the generating
minimum. Increasing the energy sampling or reducing the derivative step did
not remove that secondary solution. The former angular coordinates also had
difficulty from the equivalent physical start. This is consistent with the
known shallow/degenerate likelihood; it does not contradict the independent
spectral and geometric equivalence checks.

Following the author's direction, no further effort is spent forcing recovery
from distant starts. This historical diagnostic is not a blocker for the
coordinate implementation, and is not hidden as a successful global recovery.
Likelihood exploration with `error` and `steppar` remains necessary to assess
the scientific solution and competing minima. The routine fit check therefore
tests a reproducible constructor-like recovery and functional native profiles
and scans, rather than asserting that every remote initialization converges
to the same minimum.

Temporary references, converted sessions, build output, and real-data logs
are under `/private/tmp/wine-api-prototype/`. The small checks in the repository
are runnable without the author's photoionization tables or observations.

## Fixed-opening-angle native-link prototype

On 2026-10-27, `xspec/test_fixed_angles.py` passed in initialized PyXSPEC:
`powerlaw+windemlos*gaussian`, the same lightweight 120-channel response and
2000-bin calculation grid, and expected counts without Poisson fluctuations.
For fixed `tin=20`, `tout=60` degrees, `Cf` is frozen to their cosine difference.
The inverse `skew=A/(A+B)` is linked to `incl` through native `cosd(p6)`;
inclination has hard/soft bounds 20--60 degrees. Only inclination and continuum
normalization are free in this prototype.

- Six native link evaluations at inclinations 20, 25, 35, 45, 55, and 60 degrees
  recover both opening angles to within `1e-7` degrees.
- Starting at 45 degrees and continuum normalization 0.008, the fit recovers
  the generating inclination of 35 degrees as `35.00006366`, with
  `C=3.25e-10`. The final opening angles remain 20 and 60 degrees.
- Constants in native links must use real syntax: `1.0`, not integer `1`,
  which XSPEC treats as a parameter reference.
- At the upper inclination endpoint, degree-function rounding gives
  `skew=-8.08e-13`. The existing model adapter bounds such numerical probes;
  the strict Python comparison explicitly checks a `1e-12` coordinate tolerance
  before clipping for that reference calculation. No model changes were made.

This establishes a native-link route for this nondegenerate constrained fit,
not the general error/scan procedure. Full/zero coverage, other singular
limits, joint-model identity links, excluded/failed scan cells, session
restoration, and export remain separate implementation checks. No PDS 456
scientific refit was performed.

## Angular limits, identity links, and tiny scan

The follow-up `xspec/test_angle_scan.py` passed on 2026-10-27, using two
`windemlos*gaussian` emitters with identity-linked geometry and velocity in
one model. Both emitters' opening angles were checked at three inclinations
for seven angular configurations, including `(0,0)`, `(90,90)`, `(0,90)`,
`(0,60)`, `(20,90)`, `(20,20)`, and `(20,60)` degrees.

The exact inverse limits need no epsilon denominator:

- Full hemisphere: freeze `Cf=1`, `skew=0.5`; inclination remains free.
- Polar cap (`tin=0`, `tout<90`): freeze `skew=0`.
- Equatorial belt (`tin>0`, `tout=90`): freeze `skew=1`.
- Zero width: freeze `Cf=0`, `skew=0.5`, and `incl=tin=tout` for the formal
  LoS geometry. The vanishing emitter provides no inclination information.

A 3x3 grid with `tin={0,20,60}`, `tout={20,60,90}` returned eight finite
profiled statistics and one excluded `tin>tout` cell (`NaN`, not zero).
Both diagonal cells give the continuum-only statistic `C=357.85377`.
The generating angular pair `(20,60)` recovers `C=1.84e-10`.
Forward and reverse traversal agree within absolute statistic tolerance
`1e-4`, using the same nuisance starts at every point. Sixteen permitted
fits were performed across the two traversals.

These results support the constrained LoS scan mechanism without new model
variants or production-kernel changes. They do not yet check multiple data
groups/named models, off-wind profiling, failed-fit status recovery, caller
session preservation, public Tcl parsing, or FITS export and plotting.

## Public Tcl scan and FITS integration

`xspec/test_geometry_procedures.py` passes with native XSPEC 12.15.1. The
implemented commands `windgeom` and two-axis `windsteppar` resolve component
numbers, optional named models, and data groups, following identity ties to
their masters. A named two-group test selects a linked emitter in group 2
while fitting the full dataset. Reports use XSPEC-style indentation/columns.

The live check preserves the caller's serialized state byte for byte and its
statistic/energy grid. A disposable worker reloads the same reference at each
cell. File-backed responses are essential for reloading the test: native
`save all` does not restore a dummy response. It also omits autonomous energies
and rounds saved parameter settings to six significant digits. The new scan
therefore copies the calculation grid explicitly and refreshes unlinked
parameter settings using `tclout`'s ten-digit representation. Across multiple
spectra, a common model calculation grid is currently required.

Injected Tcl fit errors are retained as `failed` with their original message;
geometric exclusions remain `excluded`, and neither receives a measured
statistic. An XSPEC process abort retains logs and the incrementally recorded
completed points, rather than presenting an incomplete scan as success.
New minima improving the reference by more than 1.0 statistic unit are reported
and saved for inspection without altering the caller's fit.

The xspeclab generic `savesteppar result $scan` path reuses the FITS writer.
Live native and structured exports pass. Metadata are written to the table
extension expected by the reader. Excluded/failed values use FTOOLS `INDEF`
to become FITS undefined floats (not zero); physical axes have no fabricated
parameter IDs. Reader status/masking and plotting have focused portable tests;
contours do not bridge masked corners, and smoothing/interpolation of masked
surfaces are rejected. Eight xspeclab tests pass.

This is a completed scan/report/export stage, not certification of derived
one-dimensional errors, covariance reporting, global convergence, or scientific
PDS 456 confidence contours. At that stage `winderror` was still pending;
the following section records its subsequent implementation.

## Single-angle profiles and local covariance

The subsequent procedure stage implements `winderror COMPONENT QUANTITY`
and local covariance reporting in `windgeom`. It retains native-style
indentation, two-decimal confidence-limit display, full Tcl result values,
and separate WINE/native flag formats. No model kernel or normalization
was changed.

`test_geometry_procedures.py` now passes live LoS profiles for tin/tout at
delta-stat 1, direct inclination errors, fixed-Cf angular profiles,
fixed-Cf/skew monotonic angle transformations, named two-group inclination
errors, and off-wind inclination with fixed Cf/hole. Native hole/Cf links are
also checked for off-wind opening-angle trials. The fit session remains
unchanged across the public error calls. The toy tests use expected counts
without fluctuations and do not calibrate significance or global convergence.

For the two-emitter LoS toy, local propagated errors are approximately
`sigma(Cf)=0.0340`, `sigma(tin)=5.87 deg`, `sigma(tout)=2.71 deg`, and
`sigma(incl)=3.07 deg`. The asymmetric delta-stat 1 profiles are approximately
`tin=5.52--25.41 deg`, `tout=54.96--62.62 deg`, around reference 20/60 degrees.
They are illustrative functional results, not PDS 456 uncertainties or an
assertion that all optimizer basins were found. Covariance propagation agrees
with independent numerical `coverror` evaluations for Cf and a derived angular
expression, including off-diagonal terms. Pure Tcl checks compare the analytic
Jacobians to central differences in all four geometry/sector cases.

Trial initialization retains the reference's other physical opening angle
where possible. This avoids turning a nearby physical trial into Cf near 1
merely by retaining a dimensionless coordinate near tin=0. Additional live
trials at tin=1e-5 and tout=89.99999 degrees produced finite constrained
statistics. These are one-sided endpoint approaches; exact 0/90 spectra are
covered by the angular scan instead. The search records its actual bounds.

Boundary flags are not measured limits; a displayed zero is a sentinel.
Fixed, singular, and unavailable covariance uncertainties remain distinct
internally and do not become printed zero errors. Free off-wind geometry still
requires a dedicated physical-inclination constraint: this version rejects
`winderror incl` in that case rather than relabelling a view interval. Fixed
endpoint-skew charts and restricted hard ranges on an eliminated free Cf are
also rejected explicitly. Covariance is omitted if its current fit-parameter
mapping cannot be verified, for example in some mixed inactive-model sessions.
The [procedure guide](../procedures/README.md) defines the supported cases.

## Clean public example: two 20x20 profiles

The author's cleanup request replaces the example's Python test/benchmark
scripts with a short native XSPEC Bash simulation recipe and one general
plotting notebook. Only the Poisson spectrum is retained. `ftgrouppha optmin`
uses ten minimum counts; this is adaptive grouping rather than the earlier
fixed 20 eV grouping. The model retains four free parameters (continuum
normalization and shared Cf/skew/incl).

For this realization the fit gives C=193.9989892 with 188 degrees of freedom
from 192 noticed grouped bins. The native `(incl,Cf)` profile has 400 finite
points. The physical `(tout,tin)` profile has 400 stored cells: 210 finite
physical fits and 190 masked exclusions, with no failed permitted fits.
Nineteen intervals produce twenty samples on each axis over the full range.
The native run preserves the free skew state despite the expected weak-
sensitivity/zero-alpha warnings at degenerate boundaries. All native Cf=0
points and all zero-width angular points give the same continuum-only
statistic, C approximately 286.1797, rather than artificial zero statistics.

The user-facing workflow is in `example/README.md`. The general adaptation
of the author's plotting notebook is `example/results/geometry.ipynb`, using
saved FITS, joint two-parameter likelihood levels and constant-Cf curves.
No interpolation/smoothing crosses the excluded angular triangle. The plot
marks grid minima explicitly and does not equate them to continuous best fits.
These are first-look likelihood maps, not a coverage or global-minimum study.


## Skew derivatives and minimizer steps (2026-10-28)

A focused diagnostic used the compiled `windemlos*gaussian` model, beta=0.2,
a 6.7 keV line with sigma=0.025 keV, and the example response. Three interior
geometries `(Cf,skew,incl)` were tested: `(0.25,0.64,50)`,
`(0.35,0.25,45)`, and `(0.7,0.5,30)` (degrees). Centered derivatives used
skew steps 0.01, 0.003 and 0.001 on 3000, 6000 and 12000 logarithmic bins
from 0.01 to 40 keV; a 24000-bin derivative with step 0.003 was a numerical
comparison, not an exact reference. All comparisons used folded detector
spectra on the same channels.

At 6000 bins the relative L2 derivative differences against that comparison
were 2.2-6.2% for step 0.01, 8.9-16.6% for 0.003, and 9.6-17.7% for
0.001. Folded spectrum differences were 0.20-0.59%. Doubling to 12000 bins
reduced some derivative discrepancies but did not remove their dependence
on step. Small steps can resolve the sampled kernel's edge structure rather
than a smoother response. These checks establish nonzero skew sensitivity;
they do not identify the cause of any particular Levenberg pivot warning or
validate all geometries, line widths, or spectral grids.

No kernel, parameter bounds, or default deltas were changed. Do not reduce
steps merely to improve formal numerical precision. The existing proportional
`xset delta 0.01` is a reasonable starting point; it also scales the very small
emission normalization. Globally disabling proportional deltas would require
checking every parameter's stored absolute step. A general bin-integrated
kernel could merit a separate change if derivative noise limits an actual fit;
it is not justified as a quick correction by this diagnostic alone.
