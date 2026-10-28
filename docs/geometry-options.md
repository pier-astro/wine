# Geometry constraints: exploration before selecting an implementation

2026-10-27. The exploration below records the route to the agreed API. The
geometry/model implementation and first checks are now recorded in
[geometry-validation.md](geometry-validation.md). Earlier temporary experiments
remain historical evidence, not descriptions of the current public models.

## Current API direction after discussion

This section supersedes the exploratory alternatives below. The coordinate
implementation is in place. Geometry reporting, supported angular profiles,
physical-angle scans, and structured export are now described in
[geometry-procedures.md](geometry-procedures.md); the bullets below record the
agreed design and historical checks. There is no backward-compatibility requirement.
The author plans to rename the project RED WINE; no project rename is performed
as part of the API discussion.

- Use the physical hemispheric `Cf = cos(tin) - cos(tout)` as an explicit
  XSPEC parameter. The existing full-sphere emission scaling is `Cf/2`;
  changing coordinates must not silently change normalization.
- Generic geometry and off-wind geometry use `hole` for the fraction of unused
  hemispheric solid angle assigned to the polar cavity:
  `cos(tin) = 1 - hole*(1-Cf)`, `cos(tout) = (1-hole)*(1-Cf)`.
  Generic `incl` is independent. Off-wind `view` interpolates in angle from
  the lower to upper bound of the selected sector; `region` stays fixed.
- LoS geometry keeps physical `incl` and `Cf` explicit. Its remaining
  coordinate is named `skew`, as selected by the author; it differs from the
  generic/off `hole`. For `s=skew`, `c=cos(incl)`, `x=cos(tout)`, and
  `cos(tin)=x+Cf`, use the physical root
  of `(1-s)*(1-x-Cf)*(c-x) = s*x*(x+Cf-c)`. This quadratic avoids the interior
  derivative changes from the earlier min/max interpolation. Degenerate
  boundaries, numerical derivative probes, and fitting convergence remain to
  be validated.
- For non-degenerate LoS geometry, `skew=0.5` gives equal fractional coverage
  of the available solid angle on the polar and equatorial sides of `incl`.
  With polar coverage `(x+Cf-c)/(1-c)` and equatorial coverage `(c-x)/c`, both
  equal Cf at this reference. Above 0.5 coverage favors the equatorial side;
  below 0.5 it favors the polar side. This does not imply equal amounts of
  wind solid angle on both sides, symmetry in degrees about `incl`, or a
  statement about observed Doppler-weighted flux.
- Python should retain explicit-angle access through the proposed `Geometry`
  constructors rather than duplicate emission calculations for each interface.
- Helper selection uses the model component number, optional model name, and
  optional data group, not an arbitrary parameter number. Selecting a group
  identifies parameters; profiling continues to use the full joint dataset.
- The author's preferred command names are `windgeom`, `winderror`, and
  `windsteppar`. Omit an unavailable covariance error in the printed summary;
  print the value alone. Keep missing/fixed/failed states distinct internally.
- Reuse the existing profiler and its tested `toler` behavior (default 1.0 in
  absolute statistic units), separate from requested delta-stat (default
  2.706) and inner fit convergence. The existing Tcl 8.6 profile checks passed
  again; adapting the physical constraints still requires new focused tests.
- Retain `savesteppar`'s FITS schema for custom scans. Its current command reads
  native XSPEC `tcloutr steppar` state, so it cannot directly consume a custom
  result. Prefer an optional generic structured-result input to that exporter,
  sharing the FITS writer and imposing no WINE runtime dependency on xspeclab.
  The proposed call is `savesteppar result $scan outfile geometry.fits`.
  The author also explicitly accepts a dedicated RED WINE scan-saving procedure
  in xspeclab's existing XSPEC procedures if that is cleaner. Reuse the same
  FITS writing path; no separate plotting format is needed.
- Preserve every requested Cartesian grid cell, using non-finite statistics
  and a separate status for excluded/failed cells. A four-cell FITS fixture
  with physical derived axes, one excluded cell, and `P1`, `P2`, `STAT`,
  `DELTASTAT`, `STATUS` loaded through the existing `StepparData` and
  `make_matrix()` successfully. Full exporter integration and mask-safe
  interpolation/smoothing have not yet been validated.

The direct native exporter test in an isolated process did not complete:
HEASoft `ftcreate` reported an initialization/prompt error even with a writable
temporary PFILES. Do not claim a successful full Tcl export/read round trip
from the separate FITS-reader fixture. Logs are under
`/private/tmp/wind-export-check-*/`; the successful fixture is
`/private/tmp/wind-fits-contract-4d337myx/contract.fits`.

## The actual `err.xcm` failure

The author's `/Users/pierpaolo/Analysis/pds456-wine/xspec/err.xcm` was loaded in
an isolated native XSPEC 12.15.1 process. The diagnostic build changes the
exception message to print the trial parameters. With `parallel leven 1`,
`query no`, and only `fit 1 0.01`, the first iteration rejected:

| Parameter | Rejected trial value |
| --- | ---: |
| `vout` | 0.41550195885492291 |
| `uincl` | 0.84951511477526387 |
| `tout` | 52.961829156675137 deg |
| `tin` | 89.494439449518993 deg |

`uincl` is finite and inside [0, 1]. The failure is an inverted angular
interval. Multiple parameters changed together, including a large increase in
`tin` from its starting value near 6.76 deg: this is a joint optimizer trial,
not merely a small derivative probe outside the `uincl` bounds. Individual
XSPEC parameter limits do not enforce the coupled constraint `tin < tout`.
The existing combined error message obscures this distinction.

The strict diagnostic run took about 11 seconds. Its transcript is
`/private/tmp/wine-geometry-exploration/err_one_iteration.log`.
It is not a full fit or an assessment of the saved scientific solution.

## What the local-model interface can communicate

The installed C/C++ local-function interfaces return `void`, updating flux and
optional flux-error arrays. No dedicated per-evaluation "geometrically excluded,
skip this grid cell" return status is present in that interface. The installed
`Step::ParallelGrid::execute` catches an exception around the scan loop and
rethrows, rather than recording an excluded cell and continuing.

Therefore changing the wording or type of the current exception does not by
itself provide a masked rectangular `steppar` surface. Returning a spectrum is
also not equivalent to directly adding a parameter-constraint penalty to the
fit statistic. No changes to XSPEC itself were attempted.

Sources inspected in the installed HEASoft 6.36 tree:

- `Xspec/src/XSFunctions/Utilities/funcType.h`
- `Xspec/src/XSModel/Model/Component/ConComponent.cxx`
- `Xspec/src/XSFit/Fit/Step.cxx`
- `Xspec/src/XSFit/FitMethod/LevMarq/LevMarq.cxx`

## Small native XSPEC experiments

The model was `powerlaw + windemlos*gaussian`, with expected counts generated
without Poisson fluctuations from a permitted geometry. The setup used a
120-channel diagonal dummy response over 0.5--12 keV, channel offset 0.5 keV,
channel width `(12 - 0.5)/120` keV, exposure 100 ks, and a 0.1--30 keV model
grid. Reloading that same response gives a generating-model C-statistic of
about `1.4e-12`. Counts without fluctuations are an infrastructure diagnostic;
no significance or fit-quality interpretation is assigned to their statistic.

The scan was `steppar best 5 40 80 2 6 20 60 2`, containing six permitted
cells and three cells with `tout <= tin`. The temporary wrapper tested three
behaviors for those three cells:

1. Throw the current exception, augmented with parameter values.
2. Return zero emission from the affected convolution.
3. Return a non-finite emission spectrum (`NaN`).

The public code does not contain the temporary selector for these experiments.
Other invalid inputs and other model families were not generalized by the
prototype.

| Behavior outside the angular domain | Other parameters fixed | Three nuisance parameters refitted |
| --- | --- | --- |
| Exception | Scan aborts; exported statistic array has nine zeros | Same failure |
| Zero emission | Nine finite statistics | Nine finite statistics |
| `NaN` emission | Six finite statistics and three `NaN` entries | Scan aborts with "Unable to fit when starting with fit statistic = NaN"; exported array has nine zeros |

The refitted nuisance parameters were continuum normalization, `uincl`, and
Gaussian emission normalization. Geometry was fixed by `steppar` at each cell.
The zero-emission scan's four overlapping permitted cells agree, to exported
precision, with a separate entirely permitted 2x2 scan, both with fixed and
with refitted nuisance parameters. This supports post-masking on this small
example; it is not general convergence certification for `current` initialization,
parallel scans, multiple components, or real data.

### Zero emission: useful evidence and a concrete limitation

The same one-iteration `err.xcm` test succeeds with the temporary zero-emission
extension. It ends inside the physical domain at approximately `tin=64.09`,
`tout=82.95` deg, `uincl=0.6924`. The starting statistic was 1167.32 and the
one-iteration value was 931.82. These are diagnostic values, not a converged
new scientific fit. This result shows that a finite continuation can let the
native optimizer recover after the rejected geometric trial.

However, a separate toy fit intentionally started at `tout=20`, `tin=30` deg
and kept zero emission. It finished at those inverted angles with finite
C-statistic approximately 128.01 and zero-alpha warnings for the inactive
geometry/emission parameters. The zero extension therefore does **not** enforce
the physical domain. It creates a region equivalent to switching off emission,
where geometry has no derivatives to guide the optimizer back into the wind.
Starting from allowed parameters and checking the final parameters is necessary
if this extension is adopted, but has not been shown to guarantee convergence.

Logs, toy sessions, and the diagnostic runner are under
`/private/tmp/wine-geometry-exploration/`; `results.json` records the small scans.
The real-data zero-extension transcript is `err_one_iteration_zero.log`.

## Alternatives to choose together

| Route | Benefit | Cost or limitation |
| --- | --- | --- |
| Finite zero-emission continuation | Very small model change; keeps `tin` and `tout`; native rectangular scans worked in the toy check, and the real-data first iteration recovered | Permitted-domain statistics must be separated from the artificial extension; it does not prevent an inverted-angle fit or an emission-free plateau |
| Ordered angular coordinates | Guarantees ordered angles throughout fitting without a penalty spectrum; negligible numerical cost | Changes a fitted coordinate and its interpretation; physical-angle contours need coordinate conversion or a physical-domain scan; covariance/profile/export procedures need review |
| Procedure scanning only the permitted triangle | Preserves the physical model and angle names, uses native XSPEC fits at each cell, and stores exclusion/failure statuses explicitly | Additional procedure to maintain; by itself it solves the grid, not unconstrained simultaneous fitting of both opening angles |

An ordered coordinate could use

\[
\theta_{\rm out}=\theta_{\rm in}
 +q(90^\circ-\theta_{\rm in}),\qquad 0<q\leq1.
\]

The zero-volume cases `q=0` or `theta_in=90 deg` still need an explicit limit
policy; the undefined inclination fraction at zero width should not be treated
as a measured angle. An absolute width alone does not remove the coupling to
the upper 90-degree bound. An optional model entry point could preserve legacy
models, but model/parameter names and downstream support should be discussed
before implementing it. Reparameterization preserves the physical model
interior, while a flat prior in a new coordinate is generally a different
Bayesian prior.

A procedure can scan rows: fix `tin`, step `tout` only above it, and refit the
other parameters from the same accepted reference. Both opening angles are then
fixed at each cell, so the original crossing cannot occur during that cell's
fit. Geometrically excluded cells and failed permitted fits need separate
statuses. The procedure should preserve the caller session and save each
completed point immediately. A good grid reference must still be obtained;
this is not an automatic replacement for a robust simultaneous free-angle fit.

Do not silently sort or clip the existing named inner/outer angles. Large
finite spectra used as penalties are not a clean substitute for a constraint:
they alter predictions and can be weakened by a free emission normalization.
Smaller finite-difference steps or a different optimizer alone do not encode
the required angle ordering.

## Earlier decision point (superseded by the agreed API above)

Select whether preserving the present angle coordinates or guaranteeing the
physical domain during every optimizer trial takes priority. Keep the
procedure as a viable fallback for the physical-angle grid. No proposed route
has yet been selected or applied to the working model package.

## Follow-up: fitting reliability and derived-angle reporting

Two additional zero-continuation toy fits started inside the physical domain,
using the same expected-count spectrum and a budget of 20 iterations. A
moderately displaced start (`tout=65`, `tin=25`, `uincl=0.4`, continuum norm
0.011, line norm 0.004) ended at `tout=60.269`, `tin=30.359` deg with C about
`1.41e-4`. A second start increased the line normalization fourfold to 0.02
and ended at `tout=59.541`, `tin=29.523` deg with C about `5.89e-5`.
Neither produced geometry exceptions or zero-alpha warnings. These two
successes do not demonstrate that all allowed starts avoid the artificial
zero-emission plateau. No broad convergence search was performed. Logs and
`valid_start_results.json` are in the temporary exploration directory.

Ordered coordinates remove inverted angles from the accepted parameter box;
they do not remove geometric correlations or guarantee optimizer convergence.
Finite-difference probes at hard bounds, zero-width/zero-volume limits, kernel
sampling, and realistic displaced starts still need focused validation before
any public replacement. In particular, do not claim that mapping a bounded
fraction automatically makes every out-of-bounds numerical derivative probe
physical.

For `tout = tin + q*(90 - tin)`, `tin` is a native fitted angle and `tout` is
derived. The existing viewing-angle tools need an explicit extension; an
`error` interval on `q` cannot be relabelled as the `tout` interval while `tin`
varies. In the local Gaussian approximation (angles in degrees),

\[
\operatorname{Var}(\theta_{\rm out})\simeq
(1-q)^2\operatorname{Var}(\theta_{\rm in})
+(90^\circ-\theta_{\rm in})^2\operatorname{Var}(q)
+2(1-q)(90^\circ-\theta_{\rm in})
\operatorname{Cov}(\theta_{\rm in},q).
\]

For non-Gaussian/asymmetric likelihood limits, profile the physical outer angle
itself. The existing trial-link approach is applicable: holding `tout=T`
implies `q=(T-tin)/(90-tin)` with a compatible range for `tin`. For posterior
intervals, transform each joint chain sample and take the desired quantiles.
These are different uncertainty summaries and should remain labelled as such.

A joint profiled `steppar` surface in `(tin,q)` transforms pointwise into
`(tin,tout)` without any further fits: keep the statistic and transform both
coordinate arrays. The resulting grid is non-rectangular and has degenerate
boundary points. It is not correct to merely relabel the `q` axis or transform
one-dimensional confidence endpoints at fixed best-fit `tin`. Finite grid
sampling/interpolation remains an approximation. Likelihood contours are
coordinate invariant; Bayesian prior measures require the appropriate
Jacobian if the original physical prior is to be retained.

The recommended next candidate is an ordered-coordinate prototype, assessed
before changing the public models. Its implementation costs are deterministic
geometry conversion and adapted error/export/plotting support. Zero
continuation should not be adopted as the default unrestricted fitting policy
on the current evidence. This is a design recommendation, not an implemented
or convergence-certified replacement.
