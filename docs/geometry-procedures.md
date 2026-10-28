# Geometry procedures

2026-10-27. `windgeom`, the two-opening-angle `windsteppar`, and structured FITS
export are implemented and tested. `winderror` adds derived opening-angle
profiles and supported inclination cases; `windgeom` includes local covariance
uncertainties. The `wineincl*` commands still target the retired coordinates.

## Three commands, one component selector

```tcl
windgeom 2
winderror 2 tin delta 1.0 toler 1.0
windsteppar 2 tin 0 90 4 tout 0 90 4
savesteppar result $xspec_tclout outfile geometry.fits
```

`2` is the XSPEC model component number, not a parameter number. Optional
`model NAME` and `group N` select a named model and data group; defaults are
the unnamed model and group 1. Selecting a group identifies the component;
fits still use the complete joint dataset. Check component type before using
its parameter layout. Follow identity links; reject unsupported expressions
with a specific explanation rather than silently changing them.

- `windgeom COMPONENT`: report `Cf`, `tin`, `tout`, and physical `incl`.
  Reports full-covariance local Gaussian 1-sigma uncertainties when their mapping
  and analytic derivatives are valid. Fixed, singular, and unavailable cases
  remain distinct internally and do not print a fictitious zero error.
  No separate covariance command is needed.
- `winderror COMPONENT QUANTITY`: profile `tin`, `tout`, or supported `incl`.
  Retain the tested search, diagnostics, session preservation, and `toler`
  semantics. Defaults: `delta 2.706`, `toler 1.0` in absolute statistic
  units. Direct fitted `Cf` and `incl` remain accessible with native `error`;
  no dedicated WINE commands for those parameters are needed. An interval on
  `hole`, `skew`, or `view` cannot simply be relabelled as an angular interval
  while correlated parameters vary.
- `windsteppar COMPONENT AXIS MIN MAX INTERVALS ...`: profile physical axes,
  starting with the requested `(tin,tout)` use case. Like native `steppar`,
  four intervals mean five samples including both endpoints. Store a Tcl
  dictionary in `::xspec_tclout` containing axes, reference statistic, per-cell statistics and
  statuses, and metadata for export. No new plotting command.

The scan result in `::xspec_tclout` contains schema version 1, `axes`, `rows`, reference `zerostat`,
`statmethod`, `dof`, and `bins`. Axes contain `name`, `unit`, and `values`;
rows contain `values`, `stat`, `deltastat`, `status`, and `message`.
It uses shared geometry/context resolution and an isolated constrained-fit
worker, reloading the same accepted reference at every point. Completed points
are recorded incrementally; the caller's fit session is preserved. The error
command shares context resolution and full-state snapshots, reusing the tested
legacy bracket/bisection through an optional trial callback.

## Why an angular scan needs constrained fits

For generic geometry, `(Cf,hole)` determines `(tin,tout)` independently of
inclination. A joint native profile grid in these coordinates can be mapped
pointwise without refitting, although its physical grid is not rectangular.
For LoS geometry the mapping also depends on inclination. Converting a
`(Cf,skew)` surface after profiling inclination does not in general give a
profile at fixed physical opening angles.

For a fixed nondegenerate angular pair, define
`r=cos(tin)`, `x=cos(tout)`, and `c=cos(incl)`. The inverse relations are

```text
Cf   = r - x
hole = (1-r)/(1-Cf)
A    = (1-r)*(c-x)
B    = x*(r-c)
skew = A/(A+B)
```

The LoS relation requires `tin <= incl <= tout`. Off-wind geometry requires
inclination in its selected polar or equatorial sector. The worker must enforce
these physical constraints during refitting, not just at the starting point.
The inverse formulas suggest reusing parameter links and bounds, but their
singular limits still need targeted tests. The live
[`test_fixed_angles.py`](../xspec/test_fixed_angles.py) prototype confirms
native expression support and one displaced-start constrained LoS fit with
`tin=20`, `tout=60` degrees; see the [results](geometry-validation.md#fixed-opening-angle-native-link-prototype).
Do not add model variants to bypass the remaining boundary checks.

The subsequent [`test_angle_scan.py`](../xspec/test_angle_scan.py) verifies
exact LoS angular limits, two identity-linked emitters in one model, and a
tiny physical-angle scan in both orders. The boundary treatment is recorded
in the [validation report](geometry-validation.md#angular-limits-identity-links-and-tiny-scan).
The public Tcl implementation is in `procedures/wind_geometry.tcl`.

## Export and first implementation gate

Keep every requested Cartesian cell. An excluded `tin > tout` point and a
failed permitted fit must have different statuses and non-finite statistics.
Treat `Cf=0`, `Cf=1`, collapsed off-wind sectors, and undefined inverse
coordinates explicitly; a vanishing denominator does not by itself establish
that the physical spectrum is forbidden.

Reuse xspeclab's FITS writer through the proposed generic `result` input:
`P1`, `P2`, `STAT`, `DELTASTAT`, plus `STATUS`; physical axes must not receive
fabricated XSPEC parameter indices. WINE must not require xspeclab at runtime.
The reader now exposes `status` and masks excluded/failed cells. Contours use
whole-cell masking; interpolation/smoothing of masked grids is rejected.
The export/read/plot path and failed/excluded accounting have been checked.

`xspec/test_geometry_procedures.py` tests real native Tcl commands, same-model
identity ties, named models/two groups, complete caller snapshot preservation,
failed versus excluded cells, and both native and structured FITS export.
It uses a file-backed diagonal response because native snapshots cannot
restore a dummy response. The child grid is explicitly copied; `save all`
alone does not preserve an autonomous calculation grid. The current scan
requires a common grid across spectra and refreshes snapshot parameters with
`tclout`'s ten-digit representation. It does not certify optimizer convergence.

## Single-angle profiles and remaining limits

For free-Cf LoS profiles, define `c=cos(incl)` and `s=skew`. At fixed
`r=cos(tin)`, let `q=r-c`, `A=(1-s)*(1-r)` and link
`Cf=q*(A+s*r)/(A+s*q)`. At fixed `x=cos(tout)`, let `d=c-x` and link
`Cf=d*((1-s)*(1-x)+s*x)/((1-s)*d+s*x)`. Inclination is bounded to the
compatible side of the selected angle. Generic/off geometry links
`Cf=1-(1-r)/hole` or `Cf=1-x/(1-hole)`, with compatible hole bounds.
The trial starts from the reference's other physical opening angle where
possible, rather than preserving a dimensionless coordinate that becomes
ill-conditioned near the angular endpoint.

Fixed Cf uses the two-angle constraint; fixed Cf/skew LoS geometry permits a
monotonic native-inclination transformation. Direct generic/LoS inclination
uses native `error`. Off-wind inclination is supported only with fixed Cf/hole,
where view-to-angle conversion is exact. Free off-wind geometry still needs a
dedicated physical-inclination constraint; it is explicitly refused, not
approximated by relabelling a view interval.

The full Cf hard range is currently required when Cf is eliminated by a link.
Fixed endpoint-skew charts and entirely fixed requested quantities are rejected.
Derived searches approach 0/90 degrees to `1e-5` degrees; exact boundary spectra
remain available in the two-angle scan. Returned `sides`, `angular_bounds`,
`flag_format`, and (for transformations) `native_limits` distinguish boundary
sentinels from measured limits and four WINE flags from nine native flags.
Covariance reporting is local, not a substitute for the asymmetric profile.
