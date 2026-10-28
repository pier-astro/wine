# WINE XSPEC procedures

The current geometry API uses `Cf`, `hole`/`skew`, and `incl`/`view`.
Load the new commands with:

```tcl
source /absolute/path/to/wine/procedures/wind_geometry.tcl
windgeom 2
winderror 2 tin delta 1.0 toler 1.0
windsteppar 2 tin 0 90 4 tout 0 90 4
```

`2` selects a model component, not a parameter. Optional `model NAME` and
`group N` identify its context; all datasets remain included in the fits.
`windgeom` reports physical values in an indented XSPEC-style table and stores
a Tcl dictionary in `::xspec_tclout`. When a valid covariance matrix can be mapped and its order
verified against native parameter sigmas, it includes local Gaussian 1-sigma
uncertainties propagated with the full covariance, including cross terms.
Analytic derivatives use degrees. Fixed, missing, or singular uncertainties
are omitted from the printed values and distinguished in `uncertainty_state`.
Uncertainties are not printed if the free-parameter mapping is ambiguous,
including some sessions retaining inactive models.
`windsteppar` currently supports exactly the two physical opening-angle
axes, in either order. Four intervals give five samples per axis.

The scan fits in an isolated native XSPEC process and leaves the calling fit
unchanged. Every requested cell is retained as `ok`, `excluded`, or `failed`,
with statistics, delta-statistics, and diagnostic messages. Failed/excluded
statistics are non-finite. Zero-width emission is a permitted limit, subject
to any frozen LoS inclination constraints. Each fit starts from the same saved
reference. A completed finite fit does not certify global convergence.

The current version requires a common model energy grid across spectra:
set `energies` explicitly when required. The grid is copied to the worker;
native `save all` alone does not preserve it. The snapshot refreshes unlinked
parameter settings using `tclout`'s ten-digit values. Non-identity geometric
links are rejected. Aborted worker processes retain a snapshot, log, and the
completed points; failed fits retain their diagnostics. A newly found minimum
improving the reference by more than 1.0 statistic unit is reported and saved
for inspection, without replacing the caller's fit.

For FITS export, optionally load xspeclab's exporter:

```tcl
source /absolute/path/to/xspeclab/procedures/export_steppar.tcl
savesteppar result $xspec_tclout outfile geometry.fits
```

No xspeclab dependency is needed for the scan itself. The generic exporter
preserves statuses and physical axis names without fabricated parameter IDs.
Use native XSPEC `error` for directly fitted parameters, or `winderror` for the
physical quantities `tin`, `tout`, and supported `incl` cases below. The old
commands below cannot be used with the new coordinates.

## Angular confidence ranges

`winderror COMPONENT QUANTITY ?delta VALUE? ?toler VALUE? ?model NAME? ?group N?`
uses an isolated process, preserves the caller's state, and stores its dictionary
in `::xspec_tclout`.
Defaults are `delta 2.706` and `toler 1.0` in statistic units. `delta` is the
requested profile increment; `toler` controls quiet approximate crossings and
better-minimum diagnostics, not the inner fit convergence. Original XSPEC fit
settings are retained. For `tin`/`tout`, the other originally free geometry and
spectral parameters remain nuisance parameters. Each trial starts from the same
refitted reference, with a physically consistent geometry initialization.

- With free `Cf`, native links impose the selected opening angle. LoS geometry
  restricts inclination to the compatible side of that angle; generic/off
  geometry restricts `hole` instead. The eliminated `Cf` needs its full hard
  range 0--1. Fixed endpoint `skew` charts are explicitly rejected rather than
  regularized by arbitrary denominator floors.
- With fixed `Cf`, the requested angle also fixes the other opening angle;
  the existing two-angle constraint is reused. With fixed `Cf` and `skew`,
  LoS angular limits can instead be transformed monotonically from native
  inclination limits. Fixed constraints are never silently released.
- Generic/LoS `incl` uses native XSPEC `error`. Off-wind `incl` is supported
  when **both `Cf` and `hole` are fixed**: native `view` limits transform exactly
  to degrees. With free off-wind geometry this conversion is not a physical
  inclination profile, so this version explicitly refuses it.

Derived searches reuse the existing bracketing/bisection and tested tolerance
semantics. Mathematical angular endpoints are approached to `1e-5` degrees;
fixed-coordinate constraints may narrow this range. Four `flags` positions
(`flag_format=wind`) mark lower/upper boundaries and lower/upper numerical
instability. A displayed limit of 0 on a boundary is a pegged-limit sentinel,
**not a measured bound**; `sides` and `angular_bounds` retain actual search
limits. Native errors retain XSPEC's distinct nine-position flags
(`flag_format=native`); monotonic conversions retain `native_limits` too.
Neither kind of interval is a posterior interval.

The explicit derived search stops if a new minimum exceeds `toler` and retains
a candidate snapshot, calculation grid, and diagnostics for inspection. Failed
fits return an error, not confidence limits; completed trial statistics are
recorded incrementally. Native error status is preserved rather than translated
into the four WINE flags. These checks do not establish global convergence.

See [the interface and validation status](../docs/geometry-procedures.md).

## Legacy reference: retired angular coordinates

`wine_inclination.tcl` provides `wineinclcov` (local full-covariance error)
and `wineinclerror` (a wrapper around native XSPEC `error` for direct `incl`
or fixed opening angles). `wine_profile.tcl` provides `wineinclprofile`, the
XSPEC `error` counterpart for the physical viewing angle when angles are free. It
refits those angles and other nuisance parameters on the existing `los` and
`off` models.

These commands require XSPEC and the local WINE model build, but not xspeclab.
The generic xspeclab `coverror` command remains useful for derived quantities
outside WINE. See the [XSPEC guide](../docs/xspec.md#geometric-uncertainty-in-the-development-api)
for usage and limitations.

## Permanent installation

Set `USER_SCRIPT_DIRECTORY` in `~/.xspec/Xspec.init` to WINE's `procedures/`
directory:

```text
USER_SCRIPT_DIRECTORY: /absolute/path/to/wine/procedures
```

Generate the command index once in XSPEC, then start a new session:

```tcl
cd /absolute/path/to/wine/procedures
auto_mkindex .
exit
```

If `USER_SCRIPT_DIRECTORY` already points to a directory of your own, keep it
and link `wind_geometry.tcl` into that directory instead. Run `auto_mkindex .`
there. The generated `tclIndex` is ignored by Git. Regenerate it when adding or
renaming procedures; editing an indexed procedure's body only requires a new
XSPEC session. The WINE model itself still needs `lmod wine ...` in each
session unless you arrange to load it at startup.

## Load for one session

Without permanent installation, source both files in the current XSPEC
session:

```tcl
source /absolute/path/to/wine/procedures/wine_inclination.tcl
source /absolute/path/to/wine/procedures/wine_profile.tcl
```

This route needs no `tclIndex`. Source an edited file again to update its
definitions in an open session.

## Command reference

### wineinclcov

**NAME** — Local covariance uncertainty for the physical viewing angle.

**USAGE** — `wineinclcov N`

**DESCRIPTION** — Propagates XSPEC's full fit covariance through the WINE
viewing-angle expression. `N` is the XSPEC parameter number of `incl`,
`uincl`, or `uview`; opening angles may be free. The result is a local,
Gaussian 1-sigma uncertainty in degrees.

**OUTPUT** — Prints the angle and uncertainty to two decimal places.
`$xspec_tclout` and the Tcl return value contain the unrounded uncertainty
for scripts.

### wineinclerror

**NAME** — Convert native XSPEC `error` limits to viewing-angle degrees.

**USAGE** — `wineinclerror N ?delta-stat?`

**DESCRIPTION** — Calls XSPEC `error` for the selected `incl`, `uincl`, or
`uview` parameter. It profiles the physical angle only for direct `incl`, or
when both opening angles are frozen and unlinked. Unload XSPEC chains to get
fit-based rather than chain-based limits. The default `delta-stat` is `2.706`.

**OUTPUT** — Prints the converted limits. `$xspec_tclout` and the Tcl return
value hold the unrounded lower and upper limits followed by XSPEC's native
error string.

### wineinclprofile

**NAME** — Likelihood profile for the physical WINE viewing angle.

**USAGE** — `wineinclprofile ?delta-stat? N ?toler value? ?-fit-delta fraction?`

**DESCRIPTION** — Holds the physical angle fixed and refits the other free
parameters, including opening angles. `N` is the `uincl` or `uview` parameter
number. The optional leading real is the change in fit statistic, as in XSPEC
`error`: `wineinclprofile 1. N` requests Δstat = 1. Its default is `2.706`.
The calling XSPEC session remains unchanged on success or failure.

**PARAMETERS**

- `toler value` sets the absolute fit-statistic tolerance for quiet numerical
  crossings and for reporting a better minimum. Default: `1.0`. It does not
  change `delta-stat` or the inner XSPEC fit convergence.
- `-fit-delta fraction` changes the numerical derivative step in the child fit.
  Default: `0.001`. Use it when diagnosing convergence; XSPEC `error stopat`
  does not control this procedure's angle search.

**OUTPUT** — Prints the reference angle, confidence range in degrees, and
offsets from the reference, in a format similar to XSPEC `error`. The table is
rounded to two decimal places. `$xspec_tclout` and the Tcl return value contain
`lower upper flags` with unrounded limits and one four-letter diagnostic
string. The flag positions mean:

1. Lower allowed boundary reached before the requested Δstat.
2. Upper allowed boundary reached before the requested Δstat.
3. Lower crossing numerically unstable beyond `toler`.
4. Upper crossing numerically unstable beyond `toler`.

Each position is `T` when present and `F` otherwise. `FFFF` means neither
side reached an allowed boundary or exceeded the numerical tolerance. A
crossing accepted within `toler` is not flagged. For a boundary, the
corresponding reported limit is `0`, following XSPEC's pegged-limit display
convention; it is **not** a measured confidence limit. The actual allowed
boundary remains in `$wine_profile_sides`, and trial angles and statistics in
`$wine_profile_points`. A failed fit returns an error instead of limits.
These four positions are specific to WINE; XSPEC's native
[`tclout error`](https://heasarc.gsfc.nasa.gov/docs/software/xspec/manual/node86.html)
uses a different nine-position string.

If a trial improves the starting fit statistic by more than `toler`, profiling
stops and saves `better-fit.xcm` in the current working directory. Existing
files are preserved using numbered names. Inspect the candidate in a separate
XSPEC session before deciding whether to fit again or increase `toler`.

**EXAMPLES**

```tcl
wineinclcov 43
wineinclprofile 1. 43
puts $xspec_tclout
wineinclprofile 1. 43 toler 0.3
wineinclprofile -trial 36.71 43
```

The `-trial angle N` form runs one constrained diagnostic fit. Its full fit
statistic is placed in `$xspec_tclout`; the calling fit remains unchanged. On
failure, the error names a directory with the child XSPEC log and saved state.
See the [XSPEC guide](../docs/xspec.md#geometric-uncertainty-in-the-development-api) for
the geometry and model-link requirements.

The angular scan uses an isolated XSPEC session to preserve the caller's fit
state. Completed-point counts appear in the calling session; full worker
output stays in its log.

`windgeom`, `winderror`, and `windsteppar` print formatted output and store
the structured result in `::xspec_tclout`, returning an empty string so XSPEC
does not echo the dictionary. Read or copy this result before another
`tclout` or WINE command overwrites it. The result is cleared on entry so a
failed command cannot leave a previous successful result available.
