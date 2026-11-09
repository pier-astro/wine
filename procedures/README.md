# RED-WINE XSPEC procedures

The same procedures support the direct-angle and `_Cf` emission models.

## Load automatically in XSPEC

Generate the index from a terminal:

```sh
cd /absolute/path/to/red-wine/procedures
tclsh <<'EOF'
auto_mkindex . wind_geometry.tcl
EOF
```

Add the directory to Tcl's search path in `~/.xspec/xspec.rc`, then restart XSPEC:

```tcl
set auto_path [linsert $auto_path 0 /absolute/path/to/red-wine/procedures]
```

Keep the generated `tclIndex` alongside the scripts. Git ignores it. Regenerate
it when a command or script is added, removed, or renamed.

Keep `lmod redwine /absolute/path/to/red-wine/xspec` in the startup file
to load the models. For a single session, you can source the procedure file directly:

```tcl
source /absolute/path/to/red-wine/procedures/wind_geometry.tcl
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
physical quantities `tin`, `tout`, and supported `incl` cases below.

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

- Direct-angle models constrain the other edge when profiling an opening angle.
  For physical inclination, the LoS and off models link the free relative
  viewing coordinate to the target inclination as the edges vary.
- For `_Cf` models with free `Cf`, native links impose the selected opening angle. LoS geometry
  restricts inclination to the compatible side of that angle; generic/off
  geometry restricts `hole` instead. The eliminated `Cf` needs its full hard
  range 0--1. Fixed endpoint `skew` charts are explicitly rejected rather than
  regularized by arbitrary denominator floors.
- With fixed `Cf`, the requested angle also fixes the other opening angle;
  the existing two-angle constraint is reused. With fixed `Cf` and `skew`,
  LoS angular limits can instead be transformed monotonically from native
  inclination limits. Fixed constraints are never silently released.
- In the `_Cf` family, generic/LoS `incl` uses native XSPEC `error`.
  Off-wind `incl` is supported
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
into the four RED-WINE flags. These checks do not establish global convergence.

See [the interface and validation status](../docs/geometry-procedures.md).

For the direct-angle models, `windgeom` reports the derived hemispheric covering
fraction and propagates covariance in the native coordinates. In angular scans,
a cell is excluded when its requested angles violate the ordering or disagree
with a frozen opening angle. The angular inclination profiler requires a free
`uincl` or `uview` when using the relative viewing coordinates.
