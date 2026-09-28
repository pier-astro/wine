# WINE XSPEC procedures

`wine_inclination.tcl` provides `wineinclcov` (local full-covariance error)
and `wineinclerror` (a wrapper around native XSPEC `error` for direct `incl`
or fixed opening angles). `wine_profile.tcl` provides `wineinclprofile`, the
XSPEC `error` counterpart for the physical viewing angle when angles are free. It
refits those angles and other nuisance parameters on the existing `los` and
`off` models.

These commands require XSPEC and the local WINE model build, but not xspeclab.
The generic xspeclab `coverror` command remains useful for derived quantities
outside WINE. See the [XSPEC guide](../docs/xspec.md#viewing-angle-uncertainty-in-xspec)
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
and link both WINE Tcl files into that directory instead. Run `auto_mkindex .`
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
See the [XSPEC guide](../docs/xspec.md#viewing-angle-uncertainty-in-xspec) for
the geometry and model-link requirements.
