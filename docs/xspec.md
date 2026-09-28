# XSPEC models and local build

`xspec/` contains the local XSPEC implementation. `lmodel.dat` is the source of truth for registered models and parameter bounds. The [model reference](../xspec/README.md) lists every parameter, geometry transform, and example expression.

## Build on your machine

Install and initialize HEASoft/XSPEC, including its `initpackage` and `hmake` commands, and make FFTW3 available to the compiler. From a clone of this repository:

```bash
cd xspec
./build.sh
python test_geometry_models.py
```

`build.sh` removes old generated files, runs `initpackage wine lmodel.dat .`, and compiles the package with `hmake`. To remove generated files and rebuild:

```bash
./clean.sh
./build.sh
```

The scripts do not install HEASoft or modify a system installation. Generated Makefiles, XSPEC glue code, object files, and shared libraries are ignored by Git; the C++ sources, `lmodel.dat`, scripts, and tests are tracked. Build products remain available in your working tree until you run `clean.sh`.

Load the resulting package in XSPEC using an absolute path to the `xspec` directory:

```tcl
lmod wine /absolute/path/to/wine/xspec
```

## Model families

| Models | XSPEC type | Purpose |
| --- | --- | --- |
| `wind`, `zwind` | Convolution | Relativistic shift of a multiplicative transmission spectrum, with optional cosmological redshift |
| `awind`, `azwind` | Convolution | Relativistic shift of an additive spectrum, with optional cosmological redshift |
| `windem`, `windemlos`, `windemoff` | Convolution | Relativistic emission broadening for a supplied additive spectrum |
| `windline`, `windlinelos`, `windlineoff` | Additive | Single relativistic emission-line profile |

The `los` models map `uincl` in `[0, 1]` into `tin <= incl <= tout`. The `off` models map `uview` into one of two off-cone regions selected by a frozen `region` switch. At `uview=0`, the sightline lies exactly on a wind boundary. The two regions should be fitted as separate geometrical cases. See the [model reference](../xspec/README.md) for parameter order and expressions.

## Viewing-angle uncertainty in XSPEC

The viewing angle is derived from `uincl` or `uview`. Profiling either fraction with XSPEC's `error` does not generally profile the angle when `tin` or `tout` is free. The procedures in [`procedures/`](../procedures/README.md) provide two complementary estimates:

```tcl
source /absolute/path/to/wine-v2/procedures/wine_inclination.tcl
source /absolute/path/to/wine-v2/procedures/wine_profile.tcl
fit
wineinclcov 2             ;# local Gaussian 1-sigma error, degrees
wineinclprofile 2          ;# profile limits, delta-stat = 2.706
wineinclprofile 1.0 2      ;# profile limits, delta-stat = 1
puts $xspec_tclout         ;# unrounded lower upper and diagnostic flags
```

Here `2` is the XSPEC parameter number of `uincl` or `uview`; find it with `show parameters`. As in XSPEC's `error`, an optional leading real number is the **change in fit statistic**, not a sigma value; write `1.0` or `1.` to distinguish it from an integer parameter number. For a one-parameter Gaussian approximation, 1.3 sigma corresponds to a change of `1.69`. `wineinclcov` follows identity ties and propagates the full XSPEC fit covariance with analytic derivatives, including free opening angles. It needs no xspeclab installation.

`wineinclprofile` is the direct analogue of XSPEC `error` for the **physical viewing angle**: it holds that angle fixed while refitting the free opening angles and other nuisance parameters. `wineinclerror` instead calls native XSPEC `error` on an independent `incl`, `uincl`, or `uview` parameter and converts its limits to degrees. That conversion gives the physical-angle profile only when both opening angles are frozen and unlinked, because the transformation is then one-to-one with fixed coefficients. It is a convenient special case, not a substitute for `wineinclprofile` with free opening angles. Unload XSPEC chains before using `wineinclerror` for a likelihood profile: native `error` returns chain intervals when chains are loaded.

`wineinclprofile` keeps the original XSPEC session intact, including the model expression, links, fit and covariance. It uses the current fit method to refit the other parameters at each fixed physical angle. If a trial finds a substantially better fit, the procedure stops and saves `better-fit.xcm` in the current working directory for inspection; an existing file receives a numbered successor. The optional `toler value` controls the absolute fit-statistic difference used for that decision and for numerical crossing warnings (default `1.0`). It is separate from the requested Δstat. See the [procedure reference](../procedures/README.md#wineinclprofile) for diagnostics and advanced settings.

With MCMC, there is no need to hold the physical angle fixed during sampling: each chain row is a joint draw of the free parameters. Derive `incl` separately for every row using its `uincl` or `uview` and opening-angle values, then calculate quantiles of the derived-angle samples. Use fixed values for frozen angles and the sampled master values for linked angles. This retains parameter correlations even when `tin` and `tout` vary; transforming marginal `u` error bars alone does not. These are posterior credible intervals, not the likelihood-profile limits returned by `wineinclprofile`.

The printed result resembles XSPEC `error`: a reference angle, lower and upper confidence values, and their offsets. `$xspec_tclout` holds unrounded `lower upper flags`; the four `T`/`F` flags record lower or upper allowed boundaries and lower or upper numerical instability, in that order. `FFFF` indicates neither condition. A crossing accepted within `toler` needs no warning flag. A reached boundary is printed as `0`, following XSPEC's pegged-limit convention, and is not a measured confidence limit. The [procedure reference](../procedures/README.md#wineinclprofile) documents the flag positions and diagnostic variables. For a failed trial, run `wineinclprofile -trial angle N` to check that angle alone. Geometry parameters (`u`, `tin`, `tout`) may share identity links across matching components; more complex links on those parameters are unsupported. To share the *physical inclination*, link all three geometric parameters across components, not only `u`. Other parameter links are preserved. The separate process must be able to reload the fit's data and table files.

For permanent availability, see the [WINE procedure installation guide](../procedures/README.md). The default `2.706` is the usual 90% one-parameter change in fit statistic; choose a different threshold explicitly when needed.

These examples use placeholder tables generated elsewhere:

```tcl
lmod wine /absolute/path/to/wine/xspec
model windemlos*atable{emission.fits}
model (wind*mtable{absorption.fits})*powerlaw
```

`wind` and `zwind` require the table to cover the shifted energies evaluated by XSPEC. Extend the XSPEC model grid as needed and keep it within the table's energy coverage; the [model reference](../xspec/README.md#energy-coverage-for-wind-and-zwind) gives an example.
