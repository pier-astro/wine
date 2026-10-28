# XSPEC models and local build

`xspec/` contains the local XSPEC implementation. `lmodel.dat` is the source of truth for registered models and parameter bounds. The [model reference](../xspec/README.md) lists every parameter, geometry transform, and example expression.

## Build on your machine

Install and initialize HEASoft/XSPEC, including its `initpackage` and `hmake` commands, and make FFTW3 available to the compiler. From a clone of this repository:

```bash
cd xspec
./build.sh
python test_geometry_models.py
python test_energy_coverage.py
python test_geometry_fit.py
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

The geometry API is `vout, Cf, hole, incl` for generic emission,
`vout, Cf, skew, incl` for LoS emission, and `vout, Cf, hole, view, region`
for off emission. Physical angles are constructed before the shared kernel.
The [model reference](../xspec/README.md) defines the geometry and limits.

## Geometric uncertainty in the development API

`Cf` and generic/LoS `incl` are ordinary physical XSPEC fit parameters; use
native `error` or joint chains for them. `windgeom` reports physical geometry
by component, and `windsteppar` profiles the two opening angles in an isolated
session. See the [procedure guide](../procedures/README.md) for usage and
limitations. `winderror` profiles derived opening angles and supported
inclination cases; `windgeom` also reports full-covariance local uncertainties
when available. Free-geometry off-wind inclination profiles remain unsupported.
The existing `wineincl*` procedures
and their installation guide refer to the retired `uincl`/`uview` coordinates;
do not use them to interpret the new parameters. Their tested search and
`toler` behavior will be retained when adapting the constraints.

For posterior summaries, transform each joint sample with `Geometry.free`,
`Geometry.los`, or `Geometry.off` and compute quantiles of the physical
angles. A flat prior in the new coordinates is not generally the same prior
as a flat distribution in the former coordinates.

These examples use placeholder tables generated elsewhere:

```tcl
lmod wine /absolute/path/to/wine/xspec
model windemlos*atable{emission.fits}
model (wind*mtable{absorption.fits})*powerlaw
```

Velocity shifts and emission broadening need input energies beyond the fitted
band. Extend the XSPEC calculation grid over the shifts allowed by the fit,
within the supplied tables' energy coverage. The
[energy-coverage guide](../xspec/README.md#energy-coverage-for-shifts-and-emission-broadening)
explains the bounds and gives an `energies extend` example. WINE's notices are
reminders, not checks that a grid is sufficient. Verify stability in the fitted
band after extending the calculation range and checking its sampling.
