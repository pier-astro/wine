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

`build.sh` removes old generated files, runs `initpackage redwine lmodel.dat .`, and compiles the package with `hmake`. To remove generated files and rebuild:

```bash
./clean.sh
./build.sh
```

The scripts do not install HEASoft or modify a system installation. Generated Makefiles, XSPEC glue code, object files, and shared libraries are ignored by Git; the C++ sources, `lmodel.dat`, scripts, and tests are tracked. Build products remain available in your working tree until you run `clean.sh`.

Load the resulting package in XSPEC using an absolute path to the `xspec` directory:

```tcl
lmod redwine /absolute/path/to/redwine/xspec
```

## Model families

| Models | XSPEC type | Purpose |
| --- | --- | --- |
| `wshift`, `zwshift` | Convolution | Relativistic shift of a multiplicative transmission spectrum, with optional cosmological redshift |
| `awshift`, `azwshift` | Convolution | Relativistic shift of an additive spectrum, with optional cosmological redshift |
| `wem`, `wemlos`, `wemoff` and `wem_Cf`, `wemlos_Cf`, `wemoff_Cf` | Convolution | Relativistic emission broadening for a supplied additive spectrum |
| `wline`, `wlinelos`, `wlineoff` and `wline_Cf`, `wlinelos_Cf`, `wlineoff_Cf` | Additive | Single relativistic emission-line profile |

The models with suffix `_Cf` take covering fraction and a shape coordinate.
The other emission models take `tin` and `tout` directly. Both evaluate the
same physical profile. See the [model reference](../xspec/README.md).

## Geometric uncertainty

`windgeom` reports physical geometry and local covariance uncertainties.
`windsteppar` scans the two opening angles in an isolated XSPEC session.
`winderror` profiles physical opening angles and supported inclination cases.
See the [procedure guide](../procedures/README.md) for constraints.

For posterior summaries, transform each joint sample with `Geometry.free`,
`Geometry.los`, or `Geometry.off` and compute quantiles of the physical
angles. Flat priors in angular and covering-fraction coordinates generally yield different physical priors.

These examples use placeholder tables generated elsewhere:

```tcl
lmod redwine /absolute/path/to/redwine/xspec
model wemlos_Cf*atable{emission.fits}
model (wshift*mtable{absorption.fits})*powerlaw
```

Velocity shifts and emission broadening need input energies beyond the fitted
band. Extend the XSPEC calculation grid over the shifts allowed by the fit,
within the supplied tables' energy coverage. The
[energy-coverage guide](../xspec/README.md#energy-coverage-for-shifts-and-emission-broadening)
explains the bounds and gives an `energies extend` example. RED-WINE's notices are
reminders, not checks that a grid is sufficient. Verify stability in the fitted
band after extending the calculation range and checking its sampling.
