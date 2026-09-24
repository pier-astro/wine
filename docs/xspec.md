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

These examples use placeholder tables generated elsewhere:

```tcl
lmod wine /absolute/path/to/wine/xspec
model windemlos*atable{emission.fits}
model (wind*mtable{absorption.fits})*powerlaw
```

`wind` and `zwind` require the table to cover the shifted energies evaluated by XSPEC. Extend the XSPEC model grid as needed and keep it within the table's energy coverage; the [model reference](../xspec/README.md#energy-coverage-for-wind-and-zwind) gives an example.
