# WINE v2

Relativistic wind profiles, engine-independent Python spectral models (`import wine`), and local XSPEC models. WINE consumes precomputed spectra; ionization-engine runs and OGIP table production belong in `iongrid`.

## Python

From this checkout, install with `/usr/local/bin/python -m pip install -e '.[fits]'` in your chosen environment. `astropy` is needed only to read OGIP FITS tables. The model layer needs NumPy and SciPy, not CLOUDY or XSTAR.

`wine.SpectralGrid` accepts arrays of rest-frame optical depth and bin-integrated emission on `(logxi, vturb, log_nh, energy)` axes. `wine.SpectralGrid.from_tables(optical_depth={"absorption": "tau.fits"}, emission={"total": "emission.fits"})` reads matching OGIP exponential/multiplicative and additive tables. For tables that omit a singleton axis, pass its value in `fixed_parameters`, for example `{"vturb": 100.0}`. The grid exposes `optical_depth()` and `emitted()` for `absorption_slab()` and `emission_slab()`.

Spectral interpolation defaults to `interpolation="log"`: it interpolates log10 of nonnegative values, replacing exact zeros with `1e-30`. Select `interpolation="linear"` to retain exact zeros. The parameter interpolation method is linear for small grids and PCHIP when every axis has at least four points, matching the former Python `PlasmaGrid` choice. OGIP `METHOD=1` logarithmic parameter coordinates are respected. The standalone `wine.xspec.read_table_model(...).evaluate()` retains XSPEC's linear interpolation of table *spectra*; it is useful when reproducing XSPEC exactly. Consequently, WINE's default Python grid and XSPEC tables need not agree between grid nodes.

## Local XSPEC build

The canonical C++ sources are in `xwine/`. Initialize HEASoft/XSPEC on your laptop, then run:

```bash
cd xwine
./build.sh
/usr/local/bin/python test_geometry_models.py
./clean.sh
./build.sh
```

`build.sh` always starts from clean generated files. `clean.sh` removes only generated build products. Compiled libraries, objects, generated XSPEC files and local Python environments are ignored by Git; source files and build scripts are tracked. Load the built models in XSPEC with `lmod wine /absolute/path/to/wine-v2/xwine`.

## Provenance

Initial Python kernels, slab models, OGIP reader and local XSPEC sources were extracted from `pier-astro/xflow` at commit `e55c58b`. The legacy `WINE/xwine` implementation was not copied. GPL-3.0-or-later licensing is retained.
