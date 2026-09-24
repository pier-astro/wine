# WINE — Wind in the Ionised Nuclear Environment

WINE is a spectroscopic model for emission and absorption from outflows around compact sources. It combines fast analytical wind profiles with accurate special relativistic effects and viewing geometry, treating both spectral components consistently. It is designed for outflows up to mildly relativistic speeds, including ultra-fast outflows in active galactic nuclei (AGN), and can be used to explore winds in X-ray binaries and gamma-ray bursts when its assumptions apply.

## WINE 2.0

WINE 2.0 separates the wind model from the calculation of its rest-frame spectrum. It acts as a post-processing layer: users can supply their own spectra or generate them with photoionisation codes such as CLOUDY or XSTAR, then apply WINE's relativistic emission and absorption models. This makes it practical to explore individual wind shells and their geometry without rerunning a photoionisation calculation at every model evaluation. The present shell-based approach is aimed primarily at approximately Compton-thin winds (roughly $N_\mathrm{H} \lesssim 10^{24}\,\mathrm{cm}^{-2}$); its validity still depends on the supplied spectra and the physical assumptions of a particular application.

This direction is motivated in part by the detail now accessible in high-resolution X-ray spectra from XRISM.

The [`wine v1` implementation](https://baltig.infn.it/ionisation/wine) remains available as a complementary tool. It models a connected, initially homogeneous outflow through successive XSTAR calculations, retaining the coupled wind treatment useful for that approach.

## What is in this repository?

| Path | Contents |
| --- | --- |
| [`src/wine/`](src/wine/) | Python wind kernels, slab models, spectral grids, and OGIP table reader |
| [`xspec/`](xspec/) | Local XSPEC models for emission profiles and relativistic spectral shifts |
| [`tests/`](tests/) | Python tests for grids, table reading, and model evaluation |
| [`docs/`](docs/) | API, spectral input, XSPEC, and development documentation |

WINE reads precomputed spectra; it does not run CLOUDY or XSTAR. No production spectral grid is bundled with this repository.

## Getting started

For Python 3.10 or newer, install in an environment of your choice:

```bash
git clone https://github.com/pier-astro/wine.git
cd wine
python -m pip install -e '.[fits]'
```

The `fits` extra provides Astropy for OGIP FITS tables. NumPy and SciPy are the core Python dependencies. Import the model with `import wine`; see the [Python guide](docs/python.md) for array and table inputs, units, and examples.

For XSPEC, initialize a local HEASoft/XSPEC installation, then build and load the models:

```bash
cd xspec
./build.sh
# In XSPEC: lmod wine /absolute/path/to/wine/xspec
```

The [XSPEC guide](docs/xspec.md) describes the models, tests, and `./clean.sh` rebuild workflow. Generated binaries and local build files are excluded from Git.

## Documentation

- [Python models and spectral inputs](docs/python.md)
- [XSPEC models and local build](docs/xspec.md)
- [Repository structure and development](docs/development.md)

The physical model and its domain of validity will be documented in greater depth alongside the associated publication. The earlier WINE model is described in [Luminari et al. (2024)](https://arxiv.org/abs/2410.13933).

WINE 2.0 is distributed under [GPL-3.0-or-later](LICENSE).
