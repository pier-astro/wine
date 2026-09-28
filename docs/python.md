# Python models and spectral inputs

The Python package applies WINE's kinematic and geometric calculations to rest-frame spectra prepared elsewhere. `src/wine/kernels.py` contains the analytical line profile and convolution; `src/wine/models.py` combines those calculations with absorption or emission spectra; `src/wine/grid.py` interpolates the supplied spectra. The package does not call an ionisation engine.

## Spectral grid contract

`wine.SpectralGrid` accepts increasing one-dimensional axes `energy` (keV), `logxi`, `vturb` (km/s), and `log_nh = log10(N_H / cm^-2)`. Each spectral array has shape `(n_logxi, n_vturb, n_log_nh, n_energy)`. Optical-depth arrays are dimensionless and nonnegative. Emission arrays contain nonnegative **bin-integrated** photon spectra; their physical normalisation is supplied by the grid producer, and `emission_slab(..., norm=...)` scales it.

The supported dictionary keys are set by the caller. The slab functions use optical-depth modes `absorption`, `scattering`, and `total`, and emission modes `total`, `lines`, and `continuum`. Only modes present in the grid can be evaluated.

```python
import numpy as np
from wine import SpectralGrid, absorption_slab

energy = np.geomspace(0.3, 20.0, 1000)  # rest-frame bin centres, keV
tau = np.full((2, 1, 2, energy.size), 0.1)
grid = SpectralGrid(
    energy,
    logxi=[3.0, 4.0],
    vturb=[100.0],
    log_nh=[22.0, 23.0],
    optical_depth={"absorption": tau},
)
# Direct grid interpolation; no wind model is needed.
rest_frame_tau = grid.optical_depth(logxi=3.5, vturb=100.0, log_nh=22.5)
observed_energy = np.geomspace(0.5, 10.0, 500)
transmission = absorption_slab(
    observed_energy, grid,
    logxi=3.5, log_nh=22.5, vturb=100.0,
    beta=0.2, C_f=1.0,
)
```

`SpectralGrid.from_tables(...)` provides the same interface for matching OGIP XSPEC tables:

```python
grid = SpectralGrid.from_tables(
    optical_depth={"absorption": "optdepth.emod"},
    emission={"lines": "linem.amod"},
)
rest_frame_tau = grid.optical_depth(logxi=3.5, vturb=100.0, log_nh=22.5)
rest_frame_lines = grid.emitted(logxi=3.5, vturb=100.0, log_nh=22.5, mode="lines")
```

These file names match the bundled `iongrid` recipe but are placeholders; no production tables ship with WINE. Optical-depth tables may be exponential tables containing optical depth or multiplicative tables containing transmission. Emission tables must be additive and contain bin-integrated photons. All supplied tables must have identical energy bins, parameter axes, and parameter-coordinate interpolation methods. Recognised parameter names are `logxi`, `vturb`, and `lognH`/`log_nh` (case-insensitive). Supply each omitted singleton axis through `fixed_parameters`; omit `fixed_parameters` if all three axes vary in the tables. Reading FITS files requires the `fits` extra.

The general `wine.ogip.read_table_model(path)` reader can also evaluate one OGIP table directly. It supports regular interpolation grids and rejects additional-parameter spectra and XSPEC filter expressions.

## Interpolation

`SpectralGrid` defaults to `interpolation="log"`: it interpolates the logarithm of each spectrum and replaces exact zeros by the configurable positive `floor` (default `1e-30`). Use `interpolation="linear"` when exact zeros or arithmetic interpolation are required. This choice affects **spectral values**, not the grid coordinates.

Parameter coordinates follow OGIP `METHOD=0` (linear) or `METHOD=1` (logarithmic) when loaded from tables. The spatial interpolator defaults to linear if any parameter axis has fewer than four points, otherwise PCHIP; `method="linear"` or `method="pchip"` can be selected explicitly. The direct OGIP reader's `evaluate()` always interpolates spectral values linearly, as XSPEC does. Thus the default Python `SpectralGrid` and XSPEC may differ between table nodes even when both read the same files.

## Slab functions

| Function | Input spectrum | Output | Main controls |
| --- | --- | --- | --- |
| `absorption_slab(E_obs, grid, ...)` | Rest-frame optical depth | Transmission on observer-frame energy centres | `logxi`, `log_nh`, `vturb`, line-of-sight `beta`, covering fraction `C_f`, redshift `z`, optical-depth `mode` |
| `emission_slab(E_obs, grid, ...)` | Rest-frame bin-integrated emission | Observer-frame photon-flux density, scaled by `norm` | Grid coordinates, outflow `beta`, `incl`, `theta_in`, `theta_out`, `norm`, `z`, convolution `method`, emission `mode` |

`E_obs` is an increasing array of energy centres in keV. The emission angles are **radians** in Python; the XSPEC local models use **degrees**. Positive absorption `beta=v/c` denotes motion toward the observer and blueshifts features. `C_f` lies between zero and one. The slab functions clamp grid coordinates to the supplied tabulated ranges; validate the scientific suitability of any clamped evaluation in your analysis.

For emission, `norm` is a simple multiplicative factor. WINE does not infer a source luminosity, distance, or shell covering factor from the FITS header. The returned flux density inherits the grid producer's photon normalisation divided by energy-bin width and multiplied by `norm`. For a consistent joint model, emission and absorption spectra must use compatible physical conditions and conventions.

## Current limits

WINE 2.0 currently evaluates individual slabs or shells. It does not calculate ionisation balance, attenuation of one shell by another, or a self-consistent multi-shell radiation field. The approximate Compton-thin focus in the README is a modelling scope, not a hard-coded column-density cut. The supplied grid's energy coverage and resolution must support the requested Doppler shifts and observer band.
