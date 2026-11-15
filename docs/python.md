# Python models and spectral inputs

The Python package applies RED-WINE's kinematic and geometric calculations to rest-frame spectra prepared elsewhere. `src/redwine/kernels.py` contains the analytical line profile and convolution; `src/redwine/models.py` combines those calculations with absorption or emission spectra; `src/redwine/grid.py` interpolates the supplied spectra. The package does not call an ionisation engine.

## Spectral grid contract

`redwine.SpectralGrid` accepts increasing one-dimensional axes `energy` (keV), `logxi`, `vturb` (km/s), and `log_nh = log10(N_H / cm^-2)`. Each spectral array has shape `(n_logxi, n_vturb, n_log_nh, n_energy)`. Optical-depth arrays are dimensionless and nonnegative. Emission arrays contain nonnegative **bin-integrated** photon spectra; their physical normalisation is supplied by the grid producer, and `emission_slab(..., norm=...)` scales it.

The supported dictionary keys are set by the caller. The slab functions use optical-depth modes `absorption`, `scattering`, and `total`, and emission modes `total`, `lines`, and `continuum`. Only modes present in the grid can be evaluated.

```python
import numpy as np
from redwine import SpectralGrid, absorption_slab

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

These file names match the bundled `iongrid` recipe but are placeholders; no production tables ship with RED-WINE. Optical-depth tables may be exponential tables containing optical depth or multiplicative tables containing transmission. Emission tables must be additive and contain bin-integrated photons. All supplied tables must have identical energy bins, parameter axes, and parameter-coordinate interpolation methods. Recognised parameter names are `logxi`, `vturb`, and `lognH`/`log_nh` (case-insensitive). Supply each omitted singleton axis through `fixed_parameters`; omit `fixed_parameters` if all three axes vary in the tables. Reading FITS files requires the `fits` extra.

The general `redwine.ogip.read_table_model(path)` reader can also evaluate one OGIP table directly. It supports regular interpolation grids and rejects additional-parameter spectra and XSPEC filter expressions.

## Interpolation

`SpectralGrid` defaults to `interpolation="log"`: it interpolates the logarithm of each spectrum and replaces exact zeros by the configurable positive `floor` (default `1e-30`). Use `interpolation="linear"` when exact zeros or arithmetic interpolation are required. This choice affects **spectral values**, not the grid coordinates.

Parameter coordinates follow OGIP `METHOD=0` (linear) or `METHOD=1` (logarithmic) when loaded from tables. The spatial interpolator defaults to linear if any parameter axis has fewer than four points, otherwise PCHIP; `method="linear"` or `method="pchip"` can be selected explicitly. The direct OGIP reader's `evaluate()` always interpolates spectral values linearly, as XSPEC does. Thus the default Python `SpectralGrid` and XSPEC may differ between table nodes even when both read the same files.

## Slab functions

| Function | Input spectrum | Output | Main controls |
| --- | --- | --- | --- |
| `absorption_slab(E_obs, grid, ...)` | Rest-frame optical depth | Transmission on observer-frame energy centres | `logxi`, `log_nh`, `vturb`, line-of-sight `beta`, covering fraction `C_f`, redshift `z`, optical-depth `mode` |
| `emission_slab(E_obs, grid, ...)` | Rest-frame bin-integrated emission | Observer-frame photon-flux density, scaled by `norm` | Grid coordinates, outflow `beta`, `geometry`, `norm`, `z`, convolution `method`, emission `mode` |

`E_obs` is an increasing array of energy centres in keV. Public `Geometry`
angles and XSPEC inclination parameters are in **degrees**; conversion to
radians occurs at the numerical-kernel boundary. The low-level functions in
`redwine.kernels` still take radians. Positive absorption `beta=v/c` denotes motion
toward the observer and blueshifts features. Absorption `C_f` lies between zero
and one and is distinct from the wind's hemispheric `Geometry.Cf`. The slab
functions clamp grid coordinates to supplied table ranges; validate the
scientific suitability of any clamped evaluation in your analysis.

## Physical geometry

```python
from redwine import Geometry

explicit = Geometry.angles(incl=45, tin=30, tout=60)
generic = Geometry.free(Cf=0.3, hole=0.5, incl=45)
through_wind = Geometry.los(Cf=0.3, skew=0.5, incl=45)
outside = Geometry.off(Cf=0.3, hole=0.5, view=0.5, region="polar")
print(through_wind.incl, through_wind.tin, through_wind.tout, through_wind.Cf)
```

These constructors produce the same physical representation and perform no
spectral calculation. `Cf=cos(tin)-cos(tout)` is relative to one hemisphere;
the spectral full-sphere angular scaling remains `Cf/2`. They accept broadcast
arrays for converting joint chain samples or scans. A single spectrum requires
scalar geometry:

```python
flux = emission_slab(
    E_obs, grid, logxi=4.0, log_nh=23.0, vturb=2000.0,
    beta=0.25, geometry=through_wind, norm=1.0,
)
```

Fractions must be finite and inside [0, 1]; public angles must be inside
0--90 degrees. Inverted explicit boundaries warn and represent zero emission.
Zero coverage and collapsed off sectors are limiting geometries, not measured
angles. The [XSPEC model reference](../xspec/README.md) explains `hole`, `skew`,
and `view`, which have distinct meanings. Joint posterior samples retain
their correlations when transformed; a converted native-coordinate scan is
not generally a profile of derived axes that also depend on refitted nuisance
parameters.

For emission, `norm` is a simple multiplicative factor. RED-WINE does not infer a source luminosity, distance, or shell covering factor from the FITS header. The returned flux density inherits the grid producer's photon normalisation divided by energy-bin width and multiplied by `norm`. For a consistent joint model, emission and absorption spectra must use compatible physical conditions and conventions.

## Current limits

RED-WINE 2.0 currently evaluates individual slabs or shells. It does not calculate ionisation balance, attenuation of one shell by another, or a self-consistent multi-shell radiation field. The approximate Compton-thin focus in the README is a modelling scope, not a hard-coded column-density cut. The supplied grid's energy coverage and resolution must support the requested Doppler shifts and observer band.

## Experimental velocity distributions

`redwine.experimental.kernels` evaluates radial-velocity distributions for
uniform intrinsic emissivity per unit solid angle. For example:

```python
import numpy as np
from redwine.experimental.kernels import Wind, line

edges = np.linspace(5, 10, 5001)  # keV
wind = Wind("gaussian", center=0.25, spread=0.03,
            inner=20, outer=50, inclination=45)
photon_density = line(edges, 7.0, wind) / np.diff(edges)
```

`Wind` uses degrees for angles and units of c for velocities. For
`"uniform"`, `spread` is the full interval width and `center` is a speed
magnitude. Set `inflow=True` to reverse it. For `"gaussian"`, `spread` is
the parent standard deviation and `center` is signed. The Gaussian is
truncated at ±3 parent sigmas and at the physical limits −1 < β < 1.
Its PDF is normalized over the retained interval. Set `bicone=True` to include the
antipodal cone.

`line` returns bin-integrated photons per unit intrinsic all-sky line flux.
`convolve(edges, flux, wind)` redistributes bin-integrated photons assuming
constant input photon density within each bin. These kernels retain the
geometric covering fraction and the moving-emitter Doppler weight δ³ for
integrated photons. The original `redwine.kernels.lineshape` uses δ² for
integrated photons. The uniform kernel integrates constant azimuth panels
analytically. Both kernels use one-dimensional velocity quadrature on the
remaining panels and share their implementation with the XSPEC windkin models.

The three `Slider Emission Profiles` notebooks in `example/` run with
`matplotlib`, `ipympl`, `ipywidgets` and `astropy` in the notebook
environment. Run all cells to initialize the sliders.
