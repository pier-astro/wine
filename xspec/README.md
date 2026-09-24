# WINE XSPEC Models

XSPEC models for relativistic wind emission and absorption.

## Models

### windem - Emission Convolution
FFT-based relativistic emission line convolution with velocity and geometry.

**Parameters:**
1. `vout` (c): Wind velocity [0-0.99]. When `vout=0`, applies solid angle scaling only.
2. `incl` (deg): Viewing angle [0-90°]
3. `tout` (deg): Wind outer angle [0-90°], frozen by default
4. `tin` (deg): Wind inner angle [0-90°], frozen by default

**Constraint:** `tout > tin`

**Algorithm:** Conservative rebinning to a uniform log-energy grid followed by
a zero-padded linear FFT convolution. This preserves integrated additive flux
on uneven grids and prevents wrap-around between the energy boundaries.

**Complexity:** $O(n \log n)$

### windemlos - LoS-Constrained Emission Convolution

`windemlos` uses the same convolution as `windem`, but replaces `incl` with
the dimensionless parameter `uincl`:

```text
incl = tin + uincl * (tout - tin),    0 <= uincl <= 1
```

This parameterization guarantees `tin <= incl <= tout`, so it is appropriate
when the model assumes that the observed line of sight crosses the wind.
`uincl=0` and `uincl=1` select the inner and outer wind boundaries;
`uincl=0.5` selects their angular midpoint.

**Parameters:** `vout`, `uincl`, `tout`, `tin`. `uincl` is free with default
0.5; the remaining defaults match `windem`.

### windemoff - Off-Cone Emission Convolution

`windemoff` constrains the line of sight not to cross the wind. Because the
allowed inclination range has two disconnected parts, `region` selects one
part and remains frozen while `uview` varies continuously within it:

```text
region = 0:  incl = (1 - uview) * tin
region = 1:  incl = tout + uview * (90 deg - tout)
```

Here `0 <= uview <= 1`. In either region, `uview=0` is the grazing wind
boundary and `uview=1` is the point farthest from the wind (the polar axis or
the equatorial plane). `region=0` denotes the polar cavity and `region=1` the
region outside the outer wind boundary. The two regions are distinct
geometrical hypotheses and should be fitted separately rather than stepping
or fitting `region`.

**Parameters:** `vout`, `uview`, `tout`, `tin`, `region`. The default opening
angles are 30 and 60 degrees so both off-cone regions are non-degenerate.

### windline - Single Relativistic Line
Additive model for a single relativistic emission line with wind geometry.

**Parameters:**
1. `vout` (c): Wind velocity [0-0.99]. When `vout=0`, applies solid angle scaling only.
2. `incl` (deg): Viewing angle [0-90°]
3. `tout` (deg): Wind outer angle [0-90°], frozen by default
4. `tin` (deg): Wind inner angle [0-90°], frozen by default
5. `LineE` (keV): Rest-frame line energy
6. `norm`: XSPEC additive normalization (photons/cm²/s)

**Constraint:** `tout > tin`

**Complexity:** $O(n)$

### windlinelos - LoS-Constrained Single Line

`windlinelos` is the single-line counterpart of `windemlos`. It has parameters
`vout`, `uincl`, `tout`, `tin`, and `LineE`, followed by XSPEC's external
additive normalization. It uses the same `uincl` transformation and bounds.

### windlineoff - Off-Cone Single Line

`windlineoff` is the single-line counterpart of `windemoff`. Its parameters
are `vout`, `uview`, `tout`, `tin`, `region`, and `LineE`, followed by XSPEC's
external additive normalization. It uses the same two `region` definitions.

### wind - Absorption Velocity Shift
Relativistic velocity shift for multiplicative transmission spectra.

**Parameters:**
1. `vout` (c): Wind velocity [0-0.99]

**Physics:** Positive `vout` denotes an approaching outflow and blueshifts
features: E_shifted = E × sqrt((1+β)/(1-β)).

**Scope:** Apply only to multiplicative components such as `mtable`. Additive
spectra require `awind` or `azwind`.

**Complexity:** $O(n)$

### zwind - Absorption with Redshift
Velocity shift including cosmological redshift for multiplicative spectra.

**Parameters:**
1. `vout` (c): Wind velocity [0-0.99]
2. `z`: Redshift [0-10]

**Physics:** Combined factor: sqrt((1+β)/(1-β)) / (1+z)

**Complexity:** $O(n)$

### awind / azwind - Additive Velocity Shifts
Flux-conserving relativistic shifts for additive spectra. `awind` takes
`vout` in units of `v/c`; `azwind` also takes `z` and includes XSPEC's
cosmological time-dilation factor. They match the numerical roles of
`vashift` and `zashift`, but use the exact relativistic velocity factor.

Use these around additive components (`atable`, continuum, lines, or additive
sums), never around multiplicative transmission spectra.

## Installation

```bash
cd xspec
./build.sh
```

Or manually:
```bash
initpackage wine lmodel.dat .
hmake
```

## Usage in XSPEC

```tcl
# Load package
lmod wine /path/to/xspec

# Emission convolution
model windem*atable{emission.fits}
0.3     # vout = 0.3c
/*      # Accept defaults

# Emission convolution with tin <= incl <= tout by construction
model windemlos*atable{emission.fits}
0.3     # vout = 0.3c
0.5     # uincl: angular midpoint of the wind
/*      # Accept angle defaults

# Off-cone emission through the polar cavity (region=0)
model windemoff*atable{emission.fits}
0.3     # vout = 0.3c
0.5     # uview: halfway from the wind boundary to the pole
60      # tout
30      # tin
0       # region: polar cavity; keep frozen

# Single line profile
model windline
0.1 45 90 0 6.5 1.0

# LoS-constrained single line profile
model windlinelos
0.1 0.5 90 0 6.5 1.0

# Off-cone single line through the outer region
model windlineoff
0.1 0.5 60 30 1 6.5 1.0

# Absorption with velocity shift
model (wind*mtable{absorption.fits})*powerlaw

# Absorption with velocity and redshift
model (zwind*mtable{absorption.fits})*powerlaw
```
> **NOTE:** Remember to use the parentheses otherwise the whole continuum+lines will be shifted!

### Energy coverage for `wind` and `zwind`

These are convolution models: XSPEC must calculate the table over every input
energy needed after the Doppler shift. With no `energies` command, XSPEC uses
each response's native energy grid. The `.1--10 keV, 1000 linear bins` values
in the XSPEC manual are defaults for an *invoked* `energies` command with
omitted fields, not the active grid in that case.

For `wind`, positive `beta` moves a rest-frame feature to
`E_out = B E_rest`, where `B = sqrt((1+beta)/(1-beta))`. Evaluating the
detector band therefore requires the model below its low-energy boundary.
After loading all data, for example:

```tcl
energies extend low 0.01 1000 log
```

Use both low and high extensions for `zwind` when its combined Doppler and
cosmological factor requires them. Choose bounds within the FITS table's
energy coverage. Sparse UV/hard-X-ray padding is useful for this purpose, but
the fine core grid is only used if XSPEC is instructed to evaluate there.
`wind` and `zwind` issue a one-time reminder to extend the model grid beyond
the detector response. Like XSPEC's `vmshift`, they interpolate multiplicative
factors directly and therefore behave correctly on uneven XSPEC energy grids.

## Requirements

- XSPEC 12+
- FFTW3 library (for windem FFT convolution)
