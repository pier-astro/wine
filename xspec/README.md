# WINE XSPEC Models

XSPEC models for relativistic wind emission and absorption.

## Models

### Emission geometry (development API)

| Model | Parameters before the supplied additive component |
| --- | --- |
| `windem` | `vout`, `Cf`, `hole`, `incl` |
| `windemlos` | `vout`, `Cf`, `skew`, `incl` |
| `windemoff` | `vout`, `Cf`, `hole`, `view`, `region` |

`vout` is in units of c; `incl` is in degrees. `Cf`, `hole`, `skew`, and
`view` lie in [0, 1]. `region` is a frozen switch, 0 for the polar sector and
1 for the equatorial sector. Fit the sectors separately.

**This is a breaking development change.** Existing XCM sessions must be
converted before loading: the parameter counts can be unchanged while the
parameter meanings and bounds differ. No legacy model aliases are provided.

`Cf = cos(tin) - cos(tout)` is the occupied fraction of one hemisphere. The
existing full-sphere emission scaling remains **Cf/2**, including the
zero-velocity limit. `Cf=0` gives zero emission; it does not measure angles.
`Cf=1` fills the hemisphere and leaves its placement coordinate unconstrained.

For generic geometry, `incl` is independent and `hole` divides the unused
solid angle between the polar cavity and the equatorial region:

```text
cos(tin)  = 1 - hole*(1-Cf)
cos(tout) = (1-hole)*(1-Cf)
```

The LoS model retains `incl` and `Cf` as native fitted quantities. `skew`
selects the shell containing that sightline. For `c=cos(incl)` and `x=cos(tout)`,
its physical quadratic root satisfies

```text
(1-skew)*(1-x-Cf)*(c-x) = skew*x*(x+Cf-c)
cos(tin) = x + Cf
```

`skew<0.5` favors polar fractional coverage, `skew>0.5` favors equatorial
fractional coverage. At 0.5 both sides have equal fractions of their available
solid angle covered by wind. This is not symmetry in degrees or in observed
flux. The implementation uses a cancellation-resistant root, including its
linear case at `skew=0.5`.

The off model first constructs the same shell as the generic model, then sets

```text
region=0: incl = view*tin
region=1: incl = tout + view*(90 deg-tout)
```

`view=0` always means the lower inclination of the sector and `view=1` the
upper. Endpoints can graze a boundary. If a sector collapses, its view
coordinate does not identify an independently measurable inclination.

The three `windline` counterparts use the same four/five geometry parameters,
followed by `LineE` in keV and XSPEC's additive `norm`. Each computes one line;
the `windem` family convolves a supplied additive spectrum with the same
physical profile, using conservative rebinning and a linear zero-padded FFT.

XSPEC's hard ranges constrain independent fit coordinates. For finite
numerical probes outside a range, the adapter evaluates the nearest physical
boundary and issues a one-time notice per parameter/spectrum. This supplies
a one-sided numerical extension; it is not permission to fit nonphysical
coordinates or impose out-of-range parameter links. Non-finite inputs remain
errors. Public Python constructors validate inputs strictly. Explicit inverted
angles in `Geometry.angles` issue a warning and represent zero emission.

The geometry implementation and native fit tests are available now. The
consolidated `windgeom` report and two-opening-angle `windsteppar` scan are
available, with optional generic FITS export; see the
[procedure guide](../procedures/README.md). `winderror` provides derived
opening-angle profiles and supported inclination errors; `windgeom` includes
local covariance uncertainties. Free-geometry off-wind inclination profiles
remain unsupported. Older `wineincl*` procedures
were written for the retired angular coordinates and must not be used as
conversion tools for this API. Native XSPEC `error` applies directly to `Cf`
and to generic/LoS `incl`; clear loaded chains when requesting a likelihood
profile rather than a posterior interval.

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
lmod wine /path/to/wine/xspec

# Velocity, hemispheric covering fraction, LoS skew, inclination in degrees.
model windemlos*atable{emission.fits}
0.3
0.4
0.5
45
/*

# Generic shell: same Cf, independent inclination, unused-sky split=0.5.
model windem*atable{emission.fits}
0.3
0.4
0.5
45
/*

# Wind-shift only the multiplicative transmission.
model (wind*mtable{absorption.fits})*powerlaw
```

Supply the remaining table parameters appropriate to your spectra. The
geometry itself does not add a cosmological-redshift parameter.

### Energy coverage for shifts and emission broadening

The convolution models (`wind`, `zwind`, `awind`, `azwind`, and the `windem`
family) receive a spectrum already evaluated by XSPEC. They remap those values;
they cannot request the underlying table or continuum at additional energies.
The **model calculation grid** must therefore cover the input energies that
can shift into the **fitted detector band**. Missing input coverage can create
artificial zero transmission or suppress emission near a grid boundary. FFT
zero padding prevents wrap-around, but does not supply the missing spectrum.

With no `energies` command, XSPEC uses each response's native energy grid. The
`.1--10 keV, 1000 linear bins` values in the XSPEC manual are defaults for an
*invoked* `energies` command with omitted fields, not the active grid otherwise.

For a shift `E_out = g E_in`, evaluating an output band `[E_low, E_high]`
requires input coverage `[E_low/g, E_high/g]`. If broadening or fitted parameter
ranges allow `g_min <= g <= g_max`, cover at least
`[E_low/g_max, E_high/g_min]` as well as the output band itself.

- `wind` and `awind` use `g = B = sqrt((1+beta)/(1-beta))`: positive outflow
  velocity requires a low-energy extension.
- `zwind` and `azwind` use `g = B/(1+z)`. Extend the low or high boundary as
  needed over the full allowed velocity/redshift range.
- The `windem` family broadens with `g = sqrt(1-beta^2)/(1-beta*cos(alpha))`,
  where `alpha` is the angle between the local velocity and the sightline.
  `1/B <= g <= B` is a conservative range; the selected geometry can narrow
  it. These models have no redshift parameter; apply any required cosmological
  redshift in the supplied additive component. Account for other convolutions
  in the expression when choosing the full calculation range.

For example, for a response spanning 1--10 keV and emission velocities up to
0.4c, 0.5--20 keV supplies margins beyond the conservative required input
range of about 0.65--15.3 keV, before any additional shifts. After loading the
data, extend both ends while retaining the response's native interior bins:

```tcl
energies extend low 0.5 400 log
energies extend high 20.0 400 log
```

The two extensions coexist. Choose their bounds within the supplied tables'
energy coverage, accounting for table redshift parameters. A wider XSPEC grid
cannot recover information beyond a table's tabulated support. Use enough
extension bins to resolve the input spectral features; coarse padding can
smear features that shift into the fitted band. Do not infer physical zero or
unit transmission outside the tabulated range.

Check the calculation by extending the bounds further, with adequate sampling,
and verifying that the model and fit statistic in the fitted band remain
stable to the precision required by the analysis. Also check sampling
convergence: wider coverage alone does not guarantee sufficient resolution.
The one-time WINE notices remind users of this requirement; they do not test
whether the current grid is adequate.

For a model-only plot, calculate on the wider grid and restrict the displayed
range separately. `energies 1 10 1000 lin` would replace the calculation grid
with 1--10 keV and remove its exterior coverage. `energies reset` returns to
the native response grids and removes the extensions.

The shifts follow the same supplied-array approach as XSPEC's `vashift` and
`vmshift`; WINE uses its documented relativistic outflow convention.
`wind` and `zwind` interpolate multiplicative factors directly on uneven grids,
whereas `awind` and `azwind` conservatively rebin additive flux.

See the [XSPEC energies reference](https://heasarc.gsfc.nasa.gov/docs/software/xspec/manual/XSenergies.html).

## Requirements

- XSPEC 12+
- FFTW3 library (for windem FFT convolution)
