# windkin

XSPEC local models for emission from radial winds with an intrinsic, emissivity-weighted distribution of speeds. `wuline` and `wuconv` use a uniform speed distribution (B). `wgline` and `wgconv` use a truncated Gaussian (C). The additive models return photons cm⁻² s⁻¹ per bin after XSPEC applies `norm`. The convolution models take and return bin-integrated photon flux.

## XSPEC command line

The native local-model library embeds Python to share the numerical kernel with
PyXspec. It requires NumPy, SciPy and the RED-WINE checkout. From the checkout
root with HEASoft initialized, rebuild with:

```sh
WINDKIN_PYTHON_CONFIG=/usr/local/bin/python3.11-config xspec/windkin/native/build.sh
xspec
```

Choose the Python config that matches the installed NumPy and SciPy. In
XSPEC, load the library from the path printed by the build. From the
checkout root the relative form also works:

```tcl
lmod windkin ./xspec/windkin/build/native
```

The library finds the Python sources relative to its own location inside
the checkout. Load it from `xspec/windkin/build/native` using an absolute
path when running XSPEC from another directory. `lmod windkin` needs no
environment exports. XSPEC's `syscall export` changes only the child shell's
environment. Load the main `redwine` package as well if the fit uses its
components.

The library registers `wuline`, `wuconv`, `wgline`, `wgconv`, `zwgline`, the three Gaussian LOS variants, and the three constant-speed sector models below.
Load the libraries before restoring XCM files that reference them.

## PyXspec

Use PyXspec in an initialized HEASoft environment. From the root of the RED-WINE checkout:

```sh
export PYTHONPATH="$PWD/xspec/windkin:$PWD/src${PYTHONPATH:+:$PYTHONPATH}"
python
```

```python
from xspec import AllModels, Model
from windkin.xspec import register

register()
line = Model("wgline")
line.setPars(.25, .03, 20, 50, 45, 0, 6.7, 1)
# Or convolve an intrinsic spectral model, for example:
broadened = Model("wuconv*gaussian")
```

Register before restoring an XCM that uses these components. Load the
compiled RED-WINE package in the same session with
`AllModels.lmod("redwine", ...)`. If you load the native windkin library via
`lmod`, its components are already registered. Both interfaces call the binned photon kernels in
`redwine.experimental.kernels`.

| Parameter | Meaning |
|---|---|
| `vmean` | Signed center in units of c. Positive is outflow, negative is inflow |
| `width` | B full width in units of c, from `abs(vmean) - width/2` to `abs(vmean) + width/2` in speed magnitude |
| `sigma` | C parent Gaussian standard deviation in units of c, with support truncated at ±3σ and at −1 < β < 1 |
| `tin`, `tout`, `incl` | Inner and outer cone half-opening angles and inclination, all in degrees |
| `bicone` | Switch, 0 for one cone and 1 to include the antipodal cone |
| `LineE`, `norm` | Additive only. Intrinsic line energy in keV and intrinsic all-sky-equivalent photon flux |

For C, the parent Gaussian is centered at the signed `vmean`. Negative
velocities are inflow and positive velocities are outflow. A Gaussian whose
support crosses zero contributes both directions, with the relative weights
given by the same normalized PDF. The truncated mean can differ from the
parent center. For B, the PDF stays uniform in speed magnitude and the sign
of `vmean` selects the direction for the whole interval.

For example, `newpar 1 -0.25 0.001 -0.99 -0.99 0.99 0.99` enables a negative
center for a standalone windkin component.

Zero width gives a single speed. For B require
`0 <= abs(vmean)-width/2 <= abs(vmean)+width/2 < 1`. Require `tin <= tout` and `tout <= 90` for a bicone. The Python API raises an error for invalid geometry or speed parameters. XSPEC emits a warning and returns zero for invalid combinations. Set fit bounds or parameter links so trial parameters stay inside the coupled constraints.

## Geometry and Doppler weights

The model assumes isotropic comoving emission from moving elements, a radial-velocity PDF independent of position within the emitting cone, and uniform intrinsic emissivity per unit wind solid angle. The azimuth fraction is evaluated analytically in the observer frame. One hollow cone subtends `(cos(tin)-cos(tout))/2` of the sky and a symmetric bicone twice that fraction. Inflow reverses the velocity directions. It maps a single cone at inclination `incl` to a velocity cone at `180°-incl`. A symmetric bicone and a transverse cone are invariant under this reversal.

Let μ be the cosine between the outward radial direction and the sight line and δ = √(1−β²)/(1−βμ). For a moving emitter, invariance of specific energy intensity divided by energy cubed gives
`I_E(E) = δ³ I′_E(E/δ)`. Dividing by E gives the photon density
`N_E(E) = δ² N′_E(E/δ)`. Then `dE = δ dE′` gives δ³ for integrated photons and δ⁴ for bolometric energy. Equivalently, for an intrinsic line, the Dirac factor transforms as
`δ_D(E/δ−E0) = δ δ_D(E−δ E0)`. For a steady optically thin Eulerian volume integral, the transfer invariant `j_E/E²` instead gives δ² for *integrated photon emissivity* at fixed lab-frame volume. Which convention applies depends on the emitter and the definition of its intrinsic emissivity. See [Urry & Padovani, Appendix B](https://ned.ipac.caltech.edu/level5/Urry1/UrryP_appenb.html) for the distinction between moving sources and continuous flows.

For `x = E/E0`, the Jacobian is `|dμ/dE| = √(1−β²)/(|β| E0 x²)`. Thus after eliminating μ the photon kernel per keV is

```text
K(E,E0) = x/(4π E0) ∫ p(β) w_phi(μ(x,β)) √(1−β²)/|β| dβ .
```

RED-WINE's existing mono-speed photon density has the same Jacobian with integrated weight δ². The additional factor x converts that density to the moving-emitter convention adopted here. The full energy integral is the wind-solid-angle average of δ³ over the speed PDF, with the solid-angle measure divided by 4π. It equals the geometric covering fraction at zero speed. For a full sphere at fixed β it equals γ. The speed PDF has unit integral. All geometric weights and Doppler factors remain in the output.

B reuses the semi-analytic β panels and the primitive on panels with constant azimuth fraction. C splits the signed Gaussian at zero and uses the same geometric panel
boundaries on each branch. The negative branch reflects the cone
inclination and uses the PDF at negative beta. Gauss–Legendre quadrature
evaluates each remaining one-dimensional integral. The Gaussian PDF is normalized over the intersection of the ±3σ interval and the physical domain −1 < β < 1. The floating-point representation limits tabulation at the last representable speed below c. Integrated bins stay finite at the unshifted-energy cusp when the PDF reaches β = 0.

## Finite energy bins

The kernel is tabulated in log Doppler factor using adaptive panels and integrated cubic interpolants. A convolution treats the input photon density as constant within each source energy bin. It uses two kernel moments to integrate the fraction of that bin shifted below each target edge, with direct CDF quadrature for very narrow source bins. Exact geometric energy grids admit a Toeplitz convolution evaluated via FFT. Irregular grids use the same bin-overlap integrals directly. Runtime for a dense spectrum on a very fine irregular grid scales with the number of overlapping source/target bin pairs.

Doppler shifts can move photons outside the model energy grid, and flux beyond the input grid cannot shift in. Use an extended grid where those contributions matter. With PyXspec, for example, `AllModels.setEnergies("0.1 40.0 20000 log")` sets a geometric calculation grid in keV. Check convergence against a finer and wider grid for the actual spectrum and response. The kernel describes bulk motion. Intrinsic thermal or turbulent line structure belongs in the input spectrum.

## Reproduce the checks

From the RED-WINE root with the same `PYTHONPATH`:

```sh
MPLBACKEND=Agg python xspec/windkin/check.py --plot /path/to/comparison.png
python xspec/windkin/check.py --xspec --redwine "$PWD"
```

For a real Resolve response run the XSPEC check from the spectrum directory and add `--spectrum /path/to/spectrum.pi`. Add `--dense` for continuum timing on the model grid. The checks compare independent three-dimensional quasi-Monte Carlo angular/velocity integrals, exact full-sphere results, limiting cases and finite-bin convolution. They also compare the Python outputs with real XSPEC model evaluations.

`zwgline` has the same parameters as `wgline` plus cosmological redshift
`z` immediately before the XSPEC norm. It evaluates the kernel at source-frame
bin edges and applies the same integrated photon-flux factor `1/(1+z)` as
RED-WINE’s `azwshift`.

## Gaussian line-of-sight variants

`wglinelos_Cf`, `zwglinelos_Cf`, and `wgconvlos_Cf` replace `tin` and `tout`
with `Cf` and `skew`. All other parameters and the Doppler kernel are shared
with the corresponding angular models. The RED-WINE LOS transform computes
a hollow cap with `tin <= incl <= tout`, where `incl` is the physical
inclination in degrees and can be linked across datasets. `Cf` and `skew`
range from zero to one. The solid-angle fraction of the full sky is
`Cf/2` for a single cap and `Cf` for an antipodal pair. `Cf=0` gives zero
emitting area. The LOS constraint concerns angular geometry.

## Constant-speed azimuthal sectors

`wsect` convolves bin-integrated photon spectra. `wsectline` emits a
standalone line, and `zwsectline` adds explicit cosmological redshift.
All three start with `vout, tin, tout, incl, phi, dphi, bicone`. The
additive variants then take `LineE`, optional `z` for `zwsectline`, and
XSPEC's additive `norm`. As for `zwgline`, redshift maps the energy edges
to the source frame and divides the integrated flux by `1+z`.

The velocity is signed in units of c and purely radial. All angles are in
degrees. Use `0 <= tin <= tout <= 180` for a single sector and
`tout <= 90` for a bicone. Inclination ranges from 0° to 180°. The sector
spans `phi ± dphi/2` modulo 360°, with `phi=0` toward the observer's
equatorial projection. At pole-on inclinations all sector orientations
are equivalent. Reflection `phi → -phi` gives the same profile, and an
antipodal bicone is invariant under reversal of the signed velocity.
`bicone=1` adds the sector at `(180°-theta, phi+180°)`.

For each Doppler factor δ, the signed Doppler relation fixes a ring at
observer polar cosine μ. The code intersects the polar-cap arcs with
the analytical roots of the two sector meridians on this ring. Sectors
wider than 180° use the complement of the smaller wedge. Meridional and
polar tangencies and their intersections partition the remaining
one-dimensional integration in log δ. The existing adaptive bin and
convolution machinery integrates the resulting kernel.

These models use the windkin moving-emitter photon weight δ³.
`dphi=360` is exactly the zero-width axisymmetric windkin profile.
At zero velocity the integrated fraction for one sector is
`(cos(tin)-cos(tout))*dphi/720` with the angles expressed in degrees.
An antipodal pair doubles that fraction. At finite velocity the flux is
the integral of δ³ over the emitting solid angle divided by 4π.
A zero-width sector or an equal pair of polar bounds emits zero flux.

The kernel has no intrinsic line width. In XSPEC, `wsect` conservatively
rebins dense non-logarithmic grids to an internal logarithmic grid with four
times as many bins, convolves by FFT, and rebins back. It retains the native
energy limits, so input outside those limits still requires explicit energy
extension. It uses direct bin overlaps for small grids and grids starting
at zero. RED-WINE's axisymmetric geometry procedures do not report sector
solid angles or sector parameter uncertainties.

In Python use `redwine.experimental.kernels.Sector` with the same seven
parameters, then pass it to `line(edges, rest_energy, sector)` or
`convolve(edges, input_flux, sector)`. The
[sector slider notebook](../../example/Slider-Emission-Profiles-Sector.ipynb)
compares the sector to the full-azimuth profile at fixed intrinsic line
energy and normalization.
