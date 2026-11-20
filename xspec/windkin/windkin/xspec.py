"""XSPEC local-model adapters. Call register() once per PyXspec session."""

import warnings

import numpy as np

from redwine.geometry import Geometry
from redwine.experimental.kernels import Sector

from .kernel import Wind, convolve, line

_COMMON = (
    'vmean c 0.25 0 0 0.99 0.99 0.001',
    'tin deg 20 0 0 90 180 0.01',
    'tout deg 50 0 0 90 180 0.01',
    'incl deg 45 0 0 90 180 0.01',
    '$bicone 0',
)
_SPREAD = {
    "uniform": 'width c 0.10 0 0 1 1 0.001',
    "gaussian": 'sigma c 0.03 0 0 1 10 0.001',
}
_LOS_GEOMETRY = (
    'Cf "" 0.3 0 0 1 1 0.01',
    'skew "" 0.5 0 0 1 1 0.01',
    'incl deg 45 0 0 90 90 0.01',
    '$bicone 0',
)
_LINE_ENERGY = 'LineE keV 6.7 1e-6 0.01 100 1e6 0.001'
_REDSHIFT = 'z "" 0 0 0 10 10 0.001'
_SECTOR = (
    'vout c 0.25 -0.99 -0.99 0.99 0.99 0.001',
    *_COMMON[1:4],
    'phi deg 0 0 0 360 360 0.1',
    'dphi deg 90 0 0 360 360 0.1',
    '$bicone 0',
)
_registered = False


def _evaluate(law, additive, edges, params, flux, *, los=False):
    vmean, spread, inner, outer, incl, bicone = params[:6]
    try:
        if los:
            geometry = Geometry.los(Cf=inner, skew=outer, incl=incl)
            inner, outer = geometry.tin, geometry.tout
        center = vmean if law == "gaussian" else abs(vmean)
        inflow = law == "uniform" and vmean < 0
        wind = Wind(law, center, spread, inner, outer, incl, inflow, bool(bicone))
        values = line(edges, params[6], wind) if additive else convolve(edges, flux, wind)
    except ValueError as exc:
        warnings.warn(f"windkin: {exc}", RuntimeWarning, stacklevel=3)
        flux[:] = [0.] * (len(edges)-1)
    else:
        flux[:] = values.tolist()


def wuline(edges, params, flux):
    _evaluate("uniform", True, edges, params, flux)


def wuconv(edges, params, flux):
    _evaluate("uniform", False, edges, params, flux)


def wuconvlos_Cf(edges, params, flux):
    _evaluate("uniform", False, edges, params, flux, los=True)


def wgline(edges, params, flux):
    _evaluate("gaussian", True, edges, params, flux)


def _redshifted_line(edges, params, flux, *, los=False):
    redshift_factor = 1 + params[7]
    _evaluate("gaussian", True, np.asarray(edges)*redshift_factor, params, flux, los=los)
    # Same additive normalization convention as RED-WINE's azwshift.
    flux[:] = [value/redshift_factor for value in flux]


def zwgline(edges, params, flux):
    _redshifted_line(edges, params, flux)


def wglinelos_Cf(edges, params, flux):
    _evaluate("gaussian", True, edges, params, flux, los=True)


def zwglinelos_Cf(edges, params, flux):
    _redshifted_line(edges, params, flux, los=True)


def wgconvlos_Cf(edges, params, flux):
    _evaluate("gaussian", False, edges, params, flux, los=True)


def wgconv(edges, params, flux):
    _evaluate("gaussian", False, edges, params, flux)


def _evaluate_sector(edges, params, flux, *, additive, redshift_factor=1.):
    try:
        sector = Sector(*params[:7])
        if additive:
            values = line(np.asarray(edges)*redshift_factor, params[7], sector)/redshift_factor
        else:
            values = convolve(edges, flux, sector, resample=True)
    except ValueError as exc:
        warnings.warn(f"windkin: {exc}", RuntimeWarning, stacklevel=3)
        flux[:] = [0.] * (len(edges)-1)
    else:
        flux[:] = values.tolist()


def wsect(edges, params, flux):
    _evaluate_sector(edges, params, flux, additive=False)


def wsectline(edges, params, flux):
    _evaluate_sector(edges, params, flux, additive=True)


def zwsectline(edges, params, flux):
    _evaluate_sector(edges, params, flux, additive=True, redshift_factor=1+params[8])


def register():
    """Register both additive lines and convolution models in XSPEC."""
    global _registered
    if _registered:
        return
    from xspec import AllModels

    for law, additive, convolution in (("uniform", wuline, wuconv),
                                       ("gaussian", wgline, wgconv)):
        params = (_COMMON[0], _SPREAD[law], *_COMMON[1:])
        AllModels.addPyMod(additive, (*params, _LINE_ENERGY), "add")
        AllModels.addPyMod(convolution, params, "con")
    params = (_COMMON[0], _SPREAD["gaussian"], *_COMMON[1:])
    AllModels.addPyMod(zwgline, (*params, _LINE_ENERGY, _REDSHIFT), "add")
    AllModels.addPyMod(wuconvlos_Cf, (_COMMON[0], _SPREAD["uniform"], *_LOS_GEOMETRY), "con")
    params = (_COMMON[0], _SPREAD["gaussian"], *_LOS_GEOMETRY)
    AllModels.addPyMod(wglinelos_Cf, (*params, _LINE_ENERGY), "add")
    AllModels.addPyMod(zwglinelos_Cf, (*params, _LINE_ENERGY, _REDSHIFT), "add")
    AllModels.addPyMod(wgconvlos_Cf, params, "con")
    # XSPEC resolves prefixes in registration order.
    AllModels.addPyMod(wsect, _SECTOR, "con")
    AllModels.addPyMod(wsectline, (*_SECTOR, _LINE_ENERGY), "add")
    AllModels.addPyMod(zwsectline, (*_SECTOR, _LINE_ENERGY, _REDSHIFT), "add")
    _registered = True
