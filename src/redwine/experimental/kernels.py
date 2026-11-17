"""Binned photon kernels for radial velocity distributions and emitting sectors.

Energies are in keV and fluxes are integrated per bin. The intrinsic velocity
PDF has unit integral. The observed kernel retains the wind solid angle and
the moving-emitter Doppler weight (delta cubed for integrated photons).
"""

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from scipy.interpolate import PPoly
from scipy.signal import fftconvolve
from scipy.special import erf, roots_legendre

from ._uniform import Prepared, Profile
from ._sector import Sector, SectorRings

__all__ = ["Wind", "Sector", "line", "convolve"]


@dataclass(frozen=True)
class Wind:
    """Uniform (B) or truncated-normal (C) speed PDF and hollow-cone geometry.

    Center and spread are in units of c. Uniform center is a speed magnitude
    and spread is the full width. Gaussian center is a signed radial velocity
    and spread is the parent sigma. The Gaussian support is limited to three
    sigma and to subluminal velocities. Inflow reverses all velocities. Angles
    are in degrees.
    """

    law: str
    center: float
    spread: float
    inner: float
    outer: float
    inclination: float
    inflow: bool = False
    bicone: bool = False

    def __post_init__(self):
        if (self.law not in ("uniform", "gaussian")
                or not np.isfinite([self.center, self.spread, self.inner,
                                    self.outer, self.inclination]).all()
                or not abs(self.center) < 1 or self.spread < 0
                or not 0 <= self.inner <= self.outer <= (90 if self.bicone else 180)
                or not 0 <= self.inclination <= 180):
            raise ValueError("Invalid speed distribution or hollow-cone geometry.")
        if self.law == "uniform":
            if not 0 <= self.center - self.spread/2 <= self.center + self.spread/2 < 1:
                raise ValueError("Uniform speed interval must lie within [0, 1).")

    @property
    def solid_angle_fraction(self):
        inner, outer = np.deg2rad([self.inner, self.outer])
        return np.sin((outer + inner)/2)*np.sin((outer - inner)/2)*(1 + self.bicone)

    def _gaussian_limits(self):
        """Bounds of the truncated Gaussian in standard deviations."""
        return (max(-3., (-1 - self.center)/self.spread),
                min(3., (1 - self.center)/self.spread))

    def pdf(self, beta):
        """Intrinsic density per unit beta, away from the zero-width Dirac limit."""
        if self.spread == 0:
            raise ValueError("Zero width is a Dirac measure.")
        beta = np.asarray(beta, dtype=float)
        if self.law == "uniform":
            inside = abs(beta - self.center) <= self.spread/2
            return np.where(inside, 1/self.spread, 0.)
        z = (beta - self.center)/self.spread
        lower, upper = self._gaussian_limits()
        # The CDF difference straddles zero. erf avoids cancellation for broad PDFs.
        area = (erf(upper/np.sqrt(2)) - erf(lower/np.sqrt(2)))/2
        return np.where((abs(beta) < 1) & (z >= lower) & (z <= upper),
                        np.exp(-z*z/2)/(self.spread*area*np.sqrt(2*np.pi)), 0.)


def energy_edges(edges):
    edges = np.asarray(edges, dtype=float)
    if (edges.ndim != 1 or edges.size < 2 or not np.isfinite(edges).all()
            or edges[0] < 0 or np.any(np.diff(edges) <= 0)):
        raise ValueError("Energy edges must be finite, nonnegative and strictly increasing.")
    return edges


def _log_delta(beta, mu):
    return (np.log1p(-beta) + np.log1p(beta))/2 - np.log1p(-beta*mu)


# On each log-Doppler panel fit cubics at four interior Gauss nodes, and
# check both moments at eight independent nodes before integrating the cubics.
_fit_nodes, _ = roots_legendre(4)
_check_nodes, _check_weights = roots_legendre(8)
_fit_nodes = (_fit_nodes + 1)/2
_check_nodes = (_check_nodes + 1)/2
_check_weights /= 2
_fit = np.linalg.inv(np.vander(_fit_nodes, 4))
_check = np.vander(_check_nodes, 4)
_beta_nodes, _beta_weights = roots_legendre(24)
_beta_position = np.sin(np.pi*(_beta_nodes + 1)/4)**2
_beta_weights *= np.pi/4*np.sin(np.pi*(_beta_nodes + 1)/2)


@lru_cache(maxsize=8)
def _moments(wind, lower, upper, tolerance):
    """Antiderivatives of the photon kernel and of the kernel divided by delta."""
    if isinstance(wind, Sector):
        return _sector_moments(wind, lower, upper, tolerance)
    if wind.inner == wind.outer or lower >= upper:
        return None
    incl, inner, outer = np.deg2rad([wind.inclination, wind.inner, wind.outer])
    if wind.inflow:
        incl = np.pi - incl
    speed = abs(wind.center)
    mono = wind.spread <= 4*np.spacing(speed)
    uniform = wind.law == "uniform" and not mono
    if mono:
        branches = [(np.array([speed]), -1 if wind.center < 0 else 1)]
    elif uniform:
        branches = [(wind.center + wind.spread*np.array([-.5, .5]), 1)]
    else:
        # Split the signed Gaussian at zero and map each branch to b=|beta|.
        z_lower, z_upper = wind._gaussian_limits()
        signed_knots = wind.center + wind.spread*np.array([z_lower, 0., z_upper])
        branches = [(np.unique(np.clip(sign*signed_knots, 0, 1)), sign)
                    for sign in (1, -1)]
        branches = [(knots, sign) for knots, sign in branches if len(knots) > 1]

    # Prepared supplies azimuth geometry and velocity panel boundaries.
    caps = []
    knots = [lower, upper]
    for beta_knots, sign in branches:
        beta_min, beta_max = beta_knots[[0, -1]]
        angle = incl if sign > 0 else np.pi - incl
        inclinations = (angle, np.pi - angle) if wind.bicone else (angle,)
        if mono or uniform:
            center, width = speed, (wind.spread if uniform else 0.)
        else:
            limit = min(beta_max, np.nextafter(1., 0.))
            center = beta_min + (limit-beta_min)/2
            width = min(2*(center-beta_min), 2*(limit-center))
        for angle in inclinations:
            cap = Prepared(Profile(1., center, width, angle, inner, outer))
            caps.append((cap, beta_knots, sign))
            for mu in cap.tangencies:
                beta = np.r_[beta_knots, mu] if beta_min < mu < beta_max else beta_knots
                knots.extend(_log_delta(beta[beta < 1], mu))
    # All-sky envelope. Include zero for the integrable beta=0 cusp.
    max_speed = max(beta_knots[-1] for beta_knots, _ in branches)
    envelope = np.arctanh(min(max_speed, np.nextafter(1., 0.)))
    lo, hi = max(lower, -envelope), min(upper, envelope)
    if lo >= hi:
        return None
    def log_density(log_shift):
        shift = np.exp(log_shift)
        offset = -np.expm1(-log_shift)
        density = np.zeros_like(shift)
        for cap, beta_knots, sign in caps:
            if uniform:
                density += cap.evaluate(shift, tolerance=2e-9)
            elif mono:
                density += cap.density(speed, offset)/2
            else:
                row, left, right = cap.panels(offset)
                # The support is already partitioned by the angular roots.
                # Split long intervals at the Gaussian knots as well.
                for a, b in zip(beta_knots[:-1], beta_knots[1:]):
                    start, stop = np.maximum(left, a), np.minimum(right, b)
                    use = (stop > start) & (cap.arcs(cap.mu((start + stop)/2, offset[row])) > 0)
                    # A log beta coordinate resolves the 1/beta factor
                    # close to the unshifted energy.
                    idx, start, stop = row[use], start[use], stop[use]
                    logarithmic = (start > 0) & (start < stop/4)
                    span = stop - start
                    beta = start[:, None] + span[:, None]*_beta_position
                    jac = np.broadcast_to(span[:, None], beta.shape).copy()
                    log_span = np.log(stop[logarithmic]/start[logarithmic])
                    beta[logarithmic] = start[logarithmic, None]*np.exp(
                        log_span[:, None]*_beta_position)
                    jac[logarithmic] = beta[logarithmic]*log_span[:, None]
                    value = cap.density(beta, offset[idx, None])*wind.pdf(sign*beta)*jac
                    np.add.at(density, idx, (value @ _beta_weights)/2)
        # RED-WINE's density per x contains the delta^2 photon weight.
        # Multiply by x for delta^3, and again for dx/dlog(x).
        return np.stack((density*shift**2, density*shift), axis=-1)

    return _tabulate_moments(log_density, knots, lo, hi, wind.solid_angle_fraction, tolerance)


def _sector_moments(sector, lower, upper, tolerance):
    if sector.solid_angle_fraction == 0 or lower >= upper:
        return None
    if sector.dphi == 360:
        return _moments(Wind("gaussian", sector.vout, 0., sector.tin, sector.tout,
                             sector.incl, bicone=sector.bicone), lower, upper, tolerance)

    rings = SectorRings(sector)
    beta = sector.vout
    speed_factor = np.sqrt((1-beta)*(1+beta))
    knots = _log_delta(beta, rings.knots)
    envelope = np.arctanh(abs(beta))
    lo, hi = max(lower, -envelope), min(upper, envelope)
    if lo >= hi:
        return None

    def log_density(log_shift):
        # Invert the Doppler relation with the signed beta. On each ring,
        # the analytic meridian roots select the emitting azimuth intervals.
        offset = -np.expm1(-log_shift)
        mu = offset/beta + (1-offset)*beta/(1+speed_factor)
        density = rings.arcs(mu)*speed_factor/(2*abs(beta))
        shift = np.exp(log_shift)
        return np.stack((density*shift**2, density*shift), axis=-1)

    return _tabulate_moments(log_density, knots, lo, hi,
                             sector.solid_angle_fraction, tolerance)


def _tabulate_moments(log_density, knots, lo, hi, angular_fraction, tolerance):
    knots = np.unique(np.clip(np.r_[knots, lo, hi, 0.], lo, hi))
    # Remove sub-ulp gaps between equivalent tangencies.
    knots = knots[np.r_[True, np.diff(knots) > 8*np.finfo(float).eps *
                        np.maximum(abs(knots[:-1]), min(1e-12, hi-lo))]]
    # Limit initial panel size even on extended grids and broad distributions.
    knots = np.unique(np.r_[knots, np.linspace(lo, hi, 1 + int(np.ceil((hi-lo)/.12)))])
    left, right = knots[:-1], knots[1:]
    accepted_left, accepted_coefficients = [], []
    for _ in range(48):
        if not len(left):
            break
        width = right - left
        sample = left[:, None] + width[:, None]*np.r_[_fit_nodes, _check_nodes]
        values = log_density(sample.ravel()).reshape(-1, 12, 2)
        coeff = _fit @ values[:, :4]
        actual = values[:, 4:]
        error = _check_weights @ abs(_check @ coeff - actual)
        scale = _check_weights @ abs(actual)
        accept = np.all(width[:, None]*error <= tolerance*(
            width[:, None]*scale + angular_fraction*(np.sqrt(width[:, None]/(hi-lo)) + 1e-4)), axis=1)
        # Keep the exact panel width in converting the coefficients to PPoly.
        accepted_left.append(left[accept])
        accepted_coefficients.append(coeff[accept]/width[accept, None, None]**np.arange(3, -1, -1)[None, :, None])
        mid = left[~accept] + width[~accept]/2
        left, right = np.r_[left[~accept], mid], np.r_[mid, right[~accept]]
    if len(left):
        raise RuntimeError(f"Doppler kernel tabulation did not converge on {len(left)} panels.")
    left = np.concatenate(accepted_left)
    order = np.argsort(left)
    coeff = np.concatenate(accepted_coefficients)[order].transpose(1, 0, 2)
    return PPoly(coeff, np.r_[left[order], hi], extrapolate=False).antiderivative()


def _static(wind):
    if isinstance(wind, Sector):
        return wind.vout == 0
    return wind.center == 0 and wind.spread == 0


def line(edges, rest_energy, wind, *, tolerance=1e-7):
    """Photon flux per observer bin for a Wind or Sector with unit intrinsic flux."""
    edges = energy_edges(edges)
    if not np.isfinite(rest_energy) or rest_energy <= 0:
        raise ValueError("Rest energy must be positive and finite.")
    if _static(wind):
        return wind.solid_angle_fraction*((edges[:-1] <= rest_energy) &
                                          (rest_energy < edges[1:]))
    with np.errstate(divide="ignore"):
        log_ratio = np.log(edges/rest_energy)
    moments = _moments(wind, max(-360., log_ratio[0]), min(360., log_ratio[-1]), tolerance)
    if moments is None:
        return np.zeros(len(edges)-1)
    cumulative = moments(np.clip(log_ratio, *moments.x[[0, -1]]))[:, 0]
    return np.maximum(0., np.diff(cumulative))


def _partial_cdf(y, a, b, moments):
    """Cumulative output below y from uniform photon density in [a,b]."""
    lo, hi = moments.x[[0, -1]]
    with np.errstate(divide="ignore"):
        t_b = np.clip(np.log(y/b), lo, hi)
        t_a = np.clip(np.log(y/a), lo, hi)
    if a > 0 and (b-a)/a < 1e-5:
        # Near coincident source edges, subtracting antiderivatives and
        # dividing by (b-a) loses precision. Integrate the CDF instead.
        with np.errstate(divide="ignore"):
            log_y_over_b = np.log(y/b)
        cumulative = moments(np.clip(log_y_over_b[..., None]
                                     - np.log1p(-((b-a)/b)*_check_nodes), lo, hi))[..., 0] @ _check_weights
        first = np.searchsorted(moments.x, t_b, side="right")
        last = np.searchsorted(moments.x, t_a, side="left")
        for index in np.flatnonzero(last > first):
            # Split where the CDF polynomial changes. Integrate in
            # fractional source energy u=(b-E)/(b-a) with a constant PDF.
            ratio = (b-a)/b
            t0 = log_y_over_b[index]
            breaks = -np.expm1(t0 - moments.x[first[index]:last[index]])/ratio
            cuts = np.r_[0., np.clip(breaks, 0., 1.), 1.]
            widths = np.diff(cuts)
            t = t0 - np.log1p(-ratio*(cuts[:-1, None] + widths[:, None]*_check_nodes))
            cumulative[index] = widths @ (moments(np.clip(t, lo, hi))[..., 0] @ _check_weights)
        return cumulative
    lower = moments(t_b)
    upper = moments(t_a)
    partial = (y*(upper[..., 1] - lower[..., 1])
               - a*(upper[..., 0] - lower[..., 0]))/(b-a)
    return np.clip(lower[..., 0] + partial, 0., moments(hi)[0])


def convolve(edges, flux, wind, *, tolerance=1e-7, resample=False):
    """Apply the wind kernel to integrated photons, assuming constant dN/dE per bin.

    With resample=True, use a fourfold-oversampled log-energy grid for
    dense non-geometric grids. Conservative rebinning preserves photon
    integrals on the common domain, but approximates within-bin structure.
    """
    edges = energy_edges(edges)
    flux = np.asarray(flux, dtype=float)
    if flux.shape != (len(edges)-1,) or not np.isfinite(flux).all():
        raise ValueError("Supply one finite integrated flux per energy bin.")
    if _static(wind):
        return wind.solid_angle_fraction*flux
    result = np.zeros_like(flux)
    source = np.flatnonzero(flux)
    if not len(source):
        return result
    with np.errstate(divide="ignore"):
        log_edges = np.log(edges)
    n = len(flux)
    geometric = edges[0] > 0 and np.allclose(
        np.diff(log_edges), log_edges[1] - log_edges[0], rtol=1e-10, atol=1e-14
    )
    # Small grids are cheap to integrate directly. Logarithmic rebinning
    # requires a strictly positive lower edge.
    if resample and not geometric and n >= 1024 and edges[0] > 0:
        uniform_edges = np.geomspace(edges[0], edges[-1], 4*n + 1)
        input_cumulative = np.r_[0., np.cumsum(flux)]
        uniform_flux = np.diff(np.interp(uniform_edges, edges, input_cumulative))
        broadened = convolve(uniform_edges, uniform_flux, wind, tolerance=tolerance)
        output_cumulative = np.r_[0., np.cumsum(broadened)]
        return np.diff(np.interp(edges, uniform_edges, output_cumulative))

    moments = _moments(wind, max(-360., log_edges[0] - log_edges[source[-1]+1]),
                       min(360., log_edges[-1] - log_edges[source[0]]), tolerance)
    if moments is None:
        return result

    if geometric:
        # A geometric grid gives a translation-invariant bin-to-bin operator.
        log_step = (log_edges[-1] - log_edges[0])/n
        offsets = np.arange(-n+1, n+1)*log_step
        with np.errstate(over="ignore", under="ignore"):
            weights = np.diff(_partial_cdf(np.exp(np.clip(offsets, -360, 360)),
                                           1., np.exp(log_step), moments))
        result = fftconvolve(flux, np.maximum(0., weights), mode="same")
        return np.maximum(0., result) if np.all(flux >= 0) else result

    # Only pairs whose Doppler support intersects the target bin contribute.
    # Group by source bin to keep memory linear even on a wide Resolve grid.
    lo, hi = np.exp(moments.x[[0, -1]])
    for index in source:
        a, b = edges[index:index+2]
        first = max(0, np.searchsorted(edges, a*lo, side="right")-1)
        last = min(n, np.searchsorted(edges, b*hi))
        if last > first:
            y = edges[first:last+1]
            cumulative = np.zeros_like(y)
            cumulative[y > 0] = _partial_cdf(y[y > 0], a, b, moments)
            result[first:last] += flux[index]*np.maximum(0., np.diff(cumulative))
    return result
