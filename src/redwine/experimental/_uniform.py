"""Experimental RED-WINE line densities per keV, for uniform intrinsic d beta.

Angles are in radians. Only outward speeds, 0 <= beta < 1, are supported.
An interval reaching beta=0 can give an integrable infinite density at E0.
"""
from dataclasses import dataclass
from math import atan2, pi, sqrt

import numpy as np
from scipy.integrate import quad
from scipy.special import roots_legendre

from redwine.kernels import lineshape


@dataclass(frozen=True)
class Profile:
    E0: float
    beta0: float
    width: float
    inclination: float
    inner: float
    outer: float
    photon: bool = True
    norm: float = 1.0

    def __post_init__(self):
        if (not np.isfinite([self.E0, self.beta0, self.width, self.inclination,
                            self.inner, self.outer, self.norm]).all()
                or self.E0 <= 0 or self.norm < 0 or self.width < 0
                or not 0 <= self.lower <= self.upper < 1
                or not all(0 <= angle <= np.pi
                           for angle in (self.inclination, self.inner, self.outer))):
            raise ValueError("Require E0 > 0, norm >= 0, angles in [0, pi], and 0 <= beta < 1.")
        if self.width and self.upper == self.lower:
            raise ValueError("Width is below floating-point resolution.")

    @property
    def lower(self):
        return self.beta0 - self.width / 2

    @property
    def upper(self):
        return self.beta0 + self.width / 2


RULES = [roots_legendre(order) for order in (8, 16)]


class Prepared:
    def __init__(self, profile):
        self.p = profile
        self.cos_i = 0. if profile.inclination == np.pi/2 else np.cos(profile.inclination)
        self.sin_i = np.sin(profile.inclination)
        self.cos_in, self.cos_out = [0. if angle == np.pi/2 else np.cos(angle)
                                     for angle in (profile.inner, profile.outer)]
        i, a, b = profile.inclination, profile.inner, profile.outer
        limb = abs(np.pi/2-i)
        angle_tol = 4*np.finfo(float).eps
        if self.sin_i > angle_tol:
            if abs(a-limb) <= angle_tol or abs(a-(np.pi-limb)) <= angle_tol:
                self.cos_in = np.copysign(self.sin_i, self.cos_in)
            if abs(b-limb) <= angle_tol or abs(b-(np.pi-limb)) <= angle_tol:
                self.cos_out = np.copysign(self.sin_i, self.cos_out)
        tangencies = np.cos([i-a, i+a, i-b, i+b, 0, np.pi])
        self.tangencies = np.unique(np.where(np.abs(tangencies) <= angle_tol, 0., tangencies))
        self.empty = a >= b or profile.norm == 0

    @staticmethod
    def mu(beta, offset):
        # offset = 1 - E0/E. Rationalize 1 - sqrt(1-beta^2).
        speed_factor = np.sqrt((1-beta)*(1+beta))
        return offset/beta + (1-offset)*beta/(1+speed_factor)

    def polar_limits(self, mu):
        """Positive observer-ring azimuth interval between the two cap roots."""
        center = self.cos_i * mu
        sine = np.sqrt(np.maximum(0, (1-mu)*(1+mu)))
        partial = (sine*self.sin_i > 0)
        # Evaluate acos via the distances to both clipping boundaries.
        # Near mu=0, preserve the small distance at an equatorial tangency.
        correction = self.sin_i*mu**2/(1+sine)

        def cap(cos_angle):
            near = np.abs(mu) < .25
            radius = self.sin_i*sine
            numerator = cos_angle-center
            minus = np.where(near, (self.sin_i-cos_angle)+center-correction, radius-numerator)
            plus = np.where(near, (self.sin_i+cos_angle)-center-correction, radius+numerator)
            return np.arctan2(np.sqrt(np.maximum(0, minus)), np.sqrt(np.maximum(0, plus)))

        aligned = (center >= self.cos_out) & (center <= self.cos_in)
        lower = np.where(partial, 2*cap(self.cos_in), 0.)
        upper = np.where(partial, 2*cap(self.cos_out), np.pi*aligned)
        return lower, upper

    def arcs(self, mu):
        lower, upper = self.polar_limits(mu)
        return np.where(np.abs(mu) <= 1, np.maximum(0., upper-lower)/np.pi, 0.)

    def panels(self, offset):
        # u = tanh(atanh(beta)/2), r = E0/E. The unsquared Doppler
        # relation is (1+r)u² - 2 mu u + 1-r = 0.
        mu = self.tangencies[None, :]
        q = offset[:, None]
        discriminant = mu**2 - q*(2-q)
        root = np.sqrt(np.maximum(0, discriminant))
        large = mu + np.copysign(root, mu)
        small = np.divide(q, large, out=np.zeros_like(large), where=large != 0)
        with np.errstate(divide="ignore", invalid="ignore"):
            u = np.concatenate((large/(2-q), small), axis=1)
        valid = np.tile(discriminant >= 0, (1, 2)) & (u > 0) & (u < 1)
        beta = 2*np.where(valid, u, 0)/(1+np.where(valid, u, 0)**2)
        p = self.p
        edges = np.sort(np.column_stack((np.full(len(q), p.lower),
                                        np.clip(beta, p.lower, p.upper),
                                        np.full(len(q), p.upper))), axis=1)
        row, column = np.nonzero(np.diff(edges, axis=1) > 0)
        return row, edges[row, column], edges[row, column+1]

    def density(self, beta, offset):
        return self.arcs(self.mu(beta, offset))*np.sqrt((1-beta)*(1+beta))/beta

    def evaluate(self, energies, method="semi", tolerance=2e-9):
        """Evaluate at positive energy centers. Direct and semi adaptively bisect Gauss panels.

        The tolerance controls successive panel estimates in units of norm/E0
        before the energy-flux factor. Reference uses scalar adaptive QUADPACK.
        """
        if method not in ("reference", "direct", "semi"):
            raise ValueError("Unknown integration method.")
        energy = np.asarray(energies, dtype=float)
        if not np.all(np.isfinite(energy) & (energy > 0)):
            raise ValueError("Energies must be finite and positive.")
        p = self.p
        if p.width == 0:
            if p.beta0 == 0:
                # The monoenergetic limit is a delta measure.
                area = p.norm*(self.cos_in-self.cos_out)/2
                return np.where((energy == p.E0) & (area > 0), np.inf, 0.)
            return lineshape(energy, p.E0, p.norm, p.beta0, p.inclination,
                             p.inner, p.outer, is_photonflux=p.photon)
        result = np.zeros(energy.size)
        if self.empty:
            return result.reshape(energy.shape)
        energy_flat = energy.ravel()
        # The all-sky kinematic envelope also bounds every conical profile.
        gmax = np.sqrt((1+p.upper)/(1-p.upper))
        active = np.flatnonzero((energy_flat > p.E0/gmax) & (energy_flat < p.E0*gmax))
        e = energy_flat[active]
        offset = (e-p.E0)/e
        row, a, b = self.panels(offset)
        q = offset[row]
        span = b-a
        midpoint_mu = self.mu(a+span/2, q)
        w = self.arcs(midpoint_mu)
        # Midpoint identifies the geometric regime because no tangency lies
        # inside a panel. An aligned cone has only full and empty rings.
        singular = (a == 0) & (q == 0)
        if np.any(singular):
            # At E0, mu tends to zero from the positive side.
            aligned = self.sin_i <= np.finfo(float).eps
            if aligned:
                divergent = ((self.cos_out <= 0 < self.cos_in) if self.cos_i > 0
                             else (self.cos_out < 0 <= self.cos_in))
            else:
                # At mu=0 the range of wind-axis polar angles is [limb, pi-limb].
                limb = abs(np.pi/2-p.inclination)
                angle_tol = 4*np.finfo(float).eps
                divergent = p.outer-limb > angle_tol and np.pi-limb-p.inner > angle_tol
            if divergent:
                result[active[row[singular]]] = np.inf
        keep = np.isfinite(result[active[row]]) & ((w > 0) | (method == "reference"))
        row, a, b, q, w, span = (x[keep] for x in (row, a, b, q, w, span))
        integral = np.zeros(len(row))
        numerical = np.ones(len(row), dtype=bool)
        if method == "semi":
            # Constant fraction, including the aligned case. Evaluate the
            # difference of F(beta)=sqrt(1-beta²)+log(beta/(1+sqrt(1-beta²))).
            constant = (w == 1) | (self.sin_i <= np.finfo(float).eps)
            if self.cos_i == 0 and ((p.inner == 0 and p.outer == np.pi/2)
                                    or (p.inner == np.pi/2 and p.outer == np.pi)):
                constant[:] = True  # Equatorial hemisphere intersects every ring in half.
            numerical = ~constant
            lo, hi = a[constant], b[constant]
            s_lo = np.sqrt((1-lo)*(1+lo))
            s_hi = np.sqrt((1-hi)*(1+hi))
            ds = -(hi-lo)*(hi+lo)/(s_hi+s_lo)
            primitive_difference = np.log1p((hi-lo)/lo) - np.log1p(ds/(1+s_lo)) + ds
            fast = lo >= .8
            if np.any(fast):
                # F = s-atanh(s). Use atanh subtraction before evaluating
                # atanh(x)-x to avoid losing the small integral near beta=1.
                product = s_lo[fast]*s_hi[fast]
                x = -ds[fast]/(1-product)
                remainder = np.where(x < .001, x**3*(1/3+x*x*(1/5+x*x/7)),
                                     np.arctanh(x)-x)
                primitive_difference[fast] = remainder + x*product
            integral[constant] = w[constant]*primitive_difference/p.width
        if method == "reference":
            def scalar_density(beta, offset):
                speed_factor = sqrt((1-beta)*(1+beta))
                mu = offset/beta + (1-offset)*beta/(1+speed_factor)
                if abs(mu) > 1:
                    return 0.
                center = self.cos_i*mu
                radius = self.sin_i*sqrt((1-mu)*(1+mu))
                if radius == 0:
                    fraction = float(self.cos_out <= center <= self.cos_in)
                else:
                    def cap(cos_angle):
                        if abs(mu) < .25:
                            correction = self.sin_i*mu**2/(1+sqrt((1-mu)*(1+mu)))
                            minus = (self.sin_i-cos_angle)+center-correction
                            plus = (self.sin_i+cos_angle)-center-correction
                        else:
                            numerator = cos_angle-center
                            minus, plus = radius-numerator, radius+numerator
                        return atan2(sqrt(max(0., minus)), sqrt(max(0., plus)))

                    fraction = max(0., cap(self.cos_out)-cap(self.cos_in))*(2/pi)
                return fraction*speed_factor/(beta*p.width)

            for k in range(len(row)):
                integral[k] = quad(scalar_density, a[k], b[k], args=(float(q[k]),),
                                   epsabs=2e-12, epsrel=2e-12, limit=250)[0]
        else:
            pending = np.flatnonzero(numerical)
            left, right = np.zeros(len(pending)), np.ones(len(pending))
            for depth in range(90):
                if not len(pending):
                    break
                estimates = []
                for nodes, weights in RULES:
                    t = left[:, None] + (right-left)[:, None]*(nodes+1)/2
                    if method == "semi":
                        position = np.sin(np.pi*t/2)**2
                        jacobian = np.pi/2*np.sin(np.pi*t)
                    else:
                        position, jacobian = t, 1.
                    beta = a[pending, None] + span[pending, None]*position
                    value = self.density(beta, q[pending, None])*jacobian
                    estimates.append((value @ weights)*(right-left)*span[pending]/(2*p.width))
                low, high = estimates
                threshold = tolerance*(span[pending]/p.width/(depth+1)**2 + np.abs(high))
                refine = np.abs(high-low) > threshold
                np.add.at(integral, pending[~refine], high[~refine])
                pending, left, right = pending[refine], left[refine], right[refine]
                middle = (left+right)/2
                pending = np.tile(pending, 2)
                left, right = np.r_[left, middle], np.r_[middle, right]
            if len(pending):
                raise RuntimeError(f"{method} did not converge on {len(pending)} subpanels")
        averaged = np.bincount(row, weights=integral, minlength=len(e))
        result[active] += averaged*(p.norm/(2*p.E0))*(1 if p.photon else e/p.E0)
        return result.reshape(energy.shape)
