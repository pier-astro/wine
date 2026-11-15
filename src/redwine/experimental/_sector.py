"""Analytical observer-ring intersections for a constant radial velocity."""

from dataclasses import dataclass

import numpy as np

from ._uniform import Prepared, Profile

_TAU = 2*np.pi


@dataclass(frozen=True)
class Sector:
    """A signed radial velocity in units of c, and angles in degrees.

    Phi is measured from the observer's equatorial projection and is periodic.
    A bicone adds the antipodal copy of the specified sector.
    """

    vout: float
    tin: float
    tout: float
    incl: float
    phi: float
    dphi: float
    bicone: bool = False

    def __post_init__(self):
        if (not np.isfinite([self.vout, self.tin, self.tout,
                             self.incl, self.phi, self.dphi]).all()
                or not abs(self.vout) < 1
                or self.bicone not in (0, 1)
                or not 0 <= self.tin <= self.tout <= (90 if self.bicone else 180)
                or not 0 <= self.incl <= 180 or not 0 <= self.dphi <= 360):
            raise ValueError("Invalid radial velocity or sector geometry.")
        object.__setattr__(self, "phi", self.phi % 360)

    @property
    def solid_angle_fraction(self):
        inner, outer = np.deg2rad([self.tin, self.tout])
        return (np.sin((outer+inner)/2)*np.sin((outer-inner)/2)
                * (1+self.bicone)*(self.dphi/360))


def _sincos(degrees):
    # Exact cardinal directions keep coincident planes from gaining a sliver.
    quadrant = degrees % 360
    if quadrant in (0, 90, 180, 270):
        return ((0., 1.), (1., 0.), (0., -1.), (-1., 0.))[int(quadrant/90)]
    angle = np.deg2rad(degrees)
    return np.sin(angle), np.cos(angle)


class SectorRings:
    def __init__(self, sector):
        self.sector = sector
        self.cap = Prepared(Profile(1., abs(sector.vout), 0.,
                                    *np.deg2rad([sector.incl, sector.tin, sector.tout])))
        self.aligned = sector.incl in (0, 180)
        self.complement = sector.dphi > 180
        center = sector.phi + (180 if self.complement else 0)
        half = min(sector.dphi, 360-sector.dphi)/2
        sin_i, cos_i = _sincos(sector.incl)
        knots = list(self.cap.tangencies)
        self.planes = []
        for edge, sign in ((center-half, 1), (center+half, -1)):
            sin_phi, cos_phi = _sincos(edge)
            # Inward meridian normal in the wind frame: sign*(-sin phi, cos phi, 0).
            nx, ny = -sign*sin_phi, sign*cos_phi
            longitudinal = nx*sin_i
            transverse = np.hypot(nx*cos_i, ny)
            # Match the cap roots: observer azimuth zero points toward +wind z.
            phase = np.arctan2(ny, -nx*cos_i) % _TAU
            self.planes.append((longitudinal, transverse, phase))
            knots.extend((-transverse, transverse))
            for theta in (sector.tin, sector.tout):
                sin_t, cos_t = _sincos(theta)
                knots.append(cos_i*cos_t + sin_i*sin_t*cos_phi)
        if sector.bicone:
            knots.extend(-np.asarray(knots))
        self.knots = np.clip(knots, -1, 1)

    @staticmethod
    def _halfspace(mu, sine, plane):
        longitudinal, transverse, phase = plane
        # In the observer frame the boundary is
        # longitudinal*mu + transverse*sine*cos(psi-phase) = 0.
        offset = longitudinal*mu
        radius = transverse*sine
        threshold = np.divide(-offset, radius, out=np.zeros_like(mu), where=radius > 0)
        half = np.where(radius > 0, np.arccos(np.clip(threshold, -1, 1)),
                        np.pi*(offset >= 0))
        left, right = phase-half, phase+half
        return ((np.maximum(0., left), np.minimum(_TAU, right)),
                (np.where(left < 0, left+_TAU, 0.),
                 np.where(left < 0, _TAU, np.maximum(0., right-_TAU))))

    def _one_cap(self, mu):
        if self.aligned:
            return self.cap.arcs(mu)*self.sector.dphi/360

        inner, outer = self.cap.polar_limits(mu)
        polar = ((inner, outer), (_TAU-outer, _TAU-inner))
        sine = np.sqrt(np.maximum(0., (1-mu)*(1+mu)))
        first, second = [self._halfspace(mu, sine, plane) for plane in self.planes]
        length = np.zeros_like(mu)
        for a, b in polar:
            for c, d in first:
                for e, f in second:
                    length += np.maximum(0., np.minimum(np.minimum(b, d), f)
                                         - np.maximum(np.maximum(a, c), e))
        fraction = length/_TAU
        if self.complement:
            fraction = np.maximum(0., (outer-inner)/np.pi-fraction)
        return fraction

    def arcs(self, mu):
        """Fraction of each constant-Doppler ring inside the emitting region."""
        fraction = self._one_cap(mu)
        if self.sector.bicone:
            fraction += self._one_cap(-mu)
        return np.where(np.abs(mu) <= 1, fraction, 0.)
