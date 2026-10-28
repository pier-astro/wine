"""Hemispheric wind geometry; public angles are in degrees.

Cf is Omega/(2*pi), so the existing emission normalization is Cf/2 relative
to a full sphere. The constructors change coordinates, not the wind physics.
"""

from dataclasses import dataclass
import warnings

import numpy as np


def _bounded(value, low, high, name):
    value = np.asarray(value, dtype=float)
    if np.any(~np.isfinite(value) | (value < low) | (value > high)):
        raise ValueError(f"{name} must be finite and in [{low}, {high}].")
    return value


@dataclass(frozen=True)
class Geometry:
    """Physical inclination and shell boundaries in degrees; arrays broadcast.

    ``angles`` also permits inverted shell boundaries as an explicit empty
    emission convention, with a warning. Fitting coordinates use ordered shells.
    """

    incl: float | np.ndarray
    tin: float | np.ndarray
    tout: float | np.ndarray

    def __post_init__(self):
        values = np.broadcast_arrays(
            _bounded(self.incl, 0.0, 90.0, "incl"),
            _bounded(self.tin, 0.0, 90.0, "tin"),
            _bounded(self.tout, 0.0, 90.0, "tout"),
        )
        for name, value in zip(("incl", "tin", "tout"), values):
            object.__setattr__(self, name, float(value) if value.ndim == 0 else value)
        if np.any(values[1] > values[2]):
            warnings.warn("tin > tout: the angular interval is empty; emission is zero.",
                          RuntimeWarning, stacklevel=2)

    @property
    def Cf(self):
        """Occupied fraction of one hemisphere, zero for an inverted interval."""
        inner, outer = np.deg2rad(self.tin), np.deg2rad(self.tout)
        # Stable cosine difference for a thin shell.
        return np.clip(2.0 * np.sin((outer + inner) / 2.0)
                       * np.sin((outer - inner) / 2.0), 0.0, 1.0)

    @classmethod
    def angles(cls, *, incl, tin, tout):
        """Specify the physical angles directly, in degrees."""
        return cls(incl=incl, tin=tin, tout=tout)

    @classmethod
    def free(cls, *, Cf, hole, incl):
        """Split unused hemispheric solid angle with hole; incl is independent."""
        cover, cavity, inclination = np.broadcast_arrays(
            _bounded(Cf, 0.0, 1.0, "Cf"),
            _bounded(hole, 0.0, 1.0, "hole"),
            _bounded(incl, 0.0, 90.0, "incl"),
        )
        cos_inner = 1.0 - cavity * (1.0 - cover)
        cos_outer = (1.0 - cavity) * (1.0 - cover)
        return cls(incl=inclination, tin=np.rad2deg(np.arccos(cos_inner)),
                   tout=np.rad2deg(np.arccos(cos_outer)))

    @classmethod
    def los(cls, *, Cf, skew, incl):
        """Construct a shell containing incl, with polar/equatorial skew.

        For x=cos(tout), c=cos(incl), and s=skew, the physical root satisfies
        (1-s)*(1-x-Cf)*(c-x) = s*x*(x+Cf-c).
        At s=0.5 both sides of the sightline have fractional coverage Cf.
        """
        cover, skew, inclination = np.broadcast_arrays(
            _bounded(Cf, 0.0, 1.0, "Cf"),
            _bounded(skew, 0.0, 1.0, "skew"),
            _bounded(incl, 0.0, 90.0, "incl"),
        )
        c = np.where(inclination == 90.0, 0.0, np.cos(np.deg2rad(inclination)))
        a = 1.0 - 2.0 * skew
        mix = (1.0 - skew) * (1.0 - cover) + skew * cover
        b = -(mix + a * c)
        constant = (1.0 - skew) * c * (1.0 - cover)
        # Positive-term discriminant avoids cancellation at the double roots.
        discriminant = (mix - c)**2 + 4.0 * skew * (1.0 - skew) * c * (1.0 - c)
        root = np.sqrt(discriminant)
        # The rationalized root also covers the linear case skew=0.5.
        # When b>0, the direct form avoids a vanishing rational denominator.
        x = np.zeros_like(c)
        np.divide(2.0 * constant, -b + root, out=x, where=(b <= 0) & (-b + root > 0))
        np.divide(-b - root, 2.0 * a, out=x, where=b > 0)
        x = np.clip(x, np.maximum(0.0, c - cover), np.minimum(c, 1.0 - cover))
        return cls(incl=inclination,
                   tin=np.rad2deg(np.arccos(np.clip(x + cover, 0.0, 1.0))),
                   tout=np.rad2deg(np.arccos(np.clip(x, 0.0, 1.0))))

    @classmethod
    def off(cls, *, Cf, hole, view, region):
        """Place the sightline in the polar or equatorial off-wind sector.

        view=0/1 always selects the lower/upper inclination of that sector.
        At a collapsed sector the inclination is its single boundary value.
        """
        if region not in ("polar", "equatorial"):
            raise ValueError("region must be 'polar' or 'equatorial'.")
        shell = cls.free(Cf=Cf, hole=hole, incl=0.0)
        view = _bounded(view, 0.0, 1.0, "view")
        inclination = (view * shell.tin if region == "polar"
                       else shell.tout + view * (90.0 - shell.tout))
        return cls(incl=inclination, tin=shell.tin, tout=shell.tout)
