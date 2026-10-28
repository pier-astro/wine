"""Check the geometric weight alone, without Doppler factors or energy grids.

Run from the repository root: python examples/check_geometry_normalization.py
The input flux refers to a full, uniformly emitting sphere (4 pi sr).
For mu = cos(theta) in the observer frame, dOmega = dphi dmu, hence
integral w_phi(mu) dmu / (4 pi) = (cos(tin) - cos(tout)) / 2.
This is also the nonrelativistic integrated line-profile limit.
"""

import numpy as np
from scipy.integrate import quad

from wine.kernels import w_phi_hollowcap


# Angles in degrees. The full sphere is a mathematical reference; the public
# emission geometries describe at most one hemisphere (tout <= 90 degrees).
regions = [
    ("sphere", 0.0, 180.0),
    ("hemisphere", 0.0, 90.0),
    ("cap", 0.0, 60.0),
    ("hollow cap", 30.0, 70.0),
    ("empty", 45.0, 45.0),
]
inclinations = [0.0, 25.0, 90.0]

print("Geometry       tin   tout   incl     Integral / 4pi       Analytic")
for name, tin_deg, tout_deg in regions:
    tin, tout = np.deg2rad([tin_deg, tout_deg])
    expected = (np.cos(tin) - np.cos(tout)) / 2.0
    for incl_deg in inclinations:
        incl = np.deg2rad(incl_deg)

        # Ring/cap tangencies mark changes in the azimuthal weight. Supply
        # them to adaptive quadrature so narrow angular regions are resolved.
        tangencies = np.cos([incl - tin, incl + tin,
                             incl - tout, incl + tout])
        points = np.unique(tangencies[(tangencies > -1) & (tangencies < 1)])
        solid_angle, error = quad(
            lambda mu: float(w_phi_hollowcap(np.arccos(mu), incl, tin, tout)),
            -1.0, 1.0, points=points, epsabs=1e-10, epsrel=1e-10,
        )
        fraction = solid_angle / (4.0 * np.pi)
        print(f"{name:12s} {tin_deg:5.0f} {tout_deg:6.0f} {incl_deg:6.0f}"
              f" {fraction:18.10f} {expected:14.10f}")
        np.testing.assert_allclose(fraction, expected, atol=1e-9, rtol=0)

print("All geometric integrals agree with the analytic solid angle.")
