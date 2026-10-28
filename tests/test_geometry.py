import unittest

import numpy as np
from scipy.optimize import brentq

from wine import Geometry


class GeometryTest(unittest.TestCase):
    def test_free_limits_and_hemisphere(self):
        full = Geometry.free(Cf=1.0, hole=[0.0, 0.5, 1.0], incl=45.0)
        np.testing.assert_allclose(full.tin, 0.0)
        np.testing.assert_allclose(full.tout, 90.0)
        np.testing.assert_allclose(full.Cf, 1.0)
        cap = Geometry.free(Cf=0.5, hole=0.0, incl=20.0)
        np.testing.assert_allclose([cap.tin, cap.tout], [0.0, 60.0], atol=1e-12)
        belt = Geometry.free(Cf=0.5, hole=1.0, incl=20.0)
        np.testing.assert_allclose([belt.tin, belt.tout], [60.0, 90.0], atol=1e-12)

    def test_los_against_bracketed_reference(self):
        # The reference solves the geometric equation without the analytic
        # quadratic formula, including skew's linear and endpoint cases.
        for cover in (0.0, 1e-8, 0.3, 0.8, 1.0 - 1e-8, 1.0):
            for inclination in (0.0, 1.0, 45.0, 89.0, 90.0):
                c = 0.0 if inclination == 90.0 else np.cos(np.deg2rad(inclination))
                lower, upper = max(0.0, c - cover), min(c, 1.0 - cover)
                previous = upper
                for skew in (0.0, 0.1, 0.5 - 1e-9, 0.5, 0.5 + 1e-9, 0.9, 1.0):
                    geometry = Geometry.los(Cf=cover, skew=skew, incl=inclination)
                    def residual(x):
                        return ((1-skew)*(1-x-cover)*(c-x)
                                - skew*x*(x+cover-c))
                    if lower == upper:
                        expected = lower
                    elif skew == 0.0:
                        expected = upper
                    elif skew == 1.0:
                        expected = lower
                    else:
                        expected = brentq(residual, lower, upper, xtol=1e-15)
                    actual = np.cos(np.deg2rad(geometry.tout))
                    self.assertAlmostEqual(actual, expected, delta=3e-15)
                    self.assertLessEqual(actual, previous + 1e-14)
                    self.assertLessEqual(geometry.tin, inclination + 1e-10)
                    self.assertGreaterEqual(geometry.tout, inclination - 1e-10)
                    self.assertAlmostEqual(float(geometry.Cf), cover, delta=1e-14)
                    previous = actual

    def test_neutral_skew_and_array_shapes(self):
        cover = np.array([0.1, 0.3, 0.8])[:, None]
        incl = np.array([20.0, 45.0, 70.0])
        geometry = Geometry.los(Cf=cover, skew=0.5, incl=incl)
        c = np.cos(np.deg2rad(incl))
        x = np.cos(np.deg2rad(geometry.tout))
        np.testing.assert_allclose((x + cover - c) / (1.0 - c), np.broadcast_to(cover, (3, 3)))
        np.testing.assert_allclose((c - x) / c, np.broadcast_to(cover, (3, 3)))

    def test_off_sector_endpoints(self):
        for region in ("polar", "equatorial"):
            geometry = Geometry.off(Cf=0.3, hole=0.5, view=[0.0, 0.5, 1.0], region=region)
            if region == "polar":
                self.assertEqual(geometry.incl[0], 0.0)
                self.assertEqual(geometry.incl[-1], geometry.tin[-1])
            else:
                self.assertEqual(geometry.incl[0], geometry.tout[0])
                self.assertEqual(geometry.incl[-1], 90.0)

    def test_explicit_empty_geometry_and_validation(self):
        with self.assertWarnsRegex(RuntimeWarning, "tin > tout"):
            geometry = Geometry.angles(incl=45, tin=60, tout=30)
        self.assertEqual(geometry.Cf, 0.0)
        for value in (-0.01, 1.01, np.nan):
            with self.assertRaises(ValueError):
                Geometry.los(Cf=value, skew=0.5, incl=45)
        with self.assertRaises(ValueError):
            Geometry.off(Cf=0.3, hole=0.5, view=0.5, region="wind")


if __name__ == "__main__":
    unittest.main()
