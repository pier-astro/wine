import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from wine import SpectralGrid, absorption_slab, emission_slab
from test_ogip_table import write_test_table


class SpectralGridTest(unittest.TestCase):
    def setUp(self):
        energy = np.linspace(0.5, 5.0, 10)
        values = np.ones((2, 2, 2, 10))
        values[1] = 100.0
        self.args = (energy, [3.0, 4.0], [100.0, 200.0], [21.0, 22.0])
        self.values = values

    def test_log_default_and_linear_option(self):
        log_grid = SpectralGrid(*self.args, optical_depth={"absorption": self.values}, emission={"total": self.values})
        linear_grid = SpectralGrid(*self.args, optical_depth={"absorption": self.values}, interpolation="linear")
        np.testing.assert_allclose(log_grid.optical_depth(3.5, 100.0, 21.0), 10.0)
        np.testing.assert_allclose(log_grid.emitted(3.5, 100.0, 21.0), 10.0)
        np.testing.assert_allclose(linear_grid.optical_depth(3.5, 100.0, 21.0), 50.5)

        zero = self.values.copy()
        zero[0] = 0.0
        log_zero = SpectralGrid(*self.args, optical_depth={"absorption": zero})
        linear_zero = SpectralGrid(*self.args, optical_depth={"absorption": zero}, interpolation="linear")
        self.assertGreater(log_zero.optical_depth(3.5, 100.0, 21.0)[0], 0.0)
        self.assertEqual(linear_zero.optical_depth(3.0, 100.0, 21.0)[0], 0.0)

    def test_models_use_engine_independent_grid(self):
        grid = SpectralGrid(*self.args, optical_depth={"absorption": self.values}, emission={"total": self.values})
        transmission = absorption_slab(np.array([1.0, 2.0, 3.0]), grid, 3.5, 21.0, 100.0, 0.0, 1.0)
        np.testing.assert_allclose(transmission, np.exp(-10.0))
        flux = emission_slab(
            np.linspace(1.0, 4.0, 20), grid, 3.5, 21.0, 100.0,
            0.1, 0.5, 0.0, 1.0, 1.0,
        )
        self.assertEqual(flux.shape, (20,))
        self.assertTrue(np.all(np.isfinite(flux)))
        self.assertFalse(any(name.startswith("xflow.cloudy") for name in sys.modules))

    def test_ogip_tables_and_fixed_axis(self):
        from astropy.io import fits

        with tempfile.TemporaryDirectory() as tmp:
            emission = Path(tmp) / "emission.fits"
            depth = Path(tmp) / "depth.fits"
            for path, model_type in ((emission, "additive"), (depth, "exponential")):
                write_test_table(path)
                with fits.open(path, mode="update") as hdul:
                    hdul["PARAMETERS"].data["NAME"][:] = ["logxi", "vturb"]
                    hdul["SPECTRA"].data["INTPSPEC"][:] = [[1, 1], [1, 1], [100, 100], [100, 100]]
                    hdul[0].header["XFXP0001"] = model_type
            grid = SpectralGrid.from_tables(
                optical_depth={"absorption": depth}, emission={"total": emission},
                fixed_parameters={"log_nh": 21.0},
            )
            np.testing.assert_allclose(grid.optical_depth(1.0, 10.0, 21.0), 10.0)
            np.testing.assert_allclose(grid.emitted(1.0, 10.0, 21.0), 10.0)
            linear = SpectralGrid.from_tables(
                optical_depth={"absorption": depth}, fixed_parameters={"log_nh": 21.0},
                interpolation="linear",
            )
            np.testing.assert_allclose(linear.optical_depth(1.0, 10.0, 21.0), 50.5)


if __name__ == "__main__":
    unittest.main()
