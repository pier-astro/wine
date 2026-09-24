import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from wine.xspec import read_table_model
from wine import SpectralGrid


def write_test_table(path):
    from astropy.io import fits

    parameters = fits.BinTableHDU.from_columns(
        [
            fits.Column(name="NAME", format="12A", array=["linear", "logarithmic"]),
            fits.Column(name="METHOD", format="J", array=[0, 1]),
            fits.Column(name="NUMBVALS", format="J", array=[2, 2]),
            fits.Column(
                name="VALUE",
                format="2E",
                array=np.array([[0.0, 2.0], [1.0, 100.0]], dtype=np.float32),
            ),
        ],
        name="PARAMETERS",
    )
    energies = fits.BinTableHDU.from_columns(
        [
            fits.Column(name="ENERG_LO", format="E", array=[1.0, 2.0]),
            fits.Column(name="ENERG_HI", format="E", array=[2.0, 3.0]),
        ],
        name="ENERGIES",
    )
    values = np.array(
        [[0.0, 0.0], [2.0, 2.0], [10.0, 10.0], [12.0, 12.0]], dtype=np.float32
    )
    spectra = fits.BinTableHDU.from_columns(
        [
            fits.Column(
                name="PARAMVAL",
                format="2E",
                array=np.array([[0.0, 1.0], [0.0, 100.0], [2.0, 1.0], [2.0, 100.0]], dtype=np.float32),
            ),
            fits.Column(name="INTPSPEC", format="2E", array=values),
        ],
        name="SPECTRA",
    )
    primary = fits.PrimaryHDU()
    primary.header["ADDMODEL"] = True
    primary.header["REDSHIFT"] = True
    primary.header["XFXP0001"] = "additive"
    fits.HDUList([primary, parameters, energies, spectra]).writeto(path)


class XspecTableTest(unittest.TestCase):
    def test_single_interpolated_parameter(self):
        from astropy.io import fits

        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "single.fits"
            parameters = fits.BinTableHDU.from_columns(
                [
                    fits.Column(name="NAME", format="12A", array=["logxi"]),
                    fits.Column(name="METHOD", format="J", array=[0]),
                    fits.Column(name="NUMBVALS", format="J", array=[2]),
                    fits.Column(name="VALUE", format="2E", array=[[3.0, 4.0]]),
                ], name="PARAMETERS",
            )
            energies = fits.BinTableHDU.from_columns(
                [
                    fits.Column(name="ENERG_LO", format="E", array=[1.0, 2.0]),
                    fits.Column(name="ENERG_HI", format="E", array=[2.0, 3.0]),
                ], name="ENERGIES",
            )
            spectra = fits.BinTableHDU.from_columns(
                [
                    fits.Column(name="PARAMVAL", format="E", array=[3.0, 4.0]),
                    fits.Column(name="INTPSPEC", format="2E", array=[[1.0, 2.0], [3.0, 4.0]]),
                ], name="SPECTRA",
            )
            primary = fits.PrimaryHDU()
            primary.header["ADDMODEL"] = True
            fits.HDUList([primary, parameters, energies, spectra]).writeto(path)
            spectrum = read_table_model(path).evaluate({"logxi": 3.5})
            np.testing.assert_allclose(spectrum.values, [2.0, 3.0])
            grid = SpectralGrid.from_tables(
                emission={"total": path}, fixed_parameters={"vturb": 100.0, "log_nh": 21.0}
            )
            np.testing.assert_allclose(grid.emitted(3.5, 100.0, 21.0), [np.sqrt(3.0), np.sqrt(8.0)])

    def test_xspec_coordinate_and_spectrum_interpolation(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "table.fits"
            write_test_table(path)
            table = read_table_model(path)
            spectrum = table.evaluate(
                {"linear": 1.0, "logarithmic": 10.0},
                redshift=1.0,
                normalization=3.0,
            )

        self.assertTrue(np.allclose(spectrum.values, [9.0, 9.0]))
        self.assertTrue(np.allclose(spectrum.energy_lo, [0.5, 1.0]))
        self.assertTrue(np.allclose(spectrum.energy_hi, [1.0, 1.5]))

    def test_multiplicative_and_exponential_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "table.fits"
            write_test_table(path)
            for model_type, expected in [
                ("multiplicative", [12.0, 12.0]),
                ("exponential", np.exp([-12.0, -12.0])),
            ]:
                from astropy.io import fits

                with fits.open(path, mode="update") as hdul:
                    hdul[0].header["ADDMODEL"] = False
                    hdul[0].header["XFXP0001"] = model_type
                spectrum = read_table_model(path).evaluate(
                    {"linear": 2.0, "logarithmic": 100.0}
                )
                self.assertTrue(np.allclose(spectrum.transmission, expected))


if __name__ == "__main__":
    unittest.main()
