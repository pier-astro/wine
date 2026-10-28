#!/usr/bin/env python3
"""Check finite-grid effects and interior convergence without external tables.

All three log grids share the same 1--10 keV bins and log-energy spacing.
Only exterior input coverage changes. Flat transmission and power-law inputs
isolate the convolution boundary from photoionization-table interpolation.
"""

from pathlib import Path

import numpy as np
from xspec import AllModels, Model, Xset


PACKAGE_DIR = Path(__file__).resolve().parent
GRIDS = (
    ("1 10 1000 log", 0),
    ("0.1 100 3000 log", 1000),
    ("0.01 1000 5000 log", 2000),
)


def evaluate_band(model, grid, first_bin):
    """Return integrated flux in the common 1--10 keV output bins."""
    AllModels.setEnergies(grid)
    edges = np.asarray(model.energies(0))[first_bin:first_bin + 1001]
    # PyXSPEC's exported energy edges are rounded to a finite decimal precision.
    np.testing.assert_allclose(edges, np.geomspace(1.0, 10.0, 1001), rtol=1e-9)
    flux = np.asarray(model.values(0))[first_bin:first_bin + 1000]
    assert flux.shape == (1000,)
    assert np.all(np.isfinite(flux))
    return flux


def main():
    """Exercise blueshift, net redshift, and angular emission broadening."""
    Xset.chatter = 0
    AllModels.lmod("wine", str(PACKAGE_DIR))
    cases = (
        ("(wind*constant)*powerlaw", (0.3, 0.5, 2.0, 1.0)),
        ("(zwind*constant)*powerlaw", (0.3, 1.0, 0.5, 2.0, 1.0)),
        ("awind*powerlaw", (0.3, 2.0, 1.0)),
        ("azwind*powerlaw", (0.3, 1.0, 2.0, 1.0)),
        ("windemlos*powerlaw", (0.3, 0.3, 0.5, 45.0, 2.0, 1.0)),
    )

    for expression, parameters in cases:
        AllModels.clear()
        model = Model(expression)
        model.setPars(*parameters)
        restricted, extended, wider = (
            evaluate_band(model, grid, first) for grid, first in GRIDS
        )
        np.testing.assert_allclose(extended, wider, rtol=1e-9, atol=0.0)
        # A grid ending at the output band cannot supply all shifted photons.
        assert np.max(np.abs(restricted / extended - 1.0)) > 0.1, expression

        # Away from the missing input boundary, constant transmission stays
        # exactly 0.5; integral of E^-2 in each bin is 1/E_low - 1/E_high.
        if expression.startswith("(wind") or expression.startswith("(zwind"):
            edges = np.geomspace(1.0, 10.0, 1001)
            expected = 0.5 * (1.0 / edges[:-1] - 1.0 / edges[1:])
            np.testing.assert_allclose(extended, expected, rtol=1e-9, atol=0.0)

        relative_change = np.max(np.abs(extended / wider - 1.0))
        print(f"{expression}: extension change {relative_change:.2g}")

    print("Energy-coverage checks passed")


if __name__ == "__main__":
    main()
