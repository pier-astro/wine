#!/usr/bin/env python3
"""Check that geometry wrappers reproduce the corresponding base models."""

from math import isclose
from pathlib import Path

from xspec import AllModels, Model, Xset


PACKAGE_DIR = Path(__file__).resolve().parent
ENERGY_GRID = "0.1 20.0 1000 log"


def evaluate(expression, parameters):
    """Evaluate one model on the fixed test grid and return its flux array."""
    AllModels.clear()
    model = Model(expression)
    for index in range(1, model.nParameters + 1):
        if model(index).name == "region":
            assert model(index).frozen, "region must remain a frozen XSPEC switch"
    for index, value in enumerate(parameters, start=1):
        if model(index).values[0] != value:
            model(index).values = str(value)
    return tuple(model.values(0))


def assert_same(actual, expected, label):
    """Require numerical identity to round-off across the full energy grid."""
    assert len(actual) == len(expected)
    assert all(
        isclose(left, right, rel_tol=1.0e-13, abs_tol=1.0e-15)
        for left, right in zip(actual, expected)
    ), label


def main():
    """Exercise in-cone and both off-cone model parameterizations."""
    Xset.chatter = 0
    AllModels.lmod("wine", str(PACKAGE_DIR))
    AllModels.setEnergies(ENERGY_GRID)

    beta = 0.2
    uview = 0.25
    tout = 60.0
    tin = 30.0
    line_energy = 6.5

    for uincl, inclination in ((0.0, tin), (0.25, 37.5), (1.0, tout)):
        los_convolution = evaluate(
            "windemlos*powerlaw", (beta, uincl, tout, tin, 2.0, 1.0)
        )
        base_convolution = evaluate(
            "windem*powerlaw", (beta, inclination, tout, tin, 2.0, 1.0)
        )
        assert_same(los_convolution, base_convolution, f"windemlos uincl={uincl}")

        los_line = evaluate(
            "windlinelos", (beta, uincl, tout, tin, line_energy, 1.0)
        )
        base_line = evaluate(
            "windline", (beta, inclination, tout, tin, line_energy, 1.0)
        )
        assert_same(los_line, base_line, f"windlinelos uincl={uincl}")

    for region, inclination in ((0, 22.5), (1, 67.5)):
        off_convolution = evaluate(
            "windemoff*powerlaw",
            (beta, uview, tout, tin, region, 2.0, 1.0),
        )
        base_convolution = evaluate(
            "windem*powerlaw",
            (beta, inclination, tout, tin, 2.0, 1.0),
        )
        assert_same(
            off_convolution,
            base_convolution,
            f"windemoff region={region} differs from windem at incl={inclination}",
        )

        off_line = evaluate(
            "windlineoff",
            (beta, uview, tout, tin, region, line_energy, 1.0),
        )
        base_line = evaluate(
            "windline",
            (beta, inclination, tout, tin, line_energy, 1.0),
        )
        assert_same(
            off_line,
            base_line,
            f"windlineoff region={region} differs from windline at incl={inclination}",
        )

    print("windem/windline geometry checks passed")


if __name__ == "__main__":
    main()
