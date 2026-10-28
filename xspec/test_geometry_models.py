#!/usr/bin/env python3
"""Check new coordinates against generic spectra at identical physical angles.

Run after building the local package. Python Geometry also has independent
bracketed-root tests in tests/test_geometry.py. This checks the C++/XSPEC
conversion, endpoint behavior, registration, and hemispheric normalization.
"""
from pathlib import Path
import sys

import numpy as np
from xspec import AllModels, Model, Xset

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from wine import Geometry

PACKAGE_DIR = Path(__file__).resolve().parent


def evaluate(expression, parameters):
    AllModels.clear()
    model = Model(expression)
    model.setPars(*parameters)
    for index in range(1, model.nParameters + 1):
        if model(index).name == "region":
            assert model(index).frozen
    flux = np.asarray(model.values(0))
    assert np.all(np.isfinite(flux)), (expression, parameters)
    return flux


def generic_parameters(beta, geometry):
    cover = float(geometry.Cf)
    hole = ((1.0 - np.cos(np.deg2rad(geometry.tin))) / (1.0 - cover)
            if cover < 1.0 else 0.5)
    return [beta, cover, float(np.clip(hole, 0.0, 1.0)), float(geometry.incl)]


def main():
    Xset.chatter = 0
    AllModels.lmod("wine", str(PACKAGE_DIR))
    AllModels.setEnergies("0.1 30 2000 log")
    for cover in (0.0, 0.3, 1.0):
        for skew in (0.0, 0.5, 1.0):
            for inclination in (0.0, 45.0, 90.0):
                geometry = Geometry.los(Cf=cover, skew=skew, incl=inclination)
                for beta in (0.0, 0.25):
                    parameters = [beta, cover, skew, inclination]
                    physical = generic_parameters(beta, geometry)
                    for family, suffix in (("windem", "*gaussian"), ("windline", "")):
                        extra = [6.7, 0.025, 1.0] if suffix else [6.7, 1.0]
                        actual = evaluate(family + "los" + suffix, parameters + extra)
                        expected = evaluate(family + suffix, physical + extra)
                        np.testing.assert_allclose(actual, expected, rtol=1e-9, atol=1e-13)

    for region, name in ((0, "polar"), (1, "equatorial")):
        for hole in (0.0, 0.5, 1.0):
            for view in (0.0, 0.5, 1.0):
                geometry = Geometry.off(Cf=0.3, hole=hole, view=view, region=name)
                for family, suffix in (("windem", "*gaussian"), ("windline", "")):
                    extra = [6.7, 0.025, 1.0] if suffix else [6.7, 1.0]
                    actual = evaluate(family + "off" + suffix, [0.25, 0.3, hole, view, region] + extra)
                    expected = evaluate(family + suffix, generic_parameters(0.25, geometry) + extra)
                    np.testing.assert_allclose(actual, expected, rtol=1e-9, atol=1e-13)

    # No Doppler shift: a narrow line integrates to Cf/2 of its input norm.
    for family, suffix in (("windem", "*gaussian"), ("windline", "")):
        extra = [6.7, 0.025, 1.0] if suffix else [6.7, 1.0]
        flux = evaluate(family + "los" + suffix, [0.0, 0.3, 0.5, 45.0] + extra)
        np.testing.assert_allclose(flux.sum(), 0.15, rtol=1e-7)

    # Reusing a model must not leave an earlier nonzero spectrum at Cf=0.
    for expression, extra in (("windemlos*powerlaw", [2.0, 1.0]),
                              ("windlinelos", [6.7, 1.0])):
        AllModels.clear()
        model = Model(expression)
        model.setPars(*([0.25, 0.3, 0.5, 45.0] + extra))
        reference = np.asarray(model.values(0))
        model(2).values = "0"
        assert np.all(np.asarray(model.values(0)) == 0.0)
        model(2).values = "0.3"
        np.testing.assert_allclose(model.values(0), reference, rtol=1e-12, atol=1e-15)
        assert [model(i).name for i in range(1, 5)] == ["vout", "Cf", "skew", "incl"]
    print("Geometry and normalization checks passed")


if __name__ == "__main__":
    main()
