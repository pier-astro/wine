#!/usr/bin/env python3
"""Prototype a fixed-opening-angle LoS fit using native XSPEC links/bounds.

Run in a separate initialized PyXSPEC process. Expected counts, no external
tables: this checks a constrained fit, not statistical significance or the
full angular scan API. The example deliberately avoids singular geometries.
"""
from pathlib import Path
import sys

import numpy as np
from xspec import AllData, AllModels, FakeitSettings, Fit, Model, Xset

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from wine import Geometry


def main():
    Xset.chatter = 0
    Xset.allowPrompting = False
    AllModels.lmod("wine", str(Path(__file__).resolve().parent))
    AllData.dummyrsp(0.5, 12., 120, "lin", 0.5, 11.5 / 120)
    model = Model("powerlaw+windemlos*gaussian")
    AllModels.setEnergies("0.1 30 2000 log")

    tin, tout, true_incl = 20., 60., 35.
    r, x = np.cos(np.deg2rad([tin, tout]))
    cover = r - x
    c = np.cos(np.deg2rad(true_incl))
    a, b = (1-r)*(c-x), x*(r-c)
    model.setPars(2., .01, .25, cover, a/(a+b), true_incl,
                  6.7, .025, .005)
    AllData.fakeit(1, FakeitSettings(exposure=100000), applyStats=False,
                   noWrite=True)
    Fit.statMethod = "cstat"
    Fit.query = "no"
    Fit.nIterations = 40
    Fit.criticalDelta = 1e-6

    # For fixed angles, Cf is constant and skew=A/(A+B), where
    # A=(1-cos(tin))*(cos(incl)-cos(tout)),
    # B=cos(tout)*(cos(tin)-cos(incl)). XSPEC's cosd uses degrees.
    for index in range(1, 10):
        model(index).frozen = True
    model(6).values = [45., .01, tin, tin, tout, tout]
    model(6).frozen = False
    model(2).values = .008
    model(2).frozen = False
    # An integer literal in an XSPEC link denotes a parameter, so use 1.0.
    a_link = f"(1.0-{r:.17g})*(cosd(p6)-{x:.17g})"
    b_link = f"{x:.17g}*({r:.17g}-cosd(p6))"
    model(5).link = f"({a_link})/(({a_link})+({b_link}))"

    # Verify actual XSPEC link evaluation, including both allowed endpoints,
    # before fitting. This is not merely a Python inverse-formula test.
    for inclination in (tin, 25., true_incl, 45., 55., tout):
        model(6).values = inclination
        skew = model(5).values[0]
        print(f"Link at incl={inclination:g}: skew={skew:.17g}", flush=True)
        assert -1e-12 <= skew <= 1.0 + 1e-12
        # XSPEC's degree conversion may round infinitesimally beyond an
        # endpoint; the model adapter already bounds such numerical probes.
        geometry = Geometry.los(Cf=model(4).values[0],
                                skew=np.clip(skew, 0., 1.),
                                incl=model(6).values[0])
        np.testing.assert_allclose([geometry.tin, geometry.tout],
                                   [tin, tout], atol=1e-7, rtol=0)
    print("Native link preserves tin=20, tout=60 at six inclinations/endpoints")

    model(6).values = 45.
    Fit.perform()
    inclination = model(6).values[0]
    geometry = Geometry.los(Cf=model(4).values[0],
                            skew=model(5).values[0], incl=inclination)
    np.testing.assert_allclose([geometry.tin, geometry.tout],
                               [tin, tout], atol=1e-7, rtol=0)
    assert tin <= inclination <= tout
    assert abs(inclination - true_incl) < .05
    assert np.isfinite(Fit.statistic) and Fit.statistic < 1e-3
    print(f"Constrained fit: C={Fit.statistic:.8g}, incl={inclination:.8f}, "
          f"tin={geometry.tin:.8f}, tout={geometry.tout:.8f}")


if __name__ == "__main__":
    main()
