#!/usr/bin/env python3
"""Small fixed-angle scan prototype: limits and identity-linked emitters.

Run in an isolated initialized PyXSPEC process. This is not the public Tcl
procedure: it tests the constraints before adding session/export machinery.
Expected counts contain no Poisson fluctuations.
"""
from pathlib import Path
import sys

import numpy as np
from xspec import AllData, AllModels, FakeitSettings, Fit, Model, Xset

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from wine import Geometry


def constrain(model, tin, tout):
    """Fix opening angles, retaining inclination as a nuisance when defined."""
    if tin > tout:
        raise ValueError("excluded angular pair")
    model(5).untie()
    model(6).values = [35., .01, 0., 0., 90., 90.]
    model(6).frozen = True
    model(5).frozen = True
    model(4).frozen = True
    r, x = np.cos(np.deg2rad([tin, tout]))
    cover = r - x
    model(4).values = cover

    # These are exact geometric limits, not epsilon replacements for 0/0.
    # A zero-width LoS shell fixes incl to its boundary but emits no flux.
    if tin == tout:
        model(5).values = .5
        model(6).values = tin
        return
    model(6).values = [float(np.clip(45., tin, tout)), .01,
                       tin, tin, tout, tout]
    model(6).frozen = False
    if tin == 0. and tout == 90.:
        model(5).values = .5  # All skew values describe the same hemisphere.
    elif tin == 0.:
        model(5).values = 0.
    elif tout == 90.:
        model(5).values = 1.
    else:
        a = f"(1.0-{r:.17g})*(cosd(p6)-{x:.17g})"
        b = f"{x:.17g}*({r:.17g}-cosd(p6))"
        model(5).link = f"({a})/(({a})+({b}))"


def main():
    Xset.chatter = 0
    Xset.allowPrompting = False
    AllModels.lmod("wine", str(Path(__file__).resolve().parent))
    AllData.dummyrsp(.5, 12., 120, "lin", .5, 11.5/120)
    model = Model("powerlaw+windemlos*gaussian+windemlos*gaussian")
    model.setPars(2., .01, .25, .4396926207859083, .242050121690,
                  35., 6.7, .025, .005,
                  .25, .4396926207859083, .242050121690,
                  35., 6.7, .025, .0025)
    AllModels.setEnergies("0.1 30 2000 log")
    # The second component shares the first component's physical geometry.
    for slave, master in ((10, 3), (11, 4), (12, 5), (13, 6)):
        model(slave).link = f"p{master}"
    AllData.fakeit(1, FakeitSettings(exposure=100000), applyStats=False,
                   noWrite=True)
    Fit.statMethod = "cstat"
    Fit.query = "no"
    Fit.nIterations = 40
    Fit.criticalDelta = 1e-6
    for index in range(1, 17):
        if index not in (10, 11, 12, 13):
            model(index).frozen = True
    model(2).frozen = False

    # Exercise all singular angular corners and both cap/belt boundaries.
    for tin, tout in ((0., 0.), (90., 90.), (0., 90.), (0., 60.),
                      (20., 90.), (20., 20.), (20., 60.)):
        constrain(model, tin, tout)
        for inclination in np.linspace(tin, tout, 3):
            model(6).values = inclination
            for indices in ((4, 5, 6), (11, 12, 13)):
                cover, skew, incl = [model(i).values[0] for i in indices]
                assert -1e-11 <= skew <= 1. + 1e-11
                geometry = Geometry.los(Cf=cover, skew=np.clip(skew, 0., 1.),
                                        incl=incl)
                np.testing.assert_allclose([geometry.tin, geometry.tout],
                                           [tin, tout], atol=2e-6, rtol=0)
            assert np.isfinite(Fit.statistic)
    print("Seven angular limits/interior cases preserve both linked geometries")

    cells = [(tin, tout) for tin in (0., 20., 60.)
             for tout in (20., 60., 90.)]
    surfaces = []
    for order in (cells, cells[::-1]):
        results = {}
        for tin, tout in order:
            if tin > tout:
                results[tin, tout] = ("excluded", np.nan)
                continue
            # Same accepted reference for nuisance starts at every cell.
            constrain(model, tin, tout)
            model(2).values = .008
            Fit.perform()
            statistic = Fit.statistic
            assert np.isfinite(statistic)
            for indices in ((4, 5, 6), (11, 12, 13)):
                cover, skew, incl = [model(i).values[0] for i in indices]
                geometry = Geometry.los(Cf=cover, skew=np.clip(skew, 0., 1.),
                                        incl=incl)
                np.testing.assert_allclose([geometry.tin, geometry.tout],
                                           [tin, tout], atol=2e-6, rtol=0)
            results[tin, tout] = ("ok", statistic)
        surfaces.append(results)
    for cell in cells:
        forward, reverse = [surface[cell] for surface in surfaces]
        assert forward[0] == reverse[0]
        np.testing.assert_allclose(forward[1], reverse[1], atol=1e-4,
                                   rtol=0., equal_nan=True)
        print(f"tin,tout={cell}: {forward[0]}, C={forward[1]:.8g}")
    assert surfaces[0][20., 60.][1] < 1e-3
    print("3x3 scan: eight finite fits, one exclusion; reverse order agrees")


if __name__ == "__main__":
    main()
