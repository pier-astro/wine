#!/usr/bin/env python3
"""Small native fit and rectangular-grid check, without external data/tables.

Expected counts (no Poisson noise) test optimizer behavior, not significance.
The fixed line normalization makes hemispheric Cf identifiable in this toy.
The required recovery starts away from the generating values. This is not
a global-convergence check for this degenerate spectral model.
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
    AllData.dummyrsp(0.5, 12.0, 120, "lin", 0.5, (12.0 - 0.5) / 120)
    model = Model("powerlaw+windemlos*gaussian")
    model.setPars(2.0, 0.01, 0.25, 0.4, 0.65, 35.0, 6.7, 0.025, 0.005)
    AllModels.setEnergies("0.1 30 2000 log")
    AllData.fakeit(1, FakeitSettings(exposure=100000), applyStats=False, noWrite=True)
    Fit.statMethod = "cstat"
    Fit.query = "no"
    Fit.nIterations = 40
    Fit.criticalDelta = 1e-5
    for index in range(1, 10):
        model(index).frozen = index not in (2, 4, 5, 6)

    # Constructor-like starting geometry differs from the generating one.
    model.setPars({2: 0.008, 4: 0.3, 5: 0.5, 6: 45.0})
    Fit.perform()
    fitted = [model(i).values[0] for i in (4, 5, 6)]
    geometry = Geometry.los(Cf=fitted[0], skew=fitted[1], incl=fitted[2])
    assert geometry.tin <= geometry.incl <= geometry.tout
    assert np.isfinite(Fit.statistic) and Fit.statistic < 1e-3
    assert abs(fitted[0] - 0.4) < 0.003
    print(f"Known-spectrum recovery: C={Fit.statistic:.6g}, fitted={fitted}")

    Fit.error("1.0 4")
    lower, upper, status = model(4).error
    assert 0.0 <= lower < model(4).values[0] < upper <= 1.0
    print(f"Native Cf profile (delta-stat=1): {lower:.6g}, {upper:.6g}; flags={status}")

    # Include physical zero/full-coverage limits, and pole/equator sightlines.
    # The nuisance parameters are refitted at each point, not simply evaluated.
    Fit.steppar("best 4 0 1 2 6 0 90 2")
    surface = np.asarray(Fit.stepparResults("statistic"))
    assert surface.shape == (9,)
    assert np.all(np.isfinite(surface)) and np.max(surface) > 0.0
    print("Native profiled 3x3 (Cf,incl) scan completed, including the boundaries")


if __name__ == "__main__":
    main()
