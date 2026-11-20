"""Measure XSPEC emission-model latency with changing geometry and a real response.

Run in an initialized HEASoft environment. Times include table interpolation,
kernel construction, convolution and folding. The Python reference additionally
includes the XSPEC-to-Python bridge.
All models use the same single-cone geometry. All inputs are photon spectra.
"""
import argparse
from pathlib import Path
from time import perf_counter
import json
import os
import sys

import numpy as np

root = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(root / "src"), str(root / "xspec/windkin")]

from redwine.geometry import Geometry
from redwine.experimental.kernels import _moments
from windkin.xspec import wuconv as python_uniform, _COMMON, _SPREAD
from xspec import AllData, AllModels, Fit, Model, Spectrum, Xset


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spectrum", type=Path, required=True)
    parser.add_argument("--table", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--trials", type=int, default=12)
    parser.add_argument("--session", type=Path, help="Optional broadband 6abs1em XCM timing check")
    args = parser.parse_args()
    spectrum_path, table = args.spectrum.resolve(), args.table.resolve()
    output_path = args.output.resolve()
    session = args.session.resolve() if args.session else None
    assert args.trials >= 2
    Xset.allowPrompting = False
    Xset.chatter = 0
    Xset.parallel.reset()
    AllModels.lmod("redwine", str(root / "xspec"))
    AllModels.lmod("windkin", str(root / "xspec/windkin/build/native"))
    def pywuconv(edges, params, flux):
        python_uniform(edges, params, flux)
    AllModels.addPyMod(pywuconv, (_COMMON[0], _SPREAD["uniform"], *_COMMON[1:]), "con")
    os.chdir(spectrum_path.parent)
    data = Spectrum(spectrum_path.name)
    data.ignore("**-3.0 11.0-**")

    rng = np.random.default_rng(456)
    inclinations = rng.uniform(15, 65, args.trials)
    covering = rng.uniform(.45, .9, args.trials)
    skews = rng.uniform(.25, .75, args.trials)
    speeds = rng.uniform(.23, .35, args.trials)
    widths = rng.uniform(.04, .20, args.trials)
    geometries = Geometry.los(Cf=covering, skew=skews, incl=inclinations)
    rows = []
    folded = {}
    agreement = []

    for bins in (20000, 40000, 80000):
        AllModels.clear()
        AllModels.setEnergies(f"0.001 60 {bins} log")
        native_trials = []
        for name in ("wemlos_Cf", "wuconv", "pywuconv"):
            model = Model(f"{name}*atable{{{table}}}")
            if name in ("wuconv", "pywuconv"):
                values = [.28, .17, 15, 75, 45, 0, 5.5, 2500, 24, .184, 6.26739e-12]
            else:
                values = [.28, .7, .5, 45, 5.5, 2500, 24, .184, 6.26739e-12]
            model.setPars(*values)
            # Common reference point for the resolution comparison.
            if name in ("wuconv", "pywuconv"):
                geometry = Geometry.los(Cf=.7, skew=.5, incl=45)
                model.setPars({3: geometry.tin, 4: geometry.tout})
            folded[bins, name] = np.asarray(model.folded(1))
            elapsed = []
            for index in range(args.trials):
                if name in ("wuconv", "pywuconv"):
                    values = {1: speeds[index], 2: widths[index],
                              3: float(geometries.tin[index]), 4: float(geometries.tout[index]),
                              5: inclinations[index], 7: 5.4 + index*.01,
                              8: 2000 + index*25}
                else:
                    values = {1: speeds[index], 2: covering[index], 3: skews[index],
                              4: inclinations[index], 5: 5.4 + index*.01,
                              6: 2000 + index*25}
                _moments.cache_clear()
                start = perf_counter()
                model.setPars(values)
                prediction = np.asarray(model.folded(1))
                elapsed.append(perf_counter() - start)
                assert np.isfinite(prediction).all() and prediction.sum() > 0
                if name == "wuconv":
                    native_trials.append(prediction)
                elif name == "pywuconv":
                    error = float(abs(prediction-native_trials[index]).sum()/prediction.sum())
                    assert error < 1e-7, (bins, index, error)
                    agreement.append(dict(bins=bins, trial=index, native_python_relative_l1=error))
            row = dict(model=name, bins=bins, trials=args.trials,
                       median_seconds=float(np.median(elapsed)),
                       p90_seconds=float(np.quantile(elapsed, .9)))
            rows.append(row)
            print(json.dumps(row), flush=True)

    convergence = []
    for name in ("wemlos_Cf", "wuconv", "pywuconv"):
        fine = folded[80000, name]
        for bins in (20000, 40000):
            coarse = folded[bins, name]
            row = dict(model=name, bins=bins,
                       folded_relative_l1=float(abs(coarse-fine).sum()/fine.sum()),
                       expected_count_discrepancy=float(data.exposure*np.sum(
                           (coarse-fine)**2/np.maximum(fine, 1e-100))))
            convergence.append(row)
            print(json.dumps(row), flush=True)
    # Native RED-WINE uses center-bin quadrature; windkin integrates bins.
    # Compare their zero-width line profiles at a fine resolution.
    mono_limits = []
    for incl in (15., 45., 75.):
        geometry = Geometry.los(Cf=.7, skew=.5, incl=incl)
        model = Model("wlinelos_Cf")
        model.setPars(.28, .7, .5, incl, 6.7, 1.)
        reference = np.asarray(model.values(1))
        model = Model("wuline")
        model.setPars(.28, 0., geometry.tin, geometry.tout, incl, 0, 6.7, 1.)
        actual = np.asarray(model.values(1))
        error = float(abs(actual-reference).sum()/reference.sum())
        assert error < 3e-4, (incl, error)
        row = dict(inclination=incl, zero_width_relative_l1=error)
        mono_limits.append(row)
        print(json.dumps(row), flush=True)

    broadband = benchmark_broadband(session, args.trials, inclinations, speeds, widths, geometries) if session else None
    output_path.write_text(json.dumps(
        dict(xspec=list(Xset.version), spectrum=str(spectrum_path), table=str(table),
             seconds=rows, convergence=convergence, mono_limits=mono_limits,
             native_python_agreement=agreement, broadband=broadband),
        indent=2) + "\n")
    AllModels.clear()
    AllData.clear()



def benchmark_broadband(session, trials, inclinations, speeds, widths, geometries):
    """Time parameter updates and likelihood for the three-detector spectrum.

    Read the saved session without changing it. Copy its evaluated parameter
    values, replace only the emission convolution, then perturb all parameters
    that were independent and thawed in that session. This tests evaluation
    cost and native/reference agreement, not a refit or posterior sample.
    """
    AllData.clear()
    os.chdir(session.parent)
    Xset.restore(str(session))
    Xset.parallel.reset()
    snapshots = {}
    for group in range(1, AllData.nGroups + 1):
        model = AllModels(group)
        snapshots[group] = [model(i).values[0] for i in range(1, model.nParameters + 1)]
    base = AllModels(1)
    expression = base.expression
    assert "wemlos_Cf" in expression
    offset = base.wemlos_Cf.vout.index - base.startParIndex
    independent = [i for i in range(1, base.nParameters+1)
                   if not base(i).frozen and not base(i).link and not offset < i <= offset+4]
    reference_geometry = Geometry.los(Cf=snapshots[1][offset+1],
                                     skew=snapshots[1][offset+2], incl=snapshots[1][offset+3])
    outputs = {}
    rows = []
    for name in ("wuconv", "pywuconv"):
        AllModels.clear()
        AllModels.setEnergies("0.001 60 40000 log")
        Model(expression.replace("wemlos_Cf", name))
        models = {group: AllModels(group) for group in snapshots}
        values = {}
        for group, old in snapshots.items():
            parameters = old[:offset] + [old[offset], .17, reference_geometry.tin,
                                         reference_geometry.tout, old[offset+3], 0.] + old[offset+4:]
            values[group] = dict(enumerate(parameters, 1))
        AllModels.setPars(*[item for group, model in models.items() for item in (model, values[group])])
        elapsed, predictions = [], []
        for trial in range(trials):
            updates = {}
            for group in models:
                changes = {}
                for i in independent:
                    new_index = i if i <= offset else i+2
                    changes[new_index] = snapshots[group][i-1]*(1 + .0005*np.sin(trial+i))
                changes.update({offset+1: speeds[trial], offset+2: widths[trial],
                                offset+3: float(geometries.tin[trial]),
                                offset+4: float(geometries.tout[trial]),
                                offset+5: inclinations[trial]})
                updates[group] = changes
            _moments.cache_clear()
            start = perf_counter()
            AllModels.setPars(*[item for group, model in models.items() for item in (model, updates[group])])
            statistic = Fit.statistic
            prediction = np.concatenate([model.folded(group) for group, model in models.items()])
            elapsed.append(perf_counter()-start)
            assert np.isfinite(statistic) and np.isfinite(prediction).all()
            predictions.append((statistic, prediction))
        outputs[name] = predictions
        row = dict(model=name, median_seconds=float(np.median(elapsed)),
                   p90_seconds=float(np.quantile(elapsed, .9)), trials=trials)
        rows.append(row)
        print("Broadband " + json.dumps(row), flush=True)
    discrepancies, statistics = [], []
    for native, python in zip(outputs["wuconv"], outputs["pywuconv"]):
        discrepancies.append(float(abs(native[1]-python[1]).sum()/python[1].sum()))
        statistics.append(float(abs(native[0]-python[0])))
    assert max(discrepancies) < 1e-7, discrepancies
    assert max(statistics) < 1e-5, statistics
    result = dict(session=str(session), bins=40000, seconds=rows,
                  max_folded_relative_l1=max(discrepancies), max_statistic_difference=max(statistics))
    print("Broadband agreement " + json.dumps(result), flush=True)
    return result

if __name__ == "__main__":
    main()
