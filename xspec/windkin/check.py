"""Run independent physics checks, numerical benchmarks and optional XSPEC checks."""

import argparse
from dataclasses import replace
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.integrate import quad
from scipy.special import ndtr, ndtri, roots_legendre
from scipy.stats import qmc

from windkin import Wind, convolve, line


def angular_reference(wind, edges, rest_energy, power=3):
    """Integrate beta, polar cosine and azimuth directly via scrambled Sobol."""
    speed, polar, azimuth = qmc.Sobol(3, scramble=True, seed=456).random_base2(21).T
    if wind.spread == 0:
        beta = wind.center
    elif wind.law == "uniform":
        beta = wind.center + wind.spread*(speed - .5)
    else:
        bounds = np.clip(wind.center + wind.spread*np.array([-3., 3.]), -1., 1.)
        lo, hi = ndtr((bounds-wind.center)/wind.spread)
        beta = wind.center + wind.spread*ndtri(lo + (hi-lo)*speed)
    cos_out, cos_in = np.cos(np.deg2rad([wind.outer, wind.inner]))
    cos_theta = cos_out + (cos_in - cos_out)*polar
    incl = np.deg2rad(wind.inclination)
    if wind.bicone:
        azimuth *= 2
        cos_theta = np.where(azimuth < 1, cos_theta, -cos_theta)
    mu = (cos_theta*np.cos(incl)
          + np.sqrt(1-cos_theta**2)*np.sin(incl)*np.cos(2*np.pi*azimuth))
    delta = np.sqrt((1-beta)*(1+beta))/(1 + (1 if wind.inflow else -1)*beta*mu)
    weights = delta**power * (wind.solid_angle_fraction/len(speed))
    return np.histogram(delta*rest_energy, edges, weights=weights)[0], delta, weights


def analytic_sphere(edges, rest_energy, beta):
    """Exact delta^3 photon weights for a mono-speed full sphere."""
    if beta == 0:
        return ((edges[:-1] <= rest_energy) & (rest_energy < edges[1:])).astype(float)
    gamma_inverse = np.sqrt((1-beta)*(1+beta))
    max_shift = np.sqrt((1+beta)/(1-beta))
    x = np.clip(edges/rest_energy, 1/max_shift, max_shift)
    return gamma_inverse*np.diff(x*x)/(4*beta)


def validate(output):
    irregular = np.array([0., 3.7, 4.5, 5.05, 5.6, 6.2, 6.7, 6.95,
                          7.35, 7.8, 8.15, 8.5, 9.2, 10., 17., 200.])
    rest = 6.7
    for beta in (1e-7, .25, .999999):
        wind = Wind("uniform", beta, 0, 0, 180, 90)
        grid = np.geomspace(rest*1e-4, rest*1e4, 1200)
        exact = analytic_sphere(grid, rest, beta)
        actual = line(grid, rest, wind)
        np.testing.assert_allclose(actual, exact, atol=2e-8, rtol=5e-7)
        np.testing.assert_allclose(actual.sum(), 1/np.sqrt(1-beta*beta), rtol=5e-7)

    for law, center, spread in (("uniform", .05, .10), ("gaussian", .02, .11),
                                ("gaussian", .98, .12), ("gaussian", .5, 1.)):
        wind = Wind(law, center, spread, 0, 180, 0)
        lo, hi = (center-spread/2, center+spread/2) if law == "uniform" else (
            max(-1., center-3*spread), min(1., center+3*spread))
        np.testing.assert_allclose(quad(wind.pdf, lo, hi)[0], 1., rtol=2e-12)
        grid = np.geomspace(.001, 500, 351)
        expected = quad(lambda beta: wind.pdf(beta)*analytic_sphere(
            grid[[0, -1]], rest, abs(beta))[0], lo, hi, epsabs=1e-10, limit=300)[0]
        actual = line(grid, rest, wind)
        fine = line(grid, rest, wind, tolerance=1e-9)
        np.testing.assert_allclose(actual.sum(), expected, rtol=2e-6)
        np.testing.assert_allclose(actual, fine, atol=3e-7, rtol=3e-6)
        print(f"Sphere {law} beta0={center} spread={spread}: "
              f"integral relative error {abs(actual.sum()/expected-1):.2g}")

    truncated = Wind("gaussian", .25, .03, 0, 180, 0)
    assert np.all(truncated.pdf([.15, .35]) == 0)
    assert np.all(truncated.pdf([.17, .33]) > 0)

    for law, width in (("uniform", .1), ("gaussian", .03)):
        wind = Wind(law, .25, width, 20, 50, 45)
        for trial in (wind, replace(wind, inflow=True), replace(wind, bicone=True)):
            ref, delta, weights = angular_reference(trial, irregular, rest)
            actual = line(irregular, rest, trial)
            error = np.abs(actual-ref).sum()/ref.sum()
            assert error < 9e-4, (trial, error)
            print(f"Direct 3D {trial.law} inflow={trial.inflow} bicone={trial.bicone}: "
                  f"relative L1 {error:.2g}, integral {actual.sum():.9f}")

        # Convolution reference uses exact overlap for each sampled Doppler shift.
        source = np.zeros(len(irregular)-1)
        source[4:7] = .15, .6, .25
        ref = np.zeros_like(source)
        for j in np.flatnonzero(source):
            a, b = irregular[j:j+2]
            for k, (low, high) in enumerate(zip(irregular[:-1], irregular[1:])):
                overlap = np.minimum(b, high/delta) - np.maximum(a, low/delta)
                ref[k] += source[j] * (weights @ np.maximum(0., overlap))/(b-a)
        actual = convolve(irregular, source, trial)
        np.testing.assert_allclose(actual, ref, atol=1e-6, rtol=1e-5)
        np.testing.assert_allclose(actual.sum(), line(irregular, rest, trial).sum(),
                                   rtol=2e-7)
        print(f"Convolution direct 3D {law}: max abs error {abs(actual-ref).max():.2g}")

        # Speed delta limit away from the static endpoint.
        fine = np.linspace(4, 10, 1901)
        mono = line(fine, rest, replace(wind, spread=0))
        narrow = line(fine, rest, replace(wind, spread=1e-8))
        assert abs(narrow-mono).sum()/mono.sum() < 1e-6

        for incl in (0, 45, 90):
            cone = replace(wind, inclination=incl)
            inflow = replace(cone, inflow=True)
            bicone = replace(cone, bicone=True)
            near, far = [line(irregular, rest, w) for w in (cone, inflow)]
            np.testing.assert_allclose(line(irregular, rest, bicone), near+far,
                                       atol=3e-8, rtol=1e-6)
            np.testing.assert_allclose(line(irregular, rest, bicone),
                                       line(irregular, rest, replace(bicone, inflow=True)),
                                       atol=2e-10, rtol=1e-9)
            if incl == 90:
                np.testing.assert_allclose(near, far, atol=2e-10, rtol=1e-9)
            else:
                assert abs(near-far).sum() > .01

        grid = np.unique(np.r_[np.geomspace(.5, 20, 701), rest + np.array([-1e-6, 1e-6])])
        source = np.zeros(len(grid)-1)
        source[np.searchsorted(grid, rest)-1] = 1.
        np.testing.assert_allclose(convolve(grid, source, wind), line(grid, rest, wind),
                                   atol=2e-9, rtol=1e-5)
        # Removing an output edge discards shifted photons, with all source
        # photons still present on the shorter grid.
        short = grid[:np.searchsorted(grid, 7.6)+1]
        clipped = convolve(short, source[:len(short)-1], wind)
        assert clipped.sum() < .98*convolve(grid, source, wind).sum()
        np.testing.assert_allclose(clipped, convolve(grid, source, wind)[:len(clipped)],
                                   atol=2e-9, rtol=1e-5)

        geom = np.geomspace(.5, 20, 701)
        almost_geom = geom.copy()
        almost_geom[1:-1] *= 1 + 1e-10*np.sin(np.arange(699))
        source = np.sin(np.arange(700))**2
        np.testing.assert_allclose(convolve(geom, source, wind),
                                   convolve(almost_geom, source, wind), atol=1e-8, rtol=2e-7)

    # A mixed population must agree with direct signed-velocity quadrature,
    # be continuous at zero, and be invariant under reversal for a bicone.
    mixed = Wind("gaussian", .02, .08, 20, 50, 20)
    ref, delta, weights = angular_reference(mixed, irregular, rest)
    actual = line(irregular, rest, mixed)
    assert np.abs(actual-ref).sum()/ref.sum() < 5e-4
    source = np.zeros(len(irregular)-1)
    source[6] = 1.
    a, b = irregular[6:8]
    ref_conv = [(weights @ np.maximum(0., np.minimum(b, hi/delta)
                                     - np.maximum(a, lo/delta)))/(b-a)
                for lo, hi in zip(irregular[:-1], irregular[1:])]
    np.testing.assert_allclose(convolve(irregular, source, mixed), ref_conv, atol=5e-7)
    near_zero = [line(irregular, rest, replace(mixed, center=center))
                 for center in (-1e-7, 0., 1e-7)]
    for profile in (near_zero[0], near_zero[-1]):
        np.testing.assert_allclose(profile, near_zero[1], atol=2e-7, rtol=3e-6)
    bicone = replace(mixed, bicone=True)
    np.testing.assert_allclose(line(irregular, rest, bicone),
                               line(irregular, rest, replace(bicone, center=-bicone.center)),
                               atol=2e-9)
    for center in (-.25, .25):
        mono = replace(mixed, center=center, spread=0.)
        narrow = replace(mono, spread=1e-8)
        np.testing.assert_allclose(line(irregular, rest, mono),
                                   line(irregular, rest, narrow), atol=2e-8)
    print("Signed Gaussian: direct line and convolution checks, zero crossing and mono limits passed.")

    static = Wind("gaussian", 0, 0, 20, 50, 45, bicone=True)
    exact = np.zeros(len(irregular)-1)
    exact[np.searchsorted(irregular, rest, side="right")-1] = static.solid_angle_fraction
    np.testing.assert_array_equal(line(irregular, rest, static), exact)
    np.testing.assert_array_equal(convolve(irregular, exact, static),
                                  static.solid_angle_fraction*exact)

    if output:
        import matplotlib.pyplot as plt

        edges = np.linspace(4.5, 10, 1101)
        energy = (edges[:-1]+edges[1:])/2
        fig, axes = plt.subplots(2, 3, figsize=(11, 6), sharex=True, sharey=True,
                                 layout="constrained")
        for row, width in enumerate((.1, .2)):
            for col, incl in enumerate((0, 45, 90)):
                ax = axes[row, col]
                for law, spread, label in (("uniform", width, "B · uniform"),
                                           ("gaussian", width/np.sqrt(12), "C · Gaussian")):
                    wind = Wind(law, .25, spread, 20, 50, incl)
                    start = perf_counter()
                    values = line(edges, rest, wind)
                    print(f"Line {law} width={width} i={incl}: "
                          f"{1000*(perf_counter()-start):.1f} ms")
                    ax.plot(energy, values/np.diff(edges), label=label)
                ax.set_title(f"$i={incl}^\\circ$, $\\Delta\\beta={width:.1f}$")
                if row == 1:
                    ax.set_xlabel("Observed energy (keV)")
                if col == 0:
                    ax.set_ylabel("Photon kernel (keV⁻¹)")
        axes[0, 0].legend(frameon=False)
        fig.suptitle("Single hollow cone 20°–50°, $E_0=6.7$ keV, $\\beta_0=0.25$",
                     fontsize=13)
        fig.savefig(output, dpi=160)


def validate_xspec(redwine, spectrum, dense):
    from xspec import AllData, AllModels, Model, Spectrum, Xset

    from redwine.experimental.kernels import Sector
    from windkin.xspec import register

    Xset.allowPrompting = False
    Xset.chatter = 5
    if redwine:
        AllModels.lmod("redwine", str(redwine / "xspec"))
    register()
    print("PyXspec/XSPEC", Xset.version)
    if spectrum:
        data = Spectrum(str(spectrum))
        # Use the observation's RMF and ARF directly from its headers.
        group = data.dataGroup
    else:
        AllData.dummyrsp(.4, 20., 1200, "log")
        group = 0

    if redwine:
        for existing in ("wline", "wem*powerlaw"):
            assert np.sum(Model(existing).values(group)) > 0

    for law, name, conv in (("uniform", "wuline", "wuconv"),
                            ("gaussian", "wgline", "wgconv")):
        params = [.25, .10, 20., 50., 45., 1] if law == "uniform" else [
            .02, .08, 20., 50., 45., 0]
        wind = Wind(law, *params[:5], bicone=bool(params[5]))
        model = Model(name)
        model.setPars(*params, 6.7, 2.)
        edges = np.array(model.energies(group))
        expected = 2*line(edges, 6.7, wind)
        np.testing.assert_allclose(model.values(group), expected, rtol=2e-7, atol=5e-10)
        if spectrum:
            assert np.sum(model.folded(group)) > 0

        # Narrow input line probes real in-place XSPEC convolution without
        # a dense continuum calculation across the full 60,000-bin response.
        source = Model("gaussian")
        source.setPars(6.7, .0001, 1.)
        input_flux = np.array(source.values(group))
        start = perf_counter()
        model = Model(f"{conv}*gaussian")
        model.setPars(*params, 6.7, .0001, 1.)
        actual = model.values(group)
        elapsed = perf_counter()-start
        expected = convolve(edges, input_flux, wind)
        np.testing.assert_allclose(actual, expected, rtol=3e-7, atol=5e-10)
        print(f"XSPEC {name}/{conv}: {len(edges)-1} bins, {elapsed:.3f} s for "
              f"convolution setup and evaluation")
        if dense:
            source = Model("powerlaw")
            input_flux = np.array(source.values(group))
            start = perf_counter()
            expected = convolve(edges, input_flux, wind)
            elapsed = perf_counter()-start
            assert expected.sum() > 0 and np.isfinite(expected).all()
            print(f"Dense {conv}: {len(edges)-1} bins, {elapsed:.3f} s, "
                  f"output/input photon flux {expected.sum()/input_flux.sum():.6f}")
    sector_params = (-.25, 20., 50., 45., 355., 80., 1)
    sector = Sector(*sector_params)
    for name, z in (("wsectline", 0.), ("zwsectline", .184)):
        model = Model(name)
        extra = (z,) if name == "zwsectline" else ()
        model.setPars(*sector_params, 6.7, *extra, 2.)
        edges = np.asarray(model.energies(group))
        expected = 2*line(edges*(1+z), 6.7, sector)/(1+z)
        np.testing.assert_allclose(model.values(group), expected, rtol=2e-7, atol=5e-10)

    source = Model("gaussian")
    source.setPars(6.7, .0001, 1.)
    input_flux = np.asarray(source.values(group))
    model = Model("wsect*gaussian")
    model.setPars(*sector_params, 6.7, .0001, 1.)
    edges = np.asarray(model.energies(group))
    np.testing.assert_allclose(model.values(group), convolve(edges, input_flux, sector, resample=True),
                               rtol=3e-7, atol=5e-10)
    print("XSPEC wsect/wsectline/zwsectline agree with the Python sector kernel.")

    if spectrum:
        # Exercise the fast geometric-grid path with the real response
        # still loaded, including XSPEC's rebinning before response folding.
        AllModels.clear()
        folded_by_grid = {}
        for bins in (20000, 40000):
            AllModels.setEnergies(f"0.1 40.0 {bins} log")
            for conv in ("wuconv", "wgconv"):
                start = perf_counter()
                model = Model(f"{conv}*powerlaw")
                folded = np.asarray(model.folded(group))
                assert np.isfinite(folded).all() and folded.sum() > 0
                folded_by_grid[bins, conv] = folded
                print(f"XSPEC {conv} on extended {bins}-bin log grid: "
                      f"{perf_counter()-start:.3f} s including response folding")
            AllModels.clear()
        for conv in ("wuconv", "wgconv"):
            coarse, fine = [folded_by_grid[bins, conv] for bins in (20000, 40000)]
            relative_l1 = np.abs(coarse-fine).sum()/fine.sum()
            assert relative_l1 < 2e-5
            print(f"XSPEC {conv} folded spectrum, 20000 vs 40000 bins: "
                  f"relative L1 {relative_l1:.2g}")
        AllModels.setEnergies("reset")
    AllModels.clear()
    AllData.clear()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plot", type=Path)
    parser.add_argument("--xspec", action="store_true")
    parser.add_argument("--redwine", type=Path)
    parser.add_argument("--spectrum", type=Path)
    parser.add_argument("--dense", action="store_true")
    args = parser.parse_args()
    if args.xspec:
        validate_xspec(args.redwine, args.spectrum, args.dense)
    else:
        validate(args.plot)
    print("Validation passed.")
