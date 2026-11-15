"""Validate the beta mixture and optionally benchmark it on log-energy grids."""
import argparse
import csv
from pathlib import Path
from time import perf_counter
import warnings

import numpy as np
from scipy.integrate import IntegrationWarning, quad
from scipy.special import roots_legendre

from redwine.experimental._uniform import Prepared, Profile
from redwine.kernels import lineshape


CASES = [
    ("aligned", .30, .20, (0, 20, 70)),
    ("inclined", .30, .20, (37, 20, 70)),
    ("limb", .40, .45, (50, 10, 40)),
    ("edge_on", .30, .20, (90, 0, 90)),
    ("fast", .80, .35, (35, 20, 70)),
    ("thin", .40, .20, (37, 39.9, 40)),
    ("zero_endpoint", .10, .20, (35, 20, 70)),
    ("empty", .30, .20, (37, 40, 40)),
]
METHODS = ("reference", "direct", "semi")


def prepare(beta, width, degrees, photon=True):
    return Prepared(Profile(6.97, beta, width, *np.deg2rad(degrees), photon=photon))


def forward_moment(p, power):
    """Independent integration in wind-axis coordinates, with no beta roots."""
    nodes, weights = roots_legendre(112)
    lo, hi = np.cos(p.outer), np.cos(p.inner)
    axis_mu = (lo+hi)/2 + (hi-lo)/2*nodes
    phi = np.arange(256)*2*np.pi/256
    mu = (np.cos(p.inclination)*axis_mu[:, None]
          + np.sin(p.inclination)*np.sqrt(1-axis_mu**2)[:, None]*np.cos(phi))
    beta = p.beta0 + p.width/2*nodes
    delta = np.sqrt(1-beta[:, None, None]**2)/(1-beta[:, None, None]*mu)
    return (hi-lo)/8 * np.einsum("b,a,ba->", weights, weights, (delta**power).mean(axis=2))


def energy_moment(kernel, moment=0):
    p = kernel.p
    # Energy breakpoints at beta endpoints and stationary delta(beta, mu).
    mu = kernel.tangencies
    beta = np.stack([np.full_like(mu, p.lower), np.full_like(mu, p.upper),
                     np.clip(mu, p.lower, p.upper)])
    knots = np.unique(np.r_[p.E0, (p.E0*np.sqrt(1-beta**2)/(1-beta*mu)).ravel()])
    return sum(quad(lambda e: float(kernel.evaluate(e))*(e/p.E0)**moment,
                    lo, hi, epsabs=2e-10, epsrel=2e-9)[0]
               for lo, hi in zip(knots[:-1], knots[1:]))


def checks():
    for photon in (False, True):
        for beta in (.01, .3, .975):
            for angles in ((0, 20, 70), (37, 20, 70), (90, 0, 90), (37, 40, 40)):
                k = prepare(beta, 0, angles, photon)
                energy = np.geomspace(0.5, 60, 1024)
                for method in METHODS:
                    np.testing.assert_array_equal(k.evaluate(energy, method),
                                                  lineshape(energy, 6.97, 1, beta, *np.deg2rad(angles), photon))
        for beta, width, angles in ((1e-4, 1e-4, (37, 20, 70)),
                                    (.3, .2, (0, 20, 70)), (.8, .35, (35, 20, 70)),
                                    (.1, .2, (35, 20, 70))):
            k = prepare(beta, width, angles, photon)
            for moment in (0, 1):
                power = (2 if photon else 3) + moment
                actual = energy_moment(k, moment)
                expected = forward_moment(k.p, power)
                np.testing.assert_allclose(actual, expected, rtol=2e-7, atol=2e-10)
                if beta == 1e-4:
                    np.testing.assert_allclose(actual, (k.cos_in-k.cos_out)/2, atol=1e-4, rtol=0)
                print(f"Moment {moment}, delta^{power}, beta={beta:g}: {actual:.11g} "
                      f"(angular {expected:.11g})", flush=True)

    for beta, width in ((.1, .3), (.9, .2), (.3, -1), (np.nan, .1), (.3, 1e-20)):
        try:
            prepare(beta, width, (0, 0, 90))
        except ValueError:
            pass
        else:
            raise AssertionError((beta, width))
    for method in METHODS:
        assert prepare(.1, .2, (35, 20, 70)).evaluate(6.97, method) == np.inf
        assert prepare(.1, .2, (0, 20, 70)).evaluate(6.97, method) == 0
        assert prepare(.3, .2, (37, 70, 20)).evaluate(6.97, method) == 0
        tangent_to_equator = prepare(.1, .2, (30, 0, 60))
        np.testing.assert_allclose(tangent_to_equator.evaluate(6.97, method), .1340521172617,
                                   rtol=2e-9)
        assert prepare(.1, .2, (30, 120, 150)).evaluate(6.97, method) == 0
        for angles in ((0, 0, 90), (37, 0, 90), (90, 0, 180)):
            k = prepare(.3, .2, angles)
            gmax = np.sqrt((1+k.p.upper)/(1-k.p.upper))
            np.testing.assert_array_equal(k.evaluate([6.97/gmax, 6.97*gmax], method), 0)
        # At a stationary Doppler shift the two beta roots coalesce.
        k = prepare(.5, .4, (0, 60, 90))
        tangent = k.p.E0/np.sqrt(1-k.cos_in**2)
        np.testing.assert_allclose(k.evaluate([tangent, tangent*(1+1e-10)], method), 0, atol=1e-14)
        assert k.evaluate(tangent*(1-1e-8), method) > 0

    # Legacy fixed-beta kernels give an independent pointwise check on smooth
    # portions of the mixture, where unsplit integration is safe.
    for photon in (True, False):
        for energy in (6.97, 7.5):
            k = prepare(.3, .02, (37, 20, 70), photon)
            p = k.p
            expected = quad(lambda beta: float(lineshape(
                energy, p.E0, p.norm, beta, p.inclination, p.inner, p.outer,
                is_photonflux=photon))/p.width, p.lower, p.upper)[0]
            for method in METHODS:
                np.testing.assert_allclose(k.evaluate(energy, method), expected, rtol=1e-11)
            narrow = prepare(.3, 1e-6, (37, 20, 70), photon)
            mono = prepare(.3, 0., (37, 20, 70), photon)
            np.testing.assert_allclose(narrow.evaluate(energy), mono.evaluate(energy), rtol=1e-9)

    # Deliberately coarse tolerances should approach the reference on refinement.
    k = prepare(.8, .35, (35, 20, 70), False)
    energy = np.geomspace(2., 20., 201)
    reference = k.evaluate(energy, "reference")
    for method in ("direct", "semi"):
        errors = [np.max(np.abs(k.evaluate(energy, method, tol)-reference))
                  for tol in (1e-4, 1e-6, 2e-9)]
        assert errors[2] < errors[0] and errors[2] < 1e-9, errors
        print(f"{method} refinement (max absolute) {errors}", flush=True)


def median_time(call, repeats):
    call()
    elapsed = []
    for _ in range(repeats):
        start = perf_counter()
        value = call()
        elapsed.append(perf_counter()-start)
    return 1e3*np.median(elapsed), value


def benchmark(output):
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["case", "nE", "method", "prep_ms", "eval_ms", "repeats",
                         "max_abs", "max_rel_above_floor", "max_abs_below_floor",
                         "relative_floor"])
        for name, beta, width, angles in CASES:
            for n in (512, 2048, 8192):
                prep, kernel = median_time(lambda: prepare(beta, width, angles), 101)
                energy = np.geomspace(.3, 25., n)
                elapsed, reference = median_time(lambda: kernel.evaluate(energy, "reference"), 3)
                assert np.min(reference) >= 0
                floor = max(1e-14, np.max(reference)*1e-6)
                above = reference > floor
                for method in METHODS:
                    repeats = 3 if method == "reference" else 9
                    ms, profile = ((elapsed, reference) if method == "reference" else
                                   median_time(lambda: kernel.evaluate(energy, method), repeats))
                    error = np.abs(profile-reference)
                    absolute = np.max(error)
                    relative = np.max(error[above]/reference[above], initial=0)
                    near_zero = np.max(error[~above], initial=0)
                    assert np.min(profile) >= 0
                    assert absolute < 1e-8 and relative < 2e-5 and near_zero < 1e-9
                    writer.writerow([name, n, method, prep, ms, repeats,
                                     absolute, relative, near_zero, floor])
                    file.flush()
                print(f"Benchmarked {name} at nE={n}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="Run benchmark and write results to CSV.")
    args = parser.parse_args()
    warnings.filterwarnings("error", category=IntegrationWarning)
    checks()
    if args.output:
        benchmark(args.output)
