from dataclasses import replace

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal
from scipy.special import roots_legendre
from scipy.stats import qmc

from redwine.experimental.kernels import Sector, Wind, convolve, line


def doppler(sector, cos_theta, phi):
    incl = np.deg2rad(sector.incl)
    mu = (cos_theta*np.cos(incl) + np.sqrt(1-cos_theta**2)*np.sin(incl)*np.cos(phi))
    beta = sector.vout
    signs = (1, -1) if sector.bicone else (1,)
    return [np.sqrt(1-beta**2)/(1-sign*beta*mu) for sign in signs]


def angular_points(sector, u, v):
    lo, hi = np.cos(np.deg2rad([sector.tout, sector.tin]))
    return lo+(hi-lo)*u, np.deg2rad(sector.phi+(v-.5)*sector.dphi)


@pytest.mark.parametrize("beta", [-.4, 0., .25])
@pytest.mark.parametrize("bicone", [False, True])
def test_full_sector(beta, bicone):
    edges = np.geomspace(.4, 3, 301)
    sector = Sector(beta, 15, 65, 37, 79, 360, bicone)
    wind = Wind("gaussian", beta, 0, 15, 65, 37, bicone=bicone)
    assert_array_equal(line(edges, 1., sector), line(edges, 1., wind))
    flux = np.exp(-((edges[:-1]-1)/.1)**2)
    assert_array_equal(convolve(edges, flux, sector), convolve(edges, flux, wind))


@pytest.mark.parametrize("bicone", [False, True])
@pytest.mark.parametrize("beta", [0., 1e-12, -1e-12])
def test_static_half_sector(beta, bicone):
    edges = [0., .9, 1.1, 2.]
    full = Sector(beta, 15, 65, 37, 359, 360, bicone)
    half = replace(full, dphi=180)
    expected = (np.cos(np.deg2rad(15))-np.cos(np.deg2rad(65)))*(1+bicone)/4
    assert_allclose(line(edges, 1., half), [0, expected, 0], rtol=1e-7)
    flux = np.array([0., 1., 0.])
    assert_allclose(convolve(edges, flux, half), flux*expected, atol=1e-8)


def test_periodicity_orientation_and_antipode():
    edges = np.linspace(.4, 2.3, 501)
    sector = Sector(.4, 20, 65, 60, 355, 40)
    forward = line(edges, 1., sector)
    assert_array_equal(forward, line(edges, 1., replace(sector, phi=715)))
    reverse = line(edges, 1., replace(sector, vout=-.4))
    both = line(edges, 1., replace(sector, bicone=True))
    assert_allclose(both, forward+reverse, atol=3e-9, rtol=1e-6)
    centers = (edges[:-1]+edges[1:])/2
    back = line(edges, 1., replace(sector, phi=175))
    assert (centers @ forward)/sum(forward) > 1.2
    assert (centers @ back)/sum(back) < 1.
    for incl in (0., 180.):
        aligned = replace(sector, incl=incl)
        expected = line(edges, 1., replace(aligned, dphi=360))*sector.dphi/360
        assert_allclose(line(edges, 1., aligned), expected, atol=1e-12, rtol=1e-7)
    assert not line(edges, 1., replace(sector, dphi=0)).any()
    assert not line(edges, 1., replace(sector, tin=sector.tout)).any()


@pytest.mark.parametrize("beta,tin,tout,incl,phi,dphi,bicone", [
    (.3, 15, 70, 63, 355, 80, False),
    (-.6, 23, 64, 21, 180, 180, False),
    (.7, 9, 80, 126, 0, 180, False),
    (-.5, 10, 65, 90, 105, 230, True),
    (.25, 20, 55, 62, 270, .001, True),
    (.5, 10, 124, 61, 0, .001, False),
    (.85, 12, 25, 112, 0, 359.999, True),
    (-.12, 0, 90, 90, 90, 180, False),
])
def test_profile_against_direct_angular_integration(beta, tin, tout, incl, phi, dphi, bicone):
    sector = Sector(beta, tin, tout, incl, phi, dphi, bicone)
    edges = np.geomspace(.25, 4., 51)
    actual = line(edges, 1., sector)
    # Directly integrate smooth total flux in wind-frame theta and phi.
    nodes, weights = roots_legendre(120)
    angles = angular_points(sector, (nodes[:, None]+1)/2, (nodes+1)/2)
    shifts = doppler(sector, *angles)
    area = (np.cos(np.deg2rad(tin))-np.cos(np.deg2rad(tout)))*dphi/720
    expected = sum(weights @ shift**2 @ weights for shift in shifts)*area/4
    assert_allclose(actual.sum(), expected, rtol=2e-6)
    # Integrate the energy-bin indicator in the original angular coordinates.
    u, v = qmc.Sobol(2, scramble=True, seed=37).random_base2(19).T
    shifts = doppler(sector, *angular_points(sector, u, v))
    binned = sum(np.histogram(shift, edges, weights=shift**2)[0] for shift in shifts)
    binned *= area/len(u)
    assert np.sum(abs(actual-binned)) < 1.8e-3*expected


@pytest.mark.parametrize("log_grid", [False, True])
def test_convolution_against_direct_bin_overlap(log_grid):
    edges = np.geomspace(.25, 4., 81) if log_grid else np.linspace(.25, 4., 81)
    sector = Sector(-.48, 20, 67, 63, 195, 55, True)
    flux = np.zeros(len(edges)-1)
    index = np.searchsorted(edges, 1.)-1
    flux[index] = 1.
    a, b = edges[index:index+2]
    u, v = qmc.Sobol(2, scramble=True, seed=5).random_base2(17).T
    expected = np.zeros_like(flux)
    for shift in doppler(sector, *angular_points(sector, u, v)):
        for k, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
            width = np.maximum(0., np.minimum(b, hi/shift)-np.maximum(a, lo/shift))
            expected[k] += np.mean(width*shift**2)/(b-a)
    expected *= sector.solid_angle_fraction/2
    actual = convolve(edges, flux, sector)
    assert np.sum(abs(actual-expected)) < 2e-4*sum(expected)
    assert_allclose(actual.sum(), expected.sum(), rtol=1e-6)


def test_resampled_sector_convolution_on_dense_linear_grid():
    edges = np.linspace(2., 12., 2049)
    centers = (edges[:-1] + edges[1:])/2
    flux = np.exp(-((centers - 6.7)/.02)**2) * np.diff(edges)
    sector = Sector(.258, 0, 90, 45, 182, 180)
    exact = convolve(edges, flux, sector)
    fast = convolve(edges, flux, sector, resample=True)
    assert np.sum(abs(fast - exact)) < 1e-4 * exact.sum()
    assert_allclose(fast.sum(), exact.sum(), rtol=1e-7)
