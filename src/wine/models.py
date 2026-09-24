from __future__ import annotations

import numpy as np
import warnings

from .grid import SpectralGrid
from .kernels import convolve

# Speed of light in km/s
C_KM_S = 299792.458


def _check_grid_resolution(E_obs: np.ndarray, grid: SpectralGrid) -> None:
    """Warn if the observer energy grid is denser than the intrinsic simulation resolution."""
    if len(E_obs) > 1:
        # We check the minimum fractional step size of the observer grid.
        min_spacing_obs = np.min(np.diff(E_obs) / E_obs[:-1])
        if min_spacing_obs < grid.min_spacing:
            warnings.warn(
                f"The requested energy grid is denser than the intrinsic spectral grid resolution. "
                f"Observer grid min fractional spacing: {min_spacing_obs:.2e}, "
                f"Simulation grid min fractional spacing: {grid.min_spacing:.2e}.",
                UserWarning,
                stacklevel=3,
            )


def rebin_flux(E_src: np.ndarray, F_src: np.ndarray, E_dst: np.ndarray) -> np.ndarray:
    """Rebin flux density from E_src to E_dst conserving total flux via CDF interpolation.

    Parameters
    ----------
    E_src : np.ndarray
        Source energy grid centers (must be sorted).
    F_src : np.ndarray
        Source flux density (e.g. photons/s/keV).
    E_dst : np.ndarray
        Destination energy grid centers (must be sorted).

    Returns
    -------
    F_dst : np.ndarray
        Rebinned flux density on E_dst.
    """
    if len(E_dst) == 1:
        return np.interp(E_dst, E_src, F_src, left=0.0, right=0.0)

    # Source bin widths and edges
    dE_src = np.diff(E_src)
    dE_src = np.append(dE_src, dE_src[-1])
    E_src_edges = np.zeros(len(E_src) + 1)
    E_src_edges[0] = E_src[0] - dE_src[0] / 2.0
    E_src_edges[1:] = E_src + dE_src / 2.0

    # Cumulative sum (CDF) of integrated flux
    cdf_src = np.zeros(len(E_src) + 1)
    cdf_src[1:] = np.cumsum(F_src * dE_src)

    # Destination bin widths and edges
    dE_dst = np.diff(E_dst)
    dE_dst = np.append(dE_dst, dE_dst[-1])
    E_dst_edges = np.zeros(len(E_dst) + 1)
    E_dst_edges[0] = E_dst[0] - dE_dst[0] / 2.0
    E_dst_edges[1:] = E_dst + dE_dst / 2.0

    # Interpolate CDF to destination edges
    cdf_dst = np.interp(E_dst_edges, E_src_edges, cdf_src, left=0.0, right=cdf_src[-1])

    # Difference to get bin integrated flux, then divide by width to get density
    F_dst = np.diff(cdf_dst) / dE_dst
    return F_dst


def absorption_slab(
    E_obs: np.ndarray,
    grid: SpectralGrid,
    logxi: float,
    log_nh: float,
    vturb: float,
    beta: float,
    C_f: float,
    z: float = 0.0,
    mode: str = "absorption",
) -> np.ndarray:
    """Compute the transmission spectrum for a single relativistic absorbing slab.

    Parameters
    ----------
    E_obs : np.ndarray
        Observer-frame energy grid in keV.
    grid : SpectralGrid
        Engine-independent interpolation grid of rest-frame spectra.
    logxi : float
        Gas-frame ionization parameter.
    log_nh : float
        Slab column density log10(N_H / cm^-2).
    vturb : float
        Microturbulent velocity in km/s.
    beta : float
        Relativistic bulk line-of-sight velocity of the slab (v/c).
        Positive beta corresponds to motion towards the observer (blueshift).
    C_f : float
        Geometric partial covering fraction of the absorber (between 0.0 and 1.0).
    z : float, optional
        Cosmological redshift of the host galaxy. Default is 0.0.
    mode : {'absorption', 'scattering', 'total'}, optional
        The optical depth component to use. Default is 'absorption'.

    Returns
    -------
    transmission : np.ndarray
        Observer-frame transmitted flux fraction (same shape as E_obs).
    """
    # 1. Validate velocity, covering fraction, and quarantine grid query coordinates
    # _check_grid_resolution(E_obs, grid)
    if not isinstance(logxi, (int, float, np.number)):
        raise TypeError(f"logxi must be a scalar number, got {type(logxi)}")
    if not isinstance(log_nh, (int, float, np.number)):
        raise TypeError(f"log_nh must be a scalar number, got {type(log_nh)}")
    if not isinstance(vturb, (int, float, np.number)):
        raise TypeError(f"vturb must be a scalar number, got {type(vturb)}")
    if not -1.0 < beta < 1.0:
        raise ValueError("Relativistic velocity beta must be in the open interval (-1, 1).")
    if not 0.0 <= C_f <= 1.0:
        raise ValueError("Covering fraction C_f must be between 0.0 and 1.0.")

    # 2. Compute Doppler factor: D = [gamma * (1 - beta * cos_theta_los)]^-1
    # For line-of-sight absorption, cos_theta_los = 1.0
    gamma = 1.0 / np.sqrt(1.0 - beta**2)
    D = 1.0 / (gamma * (1.0 - beta))

    # 3. Map observer energies to the gas rest frame: E0 = E_obs * (1 + z) / D
    E0 = E_obs * (1.0 + z) / D

    # 4. Clamp grid parameters to the pre-computed grid boundaries
    logxi_clamped = np.clip(logxi, grid.logxi[0], grid.logxi[-1])
    vturb_clamped = np.clip(vturb, grid.vturb[0], grid.vturb[-1])
    log_nh_clamped = np.clip(log_nh, grid.log_nh[0], grid.log_nh[-1])

    # 5. Interpolate rest-frame comoving optical depth from the spectral grid
    tau_grid = grid.optical_depth(logxi_clamped, vturb_clamped, log_nh_clamped, mode=mode)

    if len(E0) < 10:
        # Fallback to point-sampling for sparse energy arrays (e.g. unit tests or single line queries)
        tau_obs = np.interp(E0, grid.energy, tau_grid, left=0.0, right=0.0)
        T_obs = np.exp(-tau_obs)
    else:
        # Use Cumulative Integral Rebinning on absorbed fraction to prevent line aliasing / jittering
        A_grid = 1.0 - np.exp(-tau_grid)
        A_obs = rebin_flux(grid.energy, A_grid, E0)
        T_obs = 1.0 - A_obs

    # 8. Apply geometric partial covering fraction: T(E) = C_f * T_obs + (1 - C_f)
    return C_f * T_obs + (1.0 - C_f)


def emission_slab(
    E_obs: np.ndarray,
    grid: SpectralGrid,
    logxi: float,
    log_nh: float,
    vturb: float,
    beta: float,
    incl: float,
    theta_in: float,
    theta_out: float,
    norm: float,
    z: float = 0.0,
    method: str = "fft",
    mode: str = "total",
) -> np.ndarray:
    """Compute the relativistic emission spectrum from an outflowing shell/slab.

    Parameters
    ----------
    E_obs : np.ndarray
        Observer-frame energy grid in keV.
    grid : SpectralGrid
        Engine-independent interpolation grid of rest-frame spectra.
    logxi : float
        Gas-frame ionization parameter.
    log_nh : float
        Slab column density log10(N_H / cm^-2).
    vturb : float
        Microturbulent velocity in km/s.
    beta : float
        Relativistic bulk outflow velocity of the wind (v/c).
    incl : float
        Inclination of the observer relative to the wind axis in radians.
    theta_in : float
        Inner opening angle of the emitting wind in radians.
    theta_out : float
        Outer opening angle of the emitting wind in radians.
    norm : float
        Global normalization scaling factor (e.g., area / 4pi d_L^2).
    z : float, optional
        Cosmological redshift of the host galaxy. Default is 0.0.
    method : str, optional
        Convolution method: 'fft' (default, fast), 'skip', or 'reference'.
    mode : {'total', 'lines', 'continuum'}, optional
        The emission component to use. Default is 'total'.

    Returns
    -------
    flux : np.ndarray
        Observer-frame specific photon flux in photons/s/cm^2/keV (same shape as E_obs).
    """
    # 1. Validate velocity, opening angles, and quarantine grid query coordinates
    # _check_grid_resolution(E_obs, grid)
    if not isinstance(logxi, (int, float, np.number)):
        raise TypeError(f"logxi must be a scalar number, got {type(logxi)}")
    if not isinstance(log_nh, (int, float, np.number)):
        raise TypeError(f"log_nh must be a scalar number, got {type(log_nh)}")
    if not isinstance(vturb, (int, float, np.number)):
        raise TypeError(f"vturb must be a scalar number, got {type(vturb)}")
    if not 0.0 <= beta < 1.0:
        raise ValueError("Relativistic outflow velocity beta must be in the half-open interval [0, 1).")
    if not 0.0 <= theta_in < theta_out <= np.pi:
        raise ValueError("Wind opening angles must satisfy 0 <= theta_in < theta_out <= pi.")

    # 2. Clamp grid parameters to the pre-computed grid boundaries
    logxi_clamped = np.clip(logxi, grid.logxi[0], grid.logxi[-1])
    vturb_clamped = np.clip(vturb, grid.vturb[0], grid.vturb[-1])
    log_nh_clamped = np.clip(log_nh, grid.log_nh[0], grid.log_nh[-1])

    # 3. Query rest-frame emergent spectral luminosity from the grid
    # (now stored in bin-integrated photons/s/cm^2 at 10 kpc)
    S_nu_0 = grid.emitted(logxi_clamped, vturb_clamped, log_nh_clamped, mode=mode)

    # 4. Convert bin-integrated photon flux to specific photon flux density (photons/s/cm^2/keV)
    # by dividing by comoving grid bin widths.
    dE_grid = np.diff(grid.energy)
    dE_grid = np.append(dE_grid, dE_grid[-1])
    grid_edges = np.zeros(len(grid.energy) + 1)
    grid_edges[0] = grid.energy[0] - dE_grid[0] / 2.0
    grid_edges[1:] = grid.energy + dE_grid / 2.0
    grid_edges = np.maximum(grid_edges, 0.0)
    grid_bin_widths = np.diff(grid_edges)

    S_phot_0 = S_nu_0 / grid_bin_widths

    # 5. Shift observer grid to host galaxy frame: E_z = E_obs * (1 + z)
    E_z = E_obs * (1.0 + z)

    # 6. Build a padded energy grid to mitigate boundary/wrap-around effects in convolution.
    # The padding factor is scaled by beta to fully cover the Doppler shift range.
    dE_z = np.diff(E_z)
    dE_z = np.append(dE_z, dE_z[-1])
    E_z_min_edge = E_z[0] - dE_z[0] / 2.0
    E_z_max_edge = E_z[-1] + dE_z[-1] / 2.0

    pad_factor = 1.0 + 3.0 * beta
    # Add a small buffer to avoid boundary cutting at the edges
    E_min_padded = E_z_min_edge / (pad_factor * 1.1)
    E_max_padded = E_z_max_edge * (pad_factor * 1.1)

    n_points = max(len(E_obs), 1000)
    E_padded = np.geomspace(E_min_padded, E_max_padded, n_points)

    # 7. Rebin the rest-frame photon spectrum onto the padded grid using Cumulative Integral Rebinning
    S_phot_padded = rebin_flux(grid.energy, S_phot_0, E_padded)

    # 8. Relativistic Convolution in the host galaxy frame
    S_conv_padded = convolve(
        E_padded,
        S_phot_padded,
        beta,
        incl,
        theta_in,
        theta_out,
        is_photonflux=True,
        method=method,
    )

    # 9. Rebin back to host frame energies E_z conserving flux, and apply global normalization
    S_conv = rebin_flux(E_padded, S_conv_padded, E_z)

    return S_conv * norm
