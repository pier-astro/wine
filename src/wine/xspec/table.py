"""Numerical evaluation of OGIP XSPEC table models."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

import numpy as np


@dataclass
class TableSpectrum:
    """One table model evaluated on its native energy bins.

    ``values`` are the interpolated table values. For additive models these
    are bin-integrated photon fluxes; for exponential models they are optical
    depths; for multiplicative models they are transmissions.
    """

    energy_lo: np.ndarray
    energy_hi: np.ndarray
    values: np.ndarray
    model_type: str

    @property
    def energy(self) -> np.ndarray:
        """Energy-bin centres in keV."""
        return 0.5 * (self.energy_lo + self.energy_hi)

    @property
    def width(self) -> np.ndarray:
        """Energy-bin widths in keV."""
        return self.energy_hi - self.energy_lo

    @property
    def transmission(self) -> np.ndarray:
        """Transmission for multiplicative or exponential tables."""
        if self.model_type == "multiplicative":
            return self.values
        if self.model_type == "exponential":
            return np.exp(-self.values)
        raise ValueError("Additive table models do not define a transmission.")

    @property
    def photon_flux(self) -> np.ndarray:
        """Bin-integrated additive photon flux."""
        if self.model_type != "additive":
            raise ValueError("Only additive table models define a photon flux.")
        return self.values

    @property
    def photon_flux_density(self) -> np.ndarray:
        """Additive photon flux density in ph cm-2 s-1 keV-1."""
        return self.photon_flux / self.width


@dataclass
class XspecTableModel:
    """An OGIP XSPEC table model stored as a regular parameter grid.

    Parameters are interpolated exactly as XSPEC's table reader does: ``METHOD
    = 0`` uses a linear coordinate and ``METHOD = 1`` uses a logarithmic
    coordinate. In both cases the spectra are interpolated linearly.
    """

    energy_lo: np.ndarray
    energy_hi: np.ndarray
    parameter_names: tuple[str, ...]
    parameter_axes: tuple[np.ndarray, ...]
    parameter_methods: tuple[int, ...]
    spectra: np.ndarray
    model_type: str
    supports_redshift: bool

    def evaluate(
        self,
        parameters: Mapping[str, float],
        *,
        redshift: float = 0.0,
        normalization: float = 1.0,
    ) -> TableSpectrum:
        """Interpolate the table at one physical parameter point.

        Parameters
        ----------
        parameters
            Mapping of every table interpolation parameter to its value.
        redshift
            XSPEC table redshift parameter. It shifts the energy grid by
            ``1 / (1 + z)`` and applies time dilation to additive models.
        normalization
            XSPEC additive-model normalization. It must be one for
            multiplicative and exponential tables.
        """
        if not np.isfinite(redshift) or redshift <= -1.0:
            raise ValueError("redshift must be finite and greater than -1.")
        if not np.isfinite(normalization):
            raise ValueError("normalization must be finite.")
        if self.model_type != "additive" and normalization != 1.0:
            raise ValueError("normalization is only defined for additive table models.")

        values = self.spectra
        for axis_index in range(len(self.parameter_axes) - 1, -1, -1):
            name = self.parameter_names[axis_index]
            coordinate = _parameter_value(parameters, name)
            low, high, fraction = _bracket(
                self.parameter_axes[axis_index],
                coordinate,
                self.parameter_methods[axis_index],
                name,
            )
            if low == high:
                values = np.take(values, low, axis=axis_index)
            else:
                values = (
                    (1.0 - fraction) * np.take(values, low, axis=axis_index)
                    + fraction * np.take(values, high, axis=axis_index)
                )

        zfactor = 1.0 + redshift if self.supports_redshift else 1.0
        if self.model_type == "additive":
            values = values * normalization / zfactor
        return TableSpectrum(
            energy_lo=self.energy_lo / zfactor,
            energy_hi=self.energy_hi / zfactor,
            values=np.asarray(values, dtype=float),
            model_type=self.model_type,
        )


def read_table_model(path: str | Path) -> XspecTableModel:
    """Read an OGIP XSPEC additive, multiplicative, or exponential table.

    This reader supports regular grids without XSPEC filter expressions or
    additional parameter spectra. Those features raise explicitly rather than
    returning scientifically incomplete results.
    """
    try:
        from astropy.io import fits
    except ImportError as error:
        raise ImportError("Reading XSPEC tables requires astropy.") from error

    path = Path(path)
    try:
        hdul = fits.open(path, memmap=True)
    except OSError as error:
        raise ValueError(f"Not a readable XSPEC table model: {path}") from error

    with hdul:
        _require_extensions(hdul, path)
        primary = hdul[0].header
        parameters = hdul["PARAMETERS"].data
        spectra = hdul["SPECTRA"].data
        energy_lo = np.asarray(hdul["ENERGIES"].data["ENERG_LO"], dtype=float)
        energy_hi = np.asarray(hdul["ENERGIES"].data["ENERG_HI"], dtype=float)

        if int(hdul["PARAMETERS"].header.get("NADDPARM", 0)) != 0:
            raise NotImplementedError("XSPEC additional-parameter spectra are not supported.")
        if "PARAMVAL" not in spectra.names or "INTPSPEC" not in spectra.names:
            raise ValueError(f"SPECTRA must contain PARAMVAL and INTPSPEC: {path}")

        names = tuple(_as_text(row["NAME"]) for row in parameters)
        methods = tuple(int(row["METHOD"]) for row in parameters)
        axes = tuple(
            np.asarray(row["VALUE"], dtype=float)[: int(row["NUMBVALS"])]
            for row in parameters
        )
        if not names or any(method not in {0, 1} for method in methods):
            raise ValueError("Only interpolated XSPEC table parameters are supported.")
        if any(len(axis) == 0 or np.any(np.diff(axis) <= 0.0) for axis in axes):
            raise ValueError("XSPEC table parameter axes must be strictly increasing.")

        shape = tuple(len(axis) for axis in axes)
        values = np.asarray(spectra["INTPSPEC"], dtype=float)
        if values.shape != (int(np.prod(shape)), len(energy_lo)):
            raise ValueError("SPECTRA dimensions do not match the parameter and energy grids.")
        _validate_record_order(np.asarray(spectra["PARAMVAL"], dtype=float), axes)

        return XspecTableModel(
            energy_lo=energy_lo,
            energy_hi=energy_hi,
            parameter_names=names,
            parameter_axes=axes,
            parameter_methods=methods,
            spectra=values.reshape(*shape, len(energy_lo)),
            model_type=_model_type(primary),
            supports_redshift=bool(primary.get("REDSHIFT", False)),
        )


def _require_extensions(hdul, path: Path) -> None:
    required = {"PARAMETERS", "ENERGIES", "SPECTRA"}
    missing = required.difference(hdu.name for hdu in hdul)
    if missing:
        raise ValueError(f"Missing XSPEC table extensions {sorted(missing)}: {path}")


def _model_type(header) -> str:
    model_type = str(
        header.get("XFXP0001", "additive" if header.get("ADDMODEL") else "multiplicative")
    ).strip().lower()
    if model_type not in {"additive", "multiplicative", "exponential"}:
        raise ValueError(f"Unsupported XSPEC table type {model_type!r}.")
    return model_type


def _as_text(value) -> str:
    return value.decode().strip() if isinstance(value, bytes) else str(value).strip()


def _parameter_value(parameters: Mapping[str, float], name: str) -> float:
    lookup = {str(key).lower(): value for key, value in parameters.items()}
    if name.lower() not in lookup:
        raise KeyError(f"Missing table parameter {name!r}.")
    value = float(lookup[name.lower()])
    if not np.isfinite(value):
        raise ValueError(f"Table parameter {name!r} must be finite.")
    return value


def _bracket(axis: np.ndarray, value: float, method: int, name: str) -> tuple[int, int, float]:
    tolerance = 1.0e-6 * max(1.0, abs(value))
    if value < axis[0] - tolerance or value > axis[-1] + tolerance:
        raise ValueError(
            f"{name}={value:g} is outside the table range [{axis[0]:g}, {axis[-1]:g}]."
        )
    exact = np.flatnonzero(np.abs(axis - value) <= tolerance)
    if len(exact):
        return int(exact[0]), int(exact[0]), 0.0

    high = int(np.searchsorted(axis, value))
    low = high - 1
    if method == 0:
        fraction = (value - axis[low]) / (axis[high] - axis[low])
    else:
        if value <= 0.0 or axis[low] <= 0.0:
            raise ValueError(f"Logarithmic interpolation for {name!r} requires positive values.")
        fraction = np.log(value / axis[low]) / np.log(axis[high] / axis[low])
    return low, high, float(fraction)


def _validate_record_order(parameter_values: np.ndarray, axes: tuple[np.ndarray, ...]) -> None:
    expected = np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1).reshape(-1, len(axes))
    if len(axes) == 1 and parameter_values.ndim == 1:
        parameter_values = parameter_values[:, np.newaxis]
    if parameter_values.shape != expected.shape or not np.allclose(parameter_values, expected):
        raise ValueError("SPECTRA rows do not follow standard XSPEC parameter-grid ordering.")
