"""Engine-independent spectral grids for the Python wind models."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from .xspec import read_table_model


class SpectralGrid:
    """Interpolate rest-frame optical depths and bin-integrated emission.

    Axes are ``logxi``, ``vturb`` (km/s), ``log_nh`` (log10 cm-2), and
    energy-bin centres (keV). ``interpolation='log'`` interpolates positive
    spectra in log10 space, replacing exact zeros with ``floor``. ``'linear'``
    preserves zeros. Parameter coordinates are linear or logarithmic as
    specified by ``axis_methods`` (OGIP METHOD 0 or 1).
    """

    def __init__(
        self,
        energy: np.ndarray,
        logxi: np.ndarray,
        vturb: np.ndarray,
        log_nh: np.ndarray,
        *,
        optical_depth: Mapping[str, np.ndarray] | None = None,
        emission: Mapping[str, np.ndarray] | None = None,
        interpolation: str = "log",
        method: str = "auto",
        axis_methods: tuple[int, int, int] = (0, 0, 0),
        floor: float = 1e-30,
    ) -> None:
        if interpolation not in {"log", "linear"}:
            raise ValueError("interpolation must be 'log' or 'linear'")
        if not np.isfinite(floor) or floor <= 0:
            raise ValueError("floor must be finite and positive")

        self.energy = np.asarray(energy, dtype=float)
        self.logxi = np.asarray(logxi, dtype=float)
        self.vturb = np.asarray(vturb, dtype=float)
        self.log_nh = np.asarray(log_nh, dtype=float)
        axes = (self.logxi, self.vturb, self.log_nh)
        for name, axis in zip(("energy", "logxi", "vturb", "log_nh"), (self.energy, *axes)):
            if axis.ndim != 1 or len(axis) == 0 or not np.all(np.isfinite(axis)) or np.any(np.diff(axis) <= 0):
                raise ValueError(f"{name} must be a nonempty, finite, increasing axis")
        if len(self.energy) < 2 or np.any(self.energy <= 0):
            raise ValueError("energy must contain at least two positive centres")
        if len(axis_methods) != 3 or any(value not in (0, 1) for value in axis_methods):
            raise ValueError("axis_methods must contain three OGIP METHOD values")
        if method == "auto":
            method = "pchip" if all(len(axis) >= 4 for axis in axes) else "linear"
        if method not in {"linear", "pchip"} or (method == "pchip" and any(len(axis) < 4 for axis in axes)):
            raise ValueError("method must be 'linear' or 'pchip' with at least four points per axis")

        self.axis_methods = axis_methods
        self.interpolation = interpolation
        self.floor = floor
        self.min_spacing = float(np.min(np.diff(self.energy) / self.energy[:-1])) if len(self.energy) > 1 else 0.0
        if any(axis_method == 1 and np.any(axis <= 0) for axis, axis_method in zip(axes, axis_methods)):
            raise ValueError("logarithmic parameter axes must be positive")
        transformed_axes = tuple(
            np.log10(axis) if axis_method == 1 else axis
            for axis, axis_method in zip(axes, axis_methods)
        )

        shape = tuple(len(axis) for axis in axes) + (len(self.energy),)
        self._optical_depth = self._interpolators(optical_depth or {}, shape, transformed_axes, method)
        self._emission = self._interpolators(emission or {}, shape, transformed_axes, method)

    def _interpolators(self, spectra, shape, axes, method):
        result = {}
        for mode, values in spectra.items():
            values = np.asarray(values, dtype=float)
            if values.shape != shape or not np.all(np.isfinite(values)) or np.any(values < 0):
                raise ValueError(f"{mode!r} must be a finite, nonnegative array of shape {shape}")
            if self.interpolation == "log":
                values = np.log10(np.maximum(values, self.floor))
            result[mode] = RegularGridInterpolator(axes, values, method=method)
        return result

    def _evaluate(self, interpolators, mode, logxi, vturb, log_nh):
        if mode not in interpolators:
            raise ValueError(f"Unknown or unavailable spectral mode: {mode!r}")
        coordinates = np.asarray((logxi, vturb, log_nh), dtype=float)
        if not np.all(np.isfinite(coordinates)):
            raise ValueError("grid coordinates must be finite")
        if any(flag == 1 and value <= 0 for flag, value in zip(self.axis_methods, coordinates)):
            raise ValueError("logarithmic grid coordinates must be positive")
        coordinates = [np.log10(value) if flag == 1 else value for flag, value in zip(self.axis_methods, coordinates)]
        values = np.squeeze(interpolators[mode](coordinates))
        return 10.0**values if self.interpolation == "log" else values

    def optical_depth(self, logxi: float, vturb: float, log_nh: float, mode: str = "absorption") -> np.ndarray:
        """Return comoving optical depth sampled on ``energy``."""
        return self._evaluate(self._optical_depth, mode, logxi, vturb, log_nh)

    def emitted(self, logxi: float, vturb: float, log_nh: float, mode: str = "total") -> np.ndarray:
        """Return bin-integrated comoving emission sampled on ``energy``."""
        return self._evaluate(self._emission, mode, logxi, vturb, log_nh)

    @classmethod
    def from_tables(
        cls,
        *,
        optical_depth: Mapping[str, str | Path] | None = None,
        emission: Mapping[str, str | Path] | None = None,
        fixed_parameters: Mapping[str, float] | None = None,
        interpolation: str = "log",
        method: str = "auto",
    ) -> SpectralGrid:
        """Load matching OGIP tables without an ionization-engine dependency.

        Optical-depth tables must be exponential or multiplicative; emission
        tables must be additive. Supply omitted singleton axes through
        ``fixed_parameters``. All tables must share bins and parameter axes.
        """
        sources = [("optical_depth", mode, path) for mode, path in (optical_depth or {}).items()]
        sources += [("emission", mode, path) for mode, path in (emission or {}).items()]
        if not sources:
            raise ValueError("at least one spectral table is required")
        aliases = {"logxi": "logxi", "vturb": "vturb", "lognh": "log_nh", "log_nh": "log_nh"}
        fixed = {}
        for key, value in (fixed_parameters or {}).items():
            name = aliases.get(key.lower())
            if name is None or name in fixed or not np.isfinite(value):
                raise ValueError(f"Invalid fixed parameter {key!r}")
            fixed[name] = float(value)
        order = ("logxi", "vturb", "log_nh")
        reference = None
        arrays = {"optical_depth": {}, "emission": {}}

        for kind, mode, path in sources:
            table = read_table_model(path)
            if kind == "emission" and table.model_type != "additive":
                raise ValueError(f"{path} must be an additive table")
            if kind == "optical_depth" and table.model_type not in {"exponential", "multiplicative"}:
                raise ValueError(f"{path} must be an exponential or multiplicative table")
            try:
                names = [aliases[name.lower()] for name in table.parameter_names]
            except KeyError as error:
                raise ValueError(f"Unsupported table parameter: {error.args[0]}") from error
            if len(set(names)) != len(names) or any(name in fixed for name in names):
                raise ValueError("table parameters must be unique and separate from fixed_parameters")
            if any(name not in names and name not in fixed for name in order):
                raise ValueError("supply omitted axes in fixed_parameters")
            axes = tuple(table.parameter_axes[names.index(name)] if name in names else np.array([fixed[name]]) for name in order)
            methods = tuple(table.parameter_methods[names.index(name)] if name in names else 0 for name in order)
            signature = (table.energy_lo, table.energy_hi, axes, methods)
            if reference is None:
                reference = signature
            elif any(not np.array_equal(left, right) for left, right in zip(reference[:2], signature[:2])) or any(
                not np.array_equal(left, right) for left, right in zip(reference[2], axes)
            ) or reference[3] != methods:
                raise ValueError("tables must share energy bins, parameter axes, and interpolation methods")

            present = [name for name in order if name in names]
            values = np.transpose(table.spectra, [names.index(name) for name in present] + [len(names)])
            for index, name in enumerate(order):
                if name not in names:
                    values = np.expand_dims(values, index)
            if kind == "optical_depth" and table.model_type == "multiplicative":
                if np.any(values > 1) or np.any(values < 0):
                    raise ValueError("transmission must lie between zero and one")
                values = -np.log(np.maximum(values, 1e-30))
            arrays[kind][mode] = values

        lo, hi, axes, methods = reference
        return cls(
            0.5 * (lo + hi), *axes,
            optical_depth=arrays["optical_depth"], emission=arrays["emission"],
            interpolation=interpolation, method=method, axis_methods=methods,
        )
