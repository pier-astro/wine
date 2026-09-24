# Repository structure and development

The repository has two independent user interfaces:

- `src/wine/`: importable Python package with `kernels.py` (analytical profiles and convolution), `models.py` (observer-frame slab functions), `grid.py` (engine-independent spectral interpolation), and `xspec/table.py` (OGIP FITS reader).
- `xwine/`: C++ local XSPEC functions registered in `lmodel.dat`. Python and XSPEC can read the same precomputed spectra, but their default interpolation of spectral values differs; see the [Python guide](python.md#interpolation).

`tests/` checks Python grids, slab evaluation, and OGIP table reading. `xwine/test_geometry_models.py` checks that the constrained geometry wrappers agree with the base XSPEC models. Run the Python tests from the repository root:

```bash
python -m unittest discover -s tests -q
```

After an XSPEC build, run `python xwine/test_geometry_models.py` with a Python interpreter that provides PyXSPEC. The local build and cleanup commands are documented in the [XSPEC guide](xspec.md).

## Spectral boundary

WINE receives rest-frame spectra as arrays or OGIP tables and does not import or execute CLOUDY or XSTAR. Spectral production workflows belong in a separate tool. This separation lets a user supply custom spectra while keeping wind kinematics and geometry in WINE. Check units, energy bins, parameter definitions, and any radiative-transfer assumptions before combining files from different producers.

## Provenance and status

The initial Python kernels, slab models, OGIP reader, and XSPEC implementation were extracted from `pier-astro/xflow` at commit `e55c58b`. The older, divergent `WINE/xwine` tree was not copied into this repository. The first independent WINE 2.0 commit is `9e0d667`. Source code is released under [GPL-3.0-or-later](../LICENSE).

The physical derivation, validation range, and references will be expanded with the accompanying publication. Current software tests check implementation behaviour; they do not by themselves establish accuracy for every astrophysical source or spectral grid.
