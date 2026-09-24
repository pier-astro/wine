# Repository structure and development

The repository has two independent user interfaces:

- `src/wine/`: importable Python package with `kernels.py` (analytical profiles and convolution), `models.py` (observer-frame slab functions), `grid.py` (engine-independent spectral interpolation), and `ogip.py` (OGIP FITS reader).
- `xspec/`: C++ local XSPEC functions registered in `lmodel.dat`. Python and XSPEC can read the same precomputed spectra, but their default interpolation of spectral values differs; see the [Python guide](python.md#interpolation).

`tests/` checks Python grids, slab evaluation, and OGIP table reading. `xspec/test_geometry_models.py` checks that the constrained geometry wrappers agree with the base XSPEC models. Run the Python tests from the repository root:

```bash
python -m unittest discover -s tests -q
```

After an XSPEC build, run `python xspec/test_geometry_models.py` with a Python interpreter that provides PyXSPEC. The local build and cleanup commands are documented in the [XSPEC guide](xspec.md).

## Spectral boundary

WINE receives rest-frame spectra as arrays or OGIP tables and does not import or execute CLOUDY or XSTAR. Spectral production workflows belong in a separate tool. This separation lets a user supply custom spectra while keeping wind kinematics and geometry in WINE. Check units, energy bins, parameter definitions, and any radiative-transfer assumptions before combining files from different producers.
