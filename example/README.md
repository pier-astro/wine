# Two-zone example

The main folder contains the source `model.xcm`, its external input tables
`habs.mmod`/`hem.amod`, one simulated spectrum (`simulated.fak` and
`simulated.grp`), and `response.rmf`/`response.arf`. All analysis products belong in `results/`.

## Reproduce the spectrum

Run `aux/build_response.ipynb` from `aux/` to regenerate the responses in this
folder. It reads `aux/epic-pn-area.txt`: EPIC-pn area with an ideal Gaussian
redistribution of **50 eV FWHM** and 10 eV detector channels.

After initializing HEASoft, from this folder run:

```bash
bash aux/simulate.sh
```

Edit `exposure`, `seed`, and `min_counts` at the top of the Bash script.
It calls native XSPEC `fakeit none` with Poisson fluctuations, then
`ftgrouppha optmin` with a minimum of **10 counts**. Unlike fixed 20 eV
grouping, this is adaptive and may retain 10 eV channels. RMF and ARF are
associated by local filename; `BACKFILE=NONE`. Load the ready spectrum with
`data simulated.grp` from this directory. Re-running replaces only the two
simulated spectra and the generation log. The input tables are supplied
separately and are not redistributed through Git.

## Fit and save both scans

In XSPEC, working from this folder:

```tcl
lmod wine ../xspec
data simulated.grp
@model.xcm
energies 0.01 40 6000 log
ignore **-1.0 10.0-**
statistic cstat
method leven 50 0.001
fit
save all results/fit.xcm

source ../../xspeclab/procedures/export_steppar.tcl
steppar best 18 0 90 20 16 0 1 20
savesteppar 18 16 outfile results/incl-cf.fits clobber yes

source ../procedures/wind_geometry.tcl
windsteppar 7 tout 0 90 20 tin 0 90 20
savesteppar result $xspec_tclout outfile results/tout-tin.fits clobber yes
```

Adjust the xspeclab exporter path if installed elsewhere. **19 intervals give
20 samples per axis**, so each scan has 400 grid cells. The native scan profiles
`skew` and continuum normalization; the angular scan profiles inclination and
continuum normalization. Forbidden angular cells retain a mask/status.
Reapply the `energies` command when reloading an ordinary saved session:
native `save all` does not preserve that calculation grid.

## Plot the saved results

Open `results/geometry.ipynb` from `results/`. It reads both FITS files, draws
joint likelihood contours and curves of constant
`Cf = cos(tin) - cos(tout)`, and writes `geometry.png`/`geometry.pdf` there.
It uses NumPy, Matplotlib and xspeclab; no live XSPEC session is needed.
There is no smoothing across the forbidden triangle. Generated fit, scan,
plot and log files are ignored by Git; rerun the commands to reproduce them.
