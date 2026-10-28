# WINE: validation and paper revision checklist

Working decisions agreed on 2026-10-27. This is a checklist for targeted model
validation and subsequent manuscript revision, not a record of completed
scientific validation. The author will personally run and supervise the new
PDS 456 fits. No manuscript changes or broad parameter scans are part of this
planning stage.

## Agreed manuscript changes

### Document the geometric normalization check

- [ ] Include the solid-angle identity for the azimuthal weight, with
  `mu = cos(theta)` measured in the observer frame:

  \[
  \frac{1}{4\pi}\int_{-1}^{1} w_\phi(\mu)\,d\mu
  =\frac{\Omega}{4\pi}
  =\frac{\cos\theta_{\rm in}-\cos\theta_{\rm out}}{2}.
  \]

  For an input spectrum normalized to a full sphere, this gives the integrated
  flux fraction in the nonrelativistic limit: 1 for a full sphere and 1/2 for
  a single hemisphere. With the model's hemispheric definition
  `Cf = Omega/(2*pi)`, the fraction is `Cf/2`.
- [ ] Report the independent solid-angle reference and numerical integration
  of the implemented weight in
  [the geometric check](../tests/check_geometry_normalization.py): five
  geometries at three inclinations agree within absolute tolerance `1e-9`.
  This corroborates geometric normalization and its inclination independence;
  it does not validate finite-velocity boosting or the complete convolution.

### State explicitly how absorption and emission are combined

- [ ] Add the observer-frame model expression:

  \[
  F_{\rm obs}(E)=T_{\rm Gal}(E)
  \left[T_{\rm LOS}(E)F_{\rm cont}(E)
  +\sum_k F_{{\rm em},k}(E)\right].
  \]

- [ ] Define the continuum, the net intrinsic line-of-sight transmission, and
  the broadened emission components. Galactic absorption applies to both
  continuum and emission. Expand the net transmission into the actual
  full/partial-covering factors only where needed.
- [ ] Explain that the current prescription omits additional attenuation of
  wind emission along its escape paths through the wind. Verify and describe
  any transfer already included in the input photoionization emission tables;
  do not equate absence of an external XSPEC multiplier with absence of all
  internal radiative transfer.
- [ ] Motivate this approximation as a balance between physical information
  available, interpretability, and complexity. Multiplying every emitter by
  the same absorption tables as the central continuum assumes a common
  foreground screen that the observations do not establish.
- [ ] Do not claim that a wide opening angle alone guarantees negligible
  self-absorption, or that full-screen and unattenuated fits necessarily bracket
  all possible inferred parameters after refitting.

### Present the two PDS 456 tests as complementary

- [ ] Retain the existing hypothesis in which each selected absorber has an
  associated emitter, with the currently adopted physical links stated
  explicitly. Six independently significant emission detections are not a
  prerequisite for testing this linked physical prescription.
- [ ] Add a second test retaining only emission components supported by the
  data under a declared model-selection criterion. For each retained emitter,
  free its column density and ionization independently of its absorption
  counterpart AND independently of the other emitters.
- [ ] Tie the emission geometry across retained emitters to limit the number
  of free parameters. Sharing `uincl` alone does not share the physical viewing
  angle unless the opening angles are also shared.
- [ ] Specify the remaining velocity and normalization links before fitting;
  their treatment has not been decided by this discussion. Check normalization
  versus angular-coverage degeneracy if emission normalizations are free.
- [ ] Count the actual independent parameters. Fewer emitting components do
  not necessarily mean fewer free parameters: the second test releases physical
  links and can have more parameters than the associated-emission prescription.
- [ ] Explain that both models make assumptions and probe different questions:
  linked absorption/emission consistency versus the emission phases most
  representative of the observed spectrum. The second test does not uniquely
  identify every physical emitting zone.
- [ ] Compare geometry, plasma properties, and the main scientific conclusions,
  as well as fit statistic and model complexity. Assess sensitivity to plausible
  component-selection choices without starting with an exhaustive subset scan.
- [ ] Frame agreement as robustness against the tested modelling assumptions,
  not proof that all systematic biases have been removed.

### Reserve geometry-dependent self-absorption for future development

- [ ] Identify this as a future extension, not a correction already implemented.
- [ ] Distinguish attenuation of emission within its own emitting material from
  attenuation by other wind phases, and check what the input tables include.
- [ ] Note that an associated transmission table supplies opacity information,
  but not the emitting locations, remaining columns along escape paths, ordering
  of phases, or relative projected velocities encountered by emitted photons.
- [ ] Explain that an analytic escape/attenuation treatment is possible under
  additional spatial and emissivity assumptions. A general treatment must apply
  the appropriate energy transformation and path-dependent attenuation within
  the angular/volume integration, rather than multiplying the final broad
  emission profile by the central sightline transmission.
- [ ] Keep the present analytic geometry. Precomputing the full effect would
  require specifying that geometry during table production or tabulating the
  additional dependencies. An on-the-fly implementation is possible in
  principle, but needs access to both opacity and emissivity and a new interface
  or combined model; the current XSPEC convolution receives only its input
  spectrum and cannot obtain a separate associated table automatically.
- [ ] Keep any effective emission covering fraction distinct from global angular
  coverage and volume filling. Do not introduce a free covering parameter as a
  substitute for an unspecified transfer model.
- [ ] Cite the distinction from WINE 2024: its emission tables already included
  an approximate self-absorption treatment, although the XSPEC expression adds
  emission outside the continuum absorption term.

  Reference: [Luminari et al. 2024, Sections 2.4 and 5](https://art.torvergata.it/retrieve/46ccaa09-4f8e-45f7-ae40-697c3ea98bcf/aa49244-24.pdf).

## Implementation priorities before revising results

### 1. Keep energy coverage explicit and the operators simple

- [x] Verify that the documentation already asks users to extend the model
  energy grid beyond the fitted detector band.
- [x] Inspect the installed XSPEC shift implementations: `vashift`/`vmshift`
  remap spectra supplied by XSPEC; they do not independently evaluate the
  underlying component outside the supplied energy grid. The WINE shift
  follows the same interface pattern, with its own relativistic convention.
- [x] Confirm the soft artifact with the two-zone example. Keeping identical
  1--10 keV output bins but adding exterior model energies restores finite
  continuum transmission at 1.1 and 1.3 keV; the restricted grid gave zero.
  Emission also loses input contributions near the boundary. This is a finite
  input-coverage effect, not evidence that the FFT needs extra wrap-around
  padding.
- [x] Prefer documented `energies extend` over hidden table loading or
  extrapolation inside generic convolution operators. An operator cannot
  recover the missing spectrum simply by allocating a wider internal array.
- [x] Clarify the distinction between the model calculation grid and fitted /
  plotted energies for absorption shifts, additive shifts, and emission
  broadening. Extend the relevant low and high boundaries for the velocities
  and redshifts actually allowed in the test, within the table energy support.
- [x] Correct the current absorption warning, which describes endpoint
  extrapolation even though uncovered transmission bins can become zero.
  Do not replace missing transmission by zero or unity as an assumed physical
  continuation. Preserve the documented finite-grid contract.
- [x] Give a concise example and a convergence check: enlarge the model grid
  further and verify stability within the science band. Avoid blindly enlarging
  the grid to extreme energies or silently changing user-selected resolution.

Validation on 2026-10-27: the rebuilt package passed
`xspec/test_energy_coverage.py` (five shift/broadening cases, including an
analytic constant-transmission reference) and `xspec/test_geometry_models.py`.
The new check preserves 1--10 keV output bins and log spacing while expanding
input coverage; its maximum relative change between extended grids was below
`7e-13`. The author's two-zone `example.xcm` also agreed within `6e-13` between
0.1--31.62 and 0.01--31.62 keV calculation grids at matched log spacing.
These are fixed-parameter spectral checks, not fit-statistic or narrow-feature
sampling certification. The change updates documentation and reminder text;
it does not change convolution mathematics or introduce extrapolation.

Reference: [XSPEC energies](https://heasarc.gsfc.nasa.gov/docs/software/xspec/manual/XSenergies.html).

### 2. Separate geometry failure from fit initialization

- [x] Confirm that `steppar best` resets non-stepped parameters, but does not
  itself catch invalid model evaluations and skip those grid points.
- [x] Run an isolated infrastructure check with `windemlos*powerlaw`, a
  60-channel dummy response, fixed nuisance parameters, and XSPEC 12.15.1.
  `steppar best 3 60 80 2 4 10 30 2` returned nine finite statistics.
  `steppar best 3 40 80 2 4 20 60 2` reached five valid evaluations, then
  failed at `tout = tin`; `tclout steppar statistic` returned nine zeros,
  including positions already evaluated. The native XSPEC session remained
  alive. This is a scan/export failure, not a physical zero-statistic surface.
  This check did not establish nuisance-parameter convergence or diagnose the
  full PDS 456 reproducer.
- [x] Reproduce the `err.xcm` fitting failure with the minimum necessary
  evaluations; distinguish non-finite `uincl`, inverted opening angles, and
  numerical derivative probes. A temporary diagnostic build reproduced it in
  the first serial fit iteration: `uincl=0.8495` is valid, but `tin=89.49` exceeds
  `tout=52.96` deg. Improve diagnostic specificity in the selected implementation.
- [ ] Adopt a coherent policy for the zero-angular-volume limit `tout = tin`
  and the inverted geometry `tout < tin`; do not make every forbidden geometry
  indistinguishable from a physically absent emitter.
- [ ] Check error propagation in both native XSPEC and PyXSPEC. An isolated
  invalid evaluation previously terminated PyXSPEC with an uncaught C++
  exception, while native XSPEC caught the error.
- [ ] Initially scan only `tin < tout`, restarting nuisance parameters from
  the same accepted fit at each point and recording the fitted statistic and
  status immediately. Store forbidden locations as masked entries, never as
  physical zero statistics. A permitted point whose fit fails needs a distinct
  status from a geometrically excluded point.
- [ ] Check a few permitted points in reverse scan order before increasing
  grid density. Plot both axes over 0--90 degrees with the excluded triangle
  masked. Profile nuisance parameters for confidence contours.
- [x] Implement the agreed direction with explicit physical `Cf`, direct
  `incl` for generic/LoS models, and `view` plus fixed `region` for off models.
  Keep explicit-angle access in Python. Do not add legacy model variants.
  Use the agreed name `skew` for the LoS cavity coordinate, distinct from
  generic/off `hole`. Values above 0.5 favor equatorial fractional coverage;
  below 0.5 favor polar fractional coverage. The neutral point means equal
  fractional coverage, not angular symmetry around the sightline.
- [x] Verify geometry against an independent numerical root and spectra against
  the former explicit-angle implementation. Check native physical-parameter
  errors, a small rectangular profiled scan, and the first reproducer iteration
  in serial and parallel. See [first validation pass](geometry-validation.md).
  Distant-start convergence is not a release criterion for this shallow model;
  explore the likelihood rather than tuning the implementation to one start.
- [ ] Consolidate helpers under the proposed names `windgeom`, `winderror`,
  and `windsteppar`, using component selection and optional model/data-group
  context. Retain the tested profile `toler` behavior, reporting conventions,
  session preservation, and better-minimum diagnostics.
  `windgeom` values and the two-opening-angle `windsteppar` are implemented;
  derived opening-angle `winderror` and local covariance reporting are now
  implemented. Free-geometry off-wind inclination profiles remain open. See the
  [interface and implementation status](geometry-procedures.md).
- [x] Add generic structured-result export to xspeclab's `savesteppar` or an
  equivalent compatible adapter, reusing its FITS schema and existing plotting.
  Preserve masks and per-cell status; do not fabricate native parameter indices
  for derived axes. Test the full export/read/plot path and prevent display
  interpolation or smoothing from bridging forbidden cells.

The [geometry exploration](geometry-options.md) records the native-interface
inspection and temporary prototype tests. Zero emission outside the domain
allowed the toy profiled grid to finish and the real reproducer's first fit
iteration to recover into the allowed domain. It also allowed a toy fit to
remain at inverted angles. `NaN` only worked for a fixed-parameter surface,
not for profiled `steppar`. The author and Francesco subsequently agreed on
the physical `Cf`/inclination API, now implemented and checked in the first
validation pass. Physical-angle scans and structured export now pass targeted
integration checks; one-dimensional opening-angle errors and covariance
reporting now have targeted checks too. The earlier
alternatives are historical evidence.

Reference: [XSPEC steppar](https://heasarc.gsfc.nasa.gov/docs/software/xspec/manual/node118.html).

### 3. Build a small, transparent synthetic-response example

- [x] Export a real EPIC-pn ARF to a readable ASCII file, retaining bin bounds:
  `ENERG_LO_keV ENERG_HI_keV SPECRESP_cm2`. It represents one energy-dependent
  effective-area curve; retaining both edges avoids ambiguity when rebuilding
  bin averages or sharp instrumental structure.
- [x] Record instrument, source response, and any scale factor in its header.
  The local candidate inspected on 2026-10-27 covers approximately 0.05--16 keV:
  `/Users/pierpaolo/Analysis/Proposals/IRASF11119/XRISM_Cycle2/xspec_XMM_IRASF11119/spec_pn_opt610keV_pn_p04_r12.arf`.
- [x] Create an autonomous, concise notebook that reads the ASCII curve,
  conservatively rebins it if required, and writes the OGIP ARF. Preserve edges;
  avoid smooth fits that invent or erase instrumental discontinuities.
- [ ] Use a synthetic Gaussian RMF with constant **50 eV FWHM**, distinct from
  the ARF choice. Start with 10 eV channels, integrate probabilities over channel
  boundaries, and store the matrix sparsely. Check incident-energy sampling
  against the intrinsic spectral features, rather than equating channel width
  with a sufficient model grid by assumption.
- [x] Prefer EPIC-pn-like effective area for the example. Choosing a Chandra
  ARF would not improve an RMF whose resolution we already prescribe, and
  would reduce Fe-K counts. Label the hybrid as a synthetic detector, not an
  actual XMM or Chandra response. HETG does not have constant energy FWHM.
- [x] Start with the original ARF amplitude; a factor of 2--3 is trivial to
  apply later, with its value recorded. Select exposure from useful counts
  rather than treating area rescaling as model-evaluation optimization.
- [x] Validate FITS compatibility, response normalization including edge losses,
  and a monochromatic-line centroid/width. Generate an expected spectrum
  without Poisson fluctuations, then a reproducible Poisson realization from
  the author's `example.xcm`.
- [x] Test recovery from modest displaced starts in the two-zone example.
  Two four-parameter expected-count fits reach C below `4e-6`; one Poisson
  fit and a full-range 5x5 angular grid complete (15 valid, 10 excluded cells).
  These are targeted checks, not global convergence or scientific calibration.

Implemented: [response notebook](../example/aux/build_response.ipynb) and
[two-zone generator](../example/README.md). The current RMF has 4418
incident bins (at most 5 eV) and 1596 channels. Row losses, XSPEC line folding,
Poisson metadata/seed reproducibility and continuum-normalization recovery
pass. A tiny physical-angle grid is a functional check, not a confidence
measurement. Incident-grid convergence and geometric recovery remain separate
scientific checks. The cleaned public example now retains one Poisson spectrum;
the earlier expected-count tests are historical checks. Recipes are in
`example/aux/`, analyses in `example/results/`. Per the author's updated choice,
native `fakeit none` and `ftgrouppha optmin` with ten minimum counts replace
the PyXSPEC generator and fixed-width grouping. This adaptive grouping can
retain 10 eV channels, so it no longer imposes the former strict oversampling
cap of three. `.grp` associations are local and background is explicitly `NONE`.
Both requested 20x20 scans pass: 400 native points, and 210 physical angular
points with 190 masked exclusions. The plotting notebook retains the original
likelihood levels without smoothing or interpolating across the mask.

Reference: [HEASoft genrsp](https://heasarc.gsfc.nasa.gov/docs/software/ftools/caldb/caltools_src/genrsp/genrsp.html)
provides an independent Gaussian-response generator;
[Chandra HETG specifications](https://cxc.harvard.edu/proposer/POG/html/chap8.html)
distinguish effective area from spectral resolution.

### 4. Review numerical modules and optimize measured costs

- [ ] Trace preparation, physical geometry, relativistic transformation,
  profile, rebinning/convolution, and normalization in readable scientific steps,
  following `sci-coder`. Retain explicit units and photon/energy conventions.
- [ ] Evaluate the requested C implementation with a thin XSPEC adapter;
  language conversion alone is not a performance result. Preserve public names,
  parameter order, and links unless a separately justified change is agreed.
- [ ] Measure representative costs before caching or restructuring. FFT plan
  recreation and invariant geometry preparation are candidates to measure.
- [ ] Validate against analytic limits and independent references as well as
  Python consistency. Test accuracy and time before/after on the small example,
  then one representative full workflow.
- [ ] Check downstream geometry, covariance, and profiling procedures if model
  bounds, coordinates, or error behavior change. Leave all changes uncommitted.

## Completion boundary

Update the draft after the targeted numerical checks and the author's two
PDS 456 fitting tests. Report separately what was assumed, measured, numerically
verified, and left for future transfer modelling. Keep existing quantitative
results until replacement fits have been supervised and accepted by the author.
