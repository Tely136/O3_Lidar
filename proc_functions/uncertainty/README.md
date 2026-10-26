# Photon-count uncertainty audit

Reference: Leblanc et al. (2016), AMT 9, 4051-4078, doi:10.5194/amt-9-4051-2016.
The numbered formulas are implemented **literally**, by explicit user choice.
This does not implement the entire uncertainty budget in the paper.

## Conventions

- Signals are **counts**, not MHz, with consistent exposure/normalization.
- Every cross-section differential passed here is **single-path**, in m^2:
  `dsigma = sigma_on - sigma_off`, matching `dial.m`.
  The paper's Eq.27 uses outgoing plus return paths. For elastic backscatter,
  its `DeltaSigma = 2*dsigma`. Thus derivative uncertainty has denominator
  `2*abs(dsigma)*dz`; the factor 2 in Eqs.47/57 cancels instead.
- Cross-section uncertainty inputs are absolute standard uncertainties in m^2.
- `uTau` is standard uncertainty in seconds; `fd` is inverse dead time in MHz.
  `1/(fd*1e6)` is the dead time itself, not its uncertainty.
- Data are altitude x profile; derivative filter matrices use column 2.
- `bg_correction` stores SEs as **[slope; intercept]**. `uP_BKG` takes
  **(intercept SE, slope SE, correlation, z)**.
- Paper merge weight `w` means LOW-range weight. `merge_profiles` returns
  HIGH-range weight; pass `1-w_high` to the uncertainty merge formulas.
- Do not use photon-count uncertainty for analog or analog/PC merged ozone.

## Function-by-function findings

| Function | Equations | Review/change |
|---|---|---|
| uDET | 28-29 | sqrt(total counts) correct; optional equal-exposure mean factor added |
| uP_DET | 30-31 | Literal expression retained; counts/exposure contract and invalid raw counts documented |
| uP_SAT | 34-35 | Literal expression retained; elementwise shots/uTau, seconds and time correlation documented |
| uP_BKG | 41-42 | Elementwise profile support; caller now reverses slope/intercept SE order |
| uNO3_P_X_NC | 32,36,43 | Fixed missing squared SDs, factor 2, profile indexing; explicit vertical correlation mode |
| uNO3_P_X_FC | 38,45 | Implemented previous empty stub, preserving printed absolute-value placement |
| uNO3_dS_O3_R | 47 | Corrected factor 2 for single-path convention; elementwise arrays; printed ozone-density sign retained |
| uNO3_dS_M_R | 57 | AIR density and OZONE differential required; former Rayleigh denominator was wrong |
| uNO3_Na1 | 81 | Formula correct; elementwise oxygen term, assumptions documented |
| uNO3_Na2_NC | 86 | Formula correct; elementwise oxygen term, independent pressure/temperature documented |
| uNO3_Na2_FC | 89 | Formula correct; elementwise oxygen term, positive full correlation documented |
| uNO3_MRG_NC | 95 | Formula correct; correlation is between ranges and weight is LOW-range |
| uNO3_MRG_FC | 97 | Formula correct; weight convention and component-wise merging documented |
| uNO3 | 32-46,99 subtotal | DET independent bins, BKG/SAT correlated bins; reversed SE order, configurable ON/OFF correlation |

## Literal-formula limitations (intentional)

The signed summation in printed Eqs.36/43 can be negative. Eqs.38/45 also
put the absolute value *inside* the sum, which can give a negative result.
These expressions are preserved, without adding an outer absolute value or
replacing them with an Eq.4 covariance/Jacobian calculation. Such component
values are signed formula outputs, not nonnegative standard uncertainties.
Eq.47 also retains the sign of the retrieved ozone density, as printed.
The quadrature subtotal squares those outputs as in Eq.99. Do not interpret
agreement with the printed equations as a proof of covariance correctness.
In particular, summing pointwise ON/OFF quadrature before filtering is not
in general the same as separately filtering independent channel sensitivities
and then combining their variances. A fitted background has a covariance
matrix determined by its two fit parameters, not merely by pointwise SDs.

Eqs.30/34 use P as printed. However, differentiating Eq.21 exactly at fixed B
would give `(Q/R)^2*sqrt(R)` and `c/(2*dz*L)*Q^2*uTau`, where `Q=P+B` is
saturation-corrected counts **before** background subtraction. The numbered
formulas are not silently replaced with those derivatives. Detection and
background estimates from the same fitted raw data also need covariance for
a fully general treatment. These are limitations of literal application.

## Acquisition averaging and integration

`photon_count_uncertainty` replays each raw acquisition's correction before
propagating Eqs.30/34. For a group of n independent detections,
`uDET_mean = sqrt(sum(uDET_i.^2))/n`. For one fixed uncertain dead time shared
across acquisitions, `uSAT_mean = mean(uSAT_i)`, without a sqrt(n) reduction.
The background is fitted once to the averaged corrected signal, as in the
current `preprocess`; its SEs must not be divided by sqrt(n) again.
The wrapper applies the group background to each acquisition, preserves
signed individual noise values and uses the numbered equations literally.
It rejects rates at/above fd in propagated DET/SAT outputs.

Dependency and BDM cross-section paths in `process.m` are resolved from the
script location so the relocated project can be tested consistently.

`process.m` now captures the fit SEs and correlations, uses C_nr for NR,
uses count-domain signals, and saves FR/NR component profiles and merged PC
uncertainty. Its old call used undefined pc_P/R_avg/uB/rb variables.
`merge_pc_uncertainty` applies Eq.95 to DET and Eq.97 to BKG/SAT before
combining; the integration follows retrieval weights and omitmissing fallback.
The original `uTau=0` setting remains: SAT is deliberately disabled until a
manufacturer/calibration standard uncertainty is supplied. No values for
cross-section, temperature, pressure or interfering-gas uncertainty have been
invented. These functions exist as building blocks but are not included in
the PC subtotal. Atmospheric uncertainty components must be merged unchanged
(Eqs.93-94) and added separately when extending this to a complete budget.
Mixing-ratio PC terms are number-density PC terms /Na *1e9. This conversion
cannot substitute for Eqs.82/84 when adding uncertainty in Na itself.

## Compatibility changes

- uNO3_dS_M_R's second argument now means **ozone differential**, not Rayleigh.
- uNO3_P_X_NC defaults to independent bins (Eq.32); use eighth argument 'full'
  for literal Eqs.36/43. FC implements full ON/OFF and vertical correlation.
- uP_BKG returns an altitude column for vector z; no caller transpose needed.
- uNO3 remains a DET/BKG/SAT subtotal; optional overrides support acquisition
  averaging. Passing nonlinear averaged signals alone is insufficient.

## Verification

Run from the repository root:

```matlab
addpath(genpath('proc_functions'));
results = runtests('tests/test_uncertainty.m');
assertSuccess(results);
```

Tests compare scalar/vector/multiple-profile results against independent
hand-calculated equations, check the literal signed cases, factor-of-two
convention, slope/intercept order, averaging, range endpoints, component
merging, and a synthetic acquisition pipeline using real preprocessing.
