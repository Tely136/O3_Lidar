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
| uDET | 28-29 | sqrt(accumulated counts); no additional averaging factor |
| uP_DET | 30-31 | Literal expression retained; counts/exposure contract and invalid raw counts documented |
| uP_SAT | 34-35 | Literal expression retained; elementwise shots/uTau, seconds and time correlation documented |
| uP_BKG | 41-42 | Elementwise profile support; caller now reverses slope/intercept SE order |
| uNO3_P_X_NC | 32,36,43 | Fixed missing squared SDs, factor 2, profile indexing; explicit vertical correlation mode |
| uNO3_P_X_FC | 38,45 | Implemented previous empty stub, preserving printed absolute-value placement |
| uNO3_dS_O3_R | 47 | Corrected factor 2 for single-path convention; elementwise arrays; printed ozone-density sign retained |
| uNO3_dS_M_S | 55 | Correlated analytical Rayleigh model; single-path ozone denominator |
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

## File accumulation and integration

Each group of `n_avg` files is treated as one longer acquisition. Preprocessing
sums RAW analog and photon-count data first, sums the corresponding shot
counts, then dead-time corrects the accumulated PC signal and fits backgrounds.
`n_shots` is read per file/channel from the Licel headers. For identical files,
`n_shots_total = n_avg*n_shots`; different exposures are also summed correctly.
Only complete file groups are used, as before. `ids_avg` retains group membership
and `times_avg` retains the mean timestamps for plotting.

All photon-count uncertainty formulas now receive accumulated counts:

```matlab
u_detection = uP_DET(P_sum,R_sum);
u_saturation = uP_SAT(P_sum,dz,n_shots_total,uTau);
```

There are no per-file uncertainty averaging loops or extra sqrt(n) factors.
Poisson noise uses sqrt(R_sum); a fixed dead-time uncertainty does not improve
in relative terms just because the exposure is longer. Background-fit SEs use
the accumulated-count normalization as well. Conversion to MHz or mV divides by
the same total shots, keeping the physical signal units consistent.

`photon_count_uncertainty` calculates all five budget components in one place:
DET, SAT, BKG, O3 and RAYLEIGH. FR and NR are separate calls with their own ozone
number density and derivative filters. The six signal uncertainties are passed
to `uNO3` for DIAL propagation; the two cross-section terms and five-component
quadrature total are calculated in the same budget function. Its options hold
dz, total shots L, uTau, dsigma, C, channel correlation, cross sections and their
relative standard uncertainties. The numbered paper formulas are unchanged.

Summing before nonlinear correction is intentionally different from correcting
each file before averaging when rates vary. This design reproduces one longer
raw accumulation; it assumes the combined interval is suitable for one effective
count-rate correction. No interpolation across missing acquisition times is added.
The output records `pc_R_sum`, corrected/background-subtracted `pc_P_sum`, and
`n_shots_total`, alongside the original raw data and grouping variables.

Dependency and BDM cross-section paths in `process.m` are resolved from the
script location so the relocated project can be tested consistently.

`process.m` now captures the fit SEs and correlations, uses C_nr for NR,
uses count-domain signals, and saves FR/NR component profiles and merged PC
uncertainty. Its old call used undefined pc_P/R_avg/uB/rb variables.
`merge_pc_uncertainty` applies Eq.95 to DET and Eq.97 to BKG/SAT before
combining; the integration follows retrieval weights and omitmissing fallback.
Dead-time uncertainty is enabled with a provisional **10% relative standard
uncertainty**: `tau = 1/(fd*1e6)` and `uTau = 0.10*tau`. With `fd=270 MHz`,
this means `tau=3.704 ns` and `uTau=0.370 ns`. The choice is an engineering
assumption, not a measured uncertainty or manufacturer tolerance. The nominal
dead time used to correct the signals is unchanged. The saved metadata records
both the nominal value and assumed relative/absolute standard uncertainty.
The same assumed magnitude is currently used for each channel; ON/OFF errors
remain independent and each channel's dead time remains correlated in time.

Licel's [PMT manual, section 6.6](https://licel.com/manuals/pmtmanual.pdf)
quotes approximately 3.7 ns / 270 MHz for older transient recorders and
3.06 ns / 326 MHz for a TR40-16bit-3U at specified settings. These nominal
values are not uncertainty estimates. Confirm the actual recorder and settings
before treating 270 MHz as calibrated.

To estimate the effective dead time, align simultaneous analog and RAW PC
profiles, select a region with linear analog response and moderate PC rates,
and fit `r = q/(1+tau*q)`, with `q=a*V+b`, in Hz. Fit gain/offset and tau
together; do not use already corrected or glued PC data as the reference.
Repeat per channel at the operational HV/discriminator settings. Estimate
standard uncertainty from the fit and repeatability, and check sensitivity to
the fitting interval and signal alignment. See Licel's
[analog/PC note](https://licel.com/manuals/analogpc.pdf), section 6.5.
A larger uTau does not make heavily saturated or model-invalid data reliable.

The processing budget also includes `O3` (Eq.47) and `RAYLEIGH` (Eq.55).
Both relative standard uncertainty settings are 0.02 in `process.m`, using
Table 3 as a working estimate, not as a calibration of this instrument.
Ozone wavelength errors are treated as independent random components, as in
Table 3's ozone example. Rayleigh wavelength errors are fully correlated
because both cross sections come from one analytical model. A common 2%
Rayleigh scale uncertainty therefore gives 2% of the Rayleigh correction,
not the quadrature of two independent wavelength errors.

Ozone uncertainty is evaluated using each range's ozone profile. Rayleigh
uncertainty is the same atmospheric profile for all ranges. All components are
calculated within the FR/NR budget function. Merging treats O3 as one shared
fractional uncertainty: its signed component follows the same linear weights
as the ozone profile, equivalent to evaluating Eq.47 on the merged density. The common Rayleigh component is carried
through unchanged. They receive no independent-noise reduction (Eqs.93-94).
The saved metadata records the input fractions and correlation assumptions.
The live profile plot shows these two additional components and the updated
subtotal; rerun processing to populate the fields in older saved files.

Air-density/temperature/pressure and interfering-gas components remain
excluded. A standard atmosphere is still treated as fixed, so the subtotal
is not the full atmospheric budget of Eq.99/100.
Mixing-ratio PC terms are number-density PC terms /Na *1e9. This conversion
cannot substitute for Eqs.82/84 when adding uncertainty in Na itself.

## Compatibility changes

- uNO3_dS_M_R's denominator is the ozone differential, not Rayleigh.
- uP_DET takes (P,R); uP_SAT takes (P,dz,L,uTau). L is total accumulated
  shots. uDET takes accumulated R directly. The old ids_avg/nAverage arguments
  are removed from the signal uncertainty functions.
- preprocess returns raw pc_R_sum ninth and n_shots_total tenth. The earlier
  per-file pc_P_unavg output is removed.
- photon_count_uncertainty takes (R,P,N,Na,uB,rb,id_on,id_off,options), with
  accumulated R/P, range-specific ozone N and air density Na. options.L is
  scalar or [1 x profile x channel]. Cross-section options are required.
- merge_pc_uncertainty now accepts all five component fields, including O3
  and RAYLEIGH, and returns their merged quadrature total.
- uNO3 still accepts the six signal uncertainty arrays explicitly; no fallback
  calculations or uncertainty overrides are used.

## Verification

Run from the repository root:

```matlab
addpath(genpath('proc_functions'));
results = runtests('tests/test_uncertainty.m');
assertSuccess(results);
```

Tests compare scalar/vector/multiple-profile results against independent
hand-calculated equations, check the literal signed cases, factor-of-two
convention, slope/intercept order, accumulated-count scaling, range endpoints, component merging, variable
shot totals, and a synthetic accumulation pipeline using real preprocessing.
