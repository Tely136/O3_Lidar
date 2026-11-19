# ozone_lidar
ozone retreival from the ozone lidar

### Aerosol correction

Pass aerosol settings as the third argument to `process` to enable the
iterative Kuang et al. (2011) correction. Set the same struct in
`o3_process.m` to enable it for batch processing. For example:

```matlab
aerosol = struct('start_km',1.2,'reference_km',6,'reference_beta',1.67e-7);
out = process(input_path,output_path,aerosol);
```

Choose the lower limit where the far-range offline photon-count signal is
fully overlapped and usable, and set the offline aerosol backscatter at the
reference altitude in m⁻¹ sr⁻¹. The inversion uses a lidar ratio of 60 sr
and an Ångström exponent of 0.5. The retrieved aerosol profile supplies the
near- and far-range corrections before ozone merging.

`out.aerosol` records the offline aerosol backscatter and the ozone
number-density adjustments. NaNs in the adjustments indicate bins that
retain the Rayleigh-only retrieval, including bins without a full derivative
window and profiles where the inversion fails to converge. The existing
photon-count uncertainty budget excludes uncertainty in the aerosol correction.
