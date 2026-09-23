# Reproducing the fid chain in Python

What the code in this directory computes, how each piece corresponds to the C tool it
follows, and where it deliberately does something else. The pipeline itself is described
on the ska-wiki pages *Aspect L1 pipeline*, *Fid light corrections* and *Predicted fid
positions*; this document is only about the correspondence.

## What runs where

| Python | C | what it does |
|---|---|---|
| `fid_processing.get_caldb_calalign` | `asp_get_calib`, `load_aca_align.c`, `load_acis_fidcorr.c` | picks the CALDB alignment file applicable to the observation date and pulls the four extensions: the instrument alignment row, the per-fid table, the ACIS thermal constants and the periscope coefficients |
| `fid_processing.nominal_fid_positions` | `aca_id_image`, `Read_CALALIGN_Write_Fidprops()` | the SIM-Z correction polynomial applied to each fid's nominal LSI position |
| `fid_processing.estimated_fid_positions` | `asp_solve`, `ACA_Solve::calcFidPos` | LSI → STT → STF → FC → ACA, giving the predicted direction cosine |
| `fid_processing.process_centroids` | `aca_filter_centr` then `aca_corr_fid` | smoothing, the ACIS thermal correction, the periscope correction, and the measured direction cosine |
| `fid_smoothing.sg_smooth_per_slot` | `asp_centr_smooth` / `asp_sg_smooth` in `asputils` | per-slot Savitzky-Golay |
| `fid_smoothing.hann_smooth` | `asp_slide_smooth` in `asputils` | the sliding Hann window used on the periscope gradients |
| `fid_deltas_mica.get_fid_deltas_mica` | — | the same predicted chain built from ACACAL and the FIDPROPS header instead of CALDB, as an independent check |

The chain starts one stage later than the pipeline does: it reads `ang_y` and `ang_z`
from ACACENT, so the pixel-to-angle conversion of `aca_corr_centr` — the field distortion
polynomial, and the CTI and first-moment adjustments that precede it — has already
happened. The appendix records that polynomial anyway, because it is not written down
anywhere else we control.

## Where this differs from the pipeline, on purpose

These are standing differences, not bugs to fix quietly. Each one will show up again for
anyone who compares this code against the C or against the archived products.

1. **DH temperature.** The temperature series driving the thermal correction is built
   here from the arithmetic mean of `counts` across fid slots at each time. The C takes a
   windowed median per fid first and then averages across fids. With a cosmic ray in the
   window the two differ by tens of counts.
2. **Periscope baseline.** Here the mean subtracted from each gradient is the mean of the
   Hann-smoothed series over the observation. The C computes the means from the raw
   gradients and freezes them before smoothing.
3. **The temperature grid.** The C evaluates the temperature series on a fixed 32.8 s
   grid and assigns each centroid the most recent grid value — a step function in time.
   This code works per unique centroid time.
4. **SIM position.** `stt0_stf` is read from the obspar. The pipeline reads the same
   quantity from the SIM coordinate product, taking the first row in the interval whose
   three axes are all flagged good, and treats bad quality as fatal.
5. **Smoothing edges.** `savgol_filter(..., mode="nearest")` stands in for the C's
   replacement of bad and end samples by the nearest good one. They agree in the interior;
   the ends have not been compared sample by sample.

## Conventions that are easy to get wrong

- The correction polynomial's coefficients are in **ascending** order, so
  `numpy.polynomial.polynomial.polyval`, never `numpy.polyval`.
- Its variable is the SIM-Z offset in mm. It is geometry, not temperature.
- The C rotation matrix about the optical axis is the transpose of the right-hand-rule
  form, so the sign of `dtheta` flips if a conventional rotation is used.
- The fid identifiers in the star catalogue run 1–14 across the observatory; the
  calibration table numbers them within the instrument, 1–6 for ACIS and 1–4 for HRC.
- ACACENT angles are degrees on disk. Everything here converts to arcsec only at output.

## Status

The two chains run and produce a per-sample comparison (`fid_processing.compare_obsid`),
but the numbers have never been checked against an independent computation end to end.
Obsid 29878 is the interesting case, because its large SIM-Z offset is what exercises the
correction polynomial. Until that comparison is done and the residuals explained, treat
the output as a reading aid.

## Appendix — the field distortion polynomial

`aca_corr_centr` converts a centroid in CCD pixels to an ACA angle with a 20-coefficient
cubic in column `C`, row `R` and CCD temperature `T`, evaluated separately for the Y and
Z angles and with separate coefficient sets for stars and for fids
(`asputils/field_distortion.c`). A temperature above 150 is taken to be kelvin and
converted to °C first. In coefficient order:

```
c0
+ c1*C     + c2*R     + c3*T
+ c4*C^2   + c5*C*R   + c6*C*T   + c7*R^2   + c8*R*T   + c9*T^2
+ c10*C^3  + c11*R*C^2 + c12*T*C^2 + c13*C*R^2 + c14*C*R*T + c15*C*T^2
+ c16*R^3  + c17*T*R^2 + c18*R*T^2 + c19*T^3
```

The coefficients come from the field distortion calibration and are copied into the
ACACAL product of the observation as `fd_y_star`, `fd_z_star`, `fd_y_fid`, `fd_z_fid`,
with the inverse transforms beside them.
