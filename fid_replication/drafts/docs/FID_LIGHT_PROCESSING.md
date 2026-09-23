# Fid Light Processing in the Aspect Pipeline

This note documents how fiducial ("fid") lights are handled end-to-end in the
Chandra DS aspect pipeline (`dstools/asp/`), from raw CCD images through the
final corrected angles that feed `asp_solve`. It complements
[ASPECT_PIPELINE_OVERVIEW.md](ASPECT_PIPELINE_OVERVIEW.md) with a fid-specific
focus.

## 1. What fid lights are

Fiducial lights are six LEDs (four for HRC instruments, six for ACIS) mounted
on each of Chandra's science instruments. They are imaged by the Aspect
Camera Assembly (ACA) via the periscope and the Fiducial Transfer System
(FTS). Their nominal positions on the ACA focal plane are known from
calibration, so the difference between "where a fid light appears" and
"where it should appear" is a direct measurement of ACA-to-SI misalignment,
thermal drift, and periscope flexure. The aspect solution uses this
information together with guide-star measurements to anchor the pointing
reconstruction.

Each fid is identified by a string like `"ACIS-5"` or `"HRC-I-2"`.

## 2. The FIDPROPS data product

FIDPROPS is a multi-extension FITS table with one row per fid light used in
the observation. The struct definition is in
[dstools/asp/aspect_lib/fidprops.h](../dstools/asp/aspect_lib/fidprops.h).

Per-fid columns of interest:

| Column | Units | Description | Origin |
|--------|-------|-------------|--------|
| `id_string` | — | e.g. `"ACIS-5"` | `aca_id_image` (from matching ACA images) |
| `id_num` | — | 1–14 | `aca_id_image` |
| `p_lsi[3]` | mm | fid position in Local Science Instrument frame | CALALIGN |
| `slot` | — | ACA image slot carrying this fid | `aca_id_image` |
| `mag_i_cmd`, `mag_i_avg/min/max` | mag | commanded and observed mags | photometry |
| **`ang_y_nom`** | deg | nominal Y-angle (yag) in ACA frame | **CALALIGN ACISFIDCORR** |
| **`ang_z_nom`** | deg | nominal Z-angle (zag) in ACA frame | **CALALIGN ACISFIDCORR** |
| `id_status` | — | good / marginal / bad | fit quality |

Header keywords on the FIDPROPS extension carry additional calibration
constants used later in the pipeline: `dht_base`, `dht_mean`, `deg_per_cnt`,
`fid_cte`, `fid_y_coe`, `fid_z_coe` (thermal model), and the periscope
coefficients `peri_{y,z}_grd{3,6}`. All of them originate in CALALIGN and
are copied through by `asp_get_calib`.

## 3. Pipeline order

From [runasp/runasp.py:20-56](../runasp/runasp.py#L20-L56), the steps that
touch fid lights, in order:

```
 5  create_props_files          aca_id_image        create FIDPROPS, GSPROPS
 6  get_calibration_data        asp_get_calib       fill in CALDB-derived columns
 9  create_ACADATA              aca_read_data       L0 telemetry -> ACA images
10  apply_ccd_corrections       aca_corr_ccd        dark / flat / responsivity
14  reapply_ccd_corrections     aca_corr_ccd        rerun with updated ACACAL
15  calculate_centroids         aca_calc_centr      images -> pixel centroids
16  sort_centroids
17  replace_acacent
18  apply_centroid_corrections  aca_corr_centr      pixel -> ACA angles
19  filter_centr                aca_filter_centr    Savitzky-Golay smoothing
20  correct_acis_fids           aca_corr_fid        thermal + periscope correction
21  correct_properties
22  ...
24  create_aspect_solution      asp_solve           attitude solution
```

## 4. FIDPROPS construction

FIDPROPS is written in two phases.

### Phase A — `aca_id_image` creates the file

Entry point:
[dstools/asp/aca_id_image/aca_id_image.cc:261](../dstools/asp/aca_id_image/aca_id_image.cc#L261)
calls `obj.Read_CALALIGN_Write_Fidprops()`. This reads the ACA telemetry to
identify which slots actually contain fid lights, cross-references against
CALALIGN, and emits a FIDPROPS file populated with `id_string`, `id_num`,
`p_lsi`, `slot`, and initial magnitude fields. Columns that will be filled
in the next phase (`ang_y_nom`, `ang_z_nom`, and the thermal/periscope
header keywords) are written as sentinel / default values here.

### Phase B — `asp_get_calib` fills in CALDB-derived columns

Entry point:
[dstools/asp/asp_get_calib/load_acis_fidcorr.c:28](../dstools/asp/asp_get_calib/load_acis_fidcorr.c#L28)
(`Load_ACIS_Fidcorr`). This tool reopens FIDPROPS READWRITE and populates the
ACIS fid-correction columns from CALALIGN's ACISFIDCORR extension.

The logic (lines 183-215):

```c
while ( FIDPROPS_GetRow( &fidpr )) {
  /* Match FIDPROPS id_string (e.g. "ACIS-5") with CALALIGN.detector. */
  for ( ii = 0; ii < nfidcorr; ii++ ) {
    if ( strncmp( fidpr.id_string, fidcorrdata[ii].detector, slen ) == 0 ) {
      fid_num = atoi( strrchr(fidpr.id_string,'-')+1 );  /* "ACIS-5" -> 5 */
      found = dmTRUE;
      break;
    }
  }
  if ( found ) {
    long ndx = fid_num - 1;                 /* fid numbers are 1-based */
    fidpr.ang_y_nom = fidcorrdata[ii].fid_y_ang_nom[ndx];
    fidpr.ang_z_nom = fidcorrdata[ii].fid_z_ang_nom[ndx];
    /* ... also dht_base, deg_per_cnt, fid_cte, fid_y_coe, fid_z_coe */
  } else {
    ds_MAKE_FNAN( fidpr.ang_y_nom );
    ds_MAKE_FNAN( fidpr.ang_z_nom );
  }
  FIDPROPS_PutRow( &fidpr, "FIDPROPS" );
}
```

No trigonometry, no frame conversion, no interpolation. For the columns
shown above, the value in CALALIGN is the value in FIDPROPS.

One field that `Load_ACIS_Fidcorr` writes is **not** a CALDB copy: the
`dht_mean` header keyword. It is set just before the row loop at
[load_acis_fidcorr.c:181](../dstools/asp/asp_get_calib/load_acis_fidcorr.c#L181)
via `calc_average_dhtemp()`, which averages good-quality `c1bat` and `c1bbt`
samples from the ACIS engineering (L0) file over `[tstart, tstop]` and
converts K→°C. It is a per-observation scalar, and it becomes the DC baseline
around which `aca_corr_fid` later builds the time-resolved DH temperature
series (see §6).

### CALALIGN provenance

`fid_y_ang_nom[6]` / `fid_z_ang_nom[6]` live in the **ACISFIDCORR** extension
of CALALIGN, one row per detector (ACIS-I, ACIS-S, HRC-I, HRC-S). They are
**externally supplied** — nothing in the `ds` repo writes them:

- Only six files in the tree reference these columns, all on the reader
  side (see [aspect_lib/calalign.h:57-58](../dstools/asp/aspect_lib/calalign.h#L57-L58)).
- `asp_calc_boresight` does produce an updated CALALIGN (intended for CALDB
  re-ingestion), but it marks `fid_{y,z}_ang_nom` as `NOACCESS` on both
  input and output
  ([asp_calc_boresight.cc:77-78, 799-800](../dstools/asp/asp_calc_boresight/asp_calc_boresight.cc#L77-L78));
  it only modifies the ACA/FTS misalignment matrices.
- `CALALIGN_Create`, `CALALIGN_InsertRow` exist in the aspect I/O library
  but have zero callers. `CALALIGN_PutRow` is called only from
  `asp_calc_boresight`.

Determination of these values (pre-launch optical calibration, on-orbit fid
observation reductions, periscope/thermal analysis) happens in the Chandra
calibration team's workflow, outside this repository. Revisions do happen:
the **PERIFIDCORR** sibling extension was added in 2008 to support the
periscope correction — see the history comment in
[load_acis_fidcorr.c:17-20](../dstools/asp/asp_get_calib/load_acis_fidcorr.c#L17-L20).

## 5. How measured fid angles are produced: the `ang_*_sm` chain

ACACENT carries two pairs of angle columns
([aspect_lib/acacent.h:74-77](../dstools/asp/aspect_lib/acacent.h#L74-L77)):

```c
float ang_y;       /* raw centroid angle, ACA y */
float ang_z;       /* raw centroid angle, ACA z */
float ang_y_sm;    /* smoothed */
float ang_z_sm;    /* smoothed */
```

The `_sm` suffix is "smoothed". They are produced by the following chain:

### Stage 1 — `aca_calc_centr` (image → pixel)

Input: CCD-corrected ACADATA (after `aca_corr_ccd`). Output: ACACENT rows
with `cent_i`, `cent_j`, `flux`, `chi`, and the algorithm tag `alg`.

Fit algorithms available per slot: first-moment, PSF fit, Gaussian,
elliptical Gaussian. At this stage `ang_y` / `ang_z` are set to 0 — they
have not yet been computed.

### Stage 2 — `aca_corr_centr` (pixel → ACA angle)

This is where `ang_y` / `ang_z` acquire meaningful values. Core loop at
[aca_corr_centr.c:243-277](../dstools/asp/aca_corr_centr/aca_corr_centr.c#L243-L277):

```c
if (guide star) {
  acacent.ang_y = Model_Distortion_P(cent_i, cent_j, T, acacal.fd_y_star);
  acacent.ang_z = Model_Distortion_P(cent_i, cent_j, T, acacal.fd_z_star);
} else if (fid) {
  acacent.ang_y = Model_Distortion_P(cent_i, cent_j, T, acacal.fd_y_fid);
  acacent.ang_z = Model_Distortion_P(cent_i, cent_j, T, acacal.fd_z_fid);
}
acacent.ang_y += acacal.color_corr_i[slot];
acacent.ang_z += acacal.color_corr_j[slot];
```

`Model_Distortion_P()` lives at
[asputils/field_distortion.c:16-56](../dstools/asp/asputils/field_distortion.c#L16-L56)
and is a 20-coefficient degree-3 polynomial in `(R = cent_i, C = cent_j,
T = CCD_temperature)` with all cross terms up to cubic. CCD temperature is
auto-converted K → °C when `T > 150`.

The polynomial is the ACA camera's optical-distortion and plate-scale model.
Same functional form for stars and fids; different coefficient arrays.

#### Distortion-coefficient provenance

The coefficient arrays `fd_y_star`, `fd_z_star`, `fd_y_fid`, `fd_z_fid`
(plus the inverse-direction arrays `fd_r_star`, `fd_c_star`, `fd_r_fid`,
`fd_c_fid`) are each 20-element floats stored in the ACACAL runtime
structure
([aspect_lib/acacal.h:54-61](../dstools/asp/aspect_lib/acacal.h#L54-L61)).

They originate in the **CALFDC** CALDB product. At pipeline start,
`asp_get_calib` runs `Load_FDC()` at
[load_field_distortion.c:76-90](../dstools/asp/asp_get_calib/load_field_distortion.c#L76-L90),
which straight-copies CALFDC's coefficient arrays into ACACAL:

```c
for (ii = 0; ii < FIELD_ORDER; ii++) {
  acacal->fd_y_star[ii] = calfdc.fd_y_star[ii];
  acacal->fd_z_star[ii] = calfdc.fd_z_star[ii];
  acacal->fd_y_fid[ii]  = calfdc.fd_y_fid[ii];
  acacal->fd_z_fid[ii]  = calfdc.fd_z_fid[ii];
  acacal->fd_r_star[ii] = calfdc.fd_r_star[ii];
  acacal->fd_c_star[ii] = calfdc.fd_c_star[ii];
  acacal->fd_r_fid[ii]  = calfdc.fd_r_fid[ii];
  acacal->fd_c_fid[ii]  = calfdc.fd_c_fid[ii];
}
```

Same pattern as CALALIGN: no repo code writes CALFDC. `CALFDC_Create` /
`CALFDC_PutRow` are defined but never called. `aca_find_hotpix` produces an
updated ACACAL but preserves the distortion arrays verbatim via `memcpy`
([aca_find_hotpix/output_file_functions.c:266-273](../dstools/asp/aca_find_hotpix/output_file_functions.c#L266-L273)).

### Stage 3 — `aca_filter_centr` (smoothing)

Driver:
[aca_filter_centr/ACA_Fidmotion.cc:257-357](../dstools/asp/aca_filter_centr/ACA_Fidmotion.cc#L257-L357).
It runs a **Savitzky-Golay** filter with iterative sigma-clipping on the
`ang_y` / `ang_z` series per slot, writes `ang_y_sm` / `ang_z_sm`, and sets
outlier status bits `SIGMA_REJECT_Y` (0x40) / `SIGMA_REJECT_Z` (0x80) from
[acacent.h:23-26](../dstools/asp/aspect_lib/acacent.h#L23-L26).

Parameters (from `aca_filter_centr.par`):

- `SG_order = 4` — polynomial order of the local fit
- `T_sm_star = T_sm_fid = 40 s` — smoothing window (converted to an odd
  sample count as `box_size = T_sm / integ_time`)
- `reject_sigma = 3.0`
- `reject_iter = 5`

Points flagged by the upstream centroid algorithm (`IMG_ALG_BAD_CENTR`) are
substituted with nearest-good-neighbor values before filtering so they do
not bias the fit. Magnitude is smoothed in the same pass (stored internally
on the `ACA_Fidlight` object; ACACENT does not carry `mag_sm`).

## 6. Thermal + periscope correction: `aca_corr_fid`

This is the tool that actually uses `ang_y_nom` / `ang_z_nom`. Source:
[dstools/asp/aca_corr_fid/tool_functions.c](../dstools/asp/aca_corr_fid/tool_functions.c).

### The thermal model

The ACA focal plane expands radially with temperature. The expansion center
in focal-plane angles, `(fid_y_coe, fid_z_coe)`, is a calibration quantity
from CALALIGN / FIDPROPS. Fids sit at known nominal angles
`(ang_y_nom, ang_z_nom)`. Over an observation, their measured positions
drift; averaging the drift across all fids localizes the *observed*
expansion center relative to the calibrated one.

Step 1 — observed center of expansion
([tool_functions.c:591-657](../dstools/asp/aca_corr_fid/tool_functions.c#L591-L657)):

```c
for each fid {
  ang_y_offset = avg(ang_y_sm) - ang_y_nom;    /* line 623 */
  ang_z_offset = avg(ang_z_sm) - ang_z_nom;    /* line 624 */
  accumulate ...
}
info->obs_y_coe = fid_y_coe + mean(ang_y_offset);   /* line 642 */
info->obs_z_coe = fid_z_coe + mean(ang_z_offset);   /* line 643 */
```

Step 2 — per-centroid thermal correction
([tool_functions.c:790-802](../dstools/asp/aca_corr_fid/tool_functions.c#L790-L802)):

```c
corr  = (dhtemps[ndx] - dht_base) * fid_cte;
ycorr = corr * (acacent.ang_y_sm - obs_y_coe);
zcorr = corr * (acacent.ang_z_sm - obs_z_coe);
acacent.ang_y_sm -= ycorr;
acacent.ang_z_sm -= zcorr;
```

The correction scales linearly with the temperature difference from the
calibration reference (`dht_base`), linearly with the thermal coefficient
`fid_cte`, and linearly with the radial distance from the observed
expansion center. Without the nominal angles, `obs_y_coe` / `obs_z_coe`
cannot be pinned, and the correction has no reference frame.

`dhtemps[]` is a per-sample ACIS detector-housing temperature series —
**computed from the fid LED brightness itself**, not from housing-temperature
telemetry. Mechanism
([aca_corr_fid/tool_functions.c:434-571](../dstools/asp/aca_corr_fid/tool_functions.c#L434-L571)):
each sample's DH temperature is the observation mean `dht_mean` plus a
fluctuation term proportional to how much the fid flux at that sample
deviates from the observation-wide average flux:

```c
/* per sample isamp, after taking the median of counts across fids */
dhtemps[isamp] = info->dht_mean
               + (avg_med_cnts[isamp] - overall_avg_cnts) * info->deg_per_cnt;
```

The conversion coefficient `deg_per_cnt` is
`"Degrees per fid light count (integrated e-)"` — the physical slope between
LED brightness and temperature
([fidprops.h:36](../dstools/asp/aspect_lib/fidprops.h#L36)). The fids dim as
they heat up and the base-plate expansion moves them; both effects share the
same underlying temperature, so the fids serve as their own thermometer.

Only two inputs to the formula come from CALDB (`dht_base` = reference
temperature at which the nominal angles were calibrated, and `deg_per_cnt` =
the brightness-to-temperature slope). The absolute baseline `dht_mean` is a
per-observation scalar computed once in `asp_get_calib/load_acis_fidcorr.c`
by `calc_average_dhtemp()` ([lines 261-340](../dstools/asp/asp_get_calib/load_acis_fidcorr.c#L261-L340))
from the ACIS engineering channels `c1bat` and `c1bbt`. Fid flux provides
the time-resolved shape of `dhtemps[]`; engineering telemetry anchors its
DC value.

### The periscope correction

Applied after the thermal step, at
[tool_functions.c:858-867](../dstools/asp/aca_corr_fid/tool_functions.c#L858-L867):

```c
ycorr = acacent.ang_y_sm;
ycorr -= (smoothGradient3[ii] - meanGradient3) * peri_y_grd3;
ycorr -= (smoothGradient6[ii] - meanGradient6) * peri_y_grd6;
acacent.ang_y_sm = ycorr;
/* same for Z with peri_z_grd3 / peri_z_grd6 */
```

The two OBC engineering channels `OOBAGRD3` and `OOBAGRD6` are solar-panel
thermal gradients; they serve as a proxy for spacecraft thermal state.
`peri_{y,z}_grd{3,6}` are bilinear coefficients from CALALIGN's
**PERIFIDCORR** extension (added in 2008). This step does not read
`ang_*_nom` directly; it operates on the already thermally-corrected
angles.

### Writing back

The corrected centroid is written back to ACACENT in place:
[tool_functions.c:878](../dstools/asp/aca_corr_fid/tool_functions.c#L878)
`ACACENT_PutRow( &acacent, "ACACENT" )`. The columns that change are
`ang_y_sm` and `ang_z_sm`. `ang_y` / `ang_z` (the raw versions) remain
untouched.

## 7. Downstream: `asp_solve`

The corrected `ang_y_sm` / `ang_z_sm` are consumed by `asp_solve`, which
fits per-time-bin SIM-frame offsets `(dy, dz, dtheta)` against a
predicted fid direction-cosine reconstructed from FIDPROPS calibration.
See [PREDICTED_FID_AND_ASP_SOLVE.md](PREDICTED_FID_AND_ASP_SOLVE.md) for
the predicted chain, the chi-square fit, and the resulting ASPSOL
columns.

## 8. End-to-end data flow

```
 CALDB                                    (external origin)
   │
   ├── CALALIGN ACISFIDCORR  ──▶ asp_get_calib / Load_ACIS_Fidcorr
   │      fid_y_ang_nom[6]                        │ copy
   │      fid_z_ang_nom[6]                        ▼
   │      dh_temp_base, deg_per_count,     FIDPROPS per-fid row
   │      fid_CTE, fid_{y,z}_center_nom    (ang_y_nom, ang_z_nom, ...)
   │
   ├── CALALIGN PERIFIDCORR  ──▶ asp_get_calib
   │      peri_{y,z}_grd{3,6}                  FIDPROPS header keywords
   │
   └── CALFDC                ──▶ asp_get_calib / Load_FDC
          fd_{y,z}_{star,fid}                  ACACAL runtime struct
          fd_{r,c}_{star,fid}

 L0 ACIS eng. ──▶ asp_get_calib / calc_average_dhtemp
   c1bat, c1bbt        avg over [tstart, tstop]
                                │
                                ▼
                         FIDPROPS.dht_mean (per-obs scalar)

 L0 ACA telemetry ──▶ aca_read_data ──▶ ACADATA raw images
                                            │
                                  aca_corr_ccd (dark, flat, responsivity)
                                            │
                                            ▼
                                     ACADATA corrected
                                            │
                             aca_calc_centr (PSF / Gaussian / FM fit)
                                            │
                                            ▼
                            ACACENT (cent_i, cent_j, flux, chi)
                                            │
                         aca_corr_centr (Model_Distortion_P + color)
                                            │
                                            ▼
                               ACACENT (+ ang_y, ang_z)
                                            │
                   aca_filter_centr (Savitzky-Golay SG4, 40 s window)
                                            │
                                            ▼
                         ACACENT (+ ang_y_sm, ang_z_sm)
                                            │
                                            │   ACACENT.flux ──┐
                                            │                  ▼
                                            │          calc_dh_temp  (fids as thermometer:
                                            │                         dht_mean + Δflux·deg_per_cnt)
                                            │                  │
                                            │                  ▼
                                            │          dhtemps[] (per-sample DH temp)
                                            │                  │
                   aca_corr_fid (thermal + periscope) ◀─────────┘
                                            │
                                            ▼
                         ACACENT (ang_y_sm, ang_z_sm updated in place)
                                            │
                                            ▼
                                        asp_solve
                                            │
                                            ▼
                                         ASPSOL
```

## 9. Key files

| Concern | File |
|---------|------|
| FIDPROPS struct | [dstools/asp/aspect_lib/fidprops.h](../dstools/asp/aspect_lib/fidprops.h) |
| FIDPROPS FITS I/O | [dstools/asp/aspect_lib/fidprops.c](../dstools/asp/aspect_lib/fidprops.c) |
| CALALIGN struct | [dstools/asp/aspect_lib/calalign.h](../dstools/asp/aspect_lib/calalign.h) |
| CALALIGN FITS I/O | [dstools/asp/aspect_lib/calalign.c](../dstools/asp/aspect_lib/calalign.c) |
| FITS column templates | [dstools/asp/vario/src/vio_data_product.temp](../dstools/asp/vario/src/vio_data_product.temp) |
| Initial FIDPROPS write | [dstools/asp/aca_id_image/aca_id_image.cc](../dstools/asp/aca_id_image/aca_id_image.cc) |
| ACISFIDCORR → FIDPROPS transfer | [dstools/asp/asp_get_calib/load_acis_fidcorr.c](../dstools/asp/asp_get_calib/load_acis_fidcorr.c) |
| CALFDC → ACACAL transfer | [dstools/asp/asp_get_calib/load_field_distortion.c](../dstools/asp/asp_get_calib/load_field_distortion.c) |
| Pixel → angle polynomial | [dstools/asp/asputils/field_distortion.c](../dstools/asp/asputils/field_distortion.c) |
| Pixel → angle application | [dstools/asp/aca_corr_centr/aca_corr_centr.c](../dstools/asp/aca_corr_centr/aca_corr_centr.c) |
| Savitzky-Golay smoothing | [dstools/asp/aca_filter_centr/ACA_Fidmotion.cc](../dstools/asp/aca_filter_centr/ACA_Fidmotion.cc) |
| Thermal + periscope correction | [dstools/asp/aca_corr_fid/tool_functions.c](../dstools/asp/aca_corr_fid/tool_functions.c) |
| Pipeline driver | [runasp/runasp.py](../runasp/runasp.py) |

## 10. What is *not* in this repository

The calibration-time determination of these quantities lives elsewhere
(Chandra calibration team products ingested into CALDB):

- CALALIGN ACISFIDCORR: `fid_{y,z}_ang_nom`, `fid_{y,z}_center_nom`,
  `fid_CTE`, `dh_temp_base`, `deg_per_count`.
- CALALIGN PERIFIDCORR: `peri_{y,z}_grd{3,6}`.
- CALFDC: `fd_{y,z}_{star,fid}` and the inverse `fd_{r,c}_{star,fid}`.

The `ds` repo is a consumer of these *reference* constants. It threads them
through ACACAL / FIDPROPS and applies them to measured centroids.

That does **not** mean the ACIS fid correction itself is externally supplied.
The thermal correction is estimated per observation inside this repository,
using the fid LEDs as their own thermometer:

- The time-resolved DH temperature series `dhtemps[]` is built from
  per-sample fid flux in `aca_corr_fid/calc_dh_temp`
  ([tool_functions.c:434-571](../dstools/asp/aca_corr_fid/tool_functions.c#L434-L571)),
  with `deg_per_cnt` converting brightness fluctuations to degrees.
- Its DC baseline `dht_mean` is set per observation in
  `asp_get_calib/calc_average_dhtemp` from the `c1bat`/`c1bbt` engineering
  channels ([load_acis_fidcorr.c:261-340](../dstools/asp/asp_get_calib/load_acis_fidcorr.c#L261-L340)).
- The observed center of expansion `obs_{y,z}_coe` is averaged from this
  observation's `(ang_*_sm - ang_*_nom)` across all fids
  ([tool_functions.c:591-657](../dstools/asp/aca_corr_fid/tool_functions.c#L591-L657)).

What CALDB supplies is the *reference frame* against which those measurements
are made: nominal per-fid angles, nominal expansion center, the CTE, and the
calibration reference temperature. Questions about how those reference
numbers themselves were measured — pre-launch optical metrology, on-orbit
fid-light observation reductions, HRMA optical modeling, thermal
characterization — need the calibration team's documentation, not this
codebase.
