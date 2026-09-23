# Fid-light correction stages — CalDB-driven reference

This document is a stage-by-stage reference for reimplementing the ACIS
fid-light centroid correction chain in Python, **reading calibration values
directly from CalDB products** (CALFDC, CALALIGN) instead of relying on the
pipeline-intermediate FIDPROPS header keywords as a calibration cache. It
complements [FID_LIGHT_PROCESSING.md](FID_LIGHT_PROCESSING.md), which covers
the same ground in narrative form and treats FIDPROPS as a first-class source.

Audience: a developer writing a standalone Python script that produces
corrected `ang_y_sm` / `ang_z_sm` from raw pixel centroids.

Guiding constraint: every calibration constant is sourced from a CalDB file
directly. FIDPROPS is used only as a per-observation **slot → fid-id map**
(one row per fid light actually identified in this observation), not as a
cache of constants.

---

## 1. Input inventory

The Python script opens seven files. Only two of them (CALFDC, CALALIGN) are
CalDB products; the others are per-observation telemetry or pipeline
intermediates.

| File | Class | Extension(s) read | What we take |
|---|---|---|---|
| CALFDC | CalDB | primary | `fd_y_fid[20]`, `fd_z_fid[20]` polynomial coefficients |
| CALALIGN | CalDB | ACISFIDCORR | `dh_temp_base`, `deg_per_count`, `fid_CTE`, `fid_y_center_nom`, `fid_z_center_nom`, `fid_y_ang_nom[6]`, `fid_z_ang_nom[6]` — one row per detector |
| CALALIGN | CalDB | PERIFIDCORR | `peri_y_oobagrd3`, `peri_y_oobagrd6`, `peri_z_oobagrd3`, `peri_z_oobagrd6` — one row per detector |
| ACACENT | L1 pipeline intermediate | ACACENT | per (slot, time): `cent_i`, `cent_j`, `flux`, `status`, `alg`, `time`, `slot`; header `t_int` |
| ACADATA | L1 pipeline intermediate | ACADATA | per slot, per sample: `temperature` (CCD) |
| FIDPROPS | L1 pipeline intermediate | FIDPROPS | per fid: `slot`, `id_string`, `id_status` only (slot→fid identity map) |
| OBCENG | L0 telemetry | OBCENG | per sample: `oobagrd3`, `oobagrd6`, `time`, `quality.oobagrd3`, `quality.oobagrd6` |
| ACISENG | L0 telemetry | ACISENG | per sample: `c1bat`, `c1bbt`, `time`, `quality.c1bat`, `quality.c1bbt` |

**ACISENG is *not* a CalDB product.** It is L0 ACIS engineering telemetry
(housekeeping). In the pipeline its values are reduced once by
`asp_get_calib/calc_average_dhtemp()` to a per-observation scalar `dht_mean`
that is written into the FIDPROPS header. In a CalDB-driven reimplementation
this reduction happens in Python (see Stage 0).

Struct references: [acacent.h:48-110](../dstools/asp/aspect_lib/acacent.h#L48-L110),
[fidprops.h:22-78](../dstools/asp/aspect_lib/fidprops.h#L22-L78),
[calalign.h:23-98](../dstools/asp/aspect_lib/calalign.h#L23-L98).

---

## 2. CALALIGN column ↔ FIDPROPS header mapping

`asp_get_calib/load_acis_fidcorr.c` is the canonical bridge between CalDB and
FIDPROPS. It copies (and in one case converts) CalDB values into FIDPROPS
header keywords. A Python reimplementation reads straight from CALALIGN and
needs this mapping to know which CalDB name each FIDPROPS-header formula
refers to:

| Downstream usage (FIDPROPS header) | Source (CALALIGN ACISFIDCORR column) | Notes |
|---|---|---|
| `dht_base` | `dh_temp_base` | **K → °C conversion:** `dht_base = dh_temp_base - 273.15`. [load_acis_fidcorr.c:207](../dstools/asp/asp_get_calib/load_acis_fidcorr.c#L207) |
| `deg_per_cnt` | `deg_per_count` | direct copy |
| `fid_cte` | `fid_CTE` | direct copy (note capitalization) |
| `fid_y_coe` | `fid_y_center_nom` | direct copy (renamed) |
| `fid_z_coe` | `fid_z_center_nom` | direct copy (renamed) |
| `ang_y_nom` (per fid) | `fid_y_ang_nom[fid_num - 1]` | indexed by 1-based fid number parsed from `id_string` |
| `ang_z_nom` (per fid) | `fid_z_ang_nom[fid_num - 1]` | same |

| Downstream usage (FIDPROPS header) | Source (CALALIGN PERIFIDCORR column) | Notes |
|---|---|---|
| `peri_y_grd3` | `peri_y_oobagrd3` | direct copy (renamed) |
| `peri_y_grd6` | `peri_y_oobagrd6` | direct copy (renamed) |
| `peri_z_grd3` | `peri_z_oobagrd3` | direct copy (renamed) |
| `peri_z_grd6` | `peri_z_oobagrd6` | direct copy (renamed) |

Detector selection: CALALIGN ACISFIDCORR and PERIFIDCORR each carry one row
per detector (ACIS-I, ACIS-S, HRC-I, HRC-S). Pick the row whose `detector`
string matches the observation's instrument. Within ACISFIDCORR, the fid
number (1–6) indexes into the `fid_{y,z}_ang_nom[6]` arrays; get it by
parsing the trailing integer from the FIDPROPS `id_string`, e.g.
`"ACIS-5"` → 5. Pattern:
[load_acis_fidcorr.c:192-198](../dstools/asp/asp_get_calib/load_acis_fidcorr.c#L192-L198).

---

## 3. Stage 0 — Calibration load

A CalDB-driven reimplementation performs this work explicitly. The pipeline
does it inside `asp_get_calib` and writes the results into FIDPROPS.

### 3.1 Load CALFDC → polynomial coefficients

Open CALFDC, read one row, copy out `fd_y_fid[20]` and `fd_z_fid[20]`
(only the fid versions are needed; the `_star` coefficients are for guide
stars, out of scope here). Pattern:
[load_field_distortion.c:76-90](../dstools/asp/asp_get_calib/load_field_distortion.c#L76-L90).

### 3.2 Load CALALIGN ACISFIDCORR → thermal constants + nominal angles

Iterate rows; for the row whose `detector` matches the observation's
instrument, copy out every field listed in §1 (dh_temp_base through
fid_z_ang_nom[6]). Pattern:
[load_acis_fidcorr.c:101-116](../dstools/asp/asp_get_calib/load_acis_fidcorr.c#L101-L116).

### 3.3 Load CALALIGN PERIFIDCORR → periscope coefficients

Same file, different extension. For the detector-matched row, copy out the
four `peri_{y,z}_oobagrd{3,6}` values.

### 3.4 Compute `dht_mean` from ACISENG

The per-observation DH-temperature DC baseline. Open ACISENG, stream rows:

```
sumTemp = 0.0
nvals   = 0
for each ACISENG row:
    if row.time < tstart: continue
    if row.time > tstop:  break
    if quality.c1bat == GOOD_NEW:
        sumTemp += row.c1bat; nvals += 1
    if quality.c1bbt == GOOD_NEW:
        sumTemp += row.c1bbt; nvals += 1
dht_mean = (sumTemp / nvals) - 273.15    # K → °C
```

Note: both channels (`c1bat`, `c1bbt`) are accumulated into the *same* sum
with `nvals` counting each good sample separately — the result is the mean
of the concatenated stream, not a pairwise mean of the two channels.
Pattern: [load_acis_fidcorr.c:300-323](../dstools/asp/asp_get_calib/load_acis_fidcorr.c#L300-L323).

### 3.5 Resolve per-fid nominal angles

For each fid identified in this observation (one row of FIDPROPS), parse
`fid_num = int(id_string.split('-')[-1])`, then:

```
ang_y_nom_for_this_fid = fid_y_ang_nom[fid_num - 1]
ang_z_nom_for_this_fid = fid_z_ang_nom[fid_num - 1]
```

The FIDPROPS row is also where the slot-to-fid assignment lives — the
downstream stages loop over ACACENT by `slot`, so we need a
`slot → (fid_id, ang_*_nom)` dictionary built here.

---

## 4. Stage 1 — Pixel → ACA angle (`aca_corr_centr`)

**Purpose:** convert pixel centroids `(cent_i, cent_j)` to ACA angles
`(ang_y, ang_z)` via the optical-distortion polynomial, per (slot, sample).

### Inputs
| Item | Source | Units | Shape |
|---|---|---|---|
| `cent_i`, `cent_j` | ACACENT | pixels | per (slot, time) |
| `temperature` (CCD) | ACADATA, **pre-smoothed** | K or °C auto-detect | per (slot, time) |
| `fd_y_fid[20]`, `fd_z_fid[20]` | CALFDC | poly coeffs | scalar arrays |

### Formula

Degree-3, 20-coefficient polynomial in `(R = cent_i, C = cent_j, T)` with
all cross-terms up to cubic:

```
ang_y = fd_y_fid[ 0]
      + fd_y_fid[ 1]*C + fd_y_fid[ 2]*R + fd_y_fid[ 3]*T
      + fd_y_fid[ 4]*C^2 + fd_y_fid[ 5]*C*R + fd_y_fid[ 6]*C*T
      + fd_y_fid[ 7]*R^2 + fd_y_fid[ 8]*R*T + fd_y_fid[ 9]*T^2
      + fd_y_fid[10]*C^3 + fd_y_fid[11]*R*C^2 + fd_y_fid[12]*T*C^2
      + fd_y_fid[13]*C*R^2 + fd_y_fid[14]*C*R*T + fd_y_fid[15]*C*T^2
      + fd_y_fid[16]*R^3 + fd_y_fid[17]*T*R^2 + fd_y_fid[18]*R*T^2
      + fd_y_fid[19]*T^3
```
(`ang_z` uses `fd_z_fid` with the same expansion.)

If `T > 150`, treat as Kelvin and subtract 273.15 before evaluating the
polynomial.

Source: [field_distortion.c:16-56](../dstools/asp/asputils/field_distortion.c#L16-L56).
Application site: [aca_corr_centr.c:254-266](../dstools/asp/aca_corr_centr/aca_corr_centr.c#L254-L266).

### Color-correction caveat

The pipeline adds a per-slot color correction after the polynomial
([aca_corr_centr.c:276-277](../dstools/asp/aca_corr_centr/aca_corr_centr.c#L276-L277)),
but **`asp_get_calib` initializes `color_corr_i/j` to 0.0 for every slot**
and no other code path populates them
([load_ccd_char.c:73-74](../dstools/asp/asp_get_calib/load_ccd_char.c#L73-L74),
[initialize_outputs.c:206-207](../dstools/asp/asp_get_calib/initialize_outputs.c#L206-L207)).
It is effectively a no-op; there is no CalDB source to read. A Python
reimplementation can skip the step.

### CCD temperature pre-smoothing

The `T` argument in the formula is *not* the raw ACADATA temperature — it
is the raw temperature passed through a Savitzky-Golay pre-smooth inside
`aca_corr_centr`
([aca_corr_centr.c:607](../dstools/asp/aca_corr_centr/aca_corr_centr.c#L607)).
A faithful Python implementation applies the same SG smoothing to ADATA
temperatures before feeding them to `Model_Distortion_P`. A simpler
implementation can use raw temperatures and document the deviation.

### Output
`ang_y`, `ang_z` per (slot, time).

---

## 5. Stage 2 — Savitzky-Golay smoothing (`aca_filter_centr`)

**Purpose:** produce `ang_y_sm`, `ang_z_sm` from `ang_y`, `ang_z` per slot
via iterative-reject Savitzky-Golay filtering.

### Inputs
- ACACENT `ang_y`, `ang_z`, `status`, `time`, plus header `t_int`
  (integration time, used to convert the smoothing time to a sample count).

### Parameters (`aca_filter_centr.par`)

| Parameter | Default | Meaning |
|---|---|---|
| `T_sm_fid` | 40.0 s | smoothing window in time; `box_size = T_sm_fid / t_int`, rounded to odd |
| `order` | 4 | SG polynomial order |
| `reject_sigma` | 3.0 | σ-clip threshold |
| `reject_iter` | 5 | max sigma-clip iterations |

### Bad-sample handling

Samples flagged by `alg == IMG_ALG_BAD_CENTR`
([acacent.h:46](../dstools/asp/aspect_lib/acacent.h#L46)) are replaced with
their nearest-good-neighbor value *before* SG runs so they do not bias the
fit. Outliers rejected by the iterative sigma clip get status-bit flags set:

| Bit | Constant | Meaning |
|---|---|---|
| 0x40 | `SIGMA_REJECT_Y` | y-series outlier |
| 0x80 | `SIGMA_REJECT_Z` | z-series outlier |
| 0x20 | `SIGMA_REJECT_M` | magnitude-series outlier |

Ref: [acacent.h:22-26](../dstools/asp/aspect_lib/acacent.h#L22-L26).

### Output
`ang_y_sm`, `ang_z_sm` per (slot, time); `status` bits updated.

### Python implementation note

`scipy.signal.savgol_filter` inside a hand-written iterative sigma-clip loop
reproduces the behaviour; the bad-sample substitution must precede the SG
call on each iteration.

Source: [ACA_Fidmotion.cc:257-357](../dstools/asp/aca_filter_centr/ACA_Fidmotion.cc#L257-L357).

---

## 6. Stage 3 — Thermal + periscope correction (`aca_corr_fid`)

Four sub-steps. Every sub-step operates in place on `ang_y_sm` / `ang_z_sm`.
Parameters come from
[aca_corr_fid.par](../dstools/asp/aca_corr_fid/aca_corr_fid.par):

| Parameter | Default | Meaning |
|---|---|---|
| `nsmooth` | 5 | half-width of fid-flux median filter (samples) |
| `tsample` | 32.8 s | sample bin width for `dhtemps[]` |
| `windowlen` | 152 | Hanning-window length for OBC gradient smoothing |
| `apply_thermalcorr` | yes | run Stage 3a–3c |
| `apply_pericorr` | yes | run Stage 3d |
| `apply_medianfilt` | no | run optional median filter on OBC gradients (see §6.4) |

### 6.1 Stage 3a — DH temperature series from fid flux

Build a per-sample DH-temperature series `dhtemps[nsamp]` on a uniform time
grid with spacing `tsample`. At each sample time `t_k`:

```
for each fid:
    ndx = ACACENT row index of this fid closest (in time) to t_k
    window = ACACENT.flux[fid, ndx-nsmooth : ndx+nsmooth+1]   # clipped at array edges
    flux_med[fid, k] = median(window)
avg_med_cnts[k] = mean_over_fids(flux_med[:, k])

overall_avg_cnts = mean_over_samples(avg_med_cnts)

dhtemps[k] = dht_mean + (avg_med_cnts[k] - overall_avg_cnts) * deg_per_cnt
```

`dht_mean` and `deg_per_cnt` come from Stage 0. `nsmooth` and `tsample` from
the parameter file.

Source:
[tool_functions.c:434-571](../dstools/asp/aca_corr_fid/tool_functions.c#L434-L571)
(formula proper at
[lines 549-553](../dstools/asp/aca_corr_fid/tool_functions.c#L549-L553)).

### 6.2 Stage 3b — Observed center of expansion

Two per-observation scalars:

```
for each fid:
    y_off[fid] = mean_over_samples(ACACENT.ang_y_sm[fid]) - ang_y_nom[fid]
    z_off[fid] = mean_over_samples(ACACENT.ang_z_sm[fid]) - ang_z_nom[fid]

obs_y_coe = fid_y_coe + mean_over_fids(y_off)
obs_z_coe = fid_z_coe + mean_over_fids(z_off)
```

`ang_{y,z}_nom[fid]` from Stage 0. `fid_{y,z}_coe` = CALALIGN
`fid_{y,z}_center_nom`.

Source:
[tool_functions.c:591-657](../dstools/asp/aca_corr_fid/tool_functions.c#L591-L657).

### 6.3 Stage 3c — Thermal correction, in place

For each ACACENT row:

```
k = step_index(acacent.time, sample_times)          # sample bin containing this time
corr = (dhtemps[k] - dht_base) * fid_cte
ang_y_sm -= corr * (ang_y_sm - obs_y_coe)
ang_z_sm -= corr * (ang_z_sm - obs_z_coe)
```

**Time alignment:** nearest-earlier (step) interpolation — no linear
interp. The pipeline searches for the bin `k` such that
`sample_times[k] <= acacent.time < sample_times[k+1]`, with clamping to the
last bin beyond the end
([tool_functions.c:774-789](../dstools/asp/aca_corr_fid/tool_functions.c#L774-L789)).

`dht_base` (already °C after Stage 0) and `fid_cte` are per-observation
scalars from CALALIGN ACISFIDCORR.

Source:
[tool_functions.c:790-802](../dstools/asp/aca_corr_fid/tool_functions.c#L790-L802).

### 6.4 Stage 3d — Periscope correction, in place

Preprocessing of OBCENG gradients (done once, not per ACACENT row):

1. Stream OBCENG rows with `t ∈ [acacent.tmin, acacent.tmax]`. Replace any
   sample with `quality.oobagrd{3,6} != GOOD_NEW` by a sentinel, then
   interpolate the sentinels away via `verifyOBCUnsmoothGradients`.
2. **Compute `meanGradient3`, `meanGradient6` from the RAW (unsmoothed)
   gradient arrays.** These means are saved and reused at application time.
   ([tool_functions.c:1242-1251](../dstools/asp/aca_corr_fid/tool_functions.c#L1242-L1251))
3. If `apply_medianfilt` is `yes` (off by default), apply
   `apply_median_filter(values, n, num_sigmas=4, half_width=windowlen/4,
   pad_algo=TREND)` to each gradient in place
   ([tool_functions.c:965-1002](../dstools/asp/aca_corr_fid/tool_functions.c#L965-L1002)).
4. Apply `asp_slide_smooth(values, n, HANNING, windowlen)` to each gradient
   in place — always, regardless of `apply_medianfilt`.
   ([tool_functions.c:1009,1020](../dstools/asp/aca_corr_fid/tool_functions.c#L1009-L1020))

The resulting in-place arrays are `smoothGradient3[]`, `smoothGradient6[]`.
The `mean` values were computed *before* step 3/4 and represent the raw
gradient baseline.

Per-ACACENT-row application:

```
ii = step_index(acacent.time, obc_times)      # nearest-earlier; fallback to closer endpoint
d3 = smoothGradient3[ii] - meanGradient3
d6 = smoothGradient6[ii] - meanGradient6
ang_y_sm -= d3 * peri_y_grd3 + d6 * peri_y_grd6
ang_z_sm -= d3 * peri_z_grd3 + d6 * peri_z_grd6
```

`peri_{y,z}_grd{3,6}` from CALALIGN PERIFIDCORR (using the mapping in §2 —
these are called `peri_{y,z}_oobagrd{3,6}` in CALALIGN).

Source:
[tool_functions.c:858-867](../dstools/asp/aca_corr_fid/tool_functions.c#L858-L867)
for application,
[tool_functions.c:815-843](../dstools/asp/aca_corr_fid/tool_functions.c#L815-L843)
for the step-interpolation with endpoint fallback.

---

## 7. Dependency map

```
CALFDC  ─── fd_{y,z}_fid[20] ─────────────┐
                                          ▼
ACACENT.cent_{i,j}, ACADATA.temperature ──► Model_Distortion_P ──► ang_{y,z}
CALALIGN.color_corr (zero, skipped) ───────┘                         │
                                                                     ▼
                                            Savitzky-Golay (order 4, 40 s) ──► ang_{y,z}_sm
                                                                     │
                                                                     ▼
ACISENG.c1bat, c1bbt ─── calc_average_dhtemp ── dht_mean ───┐
CALALIGN.deg_per_count ─────────────────────────────────────┤
ACACENT.flux ── median filter (nsmooth=5) ──────────────────┤
                                                            ▼
                                                        dhtemps[k]
                                                            │
CALALIGN.dh_temp_base ── K→°C ─ dht_base ───┐               │
CALALIGN.fid_CTE ── fid_cte ────────────────┤               │
CALALIGN.fid_{y,z}_center_nom ──┐           │               │
CALALIGN.fid_{y,z}_ang_nom[fid] ┤ centroid  │               │
                                │ nominals  │               │
                                ▼           ▼               ▼
                       obs_{y,z}_coe ── thermal correction (per sample)
                                                            │
                                                            ▼
                                                       ang_{y,z}_sm
                                                            │
OBCENG.oobagrd{3,6} + quality ── (verifyUnsmooth,           │
                                   optional median filter,  │
                                   Hanning windowlen=152) ──┤
                                  meanGradient{3,6} (raw)  ─┤
CALALIGN.peri_{y,z}_oobagrd{3,6} ───────────────────────────┤
                                                            ▼
                                              periscope correction (per sample)
                                                            │
                                                            ▼
                                                   ang_{y,z}_sm (final)
```

---

## 8. Source-file index

| Concern | File |
|---|---|
| ACACENT struct + status bits | [dstools/asp/aspect_lib/acacent.h](../dstools/asp/aspect_lib/acacent.h) |
| FIDPROPS struct (for slot→id map) | [dstools/asp/aspect_lib/fidprops.h](../dstools/asp/aspect_lib/fidprops.h) |
| CALALIGN struct (ACISFIDCORR, PERIFIDCORR members) | [dstools/asp/aspect_lib/calalign.h](../dstools/asp/aspect_lib/calalign.h) |
| ACACAL struct (runtime bundle, where CALFDC/CALALIGN land) | [dstools/asp/aspect_lib/acacal.h](../dstools/asp/aspect_lib/acacal.h) |
| CALFDC → ACACAL copy | [dstools/asp/asp_get_calib/load_field_distortion.c](../dstools/asp/asp_get_calib/load_field_distortion.c) |
| CALALIGN → FIDPROPS copy + `calc_average_dhtemp` | [dstools/asp/asp_get_calib/load_acis_fidcorr.c](../dstools/asp/asp_get_calib/load_acis_fidcorr.c) |
| Color-correction zeroing | [dstools/asp/asp_get_calib/load_ccd_char.c](../dstools/asp/asp_get_calib/load_ccd_char.c), [initialize_outputs.c](../dstools/asp/asp_get_calib/initialize_outputs.c) |
| `Model_Distortion_P` polynomial | [dstools/asp/asputils/field_distortion.c](../dstools/asp/asputils/field_distortion.c) |
| Stage 1 application + CCD-temp pre-smooth | [dstools/asp/aca_corr_centr/aca_corr_centr.c](../dstools/asp/aca_corr_centr/aca_corr_centr.c) |
| Stage 2 SG filter | [dstools/asp/aca_filter_centr/ACA_Fidmotion.cc](../dstools/asp/aca_filter_centr/ACA_Fidmotion.cc) |
| Stage 2 parameters | [dstools/asp/aca_filter_centr/aca_filter_centr.par](../dstools/asp/aca_filter_centr/aca_filter_centr.par) |
| Stage 3 thermal + periscope + OBCENG load | [dstools/asp/aca_corr_fid/tool_functions.c](../dstools/asp/aca_corr_fid/tool_functions.c) |
| Stage 3 parameters | [dstools/asp/aca_corr_fid/aca_corr_fid.par](../dstools/asp/aca_corr_fid/aca_corr_fid.par) |
