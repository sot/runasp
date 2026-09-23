# Understanding `_asol1` and `_acen1` pipeline products

## Context

Building a shared mental model of the aspect pipeline before drafting a
per-step pipeline report. This document answers two questions the user
raised about specific files they have in hand:

1. What each column means in `pcadf893066380N001_acen1.fits.gz` (`cent_i`, `cent_j`, `ang_y`, `ang_z`, `ang_y_sm`, `ang_z_sm`).
2. What processing has been applied to `ang_y` / `ang_z` (not just to the `_sm` versions).

Everything below is grounded in the pipeline shell scripts in
[pipelines/asp/](../pipelines/asp/) and the tool
source under [dstools/asp/](../dstools/asp/).

---

## Part 1 — The `_acen1` file (ACACENT)

`_acen1` is the ACACENT product — the per-slot, per-sample ACA image
measurements. A single file is opened read/write by successive pipeline
stages; each stage adds or rewrites columns. The filename pattern
`pcadf<MET>N001_acen1.fits.gz` marks it as a **CAI-level** product.

### Column reference

| Column | Units | What it is |
|---|---|---|
| `cent_i` | pixels | CCD **row** centroid (image I coordinate) measured from the raw ACA image window by [aca_calc_centr](../dstools/asp/aca_calc_centr). |
| `cent_j` | pixels | CCD **column** centroid (image J coordinate), same source. |
| `ang_y` | arcmin (ACA Y) | Centroid converted to ACA-frame Y angle. See "corrections applied" below. |
| `ang_z` | arcmin (ACA Z) | Same for ACA-frame Z angle. |
| `ang_y_sm` | arcmin | Savitzky-Golay-smoothed `ang_y`, with thermal + periscope corrections applied in place (fids only). |
| `ang_z_sm` | arcmin | Same for Z. |

Struct definition: [dstools/asp/aspect_lib/acacent.h](../dstools/asp/aspect_lib/acacent.h).

### Pipeline stages that write this file (order matters)

All stages below operate in place on `$outdir/pcad$root_acen1.fits`
— the file after the last stage is what you find on disk:

```
[asp_l1_std.ped]
 246  aca_calc_centr       → creates acen1, populates cent_i, cent_j, flux, time, slot, alg
 277  aca_corr_centr       → populates ang_y, ang_z
 288  aca_filter_centr     → populates ang_y_sm, ang_z_sm (pure SG smoothing at this point)
 299  aca_corr_fid         → rewrites ang_y_sm, ang_z_sm in place (thermal + periscope)
```

(`asp_corr_props`, `asp_forward_kalman`, `asp_solve`, `asp_make_qualint` and
`asp_vv_centroids` all **read** the acen1 file afterward; they don't
modify the six columns above.)

So the delivered `_acen1.fits.gz` you have on disk is post-Stage 3.

### What has been applied to `ang_y` / `ang_z` (the unsmoothed columns)

Written by [aca_corr_centr.c](../dstools/asp/aca_corr_centr/aca_corr_centr.c).
**`ang_y` / `ang_z` are essentially the pixel → ACA-angle conversion —
the output of the field-distortion polynomial, nothing more.** No
smoothing, no thermal correction, no periscope correction.

The per-record code path for every (slot, sample):

1. **Field-distortion polynomial** — the conversion itself. A 20-term
   cubic polynomial in `(cent_j, cent_i, T_ccd)` from CALFDC, evaluated
   by `Model_Distortion_P` ([asputils/field_distortion.c:16-56](../dstools/asp/asputils/field_distortion.c#L16-L56)),
   applied at [aca_corr_centr.c:243-266](../dstools/asp/aca_corr_centr/aca_corr_centr.c#L243-L266).
   The CCD temperature `T` fed to the polynomial is SG-pre-smoothed
   ([aca_corr_centr.c:607](../dstools/asp/aca_corr_centr/aca_corr_centr.c#L607))
   — a low-noise temperature goes into the polynomial, but the
   centroids themselves are not smoothed.

2. Two tiny in-memory pre-polynomial tweaks — **not corrections to
   `ang_y` / `ang_z` in any persistent sense**:
   - CTI shift of `cent_i`/`cent_j` ([aca_corr_centr.c:221-222](../dstools/asp/aca_corr_centr/aca_corr_centr.c#L221-L222)).
   - FM-algorithm shift of `cent_i`/`cent_j`, only when the row's
     algorithm is FM ([aca_corr_centr.c:226-232](../dstools/asp/aca_corr_centr/aca_corr_centr.c#L226-L232)).

   Both modify the in-memory struct fields that then feed
   `Model_Distortion_P`. The `cent_i`/`cent_j` columns are opened
   **READONLY** ([aca_corr_centr.c:95-96](../dstools/asp/aca_corr_centr/aca_corr_centr.c#L95-L96))
   so these tweaks are not persisted — the pixel columns on disk remain
   the raw values from `aca_calc_centr`. They nudge the polynomial
   input by small amounts; they are not a separate correction stage on
   the output angles.

3. **Color correction** — *nominally* a post-polynomial angle adjustment
   ([aca_corr_centr.c:276-277](../dstools/asp/aca_corr_centr/aca_corr_centr.c#L276-L277)),
   but **effectively a no-op**: `asp_get_calib` initializes `color_corr_i/j`
   to 0.0 and no code path populates them ([load_ccd_char.c:73-74](../dstools/asp/asp_get_calib/load_ccd_char.c#L73-L74),
   [initialize_outputs.c:206-207](../dstools/asp/asp_get_calib/initialize_outputs.c#L206-L207)).

Net summary for `ang_y` / `ang_z`:
- Distortion polynomial output.
- No smoothing.
- No thermal correction. No periscope correction. No fid-light
  correction of any kind.

The corrections the user usually cares about — thermal and periscope —
live only on `ang_y_sm` / `ang_z_sm`.

### What has been applied to `ang_y_sm` / `ang_z_sm`

Two tools write these columns, the second on top of the first:

1. **[aca_filter_centr](../dstools/asp/aca_filter_centr/ACA_Fidmotion.cc)** — Savitzky-Golay smoothing (order 4, 40 s window by default) with iterative 5σ rejection. Bad-status samples are substituted with their nearest-good-neighbor value before each SG iteration so they don't bias the fit. Sigma-rejected samples get `SIGMA_REJECT_{Y,Z,M}` bits set in the `status` column ([acacent.h:22-26](../dstools/asp/aspect_lib/acacent.h#L22-L26)).

2. **[aca_corr_fid](../dstools/asp/aca_corr_fid/tool_functions.c)** — for fid-light slots only, applies two further corrections *in place* on the SG-smoothed values:
   - **Thermal correction** from the fid-flux-derived DH-temperature series ([tool_functions.c:790-802](../dstools/asp/aca_corr_fid/tool_functions.c#L790-L802)).
   - **Periscope correction** from OBCENG `oobagrd3`/`oobagrd6` gradients ([tool_functions.c:858-867](../dstools/asp/aca_corr_fid/tool_functions.c#L858-L867)).

Details of the thermal+periscope math are already captured in
[FID_CORRECTION_STAGES.md](FID_CORRECTION_STAGES.md) §5–6.

Guide-star slots are smoothed but not touched by `aca_corr_fid` —
thermal/periscope corrections are fid-specific.

### Confirming the correction state from the file header

There are **no dedicated correction-status keywords** (no `FID_CORR`,
`THERMAL_CORR`, `PERI_CORR`, `SM_APPLIED`). The canonical audit trail is
the `HISTORY` keyword chain appended by each tool's
`FW_Application::setupOutput`:

- `HISTORY` entries naming `aca_calc_centr`, `aca_corr_centr`, `aca_filter_centr`, `aca_corr_fid` — presence of all four indicates a fully-corrected acen1.
- `aca_corr_fid` early-exits (no history entry, no correction) if both `apply_thermalcorr=no` and `apply_pericorr=no` ([aca_corr_fid.c:63-64](../dstools/asp/aca_corr_fid/aca_corr_fid.c#L63-L64)).
- `n_img`, `ang_y_nea`, `ang_z_nea` per-slot header keywords are written by `aca_filter_centr` ([ACA_Fidmotion.cc:383-385](../dstools/asp/aca_filter_centr/ACA_Fidmotion.cc#L383-L385)) — their presence confirms Stage 2 ran.

---

## Plain-English reference — what each correction does

Quick narrative descriptions of the individual corrections the
pipeline applies. Useful both for orienting readers of the report and
for sanity-checking a Python reimplementation against the *intent* of
each step, not just its math.

### Field-distortion polynomial (the pixel → ACA-angle conversion)

Not a correction in the usual sense — it's the actual conversion from
*where on the CCD the star's image landed* (pixels) to *where in the
sky, relative to ACA boresight, the star is* (arcminutes in the ACA
Y/Z frame). The ACA is a small, fast optical system with real-world
imperfections: the field is not perfectly flat, and a star near the
edge of the ACA's view has its image shifted from where a perfect lens
would put it. A 20-term cubic polynomial in (row, column, CCD
temperature) reproduces the measured distortion map from ground and
in-flight calibration. Temperature enters because thermal expansion
slightly changes the optics geometry.

### CTI correction (pixel-level)

The ACA CCD reads out by shifting accumulated charge across the chip,
row by row, into a readout register. Every shift loses a tiny fraction
of the electrons — "charge-transfer inefficiency." The lost charge
trails the star image, biasing the centroid toward the readout
direction. The CTI correction is a small per-slot nudge that pushes
the pixel centroid back to where it would be if the readout were
lossless. It's slot-dependent because each slot reads from a different
place on the chip, so the total transfer path differs.

### FM (first-moment) algorithm correction (pixel-level)

FM is the simple "center of light" centroiding algorithm — fast and
robust but known to be biased. Because the star's brightest pixel
dominates the weighted sum, the computed centroid gets pulled toward
that pixel's center whenever the true peak isn't exactly centered
there. The effect looks like centroids "snapping" toward pixel
centers, a sub-pixel systematic of ~0.1 pixel. The FM correction is a
pre-measured offset that counters this snap. It runs only on rows that
used the FM algorithm; the PSF-fit and Gaussian algorithms don't have
this bias and skip it entirely.

### Color correction (angle-level, de facto no-op)

The ACA is a small telescope with a lens, and like any lens it focuses
different wavelengths slightly differently — chromatic aberration.
Stars of different colors therefore land at slightly different CCD
positions even when they sit at the same sky location. The color
correction would use each star's B-V color index (from AGASC) to
subtract this chromatic offset from `ang_y` / `ang_z`. The code path
is in place ([aca_corr_centr.c:276-277](../dstools/asp/aca_corr_centr/aca_corr_centr.c#L276-L277))
but the coefficients in CalDB are zero and nothing ever populates them
— so the correction runs on every sample but always adds nothing.
Either the effect was measured to be below the aspect error budget or
a color calibration was never produced. Fid lights would get zero
anyway (they're monochromatic LEDs, not "colored" in the stellar
sense).

### Thermal correction (fid lights, on `ang_y_sm` / `ang_z_sm`)

The ACIS detector housing warms and cools during an observation as the
instrument cycles and as the spacecraft's attitude exposes different
surfaces to sunlight. The fid lights are small LEDs mounted on the
detector housing itself, so when the housing expands the fids
*physically move* — radially outward from a center of thermal
expansion. From the ACA's point of view this looks like a coherent
drift of the fids across the field; uncorrected, the aspect solver
mistakes it for a pointing error.

The correction needs a per-sample temperature. Rather than rely on a
single thermistor, the pipeline uses the **fid lights themselves as
thermometers**: LEDs get slightly brighter or dimmer with temperature,
and that dependence has been calibrated (`deg_per_count` °C per ADU).
So the observed fid flux is converted to a "DH-temperature" series,
anchored to a baseline from an ACIS housekeeping thermistor (`c1bat`,
`c1bbt`), and then fed into the radial-expansion formula. The
resulting shift is subtracted from each fid's smoothed `ang_y_sm`,
`ang_z_sm`.

### Periscope correction (fid lights, on `ang_y_sm` / `ang_z_sm`)

The fid lights sit at the bottom of the optical bench near the
detectors; the ACA sits at the top near the mirrors. They can't see
each other directly — a small folding relay called the **fid-light
periscope** carries the fid images up to the ACA. Like everything
else on the spacecraft, the periscope flexes when its temperature
changes. A flex shifts the apparent position of *every* fid in the
ACA image by the same amount, even though no fid has actually moved.
Uncorrected, this coherent shift gets absorbed into the aspect
solution as a spurious pointing drift.

The fix exploits an empirical correlation discovered after launch:
two on-board-computer housekeeping channels — `oobagrd3` and
`oobagrd6`, optical-bench-assembly thermal-gradient readouts — track
the periscope's flex state closely. The correction is a linear
combination of deviations of these two gradients from their
observation-mean values, with coefficients from CalDB
(`peri_{y,z}_oobagrd{3,6}`). Subtracting that combination from each
fid's `ang_y_sm` / `ang_z_sm` cancels the spurious drift. This
correction was added to the pipeline in 2008 once the bench
temperature / fid-position correlation was established; before then
it didn't exist.

---

## Summary

| Column | Applied by the time you open the file |
|---|---|
| `cent_i`, `cent_j` | Raw pixel centroids from `aca_calc_centr` — untouched by later stages. |
| `ang_y`, `ang_z` | Pixel → ACA-angle via field-distortion polynomial. No smoothing, no thermal, no periscope. |
| `ang_y_sm`, `ang_z_sm` | SG-smoothed (40 s, order 4, 5σ reject). **For fids**: thermal + periscope added in place. **For guide stars**: just smoothing. |

| File | Producer | Scope | `CONTENT` |
|---|---|---|---|
| `pcadf<MET>N001_asol1.fits.gz` | `asp_solve` | single aspect interval (CAI) | `ASPSOL` |
| `pcadf<OBSID>_<OBI>N001_asol1.fits.gz` | `asp_combine` | full observation (OBI) | `ASPSOLOBI` |

---

## Part 2 — Reproducing the fid-light corrections in Python

### Goal

Take `ang_y`, `ang_z` (pixel → ACA-angle conversion, nothing else
applied) from an `_acen1.fits.gz` file and produce `ang_y_corr`,
`ang_z_corr` by applying the thermal and periscope corrections
**directly**, without the Savitzky-Golay smoothing step. For
fid-light slots only (guide stars never get these corrections in the
pipeline).

Full per-stage formulas and source-code citations are in
[FID_CORRECTION_STAGES.md](FID_CORRECTION_STAGES.md)
§5–6. This section is a compact Python recipe; consult that document
for the exact code paths when an edge case bites.

### Files to open

| File | Class | What to read |
|---|---|---|
| `_acen1.fits.gz` (ACACENT) | L1 intermediate | per row: `time`, `slot`, `ang_y`, `ang_z`, `flux`, `status`, `alg` |
| `_fidpr1.fits.gz` (FIDPROPS) | L1 intermediate | per row: `slot`, `id_string`, `id_status`. Header: `dht_mean` (shortcut) |
| CALALIGN (CalDB) | CalDB | extensions `ACISFIDCORR` and `PERIFIDCORR`, pick the row whose `detector` matches `INSTRUME` |
| `_obc*.fits.gz` (OBCENG) | L0 telemetry | per sample: `time`, `oobagrd3`, `oobagrd6`, `quality.oobagrd3`, `quality.oobagrd6` |
| `_acis*.fits.gz` (ACISENG) *(optional — can skip if using FIDPROPS header shortcut)* | L0 telemetry | per sample: `time`, `c1bat`, `c1bbt`, `quality.c1bat`, `quality.c1bbt` |

Everything opens with `astropy.io.fits`. CALALIGN extensions can be
browsed with `fits.info(path)`.

### Stage A — Calibration pull

From **CALALIGN[ACISFIDCORR]**, the row where `detector == INSTRUME`:

| Extract | As | Note |
|---|---|---|
| `dh_temp_base` (K) | `dht_base` (°C) | subtract 273.15 |
| `deg_per_count` | `deg_per_cnt` | direct |
| `fid_CTE` | `fid_cte` | direct |
| `fid_y_center_nom`, `fid_z_center_nom` | `fid_y_coe`, `fid_z_coe` | direct |
| `fid_y_ang_nom[6]`, `fid_z_ang_nom[6]` | per-fid nominal angles | indexed 1..6 |

From **CALALIGN[PERIFIDCORR]**, same detector row:
`peri_y_oobagrd3`, `peri_y_oobagrd6`, `peri_z_oobagrd3`, `peri_z_oobagrd6`
→ `peri_y_grd3`, `peri_y_grd6`, `peri_z_grd3`, `peri_z_grd6`.

### Stage B — Slot → fid identity

Open FIDPROPS. For each row where the fid is actually used:
```
fid_num   = int(row.id_string.split("-")[-1])    # e.g. "ACIS-5" → 5
ang_y_nom = fid_y_ang_nom[fid_num - 1]
ang_z_nom = fid_z_ang_nom[fid_num - 1]
slot_to_fid[row.slot] = (fid_num, ang_y_nom, ang_z_nom)
```

### Stage C — `dht_mean` (scalar, °C)

Two options:
- **Shortcut**: read `dht_mean` from the FIDPROPS header. Already in °C.
- **Self-contained**: stream ACISENG rows with `t ∈ [TSTART, TSTOP]`,
  accumulate every `c1bat` and `c1bbt` whose quality is `GOOD_NEW`
  (both channels go into the *same* sum), divide by count, subtract
  273.15. Formula in [FID_CORRECTION_STAGES.md §3.4](FID_CORRECTION_STAGES.md).

### Stage D — DH-temperature series from fid flux

Uniform time grid `sample_times[k]` over the ACACENT range with
spacing `tsample = 32.8 s` (parameter default). For each `k`:
```
flux_med = []
for each fid slot s:
    rows_s = ACACENT rows with slot == s
    ndx    = argmin(|rows_s.time - sample_times[k]|)
    window = rows_s.flux[max(0, ndx-5) : ndx+6]      # nsmooth = 5
    flux_med.append(median(window))
avg_med_cnts[k] = mean(flux_med)

overall_avg = mean(avg_med_cnts)
dhtemps[k]  = dht_mean + (avg_med_cnts[k] - overall_avg) * deg_per_cnt
```

### Stage E — Observed center of expansion (per-observation scalars)

```
y_off, z_off = [], []
for fid slot s in slot_to_fid:
    rows = ACACENT rows with slot == s and status == GOOD
    y_off.append(mean(rows.ang_y) - ang_y_nom[s])
    z_off.append(mean(rows.ang_z) - ang_z_nom[s])

obs_y_coe = fid_y_coe + mean(y_off)
obs_z_coe = fid_z_coe + mean(z_off)
```

The pipeline uses `ang_y_sm` / `ang_z_sm` here; with smoothing skipped,
using `ang_y` / `ang_z` gives essentially the same per-observation mean
(SG smoothing removes zero-mean noise).

### Stage F — Thermal correction (per ACACENT row, fid slots only)

For each fid row, find the bin `k` such that
`sample_times[k] <= row.time < sample_times[k+1]` (nearest-earlier,
clamp at the last bin):
```
corr       = (dhtemps[k] - dht_base) * fid_cte
ang_y_tc   = row.ang_y - corr * (row.ang_y - obs_y_coe)
ang_z_tc   = row.ang_z - corr * (row.ang_z - obs_z_coe)
```

### Stage G — Prepare OBCENG gradients (once per observation)

1. Load OBCENG rows with `t ∈ [ACACENT.tmin, ACACENT.tmax]`.
2. Replace any sample with bad quality (`quality.oobagrd{3,6} != GOOD_NEW`)
   by linear interpolation over good neighbors.
3. **Compute `meanGradient3`, `meanGradient6` from these RAW gradients.**
   Save — these are used as the baseline at correction time.
4. *(optional, default off)* Apply a 4σ median filter with
   half-width `windowlen/4 = 38`.
5. Slide-smooth each gradient with a **Hanning window of length 152**.
   `np.convolve(grad, hann(152) / hann(152).sum(), mode="same")` is a
   reasonable first cut; the pipeline's `asp_slide_smooth` has careful
   edge handling — see [FID_CORRECTION_STAGES.md §6.4](FID_CORRECTION_STAGES.md).

The outputs of step 5 are `smoothGradient3[]`, `smoothGradient6[]`.
Note that the *means* were frozen before smoothing.

### Stage H — Periscope correction (per ACACENT row, fid slots only)

For each fid row, find `ii` in the OBCENG time grid by nearest-earlier
(with endpoint fallback to the closer end):
```
d3         = smoothGradient3[ii] - meanGradient3
d6         = smoothGradient6[ii] - meanGradient6
ang_y_corr = ang_y_tc - (d3 * peri_y_grd3 + d6 * peri_y_grd6)
ang_z_corr = ang_z_tc - (d3 * peri_z_grd3 + d6 * peri_z_grd6)
```

---

## Verification

1. **File-level sanity checks** (answers the Part-1/Part-2 questions):
   - `fdump pcadf893066380N001_acen1.fits.gz+1 STDOUT HISTORY -` — look
     for all four tool names (`aca_calc_centr`, `aca_corr_centr`,
     `aca_filter_centr`, `aca_corr_fid`) in the HISTORY chain.
   - Compare `CONTENT` between the two asol1 files — `ASPSOL` vs
     `ASPSOLOBI`.
   - The CAI asol1 should still have `ADY`/`ADZ`/`ADTHETA`/`*_RAW`/
     `Q_ATT_RAW`; the OBI asol1 should not.

2. **Python correction vs the pipeline's final product.** The file's
   `ang_y_sm` already has thermal + periscope applied. Apply the same
   SG-smoothing (order 4, 40 s, 5σ reject) to your `ang_y_corr` and
   check it matches `ang_y_sm` within floating-point noise. If it
   doesn't, the difference tells you which stage is off.

3. **Pipeline baseline with corrections disabled.** Re-run `aca_corr_fid`
   with `apply_thermalcorr=no apply_pericorr=no`; the resulting
   `ang_y_sm` is just SG-smoothed `ang_y` with no fid correction. The
   difference between this baseline `ang_y_sm` and the normal
   pipeline's `ang_y_sm` is exactly the thermal+periscope term — your
   Python correction (after smoothing) should reproduce that
   difference.
