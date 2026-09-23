from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table
from kadi.commands import get_starcats
from mica.archive.obspar import get_obspar

# Companion to fid_process_pred.py: starts from pipeline-computed
# centroid angles in ACACENT (post-aca_corr_centr — i.e. background-
# subtracted centroids run through the CALFDC field-distortion
# polynomial) and applies the ACIS base-plate (CALALIGN HDU 3) and
# periscope (CALALIGN HDU 4) corrections that aca_corr_fid applies in
# the production pipeline. Output keeps the full time series (no per-fid
# averaging): one row per (time, fid-slot) ACACENT sample.

obsid = 29878
obspar = get_obspar(obsid)
acas = get_starcats(obsid=obsid)[0]
fids = acas[acas["type"] == "FID"]
fid_slots = np.asarray(fids["slot"], dtype=int)
tstart, tstop = obspar["tstart"], obspar["tstop"]

# Path placeholders — fill in at run time. ACACENT carries cent_i,
# cent_j, ang_y, ang_z (post-distortion), ang_y_sm, ang_z_sm (smoothed
# + post-aca_corr_fid), and flux (counts in image; the temperature
# sensor — see below). We use the unsmoothed ang_y / ang_z so the
# corrections below aren't double-counted. FIDPROPS supplies the
# scalars dht_mean and deg_per_cnt that turn fid flux into a per-time
# DH temperature. G3 / G6 are smoothed OOBAGRD3 / OOBAGRD6 housing-
# gradient thermistors from OBCENG, 1-D arrays aligned to acen["time"].
ACACENT_PATH = "test_data/pcadf<obsid>N<vers>_acen<N>.fits"   # e.g. pcadf052138759N001_acen1.fits
FIDPROPS_PATH = "test_data/pcadf<obsid>N<vers>_fidpr<N>.fits"  # e.g. pcadf052138759N001_fidpr1.fits
G3 = None  # shape (N_rows,)
G6 = None  # shape (N_rows,)

calalign_hdus = fits.open(
    "/home/ascds/DS.release/CALDB/data/chandra/pcad/align/pcadD2021-07-02alignN0010.fits"
)
calalign_3 = Table(calalign_hdus[3].data)  # ACIS base-plate fid correction
calalign_4 = Table(calalign_hdus[4].data)  # Periscope correction

# CALALIGN HDU 3 has one row per detector. For HRC observations the
# loader (asp_get_calib/load_acis_fidcorr.c) leaves the values FNaN and
# aca_corr_fid skips the ACIS step; mirror that here.
det_mask = np.array([d.strip() == obspar["detector"] for d in calalign_3["DETECTOR"]])
acis_corr = calalign_3[det_mask][0] if det_mask.any() else None

acen = Table(fits.open(ACACENT_PATH)[1].data)
mask = (
    (acen["time"] >= tstart)
    & (acen["time"] <= tstop)
    & np.isin(acen["slot"], fid_slots)
)
acen = acen[mask]

# Post-distortion centroid angles, in degrees in the ACA frame. The
# production pipeline goes one step further: a Savitzky-Golay smooth
# produces ang_y_sm / ang_z_sm before aca_corr_fid runs. Starting from
# the unsmoothed angles here is the simplest faithful entry point — we
# trade the smoothing residual for not having to re-implement the SG
# window from scratch.
ang_y = np.asarray(acen["ang_y"], dtype=float).copy()
ang_z = np.asarray(acen["ang_z"], dtype=float).copy()
t = np.asarray(acen["time"], dtype=float)
slot = np.asarray(acen["slot"], dtype=int)

# ACIS base-plate thermal correction. See aca_corr_fid/tool_functions.c:792-802.
# DH_TEMP_BASE in CALALIGN is in K; load_acis_fidcorr.c:207 subtracts 273.15.
if acis_corr is not None:
    fidpr_hdu = fits.open(FIDPROPS_PATH)[1]
    # FITS-header aliases of the C names (Vio drops underscores to fit
    # the 8-char keyword limit): dht_mean -> DHTMEAN, deg_per_cnt -> DEGPCNT.
    dht_mean = float(fidpr_hdu.header["DHTMEAN"])
    deg_per_cnt = float(fidpr_hdu.header["DEGPCNT"])

    # Per-time DH temp from fid counts (the fids are the temperature
    # sensor). See aca_corr_fid/tool_functions.c:548-553. The on-disk
    # column is "counts" — the C struct calls it "flux" but the Vio
    # template renames it. Pipeline does a per-fid window-median then
    # averages across fids; this compact form takes the mean counts
    # across fid slots per timestamp.
    counts = np.asarray(acen["counts"], dtype=float)
    _, inv = np.unique(t, return_inverse=True)
    cnt_per_t = np.bincount(inv, weights=counts) / np.bincount(inv)
    dh_temp_C = (dht_mean + (cnt_per_t - cnt_per_t.mean()) * deg_per_cnt)[inv]

    dh_base_C = float(acis_corr["DH_TEMP_BASE"]) - 273.15
    fid_cte = float(acis_corr["FID_CTE"])
    fid_y_center = float(acis_corr["FID_Y_CENTER_NOM"])
    fid_z_center = float(acis_corr["FID_Z_CENTER_NOM"])
    corr = (dh_temp_C - dh_base_C) * fid_cte
    ang_y -= corr * (ang_y - fid_y_center)
    ang_z -= corr * (ang_z - fid_z_center)

# Periscope correction. See aca_corr_fid/tool_functions.c:858-867.
# (g - mean_g) integrates to ~0 over the obs window, so this term has
# almost no effect on per-fid averages — but it does shape the time
# series, which is why it stays per-row here.
g3_dev = G3 - G3.mean()
g6_dev = G6 - G6.mean()
ang_y -= g3_dev * float(calalign_4["PERI_Y_OOBAGRD3"][0])
ang_y -= g6_dev * float(calalign_4["PERI_Y_OOBAGRD6"][0])
ang_z -= g3_dev * float(calalign_4["PERI_Z_OOBAGRD3"][0])
ang_z -= g6_dev * float(calalign_4["PERI_Z_OOBAGRD6"][0])

# Direction cosines in the ACA frame, exact tan form (no small-angle
# approximation). Shape (N_rows, 3); [:, 0] is identically 1.
# See docs/PREDICTED_FID_AND_ASP_SOLVE.md:268-273.
d_aca_meas = np.stack(
    [
        np.ones_like(ang_y),
        np.tan(np.radians(ang_y)),
        np.tan(np.radians(ang_z)),
    ],
    axis=-1,
)

meas = Table(
    {
        "time": t,
        "slot": slot,
        "ang_y_corr": ang_y,
        "ang_z_corr": ang_z,
        "d_aca_meas_x": d_aca_meas[:, 0],
        "d_aca_meas_y": d_aca_meas[:, 1],
        "d_aca_meas_z": d_aca_meas[:, 2],
    }
)

# Composing the two scripts into the LM fit (not executed, kept here as
# a pointer). For a single time slice, group meas by slot and align with
# fid_process_pred.py's d_aca_est:
#
# from scipy.optimize import least_squares
# def residuals(params, p_lsi, lsi0_stt, stt0_stf, M, rrc0_fc_x, d_aca_meas):
#     theta_X, dy, dz = params
#     ...  # build d_aca_est for this (theta_X, dy, dz)
#     return (d_aca_meas[:, 1:] - d_aca_est[:, 1:]).ravel()
# least_squares(residuals, x0=[0.0, 0.0, 0.0],
#               args=(p_lsi, lsi0_stt, stt0_stf, M, rrc0_fc_x, d_aca_meas))
