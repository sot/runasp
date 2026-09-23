import os
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table
from cheta import fetch
from cxotime import CxoTime
from kadi.commands import get_starcats
from mica.archive.asp_l1 import get_files
from mica.archive.obspar import get_obspar
from scipy.signal.windows import hann
from scipy.signal import savgol_filter
from scipy.spatial.transform import Rotation
from scipy.interpolate import interp1d


CALDB_DIR = Path(os.environ.get("CALDB", "/home/ascds/DS.release/CALDB"))


# centroid algorithms

FM = 1           # First moment
PSF = 2          # PSF fit
GAUSS = 4        # Circular Gaussian fit
EGAUSS = 8       # Elliptical Gaussian
FLUX_SMOOTH = 16 # Smoothing flux


@dataclass
class CalAlign:
    """CALDB align-file rows the predicted-fid chain needs.

    Hand-mapped from the four HDUs of `pcadD<date>alignN<vers>.fits`:

    - `align`: ext 1, single Row matched on `INSTR_ID == obspar["detector"]`.
      Carries `LSI0_STT`, `RRC0_FC_X`, `ACA_SC_ALIGN`, `ACA_MISALIGN`,
      `FTS_MISALIGN`.
    - `fids`: ext 2, full Table — fid positions and 5-coefficient
      polynomial corrections in `FID_Y_CORR` / `FID_Z_CORR`.
    - `acis_corr`: ext 3, full Table — ACIS base-plate thermal correction
      parameters (one row per detector).
    - `peri_corr`: ext 4, full Table — periscope correction coefficients
      (single-row in current CALDB).
    """

    align: object
    fids: object
    acis_corr: object
    peri_corr: object

def get_caldb_align_files():
    """Assemble a table of all calalign files in the CALDB, with their start times.

    The table is sorted by start time.
    """
    align_files = []
    for align_file in CALDB_DIR.glob("data/chandra/pcad/align/*.fits"):
        with fits.open(align_file) as calalign_hdus:
            align_files.append({
                "tstart": CxoTime(calalign_hdus[0].header["CVSD0001"]).date,
                "path": align_file
            })
    align_files = Table(align_files)
    align_files.sort(["tstart"])
    return align_files


def get_caldb_align_file(date):
    """
    Get the CALDB align file for a given date.

    Parameters
    ----------
    date : CxoTime-like
        The date for which to retrieve the align file.

    Returns
    -------
    str
        The path to the align file.
    """
    align_files = get_caldb_align_files()
    mask = (align_files["tstart"] <= date)
    if np.sum(mask) == 0:
        raise ValueError(f"Date {date} is before the earliest align file")
    return align_files[mask][-1]["path"]


def get_caldb_calalign(obspar):
    """Get the CALDB calalign tables for a given observation.

    Parameters
    ----------
    obspar : dict
        The observation parameters, as returned by get_obspar.

    Returns
    -------
    CalAlign
        The four extensions of the align file, with HDU 1 already row-selected
        on the obspar detector.
    """
    align_file = get_caldb_align_file(obspar["date_obs"])
    with fits.open(align_file) as calalign_hdus:
        calalign_1 = Table(calalign_hdus[1].data)
        calalign_1["INSTR_ID"] = [instr.strip() for instr in calalign_1["INSTR_ID"]]
        idx = np.argwhere(calalign_1["INSTR_ID"] == obspar["detector"])[0][0]
        align = calalign_1[idx]

        fids = Table(calalign_hdus[2].data)
        fids["FID_SI"] = [si.strip() for si in fids["FID_SI"]]

        acis_corr = Table(calalign_hdus[3].data)
        acis_corr["DETECTOR"] = [det.strip() for det in acis_corr["DETECTOR"]]

        peri_corr = Table(calalign_hdus[4].data)

    return CalAlign(align=align, fids=fids, acis_corr=acis_corr, peri_corr=peri_corr)


def process_centroids(*, centroids, dh_temp_calib, fid_props, obspar, calalign, as_table=True):
    """Process the fiducial centroids to apply thermal and periscope corrections,
    and compute direction cosines in the ACA frame.

    Simplifications vs the C `aca_corr_fid` pipeline (kept on purpose —
    each is a deliberate choice in this Python re-implementation):

    - DH-temperature is built from the arithmetic mean of `counts`
      across fid slots per timestamp. The C pipeline does a per-fid
      windowed median first, then averages across fids; with cosmic-ray
      hits the median vs mean differs by tens of counts.

    Centroid SG smoothing (matches `aca_filter_centr`) and OOBAGRD3/6
    Hann smoothing (matches `aca_corr_fid`) are applied here via
    `fid_smoothing`. The thermal-correction center of expansion follows
    `tool_functions.c:615-643`: `obs_y/z_coe = FID_Y/Z_CENTER_NOM + δ`,
    where δ is the mean over good ACIS fids of
    (time_mean(ang_y_sm) − ang_y_nom). `ang_y_nom` / `ang_z_nom` and the
    `id_status == "GOOD"` filter come from FIDPROPS ext-1.
    """
    t = np.asarray(centroids["time"], dtype=float)
    slot = np.asarray(centroids["slot"], dtype=int)

    # smoothing
    ang_y, ang_z = sg_smooth_per_slot(
        time=t,
        slot=slot,
        ang_y=centroids["ang_y"],
        ang_z=centroids["ang_z"],
    )

    times, idx = np.unique(t, return_inverse=True)

    # ACIS base-plate thermal correction.
    if np.isin(obspar["detector"], calalign.acis_corr["DETECTOR"]):
        det_mask = (obspar["detector"] == calalign.acis_corr["DETECTOR"])
        acis_corr = calalign.acis_corr[det_mask][0]
        dht_mean = dh_temp_calib["dht_mean"]
        deg_per_cnt = dh_temp_calib["deg_per_cnt"]

        # Per-time DH temp from fid counts; see the third bullet in the
        # docstring and aca_corr_fid/tool_functions.c:548-553.
        counts = np.asarray(centroids["counts"], dtype=float)
        counts_per_t = np.bincount(idx, weights=counts) / np.bincount(idx)
        dh_temp_C = (dht_mean + (counts_per_t - counts_per_t.mean()) * deg_per_cnt)[idx]

        dh_base_C = float(acis_corr["DH_TEMP_BASE"]) - 273.15  # in °C
        fid_cte = float(acis_corr["FID_CTE"])

        # Per-obs center of expansion (tool_functions.c:615-643).
        # δ = mean over good ACIS fids of (time_mean(ang_y_sm) − ang_y_nom).
        good = np.asarray(
            [str(s).strip() == "GOOD" for s in fid_props["id_status"]]
        )
        good_props = fid_props[good]
        offsets_y = []
        offsets_z = []
        for row in good_props:
            sel = slot == int(row["slot"])
            if not np.any(sel):
                continue
            offsets_y.append(ang_y[sel].mean() - float(row["ang_y_nom"]))
            offsets_z.append(ang_z[sel].mean() - float(row["ang_z_nom"]))
        delta_y = float(np.mean(offsets_y))
        delta_z = float(np.mean(offsets_z))
        fid_y_center = float(acis_corr["FID_Y_CENTER_NOM"]) + delta_y
        fid_z_center = float(acis_corr["FID_Z_CENTER_NOM"]) + delta_z

        corr = (dh_temp_C - dh_base_C) * fid_cte
        ang_y -= corr * (ang_y - fid_y_center)
        ang_z -= corr * (ang_z - fid_z_center)

    # Periscope correction. ±1000s pad on the cheta fetch so the
    # interpolation has neighbours at the edges of the obs window.
    msids = ["OOBAGRD3", "OOBAGRD6"]
    telem = fetch.Msidset(msids, obspar["tstart"] - 1000, obspar["tstop"] + 1000)
    telem.interpolate(times=times)

    oobagrd3 = hann_smooth(telem["OOBAGRD3"].vals)
    oobagrd6 = hann_smooth(telem["OOBAGRD6"].vals)

    g3_dev = oobagrd3 - oobagrd3.mean()
    g6_dev = oobagrd6 - oobagrd6.mean()
    dy = (g3_dev * float(calalign.peri_corr["PERI_Y_OOBAGRD3"][0])
          + g6_dev * float(calalign.peri_corr["PERI_Y_OOBAGRD6"][0]))
    dz = (g3_dev * float(calalign.peri_corr["PERI_Z_OOBAGRD3"][0])
          + g6_dev * float(calalign.peri_corr["PERI_Z_OOBAGRD6"][0]))
    ang_y -= dy[idx]
    ang_z -= dz[idx]

    # Direction cosines in the ACA frame
    d_aca_meas = np.stack(
        [
            np.ones_like(ang_y),
            np.tan(np.radians(ang_y)),
            np.tan(np.radians(ang_z)),
        ],
        axis=-1,
    )
    d_aca_meas /= np.linalg.norm(d_aca_meas, axis=1, keepdims=True)

    # this is what gets written in centroids file?
    if as_table:
        meas = Table(
            {
                "time": t,
                "slot": slot,
                "ang_y": ang_y * 3600,
                "ang_z": ang_z * 3600,
                "d_aca": d_aca_meas,
            }
        )
        return meas
    return d_aca_meas


def nominal_fid_positions(*, obspar, calalign, fids, as_table=True):
    """Compute the nominal fiducial positions in the LSI frame, including the polynomial correction.
    """
    obs_si = obspar["detector"]

    # STT origin in STF (obspar['sim_z'] in mm)
    # the following should be the same as the SIM_X/Y/Z values in the sim$root_coor0a.fits file (need to check)
    stt0_stf = np.array([obspar['sim_x'], obspar['sim_y'], obspar['sim_z']]) #  (?)

    # LSI origin in STT frame (in mm)
    lsi0_stt = calalign.align["LSI0_STT"]

    simz_offset = lsi0_stt[2] + stt0_stf[2]  # mm

    # select fid light rows.
    # OCAT fid id is absolute (1-14); CALALIGN.fid_num_si is per-instrument
    # relative (1-6 ACIS, 1-4 HRC). See Aca_Id_Image.cc:2232-2233.
    fid_rows = []
    for fid in fids:
        fid_num = int(fid["id"])
        if fid_num > 10:
            fid_num -= 10
        elif fid_num > 6:
            fid_num -= 6
        row = np.argwhere(
            (calalign.fids["FID_NUM_SI"] == fid_num) & (calalign.fids["FID_SI"] == obs_si)
        )

        fid_rows.append(row[0][0])

    fid_cal = calalign.fids[fid_rows]

    # 4th-order LSI correction: 5 coefficients, exponents 0..4.
    # See PREDICTED_FID_AND_ASP_SOLVE.md §6 and Aca_Id_Image.cc:1939-1940.
    power = np.arange(5)
    fid_pos_dy = np.sum(fid_cal["FID_Y_CORR"].data * simz_offset ** power, axis=1)
    fid_pos_dz = np.sum(fid_cal["FID_Z_CORR"].data * simz_offset ** power, axis=1)

    # clip to limits
    fid_pos_dy = np.clip(fid_pos_dy, fid_cal["FID_Y_LIM"][:, 0], fid_cal["FID_Y_LIM"][:, 1])
    fid_pos_dz = np.clip(fid_pos_dz, fid_cal["FID_Z_LIM"][:, 0], fid_cal["FID_Z_LIM"][:, 1])

    # these are the fid position values that go into fidprops.p_lsi:
    fid_pos_lsi = fid_cal["FID_POS_LSI"].data.copy()
    fid_pos_lsi[:, 1] -= fid_pos_dy
    fid_pos_lsi[:, 2] -= fid_pos_dz

    if as_table:
        fid_pos = Table(
            {
                "slot": fids["slot"],
                "fid_id": fids["id"],
                "pos_lsi": fid_pos_lsi,
            }
        )
        return fid_pos

    return fid_pos_lsi


def estimated_fid_positions(
        *,
        obspar,
        calalign,
        fid_pos_lsi,
        theta_X = 0.0,  # initial guess for theta_X (in degrees)
        dy = 0.0,  # initial guess for dy (in mm)
        dz = 0.0,  # initial guess for dz (in mm)
        as_table=True,
    ):
    """Compute the estimated fiducial positions in the ACA frame.
     
    This computes the "estimated" positions for fid lights, applying the misalignment matrices.
    """
    # This needs values of theta_X, dy and dz, which are initially zero
    # this is an iterative process, where the values of theta_X, dy and dz are updated to minimize
    # the difference between the observed and predicted fid positions.

    # STT origin in STF (obspar['sim_z'] in mm)
    # the following should be the same as the SIM_X/Y/Z values in the sim$root_coor0a.fits file (need to check)
    stt0_stf = np.array([obspar['sim_x'], obspar['sim_y'], obspar['sim_z']]) #  (?)
    lsi0_stt = calalign.align["LSI0_STT"]  # LSI origin in STT frame (in mm)
    rrc0_fc_x = calalign.align["RRC0_FC_X"]  # FC-frame X coordinate of the RRC origin (HRMA focus)

    # note for self...
    # the previous teps are the same as (from fid_props header):
    # lsi0_stt_2 = [fid_props_header["LSI0STT%d" % x] for x in [1, 2, 3]]
    # stt0_stf_2 = [fid_props_header["STT0STF%d" % x] for x in [1, 2, 3]]
    # rrc0_fc_x_2 = fid_props_header["RRC0FCX"]

    # LSI -> STT -> STF translation
    fid_pos_stf = fid_pos_lsi + lsi0_stt + stt0_stf

    # Project on FC frame (roll and translate). theta_X is in degrees,
    # matching aspsol["adtheta"] and the C `ACA_Solve::calcFidPos` convention
    # (`sind(theta_X)` / `cosd(theta_X)`, ACA_Solve.cc:148-149).
    #
    # The C builds Rot_X with `Rot_X[1][2]=+s, Rot_X[2][1]=-s` (ACA_Solve.cc:152-155)
    # and applies `p_fc = Rot_X · p_stf`. That matrix is Rx(theta_X)^T (i.e.
    # Rx(-theta_X)) in the right-hand-rule convention. We use the same
    # transpose form here so the sign of `adtheta` matches the C.
    Rot_X = Rotation.from_euler('x', theta_X, degrees=True).as_matrix()
    if Rot_X.ndim == 3:
        fid_pos_fc = np.einsum('nji,nj->ni', Rot_X, fid_pos_stf)
    else:
        fid_pos_fc = fid_pos_stf @ Rot_X

    fid_pos_fc = fid_pos_fc + np.stack([np.zeros_like(dy), dy, dz], axis=-1)
    # rrc0_fc_x is the FC-frame X coordinate of the RRC origin (HRMA focus), also a scalar FIDPROPS keyword.
    # The leading *= -1.0 is a sign convention so positive d_fc[0] points outward from the focal plane.
    v = fid_pos_fc - np.array([rrc0_fc_x, 0.0, 0.0])
    d_fc = -v / np.linalg.norm(v, axis=-1, keepdims=True)

    # (Mis)alignment matrices:
    # - ACA_SC_ALIGN (ACA nominal alignment): rotation FC -> ACA-nominal
    # - ACA_MISALIGN (ACA misalignment): rotation ACA-nominal -> ACA-actual
    # - FTS_MISALIGN (FTS misalignment): small misalignment of the Fiducial Transfer System
    #   (the periscope/fiducial-light optical path) within the FC frame
    M = calalign.align["ACA_MISALIGN"] @ calalign.align["ACA_SC_ALIGN"] @ calalign.align["FTS_MISALIGN"]

    # d_fc is (N_fid, 3); right-multiply by M.T so each row is M @ d_fc[i].
    d_aca_est = d_fc @ M.T

    yang = np.degrees(np.arctan2(d_aca_est[:, 1], d_aca_est[:, 0])) * 3600
    zang = np.degrees(np.arctan2(d_aca_est[:, 2], d_aca_est[:, 0])) * 3600

    if as_table:
        est = Table(
            {
                "ang_y": yang,
                "ang_z": zang,
                "d_aca": d_aca_est,
            }
        )
        return est

    return d_aca_est


def _odd(n):
    n = int(np.ceil(n))
    return n + 1 if n % 2 == 0 else n


def sg_smooth_per_slot(*, time, slot, ang_y, ang_z, t_smooth=40.0, polyorder=4):
    """Savitzky-Golay smooth ang_y / ang_z per slot.

    Window length per slot is `ceil_to_odd(t_smooth / median(dt))`, capped
    so it never exceeds the slot's sample count (must also stay > polyorder).
    Edge handling is `mode="nearest"`, which matches the C tool's
    nearest-good-sample replacement at endpoints.
    """
    # note that this function splits according to slot, but not according to algorithm. There should
    # be no more than one algorithm per slot.
    time = np.asarray(time, dtype=float)
    slot = np.asarray(slot)
    y_sm = np.asarray(ang_y, dtype=float).copy()
    z_sm = np.asarray(ang_z, dtype=float).copy()

    for s in np.unique(slot):
        sel = slot == s
        t_s = time[sel]
        if t_s.size <= polyorder + 1:
            continue
        dt = np.median(np.diff(t_s))
        window = min(_odd(t_smooth / dt), _odd(t_s.size - 1))
        if window <= polyorder:
            continue
        y_sm[sel] = savgol_filter(y_sm[sel], window, polyorder, mode="nearest")
        z_sm[sel] = savgol_filter(z_sm[sel], window, polyorder, mode="nearest")

    return y_sm, z_sm


def hann_smooth(values, *, windowlen=152):
    """Hann-windowed smooth with replicated-endpoint padding."""
    values = np.asarray(values, dtype=float)
    kernel = hann(windowlen, sym=True)
    kernel /= kernel.sum()
    pad = windowlen // 2
    padded = np.pad(values, pad, mode="edge")
    return np.convolve(padded, kernel, mode="same")[pad:pad + values.size]


if __name__ == "__main__":
    # obsid = 28182
    obsid = 29622
    # obsid = 29878  # large sim-Z
    obspar = get_obspar(obsid)
    acas = get_starcats(obsid=obsid)[0]
    fids = acas[acas["type"] == "FID"]
    obs_si = obspar["detector"]

    centroid_files = get_files(obsid, content=["ACACENT"])
    fidprops_files = get_files(obsid, content=["FIDPROPS"])
    aspsol_files = get_files(obsid, content=["ASPSOL"])

    centroid_file = centroid_files[0]
    with fits.open(centroid_file) as cent_hdus:
        centroids = Table(cent_hdus[1].data)

        # discard star centroids and keep only fids with Gaussian fit (EGAUSS) algorithm.
        centroids = centroids[np.isin(centroids["slot"], fids["slot"]) & (centroids["alg"] == EGAUSS)]

    fidprops_file = fidprops_files[0]
    with fits.open(fidprops_file) as fidpr_hdus:
        # ACIS correction parameters are estimated in the pipeline and stored in the FIDPROPS file.
        # We will not re-derive them.
        dh_temp_calib = {
            "dht_mean": float(fidpr_hdus[1].header["DHTMEAN"]),
            "deg_per_cnt": float(fidpr_hdus[1].header["DEGPCNT"])
        }
        fid_props = Table(fidpr_hdus[1].data)

    aspsol_file = aspsol_files[0]
    with fits.open(aspsol_file) as aspsol_hdus:
        aspsol = Table(aspsol_hdus[1].data)
        aspsol_header = aspsol_hdus[1].header

    calalign = get_caldb_calalign(obspar)

    fid_pos_nominal, fid_pos_lsi = nominal_fid_positions(
        obspar=obspar,
        calalign=calalign,
        fids=fids
    )

    # for comparison, we will call estimated_fid_positions with two sets of arguments:
    # 1. with the initial values of theta_X, dy and dz (all zero),
    # 2. with the values of theta_X, dy and dz derived from the ASPSOL file (the final solution).

    idx = np.searchsorted(fid_props["slot"], centroids["slot"])
    cen_pos_lsi = fid_pos_lsi[idx]

    get_dy = interp1d(aspsol["time"], aspsol["ady"], bounds_error=False, fill_value="extrapolate")
    get_dz = interp1d(aspsol["time"], aspsol["adz"], bounds_error=False, fill_value="extrapolate")
    get_theta_X = interp1d(aspsol["time"], aspsol["adtheta"], bounds_error=False, fill_value="extrapolate")

    dy = get_dy(centroids["time"])
    dz = get_dz(centroids["time"])
    # aspsol["adtheta"] is in degrees; estimated_fid_positions takes degrees.
    dtheta = get_theta_X(centroids["time"])

    fid_pos_expected, _ = estimated_fid_positions(
        obspar=obspar,
        calalign=calalign,
        fid_pos_lsi=fid_pos_lsi,
        dy=dy,
        dz=dz,
        theta_X=dtheta,
    )
    fid_pos_expected_0, _ = estimated_fid_positions(
        obspar=obspar,
        calalign=calalign,
        fid_pos_lsi=fid_pos_lsi,
        dy=0,
        dz=0,
        theta_X=0,
    )

    fid_pos_measured, _ = process_centroids(
        centroids=centroids,
        dh_temp_calib=dh_temp_calib,
        fid_props=fid_props,
        obspar=obspar,
        calalign=calalign
    )