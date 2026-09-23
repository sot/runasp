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
from scipy.spatial.transform import Rotation

CALDB_DIR = Path(os.environ.get("CALDB", "/home/ascds/DS.release/CALDB"))


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


def process_centroids(*, centroids, acis_corr_params, obspar, calalign):
    """Process the fiducial centroids to apply thermal and periscope corrections,
    and compute direction cosines in the ACA frame.

    Simplifications vs the C `aca_corr_fid` pipeline (kept on purpose —
    each is a deliberate choice in this Python re-implementation):

    - Operates on raw `ang_y` / `ang_z`, NOT the Savitzky-Golay smoothed
      `ang_y_sm` / `ang_z_sm` that the C tool consumes. Skipping the SG
      smooth trades a small high-frequency residual for not having to
      re-implement the upstream `aca_corr_centr` window.
    - Periscope correction uses cheta `OOBAGRD3` / `OOBAGRD6` and
      subtracts the obs-window mean. The C tool reads the OBCENG
      `smoothGradient3/6` columns and `meanGradient3/6` keywords —
      pre-smoothed and pre-meaned upstream. Both gradients drift slowly,
      so the difference is small but non-zero.
    - DH-temperature is built from the arithmetic mean of `counts`
      across fid slots per timestamp. The C pipeline does a per-fid
      windowed median first, then averages across fids; with cosmic-ray
      hits the median vs mean differs by tens of counts.
    - Center of expansion uses CALALIGN `FID_Y/Z_CENTER_NOM` directly.
      The C tool uses `obs_y/z_coe = nominal + obs-specific offset` (the
      offset comes from this observation's own `ang_y_nom` column and
      `<ang_y_sm>`).
    """
    ang_y = np.asarray(centroids["ang_y"], dtype=float).copy()
    ang_z = np.asarray(centroids["ang_z"], dtype=float).copy()
    t = np.asarray(centroids["time"], dtype=float)
    slot = np.asarray(centroids["slot"], dtype=int)

    # ACIS base-plate thermal correction.
    if np.isin(obspar["detector"], calalign.acis_corr["DETECTOR"]):
        det_mask = (obspar["detector"] == calalign.acis_corr["DETECTOR"])
        acis_corr = calalign.acis_corr[det_mask][0]
        dht_mean = acis_corr_params["dht_mean"]
        deg_per_cnt = acis_corr_params["deg_per_cnt"]

        # Per-time DH temp from fid counts; see the third bullet in the
        # docstring and aca_corr_fid/tool_functions.c:548-553.
        counts = np.asarray(centroids["counts"], dtype=float)
        _, inv = np.unique(t, return_inverse=True)
        counts_per_t = np.bincount(inv, weights=counts) / np.bincount(inv)
        dh_temp_C = (dht_mean + (counts_per_t - counts_per_t.mean()) * deg_per_cnt)[inv]

        dh_base_C = float(acis_corr["DH_TEMP_BASE"]) - 273.15  # in °C
        fid_cte = float(acis_corr["FID_CTE"])
        fid_y_center = float(acis_corr["FID_Y_CENTER_NOM"])
        fid_z_center = float(acis_corr["FID_Z_CENTER_NOM"])
        corr = (dh_temp_C - dh_base_C) * fid_cte
        ang_y -= corr * (ang_y - fid_y_center)
        ang_z -= corr * (ang_z - fid_z_center)

    # Periscope correction. ±1000s pad on the cheta fetch so the
    # interpolation has neighbours at the edges of the obs window.
    msids = ["OOBAGRD3", "OOBAGRD6"]
    times = np.unique(centroids["time"])
    telem = fetch.Msidset(msids, obspar["tstart"] - 1000, obspar["tstop"] + 1000)
    telem.interpolate(times=times)

    oobagrd3 = telem["OOBAGRD3"].vals
    oobagrd6 = telem["OOBAGRD6"].vals
    idx = np.searchsorted(telem["OOBAGRD3"].times, centroids["time"])

    g3_dev = oobagrd3 - oobagrd3.mean()
    g6_dev = oobagrd6 - oobagrd6.mean()
    dy = (g3_dev * float(calalign.peri_corr["PERI_Y_OOBAGRD3"][0])
          + g6_dev * float(calalign.peri_corr["PERI_Y_OOBAGRD6"][0]))
    dz = (g3_dev * float(calalign.peri_corr["PERI_Z_OOBAGRD3"][0])
          + g6_dev * float(calalign.peri_corr["PERI_Z_OOBAGRD6"][0]))
    ang_y -= dy[idx]
    ang_z -= dz[idx]
    # Tangent-projection form (1, tan(ang_y), tan(ang_z)) — matches the
    # C convertMeasurements build at ACA_Solve.cc:1928. NOT a true
    # direction cosine; the [:, 0] component is fixed at 1, not the
    # vector length. The C pipeline calls .unitize() row-by-row before
    # the LM fit (ACA_Solve.cc:2269) so mrqmin sees a unit vector. This
    # script leaves it un-unitized; the per-row mismatch is
    # ~tan**3(ang)/3 ~ 1e-10 at typical fid offsets, far below pipeline
    # precision. See B8 in process_fids_review.md.
    d_aca_meas = np.stack(
        [
            np.ones_like(ang_y),
            np.tan(np.radians(ang_y)),
            np.tan(np.radians(ang_z)),
        ],
        axis=-1,
    )
    # this is what gets written in centroids file?
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
    return meas, d_aca_meas


def nominal_fid_positions(*, obspar, calalign, fids):
    """Compute the nominal fiducial positions in the LSI frame, applying the polynomial correction.
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

    return fid_pos_lsi


def estimated_fid_positions(
        *,
        obspar,
        calalign,
        fid_pos_lsi,
        theta_X = 0.0,  # initial guess for theta_X (in degrees)
        dy = 0.0,  # initial guess for dy (in mm)
        dz = 0.0,  # initial guess for dz (in mm)
    ):
    """Compute the estimated fiducial positions in the ACA frame.
     
    This computes the "estimated" positions for fid lights, applying the misalignment matrices.
    """
    # STT origin in STF (obspar['sim_z'] in mm)
    # the following should be the same as the SIM_X/Y/Z values in the sim$root_coor0a.fits file (need to check)
    stt0_stf = np.array([obspar['sim_x'], obspar['sim_y'], obspar['sim_z']]) #  (?)
    lsi0_stt = calalign.align["LSI0_STT"]  # LSI origin in STT frame (in mm)
    rrc0_fc_x = calalign.align["RRC0_FC_X"]  # FC-frame X coordinate of the RRC origin (HRMA focus)

    # LSI -> STT -> STF translation
    fid_pos_stf = fid_pos_lsi + lsi0_stt + stt0_stf
    # Project on FC frame (roll and translate)
    # This needs values of theta_X, dy and dz, which are initially zero
    # this is an iterative process, where the values of theta_X, dy and dz are updated to minimize
    # the difference between the observed and predicted fid positions.

    Rot_X = Rotation.from_euler('x', theta_X, degrees=True).as_matrix()
    fid_pos_fc = fid_pos_stf @ Rot_X.T + np.array([0.0, dy, dz])

    # rrc0_fc_x is the FC-frame X coordinate of the RRC origin (HRMA focus), a scalar FIDPROPS keyword.
    # The leading *= -1.0 is a sign convention so positive d_fc[0] points outward from the focal plane.
    v = fid_pos_fc - np.array([rrc0_fc_x, 0.0, 0.0])
    d_fc = -v / np.linalg.norm(v, axis=1, keepdims=True)

    # (Mis)alignment matrices:
    # - ACA_SC_ALIGN (ACA nominal alignment): rotation FC -> ACA-nominal
    # - ACA_MISALIGN (ACA misalignment): rotation ACA-nominal -> ACA-actual
    # - FTS_MISALIGN (FTS misalignment): small misalignment of the Fiducial Transfer System
    #   (the periscope/fiducial-light optical path) within the FC frame
    M = calalign.align["ACA_MISALIGN"] @ calalign.align["ACA_SC_ALIGN"] @ calalign.align["FTS_MISALIGN"]

    # d_fc is (N_fid, 3); right-multiply by M.T so each row is M @ d_fc[i].
    d_aca_est = d_fc @ M.T

    return d_aca_est


if __name__ == "__main__":
    # an example OBSID
    obsid = 29878
    obspar = get_obspar(obsid)
    acas = get_starcats(obsid=obsid)[0]
    fids = acas[acas["type"] == "FID"]
    obs_si = obspar["detector"]

    # centroid_files = list(Path("test_data").glob(f"pcadf*_acen1.fits.gz"))
    centroid_files = get_files(obsid, content=["ACACENT"])
    fidprops_files = get_files(obsid, content=["FIDPROPS"])

    centroid_file = centroid_files[0]
    centroids = Table(fits.open(centroid_file)[1].data)

    fidprops_file = fidprops_files[0]
    with fits.open(fidprops_file) as fidpr_hdus:
        acis_corr_params = {
            "dht_mean": float(fidpr_hdus[1].header["DHTMEAN"]),
            "deg_per_cnt": float(fidpr_hdus[1].header["DEGPCNT"])
        }

    calalign = get_caldb_calalign(obspar)
    fid_pos_lsi = nominal_fid_positions(obspar=obspar, calalign=calalign, fids=fids)
    d_aca_est = estimated_fid_positions(obspar=obspar, calalign=calalign, fid_pos_lsi=fid_pos_lsi)
    meas, d_aca_meas = process_centroids(centroids=centroids, acis_corr_params=acis_corr_params, obspar=obspar, calalign=calalign)