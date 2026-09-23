"""Python reproduction of the fid-light half of the DS aspect pipeline.

Two chains, the same two the pipeline runs and `asp_solve` compares:

- the measured chain of `aca_corr_fid` — Savitzky-Golay smoothing, the ACIS
  thermal correction and the periscope correction applied to ACACENT angles,
  giving direction cosines in the ACA frame;
- the predicted chain of `aca_id_image` — the SIM-Z correction polynomial on
  each fid's nominal LSI position, projected LSI -> STT -> STF -> FC -> ACA
  through the CALALIGN alignment matrices.

It is built on mica, cheta and kadi, which the C pipeline does not have, so it
can be run over archived observations without rerunning the pipeline. It is a
reading aid and a cross-check, not a reimplementation: the deliberate
departures from the C are listed in the docstrings and in
docs/fid-replication.md.
"""

import argparse
import os
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table
from cheta import fetch
from kadi.commands import get_starcats
from mica.archive.asp_l1 import get_files
from mica.archive.obspar import get_obspar
from scipy.interpolate import interp1d
from scipy.spatial.transform import Rotation

from fid_smoothing import hann_smooth, sg_smooth_per_slot


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

    - The periscope baseline is the mean of the Hann-smoothed gradients over
      the observation. The C freezes the mean of the *raw* gradients before
      smoothing them.

    Centroid SG smoothing (matches `aca_filter_centr`) and OOBAGRD3/6
    Hann smoothing (matches `aca_corr_fid`) come from `fid_smoothing`. The thermal-correction center of expansion follows
    `aca_corr_fid/tool_functions.c`: `obs_y/z_coe = FID_Y/Z_CENTER_NOM + δ`,
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
        # docstring and aca_corr_fid/tool_functions.c.
        counts = np.asarray(centroids["counts"], dtype=float)
        counts_per_t = np.bincount(idx, weights=counts) / np.bincount(idx)
        dh_temp_C = (dht_mean + (counts_per_t - counts_per_t.mean()) * deg_per_cnt)[idx]

        dh_base_C = float(acis_corr["DH_TEMP_BASE"]) - 273.15  # in °C
        fid_cte = float(acis_corr["FID_CTE"])

        # Per-obs center of expansion (aca_corr_fid/tool_functions.c).
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
    # Taken from the obspar. The pipeline reads the same quantity from the SIM
    # coordinate product (sim<root>_coor0a.fits), quality-gated per axis; the two
    # agree for a stable SIM but are not guaranteed to. See docs/fid-replication.md.
    stt0_stf = np.array([obspar['sim_x'], obspar['sim_y'], obspar['sim_z']])

    # LSI origin in STT frame (in mm)
    lsi0_stt = calalign.align["LSI0_STT"]

    simz_offset = lsi0_stt[2] + stt0_stf[2]  # mm

    # select fid light rows.
    # OCAT fid id is absolute (1-14); CALALIGN.fid_num_si is per-instrument
    # relative (1-6 ACIS, 1-4 HRC). See Aca_Id_Image.cc in the DS source.
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
    # See docs/fid-replication.md and Aca_Id_Image.cc in the DS source.
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
    # Taken from the obspar. The pipeline reads the same quantity from the SIM
    # coordinate product (sim<root>_coor0a.fits), quality-gated per axis; the two
    # agree for a stable SIM but are not guaranteed to. See docs/fid-replication.md.
    stt0_stf = np.array([obspar['sim_x'], obspar['sim_y'], obspar['sim_z']])
    lsi0_stt = calalign.align["LSI0_STT"]  # LSI origin in STT frame (in mm)
    rrc0_fc_x = calalign.align["RRC0_FC_X"]  # FC-frame X coordinate of the RRC origin (HRMA focus)

    # The same three quantities are also written to the FIDPROPS header, as the
    # keywords LSI0STT1..3, STT0STF1..3 and RRC0FCX; fid_deltas_mica.py reads them
    # from there instead of from the calibration file.

    # LSI -> STT -> STF translation
    fid_pos_stf = fid_pos_lsi + lsi0_stt + stt0_stf

    # Project on FC frame (roll and translate). theta_X is in degrees,
    # matching aspsol["adtheta"] and the C `ACA_Solve::calcFidPos` convention
    # (`sind(theta_X)` / `cosd(theta_X)` in ACA_Solve.cc).
    #
    # The C builds Rot_X with `Rot_X[1][2]=+s, Rot_X[2][1]=-s` (ACA_Solve.cc)
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


def _aspsol_columns(aspsol):
    """Return the (dy, dz, dtheta) column names of an aspect solution table.

    A per-interval (CAI) solution carries `ady`, `adz`, `adtheta`, written by the
    pipeline's off-axis correction step; the per-OBI product that `asp_combine`
    writes drops those and keeps `dy`, `dz`, `dtheta`.
    """
    if "ady" in aspsol.colnames:
        return "ady", "adz", "adtheta"
    return "dy", "dz", "dtheta"


def compare_obsid(obsid):
    """Run both chains over one observation and return the per-sample comparison.

    The predicted positions are evaluated at the aspect solution the pipeline
    found, so the residuals below are what is left after the fit: they say how
    well this reproduction agrees with the pipeline, not how well the pipeline
    fitted the observation.
    """
    obspar = get_obspar(obsid)
    starcat = get_starcats(obsid=obsid)[0]
    fids = starcat[starcat["type"] == "FID"]

    with fits.open(get_files(obsid, content=["ACACENT"])[0]) as hdus:
        centroids = Table(hdus[1].data)
    # fid slots only, and only the algorithm the aspect solution uses
    centroids = centroids[
        np.isin(centroids["slot"], fids["slot"]) & (centroids["alg"] == EGAUSS)
    ]

    with fits.open(get_files(obsid, content=["FIDPROPS"])[0]) as hdus:
        # The pipeline derives these two in asp_get_calib and records them in the
        # FIDPROPS header; they are read rather than re-derived here.
        dh_temp_calib = {
            "dht_mean": float(hdus[1].header["DHTMEAN"]),
            "deg_per_cnt": float(hdus[1].header["DEGPCNT"]),
        }
        fid_props = Table(hdus[1].data)

    with fits.open(get_files(obsid, content=["ASPSOL"])[0]) as hdus:
        aspsol = Table(hdus[1].data)

    calalign = get_caldb_calalign(obspar)

    nominal = nominal_fid_positions(obspar=obspar, calalign=calalign, fids=fids)

    # one nominal LSI position per centroid sample, matched by slot
    row_of_slot = {int(slot): row for row, slot in enumerate(nominal["slot"])}
    rows = np.array([row_of_slot[int(slot)] for slot in centroids["slot"]])
    pos_lsi = np.asarray(nominal["pos_lsi"])[rows]

    dy_col, dz_col, dtheta_col = _aspsol_columns(aspsol)

    def at(column):
        return interp1d(
            aspsol["time"], aspsol[column], bounds_error=False, fill_value="extrapolate"
        )(centroids["time"])

    predicted = estimated_fid_positions(
        obspar=obspar,
        calalign=calalign,
        fid_pos_lsi=pos_lsi,
        dy=at(dy_col),
        dz=at(dz_col),
        theta_X=at(dtheta_col),  # degrees, as the solution stores it
    )
    measured = process_centroids(
        centroids=centroids,
        dh_temp_calib=dh_temp_calib,
        fid_props=fid_props,
        obspar=obspar,
        calalign=calalign,
    )

    comparison = Table(
        {
            "time": measured["time"],
            "slot": measured["slot"],
            "ang_y_meas": measured["ang_y"],
            "ang_z_meas": measured["ang_z"],
            "ang_y_pred": predicted["ang_y"],
            "ang_z_pred": predicted["ang_z"],
            "d_y": measured["ang_y"] - predicted["ang_y"],
            "d_z": measured["ang_z"] - predicted["ang_z"],
        }
    )
    for column in comparison.colnames[2:]:
        comparison[column].unit = "arcsec"
        comparison[column].format = ".3f"
    return comparison


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    # 29878 has a large SIM-Z offset, which is what exercises the correction
    # polynomial of the predicted chain.
    parser.add_argument("--obsid", type=int, default=29878)
    args = parser.parse_args()

    comparison = compare_obsid(args.obsid)
    print(f"obsid {args.obsid}: {len(comparison)} fid centroid samples")
    print(f"{'slot':>5} {'n':>6} {'dy mean':>10} {'dy std':>9} {'dz mean':>10} {'dz std':>9}")
    for slot in np.unique(comparison["slot"]):
        rows = comparison[comparison["slot"] == slot]
        print(
            f"{slot:>5} {len(rows):>6} {rows['d_y'].mean():>10.3f} {rows['d_y'].std():>9.3f} "
            f"{rows['d_z'].mean():>10.3f} {rows['d_z'].std():>9.3f}"
        )


if __name__ == "__main__":
    main()
