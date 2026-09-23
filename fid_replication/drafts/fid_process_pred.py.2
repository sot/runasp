from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table
from kadi.commands import get_starcats
from mica.archive.obspar import get_obspar
from scipy.spatial.transform import Rotation

obsid = 29878
obspar = get_obspar(obsid)
acas = get_starcats(obsid=obsid)[0]
fids = acas[acas["type"] == "FID"]
obs_si = obspar["detector"]
files = list(Path("test_data").glob("*"))
# files
# Info from CalDB is stored in the acal file below (we'll check)
# cal_hdus = fits.open("test_data/pcadf893066380N001_acal1.fits.gz")

calalign_hdus = fits.open(
    "/home/ascds/DS.release/CALDB/data/chandra/pcad/align/pcadD2021-07-02alignN0010.fits"
)

calalign_1 = Table(calalign_hdus[1].data)
# print("## Misalignment matrices")
# calalign_1.pprint_all()
idx = np.argwhere(calalign_1["INSTR_ID"] == obs_si)[0][0]

# STT origin in STF (obspar['sim_z'] in mm)
# the following should be the same as the SIM_X/Y/Z values in the sim$root_coor0a.fits file (need to check)
stt0_stf = np.array([obspar["sim_x"], obspar["sim_y"], obspar["sim_z"]])  #  (?)
# LSI origin in STT frame (in mm)
lsi0_stt = calalign_1["LSI0_STT"][idx]
# FC-frame X coordinate of the RRC origin (HRMA focus)
rrc0_fc_x = calalign_1["RRC0_FC_X"][idx]
# (Mis)alignment matrices
# ACA nominal alignment
# rotation FC -> ACA-nominal
aca_align = calalign_1["ACA_SC_ALIGN"][idx]
# ACA misalignment
# rotation ACA-nominal -> ACA-actual
aca_misalign = calalign_1["ACA_MISALIGN"][idx]
# FTS misalignment
# small misalignment of the Fiducial Transfer System (the periscope/fiducial-light optical path) within the FC frame
fts_misalign = calalign_1["FTS_MISALIGN"][idx]

calalign_2 = Table(calalign_hdus[2].data)
# get rid of trailing spaces in FID_SI column
calalign_2["FID_SI"] = [si.strip() for si in calalign_2["FID_SI"]]
# print("## fid positions and 5-degree polynomial parameters for corrections:\n")
# calalign_2.pprint_all()

calalign_3 = Table(calalign_hdus[3].data)
# print("## ACIS base-plate fid correction parameters:\n")
# calalign_3.pprint_all()

calalign_4 = Table(calalign_hdus[4].data)
# print("## Periscope correction parameters:\n")
# calalign_4.pprint_all()

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
        (calalign_2["FID_NUM_SI"] == fid_num) & (calalign_2["FID_SI"] == obs_si)
    )[0][0]
    fid_rows.append(row)
fid_cal = calalign_2[fid_rows]

simz_offset = lsi0_stt[2] + stt0_stf[2]  # mm

# pol-5 LSI correction
# is it 5-degree polynomial with zero intercept?
# power = np.arange(6)[1:]
# or a 4-degree polynomial?
power = np.arange(5)
fid_pos_dy = np.sum(fid_cal["FID_Y_CORR"].data * simz_offset**power, axis=1)
fid_pos_dz = np.sum(fid_cal["FID_Z_CORR"].data * simz_offset**power, axis=1)

# clip to limits
fid_pos_dy = np.clip(fid_pos_dy, fid_cal["FID_Y_LIM"][:, 0], fid_cal["FID_Y_LIM"][:, 1])
fid_pos_dz = np.clip(fid_pos_dz, fid_cal["FID_Z_LIM"][:, 0], fid_cal["FID_Z_LIM"][:, 1])

# these are the fid position values that go into fidprops.p_lsi:
fid_pos_lsi = fid_cal["FID_POS_LSI"].data.copy()
fid_pos_lsi[:, 1] -= fid_pos_dy
fid_pos_lsi[:, 2] -= fid_pos_dz

# LSI -> STT -> STF translation
fid_pos_stf = fid_pos_lsi + lsi0_stt + stt0_stf
# Project on FC frame (roll and translate)
# This needs values of theta_X, dy and dz, which are initially zero
# this is an iterative process, where the values of theta_X, dy and dz are updated to minimize
# the difference between the observed and predicted fid positions.
theta_X = 0.0  # initial guess for theta_X (in degrees)
dy = 0.0  # initial guess for dy (in mm)
dz = 0.0  # initial guess for dz (in mm)


Rot_X = Rotation.from_euler("x", theta_X, degrees=True).as_matrix()
fid_pos_fc = fid_pos_stf @ Rot_X.T + np.array([0.0, dy, dz])

# rrc0_fc_x is the FC-frame X coordinate of the RRC origin (HRMA focus), a scalar FIDPROPS keyword.
# The leading *= -1.0 is a sign convention so positive d_fc[0] points outward from the focal plane.

v = fid_pos_fc - np.array([rrc0_fc_x, 0.0, 0.0])
d_fc = -v / np.linalg.norm(v, axis=1, keepdims=True)
# misalignment matrices
M = aca_misalign @ aca_align @ fts_misalign
d_aca_est = d_fc @ M.T
