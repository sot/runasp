import numpy as np
from astropy.io import fits
from astropy.table import Table
from kadi.commands import get_starcats
from mica.archive.asp_l1 import get_files
from mica.archive.obspar import get_obspar


def get_fid_deltas_mica(
    *,
    asol,
    h_fidpr,
    aca_misalign,
    fts_misalign,
    fidprop,
    cen,
):

    from Ska.astro import sph_dist

    R2A = 206264.81
    D2R = 0.017453293

    lsi0_stt = [h_fidpr["LSI0STT%d" % x] for x in [1, 2, 3]]
    stt0_stf = [h_fidpr["STT0STF%d" % x] for x in [1, 2, 3]]
    rrc0_fc_x = h_fidpr["RRC0FCX"]

    M = np.dot(aca_misalign, fts_misalign)

    rot_x = np.zeros([3, 3])
    rot_x[0, 0] = 1

    if "ady" in asol.colnames:
        dy_col = "ady"
        dz_col = "adz"
        dtheta_col = "adtheta"
    else:
        dy_col = "dy"
        dz_col = "dz"
        dtheta_col = "dtheta"

    dy = np.zeros(len(cen["time"]))
    dz = np.zeros(len(cen["time"]))
    dr = np.zeros(len(cen["time"]))
    for fid in fidprop:
        p_lsi = fid["p_lsi"]
        p_stf = p_lsi + lsi0_stt + stt0_stf

        ok = cen["slot"] == fid["slot"]
        ceni = cen[ok]
        asol_cen_dy = np.interp(ceni["time"], asol["time"], asol[dy_col])
        asol_cen_dz = np.interp(ceni["time"], asol["time"], asol[dz_col])
        asol_cen_dtheta = (
            np.interp(ceni["time"], asol["time"], asol[dtheta_col]) * D2R
        )

        rot_x = np.zeros([len(ceni["time"]), 3, 3])
        s_th = np.sin(asol_cen_dtheta)
        c_th = np.cos(asol_cen_dtheta)
        rot_x[:, 0, 0] = 1.0
        rot_x[:, 1, 1] = c_th
        rot_x[:, 2, 1] = s_th
        rot_x[:, 1, 2] = -s_th
        rot_x[:, 2, 2] = c_th
        p_fc = np.dot(rot_x.transpose(0, 2, 1), p_stf)
        p_fc[:, 1] = p_fc[:, 1] + asol_cen_dy
        p_fc[:, 2] = p_fc[:, 2] + asol_cen_dz
        d_fc = p_fc
        d_fc[:, 0] = d_fc[:, 0] - rrc0_fc_x
        d_fc = -d_fc
        d_aca = np.dot(d_fc, M.transpose())
        yag = np.arctan2(d_aca[:, 1], d_aca[:, 0]) * R2A
        zag = np.arctan2(d_aca[:, 2], d_aca[:, 0]) * R2A
        dy[ok] = ceni["ang_y_sm"] * 3600 - yag
        dz[ok] = ceni["ang_z_sm"] * 3600 - zag
        dr[ok] = (
            sph_dist(yag / 3600, zag / 3600, ceni["ang_y_sm"], ceni["ang_z_sm"])
            * 3600
        )
    return dy, dz, dr


if __name__ == "__main__":
    # obsid = 28182
    # obsid = 29622
    obsid = 29878  # large sim-Z
    obspar = get_obspar(obsid)
    acas = get_starcats(obsid=obsid)[0]
    fids = acas[acas["type"] == "FID"]
    obs_si = obspar["detector"]

    centroid_files = get_files(obsid, content=["ACACENT"])
    fidprops_files = get_files(obsid, content=["FIDPROPS"])
    aspsol_files = get_files(obsid, content=["ASPSOL"])
    acacal_files = get_files(obsid, content=["ACACAL"])

    centroid_file = centroid_files[0]
    with fits.open(centroid_file) as cent_hdus:
        centroids = Table(cent_hdus[1].data)

        # discard star centroids and keep only fids with Gaussian fit (EGAUSS) algorithm.
        centroids = centroids[np.isin(centroids["slot"], fids["slot"]) & (centroids["alg"] == 8)]

    fidprops_file = fidprops_files[0]
    with fits.open(fidprops_file) as fidpr_hdus:
        fid_props = Table(fidpr_hdus[1].data)
        fid_props_header = fidpr_hdus[1].header


    aspsol_file = aspsol_files[0]
    with fits.open(aspsol_file) as aspsol_hdus:
        aspsol = Table(aspsol_hdus[1].data)
        aspsol_header = aspsol_hdus[1].header

    acacal_file = acacal_files[0]
    with fits.open(acacal_file) as acacal_hdus:
        acacal = Table(acacal_hdus[1].data)

    centroids["dy_mica"], centroids["dz_mica"], centroids["dr_mica"] = get_fid_deltas_mica(
        asol = aspsol,
        h_fidpr = fid_props_header,
        aca_misalign = acacal["aca_misalign"][0],
        fts_misalign = acacal["fts_misalign"][0],
        fidprop = fid_props,
        cen = centroids,
    )