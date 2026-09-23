# Drafts — superseded, kept for provenance

These are the working notes the published material came from, archived unchanged on
2026-09-23 and **not maintained**. The knowledge in them that survived review is on the
ska-wiki — *Aspect L1 pipeline*, *ASP L1 products*, *Fid light corrections*, *Predicted
fid positions*, *Attitude chain*, *runasp*, *Running the aspect pipeline*, *HEAD
machines*, *Working with the DS source* — and in `../docs/`. Read those; read these only
to see where something came from.

They are archived rather than deleted because they carry the reading trail — the C file
and line references behind each claim — which is worth having when a page is challenged.

## Do not cite these without checking

Claims in here that were found wrong when the pages were written and reviewed against the
DS source:

| claim in the drafts | what is true |
|---|---|
| ACACENT `ang_y`/`ang_z` are arcmin (`docs/CENTROID_PRODUCTS.md`) | degrees, per the product template |
| `aca_corr_fid` runs before `aca_filter_centr` (`docs/ASPECT_PIPELINE_OVERVIEW.md`) | the other way round, per `asp_l1_std.ped` |
| Levenberg-Marquardt is the default fit (`docs/PREDICTED_FID_AND_ASP_SOLVE.md`, `docs/FID_FLOWCHART.md`) | the pipeline passes `max_iterations=0`, which selects the analytic solve |
| the periscope correction was added in 2008 (several) | 2008 is the ACIS thermal correction; the periscope one is 2011 |
| all four tool names in the ACACENT `HISTORY` mean a fully corrected file (`docs/CENTROID_PRODUCTS.md`) | `aca_corr_fid` writes no history record at all, so this test always fails |
| bad samples are substituted before each smoothing iteration (`docs/FID_CORRECTION_STAGES.md`) | once, before smoothing, for bad centroids; sigma-rejected samples take a median-filtered value |
| stage names and numbers from `runasp.py` (`docs/FID_LIGHT_PROCESSING.md`, `docs/FID_FLOWCHART.md`) | that list is stale — 12 stages missing, 6 that no longer exist |
| `aca_id_image` falls back to a text file when the OCAT database is unreachable (`docs/ASPECT_PIPELINE_OVERVIEW.md`) | it has only ever read a text file; `asp_read_ocat` is the database client |
| `asp_kalman` is one tool (`docs/ASPECT_PIPELINE_OVERVIEW.md`) | a directory building `asp_forward_kalman` and `asp_smooth_kalman` |
| ACISFIDCORR and PERIFIDCORR are CALDB products (`docs/ASPECT_PIPELINE_OVERVIEW.md`) | extensions of the CALALIGN file |
| AGASC is the "Authorized Guide Star Catalog" (`docs/ASPECT_PIPELINE_OVERVIEW.md`) | the AXAF Guide and Acquisition Star Catalog |

`process_fids_review.md` and `docs/CENTROID_PRODUCTS.md` Part 2 are superseded designs:
the first is a fix-list against a plan file that no longer exists, the second applies the
corrections without smoothing, which the published implementation does not do.

`process_fids.py`, `process_fids_v1.py`, `process_fids_v2.py`, `fid_process_meas.py` and
`fid_process_pred.py` are the lineage of `../fid_processing.py`, which supersedes all of
them.
