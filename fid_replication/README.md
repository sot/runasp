# Fid replication

A Python reproduction of the fid-light half of the DS aspect pipeline: the corrections
`aca_corr_fid` applies to the measured centroids, and the predicted fid positions that
`asp_solve` compares them against. It reads archived products through mica, cheta and
kadi, so it runs over any observation without rerunning the pipeline.

It lives beside `runasp` because it answers the same kind of question — what the aspect
pipeline did to an observation — by reading rather than by running.

## Contents

| file | what it is |
|---|---|
| `fid_processing.py` | both chains, and `compare_obsid()` which runs them over one observation |
| `fid_smoothing.py` | the two smoothers the pipeline applies |
| `fid_deltas_mica.py` | the predicted chain again, from ACACAL and the FIDPROPS header, as a cross-check |
| `docs/fid-replication.md` | how each piece maps onto the C, and where it deliberately differs |
| `docs/fid-flowchart.md` | the whole fid path as one diagram |

## Running it

Needs a Ska environment with mica, cheta and kadi configured, and `$CALDB` pointing at a
CALDB tree with the PCAD alignment files (`data/chandra/pcad/align/`).

```
python fid_processing.py --obsid 29878
```

prints, per fid slot, the mean and spread of measured minus predicted angle in arcsec,
with the prediction evaluated at the aspect solution the pipeline found.

## Status

The chains are complete and the comparison runs, but the residuals have not been checked
end to end against an independent computation; obsid 29878, with its large SIM-Z offset,
is the case that exercises the correction polynomial hardest and is the one to verify
first. The deliberate departures from the C are listed in `docs/fid-replication.md`;
they are the first thing to account for in any residual.

The pipeline itself is described on the ska-wiki: *Aspect L1 pipeline*, *ASP L1
products*, *Fid light corrections* and *Predicted fid positions*.
