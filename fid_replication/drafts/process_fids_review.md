# `process_fids.py` — review report

Companion to the planning file at
`~/.claude/plans/hi-claude-i-created-distributed-valiant.md`. The plan
contained the full inventory of findings (sections A–D); this report
records which of them were applied to the code and which remain open.

## What `process_fids.py` is

A single-file Python reproduction of two C pipeline stages:

- **Predicted side** (`aca_id_image`): CALALIGN polynomial → corrected
  `p_lsi` → STF/FC projection → `d_aca_est` (predicted ACA-frame
  direction cosines).
- **Measured side** (`aca_corr_fid`): ACIS thermal + periscope
  corrections applied to ACACENT centroid angles → `d_aca_meas`.

It uses Ska / mica / cheta / kadi tools that the production C pipeline
does not have, so it works as an independent verification scaffold for
the LM fit that lives in `asp_solve`.

---

## Fixed in this pass

### Bugs (§E.1 of the plan)

| ID | Site | Fix |
| --- | --- | --- |
| A1 | [process_fids.py:298](process_fids.py#L298) | `d_aca_est = M @ d_fc` was wrong for `(N_fid, 3)` arrays. Replaced with `d_aca_est = d_fc @ M.T` so each row is `M @ d_fc[i]`. Comment added so the next reader doesn't re-introduce the left-multiply. |
| A2 | [process_fids.py:209](process_fids.py#L209) | `nominal_fid_positions` no longer reads `obs_si` from outer scope. Now derived from `obspar["detector"]` inside the function. |
| A4 | [process_fids.py:276](process_fids.py#L276) | Stray semicolon dropped. |

### Restored explanatory comments + dead-code removal (§E.2)

The `process_centroids` docstring now records the remaining deliberate
divergences from `aca_corr_fid`:

- **B3** — cheta `OOBAGRD3/6` minus obs-window mean instead of OBCENG
  `smoothGradient3/6` minus pre-computed `meanGradient3/6`.
- **B4** — arithmetic mean of `counts` across slots vs C's per-fid
  windowed median + cross-fid mean.

**B1** (center of expansion) and **B2** (SG-smoothed centroids) have
since been closed — see "Closed in follow-up" below.

The dead branch at [process_fids.py:118-120](process_fids.py#L118-L120) of
the previous version (B6) was removed:

```python
if np.any(np.diff(times) <= 0):
    times.sort()
```

`np.unique` always returns sorted output so the body could not run.

The polynomial-degree question in `nominal_fid_positions` ("4-degree or
5-degree?") was settled in the comment at
[process_fids.py:238-239](process_fids.py#L238-L239): 4th-order, 5
coefficients (`np.arange(5)` = exponents 0..4), per
`docs/PREDICTED_FID_AND_ASP_SOLVE.md` §6 and `Aca_Id_Image.cc:1939-1940`.

A small comment was added to the periscope ±1000 s cheta fetch to make
the padding explicit.

### Structural cleanup (§E.4)

| ID | Change |
| --- | --- |
| C1 | The `[calalign_1, calalign_2, calalign_3, calalign_4]` heterogeneous list returned by `get_caldb_calalign` is replaced with a `CalAlign` dataclass at [process_fids.py:19-39](process_fids.py#L19-L39) with fields `instr_row`, `fids`, `acis_corr`, `peri_corr`. Every call site (`process_centroids`, `nominal_fid_positions`, `estimated_fid_positions`) now uses `calalign.<field>`. |
| C5 | `CALDB_DIR` honours `$CALDB` from the environment, falling back to the prod ASCDS path: `Path(os.environ.get("CALDB", "/home/ascds/DS.release/CALDB"))`. |

The `instr_row` selection (the `idx = np.argwhere(...)` step matching on
`INSTR_ID == obspar["detector"]`) is now done once inside
`get_caldb_calalign`, removing duplicate row-matching from the callers.

---

## Verification

`python -c "import ast; ast.parse(open('process_fids.py').read())"` —
syntax OK. The script was not run end-to-end here (would need cheta /
kadi / mica configured plus a real CALDB tree). The fixes are
syntactic / structural; behavioural verification still owes a run on
the example obsid `29878`.

---

## Closed in follow-up

- **B1 — center of expansion.** `process_centroids` now computes
  `obs_y/z_coe = FID_Y/Z_CENTER_NOM + δ` per
  `tool_functions.c:615-643`, with `ang_y_nom`/`ang_z_nom` and the
  `ID_STATUS == "GOOD"` fid filter pulled directly from FIDPROPS
  ext-1 (columns `ANG_Y_NOM`, `ANG_Z_NOM`, `ID_STATUS`, `SLOT` per
  `vio_data_product.temp:644-654`). `__main__` now passes `fid_props`
  through.
- **B2 — SG-smoothed centroids.** `sg_smooth_per_slot` is now applied
  at the top of `process_centroids` (matches `aca_filter_centr`),
  so all downstream arithmetic operates on the C-equivalent of
  `ang_y_sm` / `ang_z_sm`. Reading ACACENT's persisted `ang_*_sm`
  column would have been the alternative — recomputing locally is the
  chosen path. Verification owed: confirm `sg_smooth_per_slot` matches
  `aca_filter_centr` byte-for-byte (window-length rule, edge handling).
- **B5 — periscope time-index.** `np.unique(t, return_inverse=True)`
  is used end-to-end; no `searchsorted`.
- **B8 — `d_aca_meas` unitized.** Done at
  [process_fids.py:199](process_fids.py#L199).
- **C2 — `acis_corr_params` naming.** The parameter is now
  `dh_temp_calib`; CALALIGN ext-3 is accessed via `calalign.acis_corr`
  inside the function. No more misleading name.

## Open items (deliberately not addressed)

These were on the original review (sections A.3, B.7, C.3, C.4) and
were excluded from this pass either as low-priority or as
fidelity-versus-correctness trade-offs the user wanted to defer.

### Structural items

- **C3 — `estimated_fid_positions` only computes the seed.** The
  docstring describes an iterative LM fit that the function does not
  actually do. No `solve(...)` wrapper exists; the LM residual lives
  only as a commented snippet inside `fid_process_meas.py`.
- **C4 — `__main__` doesn't compare the predicted and measured
  chains.** The script returns `d_aca_est` and `d_aca_meas` but never
  feeds them into a chi-square. Without that, the file is two parallel
  computations rather than a verification scaffold.

### Lower-priority

- **A3** — `get_caldb_align_files()` will raise on any CALDB align file
  missing `CVSD0001`; convention says this never happens, so the bare
  read was kept.
- **B7** — the ±1000 s cheta-fetch margin is now commented but not
  parameterised. Fine for the example obsid; revisit if a different
  observation hits a cheta gap.

---

## Recommended next steps

If the script is to actually drive a verification fit, in order:

1. **C3 + C4** — wrap `least_squares` (the residuals snippet at the
   bottom of `fid_process_meas.py` is ready to lift verbatim, and the
   doc gives it analytically in §4.1). Even with `(theta_X, dy, dz) =
   (0, 0, 0)` it should print a chi-square at the seed; that one
   number is what most reviewers will look at first.
2. **B4** — per-fid windowed median DH-temp. Currently a deliberate
   simplification; revisit if the seed chi-square exceeds expectations.
3. The remainder (A3, B7) is cleanup and can wait.
