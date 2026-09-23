# Understanding `_asol1` pipeline products

## Context

Building a shared mental model of the aspect pipeline before drafting a
per-step pipeline report. This document answers one question the user
raised about specific files they have in hand:

1. What produced `pcadf893066380N001_asol1.fits.gz` vs `pcadf29878_000N001_asol1.fits.gz`, and how they differ.

Everything below is grounded in the pipeline shell scripts in
[pipelines/asp/](../pipelines/asp/) and the tool
source under [dstools/asp/](../dstools/asp/).

---

## Part 1 — The two `_asol1` files

The filename encodes which orchestrator wrote the file. Both carry
"aspect solution" data but at different scopes:

| File | Written by | Pipeline | Scope | Filename fields |
|---|---|---|---|---|
| `pcadf893066380N001_asol1.fits.gz` | `asp_solve` | [asp_l1_std.ped:346](../pipelines/asp/asp_l1_std.ped#L346) | **CAI** (one aspect interval) | `893066380` = TSTART in CXC MET seconds |
| `pcadf29878_000N001_asol1.fits.gz` | `asp_combine` | [asp_l1_combine.ped:24](../pipelines/asp/asp_l1_combine.ped#L24) | **OBI** (one observation) | `29878_000` = OBSID `29878`, OBI index `000` |

`N001` in both names is the processing-version field.

### How they differ in content

`asp_combine` concatenates all per-interval CAI `_asol1` files for an
observation into a single OBI product, then trims intermediate columns
and rewrites identifying header keywords. The differences that matter:

- **CONTENT keyword** — the single most reliable discriminator.
  - CAI: `CONTENT = 'ASPSOL'`
  - OBI: `CONTENT = 'ASPSOLOBI'` — set at [asp_combine.c:333](../dstools/asp/asp_combine/asp_combine.c#L333).

- **Trimmed columns in OBI.** `asp_combine` drops the raw/intermediate
  columns `ADY, ADZ, ADTHETA, RA_RAW, DEC_RAW, ROLL_RAW, Q_ATT_RAW`
  ([asp_combine.c:174](../dstools/asp/asp_combine/asp_combine.c#L174))
  — the CAI file keeps them; the OBI file doesn't.

- **OBI-only pointing keywords.** `asp_combine` adds `RA_NOM / DEC_NOM / ROLL_NOM`
  (nominal pointing, copied from OBSPAR) and `RA_PNT / DEC_PNT / ROLL_PNT`
  (actual mean pointing computed from the combined rows), and
  `CAIFILE1 …` keywords listing the input CAI files.

- **HDUCLAS\* keywords.** Only `asp_combine` writes `HDUCLASS='ASC'`,
  `HDUCLAS1='TEMPORALDATA'`, `HDUCLAS2='ASPSOL'`
  ([asp_combine.c:335-338](../dstools/asp/asp_combine/asp_combine.c#L335-L338)).

In short: the `893066380` file is the raw per-interval solver output,
and the `29878_000` file is the per-observation merged+trimmed product
delivered to users.

