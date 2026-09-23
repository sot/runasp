# Chandra DS Repository & Aspect Pipeline — Technical Overview

## 1. Repository Overview

**DS (Data Systems)** is the ground software suite for the **Chandra X-ray Observatory**, version 10.15.0. It handles the complete data lifecycle: real-time telemetry ingestion, instrument data processing, aspect (pointing) solution calculation, archive services, and observation planning tools. The repository is a large, multi-language scientific system tested on RHEL 8.

### Technology Stack

| Language | Files | Primary Role |
|----------|-------|-------------|
| C/C++ | ~7,400 | Core libraries, data processing, telemetry, aspect pipeline |
| Java | ~2,448 | Observation Cycle web apps (CDRS, XO, CPS), TAP/UWS services |
| Python | ~803 | Archive web services, testing utilities |
| Fortran | ~124 | Scientific calculations |

**Build system:** CMake 3.19+ (C++17/C11) with Ninja, Gradle for Java components, Conda/Micromamba for environment management.

**Key external dependencies:** CIAO, CFITSIO, FFTW3, GSL, Xerces-C, Sybase, ds9, XPA.

### Top-Level Structure

```
ds/
├── dslibs/          Core C/C++ libraries (time, file I/O, IPC, communication)
├── dsguilibs/       Astronomy/GUI libraries (coords, AGASC, frames, obs data)
├── dstools/         Data processing tools (ACIS, HRC, ASP, MTA, analysis)
├── ap/              Alert Pipeline — real-time flight data processing
├── telem/           Telemetry processing (ACORN, packet decompression)
├── vv/              Verification & Validation tools
├── obscycle/        Observation planning & proposal management (Java)
├── arcweb/          Archive web services (TAP, UWS)
├── arcserver/       Archive server infrastructure
├── arcclient/       Archive client tools (OCAT, proposal submission)
├── common/          Shared config and runtime setup
├── test/            Test infrastructure (BATS, pytest, CTest)
└── cmake/           CMake modules
```

### Build Workflow

```bash
conda activate ds-dev
cmake --preset conda-dev
cmake --build build
cmake --install build
```

Four presets in `CMakePresets.json`: `default`, `ninja-multi`, `dev`, `conda-dev`.

---

## 2. The Aspect Pipeline

The aspect pipeline lives in `dstools/asp/` and is responsible for one of the most critical functions in X-ray astronomy: precisely determining *where the telescope was actually pointing* at every moment during an observation. This is what makes it possible to assign sky coordinates to detected photons.

The pipeline processes **Level 0 (raw telemetry)** data from the ACA (Aspect Camera Assembly) and gyroscopes into **Level 1** calibrated aspect solutions accurate to sub-arcsecond precision.

### 2.1 Data Product Hierarchy

```
L0 Telemetry
├── TELIMG     Raw ACA image frames
├── PCADENG    PCAD engineering data (gyro, attitude, bias)
└── CCDMENG    CCD module engineering

     ↓ aca_read_data

ACADATA        L1 ACA image stack (calibrated images)

     ↓ aca_corr_ccd → aca_calc_centr → aca_corr_centr
     ↓ aca_corr_fid → aca_filter_centr → aca_calc_photom

ACACENT        Star/fid centroids with corrected positions & fluxes

     ↓ aca_id_image (+ AGASC catalog + OCAT database)

GSPROPS        Guide star properties
FIDPROPS       Fiducial light properties

     ↓ asp_solve

ASPSOL         Aspect solution: RA/Dec/Roll as function of time

     ↓ asp_make_qualint

ASPQUAL        Quality intervals (RED/YELLOW/GREEN health flags)

     ↓ asp_combine

ASPSOLOBI      Combined aspect solution per observation

Parallel gyro branch:
PCADENG → gyro_read_data → GYRODATA → gyro_corr_bias
       → gyro_process → GYRODATA (calibrated, gap-filled)
       → asp_kalman → KALMAN (smoothed attitude state)
```

---

### 2.2 Tool Catalog

#### ACA Image Processing

**aca_read_data** — Translates L0 ACA image telemetry into the L1 `ACADATA` product. Inputs: `TELIMG`, `PCADENG`, `ACACAL`. Outputs: `ACADATA`.

**aca_corr_ccd** — Applies CCD-level pixel corrections: dark current subtraction, bias removal, flat-field, and pixel responsivity. Uses dark current maps and flat fields from CALDB. Inputs: `ACADATA`, `ACACAL`. Outputs: corrected `ACADATA`.

**aca_calc_centr** — Calculates star centroids and fluxes. Supports multiple fitting algorithms: first-moment, PSF-model, Gaussian, and elliptical Gaussian. Key parameters: `algorithm`, `n_fit_par`, `n_interp` (sub-pixel), `atol` (convergence). Inputs: `ACADATA`, `ACACAL`. Outputs: `ACACENT`.

**aca_corr_centr** — Applies centroid position corrections for focal-mechanism (FM) flexure, charge transfer inefficiency (CTI), and color-dependent effects. Uses polynomial smoothing (order 0–4) over a configurable time scale. Inputs: `ACADATA`, `ACACAL` (FM_CORR, CTI_CORR, COLOR_CORR tables), `ACACENT`. Outputs: corrected `ACACENT`.

**aca_corr_fid** — Applies thermal and periscope corrections to fiducial light centroids. Configurable: `apply_thermalcorr`, `apply_pericorr`, `apply_medianfilt`, Hanning window smoothing. Inputs: `ACACENT`, `FIDPROPS`, OBC engineering. Outputs: `ACACENT`.

**aca_filter_centr** — Smooths centroid motion using Savitzky-Golay filtering with sigma rejection, separately for stars (`T_sm_star`) and fiducials (`T_sm_fid`). Inputs: `ACACENT`, `ACACAL`. Outputs: smoothed `ACACENT`.

**aca_calc_photom** — Computes photometric fluxes for ACA monitor images using Savitzky-Golay smoothing with adaptive windows tuned to a target S/N ratio. Handles 4×4 and 6×6 CCD images with corner pixel exclusion. Inputs: `ACADATA`. Outputs: `ACACENT` with flux values.

**aca_find_hotpix** — Identifies hot and warm pixels by statistical analysis across image stacks. Generates bad pixel lists for CALDB. Parameters: `percentile`, `hotthresh`, `hotfrac`, `warmthresh`, `merge`. Inputs: `ACADATA`, optional `CALBPL`. Outputs: `CALBPL`, dark current report.

**aca_make_int** — Determines PCAD mode, submode, and aspect mode as a time-dependent function from engineering telemetry. Produces the aspect interval properties product. Inputs: `PCADENG`, `CCDMENG`, `SIM`. Outputs: `AIPROPS`.

**aca_make_psf** — Creates interpolated ACA PSFs from CALDB PSF library grids, parameterized by star color (B-V) and detector position. Inputs: `ACACENT`, `ACADATA`, `ASPSOL`, `GSPROPS`, `ACAPSF`. Outputs: `ACACAL` with PSF library.

**aca_id_image** — Identifies guide stars and fiducial lights by cross-referencing observed ACA images against the OCAT (Observation Catalog) and AGASC. Applies proper motion corrections. Falls back to a text file (`ocat_data.dat`) when the OCAT database is unavailable. Inputs: `PCADENG`, `SIMCOOR`, OCAT, `CALALIGN`. Outputs: `GSPROPS`, `FIDPROPS`.

**aca_proc_dark** — Reads telemetry calibration blocks and constructs dark current calibration images. Outputs: `ACADRK`, `ACAEXP`.

#### Aspect Solution Tools

**asp_solve** — The core attitude solver. Determines the actual position of fiducial lights in the focal plane relative to their nominal positions, then computes the pointing offset (RA/Dec/Roll). Uses iterative convergence with configurable algorithm (PSF, FM, GAUSS, EGAUSS). Inputs: `ACACENT`, `FIDPROPS`, `ACACAL`, `KALMAN`. Outputs: `ASPSOL`.

**asp_calc_offsets** — Creates aspect offsets file (dy, dz, dtheta) from the aspect solution relative to nominal pointing, including SIM offsets. Inputs: `ASPSOL`, `OBSPAR`, SIM coordinates. Outputs: aspect offsets file, SIM offsets file.

**asp_calc_boresight** — Calculates boresight offsets and science instrument alignment matrices using a Nelder-Mead simplex optimizer (`amoeba()`) minimizing `D = N - R·K·M` where M=observations, N=predictions, K=scale, R=rotation (Z-Y-X Euler). Inputs: alignment data, FTS data. Outputs: updated `CALALIGN`.

**asp_obc_solve** — Calculates aspect solution directly from onboard computer (OBC) engineering data. Can output in MNC or ACA frame. Inputs: `PCADENG`, `SIMCOOR`, `CALALIGN`. Outputs: `ASPSOL`.

**asp_calc_photo** — Photometric aspect determination using guide star flux ratios. Useful for cross-checking the primary centroid-based solution. Inputs: `ACACENT`, prior `ASPSOLN`, `GSPROPS`, `GYRODATA`, `GYROCAL`. Outputs: photometric solution file.

**asp_kalman** — Applies Kalman filtering (forward and backward passes) to smooth the aspect solution, combining gyro and star tracker data. Inputs: `TKLMN`, `GYRODATA`. Outputs: `KALMAN`.

**asp_combine** — Combines CAI-level (per-aspect-interval) `ASPSOL` files into a single per-observation `ASPSOLOBI` product. Inputs: `ASPSOL` stack, `OBSPAR`. Outputs: `ASPSOLOBI`.

**asp_make_qualint** — Generates the `ASPQUAL` quality interval product by aggregating health flags across multiple input products. Flags include: image diameter/roll RMS, gyro glitches/gaps/inconsistency, fit chi-squared, ACA data quality, reaction wheel speed. Each flag assigned RED/YELLOW/GREEN based on configurable thresholds. Inputs: `ASPSOL`, `GYRODATA`, `ACACENT`, `PCADENG`, `CALRWS`. Outputs: `ASPQUAL`.

**asp_corr_props** — Corrects `GSPROPS` and `FIDPROPS` for optical distortions and temperature effects. Inputs: `GSPROPS`, `FIDPROPS`, `ACACENT`, `ASPSOL`, `ACACAL`. Outputs: corrected properties files.

**asp_vv_centroids** — Evaluates centroid quality and produces slot exclusion lists for the aspect solution. Parameters: `imprv_thresh`, `abs_rms_thresh`, `rad_tol`. Inputs: `ACACENT`, `GSPROPS`, `ASPSOL`. Outputs: slot exclusion list.

**asp_get_calib** — Extracts calibration data from CALDB and builds `ACACAL` and `GYROCAL` products. Selection is both **time-dependent** (using TSTART/TSTOP from PCADENG) and **instrument-dependent** (SIM-Z position → instrument ID). See Section 2.4 for details.

**asp_create_psflib** — Combines 2D PSF library slices into the 4D `ACAPSF` file parameterized by color and focal-plane position. Inputs: 2D PSF stack. Outputs: `ACAPSF`.

**asp_filter_startable** — Filters specific ACA slots from aspect products by slot ID list.

**asp_read_ocat** — Reads the OCAT observation database to retrieve nominal pointing and target information. Falls back to text file when DB is unavailable.

**asp_sgt2cent** — Legacy format converter: translates old `ACASIGHT` centroid data into the current `ACACENT` format.

#### Gyroscope Tools

**gyro_read_data** — Extracts L0 gyro data from PCAD engineering telemetry. Nominal sample period 0.25625 sec. Inputs: `PCADENG`, optional `GYROCAL`. Outputs: `GYRODATA`.

**gyro_corr_bias** — Applies gyro bias corrections from the `aogbias` field in PCAD engineering data. Inputs: `GYRODATA`, `PCADENG`. Outputs: bias-corrected `GYRODATA`.

**gyro_process** — Fills telemetry gaps in gyro rate data and converts raw counts to spacecraft angular rates using calibrated scale factors and alignment matrices. Parameters: sigma-edit smoothing (`sig_ed_*`), gap interpolation policy (`gap_max`, `gap_min_contig`, polynomial `order`). Inputs: `GYRODATA`, `GYROCAL`. Outputs: calibrated, gap-filled `GYRODATA`.

---

### 2.3 Shared Libraries

#### asputils/
The numerical backbone of the pipeline:
- **`quaternion.cc`** — Quaternion class for attitude representation and composition
- **`matr.cc`** — Template matrix class `Matr<double>` for linear algebra (used in boresight optimization)
- **`asp_sg_smooth()`** — Savitzky-Golay smoothing with sigma rejection
- **`field_distortion.c`** — ACA optical field distortion corrections
- **`gap_fill_poly.c`** — Polynomial interpolation for telemetry gaps
- **`aspcaldb.cc`** — CALDB interface wrapping `ciaolibs::caldb4`

#### aspnr/
Numerical Recipes routines:
- `amoeba.c` — Nelder-Mead simplex optimizer (used by boresight tool)
- `powell.c` — Powell's direction-set minimization
- `dsvdcmp.c` / `dsvbksb.c` — Singular value decomposition
- `gaussj.c`, `ludcmp.c` — Linear system solvers

#### aspect_lib/
FITS data product definitions and I/O for all aspect data types. Every tool in the pipeline reads/writes through this library. Covers 20+ data product types: `ASPSOL`, `ACACAL`, `ACACENT`, `ACADATA`, `GYROCAL`, `GYRODATA`, `GSPROPS`, `FIDPROPS`, `KALMAN`, `ASPQUAL`, and all calibration types.

---

### 2.4 Calibration Data (CALDB)

All calibration is managed through `asp_get_calib`, which queries CALDB using `ciaolibs::caldb4`. Selection is simultaneously **time-dependent** and **instrument-dependent** (based on SIM-Z position → ACIS-I, ACIS-S, HRC-I, HRC-S).

| CALDB Product | Contents | Used by |
|--------------|----------|---------|
| `CALALIGN` | ACA-to-spacecraft and FTS misalignment matrices | `asp_calc_boresight`, `asp_get_calib`, `aca_id_image` |
| `CALFDC` | Field distortion coefficients (Y/Z polynomials) | `asp_get_calib` |
| `CALCCD` | CCD characteristics and responsivity maps | `aca_corr_ccd` |
| `CALBPL` | Bad pixel lists (row/col boundary regions) | `aca_find_hotpix`, `aca_corr_ccd` |
| `ACADRK` | Dark current maps | `aca_corr_ccd` |
| `CALCTI` | Charge transfer inefficiency coefficients | `aca_corr_centr` |
| `ACISFIDCORR` | ACIS fiducial light thermal expansion | `aca_corr_fid` |
| `PERIFIDCORR` | Periscope fiducial correction (added 2008) | `aca_corr_fid` |
| `CALIRU` | IRU (Inertial Reference Unit) characteristics | `gyro_read_data` |
| `CALSFMA` | Gyro scale-factor and misalignment matrix | `gyro_process` |
| `ACAPSF` | 4D PSF library (color × position) | `aca_make_psf`, `aca_calc_centr` |

---

### 2.5 AGASC — Guide Star Catalog

The **Authorized Guide Star Catalog (AGASC)** integration lives in `dsguilibs/agasc/aspectGetAgasc.cc`. The `agascStar` class provides FITS-based access (via CFITSIO) to the catalog, keyed by the `ASCDS_AGASC` environment variable.

Each catalog entry provides 30+ properties per star:
- Sky coordinates (RA, Dec) and proper motion
- ACA-band magnitude and B-V color
- Variability flags
- Quality indicators: `ASPQ1`–`ASPQ3` (aspect quality), `ACQQ1`–`ACQQ6` (acquisition quality)
- Cross-reference IDs (Tycho, HD, etc.)

This catalog is the ground truth for `aca_id_image` — without it, the pipeline cannot associate ACA image slots with known stars.

---

### 2.6 Supporting Libraries from dsguilibs

**dsguilibs/coords/** — Fundamental coordinate transformations:
- `aca2sky.h` — ACA focal-plane mm → sky RA/Dec/Roll
- `tquat.h` — Quaternion → RA/Dec/Roll
- `offset2radec.h` — Offset angle → sky coordinate conversions

**dsguilibs/obs/** — Observatory metadata:
- `Boresight.cc` — Boresight geometry calculations
- `calc_nomroll.cc` — Nominal roll angle computation
- `roll2window.cc` — Detector roll → window mapping

**dsguilibs/frame/** — Application framework:
- `FW_Application.hh` — Base class for all C++ aspect tools (parameter handling, signal management)
- `FW_Parameter.hh` — Parameter binding and validation

**dslibs/xtime/** (`XTime.hh`) — Time system conversions between MET, TT, TAI, UTC. Used throughout the pipeline for precise TSTART/TSTOP filtering. Maintains built-in leap second tables.

---

### 2.7 Test Infrastructure

The asp tools use a **two-tier testing approach**:

**BATS (modern):** `wrap_*.bats` files in each tool directory. Uses BATS assertion library (`assert_output`, `assert_success`). Tagged with tool and library dependencies for CI filtering.

**Perl regression tests (legacy):** `*.t.in` templates configured by CMake into `*.t` scripts. Compare current output against saved baselines. Test data managed via `.lis` stack files. Use `TESTIN`/`TESTOUT`/`ASCDS_VERSION` environment variables.

**Unit tests:** Python `unittest_*.py` modules and C/C++ `test_*.cc` executables for targeted algorithm testing (e.g., `asputils/test_quat.cc` for quaternion algebra).

**CTest integration:** All test programs built via `cxcds_add_executable()` in CMakeLists.txt files and registered with CTest.

---

### 2.8 Interesting Findings

- **Photometric fall-back:** `asp_calc_photo` implements a full photometric attitude solution (from star flux ratios alone) as an alternative to the primary centroid-based method — useful for anomaly investigation when centroid tracking is degraded.

- **Periscope correction added 2008:** The `PERIFIDCORR` calibration product was added to `aca_corr_fid` in 2008, reflecting in-flight discovery that periscope thermal flexure affects fiducial light positions beyond what the earlier model accounted for.

- **Nelder-Mead boresight:** `asp_calc_boresight` uses a Nelder-Mead simplex minimizer from Numerical Recipes to solve a non-linear rotation matrix fit — the boresight can optionally float a scale factor K, or fix it at 1.0.

- **OCAT fallback:** `aca_id_image` has a database fallback path that reads `ocat_data.dat` as a plain text file when the Sybase OCAT database is unreachable — ensuring the pipeline can run offline or in test environments.

- **Gyro sample rate:** The gyroscope is sampled at exactly 0.25625 seconds — an oddly precise number driven by hardware timing constraints of the IRU.

- **Quality flags propagate forward:** `ASPQUAL` flags from `asp_make_qualint` are the primary mechanism by which downstream science processing knows which time intervals are safe to use. A RED flag effectively invalidates photon event data for that period.
