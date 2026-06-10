# CHANGELOG

## v1.3.3 — 2026-06-10

### Fix: daemons couldn't use the hyphenated work path via #define
The v1.3.2 redirect put `NIS-E-Automation` in the path, and this NIS-E build's `#define`
preprocessor *evaluates* its value — so `#define WORK_DIR ".../NIS-E-Automation/work"`
raised "Cannot Evaluate the Expression" (the `-` parsed as subtraction). Both daemons now
inline the work dir as a **direct string literal** in their path builders (direct literals
are not evaluated; only `#define` values are). `00_io_inifile.mac` likewise. Confirmed on
the rig: test 00 passes with the share path. (04/05/06 and the GUI already used direct
literals / Python strings, so were unaffected.)

## v1.3.2 — 2026-06-10

### Work directory redirected to the project share
Moved the file-flag work dir from `C:\SpheroidPA\work` to
`S:\Images\Weihao\NISeA\NIS-E-Automation\work` so the GUI and NIS-E can share it across
machines. Updated both daemons' `#define WORK_DIR`, the GUI Step 1/4 path defaults, the
self-test macros (`00`/`04`/`05`/`06`), and `verify_trigger_bridge.py`; added `work/` to
`.gitignore`.
- **Validate on a network share:** run `macro_selftest/00_io_inifile.mac` after the
  redirect — `GetPrivateProfileString` (the INI reader behind `Int_GetKeyString`) can cache
  less promptly on UNC/mapped-drive paths, so confirm the round-trip before a full run.

## v1.3.1 — 2026-06-10

### Z Intensity Correction .bin format RESOLVED (the long-standing gating unknown)
On-rig testing settled it: the panel's `Load…` filters for **Binary Files (`*.bin`)**, so
`.bin` *is* the native load format (the macro doc's ".xml" wording is generic). Our
generated bins were rejected as **"Incompatible Z Correction and Camera"** because they
carried empty metadata and the wrong device field. A known-good rig export decoded to:
- metadata `hwUnit_Name=NikonA1Grabber`, `DetectorType=10`, `ChannelBits=2`
- per item: `CH2PMTHighVoltage=5` (HV2) + `CH2LaserPower` (LP2 ramp), `iShowFlags=521`, `eItemType=3`

The codec is correct — a parse→rebuild of the reference is **byte-identical** (27,552 B).

### `generate_bin` now templates off a rig reference .bin (`spheroid_pipeline.py`)
- New `reference_bin=` parameter + `load_bin_template()` / `_detect_lp_field()` helpers.
  The generated bin inherits the reference's detector/camera metadata and per-item HV/flags,
  substituting only `ZStack` (per-spheroid Z grid) and the Beer-Lambert `CH2LaserPower` ramp.
  Item count follows the reference (sidesteps the 19-item `ROOT_TAIL_U64S` constraint).
- Legacy bare-metadata path retained as a fallback (NIS-E flags it incompatible).

### GUI (`spheroid_pa_gui.py`)
- Step 4 gains a **Reference .bin** file field (passed to `generate_bin`); warns if unset.
- Default laser field changed `CH1LaserPower` → `CH2LaserPower` (LP2).

### NIS-E daemon (`nis_macro_auto_capture.mac`)
- `#define CHECK_ZREADY 0` — `ND_ZIntensityControlIsDataReady` is A1/multilaser-only and
  returns -9 (n/a) on this rig, so the daemon no longer gates capture on it.
- `#define USE_ZCORR_RUN` switch: `ND_RunZSeriesExpWithZIntensityCorrection()` raised
  "Cannot Evaluate the Expression" on a non-A1 camera; the daemon can fall back to plain
  `ND_RunZSeriesExp()` (which runs and saves) for rigs/sessions without the A1 detector.

### Self-test macros (`macro_selftest/`) — full ladder validated on the instrument
- Build paths with `strcpy`/`strcat` string literals — a `#define` holding a dotted path
  (`…nd2`) raises "Cannot Evaluate the Expression". `05` opens the panel via `_ND_ZIntensityControl()`;
  `06` uses plain `ND_RunZSeriesExp()`.
- Rig verdicts: `00` file-flag bridge **PASS**, `04` capture/save **PASS**, `05` Z-correction
  calls respond + templated `.bin` loads, `06` Z-series runs and saves **PASS** (after taking the
  Z drive out of Escape mode). The corrected-run + `IsDataReady` are A1/confocal-only.

## v1.3 — 2026-06-10

### NIS-E macros validated on the live instrument (two dialect bugs fixed)
Ran isolated probe macros on the real NIS-E and read the interpreter errors. Both
production daemons had latent bugs that would have failed on the first real run
(the macros had never actually been executed — STATE was `not_started`):
- **`Int_GetKeyValue(...)` raised "Mismatch in dimensions"** despite matching the
  documented signature. Replaced all 7 numeric-key reads (`nis_macro_auto_capture.mac`
  ×5: stage_x/y, z_centre, z_half, z_step; `nis_macro_z_autofocus.mac` ×2: stage_x/y)
  with `Int_GetKeyString` + `atof()`.
- **`sprintf` is not C-variadic** — signature is `sprintf(buf, fmt, args)` and multiple
  args raise "Bad Number of Parameters". Args must be **one quoted comma-separated string
  of variable names**, e.g. `sprintf(buf, "%s %.2f", "name,value")`. Fixed every multi-arg
  `sprintf` (and pre-compute expressions/#defines into variables first).
- Confirmed `.mac` files must be **CRLF** (an LF-only file makes the first `//` comment
  swallow the whole file and it silently no-ops) and pure ASCII.
- The documented XY-position reader is **`StgGetPosXY(&x, &y)`**, not `StgGetPos`
  (which gives "Bad Number of Parameters"); the Z reader is `StgGetAbsPosZ(&z)`.

### Added — `macro_selftest/` on-instrument probe ladder
Seven isolated, `WaitText`-instrumented test macros (`00_io_inifile` … `06_zseries`) to
confirm NIS-E actually responds to every call the daemons make, **before** a real
experiment. Run in order; a silent/partial file pinpoints the exact call this rig rejects.
Destructive calls are made safe (XY moves to current position; tiny Z ranges). `05` is the
key check: whether `ND_ZIntensityControlLoad` accepts the pipeline's `.bin`. See its README.

### Added — `verify_trigger_bridge.py` (offline bridge proof)
Round-trips real trigger files (from the actual `spheroid_pipeline` writers) through the
Win32 `GetPrivateProfileString` API — the exact call `Int_GetKeyString`/`Int_GetKeyValue`
wrap inside NIS-E — in both directions. `RESULT: PASS` means a live macro reads identical
values. Confirmed all keys round-trip.

### Changed — robust trigger writers (`spheroid_pipeline.py`)
`write_trigger` and `trigger_autofocus_all` now write via `_atomic_write_crlf`: explicit
CRLF (not OS-dependent text-mode translation) and an `os.replace` rename so the polling
daemon never reads a half-written file.

### Removed — `nis_macro_20x_capture.mac`
Hallucinated/orphaned: built on ~10 nonexistent functions (`Stg_MoveXY`, `NdCapture`,
`Doc_SaveAs`, `InputBox`, `FileReadLine`, `StrSplit`, `ND_SetZSeriesRange`, …) and not
invoked by any Python. The real bridge is the two trigger daemons. Cleaned up the dangling
references in `CLAUDE.md` and `cross_zoom_register.py`.

### Added — `END_TO_END_TEST.md`
Full operator guide: one-time setup, a hardware-free rehearsal, and the GUI Step 1–4 ⇄
NIS-E run, with success criteria and the silent-failure troubleshooting table.

## v1.2 — 2026-06-10

### NIS-E macros — full rewrite against the official API
Both macros were previously written with **invented function names** (`Stg_MoveXY`, `NdCapture`, `Doc_SaveAs`, `ZIntCorrect_LoadFile`, `FileReadLine`, `AFC_FullRangeAutoFocus`, …) that do not exist in NIS-Elements. Rewritten and verified against the local NIS-E AR 5.42 macro reference (`Docs/nis/eng_ar`). Every NIS call now resolves to a documented function with the correct signature:
- Stage: `StgMoveXY(x,y,relative)`, `StgGetAbsPosZ(&z)`, `StgXY_GetSettleTime()`
- Autofocus: `StgFocusSetCriterion()` + `StgFocusInRangeTwoPasses(range,coarse,fine)` — single call does the coarse+fine two-pass to the spheroid midplane
- Z-series + correction: `EnableZIntensityControl()`, `ND_ZIntensityControlLoad()`, `ND_ZIntensityControlIsDataReady()`, `ND_SetZSeriesExp(...)`, `ND_RunZSeriesExpWithZIntensityCorrection()`
- Capture/save: `Capture()`, `ImageSaveAs(path,14,0)`
- INI/file: `Int_GetKeyString/Int_GetKeyValue/Int_SetKeyString`, `ExistFile`, `DeleteFile`; status via `SetCommandText()` (non-blocking)
- Correct C-dialect: `int main()` entry, top-of-body declarations, bracketed logical expressions (interpreter binds `&&`/`||`/`!` right-to-left), `long` rank counter
- Added real error handling: stage-move and Z-intensity-ready return codes are checked; `status` (ok / af_failed / move_failed / zcorr_not_ready) is written to the done file instead of always "ok"; failed autofocus skips capture so no bad Z is recorded.

### Protocol (`spheroid_pipeline.py`, `spheroid_pa_gui.py`)
- Trigger/done INI files now carry a `[spheroid]` section header (required by `Int_GetKeyString`); `parse_ini` already tolerates it, so Python round-trips unchanged.
- `poll_autofocus_done` skips ranks whose macro `status` != ok; capture queue marks a rank FAILED when the macro reports a non-ok status.

### Open hardware-verify item
- `ND_ZIntensityControlLoad` is documented to load **.xml**; the pipeline emits **.bin** (the Z Intensity Correction panel's native format). Flagged in the macro header — confirm on the real NIS-E whether the panel's Load accepts the .bin path, else have `csv_to_nis_bin.py` emit the accepted format or switch to `ND_ReuseZIntensityCorr` against a reference .nd2.

## v1.1.1 — 2026-06-10

### Pipeline GUI (`spheroid_pa_gui.py`)
- Removed the Job A–C workflow entirely: Configuration / All Spheroids / PA-Ready tabs and the "Job A → Bridge → Job B → Job C" status strip are gone. App now shows only the PA Workflow (Step 1–4) and Log tabs.
- **Step 2 anchor rank box is now editable** (was read-only). Type the rank of the spheroid you navigated to; leave blank to auto-match by NCC.
- Fixed the misleading rank display: the box previously showed a *random* `pick_anchors` suggestion (the "2") while NCC silently matched a different spheroid (the "9"). XY was always correct; only the label was wrong. After Verify, the box now updates to the **actually matched** rank, and flags a mismatch if it differs from what you typed.
- Anchor defaults now seed from the top-ranked spheroids (deterministic) instead of random picks.

### Backend (`spheroid_pipeline.py`)
- `estimate_offset_from_nd2()` now takes optional `expected_rank` and returns `(dx, dy, ncc, matched_rank)`; when a rank is given it prefers that rank's NCC match. `global_registration()` passes each capture's known rank for a constrained match.

---

## v1.1 — 2026-06-10

### Pipeline GUI (`spheroid_pa_gui.py`)
- Redesigned Pipeline tab: left sidebar with Step 1–4 cards + step-progress status dots
- Step content panels replace the old scrolling vertical list; bottom panel always shows state table + dashboard
- **Step 3 fully automated**: triggers NIS-E `nis_macro_z_autofocus.mac` per spheroid, polls done files, applies global registration across all 20X ND2s before Step 4
- Step 2: Anchor 1 now labelled "required", Anchor 2 "optional" (already skipped gracefully)
- Step 1: added editable "Top N spheroids" field (blank = all)
- Defaults corrected: `well_id` cleared (was "WellD05"), `z_half` set to 18.0 µm (was 30.0)
- Window title updated to "SpheroidPA v1.1"

### Backend (`spheroid_pipeline.py`)
- `trigger_autofocus_all()` — writes `af_trigger_XX.ini` per rank for NIS-E macro
- `poll_autofocus_done()` — reads `af_done_XX.ini` files, returns z_centre + nd2_path per rank
- `apply_autofocus_results()` — stores z_centre into records, advances status to Z_KNOWN
- `global_registration()` — NCC-matches each 20X ND2 vs mosaic, fits similarity transform (N≥2) or pure translation (N=1), updates verified_x/y on all records

### NIS-E Macros
- **New**: `nis_macro_z_autofocus.mac` — Step 3 daemon: polls `af_trigger_XX.ini`, moves stage XY, runs coarse (`AFC_FullRangeAutoFocus`) then fine (`AFC_AutoFocus`) pass to land on spheroid midplane, reads Z, captures single-plane ND2, writes `af_done_XX.ini`
- `nis_macro_auto_capture.mac` — Step 4 daemon, unchanged

### Binary format (`csv_to_nis_bin.py` / `nis_bin_to_csv.py`)
- Root tail and size formula corrected for exact byte-identical round-trip at 19 items
- Display precision matches NIS-E panel: Z → 2 dp, laser power → 1 dp

---

## v1.0 — 2026-06-08 (initial)
- Spheroid screener, NCC anchor offset, manual Z recording, Beer-Lambert bin generation, NIS-E file-flag capture daemon
