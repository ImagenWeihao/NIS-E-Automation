# CHANGELOG

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
