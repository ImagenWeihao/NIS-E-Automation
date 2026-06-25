# CHANGELOG

## v1.7.2 — 2026-06-25

- **Merged Spheroid State + Dashboard table** — the Step-3 side-panel table and the dashboard's
  status table are now one Treeview ("Spheroid State & Dashboard"), placed where the state table
  was, expanded to full width by default. Columns: **Use** (checkbox), Rank, Spheroid ID, Status,
  Mosaic XY (µm), Verified XY (µm), Z-centre (µm), Diam (µm), Score, Bin file. Styled to match the
  matplotlib dashboard (dark `bg2` field, lavender monospace headings, per-status row colors reused
  from `_ROW_BG` / `_STATUS_STYLE`). The Use toggle and Step-3 trigger/recenter exclusions are
  unchanged.
- **Dashboard is now map-only** — `PipelineDashboard`'s redundant status-table panel was removed
  (`_ax_tbl` / `_draw_table` deleted) and the 10× mosaic preview map expanded to fill the whole
  figure. The per-spheroid data now lives solely in the merged table above. Affects the GUI live
  dashboard and the standalone CLI dashboard figure alike.
- **Captured Z-Stacks moved to its own tab** — the captured-spheroid Z-stack viewer is now a
  dedicated notebook tab between **PA Workflow** and **Log**, instead of a pane inside the workflow.
  The PA Workflow tab's paned window is now two panes (steps+dashboard | merged table), so the
  table gets even more width by default.
- **Beer-Lambert compensation is now opt-in** — Step 4 gains a **"Beer-Lambert Intensity
  Compensation"** checkbox (default **OFF**). ON → `Generate All Bins` writes the depth-adaptive
  `P0*exp(depth/L)` ramp (P0 and L active). OFF → flat bins at the fixed **P0 (%)** power (default
  15%), implemented as `L → inf` so every plane gets the same power (only L is disabled; P0 stays
  editable as the base power in both modes). Bins are still always generated so the Z-corrected
  capture works in both modes; the status line / log report which mode was used.

## v1.7.1 — 2026-06-24

- **Step 3 per-spheroid selection** — the Spheroid State table gains a leading **Use**
  checkbox column (`[x]`/`[ ]`, default checked). Click a row's Use cell to include/exclude
  that spheroid. **Trigger NIS-E Z-Stack Captures** now writes `af_trigger` files only for the
  checked spheroids (unchecked are skipped, reported in the status line); **Re-center from
  captures** also honors the selection. Implements the Step-2/3 "remove unwanted spheroid IDs"
  TODO — previously the only selection was the Step-1 "Top N" count.

## v1.7.0 — 2026-06-23

- **Multi-spheroid PA macros:**
  - `nis_macro_pa_points.mac` — builds an N-point ND multipoint from the first N
    `af_trigger` files; the job's "Import Point Set from ND" then sweeps all N in one run.
  - `nis_macro_pa_validate.mac` — 1050 nm post-PA re-image (zoom 2.5, Z-stack at each
    trigger's focus Z, PFS dichroic OUT — IR imaging is incompatible with the PFS dichroic
    in path) to confirm the faded squares.
  - `nis_macro_pa_pick_current.mac` — grabs the current centred spheroid (`StgGetPos`)
    into `af_trigger_01` + a 1-point multipoint, for when the plate has moved.
- Step 2: TODO to add removal of unwanted spheroid IDs before anchoring.
- `.gitignore`: also ignore dated run folders with suffixes (e.g. `0623_2/`).

## v1.6.0 — 2026-06-23

- **Step 4 Photoactivation panel** — Job / Activation OC / Power% / Well / Loops / Zoom /
  Dichroic OUT / Remove A1 interlock inputs. "Run Photoactivation Job" writes `pa_trigger.ini`
  with PA output pinned to `<Step-1 base>/pa` (not the NIS-E default save location).
- **`nis_macro_pa_setup.mac`** (NEW) — reads `pa_trigger.ini` and preps the rig for the
  `step3_zstack_PA` JOB: `Stg_RemoveInterlock` (clear A1 interlock), `SelectOptConf` (850 nm OC),
  `Stg_PFSInsertExtractDM(0)` (dichroic OUT), optional `Stg_SetMultiLaserPower`, and
  `Confocal_SetScanArea` (centred zoom square). Validated on the rig for one spheroid.

## v1.5.1 — 2026-06-22

- **Generate All Bins falls back to the Step-3 Middle plane Z** when a record has no
  recorded z_centre (Refresh Status not run) -- the spheroids are captured centred on
  that plane, so it is the correct bin z-centre; no more silent 0/6.
- GUI header/title version label bumped to v1.5.1 (was stale at v1.2).
- Added `focus_scores()` library helper (variance / tenengrad / laplacian / brenner / fft).
  Graded all WellD05 stacks: every metric snaps to a stack *edge*, not the focus centre
  (the +/-25 um range over a ~250 um spheroid has no sharpness peak), so it is deliberately
  NOT wired into the GUI -- the PFS/home plane stays the best-focus reference.

## v1.5.0 — 2026-06-22

### Coordinate registration: flip-aware Step 2 + per-spheroid re-center (WellD05, 20X)
- **Flip-aware anchor transform** (`apply_anchor_transform`): the mosaic<->stage map is a
  ~180 deg axis flip, which a single averaged translation can't represent (only the matched
  spheroid landed). Step 2 now fits a `SimilarityTransform` from >=2 anchor correspondences
  (captures the flip+scale) -> 6/6 spheroids in frame (was 2/6). With 1 anchor it falls back to
  the known 180 deg flip + that anchor's translation (`flip1`) so even one match places all six.
- **Per-spheroid re-center from captures** (`recenter_from_captures` + GUI "Re-center from
  captures"): measures each spheroid's intensity centroid offset in its own capture and nudges
  its stage XY so the next pass lands centred -- exact per-spheroid, independent of anchor
  quality; corrects the systematic camera-centre vs stage offset the 2-anchor fit can't see.
  Flip-X/Y axis-sign toggles, "# to use" count (default all), "re-alignment in process" status,
  per-spheroid delta shown + re-centred rows highlighted in the Spheroid State table.

### Step 4: A1 confocal Z-intensity-corrected capture
- New daemon `nis_macro_capture_zcorrected.mac`: per spheroid loads its Z-correction `.bin`
  (`ND_ZIntensityControlLoad`), confirms `ND_ZIntensityControlIsDataReady`, runs
  `ND_RunZSeriesExpWithZIntensityCorrection`, saves, writes done. Reads the trigger fully then
  deletes it before capturing so the sequential queue can't be clobbered. Bins generate from a
  rig reference `.bin` (Beer-Lambert CH2 laser ramp). Z-correction is A1-confocal-ONLY -- on the
  Flash 4.0 camera `IsDataReady` returns `zcorr_not_ready` (incompatible detector).

### Workflow + viewer
- **session.ini pointer file**: daemons read `work_dir`/`nd2_dir` from `C:/SpheroidPA/session.ini`
  (written by Step 3 `trigger_autofocus_all`) -> follow the Step 1 save dir, no hardcoded paths.
- **Auto Load + filmstrip list**: one click loads all captured stacks as a vertical list of
  filmstrip rows (per spheroid), newest first, auto-locating the nd2 folder.
- **1-plane "Locate only" mode**: Step 3 writes `z_half=0` triggers; the capture macro then does a
  single `Capture()` (fast locate pass for re-centering) instead of the full Z-series.
- **Capture macro stops after the batch** (no 300 s lingering that grabbed the next pass's
  triggers); Step 3 Z-step default 5 um for 20X.
- **Removed "Apply Global Registration"** (unreliable matching 20X Z-stacks against the 10X mosaic).

## v1.4.3 — 2026-06-17

### Captured-spheroid Z-stack image viewer (its own pane)
The PA Workflow tab is now a three-pane layout: **steps + live dashboard | Z-stack viewer |
Spheroid State table** (all draggable). The dashboard stays visible at all times; the viewer
is a permanent pane between the dashboard and the table (not a Step-3 swap).
- **Spheroid selector** (dropdown) auto-scans `<out>/nd2`, `<out>/autofocus`, and the Step 4
  ND2 dir for `*.nd2`; **Refresh** rescans; **Add...** opens a file dialog to load any ND2
  manually. Reads each stack with the `nd2` lib (handles (Z,Y,X) and singleton/extra axes).
- **All planes shown at once** as a horizontally-scrollable **filmstrip** (Pillow thumbnails) —
  one image per Z plane, each captioned with its **true Z depth (µm)**. The **focus-centre plane
  is highlighted** (= the autofocus/PFS home plane, which is the best focus). Mouse-wheel scrolls.
- **Per-plane Z is read from the ND2 events log** (`Z Coord [µm]` / `Ti ZDrive [µm]`), the real
  scanned focus position. Earlier audit wrongly concluded the Z wasn't recorded — it was, in the
  events log; the field that read static (`stagePositionUm.z`) is just the coarse-stage snapshot
  (constant across the stack). The focus/centre plane is identified from the events `Z-Series == 0`.
  If a file genuinely has no per-plane events Z, the viewer falls back to reconstructing the grid
  from the GUI Middle-plane Z + the loop `stepUm` (or GUI Z step). A **[check]** line validates the
  plane count against ±half/step.
- **No focus score / auto best-focus pick.** Tried gradient-energy, Tenengrad, Laplacian, and an
  FFT high-frequency metric (with background removal) — none reliably matched the eye across stacks
  (e.g. on cell #6 they flagged the stack edges while the true best focus is the centre/home plane).
  Automatic grading was removed; the **focus-centre (PFS/home) plane is the best focus** and is
  highlighted, and the operator can read the filmstrip to confirm.
- Summary shows spheroid, centre-plane Z + source (ND2 events vs GUI), plane count + step,
  dimensions, and the geometry check. Thumbnails use 1–99.5 percentile contrast; load is off-thread.
- Default window 1280x820 → 1680x880 to fit the third pane and the filmstrip.

### New test macro: `macro_selftest/zstack_pfs_zcentre_10x.mac` (PFS-centred Z-stack at 10X)
Based on `zstack_cell06.mac`. Uses the **Perfect Focus System** to find the focused Z (the
spheroid z-centre) instead of the unreliable image autofocus:
`Stg_IsPFSPresent` → `Stg_SetPFSStatus(1)` → `Stg_WaitForPFS(8)` → `Stg_GetPFSStatus()==1` →
read Z → `Stg_SetPFSStatus(0)`, then a symmetric Z-series around that Z at **10 µm** steps
(`ND_SetZSeriesExp` type 0 → `ND_RunZSeriesExp`) and `ImageSaveAs`. Single `main()`, CRLF +
ASCII; every NIS call verified against the AR macro reference. Operator sets the PFS offset to
the spheroid plane before running.
- The z-centre is read from `StgGetPosZ(&z, 0)` — the **primary Z device (Ti ZDrive)**, the axis
  NIS-E logs per plane as `Z Coord`/`Ti ZDrive` — so the recorded per-plane Z is correct. (Not
  `stagePositionUm.z`, the constant coarse stage.) `StgGetAbsPosZ` is shown alongside in the
  confirm dialog so the rig run reveals which reading matches the captured Z Coord.

## v1.4.2 — 2026-06-17

### Step 3 Z-stack geometry is now GUI-driven (was hardcoded in the macro)
`nis_macro_capture_zstack.mac` previously hardcoded `z_centre=7680, z_half=90, z_step=10`
(19 planes) and read only XY from the trigger, so the GUI's Step 3 Z fields had no effect on
the live stack. Now:
- **GUI Step 3** exposes **Middle plane Z (um)** (new), **Z half-range (um)**, **Z step (um)**
  (defaults 7680 / 90 / 10 matching the rig). `trigger_autofocus_all(z_centre, z_half, z_step)`
  writes all three into every `af_trigger_NN.ini`.
- **The macro** reads `z_centre`/`z_half`/`z_step` per trigger (`Int_GetKeyString`+`atof`),
  computes the stack, and falls back to 7680/90/10 if a key is absent. Still a single flat
  `main()`, CRLF + ASCII.
- Trigger button relabeled "Trigger NIS-E Z-Stack Captures"; Step 3 hint updated to name
  `nis_macro_capture_zstack.mac`.

### GUI layout: Spheroid State table no longer overlaps the step content
- The step content + dashboard and the Spheroid State table are now in a horizontal
  **paned window** (draggable divider) instead of butting frames, so they can never overlap
  — previously, at Steps 3/4 the long status text sat flush against the table.
- The table pane is wider and its six columns are sized to all show by default
  (no horizontal scroll needed at the default window size); drag the divider for more.
- Default window 980x740 → 1280x820; long Step 3/4 status labels now wrap (`wraplength`)
  instead of running into the divider. The step pane keeps >= 560 px before the table shrinks.

## v1.4.1 — 2026-06-17

### Trigger safeguards: clear-stale + verify-distinct (root-cause of duplicated capture)
A capture run imaged duplicated coordinates (`sph_03 ≡ sph_01`, `sph_04 ≡ sph_02`). The
pipeline was proven correct three ways — the surviving original triggers are all distinct,
re-running `screen → apply_offset → trigger_autofocus_all` is all distinct, and writing to
the share + reading back through the Win32 profile API (the macro's read path) is all
distinct. **Root cause:** the daemon consumed a work folder holding a *mix* of fresh and
stale `af_trigger` files from several re-trigger cycles. `trigger_autofocus_all` now:
- **verify-distinct** — raises if any two records share `(stage_x, stage_y)` *before*
  writing or touching the instrument (a duplicate means a corrupted records list);
- **clear-stale** — removes existing `af_trigger_*`/`af_done_*` first, so the daemon only
  ever sees the current run's triggers.

Also added `macro_selftest/test_spacedir.mac` — probes whether the macro can use a spaced
work-dir read at runtime from a fixed config (toward pointing the work dir at the
ProgramData Jobs folder). Pending rig test.

## v1.4.0 — 2026-06-16

### First hardware run: daemons flattened to single main(); capture path validated
On the live instrument (WellD09, 12 spheroids) this NIS-E build was found to run **only a
single `main()`** — user-defined functions/procedures (even a parameterless `int bump()`)
raise "Cannot Evaluate the Expression". Both daemons were rebuilt as one flattened `main()`
(inline path building, `while` loops, nested `if`, status-flag pattern); the flattened
autofocus daemon then ran end-to-end on the rig.

### `nis_macro_z_autofocus.mac` — flattened + autofocus tuning
- Single `main()`; inline work-dir literal; `Int_GetKeyString`+`atof`; quoted-comma `sprintf`.
- Pre-positions Z with `StgMoveZ`, then `StgFocusInRangeTwoPasses` (fluorescence/noise-resistant
  criterion, tunable range/step).
- **Note:** the built-in autofocus fails on the deliberately-low fluorescence signal (returns
  −3, "total white or black scene"); exposure can't be raised (phototoxicity). Autofocus is
  shelved in favor of the Z-stack capture below + a future software focus-plane selector.

### Added `nis_macro_capture_zstack.mac` — Z-stack per spheroid (no autofocus)
Reads each `af_trigger_NN.ini`'s XY, `StgMoveXY`, runs a fixed Z-stack (7590–7770 µm, 10 µm,
19 planes), saves `<id>_zstack.nd2`, writes `af_done`, deletes the trigger. Flattened single
`main()`; reuses the `ND_SetZSeriesExp` → `ND_RunZSeriesExp` → `ImageSaveAs` path validated on
cell #6 (focus blurry→sharp→blurry on `Ti ZDrive`).

### Removed `nis_macro_auto_capture.mac`
The Z-intensity-corrected capture daemon (multi-function; correction is A1/confocal-only) is
superseded by `nis_macro_capture_zstack.mac` on this widefield rig.

### Self-test macros
- `test_func.mac` / `test_proc.mac` — confirm user-defined functions and parameterless
  procedures are unsupported on this build (both fail at the definition).
- `zstack_cell06.mac` — single-spheroid Z-stack test that validated the capture path.

### Deferred
ND2 per-frame `stagePositionUm.z` is a coarse-stage snapshot (constant ~7500), not the focus
depth; true per-plane Z is in the Z-stack loop (`ZHome` + `homeIndex`/`stepUm`). `_read_nd2_z`
and the GUI metadata reader to be patched.

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
