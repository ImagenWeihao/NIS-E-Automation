# Claude Code — SpheroidPA Project

## What this project is
Automated spheroid screening and cross-zoom registration pipeline for NIS-Elements microscopy.
Sits between NIS-E acquisition and 20X capture: detects, ranks, and validates spheroid locations.

## Repo layout (restructured 2026-07-01 / v1.8.2 — see MD/CHANGELOG.md)
- `Codes/GUI/`  — operational pipeline + Tkinter app: `spheroid_pa_gui.py`, `spheroid_pipeline.py`, `spheroid_screener.py`, `cross_zoom_v2.py`, `cross_zoom_register.py`, `csv_to_nis_bin.py`, `nis_bin_to_csv.py`
- `Codes/test/`     — harnesses: `dry_run_pipeline.py`, `dry_run_dashboard.py`, `verify_trigger_bridge.py`, `sam_nise_capture_test.py`
- `Codes/validate/` — QC/overlay: `verify_xy_offset.py`, `make_overlay.py`, `make_sub10x_overlay.py`
- `macro/` (dispatcher + 7 step macros; see `macro/README.md`), `macro_selftest/` (per-function validators)
- `MD/` (docs incl. this file, README, CHANGELOG), `Dashboard_Demo/`, `env/` (requirements.txt), `installer/`
- EDIT ONLY the S: clone `S:\Images\Weihao\NISeA\NIS-E-Automation\` (scope-PC share the GUI runs from) — never the `C:\Users\weiha\SpheroidPA\` clone; both push to the same PUBLIC origin. Working branch `NISE-dispatcher` (current v1.8.8); `master` is still the pre-restructure flat layout.

## Pipeline steps (Python in `Codes/GUI/`, macros in `macro/`)
1. `spheroid_screener.py`       — detect & rank spheroids in 10X whole-well mosaic nd2
2. `cross_zoom_register.py`     — pre-capture verification + post-capture validation (v1, with sub-10X correction)
3. `cross_zoom_v2.py`           — v2 workflow: sub-10X nd2 stage coords → NCC match → offset report
4. `nis_macro_z_autofocus.mac`  — NIS-E daemon: polls af_trigger_NN.ini, autofocuses each spheroid, writes af_done_NN.ini
5. `nis_macro_capture_zcorrected.mac` — NIS-E daemon: polls spheroid_trigger.ini, runs Z-series w/ Z-intensity (Beer-Lambert) correction from a `.bin`, writes spheroid_done.ini
6. `nis_bin_to_csv.py`          — decode NIS-E Z Intensity Correction .bin → CSV
7. `csv_to_nis_bin.py`          — encode CSV → .bin (round-trip validated)

Also in `macro/`: `nis_macro_capture_zstack.mac` (plain Z-series; optional per-trigger `oc=` → `SelectOptConf`, drives Pre-PA/Validation) and the PA workflow — `pa_setup` (activation OC + centred zoom square) → `pa_points` (ND multipoint) → **manual** `step3_zstack_PA` JOB in NIS-E → Validation re-image; `pa_validate`/`pa_pick_current` remain dispatcher-routable but their GUI cards were consolidated away in v1.8.8. `nis_macro_dispatcher.mac` routes GUI commands to all step macros by NUMERIC action_id (1 autofocus, 2 zstack, 3 zcorrected, 4 pa_setup, 5 pa_points, 6 pa_pick_current, 7 pa_validate) via `cmd.ini`/`cmd_done.ini`.

## Key facts
- 10X mosaic pixel size: ~0.644 µm/px; 20X pixel size: ~0.321 µm/px; scale factor ~0.499
- Stage coordinate frames differ between mosaic session and sub-10X/20X session (encoder re-zero)
- Confirmed stage offset (SLIM050, 2026-06-03): dx ≈ −335 µm, dy ≈ +3415 µm
- Sub-10X nd2 stage metadata IS the correct reference for 20X navigation — not the mosaic screener coords
- Watershed splitting required to separate merged blobs (min_distance = max(4, int(min_diam/2/pixel_um)))
- Screener ranks 9 spheroids (205–222 µm); 20X demo file confirmed as Rank 9 (NCC=0.991)

## Demo data paths
```
demo_data/
  2026-06-03 SLIM050 + SLIM047 QC check/
    WellD05_Channel555_spIII1_Seq0000.nd2          # 10X whole-well mosaic
    spheroid 1/sph 1 1P/
      spheroid 1 10x 1P 555nm.nd2                  # sub-10X single-position
      20x spheroid 1 1P 555nm.nd2                   # 20X capture (= Rank 9)
  screener_out/
    WellD05_Channel555_spIII1_Seq0000_screen.csv    # 9 ranked spheroids
    spheroid_screen_latest.csv
  register_out/
    verified_screen.csv                             # NCC-refined coords
    validation_report.csv                           # with correction columns
    sub10x_stage_offset.csv                         # v2 offset table
```

## Dependencies
```
pip install nd2 numpy scipy scikit-image matplotlib
pip install cellpose   # optional, for --backend cellpose
```

## Coding conventions (this project)
- No comments except non-obvious WHY (hardware quirks, format reverse-engineering)
- Unicode arrows (→, ×, µ) cause cp1252 encode errors on Windows — use ASCII in print() strings
- All stage coords in µm; pixel coords as (col, row) = (x, y)
- stage_to_pixel: col = (x - cx) / pum + W/2,  row = (y - cy) / pum + H/2
- Never use SimilarityTransform on <2 points; fall back to pure translation
- Figures saved to demo_data/register_out/ with matplotlib.use("Agg") — no display

## Known issues / TODOs
- TODO (operator-requested 2026-07-15): default the `step3_zstack_PA` JOB's ND-acquisition Save-to-File
  path to the run's `work/<run>/pa/` folder. The PA job's ND Save-to-File is currently set MANUALLY in
  NIS-E and is NOT driven by `pa_trigger.ini` (whose `save_dir=work/<run>/pa` is only a reminder) -- so on
  0715 the activation output landed in a prior session's `D:\...\Brandon` folder (stale path) and
  `work/0715/pa/` stayed empty, leaving the delivered PA dose/power un-auditable after the fact. Fix: have
  PA Setup (or the Initialization setup macro below) re-point the PA job's ND Save-to-File to
  `pa_trigger.ini`'s `save_dir` before the job runs -- candidate `ND_DefineExperiment` re-points the save
  path but also redefines T/XY/Z/L, so preserve the job's loop/z config (verify the exact save-path fn
  against the docs index at `C:\Program Files\NIS-Elements\Docs\nis\eng_ar\`). Same root cause as the ND
  "Save to File" item in the Initialization-tab TODO below.
- DONE (2026-07-17, After-PA Validation card + `phase`-parameterized `_pl_prepa_capture_ocs`/`_pl_prepa_run_thread`;
  new `_pl_postpa_oc_*` checkboxes, saves to `nd2/postPA_<tag>/`): the Validation card captures all 3 channels into `prePA_<tag>`
  subfolders (prePA_890nm_mBeRFP / prePA_940nm_PAsfGFP / prePA_1050nm_depth), but Validation is run BOTH
  before AND after PA -- so the after-PA run currently overwrites (or collides with) the before-PA data in
  the same `prePA_` folders, losing the before/after distinction (this session had to be hand-sorted into
  Before_pa/After_pa dirs for the comparison figures). Fix: add a SEPARATE "After-PA Validation" card under
  PA Points with the SAME channel-selection interface, but whose captures save to phase-named folders
  (e.g. `postPA_<tag>/` or `afterPA_<tag>/`) instead of `prePA_<tag>/`. So the flow reads: run the
  before-PA Validation card (-> prePA_*), run PA Setup/Points/Job3, then run the after-PA card (-> postPA_*),
  keeping the two phases in distinct, correctly-named folders. Implementation: the capture path is shared
  (`_pl_prepa_capture_ocs` writes `nd2/prePA_<tag>` via `_point_nd2_dir`); parameterize the phase prefix
  ("prePA"/"postPA") per card (or add a before/after toggle) so the same code serves both. This also
  auto-produces the Before_pa/After_pa layout the comparison figures already expect, removing the manual
  sorting step.
- DONE (2026-07-20, v1.13.0; Init card + `_pl_init_rig` -> `nis_macro_init_rig.mac` action 8: zoom=2 / dichroic / PFS,
  a1_on-gated; per-OC re-save + gains still rig-side): add an "Initialization" tab/card in Step 4 that puts the rig into
  a known-good default state before an experiment, via a dispatcher-run setup macro. Primary job: set the
  confocal scan ZOOM to 2 for ALL OC profiles used by the pipeline (the OC-baked-zoom inconsistency keeps
  biting -- e.g. 0715 the 890nm OC "890nm_Galvo_600nm_NDD2_BT" was baked at zoom 3.3 while 940/1050 were at
  zoom 2, producing a wrong-scale 890 pass that had to be re-captured). Two possible implementations:
  (a) persistently re-save each OC's scan area at zoom 2 (loop SelectOptConf -> Confocal_SetScanArea(2,..)
  -> save-OC; need to find the OC-save macro fn), or (b) a runtime enforcement that the capture macro
  applies Confocal_SetScanArea(2,..) after every SelectOptConf regardless of the OC's baked value (more
  robust; already tracked as the zoom-normalization macro fix). The tab should ALSO surface + default the
  other recurring global variables that have caused incidents this project, each editable but pre-set to
  its safe default:
    * confocal zoom = 2 (all OCs)
    * PA activation power = 30% (MAX_PA_ACTIVATION_POWER_PCT; 80% burned A02 sph#9)
    * ND "Save to File" = OFF / re-pointed to the pipeline nd2_dir (0715 dumped 940/890/1050 dupes into a
      prior session's D:\...\Brandon folder because Save-to-File was left checked with a stale path; no
      macro fn disables it -- ND_DefineExperiment only re-points, and redefines T/XY/Z/L, so use carefully)
    * per-OC detector channel correct (940 must route 488 not 640 -- the "..._JL" vs "..._JL2" bug)
    * detector gains sane (avoid ReflectNDDPMT overload)
    * z-stack defaults (z_centre/z_half/z_step), dichroic/PFS state
  Goal: one click at experiment start guarantees zoom/power/save/channel/gain are all at known-good values
  instead of inheriting whatever a prior manual/other-user session left in NIS-E. Pairs with the ND-setup
  macro and Save-to-File findings above.
- TODO (operator-requested 2026-07-10): add an "Abort all running macros" button in the GUI (dispatcher
  panel / Step 4) to stop the currently-running macro without alt-tabbing to NIS-E and pressing Esc.
  Non-trivial because a worker macro launched via the dispatcher's `RunMacro()` BLOCKS the dispatcher
  loop, so the dispatcher can't process a new cmd.ini "abort" while a macro is running. Design: (1) the
  GUI Abort button writes an `abort.ini` flag into work_dir (and clears the GUI-side `_pl_dispatch_busy`
  lock so the GUI unblocks immediately); (2) every long-running loop macro (`nis_macro_capture_zstack`,
  `nis_macro_pa_points`, `nis_macro_z_autofocus`, `nis_macro_capture_zcorrected`) polls
  `ExistFile(work_dir/abort.ini)` at the top of each loop iteration and returns cleanly (writing
  status=aborted to cmd_done.ini + deleting the abort flag) if present; (3) for a macro NOT in a poll
  loop (mid-capture / mid-ND-Z-series), the abort can only take effect at the next iteration boundary --
  document that a true immediate hard-stop still needs NIS-E-side Esc (verify whether a
  `StopAllMacros()`/equivalent NIS-E fn exists that could be triggered, but it would have to run from
  within NIS-E, not the GUI). Keep the abort flag one-shot (delete after honoring) so it doesn't kill the
  next dispatched macro.
- TODO (Validation OC-switch overhead): the GUI Validation capture ALREADY dispatches wavelength-outer
  (`_pl_prepa_capture_ocs`, spheroid_pa_gui.py ~3081: outer loop over checked OCs, each dispatch captures
  ALL checked spheroids in one `nis_macro_capture_zstack.mac` pass) -- so with N spheroids it images all
  N at 890, then all N at 940, then all N at 1050, NOT 3 wavelengths per spheroid (single-spheroid tests
  just make the two orders look identical). The remaining waste is INSIDE the macro: `SelectOptConf(oc)`
  sits in the per-trigger loop, so it re-selects the SAME oc once per spheroid -> (N-1) redundant OC
  switches per wavelength (3*(N-1) total). Fix: track the last-selected OC and skip `SelectOptConf` when
  the trigger's oc equals the currently-active one (or hoist the select to once per pass), so each
  wavelength's OC switch happens exactly once regardless of spheroid count. Cuts the between-channel wait
  the operator flagged 2026-07-09. Verify on the rig whether SelectOptConf is already a no-op when the OC
  is unchanged (if so this is free); if not, the guard is the win.
- TODO (operator-requested 2026-07-09): during Validation mode, DEFAULT to NOT displaying each captured
  ND2 in NIS-E. Right now `nis_macro_capture_zstack.mac` leaves every captured image open after
  `ImageSaveAs`, so a multi-spheroid x multi-wavelength x z-stack Validation run opens dozens of image
  windows -- clutters NIS-E and can bog it down / eat memory (worst with z-stacks: N spheroids x 3 OCs x
  many planes). Fix: after `ImageSaveAs(nd2_path, 14, 0)` in the capture macro, close the just-saved
  document (verify the exact NIS-E fn -- candidates `CloseCurrentDocument()` / `CloseDocument()` /
  `Close()` -- against the docs index at `C:\Program Files\NIS-Elements\Docs\nis\eng_ar\`), or capture
  into a non-displayed buffer. Keep it the DEFAULT for Validation (the GUI's Captured Z-Stacks viewer is
  where the operator reviews them anyway, so no need to also show them live in NIS-E); optionally expose a
  "keep captures open in NIS-E" checkbox for debugging. Don't close the ACTIVE/live window or anything the
  operator needs -- only the just-saved capture doc.
- DONE (2026-07-17, `_pl_regen_triggers` in spheroid_pa_gui.py, wired into `_pl_pa_run_macro` /
  `_pl_prepa_capture_ocs` / `_pl_pa_run_pipeline`; also re-syncs the stale session.ini work_dir):
  EVERY Step 4 Run button (Validation captures,
  PA Points, and Run Pipeline) should first REGENERATE the `af_trigger_*.ini` files from Step 3's
  currently-checked spheroid selection, instead of relying on whatever triggers happen to already be on
  disk. Today none of the Step 4 actions have a coordinate source of their own -- they all read the
  trigger files Step 3 last wrote, which get silently invalidated at least four ways seen this session:
  (1) a plate rescan supersedes the coords (stale, ~5 mm off, 0709 #5); (2) reusing an archived trigger
  set (`#1-5/`) that predates the rescan; (3) Validation/Pre-PA consuming them via
  `nis_macro_capture_zstack.mac`'s DeleteFile-after-capture (defensively restored in v1.9.1, but the
  fragility remains); (4) simply forgetting to re-click Step 3 Trigger after changing the checkbox
  selection. Fix: factor out Step 3's trigger-write (`trigger_autofocus_all` over `self._pl_records`
  minus `self._pl_excluded`) into a helper, and call it at the top of each Step 4 Run handler
  (`_pl_pa_run_macro` for pa_points, `_pl_prepa_run_thread`/`_pl_prepa_capture_ocs` for Validation,
  `_pl_pa_run_pipeline` for the pipeline) so the triggers are always freshly derived from the live
  checked selection at click time. This makes Step 3's checkbox table the single authoritative selection
  end to end, and eliminates the whole class of stale-trigger / wrong-location bugs that dominated the
  0709 session. (Superseded the earlier pa_points-only version of this TODO -- the operator asked to
  generalize it to all Step 4 Run buttons.)
- `SimilarityTransform.estimate()` deprecation: switch to `.from_estimate()` (skimage ≥ 0.26)
- S2 near-edge detection (col≈110/1024) gives unreliable NCC patch — clip to max_half or skip
- With only 2 anchor points, similarity transform has ~72 µm residual at Rank 6 (~900 µm from anchor)
- `csv_to_nis_bin.py` footer hardcoded for exactly 19 items; generalise for N ≠ 19
- Add a GUI button to auto-launch the `Step1_Locate_via_scan` JOB (10X whole-well mosaic) instead of
  requiring a manual run in NIS-E JOBS Explorer -- unlike `step3_zstack_PA` this has no laser-safety
  gate, so it's a reasonable auto-launch target (`_Jobs_RunJobOrWizardByName("IMAGEN",
  "Step1_Locate_via_scan")` is a valid Protocol_Commands call, confirmed in the lxapp log)
- TODO (operator-requested 2026-07-09): the Step 4 "Job3 (step3_zstack_PA -- manual JOB run)" card is
  currently NOTE-ONLY (just displays the shared Well field + a reminder to run the JOB by hand in NIS-E
  JOBS Explorer). Investigate launching step3_zstack_PA directly FROM the GUI via the Job Wizard, so the
  operator does not have to switch to NIS-E and hit Run manually. Same mechanism as the Step1 auto-launch
  note above -- a Protocol_Commands call like `_Jobs_RunJobOrWizardByName("IMAGEN", "step3_zstack_PA")`,
  dispatched through the existing cmd.ini bridge (add an action_id + dispatcher case, or a dedicated
  launch macro). CRITICAL SAFETY DIFFERENCE from Step1: step3_zstack_PA FIRES THE 850 nm PA ACTIVATION
  LASER, so unlike the harmless Step1 mosaic this MUST stay behind an explicit laser-safety gate --
  require the "A1 powered ON" confirmation AND a deliberate confirm click (never auto-fire), and confirm
  the ND multipoint import happened first (see the Step-4-Run trigger-regen TODO). Verify the exact
  wizard/job name and that _Jobs_RunJobOrWizardByName can start a JOB that itself fires the laser
  (vs only opening the wizard) before wiring it. Prove on the rig with the laser interlocked first.
- DONE (2026-07-20, v1.13.0; `nis_macro_pa_done.mac` action 9 writes `pa_done.ini`, GUI "Watch for PA done" polls it):
  add a "PA done" popup when step3_zstack_PA finishes -- currently
  the JOB completes with NO notification, so the operator has to watch the Job Execution progress bar to
  know when to run the after-PA Validation. Options: (a) a final job step that runs a command/macro doing
  `WaitText(0, "PA complete for N spheroid(s) -- run after-PA Validation next")`; (b) the JOB's
  "Execute Command after Capture" hook (seen in the ND panel) pointing at a tiny popup macro; (c) have
  that end-of-job macro also write a `pa_done.ini` flag into work_dir so the GUI can detect completion
  and surface it in the (now-timestamped) log / status line -- best tied into the existing cmd.ini bridge
  so the notification shows on the GUI side too, not just NIS-E. Pairs with the Job3 GUI-launch TODO
  above (if the GUI launches the job, it can also poll for the done-flag and pop its own toast).
  (`SelectOptConf`) instead of relying on whatever OC is already active or an explicit per-trigger
  `oc=`, matching the established anchor@555/recenter@1050 convention; add its own GUI launch button
  separate from the general z-stack capture dispatch
- DONE (2026-07-07): Validation captures now auto-load into "Captured Z-Stacks" on completion
  (`_pl_zv_load_captured`, no manual Auto Load/Refresh click needed). STILL OPEN: default view is
  one-at-a-time via the combobox, not the originally-envisioned one-row-per-channel (3 rows,
  890/940/1050nm) side-by-side layout. Also still open: crop/rescale 890nm and 940nm to match
  1050nm's FOV before display -- 890/940nm capture at ~636 um FOV vs 1050nm's ~318 um (that OC's
  own saved zoom, see 2026-07-06 finding), so a like-for-like visual comparison needs the wider
  channels center-cropped to 1050nm's physical extent (known px sizes: 890nm 1.243 um/px 512px,
  940nm 0.622 um/px 1024px, 1050nm 0.622 um/px 512px)
- Confocal scan resolution/format (512/1024/512x128/1024x256/256x256, shown as preset buttons in
  NIS-E's "A1plus Scan Area" panel) is a DIFFERENT parameter from `Confocal_SetScanArea`'s Zoom/
  Angle/XOffset/YOffset/Aspect (confirmed via the local NIS-E function reference -- that function
  has no resolution/frame-size argument). `CameraFormatSet(int CameraPropMode, char
  *CameraFormatParam)` is a candidate but UNCONFIRMED whether it applies to the A1plus confocal
  scan format specifically vs a widefield camera -- needs verification on the rig (or in NIS-E's
  macro test/command-line panel) before wiring into `nis_macro_capture_zstack.mac`. This is what
  caused Pre-PA/Validation's 890/940nm-vs-1050nm scale mismatch (2026-07-06/07 finding above).
- Off-center PA box: the faded activation square lands ~18 µm lower-right of the spheroid centre — a
  CONSTANT scanner zoom-centre offset (centre/image at zoom 2.5 vs activate at zoom 8), NOT a Step 1-3
  centring error and NOT per-spheroid. Fix = one calibrated `XOffset`/`YOffset` in `pa_setup`'s
  `Confocal_SetScanArea(zoom, 0, 0, 0, 1)` (currently 0,0), measured once on a uniform fluorescent slide
  (µm offset → scan-area units). Deferred (needs rig). Not yet a README TODO row.
- PA activation power hard-capped at 30% (`MAX_PA_ACTIVATION_POWER_PCT`, v1.8.8): 80% burned A02 sph#9
  (saturated hot spot across 890/940/1050 nm) -- corrected 2026-07-09 from an earlier "50%" attribution,
  per the user (was mistaken originally; 80% is the confirmed damaging power). Soft clamp only —
  `step3_zstack_PA` (run manually in the Job Wizard) does NOT read `pa_trigger.ini`'s `power_pct`; real
  laser power is set by hand in NIS-E.
- PA depth coverage: `step3_zstack_PA` activated only ~6 planes / ~25 µm of Z (fixed power across them),
  so a ~250 µm spheroid's core is under-dosed — two orthogonal gaps: dose-vs-depth (README TODO #8) and
  the too-narrow Z RANGE.
- SUPERSEDED (2026-07-09, v1.9.0): PA Points no longer has a Count field at all. It previously conflated
  COUNT with a rank-number ceiling (`i<=count` instead of "how many to include" -- Count=1/2 silently
  found 0 when testing a single selectively-triggered spheroid whose filename kept its original rank,
  e.g. `af_trigger_09.ini`) and separately could read a stale/not-yet-regenerated trigger set if Step 3
  Trigger wasn't (re-)clicked first. Both are moot now: `nis_macro_pa_points.mac` scans the full 1..96
  rank range and builds the multipoint from EVERY currently-active `af_trigger_NN.ini`, with no cap --
  Step 3's checked-row table + "Trigger" click (which clears stale af_trigger_*/af_done_* first, see
  `trigger_autofocus_all`, spheroid_pipeline.py:336-364) is now the ONLY selection mechanism. To PA a
  different subset: check/uncheck rows in Step 3 and click Trigger again before running PA Points.
- FIXED (2026-07-09, v1.8.9): Validation's 940nm PAsfGFP capture always came back with only channel 640
  (no green signal) -- NOT a macro bug (`nis_macro_capture_zstack.mac`'s only OC-related code is a single
  `SelectOptConf(oc)` call, confirmed by full re-read) and NOT `SelectOptConf` clobbering channel state as
  a side effect. Root cause: NIS-E's OC tree has TWO distinct, separately-configured entries with
  confusingly similar names -- `940nm_Galvo_488nm_NDD2_JL2` (640 active, 488 off -- what the dispatcher
  was actually calling) vs `940nm_Galvo_488nm_NDD2_JL` (no trailing "2" -- 488 live, sane gain, confirmed
  on-rig to capture real PAsfGFP puncta). `_pl_prepa_checked_ocs()` now points at the "..._JL" (no "2")
  name. Lesson for any future OC wiring in this pipeline: OC names in this library are NOT self-describing
  reliably -- verify the ACTUAL live channel/gain state in NIS-E's A1plus Pad for the EXACT name being
  passed to `SelectOptConf`, don't assume similarly-named configs are interchangeable.
