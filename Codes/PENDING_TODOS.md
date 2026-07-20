# SpheroidPA — Pending TODO Items

_Generated 2026-07-15. Un-implemented items only, pulled from `MD/CLAUDE.md`._
_Excludes items already marked DONE / FIXED / SUPERSEDED (940nm OC fix, CRLF macro bug, PA power cap, PA Points Count removal, Validation auto-load, per-channel crop/tint, GUI log timestamps + PA status)._

Priority key: **P0** = blocks correct/safe runs · **P1** = high-value quality-of-life · **P2** = nice-to-have / cleanup · **DEFER** = needs rig calibration or unconfirmed API.

---

## P0 — Correctness / data integrity

### 1. Regenerate triggers on every Step 4 Run — DONE 2026-07-17  *(2026-07-09, HIGH)*
**Implemented:** new `_pl_regen_triggers()` helper (delegates to `trigger_autofocus_all` → clears stale af_trigger/af_done, writes current/re-centered coords from `self._pl_records` minus `self._pl_excluded`, and re-syncs `session.ini` work_dir/nd2_dir) is called at the top of `_pl_pa_run_macro`, `_pl_prepa_capture_ocs`, and `_pl_pa_run_pipeline`. Step 3's checkbox table is now the single authoritative selection end-to-end; also fixes the stale-`session.ini` work_dir bug (0706-vs-0715).
Original problem: Step 4 Runs relied on whatever triggers were on disk, which went stale ≥4 ways (plate rescan, reused archived set, capture-macro deletion, or forgetting to re-click Step 3 Trigger).

### 2. Capture macro must own its ND save target (Save-to-File stale-path)  *(2026-07-15)*
`ND_RunZSeriesExp()` honors the ND dialog's "Save to File" path/filename, which no macro controls — on 0715 it dumped 890/940/1050 duplicates into a prior session's `D:\...\Brandon` folder. No macro fn *disables* Save-to-File; `ND_DefineExperiment(...)` can **re-point** it (but also redefines T/XY/Z/L — preserve those flags).
**Fix:** before each `ND_RunZSeriesExp()` in `nis_macro_capture_zstack.mac` (+ `nis_macro_capture_zcorrected.mac`, `nis_macro_pa_validate.mac`), re-point the ND save to the pipeline `nd2_path`. Interim workaround in use: **uncheck "Save to File" manually** before pipeline runs.

### 3. PA job Save-to-File path → `work/<run>/pa/` — DONE 2026-07-20 (needs rig verify)  *(2026-07-15)*
**Implemented:** `nis_macro_pa_setup.mac` now reads `pa_trigger.ini`'s `save_dir` and re-points the ND save target via `ND_DefineExperiment(TRUE,TRUE,TRUE,TRUE,TRUE, save_dir, "PA", 0, FALSE,FALSE,FALSE)` before the JOB runs. **Verify once on the rig** that Points/ZStack/Time dimensions are unchanged and `pa/` fills.
`step3_zstack_PA`'s ND Save-to-File is set MANUALLY and ignores `pa_trigger.ini`'s `save_dir` (reminder-only) — so on 0715 the activation output went to the stale Brandon folder and `work/0715/pa/` stayed empty, leaving the delivered PA dose un-auditable.
**Fix:** have PA Setup (or the Initialization macro, item 6) re-point the PA job's ND Save-to-File to `pa_trigger.ini`'s `save_dir` before the job runs. Same root cause as item 2.

### 4. Enforce consistent confocal zoom (OC-baked-zoom)  *(2026-07-15)*
The capture macro never sets zoom — it inherits each OC's baked scan-area zoom. On 0715 the 890nm OC (`890nm_Galvo_600nm_NDD2_BT`) was baked at **zoom 3.3** while 940/1050 were at zoom 2, producing a wrong-scale 890 pass that had to be re-captured.
**Fix (robust):** add `Confocal_SetScanArea(2.0, 0,0,0, 1.0)` after `SelectOptConf(oc)` in `nis_macro_capture_zstack.mac`, with a per-OC `zoom` trigger key (GUI emits it). Alt (rig-side): re-save each OC's scan area at zoom 2.

---

## P1 — GUI / Step 4 features

### 5. "After-PA Validation" card (postPA folders) — DONE 2026-07-17  *(2026-07-15)*
**Implemented:** `_pl_prepa_checked_ocs` / `_pl_prepa_capture_ocs` / `_pl_prepa_run_thread` / `_pl_prepa_run` now take a `phase` arg (`"prePA"`/`"postPA"`) — the folder prefix, OC checkboxes, locate var, and status labels all key off it. A new **After-PA Validation** card (own `_pl_postpa_oc_890/940/1050` + locate checkboxes) sits under PA Points and runs `_pl_prepa_run_thread("postPA")`, saving each pass to `nd2/postPA_<tag>/`. Positions/Z come from Step 3's checked rows (via item #1's regen). Auto-produces the Before_pa/After_pa split the comparison figures expect.
Original problem: before AND after both wrote to `prePA_<tag>/`, so this session had to be hand-sorted.

### 6. "Initialization" tab — set rig to known-good defaults  *(2026-07-15)*
One click at experiment start that sets a dispatcher-run setup macro to normalize the rig, so it never inherits a prior session's state. Surface + default each recurring global (editable):
- confocal **zoom = 2** for all OC profiles (persistently re-save each OC, or runtime-enforce via item 4)
- **PA activation power = 30%** (80% burned A02 sph#9)
- **ND "Save to File" = OFF / re-pointed** (items 2, 3)
- **per-OC detector channel correct** (940 → 488, the `_JL` vs `_JL2` bug)
- **detector gains sane** (avoid ReflectNDDPMT overload)
- z-stack defaults (z_centre/z_half/z_step), dichroic / PFS state

### 7. "Abort all running macros" button — DONE 2026-07-20  *(2026-07-10)*
**Implemented:** "Abort All" button in the NIS-E Macro Dispatcher panel → `_pl_abort_all()` writes a one-shot `abort.ini` into work_dir and clears `_pl_dispatch_busy`. `nis_macro_capture_zstack.mac` polls `ExistFile(work_dir/abort.ini)` once per rank *and* on the idle path, exits cleanly, and writes `status=aborted` + a message into `cmd_done.ini` so the GUI log shows "aborted" not "ok". The GUI-side Validation OC loop and Run Pipeline loop poll it between passes/steps; `_pl_abort_clear()` runs at the start of each new run so a stale flag can't kill it. Helpers: `_pl_abort_path` / `_pl_abort_requested` / `_pl_abort_clear` / `_pl_abort_all`.
**Known limit (documented in the confirm dialog):** a macro already inside a Z-series/Capture stops only at the next loop boundary — immediate hard stop still needs Esc in NIS-E.
**Still open (optional):** add the same poll to `nis_macro_z_autofocus.mac` and `nis_macro_capture_zcorrected.mac` (both are long imaging loops; only `capture_zstack` — the one Validation/PA uses — is wired so far).

### 8. Hide captured ND2s in NIS-E during Validation — DONE 2026-07-20  *(2026-07-09)*
**Implemented:** Step 4 checkbox "Keep captured images open in NIS-E" (default OFF) -> `close_after=1` in each trigger -> `CloseCurrentDocument(2)` after `ImageSaveAs` in `nis_macro_capture_zstack.mac`.
`nis_macro_capture_zstack.mac` leaves every capture open after `ImageSaveAs`, so a multi-spheroid × multi-wavelength × z-stack run opens dozens of windows (clutter + memory).
**Fix:** close the just-saved doc after `ImageSaveAs` (verify fn: `CloseCurrentDocument`/`CloseDocument`/`Close` in the docs index) or capture to a non-displayed buffer. Default ON for Validation; optional "keep open" checkbox for debug. Never close the active/live window.

---

## P1 — Job automation

### 9. Launch `step3_zstack_PA` (Job 3) from the GUI  *(2026-07-09)*
The Job3 card is note-only. Launch the JOB via `_Jobs_RunJobOrWizardByName("IMAGEN", "step3_zstack_PA")` through the cmd.ini bridge.
**CRITICAL SAFETY:** step3_zstack_PA fires the 850 nm PA laser — must stay behind an explicit laser-safety gate ("A1 powered ON" + deliberate confirm click, never auto-fire) and confirm the ND multipoint import happened first. Verify the wizard call can *start* a laser-firing job (vs only open the wizard); prove with the laser interlocked first.

### 10. "PA done" completion popup / flag  *(2026-07-09)*
`step3_zstack_PA` finishes with no notification — operator watches the progress bar to know when to run after-PA Validation.
**Fix:** final job step / "Execute Command after Capture" hook runs a macro doing `WaitText("PA complete...")` and writes a `pa_done.ini` flag into work_dir so the GUI surfaces it in the (timestamped) log. Pairs with item 9.

### 11. Auto-launch `Step1_Locate_via_scan` (10X mosaic) from the GUI
Add a GUI button; `_Jobs_RunJobOrWizardByName("IMAGEN", "Step1_Locate_via_scan")` is a valid call (confirmed in lxapp log). No laser-safety gate needed (imaging-only), unlike item 9.

---

## P2 — Cleanup / minor

### 12. Validation OC-switch overhead  *(2026-07-09)*
Dispatch is already wavelength-outer (good). Remaining waste: `SelectOptConf(oc)` is re-called once per spheroid inside a pass even though the OC is identical → (N-1) redundant OC switches per wavelength.
**Fix:** track last-selected OC and skip `SelectOptConf` when unchanged (or hoist to once per pass). Verify whether SelectOptConf is already a no-op on unchanged OC (then it's free).

### 13. Refactor `nis_macro_capture_zstack.mac` → dedicated recenter macro
Auto-select 1050nm via `SelectOptConf` instead of relying on the active OC / per-trigger `oc=`, matching the anchor@555/recenter@1050 convention; add its own GUI launch button separate from the general z-stack dispatch.

### 14. Captured Z-Stacks viewer — one-row-per-channel layout
Auto-load DONE; the crop/tint DONE (v1.9). Still open: the default view is one-at-a-time via combobox, not the envisioned 3-row (890/940/1050) side-by-side layout.

### 15. `csv_to_nis_bin.py` footer hardcoded for exactly 19 items — generalize for N ≠ 19.

### 16. `SimilarityTransform.estimate()` deprecation — switch to `.from_estimate()` (skimage ≥ 0.26).

### 17. S2 near-edge detection (col ≈ 110/1024) gives unreliable NCC patch — clip to max_half or skip.

### 18. 2-anchor similarity transform has ~72 µm residual at Rank 6 (~900 µm from anchor) — known limitation; add a 3rd anchor for precise fit.

---

## DEFER — needs rig calibration or unconfirmed API

### 19. Confocal scan resolution/format (512/1024/512×128/…)
A DIFFERENT parameter from `Confocal_SetScanArea`'s zoom. `CameraFormatSet(int, char*)` is a candidate but UNCONFIRMED for the A1plus confocal — verify on the rig / NIS-E macro test panel before wiring.

### 20. Off-center PA box (~18 µm lower-right of centre)
CONSTANT scanner zoom-centre offset (image@zoom 2.5 vs activate@zoom 8), not per-spheroid. Fix = one calibrated `XOffset`/`YOffset` in `pa_setup`'s `Confocal_SetScanArea`, measured once on a uniform fluorescent slide. Needs rig.

### 21. PA depth coverage
`step3_zstack_PA` activated ~6 planes / ~25 µm at fixed power → ~250 µm spheroid core under-dosed. Two gaps: dose-vs-depth, and too-narrow Z range.

---

_Source of truth: `MD/CLAUDE.md` "Known issues / TODOs". Update there and regenerate this file as items land._
