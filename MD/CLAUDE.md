# Claude Code — SpheroidPA Project

## What this project is
Automated spheroid screening and cross-zoom registration pipeline for NIS-Elements microscopy.
Sits between NIS-E acquisition and 20X capture: detects, ranks, and validates spheroid locations.

## Pipeline steps
1. `spheroid_screener.py`       — detect & rank spheroids in 10X whole-well mosaic nd2
2. `cross_zoom_register.py`     — pre-capture verification + post-capture validation (v1, with sub-10X correction)
3. `cross_zoom_v2.py`           — v2 workflow: sub-10X nd2 stage coords → NCC match → offset report
4. `nis_macro_z_autofocus.mac`  — NIS-E daemon: polls af_trigger_NN.ini, autofocuses each spheroid, writes af_done_NN.ini
5. `nis_macro_auto_capture.mac` — NIS-E daemon: polls spheroid_trigger.ini, runs Z-series w/ Z-intensity correction, writes spheroid_done.ini
6. `nis_bin_to_csv.py`          — decode NIS-E Z Intensity Correction .bin → CSV
7. `csv_to_nis_bin.py`          — encode CSV → .bin (round-trip validated)

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
- `SimilarityTransform.estimate()` deprecation: switch to `.from_estimate()` (skimage ≥ 0.26)
- S2 near-edge detection (col≈110/1024) gives unreliable NCC patch — clip to max_half or skip
- With only 2 anchor points, similarity transform has ~72 µm residual at Rank 6 (~900 µm from anchor)
- `csv_to_nis_bin.py` footer hardcoded for exactly 19 items; generalise for N ≠ 19
- Add a GUI button to auto-launch the `Step1_Locate_via_scan` JOB (10X whole-well mosaic) instead of
  requiring a manual run in NIS-E JOBS Explorer -- unlike `step3_zstack_PA` this has no laser-safety
  gate, so it's a reasonable auto-launch target (`_Jobs_RunJobOrWizardByName("IMAGEN",
  "Step1_Locate_via_scan")` is a valid Protocol_Commands call, confirmed in the lxapp log)
- Refactor `nis_macro_capture_zstack.mac` into a dedicated recenter macro that auto-selects 1050nm
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
