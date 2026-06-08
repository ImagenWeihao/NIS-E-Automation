# SpheroidPA — Spheroid Screening & Cross-Zoom Pipeline

Automated pipeline for detecting, ranking, and imaging spheroids across 10X and 20X zoom levels in NIS-Elements. Covers the full workflow from mosaic screening to per-spheroid Beer-Lambert Z-intensity correction bin generation and NIS-E file-flag triggered capture.

---

## Pipeline Overview

```
10X mosaic nd2
    └─ spheroid_screener.py       → ranked CSV + SpheroidRecord list
         └─ cross_zoom_v2.py      → NCC offset (sub-10X nd2 → mosaic frame)
              └─ 20X autofocus nd2 → Z-centre per spheroid
                   └─ csv_to_nis_bin.py → per-spheroid .bin
                        └─ nis_macro_auto_capture.mac → NIS-E file-flag capture
```

---

## File Index

| File | Role |
|---|---|
| `spheroid_screener.py` | Detect & rank spheroids in 10X mosaic nd2 (watershed, size filter) |
| `cross_zoom_v2.py` | Sub-10X NCC matching → stage offset between mosaic and sub-10X frames |
| `spheroid_pipeline.py` | **Main state machine**: orchestrates all steps, holds `SpheroidRecord`, `PipelineDashboard` |
| `spheroid_pa_gui.py` | Tkinter GUI with 4-step Cross-Zoom Pipeline tab + existing tools |
| `csv_to_nis_bin.py` | Encode Beer-Lambert Z-intensity ramp → NIS-E `.bin` format |
| `nis_bin_to_csv.py` | Decode `.bin` → CSV (round-trip validation) |
| `nis_macro_auto_capture.mac` | NIS-E macro: polls file-flag, loads `.bin`, captures Z-stack per spheroid |
| `dry_run_pipeline.py` | Offline test: Steps 1-4 without NIS-E (uses local demo data) |
| `dry_run_dashboard.py` | Saves dashboard PNGs across all pipeline stages for visual verification |

---

## Key Constants & Parameters

| Parameter | Value | Notes |
|---|---|---|
| 10X pixel size | 0.644 µm/px | mosaic |
| 20X pixel size | 0.321 µm/px | scale factor ~0.499 |
| Stage offset (SLIM050) | dx ≈ −335 µm, dy ≈ +3415 µm | mosaic → sub-10X frame |
| NCC accept threshold | 0.40 | `NCC_MIN_ACCEPT` in pipeline |
| Anchor strategy | best-NCC-only | single highest match per sub-10X nd2 |
| Z half-range default | 18 µm | |
| Z step default | 2 µm | gives exactly 19 items (ROOT_TAIL_U64S constraint) |
| Beer-Lambert L | 165 µm | scattering length, user-adjustable |
| Trigger filename | `spheroid_trigger.ini` | INI key=value, polled by macro |
| Done filename | `spheroid_done.ini` | written by macro when capture completes |

---

## SpheroidRecord Status States

```
DETECTED → VERIFIED → Z_KNOWN → BIN_READY → QUEUED → IMAGING → IMAGED
                                                              ↘ FAILED
```

Dashboard circles fill from 25% alpha (DETECTED) to 100% (IMAGED/FAILED).

---

## Operator Workflow (NIS-E Hardware)

### One-time setup
1. Copy `nis_macro_auto_capture.mac` into the NIS-E macros directory.
2. Verify `ZIntCorrect_LoadFile()` is callable for your NIS-E version (test in macro editor).
3. Set trigger directory to the same path in both the GUI and the macro.

### Per-well run

**Step 1 — 10X mosaic (NIS-E then Python)**
- Acquire whole-well mosaic at 10X, save as `.nd2`.
- GUI Step 1: point at nd2, run screener → ranked list appears.

**Step 2 — Sub-10X anchor capture (NIS-E then Python)**
- Navigate to 1–2 anchor spheroids at 10X in a single-position ND acquisition, save `.nd2`.
- GUI Step 2: point at sub-10X nd2, click Estimate Offset.
- Verify dx/dy are in the expected range (~−335, ~+3415 µm for SLIM050).

**Step 3 — 20X Z autofocus (NIS-E then Python, per spheroid)**
- Switch to 20X; navigate to each spheroid's verified XY.
- Run NIS-E autofocus to find equatorial plane; capture a single-position nd2 to record Z metadata.
- GUI Step 3: load each spheroid's 20X nd2, click Record Z.

**Step 4 — Automated capture (macro + Python)**
- In NIS-E: open and run `nis_macro_auto_capture.mac` (it enters a polling loop).
- GUI Step 4: click Generate Bins + Start Queue.
- Python writes trigger → macro captures → macro writes done → Python advances to next.
- Dashboard fills green as each spheroid completes.

---

## Dashboard

`PipelineDashboard` (in `spheroid_pipeline.py`) renders a 2-panel matplotlib figure:
- **Top**: 10X mosaic background with circles colour-coded by status and zig-zag capture path.
- **Bottom**: live status table (rank, ID, mosaic XY, verified XY, Z-centre, diameter, score, bin file).

Circles always use **mosaic XY** — verified XY is a different stage frame (sub-10X encoder zero) and is shown in the table only.

---

## Dry Run Results (2026-06-08, SLIM050 demo data)

- Mosaic: `WellD05_Channel555_spIII1_Seq0000.nd2`
- Sub-10X: `spheroid 1 10x 1P 555nm.nd2`
- 20X: `20x spheroid 1 1P 555nm.nd2`

| Result | Value |
|---|---|
| Spheroids detected | 9 |
| Diameter range | 205–222 µm |
| NCC best match (rank 9) | 0.991 |
| Offset dx | −334.6 µm |
| Offset dy | +3416.2 µm |
| Z-centre (rank 9) | 7686.00 µm |
| Bin size | 27 524 bytes |
| Bin items | 19 |
| LP ramp | 15.0% → 18.7% (Beer-Lambert, L=165 µm) |

Dashboard output: `demo_data/pipeline_dryrun/dash_*.png` (15 files).

---

## Known Issues & TODO

| # | Issue | Status |
|---|---|---|
| 1 | `csv_to_nis_bin.py` ROOT_TAIL_U64S hardcoded for exactly 19 items — workaround: use z_half=18, z_step=2 | Open |
| 2 | `ZIntCorrect_LoadFile()` function name unverified against actual NIS-E installation | Open |
| 3 | Near-edge NCC false positives (col < ~100/1024 px) give unreliable offset — pipeline now takes only best-NCC match, but edge clip would be cleaner | Open |
| 4 | `SimilarityTransform.estimate()` deprecation in skimage >= 0.26 — switch to `.from_estimate()` | Open |
| 5 | With only 2 anchor points, similarity transform has ~72 µm residual at ~900 µm from anchor | Known limitation |
| 6 | GUI `_btn()` defined twice (lines ~720 and ~1646) — harmless (Python uses last definition), but should be merged | Open |
| 7 | `nis_macro_auto_capture.mac` not yet tested against live NIS-E hardware | Open |

---

## Dependencies

```
pip install nd2 numpy scipy scikit-image matplotlib
pip install cellpose   # optional, for --backend cellpose
```

NIS-Elements: macro engine required; Z Intensity Correction module required for `ZIntCorrect_LoadFile()`.
