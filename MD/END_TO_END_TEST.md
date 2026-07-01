# SpheroidPA — End-to-End Test Guide

How to exercise the full pipeline: **SpheroidPA GUI (`spheroid_pa_gui.py`) ⇄ NIS-E macros**,
from a 10X mosaic to Beer-Lambert-corrected 20X Z-stacks.

The Python side and NIS-E talk through **file-flag INI triggers** in a shared work
directory. There is no socket/COM link — Python writes `*_trigger.ini`, a NIS-E macro
daemon polls for it, acts, and writes `*_done.ini` back.

```
GUI Step 3  ──af_trigger_NN.ini──►  nis_macro_z_autofocus.mac   ──af_done_NN.ini──►  GUI
GUI Step 4  ─spheroid_trigger.ini─► nis_macro_auto_capture.mac ─spheroid_done.ini─► GUI
```

---

## 0. One-time setup

### 0.1 Python environment (control PC)
```
pip install nd2 numpy scipy scikit-image matplotlib
```
Launch the GUI from the repo folder so it can import its sibling modules:
```
cd S:\Images\Weihao\NISeA\NIS-E-Automation
python spheroid_pa_gui.py
```

### 0.2 Shared work directory
Create the directory the macros expect. The macro `WORK_DIR` `#define`s are:

| Macro | `WORK_DIR` | Reads | Writes |
|---|---|---|---|
| `nis_macro_z_autofocus.mac` | `C:/SpheroidPA/work/autofocus` | `af_trigger_NN.ini` | `af_done_NN.ini`, `<id>_af.nd2` |
| `nis_macro_auto_capture.mac` | `C:/SpheroidPA/work` | `spheroid_trigger.ini` | `spheroid_done.ini`, `<id>.nd2` |

```
mkdir C:\SpheroidPA\work\autofocus
mkdir C:\SpheroidPA\work\bins
mkdir C:\SpheroidPA\work\nd2
```
**The GUI paths MUST match these** (this is the #1 cause of a silent no-op):
- Step 1 **Output directory** → `C:\SpheroidPA\work`  (so Step 3 writes to `...\work\autofocus`)
- Step 4 **Trigger dir** → `C:\SpheroidPA\work`
- Step 4 **ND2 output dir** → `C:\SpheroidPA\work\nd2`

If you want different paths, edit the macro `#define WORK_DIR` lines to match — and keep them CRLF (see 0.4).

### 0.3 NIS-E hardware / acquisition config (do BEFORE running any macro)
- Switch the turret to the **20X objective**.
- **XY + Z stages** connected, initialized, in automatic mode.
- **Z drive must NOT be in Escape mode** — otherwise the Z-series warns *"Z Drive is in
  Escape mode, Z movement is not available"* and only a single plane is captured. Press
  Escape/Return to send the focus drive back to its working position first.
- Camera **not in Live** (capture/run functions freeze if Live is on).
- Open **Acquire ▸ ND Acquisition** and configure channels.
- For corrected capture (Step 4): enable **Z Intensity Correction** in the ND Acquisition
  dialog, and select the **A1/confocal detector** — the Z-intensity-*corrected* run
  (`ND_RunZSeriesExpWithZIntensityCorrection`) and `ND_ZIntensityControlIsDataReady` are
  A1/multilaser-only. On a non-A1 camera the panel shows "Incompatible Z Correction and
  Camera"; set `#define USE_ZCORR_RUN 0` in `nis_macro_auto_capture.mac` for a plain
  (uncorrected) Z-series.
- Autofocus (Step 3) assumes brightfield (criterion 0). For fluorescence, set
  `#define AF_CRITERION 2` in `nis_macro_z_autofocus.mac`.

### 0.4 Install the macros in NIS-E
Copy `nis_macro_z_autofocus.mac` and `nis_macro_auto_capture.mac` to the NIS-E PC.
They are already **CRLF-encoded** — keep it that way. (A `.mac` saved with Unix `LF`
line endings makes every `//` comment swallow the rest of the file and the macro
silently does nothing. Verify with `file nis_macro_auto_capture.mac` → "CRLF line terminators".)
Open each in **Macro ▸ Macro Editor**; you run them via **Macro ▸ Run** (entry point `main()`).

---

## 1. Phase 1 — Hardware-free rehearsal (no microscope)

Prove the bridge and the Python pipeline before tying up the scope.

### 1.1 Bridge round-trip (the key check)
```
python verify_trigger_bridge.py
```
This drives the **real** `write_trigger` / `trigger_autofocus_all` writers, then reads
the files back through the **Win32 profile API** — the exact `GetPrivateProfileString`
that NIS-E's `Int_GetKeyString` / `Int_GetKeyValue` wrap. It also drives the reverse
(done-file) direction. Expect:
```
RESULT: PASS -- every key the NIS-E macros read round-trips correctly.
```
A PASS means a live macro will read identical values. If it FAILs, the bridge is broken
and there is no point going to the scope.

### 1.2 GUI manual half-loop (optional, tests the GUI's polling UI)
You can confirm the GUI's trigger/await UI without NIS-E by playing the macro's part by hand:
1. Run Steps 1–4 far enough to click **Start Capture Queue** (or hand-make a `BIN_READY` record).
2. The GUI writes `C:\SpheroidPA\work\spheroid_trigger.ini` and the row goes **Imaging…**.
3. Simulate the macro: delete the trigger and write a done file so the GUI unblocks:
   ```
   python -c "import spheroid_pipeline as p, pathlib as L; w=L.Path(r'C:\SpheroidPA\work'); \
   (w/'spheroid_trigger.ini').unlink(missing_ok=True); \
   open(w/'spheroid_done.ini','w',newline='').write('[spheroid]\r\nstatus=ok\r\nnd2_path=C:/tmp/fake.nd2\r\n')"
   ```
4. The row should flip to **IMAGED**. This validates `wait_for_done()` polling and the table refresh.

---

## 2. Phase 2 — Full end-to-end on NIS-E

Open the GUI on the **PA Workflow** tab. The left sidebar has the four steps; the right
panel shows live **Spheroid State**; the **Log** tab has the full trace.

### Step 1 — Screen Mosaic
1. **10X mosaic ND2** → browse to the whole-well 10X scan.
2. **Well ID** → e.g. `D05`.  **Output directory** → `C:\SpheroidPA\work`.
3. **Top N spheroids** → e.g. `5` (blank = all).
4. Click **Run Screener**.
- **Success:** "N spheroid(s) detected and ranked", the state table fills with ranks, and
  the Live Dashboard shows the mosaic with detections.

### Step 2 — Anchor Offset  (mosaic frame → live stage frame)
The 10X mosaic and the live 20X session use different stage origins; this measures the offset.
1. In NIS-E, drive the stage to a recognizable spheroid and capture a **sub-10X single-position ND2**.
2. In the GUI **Anchor 1**: type its **Rank** (or leave blank to auto-match) and **Browse ND2** to that file.
3. (Optional) repeat for **Anchor 2** — two anchors reduce residual.
4. Click **Verify Anchors + Apply Offset**.
- **Success:** each anchor shows `matched rank R  NCC=0.9x  dx=… dy=…`, and
  "Offset applied: dx=… dy=…". Records now carry `verified_x/y_um` (what the macros move to).
- NCC below ~0.4 = bad match; re-capture the anchor or fix the typed rank.

### Step 3 — Autofocus + Registration
1. Set **Z half-range (um)** (e.g. `18`) and **Z step (um)** (e.g. `2`).
2. In **NIS-E**: **Macro ▸ Run** `nis_macro_z_autofocus.mac`. Status bar shows
   *"AF daemon: idle, waiting for triggers"* (it polls for 300 s of idle, then exits).
3. In the GUI: click **Trigger NIS-E Autofocus Captures**. It writes `af_trigger_NN.ini`
   for every spheroid into `...\work\autofocus`.
4. Watch NIS-E move to each spheroid, two-pass autofocus, snap a plane, and write
   `af_done_NN.ini`.
5. In the GUI: click **Refresh Status** until it reads `M/N autofocus captures complete`.
6. Click **Apply Global Registration** to refine all coordinates from the captured ND2s.
- **Success:** table `Z (um)` column fills (status → `Z_KNOWN`); registration reports
  `N matches, mean residual=… um`.
- **Order matters:** start the macro *before* clicking Trigger, or do it within the 300 s window.

### Step 4 — Generate Bins + Capture
1. Set **Trigger dir** = `C:\SpheroidPA\work`, **ND2 output dir** = `C:\SpheroidPA\work\nd2`.
2. Set the **Reference .bin** to a rig-exported Z-correction `.bin` (Z Intensity Correction
   panel ▸ `Save…`). **Required** — it supplies the detector/camera identity
   (`hwUnit_Name`, `DetectorType`, `ChannelBits`) and the `HV2`/flags fields. Without it,
   generated bins are rejected as *"Incompatible Z Correction and Camera"* and the panel
   stays empty. Generated bins inherit that identity and substitute only Z + the LP2 ramp.
3. Set **Laser channel field** (`CH2LaserPower` = LP2 on the A1), **P0 (%)**, **L (um)** for
   the Beer-Lambert ramp.
4. Click **Generate All Bins** → per-spheroid `.bin` profiles (status → `BIN_READY`). The item
   count follows the reference bin (e.g. 19), so `z_half`/`z_step` no longer need to hit 19 exactly.
5. In **NIS-E**: **Macro ▸ Run** `nis_macro_auto_capture.mac`. Status bar shows
   *"…polling for spheroid_trigger.ini"* (600 s idle timeout). For corrected capture, be on the
   A1 detector; otherwise set `#define USE_ZCORR_RUN 0` (plain run).
6. In the GUI: click **Start Capture Queue**. It triggers one spheroid at a time, each row
   going **Imaging…**, and blocks up to 10 min per spheroid for the macro's `done`.
7. Per spheroid the macro: loads the `.bin` correction, programs an absolute bottom-top
   Z-series around `z_centre`, moves XY, runs the Z-series (corrected on A1, else plain), saves the ND2.
- **Success:** every row reaches **IMAGED**, "Capture queue done: N/N imaged", and
  `pipeline_state.csv` is written. ND2s land in `...\work\nd2`.

---

## 3. Success criteria at a glance

| Step | GUI says | NIS-E does | Files appear |
|---|---|---|---|
| 1 | "N detected and ranked" | — | `spheroid_screen_latest.csv` |
| 2 | "Offset applied dx/dy" | (manual sub-10X capture) | `verified_screen.csv` |
| 3 | "M/N autofocus complete" | moves + autofocus + snaps | `autofocus/af_done_NN.ini`, `<id>_af.nd2` |
| 4 | "N/N imaged" | Z-series w/ Z-correction | `nd2/<id>.nd2`, `pipeline_state.csv` |

## 4. Troubleshooting (silent-failure modes)

- **Macro runs but nothing happens / no triggers consumed** → GUI path ≠ macro `WORK_DIR`.
  Confirm Step 1 Output dir and Step 4 Trigger dir both equal `C:\SpheroidPA\work`.
- **Macro "finished" instantly, did nothing** → the `.mac` was saved with LF endings; the
  first `//` comment ate the file. Re-save as CRLF.
- **Macro moves to (0,0) or reads blanks** → trigger INI unreadable by `GetPrivateProfileString`
  (wrong line endings / missing `[spheroid]` header). Run `verify_trigger_bridge.py`; it must PASS.
- **GUI capture queue times out (no `done.ini` in 10 min)** → macro not running, or it hit its
  idle timeout and exited. Re-run the macro, then click the GUI button.
- **`zcorr_not_ready` in done status** → `ND_ZIntensityControlIsDataReady` is A1/multilaser-only
  and returns -9 on other rigs. `CHECK_ZREADY` is already `0` so the daemon won't gate on it.
- **Generated bin shows "Incompatible Z Correction and Camera" / panel stays empty** → the bin
  lacks the rig detector identity. Set the **Reference .bin** (Step 4) to a rig `Save…` export;
  generated bins then inherit `NikonA1Grabber`/`DetectorType`/`ChannelBits` and load.
- **"Z Drive is in Escape mode, Z movement is not available"** → the focus drive is parked in
  Escape; press Escape/Return to its working position, then re-run. (Not a macro error.)
- **"Cannot Evaluate the Expression" on `ND_RunZSeriesExpWithZIntensityCorrection`** → corrected
  run needs the A1 detector. Use the A1 camera, or set `#define USE_ZCORR_RUN 0` for a plain run.
- Camera in **Live** mode → capture/run freezes. Take it out of Live.
