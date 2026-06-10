# NIS-E macro self-test ladder

Isolated probes to confirm **NIS-E actually responds to every function call** used by
`nis_macro_z_autofocus.mac` and `nis_macro_auto_capture.mac` — run these on the real
instrument **before** the first real experiment.

## Why a ladder
NIS-E **compiles the whole macro before running it**. One unknown symbol anywhere makes
the *entire* file fail silently (no output at all). By isolating each function group in
its own small file, a silent failure pinpoints exactly which call this rig rejects.
Every file is instrumented with `WaitText` checkpoints ("about to call X" -> call ->
"X returned = …") so you can see precisely how far it got.

**Reading the result of any file:**
- No dialog appears at all -> a function in that file is unknown to this build (compile fail).
- "ENTER" shows but a later checkpoint never does -> that call hung or aborted at runtime.
- All checkpoints show with sensible values -> NIS-E responds to those calls. PASS.

## Run order (safe -> involved)

| File | Validates | Covers production macro | Hardware action |
|---|---|---|---|
| `00_io_inifile.mac` | Int_SetKeyString/GetKeyString/GetKeyValue, ExistFile, DeleteFile | **both** (trigger handshake) | none |
| `01_stage_read.mac` | StgGetPos, StgGetAbsPosZ | both | read-only |
| `02_stage_move.mac` | StgMoveXY, StgXY_GetSettleTime, Wait | both | moves XY to current pos (no net motion) |
| `03_autofocus.mac` | StgFocusSetCriterion, StgFocusInRangeTwoPasses | z_autofocus | sweeps Z +/-50 um |
| `04_capture_save.mac` | Capture, ImageSaveAs | z_autofocus | acquires 1 frame |
| `05_zintensity.mac` | EnableZIntensityControl, ND_ZIntensityControlLoad, ND_ZIntensityControlIsDataReady | auto_capture | loads a .bin profile |
| `06_zseries.mac` | ND_SetZSeriesExp, ND_RunZSeriesExpWithZIntensityCorrection, ImageSaveAs | auto_capture | acquires a small Z-stack |

`00`–`04` fully cover `nis_macro_z_autofocus.mac`. `00`,`02`,`05`,`06` (+`04`) cover
`nis_macro_auto_capture.mac`. The highest-risk call is in **`05`**: the docs loaders are
`.xml`/`.nd2`, not `.bin`, so `05` is where you confirm whether `ND_ZIntensityControlLoad`
accepts the pipeline's `.bin`. If it rejects it, either make `csv_to_nis_bin.py` emit the
accepted format or switch the production macro to `ND_ReuseZIntensityCorr` against a
reference `.nd2`.

## How to run
1. Copy this folder to the NIS-E PC (keep the files **CRLF-encoded** — an LF-only `.mac`
   makes the first `//` comment swallow the file and it does nothing).
2. NIS-E prerequisites: 20X objective, XY+Z stages initialized/automatic, camera **out of
   Live**. For `05`/`06`, enable Z Intensity Correction in ND Acquisition and set a real
   `BIN_PATH` in `05_zintensity.mac`.
3. Open a file in **Macro ▸ Macro Editor**, then **Macro ▸ Run** (entry point `main()`).
   Click **Continue** through each checkpoint and note where (if anywhere) it stops.
4. Work down the table. The first file that shows nothing (or stops mid-way) is the exact
   call to fix before running the GUI pipeline end-to-end.

Edit the `#define` paths at the top of each file if your `WORK_DIR` differs from
`C:/SpheroidPA/work`.
