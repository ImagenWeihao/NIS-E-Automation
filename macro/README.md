# SpheroidPA NIS-Elements macros

The primary NIS-Elements AR macros that drive the microscope for the SpheroidPA
pipeline. The Python GUI (`../Codes/GUI/spheroid_pa_gui.py`) never talks to NIS-E
directly — it writes small INI **trigger** files, and these macros (run inside
NIS-E) act on them and write **done** files back. `spheroid_pipeline.py` (in
`../Codes/GUI/`) points `MACRO_DIR` at this folder, so the dispatcher resolves
every macro here with no hardcoded path.

## How they're launched

Two ways to run a macro in NIS-E:

1. **Dispatcher (preferred):** run **`nis_macro_dispatcher.mac`** once (Macro > Run,
   or flag it Run-on-Startup). It loops, reads `<work_dir>/cmd.ini` that the GUI
   writes, and `RunMacro()`s the matching macro below **in-process** — so you never
   hand-load a `.mac` per step again.
2. **Manual:** open the specific `.mac` in the NIS-E Macro editor and click Run.

## The macros

| Macro | Dispatcher `action_id` | Role | Function | Reads → Writes |
|---|---|---|---|---|
| `nis_macro_dispatcher.mac` | — (the router) | GUI↔NIS-E bridge | Started once; polls `cmd.ini` and `RunMacro()`s the matching macro, then writes `cmd_done.ini`. Preserves a step macro's own status (e.g. an abort). | `cmd.ini` → `cmd_done.ini` |
| `nis_macro_z_autofocus.mac` | 1 | Step 3 – autofocus | Per spheroid: move XY, two-pass autofocus (`StgFocusInRangeTwoPasses`), capture one plane, record the focus Z. | `af_trigger_NN.ini` → `af_done_NN.ini` + `<id>_af.nd2` |
| `nis_macro_capture_zstack.mac` | 2 | Step 3.5 – plain Z-stack | Per spheroid: move XY, run a fixed `ND` Z-series (no AF, no correction), save the stack. | `af_trigger_NN.ini` → `af_done_NN.ini` + `<id>_zstack.nd2` |
| `nis_macro_capture_zcorrected.mac` | 3 | Step 4 – Z-corrected capture | One at a time: load the spheroid's `.bin` (`ND_ZIntensityControlLoad`), run a Z-intensity-corrected Z-series (Beer-Lambert depth comp), save. | `spheroid_trigger.ini` (+`bin_path`) → `spheroid_done.ini` + nd2 |
| `nis_macro_pa_setup.mac` | 4 | Step 4 – PA prep (one-shot) | Clears the A1 interlock, selects the activation OC, dichroic OUT, sets the centred zoom square. Preps the rig for the `step3_zstack_PA` JOB. | `pa_trigger.ini` → (rig state) |
| `nis_macro_pa_points.mac` | 5 | Step 4 – PA batch points | Builds an ND multipoint (`ND_AppendMultipointPoint`) from EVERY currently-active `af_trigger_NN.ini` (i.e. Step 3's checked-row selection, whatever ranks it contains), for the job's "Import Point Set from ND." | `af_trigger_01..96` → ND multipoint |
| `nis_macro_pa_pick_current.mac` | 6 | Step 4 – PA single pick | Grabs the current live-centred stage X/Y/Z (`StgGetPos`) → `af_trigger_01` (`sph_pick_01`) + a 1-point multipoint, for manual single-spheroid PA. | (live stage) → `af_trigger_01.ini` |
| `nis_macro_pa_validate.mac` | 7 | Step 4 – PA validation | Post-PA 1050 nm re-image: select the viz OC (**laser-safety guard** — aborts if the active laser is still ~850 nm), Z-stack per spheroid to confirm the faded square. | `af_trigger_01..N` + `pa_trigger.ini` → `<id>_pacheck.nd2` |

## File-flag bridge

All paths come at runtime from `C:/SpheroidPA/session.ini` `[paths]`:
`work_dir` (trigger/done + `cmd.ini`), `nd2_dir` (captures), `macro_dir` (this folder).
The GUI writes triggers atomically (temp file + rename) with CRLF endings; macros
read INI keys with `Int_GetKeyString` + `atof`, delete the trigger after reading,
then write the done file — so the sequential queue is never clobbered.

## Rig dialect constraints (this NIS-E build)

- **single `main()` only** — no user-defined functions or globals (see
  `../macro_selftest/test_func.mac` / `test_proc.mac`).
- **CRLF + ASCII** files only — an LF-only or non-ASCII `.mac` fails silently.
- `Int_GetKeyString` + `atof` for INI values (`Int_GetKeyValue` is unreliable);
  the dispatcher selects by **numeric `action_id`** (avoids `strcmp`).
- `sprintf` takes its args as **one quoted comma-list of variable names**.
- `RunMacro(file)` runs a `.mac` in-process and returns when it completes.

## Related

- Per-function validators (stage, autofocus, capture, Z-intensity, Z-series) live in
  `../macro_selftest/`.
- Full change history: `../MD/CHANGELOG.md`. Pipeline overview: `../MD/README.md`.
