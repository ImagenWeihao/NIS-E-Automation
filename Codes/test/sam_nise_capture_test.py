#!/usr/bin/env python3
"""
sam_nise_capture_test.py -- First SAM-driven NIS-E capture test  [NISE017]

Drive the Nikon NIS-Elements rig to MOVE to a target stage XY and CAPTURE one
image, through the existing SpheroidPA file-trigger bridge:

    session.ini + af_trigger_01.ini  --> nis_macro_capture_zstack.mac
                                     <-- af_done_01.ini + <id>_zstack.nd2

This is the foundation layer the SAM<->NIS-E hardware adapter will wrap: the
move+capture primitive here is exactly what SAM's move_stage / capture modules
will call on this rig (see [NISE017] plan, Phase 2).

MODES
  --simulate   Hardware-free round-trip: write the trigger, then play the
               macro's part (write af_done with a placeholder .nd2). Proves the
               Python protocol end-to-end with NO microscope. Run this first,
               anywhere -- it must PASS before going to the scope.
  (live)       Write session.ini + af_trigger_01.ini, wait for you to run
               nis_macro_capture_zstack.mac in NIS-E, poll af_done_01.ini,
               confirm the .nd2 was written.

USAGE
  python sam_nise_capture_test.py --simulate
  python sam_nise_capture_test.py --x 64910 --y 14710 --z 7680

LIVE RIG PREREQS (END_TO_END_TEST.md sec 0.3)
  20X objective on turret; XY+Z stages initialized; Z drive NOT in Escape mode;
  camera NOT in Live; PFS OFF (it holds focus and defeats single-plane capture).
  nis_macro_capture_zstack.mac must be CRLF + ASCII (an LF-only .mac silently
  no-ops). Get --x/--y by centering the target well's spheroid in NIS-E and
  reading the stage position, or from the SpheroidPA screener's verified_x/y_um.
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

# The NIS-E macros read their working paths from this fixed pointer file.
DEFAULT_SESSION_INI = Path("C:/SpheroidPA/session.ini")
SECTION = "spheroid"
TRIGGER_NAME = "af_trigger_01.ini"   # rank 01 (matches nis_macro_capture_zstack.mac NN loop)
DONE_NAME = "af_done_01.ini"
MACRO_NAME = "nis_macro_capture_zstack.mac"
# abspath (NOT resolve): keep a mapped drive as "S:/..." for NIS-E RunMacro --
# resolve() would expand it to a \\host\share UNC that RunMacro may reject.
# This script lives in Codes/test/, so go up two levels to the repo root, then macro/.
MACRO_DIR = Path(os.path.abspath(__file__)).parent.parent.parent / "macro"  # repo-root macro/ folder


def _atomic_write_crlf(path: Path, lines: list[str]) -> None:
    """Write an INI with explicit CRLF, atomically (mirrors spheroid_pipeline).

    CRLF: NIS-E's Int_GetKeyString wraps Win32 GetPrivateProfileString, which
    only parses Windows line endings reliably. ASCII: the keys/paths are ASCII;
    enforce it so no stray byte breaks the parse. Atomic: the macro polls by
    name, so os.replace guarantees it never sees a half-written file.
    """
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text("\r\n".join(lines) + "\r\n", encoding="ascii", newline="")
    os.replace(tmp, path)


def _parse_ini(text: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for ln in text.splitlines():
        ln = ln.strip()
        if "=" in ln and not ln.startswith("["):
            k, v = ln.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def write_session(session_ini: Path, work_dir: Path, nd2_dir: Path) -> None:
    session_ini.parent.mkdir(parents=True, exist_ok=True)
    _atomic_write_crlf(session_ini, [
        "[paths]",
        f"work_dir={work_dir.as_posix()}",
        f"nd2_dir={nd2_dir.as_posix()}",
        f"macro_dir={MACRO_DIR.as_posix()}",
    ])


def write_trigger(work_dir: Path, sphid: str, x: float, y: float, z: float,
                  z_half: float = 0.0, z_step: float = 10.0) -> Path:
    work_dir.mkdir(parents=True, exist_ok=True)
    # clear any stale rank-01 trigger/done so the macro only sees this run
    (work_dir / TRIGGER_NAME).unlink(missing_ok=True)
    (work_dir / DONE_NAME).unlink(missing_ok=True)
    trig = work_dir / TRIGGER_NAME
    _atomic_write_crlf(trig, [
        f"[{SECTION}]",
        "rank=1",
        f"spheroid_id={sphid}",
        f"stage_x={x}",
        f"stage_y={y}",
        f"z_centre={z}",
        f"z_half={z_half}",   # 0 -> single-plane locate (StgMoveZ + Capture); >0 -> Z-stack
        f"z_step={z_step}",
    ])
    return trig


def poll_done(work_dir: Path, timeout_s: float = 300.0, interval: float = 2.0) -> dict | None:
    done = work_dir / DONE_NAME
    deadline = time.time() + timeout_s
    while time.time() < deadline:
        if done.exists():
            d = _parse_ini(done.read_text(encoding="ascii", errors="replace"))
            done.unlink(missing_ok=True)
            return d
        time.sleep(interval)
    return None


def simulate_macro(work_dir: Path, nd2_dir: Path, sphid: str) -> None:
    """Play nis_macro_capture_zstack.mac's part: consume trigger, write done + placeholder nd2."""
    trig = work_dir / TRIGGER_NAME
    if not trig.exists():
        raise RuntimeError("simulate: no trigger present to consume")
    d = _parse_ini(trig.read_text(encoding="ascii", errors="replace"))
    nd2_dir.mkdir(parents=True, exist_ok=True)
    fake_nd2 = nd2_dir / f"{sphid}_zstack.nd2"
    fake_nd2.write_bytes(b"")                 # placeholder, mimics ImageSaveAs(path,14,0)
    trig.unlink(missing_ok=True)             # the real macro deletes the trigger after handling
    _atomic_write_crlf(work_dir / DONE_NAME, [
        f"[{SECTION}]",
        "rank=1",
        f"spheroid_id={sphid}",
        f"z_centre={d.get('z_centre', '0')}",
        f"nd2_path={fake_nd2.as_posix()}",
        "status=ok",
    ])


def main() -> int:
    ap = argparse.ArgumentParser(description="First SAM-driven NIS-E move+capture test [NISE017]")
    ap.add_argument("--x", type=float, help="target stage X (um, absolute)")
    ap.add_argument("--y", type=float, help="target stage Y (um, absolute)")
    ap.add_argument("--z", type=float, default=7680.0, help="focus Z (um); macro default 7680")
    ap.add_argument("--id", default="nise017_test", help="capture id / nd2 stem")
    ap.add_argument("--z-half", type=float, default=0.0,
                    help="0 = single plane (default); >0 = Z-stack half-range um")
    ap.add_argument("--session-ini", default=str(DEFAULT_SESSION_INI),
                    help="pointer file the macros read (default C:/SpheroidPA/session.ini)")
    ap.add_argument("--work-dir", default="C:/SpheroidPA/work/autofocus",
                    help="folder where trigger/done INIs live (must match GUI Step 1 save dir)")
    ap.add_argument("--nd2-dir", default="C:/SpheroidPA/work/nd2",
                    help="folder where the .nd2 is saved")
    ap.add_argument("--timeout", type=float, default=300.0, help="seconds to wait for af_done_01.ini")
    ap.add_argument("--simulate", action="store_true", help="hardware-free round-trip (no scope)")
    a = ap.parse_args()

    if not a.simulate and (a.x is None or a.y is None):
        ap.error("--x and --y (target stage coords in um) are required for a live run "
                 "(or use --simulate).")

    session_ini = Path(a.session_ini)
    work_dir = Path(a.work_dir)
    nd2_dir = Path(a.nd2_dir)
    x = a.x if a.x is not None else 64910.0   # SLIM050 example coord
    y = a.y if a.y is not None else 14710.0

    print(f"[NISE017] macro_dir = {MACRO_DIR}")
    print(f"[NISE017] session.ini -> {session_ini}")
    write_session(session_ini, work_dir, nd2_dir)
    print(f"[NISE017] trigger -> {work_dir / TRIGGER_NAME}")
    print(f"          target x={x} y={y} z={a.z} z_half={a.z_half} ({'Z-stack' if a.z_half > 0 else 'single plane'})")
    write_trigger(work_dir, a.id, x, y, a.z, z_half=a.z_half)

    if a.simulate:
        print("[NISE017] SIMULATE: playing the NIS-E macro's part (no microscope)...")
        simulate_macro(work_dir, nd2_dir, a.id)
        wait = min(10.0, a.timeout)
    else:
        print("")
        print(f"  >>> In NIS-E now:  Macro -> Run  {MACRO_NAME}")
        print("  >>> (20X obj; Z NOT in Escape; camera NOT in Live; PFS OFF)")
        print("")
        wait = a.timeout

    print(f"[NISE017] polling for {DONE_NAME} (timeout {wait:.0f}s)...")
    d = poll_done(work_dir, timeout_s=wait)
    if d is None:
        print("FAIL: no af_done_01.ini within timeout.")
        print("      Causes: macro not running / hit idle timeout; GUI path != macro work_dir;")
        print("      .mac saved with LF endings; Z in Escape mode. See END_TO_END_TEST.md sec 4.")
        return 2

    status = d.get("status", "?")
    nd2 = d.get("nd2_path", "")
    print(f"[NISE017] done: status={status}  nd2_path={nd2}")
    if status != "ok":
        print(f"FAIL: macro reported status={status} (e.g. move_failed -> stage move rejected).")
        return 3
    if nd2 and Path(nd2).exists():
        size = Path(nd2).stat().st_size
        print(f"PASS: image captured -> {nd2} ({size} bytes)")
        return 0
    print("WARN: status=ok but the .nd2 was not found at the reported path "
          "(check the control-PC vs rig shared-path mapping).")
    return 1


if __name__ == "__main__":
    sys.exit(main())
