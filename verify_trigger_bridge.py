"""
verify_trigger_bridge.py — confirm the Python <-> NIS-E file-flag bridge works.

The NIS-E macros read trigger INIs with Int_GetKeyString / Int_GetKeyValue and
write done INIs with Int_SetKeyString. Those functions wrap the Windows profile
API (GetPrivateProfileString / WritePrivateProfileString). This script exercises
the REAL spheroid_pipeline writers, then reads the files back through that exact
Win32 API — so a PASS here means the actual NIS-E macros will read the same
values. It also drives the reverse direction (Win32 writes the done file, the
Python pollers read it).

No NIS-E required. Run:  python verify_trigger_bridge.py
Exit code 0 = bridge OK, 1 = mismatch found.
"""
from __future__ import annotations

import ctypes
import sys
import tempfile
from ctypes import wintypes
from pathlib import Path

import spheroid_pipeline as pl

# ── Win32 profile API (exactly what Int_Get/SetKey* wrap inside NIS-E) ─────────
_k32 = ctypes.WinDLL("kernel32", use_last_error=True)

_GetStr = _k32.GetPrivateProfileStringW
_GetStr.argtypes = [wintypes.LPCWSTR, wintypes.LPCWSTR, wintypes.LPCWSTR,
                    wintypes.LPWSTR, wintypes.DWORD, wintypes.LPCWSTR]
_GetStr.restype = wintypes.DWORD

_SetStr = _k32.WritePrivateProfileStringW
_SetStr.argtypes = [wintypes.LPCWSTR, wintypes.LPCWSTR, wintypes.LPCWSTR,
                    wintypes.LPCWSTR]
_SetStr.restype = wintypes.BOOL


def win_get_string(path: Path, section: str, key: str, default: str = "") -> str:
    """Mirror of NIS-E Int_GetKeyString(path, section, key, buf, len)."""
    buf = ctypes.create_unicode_buffer(2048)
    # GetPrivateProfileString needs an ABSOLUTE path or it searches %WINDIR%.
    _GetStr(section, key, default, buf, len(buf), str(path.resolve()))
    return buf.value


def win_get_value(path: Path, section: str, key: str, default: str = "0") -> float:
    """Mirror of NIS-E Int_GetKeyValue(path, section, key, default) -> number.

    NIS-E returns a numeric (atof of the stored string); floats like 123.4 rule
    out GetPrivateProfileInt, so the value path is string-fetch + atof.
    """
    return float(win_get_string(path, section, key, default))


def win_set_string(path: Path, section: str, key: str, value: str) -> None:
    """Mirror of NIS-E Int_SetKeyString(path, section, key, value)."""
    if not _SetStr(section, key, value, str(path.resolve())):
        raise OSError(f"WritePrivateProfileString failed: {ctypes.get_last_error()}")


# ── reporting helpers ─────────────────────────────────────────────────────────
_fails: list[str] = []


def check(label: str, got, want, tol: float | None = None) -> None:
    if tol is not None:
        ok = abs(float(got) - float(want)) <= tol
    else:
        ok = str(got) == str(want)
    mark = "OK  " if ok else "FAIL"
    print(f"  [{mark}] {label}: got={got!r} want={want!r}")
    if not ok:
        _fails.append(label)


def check_crlf(path: Path) -> None:
    raw = path.read_bytes()
    has_crlf = b"\r\n" in raw
    lone_lf = raw.replace(b"\r\n", b"").count(b"\n")
    ok = has_crlf and (lone_lf == 0)
    mark = "OK  " if ok else "FAIL"
    print(f"  [{mark}] {path.name} line endings: CRLF={has_crlf} lone_LF={lone_lf}")
    if not ok:
        _fails.append(f"{path.name} line endings")


def make_record(rank: int = 7) -> pl.SpheroidRecord:
    r = pl.SpheroidRecord(
        spheroid_id=f"WellD05_sph{rank:02d}", rank=rank, well="D05",
        mosaic_cx_px=1234.5, mosaic_cy_px=678.9, mosaic_diam_um=212.0,
        mosaic_score=0.947, mosaic_x_um=15234.5, mosaic_y_um=-3415.0,
        mosaic_z_um=380.0,
    )
    # Values the bridge actually transmits (set by Steps 2-3 in the GUI):
    r.verified_x_um = 15234.56
    r.verified_y_um = -3415.07
    r.z_centre_um   = 384.25
    r.z_half_um     = 18.0
    r.z_step_um     = 2.0
    r.bin_path      = r"S:\Images\Weihao\NISeA\NIS-E-Automation\work\bins\WellD05_sph07.bin"
    return r


# ── Bridge 1: Step 4 capture (spheroid_trigger.ini / spheroid_done.ini) ───────
def test_capture_bridge(work: Path) -> None:
    print("\n[Bridge 1] Step 4 capture  (nis_macro_auto_capture.mac)")
    rec = make_record()
    nd2_out_dir = work / "nd2"
    trigger = pl.write_trigger(rec, work, nd2_out_dir)
    check_crlf(trigger)

    print("  -- macro-side reads (Int_GetKeyString / Int_GetKeyValue) --")
    SEC = "spheroid"
    check("rank (str)",        win_get_string(trigger, SEC, "rank"),        str(rec.rank))
    check("spheroid_id (str)", win_get_string(trigger, SEC, "spheroid_id"), rec.spheroid_id)
    check("bin_path (str)",    win_get_string(trigger, SEC, "bin_path"),    rec.bin_path)
    check("nd2_out (str)",     win_get_string(trigger, SEC, "nd2_out"),     str(nd2_out_dir / f"{rec.spheroid_id}.nd2"))
    check("stage_x (val)",     win_get_value(trigger, SEC, "stage_x"),  rec.verified_x_um, tol=1e-6)
    check("stage_y (val)",     win_get_value(trigger, SEC, "stage_y"),  rec.verified_y_um, tol=1e-6)
    check("z_centre (val)",    win_get_value(trigger, SEC, "z_centre"), rec.z_centre_um,   tol=1e-6)
    check("z_half (val)",      win_get_value(trigger, SEC, "z_half"),   rec.z_half_um,     tol=1e-6)
    check("z_step (val)",      win_get_value(trigger, SEC, "z_step"),   rec.z_step_um,     tol=1e-6)

    # Reverse: macro writes spheroid_done.ini (via Int_SetKeyString), Python reads it.
    print("  -- reverse: macro writes done.ini, Python wait_for_done() reads --")
    done = work / pl.DONE_FILENAME
    out_nd2 = str(nd2_out_dir / f"{rec.spheroid_id}.nd2")
    win_set_string(done, SEC, "rank", str(rec.rank))
    win_set_string(done, SEC, "spheroid_id", rec.spheroid_id)
    win_set_string(done, SEC, "nd2_path", out_nd2)
    win_set_string(done, SEC, "status", "ok")
    parsed = pl.wait_for_done(work, timeout_s=2.0, poll_s=0.1)
    check("done status",   (parsed or {}).get("status"),   "ok")
    check("done nd2_path", (parsed or {}).get("nd2_path"), out_nd2)


# ── Bridge 2: Step 3 autofocus (af_trigger_NN.ini / af_done_NN.ini) ───────────
def test_autofocus_bridge(work: Path) -> None:
    print("\n[Bridge 2] Step 3 autofocus  (nis_macro_z_autofocus.mac)")
    recs = [make_record(rank=3), make_record(rank=12)]
    n = pl.trigger_autofocus_all(recs, work)
    check("trigger count", n, len(recs))

    SEC = "spheroid"
    for rec in recs:
        trig = work / f"{pl.AF_TRIGGER_PREFIX}{rec.rank:02d}.ini"
        print(f"  -- af_trigger_{rec.rank:02d}.ini macro-side reads --")
        check_crlf(trig)
        check(f"r{rec.rank} rank",    win_get_string(trig, SEC, "rank"),        str(rec.rank))
        check(f"r{rec.rank} id",      win_get_string(trig, SEC, "spheroid_id"), rec.spheroid_id)
        check(f"r{rec.rank} stage_x", win_get_value(trig, SEC, "stage_x"), rec.verified_x_um, tol=1e-6)
        check(f"r{rec.rank} stage_y", win_get_value(trig, SEC, "stage_y"), rec.verified_y_um, tol=1e-6)

    # Reverse: macro writes af_done_NN.ini, Python poll_autofocus_done() reads it.
    print("  -- reverse: macro writes af_done, Python poll_autofocus_done() reads --")
    for rec in recs:
        done = work / f"{pl.AF_DONE_PREFIX}{rec.rank:02d}.ini"
        win_set_string(done, SEC, "rank", str(rec.rank))
        win_set_string(done, SEC, "spheroid_id", rec.spheroid_id)
        win_set_string(done, SEC, "z_centre", f"{rec.z_centre_um:.3f}")
        win_set_string(done, SEC, "nd2_path", f"{rec.spheroid_id}_af.nd2")
        win_set_string(done, SEC, "status", "ok")
    got = pl.poll_autofocus_done(recs, work)
    for rec in recs:
        check(f"r{rec.rank} af z_centre", got.get(rec.rank, {}).get("z_centre_um"),
              rec.z_centre_um, tol=1e-3)


def main() -> int:
    print("=" * 70)
    print("NIS-E trigger bridge verification (Win32 profile API == macro reads)")
    print("=" * 70)
    work = Path(tempfile.mkdtemp(prefix="spa_bridge_"))
    print(f"scratch dir: {work}")
    try:
        test_capture_bridge(work)
        test_autofocus_bridge(work)
    finally:
        for p in work.rglob("*"):
            if p.is_file():
                p.unlink()
        # leave the (now-empty) tree for inspection if anything failed
        if not _fails:
            import shutil
            shutil.rmtree(work, ignore_errors=True)

    print("\n" + "=" * 70)
    if _fails:
        print(f"RESULT: FAIL ({len(_fails)} mismatch) -> {', '.join(_fails)}")
        print(f"scratch kept for inspection: {work}")
        return 1
    print("RESULT: PASS -- every key the NIS-E macros read round-trips correctly.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
