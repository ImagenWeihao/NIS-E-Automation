r"""Archive one session's logs into <run_dir>/Log/, sliced to the session.

    python collect_session_logs.py <run_dir>

Reads <run_dir>/Log/session_state/session_start.ini for the byte offsets recorded when
the session opened, and copies only what was written after them. That is the point of
the offsets: spheroidpa_gui.log and the NIS-E lxapp_*.log are BOTH cumulative -- one
lxapp file covers a whole day across every run -- so copying them wholesale produces a
120 MB archive in which this session is indistinguishable from the six before it.

Also snapshots the handshake/config files as they stand now, and the macro readback
files from C:\SpheroidPA\debug\.

Safe to run repeatedly; it overwrites its own outputs and touches nothing else.
"""
import configparser
import shutil
import sys
from pathlib import Path

RUN = Path(sys.argv[1])
LOGD = RUN / "Log"
START = LOGD / "session_state" / "session_start.ini"

if not START.exists():
    sys.exit(f"No {START} -- this session was not opened with a start marker.")

cp = configparser.ConfigParser()
cp.read_string("[s]\n" + START.read_text(encoding="utf-8", errors="replace"))
s = cp["s"]

DEBUG_DIR = Path(r"C:/SpheroidPA/debug")
WORKDIRS = [RUN / "autofocus", Path(r"C:/SpheroidPA")]
STATE_FILES = ["session.ini", "plate_calib.ini", "cmd.ini", "cmd_done.ini",
               "pa_trigger.ini", "job_trigger.ini", "init_trigger.ini",
               "lightpath_trigger.ini", "lightpath_state.ini", "abort.ini"]


def slice_after(src, nbytes, dest):
    """Copy src[nbytes:] to dest. If the file shrank (rotated), copy the whole thing."""
    src = Path(src)
    if not src.exists():
        return f"MISSING {src}"
    size = src.stat().st_size
    if size < nbytes:                      # rotated out from under us
        nbytes = 0
        note = " (log rotated -- copied whole)"
    else:
        note = ""
    with src.open("rb") as f:
        f.seek(nbytes)
        data = f.read()
    dest.write_bytes(data)
    return f"{dest.name}: {len(data)/1e6:.2f} MB of {size/1e6:.2f} MB{note}"


print(f"session opened {s.get('session_start_local', '?')}\n")

print("pipeline log")
print("  " + slice_after(s["spheroidpa_log"], int(s["spheroidpa_bytes_at_start"]),
                         LOGD / "spheroidpa" / "spheroidpa_gui_session.log"))
shutil.copy2(s["spheroidpa_log"], LOGD / "spheroidpa" / "spheroidpa_gui_FULL.log")
print(f"  spheroidpa_gui_FULL.log: complete history, for cross-run comparison")

print("\nNIS-E log")
print("  " + slice_after(s["nise_log"], int(s["nise_bytes_at_start"]),
                         LOGD / "nise" / "lxapp_session.log"))
# Any lxapp file touched after the session opened -- NIS-E rolls over on restart, so the
# session can span several, and the numbered ones are near-empty restart markers whose
# timestamps are the useful part.
ni = Path(s["nise_log"])
for f in sorted(ni.parent.glob("lxapp_*.log")):
    if f != ni and f.stat().st_mtime >= ni.stat().st_mtime - 86400:
        shutil.copy2(f, LOGD / "nise" / f.name)
        print(f"  also copied {f.name} ({f.stat().st_size/1e6:.2f} MB)")

print("\nmacro readbacks")
if DEBUG_DIR.is_dir():
    for f in DEBUG_DIR.glob("*.txt"):
        shutil.copy2(f, LOGD / "spheroidpa" / f.name)
        print(f"  {f.name}")
else:
    print(f"  {DEBUG_DIR} not present")

print("\nsession state")
for d in WORKDIRS:
    if not d.is_dir():
        continue
    for name in STATE_FILES:
        p = d / name
        if p.exists():
            shutil.copy2(p, LOGD / "session_state" / p.name)
            print(f"  {p}")
    for p in sorted(d.glob("af_trigger_*.ini")) + sorted(d.glob("af_done_*.ini")):
        shutil.copy2(p, LOGD / "session_state" / p.name)
        print(f"  {p.name}")
for extra in ("pipeline_state.json", "spheroid_screen_latest.csv"):
    p = RUN / extra
    if p.exists():
        shutil.copy2(p, LOGD / "session_state" / p.name)
        print(f"  {p.name}")

total = sum(f.stat().st_size for f in LOGD.rglob("*") if f.is_file())
n = sum(1 for f in LOGD.rglob("*") if f.is_file())
print(f"\n{LOGD}\n  {n} files, {total/1e6:.1f} MB")
