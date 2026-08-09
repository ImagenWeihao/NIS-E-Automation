r"""Did the PA fire where we told it to?

    python pa_verify_points.py <run_dir>

Reads the commanded points out of the pipeline log and the ACTUAL fired positions out of
the PA nd2 files' own stagePositionUm, then matches them up.

This is the only check that catches a silent fallback to the job's STORED PredefinedPoints
-- the failure mode where Jobs_RunJobInitParam rejects the whole parameter set, the job
runs on stale coordinates, and every log line still reports success. It would have caught
the 2026-08-07 PA that landed on the wrong spheroid at 17:21, at the time rather than days
later.

For a multi-point run it also proves the POINT LOOP visited every point rather than firing
N times at one of them -- distinct file count and distinct positions are different claims.
"""
import glob
import re
import sys
from pathlib import Path

import nd2

RUN = Path(sys.argv[1])
LOG = Path(r"C:/SpheroidPA/logs/spheroidpa_gui.log")
TOL_UM = 5.0        # commanded-vs-fired tolerance; the stage repeats to ~2 um

# ── commanded: the last "PA point params" block in the log ──────────────────
cmd = []
block = []
for ln in LOG.read_text(encoding="utf-8", errors="replace").splitlines():
    if "PA point params:" in ln:
        block = []
    m = re.search(r"\[(\d+)\]\s+(\S+):\s+Stage\.x/y=\(([-\d.]+),\s*([-\d.]+)\)\s*um\s+"
                  r"z=([-\d.]+)", ln)
    if m:
        block.append((int(m.group(1)), m.group(2), float(m.group(3)),
                      float(m.group(4)), float(m.group(5))))
    if "dispatched run_job" in ln and block:
        cmd = block

# ── fired: stagePositionUm out of the PA series itself ─────────────────────
#
# GROUPED BY Point####, not one row per file. A multi-plane PA writes one nd2 per
# z-plane, all at the same XY, so a file-per-point assumption reports a correct
# 2-point x 37-plane run as "72 files matched no commanded point" and "2 distinct
# positions of 74 files" -- both alarming, both wrong. What the job emits is
# Point0000_ZStack0000..36 / Point0001_ZStack0000..36. The claim to test is ONE XY per
# POINT, with z sweeping by design.
_by_point = {}
for p in sorted(glob.glob(str(RUN / "pa" / (sys.argv[2] if len(sys.argv) > 2 else "*ZStack*.nd2")))):
    m = re.search(r"(Point\d+)", Path(p).name)
    key = m.group(1) if m else Path(p).name
    with nd2.ND2File(p) as f:
        sp = f.frame_metadata(0).channels[0].position.stagePositionUm
        rec = (sp.x, sp.y, sp.z, f.sizes["X"] * f.voxel_size().x, dict(f.sizes))
    _by_point.setdefault(key, []).append(rec)

fired = []
for key, recs in sorted(_by_point.items()):
    xs = {round(r[0], 1) for r in recs}
    ys = {round(r[1], 1) for r in recs}
    zs = sorted(r[2] for r in recs)
    label = f"{key} ({len(recs)} plane{'s' if len(recs) > 1 else ''}, z {zs[0]:.0f}-{zs[-1]:.0f})"
    # XY must be constant within a point; z sweeps. XY moving mid-stack means the stage
    # drifted or the job re-homed -- surface it rather than averaging it away.
    if len(xs) > 1 or len(ys) > 1:
        label += f"  *** XY NOT CONSTANT ({len(xs)} x, {len(ys)} y) ***"
    # match on the MIDDLE plane's z: a symmetric stack straddles the commanded z, so the
    # first plane is a half-range low and would report a spurious dZ of -half.
    fired.append((label, recs[0][0], recs[0][1], zs[len(zs) // 2], recs[0][3], recs[0][4]))

print(f"commanded : {len(cmd)} point(s)")
for i, n, x, y, z in cmd:
    print(f"   [{i}] {n:24s} ({x:9.1f}, {y:9.1f})  z={z:8.1f}")
print(f"\nfired     : {len(fired)} PA point(s), {sum(len(v) for v in _by_point.values())} file(s)")
for n, x, y, z, side, sizes in fired:
    print(f"   {n[:58]:58s} ({x:9.1f}, {y:9.1f})  z_mid={z:8.1f}  {side:.1f} um sq")

if not cmd or not fired:
    print("\nCannot verify: missing commanded points or PA output.")
    sys.exit(1)

# ── match ──────────────────────────────────────────────────────────────────
print(f"\nmatching (tolerance {TOL_UM} um):")
used, ok = set(), True
for i, n, x, y, z in cmd:
    best, bd = None, 9e9
    for j, (fn, fx, fy, fz, _, _) in enumerate(fired):
        if j in used:
            continue
        d = ((fx - x) ** 2 + (fy - y) ** 2) ** 0.5
        if d < bd:
            best, bd = j, d
    if best is None:
        print(f"   [{i}] {n}: NO PA FILE LEFT TO MATCH  <-- point never fired")
        ok = False
        continue
    used.add(best)
    fn, fx, fy, fz, _, _ = fired[best]
    dz = fz - z
    flag = "ok" if bd <= TOL_UM else "*** MISMATCH ***"
    print(f"   [{i}] {n:22s} -> {fn[:44]:44s} dXY={bd:7.1f} um  dZ={dz:+7.1f} um  {flag}")
    if bd > TOL_UM:
        ok = False

extra = [fired[j][0] for j in range(len(fired)) if j not in used]
if extra:
    print(f"\n   {len(extra)} PA file(s) matched no commanded point:")
    for e in extra:
        print(f"      {e}")
    print("   -> the job fired somewhere it was not told to. Check its stored points.")
    ok = False

# distinct positions, not just distinct files
if len(fired) > 1:
    pos = {(round(f[1], 1), round(f[2], 1)) for f in fired}
    print(f"\n   distinct fired XY positions: {len(pos)} of {len(fired)} file(s)")
    if len(pos) < len(fired):
        print("   -> the POINT LOOP fired repeatedly at the SAME place, not at N points.")
        ok = False

print(f"\n{'VERIFIED' if ok else 'FAILED'}: commanded points vs where the laser actually went")
sys.exit(0 if ok else 1)
