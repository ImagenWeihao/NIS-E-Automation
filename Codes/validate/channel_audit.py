r"""Is the post-PA rise specific to the 940 nm PA readout, or does everything brighten?

For every run with a complete 3-channel pre+post set, report per channel:
  - acquisition settings pre vs post (laser power, PMT HV, offset, pinhole, zoom,
    line-average, scan speed) -- any difference invalidates the comparison outright
  - whole-stack mean, median, p99.9, max, and voxels > 1200 on the 0-4000 LUT
  - the same restricted to inside the PA square vs outside, so a global gain-like shift
    is separable from a localised one

The discriminator: 850 nm photoactivation converts PAsfGFP (940 readout). It does NOT
convert mBeRFP (890, the identity marker). If 890 rises as much as or more than 940, the
change is not photoconversion.
"""
import glob
import re
import sys
from pathlib import Path

import numpy as np
import nd2

BASE = Path(r"S:/Images/Weihao/NISeA/NIS-E-Automation/work")
TAGS = ["890nm_mBeRFP", "940nm_PAsfGFP", "1050nm_depth"]
ROLE = {"890nm_mBeRFP": "identity  (must NOT convert)",
        "940nm_PAsfGFP": "PA readout (should convert)",
        "1050nm_depth": "structure  (must NOT convert)"}
SETTING_RE = re.compile(
    r"\{(Laser Power|PMT HV|PMT Offset|Pinhole Size\(um\)|Scanner Zoom|Scan Speed|"
    r"Line Average/Integrate Count|Line Average Mode|Detector Selection)\}:\s*([^\t\r\n]+)")
PA_SIDE_DEFAULT = 79.5


def read(run, phase, tag):
    fs = glob.glob(str(BASE / run / "nd2" / f"{phase}_{tag}" / "*.nd2"))
    if not fs:
        return None
    with nd2.ND2File(fs[0]) as f:
        a = f.asarray().astype(np.float32)
        px = f.voxel_size().x
        desc = f.text_info.get("description", "")
        nz = a.shape[0] if a.ndim == 3 else 1
    settings = {}
    for k, v in SETTING_RE.findall(desc):
        settings.setdefault(k, v.strip())
    if a.ndim == 2:
        a = a[None]
    return a, px, nz, settings


def masks(shape, px, side=PA_SIDE_DEFAULT):
    _, ny, nx = shape
    Y, X = np.mgrid[:ny, :nx]
    S = side / px
    inside = np.maximum(np.abs(X - nx / 2), np.abs(Y - ny / 2)) <= S / 2
    return inside, ~inside


runs = sys.argv[1:] or ["0807", "0807_2", "0807_3"]
for run in runs:
    print("=" * 96)
    print(f"RUN {run}")
    got = {}
    for tag in TAGS:
        a, b = read(run, "prePA", tag), read(run, "postPA", tag)
        if a and b:
            got[tag] = (a, b)
    if not got:
        print("  no complete pre+post set")
        continue

    print("\n  acquisition settings (a mismatch here invalidates the comparison)")
    settings_ok = True
    for tag, ((_, _, nzp, sp), (_, _, nzq, sq)) in got.items():
        diffs = [f"{k}: {sp.get(k)} -> {sq.get(k)}"
                 for k in set(sp) | set(sq) if sp.get(k) != sq.get(k)]
        planes = f"{nzp}p -> {nzq}p" + ("  *** PLANE COUNT DIFFERS ***" if nzp != nzq else "")
        print(f"    {tag:16s} {planes:<34s} " +
              ("identical" if not diffs else "DIFFERS: " + "; ".join(diffs)))
        if diffs or nzp != nzq:
            settings_ok = False

    print(f"\n  {'channel':16s} {'role':28s} {'mean':>16s} {'>1200 voxels':>22s}")
    res = {}
    for tag, ((pa_, px, _, _), (pb, _, _, _)) in got.items():
        n = min(pa_.shape[0], pb.shape[0])
        pa_, pb = pa_[:n], pb[:n]
        ins, out = masks(pa_.shape, px)
        hi_a, hi_b = int((pa_ > 1200).sum()), int((pb > 1200).sum())
        res[tag] = dict(
            mean=(pa_.mean(), pb.mean()), hi=(hi_a, hi_b),
            ins=(pa_[:, ins].mean(), pb[:, ins].mean()),
            out=(pa_[:, out].mean(), pb[:, out].mean()))
        r = res[tag]
        print(f"    {tag:16s} {ROLE[tag]:28s} "
              f"{r['mean'][0]:7.2f}->{r['mean'][1]:7.2f} "
              f"{hi_a:9d}->{hi_b:9d} ({hi_b/max(hi_a,1):6.1f}x)")

    print(f"\n  {'channel':16s} {'inside PA square':>26s} {'outside (control)':>26s} {'in/out':>9s}")
    for tag, r in res.items():
        gi = r["ins"][1] / max(r["ins"][0], 1e-6)
        go = r["out"][1] / max(r["out"][0], 1e-6)
        print(f"    {tag:16s} {r['ins'][0]:8.2f}->{r['ins'][1]:8.2f} ({gi:5.2f}x) "
              f"{r['out'][0]:8.2f}->{r['out'][1]:8.2f} ({go:5.2f}x) {gi/max(go,1e-6):8.2f}")

    if "940nm_PAsfGFP" in res and "890nm_mBeRFP" in res:
        g940 = res["940nm_PAsfGFP"]["mean"][1] / max(res["940nm_PAsfGFP"]["mean"][0], 1e-6)
        g890 = res["890nm_mBeRFP"]["mean"][1] / max(res["890nm_mBeRFP"]["mean"][0], 1e-6)
        verdict = ("SPECIFIC -- 940 rose more than the identity channel"
                   if g940 > g890 * 1.5 else
                   "NOT SPECIFIC -- the identity channel rose as much or more")
        print(f"\n  whole-stack gain  940={g940:.2f}x  890={g890:.2f}x  ->  {verdict}")
    if not settings_ok:
        print("  NOTE: settings/plane-count mismatch above -- treat the numbers as unusable.")
    print()
