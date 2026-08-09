r"""Native-resolution depth montages: one figure per spheroid per channel.

    python pa_depth_montage.py <run_dir> [pa_glob] [vmax]

Three rows -- Pre-PA / Pre-PA + PA zone / Post-PA -- and one column per z-plane, written
at SOURCE RESOLUTION: one image pixel per output pixel, no resampling anywhere. Zooming
in resolves data rather than interpolation, which is the whole point of these sheets.
PIL is used rather than matplotlib for exactly that reason: an imshow would resample.

DISPLAY RANGE is ONE value per figure, applied identically to every panel -- all three
rows, all N z-planes, pre and post alike. That is the property that matters: neither side
of the comparison is ever favoured, which is how a retracted result once got made.

The value is chosen PER CHANNEL, by the same rule the 0807 sheets used: snap the 99.9th
percentile over pre AND post together up to the next round ceiling. It is emphatically
not a per-panel stretch. A single ceiling shared across channels does not work on this
sample -- 940 nm PAsfGFP sits at a mean of ~2 counts while 890 nm runs to hundreds, so a
range that suits the identity channel renders the PA readout black and hides the result.
Reference values, 0807_6:  890 -> 0-4000,  940 -> 0-400,  1050 -> 0-2000.

Pass a number as the third argument to pin the range instead, when a figure has to be
compared against another run pixel-for-pixel. Either way the LUT key on each figure
states the range, so the level is never implied.

The PA square and its z range come from the PA series' OWN geometry and stage position --
side from sizes['X'] * voxel_size().x, centre matched to this spheroid's readout position.
Never from a GUI field or a config value: on run 0807_2 the pa/*.nd2 turned out to be byte
copies of the readout stacks, and the square used in that analysis had come from a config
number, which is what made its in/out statistic meaningless. If no PA file matches a
spheroid, the zone row says so instead of drawing an assumed box.
"""
import glob
import re
import sys
from pathlib import Path

import numpy as np
import nd2
from PIL import Image, ImageDraw, ImageFont

RUN = Path(sys.argv[1])
PA_GLOB = sys.argv[2] if len(sys.argv) > 2 else "*ZStack*.nd2"
# "auto" (default) reproduces the rule the 0807 sheets used: per channel, snap the 99.9th
# percentile over pre AND post together up to the next round number. A single fixed
# ceiling across all three channels does not work here -- 940 nm PAsfGFP sits at a mean of
# ~2 counts against 890 nm's hundreds, so a range that suits the identity channel renders
# the PA readout black, which hides the result rather than reporting it.
# Pass a number instead to pin the range when comparing across runs.
VMAX_ARG = sys.argv[3] if len(sys.argv) > 3 else "auto"

BG, INK, MUTED, CYAN, RED = (30, 30, 46), (232, 232, 240), (154, 160, 181), (34, 229, 229), (243, 139, 168)
GAP = 3
L_M, T_M, B_M, R_M = 430, 330, 150, 50
CHANNELS = [("890nm_mBeRFP", "890 nm  mBeRFP  -- T-cell identity", (255, 32, 32)),
            ("940nm_PAsfGFP", "940 nm  PAsfGFP  -- PA readout", (32, 255, 32)),
            ("1050nm_depth", "1050 nm  depth  -- spheroid structure", (255, 32, 32))]
ROW_LABELS = ["Pre-PA", "Pre-PA + PA zone", "Post-PA"]


def nice_max(v):
    """Snap up to the next round ceiling -- the 0807 rule, kept identical so sheets from
    the two runs are read the same way."""
    for c in (25, 50, 100, 200, 400, 800, 1200, 2000, 4000, 8000, 16000, 65535):
        if v <= c:
            return float(c)
    return 65535.0


def font(sz):
    for n in ("arial.ttf", "segoeui.ttf", "DejaVuSans.ttf"):
        try:
            return ImageFont.truetype(n, sz)
        except OSError:
            pass
    return ImageFont.load_default()


def load(phase, tag, sph):
    fs = glob.glob(str(RUN / "nd2" / f"{phase}_{tag}" / f"{sph}_zstack.nd2"))
    if not fs:
        return None
    with nd2.ND2File(fs[0]) as f:
        a = f.asarray().astype(np.float32)
        px = f.voxel_size().x
        sp = f.frame_metadata(0).channels[0].position.stagePositionUm
        lp = f.unstructured_metadata()["ImageMetadataLV"]["SLxExperiment"]["LoopPars"]
        n = a.shape[0] if a.ndim == 3 else 1
        z = float(lp["ZLow"]) + float(lp["ZStep"]) * np.arange(n)
    if a.ndim == 2:
        a = a[None]
    return a, px, z, (sp.x, sp.y)


def lut_key(d, x, y, colour, vmax, f_lbl):
    """Horizontal 0..vmax gradient in the channel's own LUT, with end labels.

    Drawn in the header margin, never over a panel -- an inset laid on the image would
    sit on tissue at some z even if it misses it at the plane you happened to check."""
    w, h = 260, 22
    for i in range(w):
        v = i / (w - 1)
        d.line([(x + i, y), (x + i, y + h)],
               fill=tuple(int(v * c + 0.5) for c in colour))
    d.rectangle([x, y, x + w, y + h], outline=MUTED, width=1)
    d.text((x, y + h + 5), "0", font=f_lbl, fill=MUTED, anchor="lt")
    d.text((x + w, y + h + 5), f"{vmax:.0f}", font=f_lbl, fill=MUTED, anchor="rt")
    d.text((x + w / 2, y + h + 5), "display range (identical every panel, every spheroid)",
           font=f_lbl, fill=MUTED, anchor="mt")


# ── PA stimulus geometry, grouped per point ────────────────────────────────
pa = {}
for p in sorted(glob.glob(str(RUN / "pa" / PA_GLOB))):
    key = (re.search(r"(Point\d+)", Path(p).name) or [None, Path(p).name])[1]
    with nd2.ND2File(p) as f:
        sp = f.frame_metadata(0).channels[0].position.stagePositionUm
        pa.setdefault(key, []).append((sp.x, sp.y, sp.z, f.sizes["X"] * f.voxel_size().x))
pa_pts = [(k, v[0][0], v[0][1], min(r[2] for r in v), max(r[2] for r in v), v[0][3])
          for k, v in sorted(pa.items())]
print(f"PA points: {len(pa_pts)}")
for k, x, y, z0, z1, s in pa_pts:
    print(f"   {k}  ({x:.1f}, {y:.1f})  z {z0:.0f}-{z1:.0f}  {s:.1f} um square")

spheroids = sorted({Path(p).stem.replace("_zstack", "")
                    for p in glob.glob(str(RUN / "nd2" / "prePA_*" / "*_zstack.nd2"))})
print(f"spheroids: {spheroids}\n")

# ── one display range per CHANNEL, over EVERY spheroid and both phases ──────
# Per-spheroid ranges would rescale each spheroid to its own brightness, so a spheroid
# that lost 85% of its signal renders as bright as one that did not -- the two become
# incomparable at a glance, which is the opposite of what these sheets are for.
CH_VMAX = {}
if VMAX_ARG == "auto":
    for _tag, _lbl, _c in CHANNELS:
        _hi = 0.0
        for _s in spheroids:
            for _ph in ("prePA", "postPA"):
                _r = load(_ph, _tag, _s)
                if _r:
                    _hi = max(_hi, float(np.percentile(_r[0], 99.9)))
        CH_VMAX[_tag] = nice_max(_hi)
        print(f"  display range {_tag:16s} 0-{CH_VMAX[_tag]:.0f}"
              f"   (max p99.9 over {len(spheroids)} spheroid(s), both phases: {_hi:.0f})")
    print()

for sph in spheroids:
    for tag, label, colour in CHANNELS:
        pre, post = load("prePA", tag, sph), load("postPA", tag, sph)
        if not pre or not post:
            print(f"  {sph} {tag}: incomplete pre/post, skipped")
            continue
        (A, px, zs, xy), (B, _, _, _) = pre, post
        nz = min(A.shape[0], B.shape[0])
        A, B, zs = A[:nz], B[:nz], zs[:nz]
        ny, nx = A.shape[1], A.shape[2]

        VMAX = CH_VMAX[tag] if VMAX_ARG == "auto" else float(VMAX_ARG)

        mine = min(pa_pts, key=lambda r: (r[1] - xy[0]) ** 2 + (r[2] - xy[1]) ** 2) if pa_pts else None
        if mine and ((mine[1] - xy[0]) ** 2 + (mine[2] - xy[1]) ** 2) ** 0.5 > 50:
            mine = None
        side, PA_LO, PA_HI = (mine[5], mine[3], mine[4]) if mine else (None, None, None)

        W = L_M + nz * (nx + GAP) - GAP + R_M
        H = T_M + 3 * (ny + GAP) - GAP + B_M
        canvas = np.empty((H, W, 3), np.uint8)
        canvas[:, :] = BG
        col = np.array(colour, np.float32) / 255.0
        for r in range(3):
            stack = B if r == 2 else A
            y = T_M + r * (ny + GAP)
            for i in range(nz):
                x = L_M + i * (nx + GAP)
                canvas[y:y + ny, x:x + nx] = (
                    (np.clip(stack[i] / VMAX, 0, 1)[..., None] * col) * 255 + 0.5).astype(np.uint8)

        img = Image.fromarray(canvas)
        d = ImageDraw.Draw(img)
        f_t, f_s, f_r, f_z, f_l = font(64), font(30), font(38), font(26), font(24)

        ry = T_M + (ny + GAP)
        if side:
            s_px = side / px
            for i in range(nz):
                x = L_M + i * (nx + GAP)
                d.rectangle([x + (nx - s_px) / 2, ry + (ny - s_px) / 2,
                             x + (nx + s_px) / 2, ry + (ny + s_px) / 2],
                            outline=CYAN, width=3)
        else:
            d.text((L_M + 12, ry + 12), "no PA file matches this spheroid",
                   font=f_s, fill=RED, anchor="lt")

        for i in range(nz):
            x = L_M + i * (nx + GAP)
            inpa = side is not None and PA_LO - 1 <= zs[i] <= PA_HI + 1
            if inpa:
                d.rectangle([x, T_M - 26, x + nx, T_M - 18], fill=CYAN)
            if i % 2 == 0 or i == nz - 1:
                d.text((x + nx / 2, T_M - 36), f"{zs[i]:.0f}", font=f_z,
                       fill=CYAN if inpa else MUTED, anchor="mb")
        for r, lab in enumerate(ROW_LABELS):
            d.text((L_M - 22, T_M + r * (ny + GAP) + ny / 2), lab, font=f_r, fill=INK, anchor="rm")

        d.text((L_M, 46), f"{sph}   {label}", font=f_t, fill=INK, anchor="lt")
        pa_txt = (f"cyan = 850 nm PA, {side:.0f} um square, z {PA_LO:.0f}-{PA_HI:.0f} um"
                  if side else "PA geometry not identifiable in pa/")
        d.text((L_M, 130), f"{nz} z-planes {zs[0]:.0f}-{zs[-1]:.0f} um, {zs[1]-zs[0]:.0f} um step"
               f"   |   {pa_txt}", font=f_s, fill=MUTED, anchor="lt")
        d.text((L_M, 172), f"native 1:1 -- every panel is the full {nx}x{ny} source, no resampling "
               f"({px:.5f} um/px)", font=f_s, fill=MUTED, anchor="lt")
        lut_key(d, L_M, 214, colour, VMAX, f_l)

        bar = 50.0 / px
        by = T_M + 3 * (ny + GAP) - GAP + 46
        d.rectangle([L_M, by, L_M + bar, by + 9], fill=INK)
        d.text((L_M, by + 20), "50 um", font=f_s, fill=INK, anchor="lt")

        out = RUN / f"{sph}_depth_montage_{tag}_native.png"
        img.save(out, compress_level=6)
        print(f"  wrote {out.name}   {W}x{H} px   display 0-{VMAX:.0f}"
              f"   (p99.9 pre {np.percentile(A,99.9):.0f} / post {np.percentile(B,99.9):.0f})")
