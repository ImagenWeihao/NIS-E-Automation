r"""Pre-PA | Pre-PA + PA zone | Post-PA, one figure per spheroid, all channels.

    python pa_fig_prepost.py <run_dir> [pa_glob]

Works for single-plane and z-stack runs alike (a stack collapses to the plane nearest the
PA depth, which is stated on the figure).

CONTRAST: one fixed display range per ROW, computed once over the pre AND post panels
together and applied to both. Never a per-panel percentile stretch -- that renormalises
each image to its own noise floor and flatters whichever side is dimmer. The range is
printed on every row so the level is never implied.

The PA square is drawn from the PA series' OWN geometry and position -- side from
sizes['X'] * voxel_size().x, centre from each PA file's stagePositionUm matched to this
spheroid's readout position. It is never taken from a GUI field or config file: on run
0807_2 the pa/*.nd2 turned out to be byte copies of the readout stacks, and the square
used in that analysis had come from a config value rather than the data. If no PA file
matches this spheroid, the zone panel says so instead of drawing an assumed box.
"""
import glob
import sys
from pathlib import Path

import numpy as np
import nd2
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap

RUN = Path(sys.argv[1])
PA_GLOB = sys.argv[2] if len(sys.argv) > 2 else "*ZStack*.nd2"
BG, INK, MUTED, CYAN, RED = "#1e1e2e", "#e8e8f0", "#9aa0b5", "#22e5e5", "#f38ba8"
ROWS = [("890nm_mBeRFP", "890 nm  mBeRFP\nT-cell identity", "#ff2020"),
        ("940nm_PAsfGFP", "940 nm  PAsfGFP\nPA readout", "#20ff20"),
        ("1050nm_depth", "1050 nm\nspheroid structure", "#ff2020")]


def lut(hi):
    return LinearSegmentedColormap.from_list("l", ["#000000", hi])


def nice_max(v):
    for s in (25, 50, 100, 200, 400, 800, 1200, 2000, 4000, 8000, 16000, 65535):
        if v <= s:
            return s
    return 65535


def load(phase, tag, sph):
    fs = glob.glob(str(RUN / "nd2" / f"{phase}_{tag}" / f"{sph}_zstack.nd2"))
    if not fs:
        return None
    with nd2.ND2File(fs[0]) as f:
        a = f.asarray().astype(np.float32)
        px = f.voxel_size().x
        sp = f.frame_metadata(0).channels[0].position.stagePositionUm
        lp = f.unstructured_metadata()["ImageMetadataLV"]["SLxExperiment"]["LoopPars"]
        z = float(lp["ZLow"]) + float(lp["ZStep"]) * np.arange(a.shape[0] if a.ndim == 3 else 1)
    if a.ndim == 2:
        a = a[None]
    return a, px, z, (sp.x, sp.y)


# ── PA stimulus geometry, from the PA files themselves ─────────────────────
pa = []
for p in sorted(glob.glob(str(RUN / "pa" / PA_GLOB))):
    with nd2.ND2File(p) as f:
        sp = f.frame_metadata(0).channels[0].position.stagePositionUm
        pa.append((Path(p).name, sp.x, sp.y, sp.z, f.sizes["X"] * f.voxel_size().x))
print(f"PA files: {len(pa)}")
for n, x, y, z, s in pa:
    print(f"   {n[:56]:56s} ({x:.1f}, {y:.1f}) z={z:.1f}  {s:.1f} um sq")

spheroids = sorted({Path(p).stem.replace("_zstack", "")
                    for p in glob.glob(str(RUN / "nd2" / "prePA_*" / "*_zstack.nd2"))})
print(f"\nspheroids with a pre-PA set: {spheroids}")

for sph in spheroids:
    got = {t: (load("prePA", t, sph), load("postPA", t, sph)) for t, _, _ in ROWS}
    got = {t: v for t, v in got.items() if v[0] and v[1]}
    if not got:
        print(f"  {sph}: no complete pre+post set, skipped")
        continue
    (_, px, _, xy) = list(got.values())[0][0]

    # Match this spheroid to its own PA exposure by stage position.
    mine = min(pa, key=lambda r: (r[1] - xy[0]) ** 2 + (r[2] - xy[1]) ** 2) if pa else None
    if mine and ((mine[1] - xy[0]) ** 2 + (mine[2] - xy[1]) ** 2) ** 0.5 > 50:
        mine = None
    side = mine[4] if mine else None

    rows = [(t, lab, col) for t, lab, col in ROWS if t in got]
    fig, axes = plt.subplots(len(rows), 3, figsize=(11.6, 4.0 * len(rows)), facecolor=BG,
                             squeeze=False)
    fig.subplots_adjust(left=0.135, right=0.995, top=1 - 0.30 / len(rows),
                        bottom=0.045, wspace=0.02, hspace=0.02)

    zshow = None
    for r, (tag, label, colour) in enumerate(rows):
        (pre, px, zs, _), (post, _, _, _) = got[tag]
        k = 0 if pre.shape[0] == 1 else int(np.argmin(np.abs(zs - (mine[3] if mine else zs.mean()))))
        zshow = zs[k]
        a_pre, a_post = pre[k], post[k]
        vmax = nice_max(max(np.percentile(a_pre, 99.9), np.percentile(a_post, 99.9)))
        h, w = a_pre.shape
        for c in range(3):
            ax = axes[r][c]
            ax.set_xticks([]); ax.set_yticks([]); ax.set_facecolor("black")
            for s_ in ax.spines.values():
                s_.set_visible(False)
            ax.imshow(a_post if c == 2 else a_pre, cmap=lut(colour), vmin=0, vmax=vmax,
                      interpolation="nearest")
            # Square on the zone column AND the post-PA column. Without it on post-PA the
            # reader has to carry the position across from the middle panel by eye to judge
            # whether a change sits inside the stimulated volume -- the whole question.
            if c in (1, 2):
                if side:
                    S = side / px
                    ax.add_patch(Rectangle((w / 2 - S / 2, h / 2 - S / 2), S, S,
                                           fill=False, edgecolor=CYAN, lw=1.8))
                elif c == 1:
                    ax.text(0.5, 0.06, "no PA file matches this spheroid",
                            transform=ax.transAxes, color=RED, fontsize=9,
                            ha="center", va="bottom")
            if c == 0:
                ax.set_ylabel(f"{label}\ndisplay 0-{vmax}", color=INK, fontsize=10, labelpad=8)
            if r == 0:
                ax.set_title(["Pre-PA",
                              f"Pre-PA + PA zone ({side:.0f} um sq)" if side
                              else "Pre-PA + PA zone (unknown)",
                              "Post-PA + PA zone" if side else "Post-PA"][c],
                             color=INK, fontsize=12, pad=8)
        if r == len(rows) - 1:
            ax = axes[r][0]
            ax.plot([12, 12 + 50.0 / px], [h - 18] * 2, color=INK, lw=3, solid_capstyle="butt")
            ax.text(12, h - 26, "50 um", color=INK, fontsize=10, va="bottom")

    pa_txt = (f"850 nm PA at ({mine[1]:.0f}, {mine[2]:.0f}) um, z={mine[3]:.0f} um, "
              f"{side:.0f} um square" if mine else "PA exposure not identifiable in pa/")
    fig.suptitle(f"{sph}   pre vs post photoactivation   readout z = {zshow:.0f} um\n"
                 f"{pa_txt}\n"
                 f"each row: ONE fixed display range, identical for pre and post",
                 color=INK, fontsize=11.5, y=0.995, va="top")
    out = RUN / f"{sph}_PA_prepost.png"
    fig.savefig(out, dpi=150, facecolor=BG); plt.close(fig)
    print(f"  wrote {out.name}")
