r"""Depth-resolved PA efficiency for one run.

    python pa_efficiency.py <run_dir> <pa_glob> <sph_id>

WHY COUNTS AND NOT A post/pre RATIO
    On this sample the pre-PA 940 nm channel sits on the detector floor -- median 1
    count, 99.9th percentile 12. A post/pre ratio therefore divides by ~1 and explodes
    (peaks over 500x) without meaning anything: the same blob over a floor of 2 instead
    of 1 would read 275x. So the primary panel is mean counts per plane, with the
    spatial control carrying the comparison instead:

    inside(z)  = mean 940 counts inside the PA square
    outside(z) = mean 940 counts outside it, same plane
    contrast   = inside / outside, POST only. Both regions saw the same readout laser,
                 the same depth and the same bleaching, so a contrast above the pre-PA
                 baseline contrast is 850 nm effect and nothing else.

The second panel is the level-based read: voxels above 1200 on the 0-4000 LUT NIS-E
displays, per plane, pre and post. That is the operator's own true-signal rule applied
depth by depth.
"""
import sys
import glob
from pathlib import Path

import numpy as np
import nd2
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RUN, PA_GLOB, SPH = Path(sys.argv[1]), sys.argv[2], sys.argv[3]
BG, INK, MUTED, GRID = "#1e1e2e", "#e8e8f0", "#9aa0b5", "#313244"
CYAN, GREEN, RED, PINK = "#22e5e5", "#a6e3a1", "#f38ba8", "#cba6f7"
TRUE_SIGNAL = 1200.0     # NIS-E displays these channels on a 0-4000 LUT


def load(phase, tag):
    p = glob.glob(str(RUN / "nd2" / f"{phase}_{tag}" / "*.nd2"))[0]
    with nd2.ND2File(p) as f:
        a = f.asarray().astype(np.float32)
        px = f.voxel_size().x
        lp = f.unstructured_metadata()["ImageMetadataLV"]["SLxExperiment"]["LoopPars"]
        z = float(lp["ZLow"]) + float(lp["ZStep"]) * np.arange(a.shape[0])
    return a, px, z


pa = sorted(glob.glob(PA_GLOB))
with nd2.ND2File(pa[0]) as f:
    PA_SIDE = f.sizes["X"] * f.voxel_size().x
pz = sorted(nd2.ND2File(p).frame_metadata(0).channels[0].position.stagePositionUm.z
            for p in pa)
PA_LO, PA_HI = pz[0], pz[-1]

pre9, px, zs = load("prePA", "940nm_PAsfGFP")
post9, _, _ = load("postPA", "940nm_PAsfGFP")
nz, ny, nx = pre9.shape

S = PA_SIDE / px
Y, X = np.mgrid[:ny, :nx]
inside = np.maximum(np.abs(X - nx / 2), np.abs(Y - ny / 2)) <= S / 2
outside = ~inside

pre_in = np.array([pre9[i][inside].mean() for i in range(nz)])
pre_out = np.array([pre9[i][outside].mean() for i in range(nz)])
post_in = np.array([post9[i][inside].mean() for i in range(nz)])
post_out = np.array([post9[i][outside].mean() for i in range(nz)])

FLOOR = 0.2                                    # keeps the contrast finite at the floor
c_pre = pre_in / np.maximum(pre_out, FLOOR)
c_post = post_in / np.maximum(post_out, FLOOR)

hi_pre = np.array([(pre9[i] > TRUE_SIGNAL).sum() for i in range(nz)], float)
hi_post = np.array([(post9[i] > TRUE_SIGNAL).sum() for i in range(nz)], float)

band = (zs >= PA_LO) & (zs <= PA_HI)
peak = int(np.argmax(post_in))

print(f"{SPH}   PA {len(pa)} planes  z {PA_LO:.0f}-{PA_HI:.0f}  square {PA_SIDE:.0f} um")
print(f"  readout {nz} planes {zs[0]:.0f}-{zs[-1]:.0f} @ {zs[1]-zs[0]:.0f} um\n")
print(f"  940 mean counts INSIDE the square    pre {pre_in.mean():7.2f}   "
      f"post {post_in.mean():7.2f}   ({post_in.mean()/max(pre_in.mean(),1e-6):.1f}x)")
print(f"  940 mean counts OUTSIDE (control)    pre {pre_out.mean():7.2f}   "
      f"post {post_out.mean():7.2f}   ({post_out.mean()/max(pre_out.mean(),1e-6):.1f}x)")
print(f"  inside/outside contrast              pre {np.median(c_pre):7.2f}   "
      f"post {np.median(c_post):7.2f}   (median over planes)")
print(f"  peak post inside-square signal at z={zs[peak]:.0f} um "
      f"({post_in[peak]:.1f} counts), PA centre {(PA_LO+PA_HI)/2:.0f}")
print(f"  voxels >{TRUE_SIGNAL:.0f}/4000        pre {hi_pre.sum():.0f}   "
      f"post {hi_post.sum():.0f}   ({hi_post.sum()/max(hi_pre.sum(),1):.0f}x)")
inb = hi_post[band].sum()
print(f"     of the post voxels, {inb:.0f}/{hi_post.sum():.0f} "
      f"({100*inb/max(hi_post.sum(),1):.0f}%) lie inside the PA depth range")

fig, ax = plt.subplots(2, 1, figsize=(11, 9), facecolor=BG, sharex=True,
                       gridspec_kw={"height_ratios": [1.35, 1]})
fig.subplots_adjust(left=0.115, right=0.975, top=0.855, bottom=0.082, hspace=0.09)

for a_ in ax:
    a_.set_facecolor(BG)
    a_.grid(True, color=GRID, lw=0.7, which="both")
    a_.tick_params(colors=MUTED, labelsize=10)
    for s in a_.spines.values():
        s.set_color(GRID)
    a_.axvspan(PA_LO, PA_HI, color=CYAN, alpha=0.10, lw=0)
    a_.axvline(zs[peak], color=PINK, lw=1.2, ls=":")

ax[0].semilogy(zs, np.maximum(post_in, 0.05), color=GREEN, lw=2.4, marker="o", ms=4.5,
               label="post-PA, inside the PA square")
ax[0].semilogy(zs, np.maximum(post_out, 0.05), color=RED, lw=1.7, marker="o", ms=3.5,
               label="post-PA, outside (internal control)")
ax[0].semilogy(zs, np.maximum(pre_in, 0.05), color=GREEN, lw=1.4, ls="--", alpha=0.75,
               label="pre-PA, inside")
ax[0].semilogy(zs, np.maximum(pre_out, 0.05), color=RED, lw=1.2, ls="--", alpha=0.6,
               label="pre-PA, outside")
ax[0].set_ylabel("940 nm PAsfGFP\nmean counts per plane (log)", color=INK, fontsize=11)
lg = ax[0].legend(facecolor="#181825", edgecolor=GRID, fontsize=9.5, loc="upper left")
for t in lg.get_texts():
    t.set_color(INK)

ax[1].semilogy(zs, np.maximum(hi_post, 0.5), color=GREEN, lw=2.4, marker="o", ms=4.5,
               label="post-PA")
ax[1].semilogy(zs, np.maximum(hi_pre, 0.5), color=MUTED, lw=1.7, marker="o", ms=3.5,
               label="pre-PA")
ax[1].set_ylabel(f"940 nm voxels > {TRUE_SIGNAL:.0f} / 4000\nper plane (log)",
                 color=INK, fontsize=11)
ax[1].set_xlabel("stage z (um)", color=INK, fontsize=11)
lg = ax[1].legend(facecolor="#181825", edgecolor=GRID, fontsize=9.5, loc="upper left")
for t in lg.get_texts():
    t.set_color(INK)

fig.suptitle(
    f"{SPH}   depth-resolved photoactivation efficiency\n"
    f"850 nm PA, z {PA_LO:.0f}-{PA_HI:.0f} um ({len(pa)} planes, {PA_SIDE:.0f} um square) = "
    f"shaded band\nreadout 940 nm PAsfGFP, identical settings pre and post",
    color=INK, fontsize=11.5, y=0.985, va="top")
fig.text(0.115, 0.019,
         f"mean counts inside the square {pre_in.mean():.2f} -> {post_in.mean():.2f} "
         f"({post_in.mean()/max(pre_in.mean(),1e-6):.0f}x);  outside "
         f"{pre_out.mean():.2f} -> {post_out.mean():.2f} "
         f"({post_out.mean()/max(pre_out.mean(),1e-6):.1f}x).   "
         f"Voxels over {TRUE_SIGNAL:.0f}: {hi_pre.sum():.0f} -> {hi_post.sum():.0f}, "
         f"{100*inb/max(hi_post.sum(),1):.0f}% inside the PA depth range.",
         color=MUTED, fontsize=9.5, ha="left")

out = RUN / f"{SPH}_PA_efficiency_depth.png"
fig.savefig(out, dpi=160, facecolor=BG)
print(f"\nwrote {out}")
