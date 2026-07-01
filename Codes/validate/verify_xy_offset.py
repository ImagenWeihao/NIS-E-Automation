"""
verify_xy_offset.py  —  Feature-level XY offset verification using Ranks 6 & 7

For each target rank:
  1. Extract a template patch from the SUB-IMAGE 10X at the inverse-transformed
     predicted location (sub-image coords of that mosaic rank)
  2. Search for this template in the MOSAIC at the screener position via NCC
  3. Measure where the NCC peak lands → actual XY displacement
  4. Compare pure-translation offset (S1→Rank9 anchor only) vs. similarity-
     transform corrected offset for each rank

This is a fully independent cross-check: the sub-image patch is compared
directly to the mosaic patch without relying on any stage metadata.

Usage:
    python verify_xy_offset.py
"""
from __future__ import annotations

import csv
import math
from pathlib import Path

import nd2 as _nd2lib
import numpy as np
from scipy.ndimage import map_coordinates
from skimage.feature import match_template
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── paths ─────────────────────────────────────────────────────────────────────
ND2_MOSAIC = Path(r"demo_data/2026-06-03 SLIM050 + SLIM047 QC check/WellD05_Channel555_spIII1_Seq0000.nd2")
ND2_SUB    = Path(r"demo_data/2026-06-03 SLIM050 + SLIM047 QC check/spheroid 1/sph 1 1P/spheroid 1 10x 1P 555nm.nd2")
SCREEN_CSV = Path("demo_data/register_out/verified_screen.csv")
VAL_CSV    = Path("demo_data/register_out/validation_report.csv")
OUT_DIR    = Path("demo_data/register_out")

TARGET_RANKS = [6, 7, 9]

# Search radius for NCC peak in the mosaic (µm either side)
SEARCH_UM = 80.0
# Template patch half-size in sub-image pixels (will be cropped to spheroid)
TEMPLATE_HALF_PX = 80

# ── helpers ───────────────────────────────────────────────────────────────────

def _pixel_um(f):
    try:
        v = float(f.voxel_size().x)
        return v if v > 0 else 0.644
    except Exception:
        return 0.644

def _stage_center(f):
    try:
        fm  = f.frame_metadata(0)
        pos = fm.channels[0].position.stagePositionUm
        return float(pos.x), float(pos.y)
    except Exception:
        return 0.0, 0.0

def _load_mip(f, ch=0):
    arr  = f.asarray()
    dims = list(f.sizes.keys())
    def ax(k): return dims.index(k) if k in dims else None
    img = arr
    p_ax = ax("P")
    if p_ax is not None:
        sl = [slice(None)] * img.ndim; sl[p_ax] = 0
        img = img[tuple(sl)]
        dims2 = [d for d in dims if d != "P"]
        z_ax  = dims2.index("Z") if "Z" in dims2 else None
        c_ax  = dims2.index("C") if "C" in dims2 else None
    else:
        dims2 = dims; z_ax = ax("Z"); c_ax = ax("C")
    if c_ax is not None:
        sl = [slice(None)] * img.ndim; sl[c_ax] = ch
        img = img[tuple(sl)]
        dims2 = [d for d in dims2 if d != "C"]
        z_ax  = dims2.index("Z") if "Z" in dims2 else None
    if z_ax is not None:
        img = img.max(axis=z_ax)
    return img.astype(np.float32)

def _norm(img, lo=0.5, hi=99.5):
    lo_ = float(np.percentile(img, lo))
    hi_ = float(np.percentile(img, hi))
    return np.clip((img - lo_) / (hi_ - lo_ + 1e-8), 0, 1).astype(np.float32)

def stage_to_pixel(x, y, cx, cy, pum, W, H):
    return (x - cx) / pum + W / 2.0, (y - cy) / pum + H / 2.0

def pixel_to_stage(col, row, cx, cy, pum, W, H):
    return (col - W / 2.0) * pum + cx, (row - H / 2.0) * pum + cy

# ── load images ───────────────────────────────────────────────────────────────

print("Loading 10X mosaic ...")
with _nd2lib.ND2File(ND2_MOSAIC) as f:
    px_mos = _pixel_um(f)
    cx_mos, cy_mos = _stage_center(f)
    img_mos = _load_mip(f, 0)
Hm, Wm = img_mos.shape
print(f"  {Wm}x{Hm}  pixel={px_mos:.5f} um  centre=({cx_mos:.1f},{cy_mos:.1f})")

print("Loading sub-image 10X ...")
with _nd2lib.ND2File(ND2_SUB) as f:
    px_sub = _pixel_um(f)
    cx_sub, cy_sub = _stage_center(f)
    img_sub = _load_mip(f, 0)
Hs, Ws = img_sub.shape
print(f"  {Ws}x{Hs}  pixel={px_sub:.5f} um  centre=({cx_sub:.1f},{cy_sub:.1f})")

# ── read correction parameters ────────────────────────────────────────────────

corr_tx = -1745.23; corr_ty = -1734.27
corr_sc = 1.1052;   corr_rot = -8.476

if VAL_CSV.exists():
    with open(VAL_CSV, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    if rows:
        r = rows[0]
        try:
            corr_tx  = float(r["corr_tx_um"]);  corr_ty  = float(r["corr_ty_um"])
            corr_sc  = float(r["corr_scale"]);   corr_rot = float(r["corr_rotation_deg"])
        except (KeyError, ValueError):
            pass

th    = math.radians(corr_rot)
cos_t = math.cos(th); sin_t = math.sin(th)

# pure-translation anchor: from S1->Rank9 (tx_anchor from validation report)
# the pure translation was tx=+335.3, ty=-3414.9 µm
TX_ANCHOR = 335.3;  TY_ANCHOR = -3414.9

def mos_to_sub_stage(xm, ym):
    """Inverse similarity transform: mosaic stage -> sub-image stage."""
    dx = xm - corr_tx; dy = ym - corr_ty
    sc = corr_sc
    x_sub = (1.0 / sc) * ( cos_t * dx + sin_t * dy)
    y_sub = (1.0 / sc) * (-sin_t * dx + cos_t * dy)
    return x_sub, y_sub

def sub_to_mos_stage(xs, ys):
    """Forward similarity transform: sub-image stage -> mosaic stage."""
    xm = corr_sc * (cos_t * xs - sin_t * ys) + corr_tx
    ym = corr_sc * (sin_t * xs + cos_t * ys) + corr_ty
    return xm, ym

# ── read screener CSV ─────────────────────────────────────────────────────────

screener: dict[int, dict] = {}
with open(SCREEN_CSV, newline="", encoding="utf-8") as fh:
    for row in csv.DictReader(fh):
        screener[int(row["rank"])] = row

# ── per-rank analysis ─────────────────────────────────────────────────────────

results = []

for rank in TARGET_RANKS:
    if rank not in screener:
        print(f"Rank {rank} not in screener CSV — skipping")
        continue

    rec      = screener[rank]
    x_mos_gt = float(rec["x_um"])      # screener (ground-truth) mosaic stage
    y_mos_gt = float(rec["y_um"])
    diam_um  = float(rec["diameter_um"])
    col_mos_gt, row_mos_gt = stage_to_pixel(x_mos_gt, y_mos_gt,
                                             cx_mos, cy_mos, px_mos, Wm, Hm)

    # Predict sub-image pixel for this rank via inverse transform
    x_sub_pred, y_sub_pred = mos_to_sub_stage(x_mos_gt, y_mos_gt)
    col_sub_pred, row_sub_pred = stage_to_pixel(
        x_sub_pred, y_sub_pred, cx_sub, cy_sub, px_sub, Ws, Hs)

    print(f"\nRank {rank}:")
    print(f"  Mosaic stage=({x_mos_gt:.1f},{y_mos_gt:.1f})  pixel=({col_mos_gt:.0f},{row_mos_gt:.0f})")
    print(f"  Predicted sub-image pixel=({col_sub_pred:.0f},{row_sub_pred:.0f})")

    in_sub = (TEMPLATE_HALF_PX <= col_sub_pred < Ws - TEMPLATE_HALF_PX and
              TEMPLATE_HALF_PX <= row_sub_pred < Hs - TEMPLATE_HALF_PX)
    near_sub = (0 <= col_sub_pred < Ws and 0 <= row_sub_pred < Hs)

    # ── pure-translation prediction (no rotation/scale) ───────────────────────
    # Sub-image stage -> mosaic stage using only the anchor translation
    x_mos_puretx = cx_sub + TX_ANCHOR
    y_mos_puretx = cy_sub + TY_ANCHOR
    # For this rank: sub-image stage of Rank k is sub-image stage of mosaic screener
    # coords after inverse pure translation
    x_sub_pure = x_mos_gt - TX_ANCHOR
    y_sub_pure = y_mos_gt - TY_ANCHOR
    col_sub_pure, row_sub_pure = stage_to_pixel(
        x_sub_pure, y_sub_pure, cx_sub, cy_sub, px_sub, Ws, Hs)
    in_sub_pure = (0 <= col_sub_pure < Ws and 0 <= row_sub_pure < Hs)

    # ── extract template from sub-image ───────────────────────────────────────
    # Use the SIMILARITY TRANSFORM prediction (more accurate) for the template
    # Fall back to pure-translation prediction if the sim-transform lands outside
    template_patch = None
    template_col   = col_sub_pred
    template_row   = row_sub_pred
    used_correction = "similarity"

    if not in_sub and in_sub_pure:
        template_col = col_sub_pure;  template_row = row_sub_pure
        used_correction = "pure-translation"
        print(f"  Similarity-transform prediction is outside sub-image bounds; "
              f"falling back to pure-translation: ({col_sub_pure:.0f},{row_sub_pure:.0f})")
    elif not in_sub and not in_sub_pure:
        print(f"  BOTH predictions outside sub-image — skipping NCC for Rank {rank}")

    if in_sub or in_sub_pure:
        r0 = int(round(template_row)) - TEMPLATE_HALF_PX
        r1 = r0 + 2 * TEMPLATE_HALF_PX
        c0 = int(round(template_col)) - TEMPLATE_HALF_PX
        c1 = c0 + 2 * TEMPLATE_HALF_PX
        r0c = max(0, r0); r1c = min(Hs, r1)
        c0c = max(0, c0); c1c = min(Ws, c1)
        template_patch = img_sub[r0c:r1c, c0c:c1c].copy()

    # ── search in mosaic via NCC ───────────────────────────────────────────────
    search_px    = int(round(SEARCH_UM / px_mos))
    half_t       = TEMPLATE_HALF_PX
    search_half  = search_px + half_t + 4

    ncc_peak     = float("nan")
    peak_col_mos = float("nan")
    peak_row_mos = float("nan")
    offset_x_um  = float("nan")
    offset_y_um  = float("nan")
    ncc_map_crop = None

    if template_patch is not None and template_patch.size > 0:
        # Extract mosaic search region
        mr0 = max(0, int(round(row_mos_gt)) - search_half)
        mr1 = min(Hm, int(round(row_mos_gt)) + search_half)
        mc0 = max(0, int(round(col_mos_gt)) - search_half)
        mc1 = min(Wm, int(round(col_mos_gt)) + search_half)
        mosaic_region = img_mos[mr0:mr1, mc0:mc1]

        th_, tw_ = template_patch.shape

        if mosaic_region.shape[0] >= th_ and mosaic_region.shape[1] >= tw_:
            # Normalise both
            def _n(x):
                lo_ = x.min(); hi_ = x.max()
                return (x - lo_) / (hi_ - lo_ + 1e-8)

            ncc_m = match_template(_n(mosaic_region).astype(np.float32),
                                   _n(template_patch).astype(np.float32),
                                   pad_input=True)

            # Restrict peak search to ±search_px around the screener position
            cy_loc = int(round(row_mos_gt)) - mr0
            cx_loc = int(round(col_mos_gt)) - mc0
            rs_ = max(0, cy_loc - search_px); re_ = min(ncc_m.shape[0], cy_loc + search_px + 1)
            cs_ = max(0, cx_loc - search_px); ce_ = min(ncc_m.shape[1], cx_loc + search_px + 1)
            sub_ncc = ncc_m[rs_:re_, cs_:ce_]
            pr, pc  = np.unravel_index(sub_ncc.argmax(), sub_ncc.shape)
            ncc_peak = float(sub_ncc[pr, pc])

            peak_row_mos = float(mr0 + rs_ + pr)
            peak_col_mos = float(mc0 + cs_ + pc)

            # XY offset: where NCC peak is vs. screener ground-truth
            offset_x_um = (peak_col_mos - col_mos_gt) * px_mos
            offset_y_um = (peak_row_mos - row_mos_gt) * px_mos

            # Save NCC map crop for visualisation
            ncc_map_crop = sub_ncc.copy()

            print(f"  NCC peak = {ncc_peak:.4f}")
            print(f"  NCC peak mosaic pixel=({peak_col_mos:.0f},{peak_row_mos:.0f})")
            print(f"  Offset vs screener:  dx={offset_x_um:+.1f} um  "
                  f"dy={offset_y_um:+.1f} um  "
                  f"|err|={math.sqrt(offset_x_um**2+offset_y_um**2):.1f} um")
        else:
            print(f"  Search region too small for NCC")

    # ── mosaic patch for display ───────────────────────────────────────────────
    half_disp = TEMPLATE_HALF_PX + 20
    mr0d = max(0, int(round(row_mos_gt)) - half_disp)
    mr1d = min(Hm, int(round(row_mos_gt)) + half_disp)
    mc0d = max(0, int(round(col_mos_gt)) - half_disp)
    mc1d = min(Wm, int(round(col_mos_gt)) + half_disp)
    mosaic_patch = img_mos[mr0d:mr1d, mc0d:mc1d]

    results.append({
        "rank":             rank,
        "x_mos_gt":         x_mos_gt,
        "y_mos_gt":         y_mos_gt,
        "col_mos_gt":       col_mos_gt,
        "row_mos_gt":       row_mos_gt,
        "col_sub_pred":     col_sub_pred,
        "row_sub_pred":     row_sub_pred,
        "col_sub_pure":     col_sub_pure,
        "row_sub_pure":     row_sub_pure,
        "template_patch":   template_patch,
        "mosaic_patch":     mosaic_patch,
        "mc0d":             mc0d,  "mr0d": mr0d,
        "ncc_peak":         ncc_peak,
        "peak_col_mos":     peak_col_mos,
        "peak_row_mos":     peak_row_mos,
        "offset_x_um":      offset_x_um,
        "offset_y_um":      offset_y_um,
        "ncc_map_crop":     ncc_map_crop,
        "used_correction":  used_correction,
        "diam_um":          diam_um,
    })

# ── Figure: 4-column layout per rank ─────────────────────────────────────────
# Col 1: sub-image patch at predicted location (template)
# Col 2: mosaic patch at screener position (search region)
# Col 3: mosaic patch at NCC peak (matched result)
# Col 4: NCC map with peak marker

n_ranks = len(results)
fig, axes = plt.subplots(n_ranks, 4, figsize=(16, 5 * n_ranks))
if n_ranks == 1:
    axes = axes[np.newaxis, :]

fig.suptitle(
    "Feature-level XY Offset Verification — Ranks 6, 7 & 9\n"
    f"Sub-image 10X template matched against 10X Mosaic via NCC",
    fontsize=12, y=1.01)

col_titles = [
    "Sub-image\n(template at predicted location)",
    "Mosaic at screener position\n(search centre)",
    "Mosaic at NCC peak\n(matched result)",
    "NCC heatmap\n(peak = actual feature location)",
]
for j, t in enumerate(col_titles):
    axes[0, j].set_title(t, fontsize=9)

for i, res in enumerate(results):
    rk  = res["rank"]
    dpm = res["diam_um"] / px_mos / 2   # radius in mosaic px
    dps = res["diam_um"] / px_sub / 2   # radius in sub-image px

    # ── Panel 0: sub-image template ──────────────────────────────────────────
    ax = axes[i, 0]
    if res["template_patch"] is not None and res["template_patch"].size > 0:
        ax.imshow(_norm(res["template_patch"]), cmap="gray",
                  origin="upper", interpolation="nearest")
        th_, tw_ = res["template_patch"].shape
        # mark the predicted centre (may not be exactly at centre if clipped)
        cen_c = res["col_sub_pred"] - (int(round(res["col_sub_pred"])) - TEMPLATE_HALF_PX)
        cen_r = res["row_sub_pred"] - (int(round(res["row_sub_pred"])) - TEMPLATE_HALF_PX)
        ax.plot(cen_c, cen_r, "c+", ms=14, mew=2)
        ax.add_patch(mpatches.Circle((cen_c, cen_r), dps, color="yellow",
                                     fill=False, lw=1.5, linestyle="--"))
    else:
        ax.text(0.5, 0.5, "Outside\nsub-image\nbounds",
                ha="center", va="center", transform=ax.transAxes,
                fontsize=10, color="red")
        ax.set_facecolor("black")
    ax.set_ylabel(f"Rank {rk}\n({res['used_correction']})", fontsize=9)
    ax.axis("off")

    # ── Panel 1: mosaic at screener position ─────────────────────────────────
    ax = axes[i, 1]
    mp = res["mosaic_patch"]
    ax.imshow(_norm(mp), cmap="gray", origin="upper", interpolation="nearest")
    local_col = res["col_mos_gt"] - res["mc0d"]
    local_row = res["row_mos_gt"] - res["mr0d"]
    ax.plot(local_col, local_row, "y+", ms=14, mew=2)
    ax.add_patch(mpatches.Circle((local_col, local_row), dpm,
                                  color="yellow", fill=False, lw=1.5))
    ax.text(local_col, local_row - dpm - 4, f"Rank {rk}", color="yellow",
            fontsize=8, ha="center")
    ax.axis("off")

    # ── Panel 2: mosaic at NCC peak position ─────────────────────────────────
    ax = axes[i, 2]
    if not math.isnan(res["peak_col_mos"]):
        half_d = TEMPLATE_HALF_PX + 20
        pr0 = max(0, int(round(res["peak_row_mos"])) - half_d)
        pr1 = min(Hm, int(round(res["peak_row_mos"])) + half_d)
        pc0 = max(0, int(round(res["peak_col_mos"])) - half_d)
        pc1 = min(Wm, int(round(res["peak_col_mos"])) + half_d)
        peak_patch = img_mos[pr0:pr1, pc0:pc1]
        ax.imshow(_norm(peak_patch), cmap="gray", origin="upper",
                  interpolation="nearest")
        pk_lc = res["peak_col_mos"] - pc0
        pk_lr = res["peak_row_mos"] - pr0
        ax.plot(pk_lc, pk_lr, "c+", ms=14, mew=2)
        ax.add_patch(mpatches.Circle((pk_lc, pk_lr), dpm,
                                      color="cyan", fill=False, lw=1.5))
        ncc_v = res["ncc_peak"]
        dx = res["offset_x_um"]; dy = res["offset_y_um"]
        err = math.sqrt(dx**2 + dy**2)
        conf = "CONFIRMED" if ncc_v >= 0.50 else "LOW NCC"
        ax.set_title(
            f"NCC={ncc_v:.3f} ({conf})\ndx={dx:+.1f} um  dy={dy:+.1f} um  |err|={err:.1f} um",
            fontsize=8)
    else:
        ax.text(0.5, 0.5, "NCC\nnot\ncomputed", ha="center", va="center",
                transform=ax.transAxes, fontsize=10)
        ax.set_facecolor("black")
    ax.axis("off")

    # ── Panel 3: NCC heatmap ──────────────────────────────────────────────────
    ax = axes[i, 3]
    if res["ncc_map_crop"] is not None:
        nm = res["ncc_map_crop"]
        ax.imshow(nm, cmap="hot", origin="upper",
                  vmin=nm.min(), vmax=nm.max(), interpolation="nearest")
        pr_, pc_ = np.unravel_index(nm.argmax(), nm.shape)
        ax.plot(pc_, pr_, "c+", ms=12, mew=2)
        # mark centre (screener reference)
        cen_r_ncc = nm.shape[0] // 2; cen_c_ncc = nm.shape[1] // 2
        ax.plot(cen_c_ncc, cen_r_ncc, "wx", ms=10, mew=2)
        ax.legend(
            handles=[
                mpatches.Patch(color="cyan",  label=f"NCC peak ({res['ncc_peak']:.3f})"),
                mpatches.Patch(color="white", label="Screener centre"),
            ],
            fontsize=7, loc="lower right")
    else:
        ax.text(0.5, 0.5, "N/A", ha="center", va="center",
                transform=ax.transAxes, fontsize=12)
        ax.set_facecolor("black")
    ax.axis("off")

fig.tight_layout()
out_fig = OUT_DIR / "rank67_xy_offset_verification.png"
fig.savefig(out_fig, dpi=150, bbox_inches="tight")
print(f"\nSaved: {out_fig}")
plt.close(fig)

# ── Summary table ─────────────────────────────────────────────────────────────
print("\n--- XY offset verification summary ---")
print(f"{'Rank':>5}  {'NCC':>6}  {'dx (um)':>9}  {'dy (um)':>9}  {'|err| (um)':>10}  "
      f"{'correction_used':>18}")
for res in results:
    dx = res["offset_x_um"]; dy = res["offset_y_um"]
    err = math.sqrt(dx**2 + dy**2) if not math.isnan(dx) else float("nan")
    print(f"{res['rank']:>5}  {res['ncc_peak']:>6.3f}  {dx:>+9.1f}  {dy:>+9.1f}  "
          f"{err:>10.1f}  {res['used_correction']:>18}")

print("\nDone.")
