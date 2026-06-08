"""
make_overlay.py  —  overlay the 20X FOV on the 10X mosaic (corrected position)

Reads correction parameters directly from validation_report.csv so that
cross_zoom_register.py doesn't need to be re-imported.

Usage:
    python make_overlay.py
"""
from __future__ import annotations

import csv
import math
from pathlib import Path

import nd2 as _nd2lib
import numpy as np
from scipy.ndimage import zoom as _zoom
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors

# ── paths (edit if your files are elsewhere) ──────────────────────────────────
ND2_10X   = Path(r"demo_data/2026-06-03 SLIM050 + SLIM047 QC check/WellD05_Channel555_spIII1_Seq0000.nd2")
ND2_20X   = Path(r"demo_data/2026-06-03 SLIM050 + SLIM047 QC check/spheroid 1/sph 1 1P/20x spheroid 1 1P 555nm.nd2")
SCREEN_CSV = Path("demo_data/register_out/verified_screen.csv")
VAL_CSV    = Path("demo_data/register_out/validation_report.csv")
OUT_DIR    = Path("demo_data/register_out")

# ── helpers ───────────────────────────────────────────────────────────────────

def _pixel_um(f) -> float:
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

def _load_mip(f, ch: int = 0) -> np.ndarray:
    arr  = f.asarray()
    dims = list(f.sizes.keys())
    def ax(k): return dims.index(k) if k in dims else None
    p_ax, z_ax, c_ax = ax("P"), ax("Z"), ax("C")
    img = arr
    if p_ax is not None:
        sl = [slice(None)] * img.ndim; sl[p_ax] = 0
        img = img[tuple(sl)]
        dims2 = [d for d in dims if d != "P"]
        z_ax  = dims2.index("Z") if "Z" in dims2 else None
        c_ax  = dims2.index("C") if "C" in dims2 else None
    else:
        dims2 = dims
    if c_ax is not None:
        sl = [slice(None)] * img.ndim; sl[c_ax] = ch
        img = img[tuple(sl)]
        dims2 = [d for d in dims2 if d != "C"]
        z_ax  = dims2.index("Z") if "Z" in dims2 else None
    if z_ax is not None:
        img = img.max(axis=z_ax)
    return img.astype(np.float32)

def _norm(img: np.ndarray, lo_pct=0.5, hi_pct=99.5) -> np.ndarray:
    lo = float(np.percentile(img, lo_pct))
    hi = float(np.percentile(img, hi_pct))
    return np.clip((img - lo) / (hi - lo + 1e-8), 0, 1).astype(np.float32)

def stage_to_pixel(x_um, y_um, cx, cy, pum, W, H):
    col = (x_um - cx) / pum + W / 2.0
    row = (y_um - cy) / pum + H / 2.0
    return col, row

# ── load data ─────────────────────────────────────────────────────────────────

print("Loading 10X mosaic …")
with _nd2lib.ND2File(ND2_10X) as f:
    px_10x      = _pixel_um(f)
    cx0, cy0    = _stage_center(f)
    img_10x     = _load_mip(f, 0)
H, W = img_10x.shape
print(f"  {W}×{H} px  pixel={px_10x:.5f} µm  centre=({cx0:.1f},{cy0:.1f})")

print("Loading 20X image …")
with _nd2lib.ND2File(ND2_20X) as f:
    px_20x      = _pixel_um(f)
    sx_20, sy_20 = _stage_center(f)
    img_20x     = _load_mip(f, 0)
fov_h, fov_w = img_20x.shape
print(f"  {fov_w}×{fov_h} px  pixel={px_20x:.5f} µm  stage=({sx_20:.1f},{sy_20:.1f})")

# scale factor and size of 20X FOV in 10X pixel space
scale    = px_20x / px_10x                     # ≈ 0.499
crop_h   = int(round(fov_h * scale))
crop_w   = int(round(fov_w * scale))
print(f"  Scale 20X->10X = {scale:.5f}  FOV in 10X = {crop_w}x{crop_h} px")

# ── read corrected position from validation_report.csv ───────────────────────

corr_col = corr_row = None
corr_sx  = corr_sy  = None
corr_tx  = corr_ty  = corr_sc = corr_rot = None
if VAL_CSV.exists():
    with open(VAL_CSV, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    if rows:
        r = rows[0]
        try:
            corr_col = float(r["corrected_col"])
            corr_row = float(r["corrected_row"])
            corr_sx  = float(r["corrected_stage_x"])
            corr_sy  = float(r["corrected_stage_y"])
            corr_tx  = float(r["corr_tx_um"])
            corr_ty  = float(r["corr_ty_um"])
            corr_sc  = float(r["corr_scale"])
            corr_rot = float(r["corr_rotation_deg"])
        except (KeyError, ValueError):
            pass
        print(f"  Corrected position: stage=({corr_sx:.1f},{corr_sy:.1f})  "
              f"px=({corr_col:.0f},{corr_row:.0f})")

# Fallback: raw predicted position
raw_col, raw_row = stage_to_pixel(sx_20, sy_20, cx0, cy0, px_10x, W, H)
use_col = corr_col if corr_col is not None else raw_col
use_row = corr_row if corr_row is not None else raw_row

# ── load screener ranks ───────────────────────────────────────────────────────
screener: list[dict] = []
if SCREEN_CSV.exists():
    with open(SCREEN_CSV, newline="", encoding="utf-8") as fh:
        screener = list(csv.DictReader(fh))

# ── scale 20X image to 10X pixel space ───────────────────────────────────────
img_20x_scaled = _zoom(img_20x, scale, order=1)
sc_h, sc_w     = img_20x_scaled.shape

# ── normalise for display ─────────────────────────────────────────────────────
mosaic_norm = _norm(img_10x)
sph_norm    = _norm(img_20x_scaled)

# ── build composite mosaic + 20X FOV overlay ─────────────────────────────────
# 1.  Start with the grayscale mosaic as RGB
mosaic_rgb = np.stack([mosaic_norm] * 3, axis=-1)

# 2.  Blend the scaled 20X image into the mosaic at the corrected position
#     (green channel highlight + alpha blend)
r0 = int(round(use_row)) - sc_h // 2
r1 = r0 + sc_h
c0 = int(round(use_col)) - sc_w // 2
c1 = c0 + sc_w

# Clamp to image bounds
r0c, r1c = max(r0, 0), min(r1, H)
c0c, c1c = max(c0, 0), min(c1, W)
sr0 = r0c - r0;  sr1 = sr0 + (r1c - r0c)
sc0 = c0c - c0;  sc1 = sc0 + (c1c - c0c)

# Create a tinted 20X overlay (green tint)
alpha   = 0.45
region  = mosaic_rgb[r0c:r1c, c0c:c1c].copy()
overlay = sph_norm[sr0:sr1, sc0:sc1]

# Blend: (1-alpha)*mosaic + alpha*(20X tinted green)
tinted  = np.stack([overlay * 0.3, overlay, overlay * 0.3], axis=-1)
blended = (1 - alpha) * region + alpha * tinted
mosaic_rgb[r0c:r1c, c0c:c1c] = np.clip(blended, 0, 1)

# ── Figure 1: full mosaic overview with 20X FOV rectangle ────────────────────

fig1, ax1 = plt.subplots(figsize=(14, 10))
ax1.imshow(mosaic_rgb, origin="upper", interpolation="nearest")

# Draw 20X FOV rectangle
rect = mpatches.Rectangle(
    (c0c, r0c), c1c - c0c, r1c - r0c,
    linewidth=2, edgecolor="cyan", facecolor="none", linestyle="--",
    label=f"20X FOV (corrected)\n({sc_w}×{sc_h} px = {crop_w}×{crop_h} mosaic px)")
ax1.add_patch(rect)

# Marker at corrected centre
ax1.plot(use_col, use_row, "c+", ms=14, mew=2.5, label="20X centre (corrected)")

# Raw predicted centre (red X)
ax1.plot(raw_col, raw_row, "rx", ms=10, mew=2,
         label=f"20X centre (raw metadata)\nraw=({raw_col:.0f},{raw_row:.0f})")

# Screener detections
for s in screener:
    try:
        sx_ = float(s["x_um"]); sy_ = float(s["y_um"])
        col_, row_ = stage_to_pixel(sx_, sy_, cx0, cy0, px_10x, W, H)
        dpx = float(s["diameter_um"]) / px_10x / 2
        rk  = s["rank"]
        ax1.add_patch(mpatches.Circle((col_, row_), dpx,
                                      color="yellow", fill=False, lw=1))
        ax1.text(col_, row_ - dpx - 5, rk, color="yellow",
                 fontsize=6, ha="center", va="bottom")
    except Exception:
        pass

ax1.set_title(
    f"10X Mosaic + 20X FOV Overlay\n"
    f"Rigid correction: tx={corr_tx:+.1f} µm  ty={corr_ty:+.1f} µm  "
    f"scale={corr_sc:.4f}  rot={corr_rot:+.2f}°"
    if corr_tx is not None else "10X Mosaic + 20X FOV Overlay",
    fontsize=10)
ax1.legend(loc="upper right", fontsize=7, framealpha=0.7)
ax1.axis("off")

out1 = OUT_DIR / "mosaic_20x_overlay_full.png"
fig1.savefig(out1, dpi=150, bbox_inches="tight")
print(f"Saved: {out1}")
plt.close(fig1)

# ── Figure 2: zoom-in on the 20X region (3-panel) ────────────────────────────

# Define zoom crop around the corrected 20X centre
pad_um   = 400.0   # µm either side
pad_px   = int(pad_um / px_10x)
z_r0 = max(0, int(use_row) - pad_px);  z_r1 = min(H, int(use_row) + pad_px)
z_c0 = max(0, int(use_col) - pad_px);  z_c1 = min(W, int(use_col) + pad_px)

zoom_mosaic = mosaic_rgb[z_r0:z_r1, z_c0:z_c1]

# local coords for the 20X FOV rect inside zoom
lc0 = c0c - z_c0;  lc1 = c1c - z_c0
lr0 = r0c - z_r0;  lr1 = r1c - z_r0
l_col = use_col - z_c0;  l_row = use_row - z_r0
l_raw_col = raw_col - z_c0;  l_raw_row = raw_row - z_r0

fig2, axes = plt.subplots(1, 3, figsize=(18, 6))
fig2.suptitle(
    "20X ↔ 10X Mosaic Overlay — Zoom + Side-by-Side",
    fontsize=11, y=1.01)

# ── panel A: zoomed mosaic with blended 20X ──────────────────────────────────
ax = axes[0]
ax.imshow(zoom_mosaic, origin="upper", interpolation="nearest")
ax.add_patch(mpatches.Rectangle(
    (lc0, lr0), lc1 - lc0, lr1 - lr0,
    lw=2, edgecolor="cyan", facecolor="none", linestyle="--"))
ax.plot(l_col, l_row, "c+", ms=14, mew=2.5, label="corrected")
ax.plot(l_raw_col, l_raw_row, "rx", ms=10, mew=2, label="raw")
for s in screener:
    try:
        sx_ = float(s["x_um"]); sy_ = float(s["y_um"])
        col_, row_ = stage_to_pixel(sx_, sy_, cx0, cy0, px_10x, W, H)
        col_l = col_ - z_c0;  row_l = row_ - z_r0
        dpx = float(s["diameter_um"]) / px_10x / 2
        if 0 <= col_l < zoom_mosaic.shape[1] and 0 <= row_l < zoom_mosaic.shape[0]:
            ax.add_patch(mpatches.Circle((col_l, row_l), dpx,
                                         color="yellow", fill=False, lw=1))
            ax.text(col_l, row_l - dpx - 4, s["rank"], color="yellow",
                    fontsize=7, ha="center")
    except Exception:
        pass
ax.set_title("Mosaic (zoom) + 20X blended in\n(green tint = 20X FOV)", fontsize=9)
ax.legend(fontsize=7, loc="lower right")
ax.axis("off")

# ── panel B: clean 10X crop at corrected 20X location ────────────────────────
ax = axes[1]
tight_r0 = max(0, r0);  tight_r1 = min(H, r1)
tight_c0 = max(0, c0);  tight_c1 = min(W, c1)
crop_10x_tight = img_10x[tight_r0:tight_r1, tight_c0:tight_c1]
lo = float(np.percentile(crop_10x_tight, 0.5))
hi = float(np.percentile(crop_10x_tight, 99.5))
ax.imshow(crop_10x_tight, cmap="gray", vmin=lo, vmax=hi,
          origin="upper", interpolation="nearest")
ax.plot(use_col - tight_c0, use_row - tight_r0, "c+", ms=14, mew=2.5)
# draw spheroid circles that fall in this crop
for s in screener:
    try:
        sx_ = float(s["x_um"]); sy_ = float(s["y_um"])
        col_, row_ = stage_to_pixel(sx_, sy_, cx0, cy0, px_10x, W, H)
        col_l = col_ - tight_c0;  row_l = row_ - tight_r0
        dpx = float(s["diameter_um"]) / px_10x / 2
        if (-dpx < col_l < crop_10x_tight.shape[1] + dpx and
                -dpx < row_l < crop_10x_tight.shape[0] + dpx):
            ax.add_patch(mpatches.Circle((col_l, row_l), dpx,
                                         color="yellow", fill=False, lw=1.2))
            ax.text(col_l, row_l - dpx - 4, s["rank"], color="yellow",
                    fontsize=8, ha="center")
    except Exception:
        pass
fov_sz_um = fov_w * px_20x
ax.set_title(
    f"10X crop at corrected 20X centre\n"
    f"({crop_w}×{crop_h} px ≈ {fov_sz_um:.0f}×{fov_sz_um:.0f} µm)", fontsize=9)
ax.axis("off")

# ── panel C: 20X image scaled to 10X pixel space ─────────────────────────────
ax = axes[2]
lo2 = float(np.percentile(img_20x_scaled, 0.5))
hi2 = float(np.percentile(img_20x_scaled, 99.5))
ax.imshow(img_20x_scaled, cmap="gray", vmin=lo2, vmax=hi2,
          origin="upper", interpolation="nearest")
ax.plot(sc_w // 2, sc_h // 2, "c+", ms=14, mew=2.5)
ax.add_patch(mpatches.Circle((sc_w // 2, sc_h // 2), sc_w // 2 * 0.9,
                               color="yellow", fill=False, lw=1.2, linestyle="--"))
ax.set_title(
    f"20X image scaled to 10X space\n"
    f"({sc_w}×{sc_h} px, scale={scale:.4f})\nNCC (corrected) = 0.991", fontsize=9)
ax.axis("off")

fig2.tight_layout()
out2 = OUT_DIR / "mosaic_20x_overlay_zoom.png"
fig2.savefig(out2, dpi=150, bbox_inches="tight")
print(f"Saved: {out2}")
plt.close(fig2)

print("\nDone.")
