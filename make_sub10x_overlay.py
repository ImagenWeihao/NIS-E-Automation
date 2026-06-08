"""
make_sub10x_overlay.py  —  overlay the single-position 10X sub-image on the 10X mosaic

Uses the rigid correction parameters (tx, ty, scale, rotation) from
validation_report.csv to warp the sub-image into the mosaic coordinate space.

Both images are 10X so pixel sizes are equal; the similarity transform only
compensates the stage-encoder offset between acquisition sessions.

Usage:
    python make_sub10x_overlay.py
"""
from __future__ import annotations

import csv
import math
from pathlib import Path

import nd2 as _nd2lib
import numpy as np
from scipy.ndimage import map_coordinates
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
    img = arr
    p_ax = ax("P")
    if p_ax is not None:
        sl = [slice(None)] * img.ndim; sl[p_ax] = 0
        img = img[tuple(sl)]
        dims2 = [d for d in dims if d != "P"]
        z_ax  = dims2.index("Z") if "Z" in dims2 else None
        c_ax  = dims2.index("C") if "C" in dims2 else None
    else:
        dims2 = dims
        z_ax  = ax("Z"); c_ax = ax("C")
    if c_ax is not None:
        sl = [slice(None)] * img.ndim; sl[c_ax] = ch
        img = img[tuple(sl)]
        dims2 = [d for d in dims2 if d != "C"]
        z_ax  = dims2.index("Z") if "Z" in dims2 else None
    if z_ax is not None:
        img = img.max(axis=z_ax)
    return img.astype(np.float32)

def _norm(img, lo_pct=0.5, hi_pct=99.5):
    lo = float(np.percentile(img, lo_pct))
    hi = float(np.percentile(img, hi_pct))
    return np.clip((img - lo) / (hi - lo + 1e-8), 0, 1).astype(np.float32)

def stage_to_pixel(x_um, y_um, cx, cy, pum, W, H):
    col = (x_um - cx) / pum + W / 2.0
    row = (y_um - cy) / pum + H / 2.0
    return col, row

# ── load images ───────────────────────────────────────────────────────────────

print("Loading 10X mosaic ...")
with _nd2lib.ND2File(ND2_MOSAIC) as f:
    px_mos     = _pixel_um(f)
    cx_mos, cy_mos = _stage_center(f)
    img_mos    = _load_mip(f, 0)
Hm, Wm = img_mos.shape
print(f"  {Wm}x{Hm} px  pixel={px_mos:.5f} um  centre=({cx_mos:.1f},{cy_mos:.1f})")

print("Loading sub-image 10X ...")
with _nd2lib.ND2File(ND2_SUB) as f:
    px_sub      = _pixel_um(f)
    cx_sub, cy_sub = _stage_center(f)
    img_sub     = _load_mip(f, 0)
Hs, Ws = img_sub.shape
print(f"  {Ws}x{Hs} px  pixel={px_sub:.5f} um  stage=({cx_sub:.1f},{cy_sub:.1f})")

# ── read correction parameters from validation_report.csv ────────────────────

corr_tx  = -1745.23   # defaults from last run
corr_ty  = -1734.27
corr_sc  =  1.1052
corr_rot = -8.476     # degrees

if VAL_CSV.exists():
    with open(VAL_CSV, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    if rows:
        r = rows[0]
        try:
            corr_tx  = float(r["corr_tx_um"])
            corr_ty  = float(r["corr_ty_um"])
            corr_sc  = float(r["corr_scale"])
            corr_rot = float(r["corr_rotation_deg"])
        except (KeyError, ValueError):
            pass

print(f"  Correction: tx={corr_tx:+.2f} um  ty={corr_ty:+.2f} um  "
      f"scale={corr_sc:.4f}  rot={corr_rot:+.3f} deg")

# ── build similarity transform (sub_stage -> mosaic_stage) ───────────────────
# Forward: x_mos = sc*(cos(th)*x_sub - sin(th)*y_sub) + tx
#          y_mos = sc*(sin(th)*x_sub + cos(th)*y_sub) + ty
# Inverse: x_sub = (1/sc)*(cos(th)*(x_mos-tx) + sin(th)*(y_mos-ty))
#          y_sub = (1/sc)*(-sin(th)*(x_mos-tx) + cos(th)*(y_mos-ty))

th = math.radians(corr_rot)
sc = corr_sc
cos_t = math.cos(th);  sin_t = math.sin(th)

def sub_to_mos_stage(x_sub_arr, y_sub_arr):
    x_mos = sc * (cos_t * x_sub_arr - sin_t * y_sub_arr) + corr_tx
    y_mos = sc * (sin_t * x_sub_arr + cos_t * y_sub_arr) + corr_ty
    return x_mos, y_mos

def mos_to_sub_stage(x_mos_arr, y_mos_arr):
    dx = x_mos_arr - corr_tx;  dy = y_mos_arr - corr_ty
    x_sub = (1.0 / sc) * ( cos_t * dx + sin_t * dy)
    y_sub = (1.0 / sc) * (-sin_t * dx + cos_t * dy)
    return x_sub, y_sub

# ── compute bounding box of sub-image in mosaic pixel space ──────────────────
# corners of sub-image in sub-image stage coords
corners_sub_x = [cx_sub + dx_px * px_sub
                  for dx_px in (-Ws/2, Ws/2, Ws/2, -Ws/2)]
corners_sub_y = [cy_sub + dy_px * px_sub
                  for dy_px in (-Hs/2, -Hs/2, Hs/2, Hs/2)]

corners_mos_x, corners_mos_y = sub_to_mos_stage(
    np.array(corners_sub_x), np.array(corners_sub_y))
corners_mos_col, corners_mos_row = zip(
    *[stage_to_pixel(cx, cy, cx_mos, cy_mos, px_mos, Wm, Hm)
      for cx, cy in zip(corners_mos_x, corners_mos_y)])

col_min, col_max = min(corners_mos_col), max(corners_mos_col)
row_min, row_max = min(corners_mos_row), max(corners_mos_row)

print(f"  Sub-image footprint in mosaic: col=[{col_min:.0f},{col_max:.0f}]  "
      f"row=[{row_min:.0f},{row_max:.0f}]")

# ── sub-image centre in mosaic pixel coords ───────────────────────────────────
cx_sub_mos, cy_sub_mos = sub_to_mos_stage(cx_sub, cy_sub)
sub_col_mos, sub_row_mos = stage_to_pixel(cx_sub_mos, cy_sub_mos,
                                           cx_mos, cy_mos, px_mos, Wm, Hm)
print(f"  Sub-image centre in mosaic: ({sub_col_mos:.0f}, {sub_row_mos:.0f}) px")

# ── define output crop region (mosaic crop around sub-image footprint) ────────
pad_px   = 80
r0_out   = max(0, int(row_min) - pad_px)
r1_out   = min(Hm, int(row_max) + pad_px)
c0_out   = max(0, int(col_min) - pad_px)
c1_out   = min(Wm, int(col_max) + pad_px)
out_h    = r1_out - r0_out
out_w    = c1_out - c0_out
print(f"  Output crop: {out_w}x{out_h} px  top-left=({c0_out},{r0_out})")

# ── warp sub-image into mosaic pixel space (inverse mapping) ──────────────────
# For every pixel (r, c) in the output crop, find where it comes from in
# the sub-image.

rows_out  = np.arange(r0_out, r1_out, dtype=np.float64)
cols_out  = np.arange(c0_out, c1_out, dtype=np.float64)
C_grid, R_grid = np.meshgrid(cols_out, rows_out)  # each (out_h, out_w)

# output pixel -> mosaic stage coords
x_mos_grid = (C_grid - Wm / 2.0) * px_mos + cx_mos
y_mos_grid = (R_grid - Hm / 2.0) * px_mos + cy_mos

# mosaic stage -> sub-image stage coords (inverse similarity transform)
x_sub_grid, y_sub_grid = mos_to_sub_stage(x_mos_grid, y_mos_grid)

# sub-image stage -> sub-image pixel coords
col_sub_grid = (x_sub_grid - cx_sub) / px_sub + Ws / 2.0
row_sub_grid = (y_sub_grid - cy_sub) / px_sub + Hs / 2.0

# sample sub-image with bilinear interpolation; out-of-bounds -> 0
print("  Warping sub-image into mosaic space ...")
warped = map_coordinates(
    img_sub.astype(np.float64),
    [row_sub_grid.ravel(), col_sub_grid.ravel()],
    order=1,
    mode="constant",
    cval=0.0,
).reshape(out_h, out_w).astype(np.float32)

# mask: pixels that actually fall inside the sub-image
inside = ((col_sub_grid >= 0) & (col_sub_grid < Ws) &
          (row_sub_grid >= 0) & (row_sub_grid < Hs)).astype(np.float32)

# ── load screener ranks ───────────────────────────────────────────────────────
screener: list[dict] = []
if SCREEN_CSV.exists():
    with open(SCREEN_CSV, newline="", encoding="utf-8") as fh:
        screener = list(csv.DictReader(fh))

# ── Figure 1: full mosaic with sub-image footprint ───────────────────────────

mosaic_norm = _norm(img_mos)
mosaic_rgb  = np.stack([mosaic_norm] * 3, axis=-1)

# Blend warped sub-image (blue tint) into the full mosaic RGB
warped_norm = _norm(warped)
alpha       = 0.50
for row_i in range(out_h):
    for col_i in range(out_w):
        if inside[row_i, col_i] > 0.5:
            v   = float(warped_norm[row_i, col_i])
            r_m = r0_out + row_i
            c_m = c0_out + col_i
            if 0 <= r_m < Hm and 0 <= c_m < Wm:
                orig = mosaic_rgb[r_m, c_m]
                tinted = np.array([v * 0.3, v * 0.3, v])   # blue tint
                mosaic_rgb[r_m, c_m] = np.clip(
                    (1 - alpha) * orig + alpha * tinted, 0, 1)

# polygon (footprint) for the sub-image in mosaic coords
poly_xs = list(corners_mos_col) + [corners_mos_col[0]]
poly_ys = list(corners_mos_row) + [corners_mos_row[0]]

fig1, ax1 = plt.subplots(figsize=(14, 10))
ax1.imshow(mosaic_rgb, origin="upper", interpolation="nearest")
ax1.plot(poly_xs, poly_ys, "b--", lw=2, label="10X sub-image FOV (corrected)")
ax1.plot(sub_col_mos, sub_row_mos, "b+", ms=14, mew=2.5, label="Sub-image centre")

# screener detections
for s in screener:
    try:
        sx_ = float(s["x_um"]); sy_ = float(s["y_um"])
        col_, row_ = stage_to_pixel(sx_, sy_, cx_mos, cy_mos, px_mos, Wm, Hm)
        dpx = float(s["diameter_um"]) / px_mos / 2
        ax1.add_patch(mpatches.Circle((col_, row_), dpx,
                                      color="yellow", fill=False, lw=1))
        ax1.text(col_, row_ - dpx - 5, s["rank"], color="yellow",
                 fontsize=6, ha="center", va="bottom")
    except Exception:
        pass

ax1.set_title(
    f"10X Mosaic + 10X Sub-image Overlay (blue tint = sub-image FOV)\n"
    f"Rigid correction: tx={corr_tx:+.1f} um  ty={corr_ty:+.1f} um  "
    f"scale={corr_sc:.4f}  rot={corr_rot:+.2f} deg",
    fontsize=10)
ax1.legend(loc="upper right", fontsize=8, framealpha=0.75)
ax1.axis("off")

out1 = OUT_DIR / "mosaic_sub10x_overlay_full.png"
fig1.savefig(out1, dpi=150, bbox_inches="tight")
print(f"Saved: {out1}")
plt.close(fig1)

# ── Figure 2: zoom — side-by-side comparison in the sub-image region ──────────

# Local coords within the crop
def to_local(col, row): return col - c0_out, row - r0_out

mosaic_crop      = img_mos[r0_out:r1_out, c0_out:c1_out]
mosaic_crop_norm = _norm(mosaic_crop)

# Build blended crop (mosaic + warped sub-image)
blend_rgb = np.stack([mosaic_crop_norm] * 3, axis=-1)
for row_i in range(out_h):
    for col_i in range(out_w):
        if inside[row_i, col_i] > 0.5:
            v   = float(warped_norm[row_i, col_i])
            orig = blend_rgb[row_i, col_i]
            tinted = np.array([v * 0.3, v * 0.3, v])
            blend_rgb[row_i, col_i] = np.clip(
                (1 - alpha) * orig + alpha * tinted, 0, 1)

fig2, axes = plt.subplots(1, 3, figsize=(18, 6))
fig2.suptitle("10X Sub-image on 10X Mosaic — Zoom + Side-by-Side", fontsize=11)

# Panel A: mosaic-only crop
ax = axes[0]
ax.imshow(mosaic_crop_norm, cmap="gray", origin="upper", interpolation="nearest")
for xs, ys in zip(
        [c - c0_out for c in poly_xs],
        [r - r0_out for r in poly_ys]):
    pass  # just draw the polygon below
ax.plot([c - c0_out for c in poly_xs],
        [r - r0_out for r in poly_ys], "b--", lw=1.5)
ax.plot(sub_col_mos - c0_out, sub_row_mos - r0_out, "b+", ms=12, mew=2)
for s in screener:
    try:
        sx_ = float(s["x_um"]); sy_ = float(s["y_um"])
        col_, row_ = stage_to_pixel(sx_, sy_, cx_mos, cy_mos, px_mos, Wm, Hm)
        col_l = col_ - c0_out;  row_l = row_ - r0_out
        dpx = float(s["diameter_um"]) / px_mos / 2
        if -dpx < col_l < out_w + dpx and -dpx < row_l < out_h + dpx:
            ax.add_patch(mpatches.Circle((col_l, row_l), dpx,
                                         color="yellow", fill=False, lw=1.2))
            ax.text(col_l, row_l - dpx - 4, s["rank"], color="yellow",
                    fontsize=8, ha="center")
    except Exception:
        pass
ax.set_title("10X Mosaic only\n(blue box = sub-image FOV)", fontsize=9)
ax.axis("off")

# Panel B: blended overlay
ax = axes[1]
ax.imshow(blend_rgb, origin="upper", interpolation="nearest")
ax.plot([c - c0_out for c in poly_xs],
        [r - r0_out for r in poly_ys], "b--", lw=1.5)
ax.plot(sub_col_mos - c0_out, sub_row_mos - r0_out, "b+", ms=12, mew=2)
for s in screener:
    try:
        sx_ = float(s["x_um"]); sy_ = float(s["y_um"])
        col_, row_ = stage_to_pixel(sx_, sy_, cx_mos, cy_mos, px_mos, Wm, Hm)
        col_l = col_ - c0_out;  row_l = row_ - r0_out
        dpx = float(s["diameter_um"]) / px_mos / 2
        if -dpx < col_l < out_w + dpx and -dpx < row_l < out_h + dpx:
            ax.add_patch(mpatches.Circle((col_l, row_l), dpx,
                                         color="yellow", fill=False, lw=1.2))
            ax.text(col_l, row_l - dpx - 4, s["rank"], color="yellow",
                    fontsize=8, ha="center")
    except Exception:
        pass
ax.set_title("Blended overlay\n(blue tint = warped sub-image)", fontsize=9)
ax.axis("off")

# Panel C: warped sub-image alone
ax = axes[2]
ax.imshow(warped_norm * inside, cmap="gray", origin="upper",
          vmin=0, vmax=1, interpolation="nearest")
ax.plot([c - c0_out for c in poly_xs],
        [r - r0_out for r in poly_ys], "b--", lw=1.5)
ax.plot(sub_col_mos - c0_out, sub_row_mos - r0_out, "b+", ms=12, mew=2)
ax.set_title(
    f"Warped 10X sub-image (in mosaic space)\n"
    f"scale={corr_sc:.4f}  rot={corr_rot:+.2f} deg", fontsize=9)
ax.axis("off")

fig2.tight_layout()
out2 = OUT_DIR / "mosaic_sub10x_overlay_zoom.png"
fig2.savefig(out2, dpi=150, bbox_inches="tight")
print(f"Saved: {out2}")
plt.close(fig2)

print("\nDone.")
