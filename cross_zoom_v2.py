"""
cross_zoom_v2.py  —  Sub-10X Stage Coordinate Verification (v2 workflow)

REAL EXPERIMENTAL WORKFLOW
--------------------------
  Step 1  10X whole-well mosaic scan
            → spheroid_screener.py produces ranked CSV with stage coords
              (mosaic session stage frame)

  Step 2  Operator navigates stage to approximate Rank-N area
            → acquires sub-10X single-position snapshot
              (sub-image session stage frame — may differ from mosaic frame)

  Step 3  THIS SCRIPT
            a. Detect spheroids in sub-10X image
            b. Compute their stage coordinates from sub-10X nd2 metadata
               (pixel_offset × pixel_um + nd2_stage_centre)
            c. Match each sub-10X spheroid to a mosaic rank by NCC
               (image patches at same 10X magnification → direct pixel comparison)
            d. Report, per matched pair:
                 sub-10X stage coords  (what the microscope reads from the nd2)
                 mosaic stage coords   (what the screener recorded)
                 Δx, Δy               (session-to-session stage offset)

  Step 4  Operator uses the sub-10X stage coordinates to drive 20X capture
            (NOT the mosaic screener coords — those are in a different frame)

Key verification question
--------------------------
  "Does the nd2 stage metadata of the sub-10X image agree with the mosaic
   screener stage coordinates for the same physical spheroid?"

  If yes  → Δ ≈ 0 → screener coords can drive 20X directly
  If no   → Δ reported → sub-10X nd2 coords must be used for 20X

Usage
-----
  python cross_zoom_v2.py \\
      --mosaic   "demo_data/.../WellD05_Channel555_spIII1_Seq0000.nd2" \\
      --sub10x   "demo_data/.../spheroid 1 10x 1P 555nm.nd2" \\
      --screen   "demo_data/register_out/verified_screen.csv" \\
      --out-dir  "demo_data/register_out"

Dependencies: nd2 numpy scipy scikit-image matplotlib
"""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

_MISSING: list[str] = []
try:
    import nd2 as _nd2lib
except ImportError:
    _nd2lib = None; _MISSING.append("nd2")
try:
    import numpy as np
except ImportError:
    np = None; _MISSING.append("numpy")
try:
    from scipy import ndimage as _ndi
except ImportError:
    _ndi = None; _MISSING.append("scipy")
try:
    from skimage.feature import match_template, peak_local_max as _plm
    from skimage.measure import regionprops
    from skimage.segmentation import watershed as _watershed
    _SKI_OK = True
except ImportError:
    _SKI_OK = False; _MISSING.append("scikit-image")
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    _MPL_OK = True
except ImportError:
    _MPL_OK = False


# ── ND2 helpers ───────────────────────────────────────────────────────────────

def _pixel_um(f) -> float:
    try:
        v = float(f.voxel_size().x)
        return v if v > 0 else 0.644
    except Exception:
        return 0.644


def _stage_center(f) -> tuple[float, float]:
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
        dims2 = dims; z_ax = ax("Z"); c_ax = ax("C")
    if c_ax is not None:
        sl = [slice(None)] * img.ndim; sl[c_ax] = ch
        img = img[tuple(sl)]
        dims2 = [d for d in dims2 if d != "C"]
        z_ax  = dims2.index("Z") if "Z" in dims2 else None
    if z_ax is not None:
        img = img.max(axis=z_ax)
    return img.astype(np.float32)


# ── coordinate helpers ────────────────────────────────────────────────────────

def pixel_to_stage(col: float, row: float,
                   cx: float, cy: float,
                   pum: float, W: int, H: int) -> tuple[float, float]:
    """Pixel position → stage coordinates using nd2 image-centre metadata."""
    return cx + (col - W / 2.0) * pum, cy + (row - H / 2.0) * pum


def stage_to_pixel(x: float, y: float,
                   cx: float, cy: float,
                   pum: float, W: int, H: int) -> tuple[float, float]:
    return (x - cx) / pum + W / 2.0, (y - cy) / pum + H / 2.0


# ── spheroid detector (same logic as screener) ────────────────────────────────

def detect_centroids(img: np.ndarray, pum: float,
                     min_diam_um: float = 50.0,
                     max_diam_um: float = 400.0,
                     sigma: float = 2.0,
                     thresh_k: float = 0.5) -> list[dict]:
    """
    Detect spheroids and return list of dicts with pixel centroid + stage coords.
    """
    if not _SKI_OK or _ndi is None:
        return []
    H, W = img.shape
    blurred = _ndi.gaussian_filter(img, sigma)
    thr     = blurred.mean() + thresh_k * blurred.std()
    filled  = _ndi.binary_fill_holes(blurred > thr)
    dist    = _ndi.distance_transform_edt(filled)
    min_d   = max(4, int(min_diam_um * 0.5 / pum))
    coords  = _plm(dist, min_distance=min_d, labels=filled)
    pm = np.zeros(dist.shape, dtype=bool)
    if coords.size:
        pm[tuple(coords.T)] = True
    markers, _ = _ndi.label(pm)
    labeled = _watershed(-dist, markers, mask=filled)
    min_a = math.pi * (min_diam_um / pum / 2) ** 2
    max_a = math.pi * (max_diam_um / pum / 2) ** 2
    dets: list[dict] = []
    for rp in regionprops(labeled):
        if min_a <= rp.area <= max_a:
            cy_px, cx_px = rp.centroid
            diam = 2.0 * math.sqrt(rp.area / math.pi) * pum
            dets.append({
                "cx_px": cx_px, "cy_px": cy_px,
                "diam_um": round(diam, 2),
                "area_px": int(rp.area),
            })
    return dets


# ── NCC patch comparison (both images are 10X → same pixel size) ──────────────

def _norm(img: np.ndarray) -> np.ndarray:
    lo = float(img.min()); hi = float(img.max())
    return ((img - lo) / (hi - lo + 1e-8)).astype(np.float32)


def _safe_patch(img: np.ndarray, col: float, row: float, half: int) -> np.ndarray | None:
    """Extract a square patch; clips half-size to fit image bounds (min 20 px)."""
    H, W = img.shape
    max_half = min(int(col), W - int(col) - 1,
                   int(row), H - int(row) - 1, half)
    if max_half < 20:
        return None
    h = max_half
    r0 = int(round(row)) - h; r1 = r0 + 2 * h
    c0 = int(round(col)) - h; c1 = c0 + 2 * h
    return img[r0:r1, c0:c1].copy()


def ncc_patch_score(patch_a: np.ndarray, patch_b: np.ndarray) -> float:
    """Normalised cross-correlation between two same-size patches."""
    a = _norm(patch_a).ravel().astype(np.float64)
    b = _norm(patch_b).ravel().astype(np.float64)
    a -= a.mean(); b -= b.mean()
    denom = (np.linalg.norm(a) * np.linalg.norm(b))
    return float(np.dot(a, b) / denom) if denom > 1e-8 else 0.0


# ── main verification function ────────────────────────────────────────────────

def verify(nd2_mosaic: Path, nd2_sub10x: Path,
           screener_csv: Path, out_dir: Path,
           ch: int = 0, patch_half: int = 100) -> list[dict]:
    """
    Match sub-10X detections to mosaic screener ranks by NCC (same-magnification
    image comparison), then report the stage offset per matched pair.
    """
    if _MISSING:
        raise ImportError(f"Missing: {', '.join(_MISSING)}")

    out_dir.mkdir(parents=True, exist_ok=True)

    # ── load mosaic ───────────────────────────────────────────────────────────
    print(f"Loading mosaic: {nd2_mosaic.name}")
    with _nd2lib.ND2File(nd2_mosaic) as f:
        px_mos        = _pixel_um(f)
        cx_mos, cy_mos = _stage_center(f)
        img_mos       = _load_mip(f, ch)
    Hm, Wm = img_mos.shape
    print(f"  {Wm}x{Hm} px  pixel={px_mos:.5f} um  centre=({cx_mos:.1f},{cy_mos:.1f})")

    # ── load screener CSV ─────────────────────────────────────────────────────
    screener: list[dict] = []
    with open(screener_csv, newline="", encoding="utf-8") as fh:
        screener = list(csv.DictReader(fh))
    print(f"Screener: {len(screener)} ranks")

    # Pre-compute mosaic patches for each screener rank
    mos_patches: dict[int, np.ndarray | None] = {}
    for s in screener:
        rk   = int(s["rank"])
        col_ = float(s["canvas_col"])
        row_ = float(s["canvas_row"])
        mos_patches[rk] = _safe_patch(img_mos, col_, row_, patch_half)

    # ── load sub-10X ─────────────────────────────────────────────────────────
    print(f"\nLoading sub-10X: {nd2_sub10x.name}")
    with _nd2lib.ND2File(nd2_sub10x) as f:
        px_sub          = _pixel_um(f)
        cx_sub, cy_sub  = _stage_center(f)
        img_sub         = _load_mip(f, ch)
    Hs, Ws = img_sub.shape
    print(f"  {Ws}x{Hs} px  pixel={px_sub:.5f} um  centre=({cx_sub:.1f},{cy_sub:.1f})")
    print(f"  Sub-10X nd2 stage centre: ({cx_sub:.2f}, {cy_sub:.2f}) um")

    # ── detect spheroids in sub-10X ───────────────────────────────────────────
    print("\nDetecting spheroids in sub-10X image ...")
    dets = detect_centroids(img_sub, px_sub)
    print(f"  Found {len(dets)} spheroid(s)")
    for i, d in enumerate(dets):
        # Stage coords from sub-10X nd2 metadata
        sx, sy = pixel_to_stage(d["cx_px"], d["cy_px"],
                                 cx_sub, cy_sub, px_sub, Ws, Hs)
        d["stage_x"] = round(sx, 2)
        d["stage_y"] = round(sy, 2)
        print(f"  S{i+1}: pixel=({d['cx_px']:.0f},{d['cy_px']:.0f})  "
              f"stage=({sx:.1f},{sy:.1f})  diam={d['diam_um']:.1f} um")

    # ── NCC matching: sub-10X detections → mosaic ranks ───────────────────────
    # Both images are 10X → same pixel size → patches are directly comparable
    print("\nNCC patch matching (sub-10X vs mosaic, same 10X magnification) ...")
    ncc_matrix: dict[tuple[int,int], float] = {}   # (sub_idx, rank) -> ncc

    for si, d in enumerate(dets):
        sub_patch = _safe_patch(img_sub, d["cx_px"], d["cy_px"], patch_half)
        if sub_patch is None:
            continue
        for s in screener:
            rk = int(s["rank"])
            mp = mos_patches[rk]
            if mp is None:
                continue
            # Resize to same shape if slightly different (shouldn't happen here)
            min_h = min(sub_patch.shape[0], mp.shape[0])
            min_w = min(sub_patch.shape[1], mp.shape[1])
            score = ncc_patch_score(sub_patch[:min_h,:min_w], mp[:min_h,:min_w])
            ncc_matrix[(si, rk)] = score

    NCC_MIN_ACCEPT = 0.40   # discard matches below this threshold

    # Greedy best-match assignment (each sub-det → best unmatched rank)
    used_ranks: set[int] = set()
    matches: list[dict] = []

    for si, d in enumerate(dets):
        scores = {rk: ncc_matrix.get((si, rk), -1.0)
                  for s in screener
                  for rk in [int(s["rank"])]
                  if rk not in used_ranks}
        if not scores:
            continue
        best_rank = max(scores, key=lambda r: scores[r])
        best_ncc  = scores[best_rank]

        if best_ncc < NCC_MIN_ACCEPT:
            print(f"  S{si+1}: best NCC={best_ncc:.4f} < {NCC_MIN_ACCEPT} "
                  f"-- no confident match (near-edge or small patch)")
            continue

        # Mosaic stage coords for this rank (from screener)
        sc_rec   = next(s for s in screener if int(s["rank"]) == best_rank)
        x_mos_rk = float(sc_rec["x_um"])
        y_mos_rk = float(sc_rec["y_um"])

        # Stage offset: sub-10X nd2 stage − mosaic screener stage
        dx = d["stage_x"] - x_mos_rk
        dy = d["stage_y"] - y_mos_rk
        err = math.sqrt(dx**2 + dy**2)

        used_ranks.add(best_rank)
        matches.append({
            "sub_idx":          si + 1,
            "sub_cx_px":        round(d["cx_px"], 1),
            "sub_cy_px":        round(d["cy_px"], 1),
            "sub_stage_x":      d["stage_x"],
            "sub_stage_y":      d["stage_y"],
            "matched_rank":     best_rank,
            "ncc":              round(best_ncc, 4),
            "mosaic_stage_x":   round(x_mos_rk, 2),
            "mosaic_stage_y":   round(y_mos_rk, 2),
            "offset_x_um":      round(dx, 2),
            "offset_y_um":      round(dy, 2),
            "offset_mag_um":    round(err, 2),
            "diam_um":          d["diam_um"],
            "mosaic_col":       float(sc_rec["canvas_col"]),
            "mosaic_row":       float(sc_rec["canvas_row"]),
        })

        print(f"  S{si+1} -> Rank {best_rank}  NCC={best_ncc:.4f}")
        print(f"    Sub-10X  stage = ({d['stage_x']:.2f}, {d['stage_y']:.2f}) um  "
              f"[from nd2 metadata]")
        print(f"    Mosaic   stage = ({x_mos_rk:.2f}, {y_mos_rk:.2f}) um  "
              f"[from screener CSV]")
        print(f"    Offset:  dx={dx:+.1f} um  dy={dy:+.1f} um  "
              f"|offset|={err:.1f} um")

    # Summary
    print("\n--- Stage offset summary ---")
    print(f"  {'Sub':>4}  {'Rank':>5}  {'NCC':>6}  "
          f"{'dx (um)':>9}  {'dy (um)':>9}  {'|offset| (um)':>14}")
    for m in matches:
        print(f"  {m['sub_idx']:>4}  {m['matched_rank']:>5}  {m['ncc']:>6.4f}  "
              f"{m['offset_x_um']:>+9.1f}  {m['offset_y_um']:>+9.1f}  "
              f"{m['offset_mag_um']:>14.1f}")

    if matches:
        dxs = [m["offset_x_um"] for m in matches]
        dys = [m["offset_y_um"] for m in matches]
        print(f"\n  Mean offset:  dx={sum(dxs)/len(dxs):+.1f} um  "
              f"dy={sum(dys)/len(dys):+.1f} um")
        print(f"  Interpretation:")
        print(f"    When the 20X capture uses the sub-10X nd2 stage coordinates,")
        print(f"    the physical position will be offset from the mosaic screener")
        print(f"    prediction by approximately ({sum(dxs)/len(dxs):+.0f}, "
              f"{sum(dys)/len(dys):+.0f}) um.")
        print(f"    This offset is the stage encoder re-zero between sessions.")

    # ── write CSV ─────────────────────────────────────────────────────────────
    out_csv = out_dir / "sub10x_stage_offset.csv"
    fieldnames = [
        "sub_idx", "matched_rank", "ncc",
        "sub_stage_x", "sub_stage_y",
        "mosaic_stage_x", "mosaic_stage_y",
        "offset_x_um", "offset_y_um", "offset_mag_um",
        "diam_um",
    ]
    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames,
                           extrasaction="ignore")
        w.writeheader(); w.writerows(matches)
    print(f"\nWrote: {out_csv}")

    # ── figure ─────────────────────────────────────────────────────────────────
    if _MPL_OK and matches:
        _plot_verification(
            img_sub, img_mos, px_sub, px_mos,
            cx_sub, cy_sub, cx_mos, cy_mos,
            Ws, Hs, Wm, Hm,
            dets, matches, screener, mos_patches, patch_half,
            out_dir)

    return matches


def _norm_display(img, lo=0.5, hi=99.5):
    lo_ = float(np.percentile(img, lo))
    hi_ = float(np.percentile(img, hi))
    return np.clip((img - lo_) / (hi_ - lo_ + 1e-8), 0, 1).astype(np.float32)


def _plot_verification(img_sub, img_mos, px_sub, px_mos,
                        cx_sub, cy_sub, cx_mos, cy_mos,
                        Ws, Hs, Wm, Hm,
                        dets, matches, screener, mos_patches, patch_half,
                        out_dir: Path):

    n = len(matches)
    # Layout: overview row + n detail rows, each with 4 panels
    fig = plt.figure(figsize=(18, 5 + 5 * n))
    fig.suptitle(
        "Cross-Zoom V2: Sub-10X Stage Coordinate Verification\n"
        "Sub-10X nd2 stage coords vs Mosaic screener stage coords (same spheroid)",
        fontsize=12, y=0.98)

    # ── Row 0: overview panels ────────────────────────────────────────────────
    ax_sub = fig.add_subplot(n + 1, 4, 1)
    ax_sub.imshow(_norm_display(img_sub), cmap="gray",
                  origin="upper", interpolation="nearest")
    for i, d in enumerate(dets):
        r_px = d["diam_um"] / px_sub / 2
        ax_sub.add_patch(mpatches.Circle((d["cx_px"], d["cy_px"]), r_px,
                                          color="cyan", fill=False, lw=1.5))
        m = next((m for m in matches if m["sub_idx"] == i+1), None)
        lbl = f"S{i+1}\nRk{m['matched_rank']}" if m else f"S{i+1}"
        ax_sub.text(d["cx_px"], d["cy_px"] - r_px - 6, lbl,
                    color="cyan", fontsize=7, ha="center")
        # Stage coord annotation
        ax_sub.text(d["cx_px"], d["cy_px"] + r_px + 10,
                    f"({d['stage_x']:.0f},{d['stage_y']:.0f})",
                    color="yellow", fontsize=6, ha="center")
    ax_sub.set_title("Sub-10X image\n(cyan circle = detected; "
                     "yellow = nd2 stage coords)", fontsize=8)
    ax_sub.axis("off")

    # Mosaic overview (crop around the matched ranks)
    ax_mos = fig.add_subplot(n + 1, 4, 2)
    all_cols = [m["mosaic_col"] for m in matches]
    all_rows = [m["mosaic_row"] for m in matches]
    # Ensure minimum extent so we always see something
    span = max(max(all_cols) - min(all_cols), max(all_rows) - min(all_rows), 1)
    pad_ov = max(int(300 / px_mos), int(span * 0.5))
    r0m = max(0, int(min(all_rows)) - pad_ov)
    r1m = min(Hm, int(max(all_rows)) + pad_ov)
    c0m = max(0, int(min(all_cols)) - pad_ov)
    c1m = min(Wm, int(max(all_cols)) + pad_ov)
    mos_crop_ov = img_mos[r0m:r1m, c0m:c1m]
    ax_mos.imshow(_norm_display(mos_crop_ov), cmap="gray",
                  origin="upper", interpolation="nearest")
    for s in screener:
        col_ = float(s["canvas_col"]) - c0m
        row_ = float(s["canvas_row"]) - r0m
        r_px = float(s["diameter_um"]) / px_mos / 2
        if -r_px < col_ < mos_crop_ov.shape[1]+r_px and \
           -r_px < row_ < mos_crop_ov.shape[0]+r_px:
            clr = "lime" if int(s["rank"]) in {m["matched_rank"] for m in matches} \
                  else "yellow"
            ax_mos.add_patch(mpatches.Circle((col_, row_), r_px,
                                              color=clr, fill=False, lw=1.2))
            ax_mos.text(col_, row_ - r_px - 4, s["rank"],
                        color=clr, fontsize=7, ha="center")
    ax_mos.set_title("Mosaic 10X (zoom)\n(green = matched ranks)", fontsize=8)
    ax_mos.axis("off")

    # Stage offset summary text panel
    ax_txt = fig.add_subplot(n + 1, 4, (3, 4))
    ax_txt.axis("off")
    lines = [
        "STAGE OFFSET SUMMARY",
        "(sub-10X nd2 stage  −  mosaic screener stage)",
        "",
        f"{'Sub':>4}  {'Rank':>5}  {'NCC':>6}  {'dx (um)':>9}  "
        f"{'dy (um)':>9}  {'|err| (um)':>11}",
        "─" * 52,
    ]
    for m in matches:
        lines.append(
            f"{'S'+str(m['sub_idx']):>4}  {m['matched_rank']:>5}  "
            f"{m['ncc']:>6.4f}  {m['offset_x_um']:>+9.1f}  "
            f"{m['offset_y_um']:>+9.1f}  {m['offset_mag_um']:>11.1f}")
    lines += ["", "─" * 52]
    if matches:
        dxs = [m["offset_x_um"] for m in matches]
        dys = [m["offset_y_um"] for m in matches]
        lines.append(f"  Mean offset:  dx={sum(dxs)/len(dxs):+.1f} um  "
                     f"dy={sum(dys)/len(dys):+.1f} um")
        lines.append("")
        lines.append("INTERPRETATION")
        lines.append("  Stage coords from sub-10X nd2 metadata differ from")
        lines.append(f"  mosaic screener by ({sum(dxs)/len(dxs):+.0f}, "
                     f"{sum(dys)/len(dys):+.0f}) um.")
        lines.append("  The sub-10X nd2 coords ARE the correct coords")
        lines.append("  to use for 20X capture (Step 4).")
        lines.append("  The mosaic screener coords require this offset")
        lines.append("  applied BEFORE driving the stage to 20X position.")
    ax_txt.text(0.02, 0.98, "\n".join(lines), transform=ax_txt.transAxes,
                va="top", ha="left", fontsize=9,
                fontfamily="monospace",
                bbox=dict(boxstyle="round", fc="black", ec="gray", alpha=0.6))

    # ── Rows 1..n: per-match detail ───────────────────────────────────────────
    for i, m in enumerate(matches):
        row_base = (i + 1) * 4 + 1   # subplot indices start at 1

        si_idx  = m["sub_idx"] - 1
        d       = dets[si_idx]
        rk      = m["matched_rank"]
        sc_rec  = next(s for s in screener if int(s["rank"]) == rk)
        diam_px_sub = d["diam_um"] / px_sub / 2
        diam_px_mos = d["diam_um"] / px_mos / 2

        # Sub-10X patch
        ax = fig.add_subplot(n + 1, 4, row_base)
        sp = _safe_patch(img_sub, d["cx_px"], d["cy_px"], patch_half)
        if sp is not None:
            ax.imshow(_norm_display(sp), cmap="gray",
                      origin="upper", interpolation="nearest")
            ax.plot(patch_half, patch_half, "c+", ms=14, mew=2)
            ax.add_patch(mpatches.Circle((patch_half, patch_half), diam_px_sub,
                                          color="cyan", fill=False, lw=1.5))
        ax.set_title(
            f"S{m['sub_idx']} in sub-10X\n"
            f"nd2 stage: ({m['sub_stage_x']:.1f}, {m['sub_stage_y']:.1f}) um",
            fontsize=8)
        ax.set_ylabel(f"Matched to Rank {rk}\nNCC={m['ncc']:.4f}", fontsize=8)
        ax.axis("off")

        # Mosaic patch at screener position
        ax = fig.add_subplot(n + 1, 4, row_base + 1)
        mp = mos_patches[rk]
        if mp is not None:
            ax.imshow(_norm_display(mp), cmap="gray",
                      origin="upper", interpolation="nearest")
            ax.plot(patch_half, patch_half, "y+", ms=14, mew=2)
            ax.add_patch(mpatches.Circle((patch_half, patch_half), diam_px_mos,
                                          color="yellow", fill=False, lw=1.5))
        ax.set_title(
            f"Rank {rk} in mosaic\n"
            f"screener stage: ({m['mosaic_stage_x']:.1f}, "
            f"{m['mosaic_stage_y']:.1f}) um",
            fontsize=8)
        ax.axis("off")

        # Side-by-side comparison
        ax = fig.add_subplot(n + 1, 4, row_base + 2)
        if sp is not None and mp is not None:
            min_h = min(sp.shape[0], mp.shape[0])
            min_w = min(sp.shape[1], mp.shape[1])
            combined = np.concatenate([
                _norm_display(sp[:min_h,:min_w]),
                np.ones((min_h, 4), dtype=np.float32) * 0.5,   # divider
                _norm_display(mp[:min_h,:min_w]),
            ], axis=1)
            ax.imshow(combined, cmap="gray", origin="upper",
                      interpolation="nearest")
            ax.axvline(min_w + 2, color="white", lw=1, linestyle="--")
            ax.text(min_w // 2, -8, "Sub-10X", color="cyan",
                    ha="center", fontsize=8)
            ax.text(min_w + 4 + min_w // 2, -8, "Mosaic", color="yellow",
                    ha="center", fontsize=8)
        ax.set_title(f"Side-by-side comparison\n"
                     f"NCC = {m['ncc']:.4f}", fontsize=8)
        ax.axis("off")

        # Stage offset arrow
        ax = fig.add_subplot(n + 1, 4, row_base + 3)
        ax.axis("off")
        ax.set_facecolor("black")
        dx = m["offset_x_um"]; dy = m["offset_y_um"]
        err = m["offset_mag_um"]
        # Draw offset arrow in a schematic
        ax.set_xlim(-1.2, 1.2); ax.set_ylim(-1.2, 1.2)
        ax.plot(0, 0, "o", color="yellow", ms=12, label="Mosaic stage")
        ex = dx / (max(abs(dx), abs(dy), 1e-3)) * 0.7
        ey = dy / (max(abs(dx), abs(dy), 1e-3)) * 0.7
        ax.annotate("", xy=(ex, ey), xytext=(0, 0),
                    arrowprops=dict(arrowstyle="->", color="cyan", lw=2))
        ax.plot(ex, ey, "o", color="cyan", ms=12, label="Sub-10X nd2 stage")
        ax.legend(fontsize=7, loc="upper left")
        ax.set_title(
            f"Stage offset (nd2 - screener)\n"
            f"dx = {dx:+.1f} um\n"
            f"dy = {dy:+.1f} um\n"
            f"|offset| = {err:.1f} um",
            fontsize=8, color="white")
        ax.set_facecolor("#111111")
        for spine in ax.spines.values():
            spine.set_edgecolor("gray")

    fig.tight_layout(rect=[0, 0, 1, 0.97])
    out_png = out_dir / "cross_zoom_v2_verification.png"
    fig.savefig(out_png, dpi=150, bbox_inches="tight", facecolor="black")
    print(f"Saved figure: {out_png}")
    plt.close(fig)


# ── CLI ────────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description="Sub-10X stage coordinate verification (v2 workflow)")
    ap.add_argument("--mosaic",  type=Path, required=True,
                    help="10X whole-well mosaic nd2")
    ap.add_argument("--sub10x",  type=Path, required=True,
                    help="Single-position sub-10X nd2")
    ap.add_argument("--screen",  type=Path, required=True,
                    help="Screener CSV (verified_screen.csv)")
    ap.add_argument("--out-dir", type=Path,
                    default=Path("demo_data/register_out"),
                    help="Output directory for CSV and figures")
    ap.add_argument("--ch",      type=int, default=0,
                    help="Channel index (default 0)")
    ap.add_argument("--patch",   type=int, default=100,
                    help="Patch half-size in pixels for NCC matching (default 100)")
    args = ap.parse_args()

    verify(
        nd2_mosaic   = args.mosaic,
        nd2_sub10x   = args.sub10x,
        screener_csv = args.screen,
        out_dir      = args.out_dir,
        ch           = args.ch,
        patch_half   = args.patch,
    )


if __name__ == "__main__":
    main()
