"""
cross_zoom_register.py
Cross-Zoom Coordinate Verification: 10X Whole-Well → 20X Single-Position

Sits between spheroid_screener.py and the NIS-E capture stage
(nis_macro_z_autofocus.mac / nis_macro_auto_capture.mac).

Background
----------
NIS-E outputs the 10X scan as a single stitched whole-well nd2
(e.g. WellD05_Channel555_spIII1_Seq0000.nd2, ~10 k×10 k px, P=1).
The screener detects spheroids in that image and stores their stage
coordinates (µm) in spheroid_screen_latest.csv.

This script answers the question:
  "Does the screener's (x_um, y_um) actually land on the same
   spheroid in the 20X image?"

Stage → pixel conversion for a P=1 stitched whole-well image:
  col = (x_um - center_x_um) / pixel_um  +  width  / 2
  row = (y_um - center_y_um) / pixel_um  +  height / 2
where (center_x_um, center_y_um) is the stage position stored in
frame_metadata(0) — the physical centre of the stitched image.

Two modes
---------
pre   (default) — run BEFORE 20X capture
  Load the whole-well 10X and the screener CSV.
  For every detected spheroid, extract a crop from the whole-well
  image at the screener coordinates and, if requested, refine the
  centroid with intra-image NCC (disk template vs. crop).
  Output: verified_screen.csv with NCC-refined coords + QC flags.
  --gui shows the whole-well mosaic with detections and per-spheroid
  zoom-in crops; click to override centroid, a/r to accept/reject.

validate — run AFTER 20X capture
  For each 20X nd2 in --nd2-20x (folder or single file):
    1. Read 20X stage position from metadata
    2. Predict where it falls in the whole-well 10X (same transform)
    3. Extract 10X crop matching the physical extent of the 20X FOV
    4. Scale 20X MIP to 10X pixel space
    5. NCC between the crops → match quality score
  A high NCC score confirms the screener coordinates were accurate;
  a low score flags potential stage drift or mis-detection.
  Output: validation_report.csv + optional --gui side-by-side figures.

Usage
-----
  # Pre-capture verification (refine + QC screen detections)
  python cross_zoom_register.py --mode pre \\
      --nd2-10x  D:\\scan\\WellD05_Channel555_spIII1_Seq0000.nd2 \\
      --screen   D:\\scan\\spheroid_screen_latest.csv

  # Post-capture validation (confirm 20X images match screener coords)
  python cross_zoom_register.py --mode validate \\
      --nd2-10x  D:\\scan\\WellD05_Channel555_spIII1_Seq0000.nd2 \\
      --nd2-20x  D:\\scan\\20X_captures\\

  # Either mode with interactive GUI
  python cross_zoom_register.py --mode pre --gui ...

Dependencies:
  pip install nd2 numpy scipy scikit-image matplotlib
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

# ── optional stack ────────────────────────────────────────────────────────────

_MISSING: list[str] = []

try:
    import nd2 as _nd2lib
except ImportError:
    _nd2lib = None; _MISSING.append("nd2")

try:
    import numpy as np
except ImportError:
    np = None; _MISSING.append("numpy")           # type: ignore[assignment]

try:
    from scipy import ndimage as _ndi
    from scipy.ndimage import zoom as _zoom
except ImportError:
    _ndi = None; _zoom = None; _MISSING.append("scipy")

try:
    from skimage.feature import match_template, peak_local_max as _plm
    from skimage.measure import regionprops, label as _ski_label
    from skimage.segmentation import watershed as _watershed
    from skimage.transform import SimilarityTransform
    _SKI_OK = True
except ImportError:
    _SKI_OK = False; _MISSING.append("scikit-image")

try:
    import matplotlib
    matplotlib.use("TkAgg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    _MPL_OK = True
except ImportError:
    _MPL_OK = False


# ── constants ─────────────────────────────────────────────────────────────────

NCC_SEARCH_UM   = 30.0  # ± search radius for intra-image centroid refinement
NCC_MIN_PEAK    = 0.30  # NCC peak below this → flag "low-ncc"
EDGE_MARGIN_PX  = 30    # pixels within FOV edge → flag "near-edge"
QC_MIN_CIRC     = 0.60
QC_MIN_SOL      = 0.70
MATCH_GOOD_NCC  = 0.50  # validate mode: threshold for "confirmed" match


# ── ND2 helpers ───────────────────────────────────────────────────────────────

def _pixel_um(f) -> float:
    try:
        v = float(f.voxel_size().x)
        return v if v and v > 0 else 0.644
    except Exception:
        return 0.644


def _stage_center(f) -> tuple[float, float, float]:
    """Stage position of the image centre (frame 0, which for P=1 is the only frame)."""
    try:
        fm  = f.frame_metadata(0)
        pos = fm.channels[0].position.stagePositionUm
        return float(pos.x), float(pos.y), float(pos.z)
    except Exception:
        return 0.0, 0.0, 0.0


def _load_mip(f, ch_idx: int = 0) -> np.ndarray:
    """
    Return a 2-D float32 image.
    Handles any dimension order; max-projects over Z if present.
    For a flat P=1 whole-well image (sizes = {Y, X}) this is just arr[].astype.
    """
    arr  = f.asarray()
    dims = list(f.sizes.keys())   # e.g. ['Y','X'] or ['Z','C','Y','X']

    def ax(k): return dims.index(k) if k in dims else None

    p_ax = ax("P"); z_ax = ax("Z"); c_ax = ax("C")
    img  = arr

    if p_ax is not None:
        sl = [slice(None)] * img.ndim; sl[p_ax] = 0
        img = img[tuple(sl)]
        dims2 = [d for d in dims if d != "P"]
        z_ax  = dims2.index("Z") if "Z" in dims2 else None
        c_ax  = dims2.index("C") if "C" in dims2 else None
    else:
        dims2 = dims

    if c_ax is not None:
        sl = [slice(None)] * img.ndim; sl[c_ax] = ch_idx
        img = img[tuple(sl)]
        dims2 = [d for d in dims2 if d != "C"]
        z_ax  = dims2.index("Z") if "Z" in dims2 else None

    if z_ax is not None:
        img = img.max(axis=z_ax)

    return img.astype(np.float32)


# ── coordinate helpers ────────────────────────────────────────────────────────

def stage_to_pixel(x_um: float, y_um: float,
                   center_x: float, center_y: float,
                   pixel_um: float,
                   img_w: int, img_h: int) -> tuple[float, float]:
    """
    Stage µm → pixel (col, row) in a P=1 stitched whole-well image.
    (center_x, center_y) is the stage position of the image centre.
    """
    col = (x_um - center_x) / pixel_um + img_w / 2.0
    row = (y_um - center_y) / pixel_um + img_h / 2.0
    return col, row


def pixel_to_stage(col: float, row: float,
                   center_x: float, center_y: float,
                   pixel_um: float,
                   img_w: int, img_h: int) -> tuple[float, float]:
    x_um = (col - img_w / 2.0) * pixel_um + center_x
    y_um = (row - img_h / 2.0) * pixel_um + center_y
    return x_um, y_um


def _safe_crop(img: np.ndarray, row: float, col: float,
               half_h: int, half_w: int) -> np.ndarray | None:
    """Extract a rectangular crop; return None if it goes out of bounds."""
    r0 = int(round(row)) - half_h
    r1 = int(round(row)) + half_h
    c0 = int(round(col)) - half_w
    c1 = int(round(col)) + half_w
    if r0 < 0 or r1 > img.shape[0] or c0 < 0 or c1 > img.shape[1]:
        return None
    return img[r0:r1, c0:c1].copy()


# ── NCC disk-template refinement (pre mode) ──────────────────────────────────

def _disk_template(radius_px: float, sigma: float = 2.0) -> np.ndarray:
    r  = int(math.ceil(radius_px + 3 * sigma)) + 2
    sz = 2 * r + 1
    Y, X = np.ogrid[:sz, :sz]
    d    = np.sqrt((Y - r) ** 2 + (X - r) ** 2)
    mask = (d <= radius_px).astype(np.float32)
    blr  = _ndi.gaussian_filter(mask, sigma=sigma)
    return blr / (blr.max() + 1e-8)


def ncc_refine_centroid(img_10x: np.ndarray,
                        col0: float, row0: float,
                        diameter_um: float, pixel_um: float,
                        search_um: float = NCC_SEARCH_UM
                        ) -> tuple[float, float, float]:
    """
    Refine a centroid guess (col0, row0) by sliding a soft-disk template
    around it in the whole-well 10X image.

    Returns (refined_col, refined_row, ncc_peak).
    Falls back to the input centroid if the search window is too small.
    """
    if not _SKI_OK:
        return col0, row0, 0.0

    radius_px  = diameter_um / 2.0 / pixel_um
    search_px  = int(round(search_um / pixel_um))
    template   = _disk_template(radius_px)
    th, tw     = template.shape
    pad        = search_px + th // 2 + 4

    r0i = max(0, int(round(row0)) - pad)
    r1i = min(img_10x.shape[0], int(round(row0)) + pad + 1)
    c0i = max(0, int(round(col0)) - pad)
    c1i = min(img_10x.shape[1], int(round(col0)) + pad + 1)
    region = img_10x[r0i:r1i, c0i:c1i]

    if region.shape[0] < th or region.shape[1] < tw:
        return col0, row0, 0.0

    region_n  = (region - region.min()) / ((region.max() - region.min()) + 1e-8)
    ncc_map   = match_template(region_n, template, pad_input=True)

    # Restrict search to ±search_px around original centroid in the crop
    cy = int(round(row0)) - r0i
    cx = int(round(col0)) - c0i
    rs = max(0, cy - search_px); re = min(ncc_map.shape[0], cy + search_px + 1)
    cs = max(0, cx - search_px); ce = min(ncc_map.shape[1], cx + search_px + 1)

    if re <= rs or ce <= cs:
        return col0, row0, 0.0

    sub = ncc_map[rs:re, cs:ce]
    pr, pc    = np.unravel_index(sub.argmax(), sub.shape)
    ncc_peak  = float(sub[pr, pc])
    ref_row   = float(r0i + rs + pr)
    ref_col   = float(c0i + cs + pc)

    return ref_col, ref_row, ncc_peak


# ── QC flag helper ────────────────────────────────────────────────────────────

def _qc_flags(rec: dict, img_h: int, img_w: int) -> str:
    flags: list[str] = []
    col = float(rec.get("canvas_col", 0))
    row = float(rec.get("canvas_row", 0))
    if (row < EDGE_MARGIN_PX or row > img_h - EDGE_MARGIN_PX or
            col < EDGE_MARGIN_PX or col > img_w - EDGE_MARGIN_PX):
        flags.append("near-edge")
    if float(rec.get("circularity", 1)) < QC_MIN_CIRC:
        flags.append("low-circularity")
    if float(rec.get("solidity", 1)) < QC_MIN_SOL:
        flags.append("low-solidity")
    if float(rec.get("ncc_peak", 1)) < NCC_MIN_PEAK:
        flags.append("low-ncc")
    return "|".join(flags)


# ── PRE mode: verify screener output against whole-well 10X ──────────────────

_VERIFIED_COLS = [
    "rank", "position_idx",
    "x_um", "y_um", "z_surface_um",
    "x_um_raw", "y_um_raw",
    "canvas_col", "canvas_row",
    "shift_um",
    "diameter_um", "circularity", "solidity", "score",
    "ncc_peak", "qc_flags", "verified",
]


def pre_verify(nd2_10x: Path, screener_csv: Path, ch_idx: int,
               search_um: float, out_dir: Path, gui: bool = False
               ) -> tuple[list[dict], Path]:
    """
    Load the whole-well 10X stitched nd2 and the screener CSV.
    For each detection, refine the centroid with intra-image NCC and
    apply QC flags.  Optionally open an interactive Matplotlib window.

    Returns (verified_records, verified_csv_path).
    """
    print(f"[pre]  Loading 10X image: {nd2_10x.name}")
    with _nd2lib.ND2File(nd2_10x) as f:
        px_um    = _pixel_um(f)
        cx0, cy0, _ = _stage_center(f)
        img      = _load_mip(f, ch_idx)
    img_h, img_w = img.shape
    print(f"[pre]  Image: {img_w}x{img_h} px  "
          f"({img_w*px_um/1000:.2f}x{img_h*px_um/1000:.2f} mm)  "
          f"pixel={px_um:.5f} um  stage_centre=({cx0:.1f},{cy0:.1f})")

    with open(screener_csv, newline="", encoding="utf-8") as f:
        screener_recs = list(csv.DictReader(f))
    print(f"[pre]  {len(screener_recs)} detection(s) from {screener_csv.name}")

    verified: list[dict] = []
    for rec in screener_recs:
        x_raw = float(rec["x_um"])
        y_raw = float(rec["y_um"])
        diam  = float(rec["diameter_um"])

        col_raw, row_raw = stage_to_pixel(x_raw, y_raw, cx0, cy0,
                                          px_um, img_w, img_h)

        # NCC refinement within the whole-well image
        ref_col, ref_row, ncc_peak = ncc_refine_centroid(
            img, col_raw, row_raw, diam, px_um, search_um)

        x_ref, y_ref = pixel_to_stage(ref_col, ref_row, cx0, cy0,
                                      px_um, img_w, img_h)
        shift_um = math.sqrt((x_ref - x_raw)**2 + (y_ref - y_raw)**2)

        v: dict = {
            "rank":          int(rec["rank"]),
            "position_idx":  int(rec["position_idx"]),
            "x_um":          round(x_ref, 2),
            "y_um":          round(y_ref, 2),
            "z_surface_um":  float(rec["z_surface_um"]),
            "x_um_raw":      round(x_raw, 2),
            "y_um_raw":      round(y_raw, 2),
            "canvas_col":    round(ref_col, 1),
            "canvas_row":    round(ref_row, 1),
            "shift_um":      round(shift_um, 2),
            "diameter_um":   diam,
            "circularity":   float(rec.get("circularity", 0)),
            "solidity":      float(rec.get("solidity", 0)),
            "score":         float(rec.get("score", 0)),
            "ncc_peak":      round(ncc_peak, 4),
            "qc_flags":      "",
            "verified":      True,
        }
        v["qc_flags"] = _qc_flags(v, img_h, img_w)

        print(f"  Rank {v['rank']:2d}  "
              f"({x_ref:8.1f},{y_ref:8.1f}) um  "
              f"pixel=({ref_col:6.0f},{ref_row:6.0f})  "
              f"shift={shift_um:4.1f} um  "
              f"NCC={ncc_peak:.3f}  [{v['qc_flags'] or 'OK'}]")
        verified.append(v)

    if gui:
        if _MPL_OK:
            verified = _gui_review_pre(img, px_um, verified)
        else:
            print("[pre]  matplotlib not installed — skipping GUI")

    out_dir.mkdir(parents=True, exist_ok=True)
    out_csv = out_dir / "verified_screen.csv"
    _write_csv(verified, out_csv, _VERIFIED_COLS)
    n_ok = sum(1 for v in verified if v["verified"])
    print(f"[pre]  {n_ok}/{len(verified)} verified.  Wrote: {out_csv}")
    return verified, out_csv


# ── 2D rigid stage-correction between mosaic and sub-image acquisitions ─────────
#
# NIS-E sometimes uses a different stage origin between the whole-well overview
# scan and individual-position acquisitions (e.g. after an objective change or
# encoder re-zero).  The function below detects spheroid centroids in the sub-
# image 10X, finds the best-fit 2D similarity transform (translation + rotation
# + uniform scale) that maps sub-image stage coords to mosaic stage coords, and
# returns that transform so the 20X metadata coordinates can be corrected before
# NCC validation.

def _detect_sub_centroids(img: np.ndarray, pum: float,
                           min_diam_um: float = 50.0,
                           max_diam_um: float = 400.0,
                           sigma: float = 2.0,
                           thresh_k: float = 0.5) -> list[dict]:
    """Detect spheroid centroids in the sub-image (same logic as screener)."""
    if not _SKI_OK or _ndi is None:
        return []
    blurred = _ndi.gaussian_filter(img, sigma)
    filled  = _ndi.binary_fill_holes(blurred > (blurred.mean() + thresh_k * blurred.std()))
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
    pts = []
    for rp in regionprops(labeled):
        if min_a <= rp.area <= max_a:
            pts.append({
                "cx_px": rp.centroid[1],
                "cy_px": rp.centroid[0],
                "diam":  2.0 * math.sqrt(rp.area / math.pi) * pum,
            })
    return pts


def _fit_stage_correction(
        nd2_sub10x: Path,
        screener_csv: Path,
        ch_idx: int = 0,
) -> dict:
    """
    Fit a 2D similarity transform (tx, ty, scale, rotation) that maps sub-image
    stage coordinates to mosaic stage coordinates.

    Algorithm
    ---------
    1. Detect spheroid centroids in the sub-image 10X.
    2. Exhaustively try every pairing of (sub-image detection → screener rank).
    3. For each anchor pair, compute the implied pure translation, apply it to all
       sub-image detections, and measure the total RMS distance to the nearest
       mosaic rank.  The pair with the lowest RMS is the anchor correspondence.
    4. With the best two corresponding point pairs identified (anchor + its closest
       second point), fit a SimilarityTransform (4 DOF: tx, ty, scale, rotation).
    5. Return the transform parameters and per-detection residuals.

    Returns dict with keys:
      tform          — skimage SimilarityTransform (apply with tform(pts))
      tx, ty         — translation (µm)
      scale          — uniform scale factor
      rotation_deg   — rotation angle in degrees
      correspondences— list of (sub_rank_idx, mosaic_rank, residual_um)
      sub_dets       — detected sub-image centroids (stage coords)
      anchor_rank    — mosaic rank used as primary anchor
    """
    print(f"[corr]  Loading sub-image 10X: {Path(nd2_sub10x).name}")
    with _nd2lib.ND2File(nd2_sub10x) as f:
        pum_s  = _pixel_um(f)
        cxs, cys, _ = _stage_center(f)
        img_s  = _load_mip(f, ch_idx)
    hs, ws = img_s.shape

    # Load mosaic screener ranks
    with open(screener_csv, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    mosaic_pts = np.array([[float(r["x_um"]), float(r["y_um"])] for r in rows])
    mosaic_ranks = [int(r["rank"]) for r in rows]

    # Detect centroids in sub-image, convert to stage coords
    sub_dets = _detect_sub_centroids(img_s, pum_s)
    if not sub_dets:
        print("[corr]  No spheroids detected in sub-image — correction skipped.")
        return {}

    sub_stage = np.array([
        [cxs + (d["cx_px"] - ws / 2) * pum_s,
         cys + (d["cy_px"] - hs / 2) * pum_s]
        for d in sub_dets
    ])

    print(f"[corr]  Sub-image: {len(sub_dets)} detection(s)  "
          f"Mosaic: {len(mosaic_pts)} rank(s)")
    for i, (sx, sy) in enumerate(sub_stage):
        print(f"  S{i+1}: stage=({sx:.1f},{sy:.1f})  diam={sub_dets[i]['diam']:.0f}um")

    # ── exhaustive anchor search (pure translation consistency) ──────────────
    best_anchor = None
    for si_idx, si in enumerate(sub_stage):
        for mj_idx, mj in enumerate(mosaic_pts):
            Tx = mj[0] - si[0];  Ty = mj[1] - si[1]
            # apply translation to ALL sub-image points → distance to nearest rank
            shifted = sub_stage + np.array([Tx, Ty])
            rms = 0.0
            for sk in shifted:
                d_all = np.linalg.norm(mosaic_pts - sk, axis=1)
                rms  += float(d_all.min() ** 2)
            rms = math.sqrt(rms / len(shifted))
            if best_anchor is None or rms < best_anchor["rms"]:
                best_anchor = {"si": si_idx, "mj": mj_idx, "Tx": Tx, "Ty": Ty,
                               "rms": rms, "rank": mosaic_ranks[mj_idx]}

    Tx0 = best_anchor["Tx"];  Ty0 = best_anchor["Ty"]
    shifted0 = sub_stage + np.array([Tx0, Ty0])
    print(f"\n[corr]  Best anchor: S{best_anchor['si']+1} -> Rank {best_anchor['rank']}  "
          f"tx={Tx0:+.1f}  ty={Ty0:+.1f}  RMS={best_anchor['rms']:.1f}um")

    # ── collect all correspondences under this translation ───────────────────
    src_pts = []   # sub-image stage coords
    dst_pts = []   # mosaic stage coords
    corr_info = []
    used_mosaic = set()
    for si_idx, sk in enumerate(shifted0):
        d_all = np.linalg.norm(mosaic_pts - sk, axis=1)
        best_m = int(d_all.argmin())
        if best_m not in used_mosaic:
            used_mosaic.add(best_m)
            src_pts.append(sub_stage[si_idx])
            dst_pts.append(mosaic_pts[best_m])
            corr_info.append((si_idx + 1, mosaic_ranks[best_m], float(d_all[best_m])))

    src_arr = np.array(src_pts)
    dst_arr = np.array(dst_pts)

    # ── fit similarity transform ─────────────────────────────────────────────
    tform = SimilarityTransform()
    if len(src_arr) >= 2:
        tform.estimate(src=src_arr, dst=dst_arr)
        scale  = float(tform.scale)
        rot_d  = float(math.degrees(tform.rotation))
        tx, ty = float(tform.translation[0]), float(tform.translation[1])
    else:
        # only 1 point → pure translation
        tx, ty = Tx0, Ty0;  scale = 1.0;  rot_d = 0.0
        tform.params = np.array([[1, 0, tx], [0, 1, ty], [0, 0, 1]])

    # ── compute residuals with similarity transform ──────────────────────────
    if len(src_arr) >= 2:
        pred = tform(src_arr)          # shape (N, 2)
        residuals = np.linalg.norm(pred - dst_arr, axis=1).tolist()
    else:
        residuals = [0.0]

    print(f"[corr]  Similarity transform fit:")
    print(f"  tx={tx:+.1f}um  ty={ty:+.1f}um  scale={scale:.4f}  "
          f"rotation={rot_d:+.2f}deg")
    print(f"  Correspondences:")
    for (si, mr, _), res in zip(corr_info, residuals):
        print(f"    S{si} -> Rank {mr}  residual={res:.1f}um")

    return {
        "tform":           tform,
        "tx":              tx,
        "ty":              ty,
        "scale":           scale,
        "rotation_deg":    rot_d,
        "correspondences": corr_info,
        "residuals":       residuals,
        "sub_dets":        sub_dets,
        "sub_stage":       sub_stage,
        "anchor_rank":     best_anchor["rank"],
    }


def _apply_correction(sx: float, sy: float, corr: dict) -> tuple[float, float]:
    """Apply the similarity transform correction to a single stage point."""
    if not corr:
        return sx, sy
    tform = corr["tform"]
    pt    = tform(np.array([[sx, sy]]))[0]
    return float(pt[0]), float(pt[1])


# ── VALIDATE mode: confirm 20X images match screener coordinates ──────────────

_REPORT_COLS = [
    "nd2_20x", "stage_x_20x", "stage_y_20x",
    "col_predicted", "row_predicted",
    "scale_20x_to_10x", "crop_size_px",
    "ncc_peak", "confirmed",
    "shift_x_um", "shift_y_um",
    "closest_screener_rank", "dist_to_screener_um",
    # rigid-correction columns (populated when --sub-10x is supplied)
    "corrected_stage_x", "corrected_stage_y",
    "corrected_col", "corrected_row",
    "corrected_nearest_rank", "corrected_dist_um",
    "ncc_corrected", "confirmed_corrected",
    "corr_tx_um", "corr_ty_um", "corr_scale", "corr_rotation_deg",
]


def validate(nd2_10x: Path, nd2_20x_list: list[Path],
             ch_10x: int, ch_20x: int,
             out_dir: Path, screener_csv: Path | None = None,
             sub_10x: Path | None = None,
             gui: bool = False) -> tuple[list[dict], Path]:
    """
    For each 20X nd2:
      1. Predict its location in the 10X whole-well image from stage coords
      2. Extract matching 10X crop (same physical size as 20X FOV)
      3. Scale 20X MIP to 10X pixel space
      4. NCC between them → match quality
      5. Flag as confirmed / drifted / mis-matched

    Returns (report_records, report_csv_path).
    """
    print(f"[val]  Loading 10X: {nd2_10x.name}")
    with _nd2lib.ND2File(nd2_10x) as f:
        px_10x       = _pixel_um(f)
        cx0, cy0, _  = _stage_center(f)
        img_10x      = _load_mip(f, ch_10x)
    img_h, img_w = img_10x.shape
    print(f"[val]  Image {img_w}x{img_h} px  pixel={px_10x:.5f} um  "
          f"centre=({cx0:.1f},{cy0:.1f})")

    # Optional: load screener CSV to find closest detection per 20X file
    screener_list: list[dict] = []
    if screener_csv and screener_csv.exists():
        with open(screener_csv, newline="", encoding="utf-8") as fh:
            screener_list = list(csv.DictReader(fh))

    # Optional: fit 2D rigid correction from sub-image 10X
    corr: dict = {}
    if sub_10x and Path(sub_10x).exists() and screener_csv and screener_csv.exists():
        corr = _fit_stage_correction(Path(sub_10x), screener_csv, ch_10x)

    reports: list[dict] = []
    pairs:   list[tuple] = []   # for GUI

    for nd2_20 in nd2_20x_list:
        print(f"[val]  {nd2_20.name}")
        with _nd2lib.ND2File(nd2_20) as f:
            px_20x  = _pixel_um(f)
            sx, sy, _ = _stage_center(f)
            img_20x = _load_mip(f, ch_20x)
        fov_h_20, fov_w_20 = img_20x.shape

        # Scale factor: shrink the 20X image to 10X pixel space
        scale   = px_20x / px_10x                  # ~ 0.5 for a 2× zoom step
        crop_h  = int(round(fov_h_20 * scale))
        crop_w  = int(round(fov_w_20 * scale))

        # Predicted centre in the 10X whole-well image
        col_pred, row_pred = stage_to_pixel(sx, sy, cx0, cy0,
                                            px_10x, img_w, img_h)

        # Extract 10X crop at the predicted location
        crop_10x = _safe_crop(img_10x, row_pred, col_pred,
                               crop_h // 2, crop_w // 2)

        # Rescale 20X to 10X pixel space
        img_20x_scaled = _zoom(img_20x, scale, order=1)

        ncc_peak  = 0.0
        shift_x   = 0.0
        shift_y   = 0.0

        if crop_10x is not None and _SKI_OK:
            # Crop both to the same size (rounding in zoom/safe_crop can differ by 1 px)
            min_h = min(crop_10x.shape[0], img_20x_scaled.shape[0])
            min_w = min(crop_10x.shape[1], img_20x_scaled.shape[1])
            crop_10x      = crop_10x     [:min_h, :min_w]
            img_20x_scaled = img_20x_scaled[:min_h, :min_w]

            # Normalise both to [0,1] before NCC
            c_n = (crop_10x - crop_10x.min()) / ((crop_10x.max() - crop_10x.min()) + 1e-8)
            t_n = (img_20x_scaled - img_20x_scaled.min()) / ((img_20x_scaled.max() - img_20x_scaled.min()) + 1e-8)

            # Use the smaller image as template; search in the larger
            # Here both are the same size; pad_input=True centres the search
            ncc_map  = match_template(c_n, t_n, pad_input=True)
            pr, pc   = np.unravel_index(ncc_map.argmax(), ncc_map.shape)
            ncc_peak = float(ncc_map[pr, pc])

            # Shift = offset from centre of NCC map (which corresponds to zero shift)
            shift_r  = pr - ncc_map.shape[0] // 2
            shift_c  = pc - ncc_map.shape[1] // 2
            shift_x  = round(shift_c * px_10x, 2)
            shift_y  = round(shift_r * px_10x, 2)

        elif crop_10x is None:
            print(f"  WARNING: predicted location ({col_pred:.0f},{row_pred:.0f}) "
                  f"is out of the 10X image bounds — check stage coordinates")

        confirmed = ncc_peak >= MATCH_GOOD_NCC

        # Find the closest screener detection (raw coordinates)
        closest_rank = ""
        min_dist     = ""
        if screener_list:
            dists = [math.sqrt((float(r["x_um"]) - sx)**2 +
                               (float(r["y_um"]) - sy)**2)
                     for r in screener_list]
            idx = int(np.argmin(dists))
            closest_rank = screener_list[idx]["rank"]
            min_dist     = round(dists[idx], 1)

        status = "CONFIRMED" if confirmed else f"MISMATCH (NCC={ncc_peak:.3f})"
        print(f"  stage=({sx:.1f},{sy:.1f})  pred_px=({col_pred:.0f},{row_pred:.0f})  "
              f"NCC={ncc_peak:.3f}  shift=({shift_x:+.1f},{shift_y:+.1f}) um  {status}")

        # ── rigid-correction block ────────────────────────────────────────────
        corr_sx, corr_sy   = _apply_correction(sx, sy, corr)
        corr_col, corr_row = stage_to_pixel(corr_sx, corr_sy, cx0, cy0,
                                             px_10x, img_w, img_h)
        corr_rank  = "";  corr_dist = "";  ncc_corr = 0.0;  conf_corr = False

        if corr and screener_list:
            c_dists = [math.sqrt((float(r["x_um"]) - corr_sx)**2 +
                                  (float(r["y_um"]) - corr_sy)**2)
                       for r in screener_list]
            ci = int(np.argmin(c_dists))
            corr_rank = screener_list[ci]["rank"]
            corr_dist = round(c_dists[ci], 1)

            # NCC at corrected location
            corr_crop = _safe_crop(img_10x, corr_row, corr_col,
                                   crop_h // 2, crop_w // 2)
            if corr_crop is not None and _SKI_OK:
                min_h2 = min(corr_crop.shape[0], img_20x_scaled.shape[0])
                min_w2 = min(corr_crop.shape[1], img_20x_scaled.shape[1])
                cc = corr_crop[:min_h2, :min_w2]
                tt = img_20x_scaled[:min_h2, :min_w2]
                cc_n = (cc - cc.min()) / ((cc.max() - cc.min()) + 1e-8)
                tt_n = (tt - tt.min()) / ((tt.max() - tt.min()) + 1e-8)
                nm   = match_template(cc_n, tt_n, pad_input=True)
                ncc_corr  = float(nm.max())
                conf_corr = ncc_corr >= MATCH_GOOD_NCC

            corr_status = "CONFIRMED" if conf_corr else f"MISMATCH (NCC={ncc_corr:.3f})"
            delta_x = corr_sx - float(screener_list[ci]["x_um"])
            delta_y = corr_sy - float(screener_list[ci]["y_um"])
            print(f"  [corrected] stage=({corr_sx:.1f},{corr_sy:.1f})  "
                  f"px=({corr_col:.0f},{corr_row:.0f})  "
                  f"-> Rank {corr_rank}  dist={corr_dist}um  "
                  f"NCC={ncc_corr:.3f}  {corr_status}")
            print(f"  [xy-error]  dx={delta_x:+.1f}um  dy={delta_y:+.1f}um  "
                  f"|err|={math.sqrt(delta_x**2+delta_y**2):.1f}um")

        r_dict = {
            "nd2_20x":               nd2_20.name,
            "stage_x_20x":           round(sx, 2),
            "stage_y_20x":           round(sy, 2),
            "col_predicted":         round(col_pred, 1),
            "row_predicted":         round(row_pred, 1),
            "scale_20x_to_10x":      round(scale, 5),
            "crop_size_px":          f"{crop_w}×{crop_h}",
            "ncc_peak":              round(ncc_peak, 4),
            "confirmed":             confirmed,
            "shift_x_um":            shift_x,
            "shift_y_um":            shift_y,
            "closest_screener_rank": closest_rank,
            "dist_to_screener_um":   min_dist,
            # rigid-correction
            "corrected_stage_x":      round(corr_sx, 2) if corr else "",
            "corrected_stage_y":      round(corr_sy, 2) if corr else "",
            "corrected_col":          round(corr_col, 1) if corr else "",
            "corrected_row":          round(corr_row, 1) if corr else "",
            "corrected_nearest_rank": corr_rank,
            "corrected_dist_um":      corr_dist,
            "ncc_corrected":          round(ncc_corr, 4) if corr else "",
            "confirmed_corrected":    conf_corr if corr else "",
            "corr_tx_um":             round(corr["tx"], 2) if corr else "",
            "corr_ty_um":             round(corr["ty"], 2) if corr else "",
            "corr_scale":             round(corr["scale"], 5) if corr else "",
            "corr_rotation_deg":      round(corr["rotation_deg"], 3) if corr else "",
        }
        reports.append(r_dict)
        pairs.append((nd2_20.stem, crop_10x, img_20x_scaled, ncc_peak, confirmed))

    if gui and _MPL_OK:
        _gui_validate(img_10x, px_10x, cx0, cy0, reports, pairs, screener_list)

    out_dir.mkdir(parents=True, exist_ok=True)
    out_csv = out_dir / "validation_report.csv"
    _write_csv(reports, out_csv, _REPORT_COLS)
    n_conf = sum(1 for r in reports if r["confirmed"])
    print(f"[val]  {n_conf}/{len(reports)} confirmed.  Wrote: {out_csv}")
    return reports, out_csv


# ── GUI helpers ───────────────────────────────────────────────────────────────

def _gui_review_pre(img_10x: np.ndarray, px_um: float,
                    records: list[dict]) -> list[dict]:
    """
    Interactive pre-capture review.
    Left: whole-well 10X overview with detection circles.
    Right: zoom crop around the currently selected spheroid.
    Keys: a=Accept  r=Reject  n/→=Next  p/←=Prev  q=Quit
    Left-click on left panel to override centroid.
    """
    recs  = [dict(r) for r in records]
    state = {"idx": 0}

    fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(17, 8))
    fig.suptitle(
        "Pre-capture Verification  —  a=Accept  r=Reject  n=Next  p=Prev  q=Quit\n"
        "Left-click on overview to relocate current spheroid centroid",
        fontsize=9)

    def draw():
        idx = state["idx"]
        r   = recs[idx]

        # ── overview (left) ───────────────────────────────────────────────
        ax_l.cla()
        vlo = float(np.percentile(img_10x, 0.5))
        vhi = float(np.percentile(img_10x, 99.5))
        ax_l.imshow(img_10x, cmap="gray", vmin=vlo, vmax=vhi,
                    origin="upper", interpolation="nearest")
        for i, rv in enumerate(recs):
            col  = rv["canvas_col"]; row = rv["canvas_row"]
            dpx  = rv["diameter_um"] / px_um / 2
            clr  = "lime" if rv["verified"] else "red"
            lw   = 3.0 if i == idx else 1.0
            ax_l.add_patch(mpatches.Circle((col, row), dpx,
                                           color=clr, fill=False, lw=lw))
            ax_l.text(col, row - dpx - 6, str(rv["rank"]),
                      color=clr, fontsize=6, ha="center")
        ax_l.set_title(f"10X whole-well  ({img_10x.shape[1]}×{img_10x.shape[0]} px)",
                       fontsize=8)
        ax_l.axis("off")

        # ── crop (right) ─────────────────────────────────────────────────
        ax_r.cla()
        col = r["canvas_col"]; row = r["canvas_row"]
        dpx = r["diameter_um"] / px_um
        pad = int(dpx * 1.5) + 30
        r0  = max(0, int(row - pad)); r1 = min(img_10x.shape[0], int(row + pad))
        c0  = max(0, int(col - pad)); c1 = min(img_10x.shape[1], int(col + pad))
        crop = img_10x[r0:r1, c0:c1]
        lo2  = float(np.percentile(crop, 0.5))
        hi2  = float(np.percentile(crop, 99.5))
        ax_r.imshow(crop, cmap="gray", vmin=lo2, vmax=hi2,
                    origin="upper", interpolation="nearest")
        ax_r.plot(col - c0, row - r0, "g+", ms=14, mew=2)
        ax_r.add_patch(mpatches.Circle((col-c0, row-r0), dpx/2,
                                       color="yellow", fill=False,
                                       lw=1.5, ls="--"))
        st = "✓ ACCEPTED" if r["verified"] else "✗ REJECTED"
        ax_r.set_title(
            f"Rank {r['rank']}  ∅={r['diameter_um']:.0f} µm  "
            f"NCC={r['ncc_peak']:.3f}  shift={r['shift_um']:.1f} µm\n"
            f"flags=[{r['qc_flags'] or 'OK'}]  {st}",
            fontsize=8)
        ax_r.axis("off")
        fig.canvas.draw_idle()

    def on_key(event):
        idx = state["idx"]
        if event.key == "a":
            recs[idx]["verified"] = True
        elif event.key == "r":
            recs[idx]["verified"] = False
        elif event.key in ("n", "right"):
            state["idx"] = min(idx + 1, len(recs) - 1)
        elif event.key in ("p", "left"):
            state["idx"] = max(idx - 1, 0)
        elif event.key == "q":
            plt.close(fig); return
        draw()

    def on_click(event):
        if event.inaxes is not ax_l or event.button != 1:
            return
        idx   = state["idx"]
        new_col, new_row = event.xdata, event.ydata
        recs[idx]["canvas_col"] = round(new_col, 1)
        recs[idx]["canvas_row"] = round(new_row, 1)
        fl = recs[idx]["qc_flags"]
        recs[idx]["qc_flags"] = (fl + "|manual" if fl else "manual")
        draw()

    fig.canvas.mpl_connect("key_press_event",    on_key)
    fig.canvas.mpl_connect("button_press_event", on_click)
    draw()
    plt.tight_layout()
    plt.show()
    return recs


def _gui_validate(img_10x: np.ndarray, px_10x: float,
                  cx0: float, cy0: float,
                  reports: list[dict],
                  pairs: list[tuple],
                  screener_list: list[dict]) -> None:
    """
    Validation overview:
    - Left panel: whole-well 10X with predicted 20X FOV boxes and screener dots
    - Right panel: side-by-side 10X crop vs. scaled 20X for current file
    Keys: n/p navigate files; q quits.
    """
    img_h, img_w = img_10x.shape
    state = {"idx": 0}

    fig, axes = plt.subplots(1, 3, figsize=(20, 7))
    ax_ov, ax_10c, ax_20c = axes
    fig.suptitle("Post-Capture Validation  —  n=Next  p=Prev  q=Quit", fontsize=9)

    def draw():
        idx = state["idx"]
        rpt = reports[idx]
        nm, crop_10x, img_20x_sc, ncc, confirmed = pairs[idx]

        # ── overview ─────────────────────────────────────────────────────
        ax_ov.cla()
        lo = float(np.percentile(img_10x, 0.5))
        hi = float(np.percentile(img_10x, 99.5))
        ax_ov.imshow(img_10x, cmap="gray", vmin=lo, vmax=hi,
                     origin="upper", interpolation="nearest")

        # Screener dots
        for sr in screener_list:
            sc, sr_row = stage_to_pixel(float(sr["x_um"]), float(sr["y_um"]),
                                        cx0, cy0, px_10x, img_w, img_h)
            ax_ov.plot(sc, sr_row, "c.", ms=4)

        # 20X predicted FOV boxes
        for i, rp in enumerate(reports):
            col_p = rp["col_predicted"]; row_p = rp["row_predicted"]
            crop_sz = rp["crop_size_px"].split("×")
            cw, ch = int(crop_sz[0]), int(crop_sz[1])
            ec  = "lime" if rp["confirmed"] else "red"
            lw  = 2.5 if i == idx else 1.0
            ax_ov.add_patch(mpatches.Rectangle(
                (col_p - cw/2, row_p - ch/2), cw, ch,
                ec=ec, fill=False, lw=lw))
            ax_ov.text(col_p, row_p - ch/2 - 5,
                       f"20X-{i+1}", color=ec, fontsize=6, ha="center")

        ax_ov.set_title("Overview (cyan=screener, box=20X FOV)", fontsize=8)
        ax_ov.axis("off")

        # ── 10X crop (centre) ─────────────────────────────────────────────
        ax_10c.cla()
        if crop_10x is not None:
            lo2 = float(np.percentile(crop_10x, 0.5))
            hi2 = float(np.percentile(crop_10x, 99.5))
            ax_10c.imshow(crop_10x, cmap="gray", vmin=lo2, vmax=hi2,
                          origin="upper", interpolation="nearest")
        else:
            ax_10c.text(0.5, 0.5, "Out of bounds", transform=ax_10c.transAxes,
                        ha="center", va="center", color="red")
        ax_10c.set_title(f"10X crop at predicted location\n"
                         f"({rpt['col_predicted']:.0f},{rpt['row_predicted']:.0f}) px",
                         fontsize=8)
        ax_10c.axis("off")

        # ── 20X scaled (right) ────────────────────────────────────────────
        ax_20c.cla()
        if img_20x_sc is not None:
            lo3 = float(np.percentile(img_20x_sc, 0.5))
            hi3 = float(np.percentile(img_20x_sc, 99.5))
            ax_20c.imshow(img_20x_sc, cmap="gray", vmin=lo3, vmax=hi3,
                          origin="upper", interpolation="nearest")
        st    = "✓ CONFIRMED" if confirmed else "✗ MISMATCH"
        color = "green" if confirmed else "red"
        ax_20c.set_title(
            f"20X MIP (scaled to 10X pixels)\n"
            f"NCC={ncc:.3f}  {st}  shift=({rpt['shift_x_um']:+.1f},{rpt['shift_y_um']:+.1f}) µm",
            fontsize=8, color=color)
        ax_20c.axis("off")

        fig.canvas.draw_idle()

    def on_key(event):
        idx = state["idx"]
        if event.key in ("n", "right"):
            state["idx"] = min(idx + 1, len(reports) - 1)
        elif event.key in ("p", "left"):
            state["idx"] = max(idx - 1, 0)
        elif event.key == "q":
            plt.close(fig); return
        draw()

    fig.canvas.mpl_connect("key_press_event", on_key)
    draw()
    plt.tight_layout()
    plt.show()


# ── CSV helper ────────────────────────────────────────────────────────────────

def _write_csv(records: list[dict], out_path: Path, cols: list[str]) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(records)


# ── Overlay renderer ─────────────────────────────────────────────────────────

def save_overlay(nd2_10x: Path, screen_csv: Path, out_path: Path,
                 ch_idx: int = 0, dpi: int = 150) -> Path:
    """
    Render the 10X whole-well image with a red-circle segmentation overlay.

    For every row in screen_csv:
      • Red circle whose radius = detected diameter / 2 (in image pixels)
      • Rank number centred inside the circle
      • Score shown below the rank in smaller text

    The output is a publication-ready PNG saved to out_path.
    Returns out_path.
    """
    if not _MPL_OK:
        raise RuntimeError("matplotlib is required for overlay rendering.  "
                           "pip install matplotlib")

    print(f"[overlay]  Loading 10X image: {nd2_10x.name}")
    with _nd2lib.ND2File(nd2_10x) as f:
        px_um       = _pixel_um(f)
        cx0, cy0, _ = _stage_center(f)
        img         = _load_mip(f, ch_idx)
    img_h, img_w = img.shape
    print(f"[overlay]  {img_w}x{img_h} px  pixel={px_um:.5f} um")

    with open(screen_csv, newline="", encoding="utf-8") as fh:
        records = list(csv.DictReader(fh))
    print(f"[overlay]  Annotating {len(records)} spheroid(s) from {screen_csv.name}")

    # ── figure ─────────────────────────────────────────────────────────────
    # Fix the figure at 16×16 inches; dpi controls the saved PNG resolution.
    fig, ax = plt.subplots(1, 1, figsize=(16, 16), facecolor="black")
    ax.set_facecolor("black")

    lo = float(np.percentile(img, 0.5))
    hi = float(np.percentile(img, 99.8))
    ax.imshow(img, cmap="gray", vmin=lo, vmax=hi,
              origin="upper", interpolation="nearest",
              extent=[0, img_w, img_h, 0])   # keeps data coords = pixel coords

    # Font size scales with the figure→image ratio so labels stay legible
    # regardless of the total image resolution.
    # A 200 µm spheroid ≈ 310 px; we want rank text ≈ 20% of that.
    ref_diam_px  = 200.0 / px_um          # reference circle size (200 µm in px)
    label_pts    = max(6.0, ref_diam_px * 16.0 / img_w * 0.25 * 72)
    sub_pts      = label_pts * 0.65

    for rec in records:
        x_um = float(rec.get("x_um",     rec.get("x_um_raw",    0)))
        y_um = float(rec.get("y_um",     rec.get("y_um_raw",    0)))
        diam = float(rec["diameter_um"])
        rank = int  (rec["rank"])
        score = float(rec.get("score", 0))

        col, row   = stage_to_pixel(x_um, y_um, cx0, cy0, px_um, img_w, img_h)
        radius_px  = diam / 2.0 / px_um

        # ── red detection circle ──────────────────────────────────────────
        circle = mpatches.Circle(
            (col, row), radius_px,
            edgecolor="red", facecolor="none",
            linewidth=max(1.5, radius_px * 0.025),   # line ~2.5% of radius
            zorder=3)
        ax.add_patch(circle)

        # ── rank label (inside circle, centred) ──────────────────────────
        ax.text(col, row, str(rank),
                color="red", fontsize=label_pts, fontweight="bold",
                ha="center", va="center", zorder=4,
                bbox=dict(boxstyle="round,pad=0.15",
                          facecolor="black", alpha=0.55, edgecolor="none"))

        # ── score below the circle ────────────────────────────────────────
        ax.text(col, row + radius_px + label_pts * 0.6,
                f"score={score:.3f}",
                color=(1.0, 0.5, 0.5),   # light red
                fontsize=sub_pts, ha="center", va="top", zorder=4,
                bbox=dict(boxstyle="round,pad=0.1",
                          facecolor="black", alpha=0.45, edgecolor="none"))

    # ── title bar ──────────────────────────────────────────────────────────
    n  = len(records)
    ax.set_title(
        f"10X Whole-Well Segmentation — {n} spheroid{'s' if n != 1 else ''} detected"
        f"  |  {nd2_10x.stem}",
        fontsize=10, color="white", pad=6)

    ax.set_xlim(0, img_w)
    ax.set_ylim(img_h, 0)    # origin upper
    ax.axis("off")

    # ── save ───────────────────────────────────────────────────────────────
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight",
                facecolor="black", edgecolor="none")
    plt.close(fig)

    mb = out_path.stat().st_size / 1e6
    print(f"[overlay]  Saved: {out_path}  ({mb:.1f} MB)")
    return out_path


# ── CLI ───────────────────────────────────────────────────────────────────────

def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(
        description="Cross-zoom coordinate verification for 10X stitched "
                    "whole-well nd2 + 20X spheroid nd2 files")
    ap.add_argument("--mode",      choices=["pre", "validate", "overlay"],
                    default="pre",
                    help="pre=verify before 20X capture  "
                         "validate=confirm after  "
                         "overlay=save annotated PNG of whole-well image")
    ap.add_argument("--nd2-10x",   type=Path, required=True,
                    help="stitched whole-well 10X nd2 (P=1, e.g. WellD05...nd2)")
    ap.add_argument("--screen",    type=Path, default=None,
                    help="screener CSV (spheroid_screen_latest.csv or "
                         "verified_screen.csv)")
    ap.add_argument("--nd2-20x",   type=Path, default=None,
                    help="[validate] folder of 20X nd2 files, or single nd2 path")
    ap.add_argument("--sub-10x",   type=Path, default=None,
                    help="[validate] single-position 10X sub-image taken at the same "
                         "location as the 20X capture; used to fit a 2D rigid "
                         "correction for stage-coordinate offset between sessions")
    ap.add_argument("--out-dir",   type=Path, default=None,
                    help="output directory (default: same folder as nd2-10x)")
    ap.add_argument("--out-png",   type=Path, default=None,
                    help="[overlay] explicit output PNG path")
    ap.add_argument("--dpi",       type=int,  default=150,
                    help="[overlay] output PNG DPI (default: 150)")
    ap.add_argument("--ch-10x",    type=int,  default=0, help="10X channel index")
    ap.add_argument("--ch-20x",    type=int,  default=0, help="20X channel index")
    ap.add_argument("--search-um", type=float, default=NCC_SEARCH_UM,
                    help=f"NCC search radius um for pre-mode refinement "
                         f"(default {NCC_SEARCH_UM})")
    ap.add_argument("--gui",       action="store_true",
                    help="open interactive Matplotlib review window")
    args = ap.parse_args(argv)

    if _MISSING:
        print("Missing packages:  " + "  ".join(_MISSING))
        print("Install:  pip install nd2 numpy scipy scikit-image matplotlib")
        sys.exit(1)

    out_dir = args.out_dir or args.nd2_10x.parent

    if args.mode == "pre":
        screen_csv = (args.screen
                      or args.nd2_10x.parent / "spheroid_screen_latest.csv")
        if not screen_csv.exists():
            print(f"Screener CSV not found: {screen_csv}\n"
                  "Run spheroid_screener.py first, or pass --screen explicitly.")
            sys.exit(1)
        pre_verify(args.nd2_10x, screen_csv,
                   args.ch_10x, args.search_um, out_dir, gui=args.gui)

    elif args.mode == "validate":
        if args.nd2_20x is None:
            print("--nd2-20x required for validate mode")
            sys.exit(1)
        p20 = args.nd2_20x
        nd2_20x_list = sorted(p20.glob("*.nd2")) if p20.is_dir() else [p20]
        if not nd2_20x_list:
            print(f"No .nd2 files in {p20}")
            sys.exit(1)
        validate(args.nd2_10x, nd2_20x_list,
                 args.ch_10x, args.ch_20x,
                 out_dir,
                 screener_csv=args.screen,
                 sub_10x=args.sub_10x,
                 gui=args.gui)

    else:   # overlay
        screen_csv = (args.screen
                      or out_dir / "spheroid_screen_latest.csv"
                      or out_dir / "verified_screen.csv")
        if not screen_csv.exists():
            # try verified_screen.csv as fallback
            fb = out_dir / "verified_screen.csv"
            if fb.exists():
                screen_csv = fb
            else:
                print(f"No screener CSV found at {screen_csv}\n"
                      "Run spheroid_screener.py first, then pass --screen.")
                sys.exit(1)
        stem    = args.nd2_10x.stem
        out_png = (args.out_png
                   or out_dir / f"{stem}_segmentation_overlay.png")
        save_overlay(args.nd2_10x, screen_csv, out_png,
                     ch_idx=args.ch_10x, dpi=args.dpi)


if __name__ == "__main__":
    main()
