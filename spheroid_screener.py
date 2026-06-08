"""
spheroid_screener.py
10X Spheroid Auto-Screen — Detect, Rank, and Export for NIS-E 20X Follow-Up

Pipeline:
  1. Find the latest .nd2 in a NIS-E output folder (or watch for new arrivals)
  2. Detect all spheroids across every FOV position in the 10X scan
  3. Score each detection by circularity (4πA/P²) and solidity (A/A_convex)
  4. Rank globally; export a stage-coordinate CSV for the NIS-E 20X macro

Usage:
  # Process the single latest .nd2 in a folder, then exit
  python spheroid_screener.py --watch-dir D:\\NIS_Output

  # Process one specific file
  python spheroid_screener.py --nd2 D:\\NIS_Output\\WellD05_Channel555.nd2

  # Watch for new files continuously (daemon mode; requires watchdog)
  python spheroid_screener.py --watch-dir D:\\NIS_Output --daemon

Output CSV columns:
  rank, position_idx, x_um, y_um, z_surface_um,
  diameter_um, circularity, solidity, area_px, score, cx_px, cy_px, pixel_um

Dependencies:
  pip install nd2 numpy scipy scikit-image
  pip install watchdog          # only needed for --daemon mode
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
import time
from datetime import datetime
from pathlib import Path

# ── optional scientific stack ─────────────────────────────────────────────────

_MISSING: list[str] = []

try:
    import nd2 as _nd2lib
except ImportError:
    _nd2lib = None
    _MISSING.append("nd2")

try:
    import numpy as np
except ImportError:
    np = None                    # type: ignore[assignment]
    _MISSING.append("numpy")

try:
    from scipy import ndimage as _ndi
except ImportError:
    _ndi = None
    _MISSING.append("scipy")

try:
    from skimage.measure import regionprops, label as _ski_label
    from skimage.feature import peak_local_max as _peak_local_max
    from skimage.segmentation import watershed as _watershed
    _SKI_OK = True
    _WATERSHED_OK = True
except ImportError:
    try:
        from skimage.measure import regionprops, label as _ski_label
        _SKI_OK = True
    except ImportError:
        _SKI_OK = False
        _MISSING.append("scikit-image  →  pip install scikit-image")
    _WATERSHED_OK = False

try:
    from cellpose import models as _cp_models
    _CELLPOSE_OK = True
except ImportError:
    _CELLPOSE_OK = False

try:
    from watchdog.observers import Observer
    from watchdog.events import FileSystemEventHandler as _FSH
    _WATCHDOG_OK = True
except ImportError:
    _WATCHDOG_OK = False
    _FSH = object                # type: ignore[assignment,misc]


# ── tuneable defaults ─────────────────────────────────────────────────────────

DEFAULT_MIN_DIAM_UM  = 50.0
DEFAULT_MAX_DIAM_UM  = 400.0
DEFAULT_DETECT_SIGMA = 2.0    # Gaussian pre-blur (px)
DEFAULT_DETECT_K     = 1.5    # threshold = mean + K * std
DEFAULT_W_CIRC       = 0.6    # circularity weight in composite score
DEFAULT_W_SOL        = 0.4    # solidity weight in composite score
DEFAULT_CH_IDX       = 0      # 0-based channel index for detection
DEFAULT_BACKEND      = "threshold"


# ── ND2 helpers ───────────────────────────────────────────────────────────────

def _extract_pixel_um(f) -> float:
    try:
        v = float(f.voxel_size().x)
        return v if v and v > 0 else 0.65
    except Exception:
        return 0.65   # 10x fallback


def _stage_xy_z(f, seq_idx: int) -> tuple[float, float, float]:
    try:
        fm  = f.frame_metadata(seq_idx)
        pos = fm.channels[0].position.stagePositionUm
        return float(pos.x), float(pos.y), float(pos.z)
    except Exception:
        return 0.0, 0.0, 0.0


def _max_projection(f, pos_idx: int, ch_idx: int) -> np.ndarray:
    """Return a 2-D float32 max-intensity projection for (position, channel)."""
    arr  = f.asarray()
    dims = list(f.sizes.keys())

    def _ax(k: str) -> int | None:
        return dims.index(k) if k in dims else None

    p_ax = _ax("P"); z_ax = _ax("Z"); c_ax = _ax("C")
    img  = arr

    if p_ax is not None:
        sl = [slice(None)] * img.ndim; sl[p_ax] = pos_idx
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


# ── morphological detection & ranking ────────────────────────────────────────

def _circularity(area: float, perimeter: float) -> float:
    """4π·A / P² — equals 1.0 for a perfect circle, <1 for non-circular shapes."""
    return min(4.0 * math.pi * area / (perimeter ** 2), 1.0) if perimeter > 1 else 0.0


def _watershed_label(binary: np.ndarray, pixel_um: float,
                     min_diam_um: float) -> np.ndarray:
    """
    Label binary mask using distance-transform watershed.
    Peaks are separated by at least half a minimum-diameter in pixels,
    so a merged pair that is ~2× the min diameter will be correctly split.
    """
    dist = _ndi.distance_transform_edt(binary)
    min_dist_px = max(4, int(min_diam_um * 0.5 / pixel_um))
    coords = _peak_local_max(dist, min_distance=min_dist_px, labels=binary)
    peak_mask = np.zeros(dist.shape, dtype=bool)
    if coords.size:
        peak_mask[tuple(coords.T)] = True
    markers, _ = _ndi.label(peak_mask)
    return _watershed(-dist, markers, mask=binary)


def _cellpose_detect(image: np.ndarray, pixel_um: float,
                     min_diam_um: float, max_diam_um: float) -> np.ndarray:
    """
    Run Cellpose cyto3 model and return a labeled mask (int32 array).
    Diameter hint is the midpoint of the allowed range.
    """
    typical_diam_px = ((min_diam_um + max_diam_um) / 2.0) / pixel_um
    model = _cp_models.Cellpose(model_type="cyto3", gpu=False)
    img_f = image.astype(np.float32)
    lo, hi = img_f.min(), img_f.max()
    if hi > lo:
        img_f = (img_f - lo) / (hi - lo) * 255.0
    img_u8 = img_f.astype(np.uint8)
    masks, _, _, _ = model.eval(img_u8, diameter=typical_diam_px, channels=[0, 0])
    return masks.astype(np.int32)


def _measure_labeled(labeled: np.ndarray, min_area: float, max_area: float,
                     pixel_um: float, w_circ: float, w_sol: float) -> list[dict]:
    """Run regionprops on a labeled array and build detection records."""
    results = []
    for rp in regionprops(labeled):
        area = rp.area
        if not (min_area <= area <= max_area):
            continue
        circ  = _circularity(area, rp.perimeter)
        sol   = float(rp.solidity)
        score = w_circ * circ + w_sol * sol
        cy, cx = rp.centroid
        diam   = 2.0 * math.sqrt(area / math.pi) * pixel_um
        results.append({
            "cx_px":        cx,
            "cy_px":        cy,
            "diameter_um":  round(diam, 2),
            "area_px":      int(area),
            "circularity":  round(circ, 4),
            "solidity":     round(sol, 4),
            "perimeter_px": round(rp.perimeter, 1),
            "score":        round(score, 4),
        })
    return results


def detect_and_rank(image: np.ndarray, pixel_um: float,
                    min_diam_um: float = DEFAULT_MIN_DIAM_UM,
                    max_diam_um: float = DEFAULT_MAX_DIAM_UM,
                    sigma:       float = DEFAULT_DETECT_SIGMA,
                    thresh_k:    float = DEFAULT_DETECT_K,
                    w_circ:      float = DEFAULT_W_CIRC,
                    w_sol:       float = DEFAULT_W_SOL,
                    backend:     str   = DEFAULT_BACKEND) -> list[dict]:
    """
    Detect spheroids in a 2-D image and rank by composite quality score.

    Parameters
    ----------
    backend : "threshold" | "cellpose"
        "threshold" — Gaussian-threshold + distance-transform watershed (default).
                      Watershed automatically splits merged blobs that are wider
                      than min_diam_um.
        "cellpose"  — Cellpose cyto3 model (requires: pip install cellpose).
                      Falls back to "threshold" if cellpose is not installed.

    Returns list of dicts (sorted best-first):
      cx_px, cy_px, diameter_um, area_px,
      circularity, solidity, perimeter_px, score
    """
    min_area = math.pi * (min_diam_um / pixel_um / 2) ** 2
    max_area = math.pi * (max_diam_um / pixel_um / 2) ** 2
    results:  list[dict] = []

    # ── cellpose path ─────────────────────────────────────────────────────────
    if backend == "cellpose":
        if not _CELLPOSE_OK:
            print("[WARN]  cellpose not installed — pip install cellpose. "
                  "Falling back to threshold backend.")
            backend = "threshold"
        elif not _SKI_OK:
            print("[WARN]  scikit-image required for regionprops with cellpose. "
                  "Falling back to threshold backend.")
            backend = "threshold"
        else:
            labeled = _cellpose_detect(image, pixel_um, min_diam_um, max_diam_um)
            results = _measure_labeled(labeled, min_area, max_area,
                                       pixel_um, w_circ, w_sol)
            results.sort(key=lambda r: r["score"], reverse=True)
            return results

    # ── threshold + watershed path ────────────────────────────────────────────
    blurred = _ndi.gaussian_filter(image, sigma=sigma)
    mu, std  = float(blurred.mean()), float(blurred.std())
    binary   = blurred > (mu + thresh_k * std)
    filled   = _ndi.binary_fill_holes(binary)

    if _SKI_OK:
        if _WATERSHED_OK:
            labeled = _watershed_label(filled, pixel_um, min_diam_um)
        else:
            labeled = _ski_label(filled)
        results = _measure_labeled(labeled, min_area, max_area,
                                   pixel_um, w_circ, w_sol)
    else:
        # scipy-only fallback: area + centroid, no shape metrics
        labeled_sci, n_obj = _ndi.label(filled)
        for lbl in range(1, n_obj + 1):
            mask = labeled_sci == lbl
            area = int(mask.sum())
            if not (min_area <= area <= max_area):
                continue
            cy, cx = _ndi.center_of_mass(mask)
            diam   = 2.0 * math.sqrt(area / math.pi) * pixel_um
            results.append({
                "cx_px":        cx,
                "cy_px":        cy,
                "diameter_um":  round(diam, 2),
                "area_px":      area,
                "circularity":  0.0,
                "solidity":     0.0,
                "perimeter_px": 0.0,
                "score":        0.0,
            })

    results.sort(key=lambda r: r["score"], reverse=True)
    return results


# ── ND2 processing ────────────────────────────────────────────────────────────

def process_nd2(nd2_path: Path, cfg: dict) -> list[dict]:
    """
    Load a 10X ND2, detect + rank all spheroids across all positions.
    Returns records sorted best-first, each with an overall 'rank' field.
    """
    ch_idx   = cfg.get("ch_idx",       DEFAULT_CH_IDX)
    min_d    = cfg.get("min_diam_um",  DEFAULT_MIN_DIAM_UM)
    max_d    = cfg.get("max_diam_um",  DEFAULT_MAX_DIAM_UM)
    sigma    = cfg.get("sigma",        DEFAULT_DETECT_SIGMA)
    thresh_k = cfg.get("thresh_k",     DEFAULT_DETECT_K)
    w_circ   = cfg.get("w_circ",       DEFAULT_W_CIRC)
    w_sol    = cfg.get("w_sol",        DEFAULT_W_SOL)
    backend  = cfg.get("backend",      DEFAULT_BACKEND)

    all_records: list[dict] = []

    with _nd2lib.ND2File(nd2_path) as f:
        pixel_um = _extract_pixel_um(f)
        sizes    = dict(f.sizes)
        n_pos    = sizes.get("P", 1)
        n_z      = sizes.get("Z", 1)
        n_ch     = sizes.get("C", 1)

        for p in range(n_pos):
            seq_idx       = p * n_z * n_ch
            sx, sy, sz    = _stage_xy_z(f, seq_idx)
            img           = _max_projection(f, p, ch_idx)
            img_h, img_w  = img.shape

            detections = detect_and_rank(
                img, pixel_um, min_d, max_d, sigma, thresh_k, w_circ, w_sol,
                backend=backend)

            for obj in detections:
                # Convert pixel centroid → absolute stage coordinate.
                # sx/sy is the stage position of the FOV centre.
                dx = (obj["cx_px"] - img_w / 2.0) * pixel_um
                dy = (obj["cy_px"] - img_h / 2.0) * pixel_um
                all_records.append({
                    "position_idx":  p,
                    "x_um":          round(sx + dx, 2),
                    "y_um":          round(sy + dy, 2),
                    "z_surface_um":  round(sz, 2),
                    "diameter_um":   obj["diameter_um"],
                    "circularity":   obj["circularity"],
                    "solidity":      obj["solidity"],
                    "area_px":       obj["area_px"],
                    "score":         obj["score"],
                    "cx_px":         round(obj["cx_px"], 1),
                    "cy_px":         round(obj["cy_px"], 1),
                    "pixel_um":      pixel_um,
                })

    # Global re-rank by score (best across all positions first)
    all_records.sort(key=lambda r: r["score"], reverse=True)
    for i, r in enumerate(all_records):
        r["rank"] = i + 1

    return all_records


# ── file discovery ────────────────────────────────────────────────────────────

def latest_nd2(watch_dir: Path) -> Path | None:
    """Return the most recently modified .nd2 file in watch_dir."""
    nd2s = sorted(watch_dir.glob("*.nd2"), key=lambda p: p.stat().st_mtime)
    return nd2s[-1] if nd2s else None


# ── CSV export ────────────────────────────────────────────────────────────────

_CSV_COLS = [
    "rank", "position_idx",
    "x_um", "y_um", "z_surface_um",
    "diameter_um", "circularity", "solidity", "area_px", "score",
    "cx_px", "cy_px", "pixel_um",
]


def export_csv(records: list[dict], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=_CSV_COLS, extrasaction="ignore")
        w.writeheader()
        w.writerows(records)


# ── one-shot runner ───────────────────────────────────────────────────────────

def run_once(nd2_path: Path, out_dir: Path, cfg: dict) -> tuple[list[dict], Path]:
    ts = datetime.now().strftime("%H:%M:%S")
    print(f"[{ts}]  Screening: {nd2_path.name}")

    records = process_nd2(nd2_path, cfg)

    # Per-run CSV (timestamped stem) + a fixed "latest" path for the NIS macro
    out_csv    = out_dir / f"{nd2_path.stem}_screen.csv"
    latest_csv = out_dir / "spheroid_screen_latest.csv"
    export_csv(records, out_csv)
    export_csv(records, latest_csv)

    ts2 = datetime.now().strftime("%H:%M:%S")
    print(f"[{ts2}]  {len(records)} spheroid(s) detected and ranked:")
    for r in records[:10]:
        print(f"  Rank {r['rank']:2d}  "
              f"XY=({r['x_um']:9.1f},{r['y_um']:9.1f})  "
              f"Z={r['z_surface_um']:7.1f}  "
              f"diam={r['diameter_um']:6.1f}um  "
              f"circ={r['circularity']:.3f}  "
              f"sol={r['solidity']:.3f}  "
              f"score={r['score']:.3f}")
    if len(records) > 10:
        print(f"  ... (+{len(records)-10} more, see CSV)")
    print(f"[{ts2}]  Wrote: {out_csv}")
    print(f"[{ts2}]  Also:  {latest_csv}  <- NIS-E macro reads this fixed path")

    return records, out_csv


# ── watchdog daemon ───────────────────────────────────────────────────────────

class _ND2Handler(_FSH):
    def __init__(self, out_dir: Path, cfg: dict, processed: set[Path]):
        self.out_dir   = out_dir
        self.cfg       = cfg
        self.processed = processed

    def on_created(self, event):
        if event.is_directory:
            return
        p = Path(event.src_path)
        if p.suffix.lower() == ".nd2" and p not in self.processed:
            # Brief delay so NIS-E finishes writing the file before we open it
            time.sleep(2.0)
            try:
                run_once(p, self.out_dir, self.cfg)
                self.processed.add(p)
            except Exception as exc:
                print(f"[ERROR]  {p.name}: {exc}")


# ── CLI ───────────────────────────────────────────────────────────────────────

def _build_cfg(args) -> dict:
    return {
        "ch_idx":      args.ch,
        "min_diam_um": args.min_diam,
        "max_diam_um": args.max_diam,
        "sigma":       args.sigma,
        "thresh_k":    args.thresh_k,
        "w_circ":      args.w_circ,
        "w_sol":       args.w_sol,
        "backend":     args.backend,
    }


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(
        description="10X spheroid auto-screen: detect, rank, export for 20X capture")
    ap.add_argument("--nd2",       type=Path,  help="explicit ND2 file to process")
    ap.add_argument("--watch-dir", type=Path,  default=Path.cwd(),
                    help="folder to monitor for new 10X ND2 files")
    ap.add_argument("--out-dir",   type=Path,  default=None,
                    help="output directory for CSV (default: same as nd2 / watch-dir)")
    ap.add_argument("--daemon",    action="store_true",
                    help="watch continuously for new nd2 files")
    ap.add_argument("--ch",        type=int,   default=DEFAULT_CH_IDX,
                    help="0-based channel index for detection")
    ap.add_argument("--min-diam",  type=float, default=DEFAULT_MIN_DIAM_UM,
                    metavar="UM",  help=f"minimum spheroid diameter µm (default {DEFAULT_MIN_DIAM_UM})")
    ap.add_argument("--max-diam",  type=float, default=DEFAULT_MAX_DIAM_UM,
                    metavar="UM",  help=f"maximum spheroid diameter µm (default {DEFAULT_MAX_DIAM_UM})")
    ap.add_argument("--sigma",     type=float, default=DEFAULT_DETECT_SIGMA,
                    help=f"Gaussian pre-blur σ in pixels (default {DEFAULT_DETECT_SIGMA})")
    ap.add_argument("--thresh-k",  type=float, default=DEFAULT_DETECT_K,
                    help=f"threshold = mean + k·std (default {DEFAULT_DETECT_K})")
    ap.add_argument("--w-circ",    type=float, default=DEFAULT_W_CIRC,
                    help=f"circularity weight in composite score (default {DEFAULT_W_CIRC})")
    ap.add_argument("--w-sol",     type=float, default=DEFAULT_W_SOL,
                    help=f"solidity weight in composite score (default {DEFAULT_W_SOL})")
    ap.add_argument("--backend",   choices=["threshold", "cellpose"],
                    default=DEFAULT_BACKEND,
                    help="segmentation backend: 'threshold' (Gaussian+watershed, default) "
                         "or 'cellpose' (Cellpose cyto3, requires: pip install cellpose)")
    args = ap.parse_args(argv)

    if _MISSING:
        print("Missing required packages:")
        for m in _MISSING:
            print(f"  {m}")
        print("Install with: pip install nd2 numpy scipy scikit-image")
        sys.exit(1)

    cfg = _build_cfg(args)

    if args.nd2:
        out_dir = args.out_dir or args.nd2.parent
        run_once(args.nd2, out_dir, cfg)

    elif args.daemon:
        if not _WATCHDOG_OK:
            print("watchdog is not installed — daemon mode unavailable.")
            print("  pip install watchdog")
            sys.exit(1)
        watch_dir = args.watch_dir
        out_dir   = args.out_dir or watch_dir
        processed: set[Path] = set()

        # Screen any existing latest file before blocking
        existing = latest_nd2(watch_dir)
        if existing:
            run_once(existing, out_dir, cfg)
            processed.add(existing)

        handler  = _ND2Handler(out_dir, cfg, processed)
        observer = Observer()
        observer.schedule(handler, str(watch_dir), recursive=False)
        observer.start()
        print(f"Watching {watch_dir} for new .nd2 files … (Ctrl+C to stop)")
        try:
            while True:
                time.sleep(1)
        except KeyboardInterrupt:
            observer.stop()
        observer.join()

    else:
        # Default: process the latest nd2 in watch_dir and exit
        watch_dir = args.watch_dir
        out_dir   = args.out_dir or watch_dir
        nd2_path  = latest_nd2(watch_dir)
        if nd2_path is None:
            print(f"No .nd2 files found in {watch_dir}")
            sys.exit(1)
        run_once(nd2_path, out_dir, cfg)


if __name__ == "__main__":
    main()
