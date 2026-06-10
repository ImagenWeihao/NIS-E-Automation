"""
spheroid_pipeline.py — Cross-Zoom + CSV-to-Bin Pipeline Orchestrator

Drives the 4-step spheroid 20X imaging pipeline:
  1. Screen 10X mosaic  → ranked CSV (spheroid_screener.process_nd2)
  2. 2-anchor NCC       → stage-frame offset (cross_zoom_v2.verify)
  3. 20X autofocus      → z_centre per rank (read from captured nd2 metadata)
  4. Per-spheroid bin   → Beer-Lambert LP ramp → NIS-E file-flag trigger

Trigger/done protocol uses simple INI-style text files (NIS-E macro
language has no JSON parser):
  trigger.ini   written by Python, consumed by NIS-E macro
  done.ini      written by NIS-E macro, consumed by Python
"""
from __future__ import annotations

import csv
import json
import math
import os
import random
import time
from datetime import datetime

# ── Optional matplotlib stack (dashboard only) ────────────────────────────────
try:
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import matplotlib.patheffects as mpe
    import matplotlib.gridspec as mgridspec
    import numpy as np
    _MPL_OK = True
except ImportError:
    _MPL_OK = False
from dataclasses import asdict, dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Optional

import spheroid_screener as _screener
import cross_zoom_v2 as _czv2

# ── Constants ─────────────────────────────────────────────────────────────────

TRIGGER_FILENAME = "spheroid_trigger.ini"
DONE_FILENAME    = "spheroid_done.ini"
STATE_FILENAME   = "pipeline_state.json"
NCC_MIN_ACCEPT   = 0.40
MAX_ANCHOR_TRIES = 4      # max anchor attempts before aborting offset step
Z_HALF_DEFAULT   = 30.0   # ±µm around z_centre for z-stack range
Z_STEP_DEFAULT   = 2.0    # µm per slice


# ── Status enum ───────────────────────────────────────────────────────────────

class Status:
    DETECTED  = "DETECTED"
    VERIFIED  = "VERIFIED"   # stage offset applied
    Z_KNOWN   = "Z_KNOWN"    # z_centre recorded
    BIN_READY = "BIN_READY"  # .bin written
    QUEUED    = "QUEUED"     # trigger.ini written
    IMAGING   = "IMAGING"    # waiting for done.ini
    IMAGED    = "IMAGED"     # nd2 received and validated
    FAILED    = "FAILED"


# ── Data classes ──────────────────────────────────────────────────────────────

@dataclass
class SpheroidRecord:
    spheroid_id:      str
    rank:             int
    well:             str
    mosaic_cx_px:     float
    mosaic_cy_px:     float
    mosaic_diam_um:   float
    mosaic_score:     float
    mosaic_x_um:      float
    mosaic_y_um:      float
    mosaic_z_um:      float
    verified_x_um:    float = 0.0
    verified_y_um:    float = 0.0
    z_centre_um:      float = 0.0
    z_half_um:        float = Z_HALF_DEFAULT
    z_step_um:        float = Z_STEP_DEFAULT
    bin_path:         str   = ""
    nd2_out_path:     str   = ""
    status:           str   = Status.DETECTED


@dataclass
class PipelineState:
    well_id:        str
    mosaic_nd2:     str
    out_dir:        str
    trigger_dir:    str
    nd2_out_dir:    str
    records:        list[SpheroidRecord] = field(default_factory=list)
    offset_x:       float = 0.0
    offset_y:       float = 0.0
    offset_anchors: list[dict] = field(default_factory=list)


# ── Step 1: Screen mosaic ─────────────────────────────────────────────────────

def screen_mosaic(nd2_path: Path, out_dir: Path, well_id: str,
                  cfg: dict | None = None) -> list[SpheroidRecord]:
    """
    Run spheroid_screener on the 10X mosaic nd2.
    Returns SpheroidRecord list with Status.DETECTED.
    """
    cfg = cfg or {}
    out_dir.mkdir(parents=True, exist_ok=True)
    raw = _screener.process_nd2(nd2_path, cfg)
    _screener.export_csv(raw, out_dir / f"{nd2_path.stem}_screen.csv")
    _screener.export_csv(raw, out_dir / "spheroid_screen_latest.csv")

    records: list[SpheroidRecord] = []
    for r in raw:
        sph_id = f"sph_{well_id}_{r['rank']:02d}"
        records.append(SpheroidRecord(
            spheroid_id   = sph_id,
            rank          = r["rank"],
            well          = well_id,
            mosaic_cx_px  = float(r["cx_px"]),
            mosaic_cy_px  = float(r["cy_px"]),
            mosaic_diam_um= float(r["diameter_um"]),
            mosaic_score  = float(r["score"]),
            mosaic_x_um   = float(r["x_um"]),
            mosaic_y_um   = float(r["y_um"]),
            mosaic_z_um   = float(r["z_surface_um"]),
        ))
    return records


# ── Step 2: Anchor offset estimation ─────────────────────────────────────────

def pick_anchors(records: list[SpheroidRecord], n: int = 2,
                 seed: int | None = None) -> list[SpheroidRecord]:
    """Return n randomly chosen records as anchor candidates."""
    rng = random.Random(seed)
    return rng.sample(records, min(n, len(records)))


def _write_czv2_csv(records: list[SpheroidRecord], out_path: Path) -> None:
    """
    Write a screener CSV compatible with cross_zoom_v2.verify().
    That function expects canvas_col/canvas_row; screener exports cx_px/cy_px.
    For single-position mosaics they are identical; for multi-tile they differ.
    """
    fields = ["rank", "canvas_col", "canvas_row",
              "x_um", "y_um", "diameter_um", "score"]
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in records:
            w.writerow({
                "rank":       r.rank,
                "canvas_col": r.mosaic_cx_px,
                "canvas_row": r.mosaic_cy_px,
                "x_um":       r.mosaic_x_um,
                "y_um":       r.mosaic_y_um,
                "diameter_um": r.mosaic_diam_um,
                "score":      r.mosaic_score,
            })


def estimate_offset_from_nd2(
    mosaic_nd2: Path,
    sub10x_nd2: Path,
    records: list[SpheroidRecord],
    out_dir: Path,
    expected_rank: int | None = None,
) -> tuple[float, float, float, int]:
    """
    Run cross_zoom_v2.verify() on one sub-10X nd2.
    Returns (dx_um, dy_um, ncc_score, matched_rank).

    If expected_rank is given, prefer the NCC match for that rank (the spheroid
    the operator says they navigated to); otherwise use the highest-NCC match
    across all ranks. Raises ValueError if no match above threshold.
    """
    compat_csv = out_dir / "_czv2_compat.csv"
    _write_czv2_csv(records, compat_csv)

    matches = _czv2.verify(
        nd2_mosaic   = mosaic_nd2,
        nd2_sub10x   = sub10x_nd2,
        screener_csv = compat_csv,
        out_dir      = out_dir,
    )
    # Use only the best single match per anchor nd2.
    # Operator aims at one specific spheroid; secondary detections (e.g. near-edge
    # objects at pixel col < 120/1024) are noise and produce wrong offsets when averaged.
    good = sorted(matches, key=lambda m: m["ncc"], reverse=True)
    good = [m for m in good if m["ncc"] >= NCC_MIN_ACCEPT]
    if not good:
        raise ValueError(
            f"No NCC match >= {NCC_MIN_ACCEPT} in {sub10x_nd2.name}"
        )
    best = None
    if expected_rank is not None:
        best = next((m for m in good if int(m["matched_rank"]) == expected_rank), None)
    if best is None:
        best = good[0]
    return best["offset_x_um"], best["offset_y_um"], best["ncc"], int(best["matched_rank"])


def apply_offset(records: list[SpheroidRecord], dx: float, dy: float) -> None:
    """Apply stage offset to all records and advance status to VERIFIED."""
    for r in records:
        r.verified_x_um = round(r.mosaic_x_um + dx, 2)
        r.verified_y_um = round(r.mosaic_y_um + dy, 2)
        if r.status == Status.DETECTED:
            r.status = Status.VERIFIED


# ── Step 3: Automated autofocus trigger / global registration ─────────────────

AF_TRIGGER_PREFIX = "af_trigger_"
AF_DONE_PREFIX    = "af_done_"

# NIS-E macros read these files with Int_GetKeyString/Int_GetKeyValue, which
# require a Windows-style [section] header. parse_ini() ignores the header line
# (no '='), so the same files round-trip on the Python side.
INI_SECTION = "spheroid"


def _atomic_write_crlf(path: Path, lines: list[str]) -> None:
    """Write an INI file with explicit CRLF endings, atomically.

    CRLF: GetPrivateProfileString (which NIS-E Int_GetKeyString/Value wrap)
    only parses Windows-style line endings reliably -- don't rely on the OS
    text-mode translation in Path.write_text, which is a no-op off Windows.
    Atomic: the NIS-E daemon polls for the file by name, so an os.replace
    rename guarantees it only ever sees a fully-written file, never a partial.
    """
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text("\r\n".join(lines) + "\r\n", encoding="utf-8", newline="")
    os.replace(tmp, path)


def trigger_autofocus_all(records: list[SpheroidRecord], work_dir: Path) -> int:
    """Write per-rank autofocus trigger files for NIS-E macro. Returns count written."""
    work_dir.mkdir(parents=True, exist_ok=True)
    n = 0
    for r in records:
        p = work_dir / f"{AF_TRIGGER_PREFIX}{r.rank:02d}.ini"
        lines = [
            f"[{INI_SECTION}]",
            f"rank={r.rank}",
            f"spheroid_id={r.spheroid_id}",
            f"stage_x={r.verified_x_um}",
            f"stage_y={r.verified_y_um}",
        ]
        _atomic_write_crlf(p, lines)
        n += 1
    return n


def poll_autofocus_done(records: list[SpheroidRecord], work_dir: Path) -> dict[int, dict]:
    """Check for af_done_XX.ini files. Returns {rank: {z_centre_um, nd2_path}}.

    Only ranks the macro reports with status=ok are returned; af_failed /
    move_failed ranks are skipped so a bad autofocus never sets a z_centre.
    """
    done: dict[int, dict] = {}
    for r in records:
        p = work_dir / f"{AF_DONE_PREFIX}{r.rank:02d}.ini"
        if p.exists():
            d = parse_ini(p.read_text(encoding="utf-8"))
            if d.get("status", "ok") != "ok":
                continue
            try:
                done[r.rank] = {
                    "z_centre_um": float(d.get("z_centre", "0")),
                    "nd2_path":    d.get("nd2_path", ""),
                }
            except ValueError:
                pass
    return done


def apply_autofocus_results(
    records: list[SpheroidRecord],
    done_results: dict[int, dict],
    z_half_um: float,
    z_step_um: float,
) -> None:
    """Store z_centre from autofocus done results into records."""
    for r in records:
        if r.rank in done_results:
            r.z_centre_um = round(done_results[r.rank]["z_centre_um"], 2)
            r.z_half_um   = z_half_um
            r.z_step_um   = z_step_um
            if r.status in (Status.VERIFIED, Status.DETECTED):
                r.status = Status.Z_KNOWN


def global_registration(
    records: list[SpheroidRecord],
    mosaic_nd2: Path,
    out_dir: Path,
    done_nd2s: dict[int, str],
) -> dict:
    """
    NCC-match each captured 20X ND2 vs mosaic to get precise XY per spheroid.
    Fit similarity transform across all matches; update verified_x/y on all records.
    Returns summary dict with n_matched, residuals_um, mean_residual_um, transform.
    """
    import math as _math2

    src: list[tuple[float, float]] = []
    dst: list[tuple[float, float]] = []

    for rank, nd2_path in done_nd2s.items():
        if not nd2_path:
            continue
        rec = next((r for r in records if r.rank == rank), None)
        if rec is None:
            continue
        try:
            dx, dy, _ncc, _mr = estimate_offset_from_nd2(
                mosaic_nd2, Path(nd2_path), records, out_dir, expected_rank=rank)
            src.append((rec.verified_x_um, rec.verified_y_um))
            dst.append((rec.mosaic_x_um + dx, rec.mosaic_y_um + dy))
        except (ValueError, Exception):
            pass

    if not src:
        raise ValueError("No successful NCC matches for global registration")

    n = len(src)

    try:
        import numpy as np
        src_a = np.array(src)
        dst_a = np.array(dst)

        if n >= 2:
            from skimage.transform import SimilarityTransform
            tform = SimilarityTransform()
            tform.estimate(src_a, dst_a)
            pred = tform(src_a)
            residuals = [float(np.linalg.norm(pred[i] - dst_a[i])) for i in range(n)]
            for r in records:
                xy_t = tform(np.array([[r.verified_x_um, r.verified_y_um]]))[0]
                r.verified_x_um = round(float(xy_t[0]), 2)
                r.verified_y_um = round(float(xy_t[1]), 2)
            return {
                "n_matched": n,
                "residuals_um": residuals,
                "mean_residual_um": float(np.mean(residuals)),
                "transform": {
                    "dx": float(tform.translation[0]),
                    "dy": float(tform.translation[1]),
                    "scale": float(tform.scale),
                    "angle_deg": float(_math2.degrees(tform.rotation)),
                },
            }
        else:
            ddx = float(dst_a[0, 0] - src_a[0, 0])
            ddy = float(dst_a[0, 1] - src_a[0, 1])
            for r in records:
                r.verified_x_um = round(r.verified_x_um + ddx, 2)
                r.verified_y_um = round(r.verified_y_um + ddy, 2)
            return {
                "n_matched": 1, "residuals_um": [0.0],
                "mean_residual_um": 0.0,
                "transform": {"dx": ddx, "dy": ddy},
            }
    except ImportError:
        ddx = sum(d[0] - s[0] for s, d in zip(src, dst)) / n
        ddy = sum(d[1] - s[1] for s, d in zip(src, dst)) / n
        for r in records:
            r.verified_x_um = round(r.verified_x_um + ddx, 2)
            r.verified_y_um = round(r.verified_y_um + ddy, 2)
        residuals = [
            _math2.sqrt((dst[i][0] - src[i][0] - ddx) ** 2 + (dst[i][1] - src[i][1] - ddy) ** 2)
            for i in range(n)
        ]
        return {
            "n_matched": n, "residuals_um": residuals,
            "mean_residual_um": sum(residuals) / n,
            "transform": {"dx": ddx, "dy": ddy},
        }


# ── Step 3: Z-centre from 20X autofocus nd2 (manual fallback) ──────────────────

def _read_nd2_z(nd2_path: Path) -> float:
    """Extract stage Z (µm) from frame 0 of an nd2 file."""
    try:
        import nd2 as _nd2lib
        with _nd2lib.ND2File(nd2_path) as f:
            fm = f.frame_metadata(0)
            return float(fm.channels[0].position.stagePositionUm.z)
    except Exception as exc:
        raise ValueError(f"Cannot read Z from {nd2_path.name}: {exc}") from exc


def record_z_centre(records: list[SpheroidRecord], rank: int,
                    nd2_path: Path,
                    z_half_um: float = Z_HALF_DEFAULT,
                    z_step_um: float = Z_STEP_DEFAULT) -> None:
    """
    Read Z from a 20X single-plane nd2 (acquired after autofocus) and
    store it as z_centre for the given rank.
    """
    z = _read_nd2_z(nd2_path)
    for r in records:
        if r.rank == rank:
            r.z_centre_um = round(z, 2)
            r.z_half_um   = z_half_um
            r.z_step_um   = z_step_um
            if r.status in (Status.VERIFIED, Status.DETECTED):
                r.status = Status.Z_KNOWN
            return
    raise ValueError(f"No record with rank={rank}")


# ── Step 4a: Generate per-spheroid .bin ──────────────────────────────────────

def _beer_lambert(P0: float, L: float, depth_um: float) -> float:
    return min(P0 * math.exp(depth_um / L), 100.0)


def _build_bin_items(z_start: float, z_end: float, z_step: float,
                     P0: float, L: float, ch_field: str) -> list[dict]:
    """
    Build Z Intensity Correction items with a Beer-Lambert LP ramp.
    z_start is the top of the spheroid (lowest depth); depth = z - z_start.
    """
    z = z_start
    items: list[dict] = []
    while z <= z_end + 1e-6:
        depth = max(0.0, z - z_start)
        lp    = _beer_lambert(P0, L, depth)
        items.append({
            "ZStack":   round(z, 3),
            ch_field:   round(lp, 3),
        })
        z = round(z + z_step, 4)
    return items


def _write_bin(metadata: dict, items: list[dict], out_path: Path) -> None:
    """
    Write a .bin via csv_to_nis_bin.build_bin().
    NOTE: ROOT_TAIL_U64S in csv_to_nis_bin.py is hardcoded for 19 items;
    bins with n_items != 19 may be malformed — see the TODO in CLAUDE.md.
    For the test workflow, configure z_half/z_step so that n_items == 19
    (e.g., z_half=18 um, z_step=2 um → 19 items).
    """
    from csv_to_nis_bin import build_bin
    blob = build_bin(metadata, items)
    out_path.write_bytes(blob)


def generate_bin(record: SpheroidRecord,
                 P0_pct: float, L_um: float,
                 ch_field: str,
                 bin_dir: Path,
                 bin_metadata: dict | None = None) -> Path:
    """
    Generate a per-spheroid Z Intensity Correction .bin file.
    ch_field: internal field name, e.g. "CH1LaserPower"
    Returns the path to the written bin file.
    """
    if record.z_centre_um == 0.0:
        raise ValueError(f"Rank {record.rank}: z_centre_um not set — run Step 3 first")

    z_start = round(record.z_centre_um - record.z_half_um, 3)
    z_end   = round(record.z_centre_um + record.z_half_um, 3)
    items   = _build_bin_items(z_start, z_end, record.z_step_um, P0_pct, L_um, ch_field)

    bin_dir.mkdir(parents=True, exist_ok=True)
    out_path = bin_dir / f"{record.spheroid_id}.bin"
    _write_bin(bin_metadata or {}, items, out_path)

    record.bin_path = str(out_path)
    record.status   = Status.BIN_READY
    return out_path


# ── Step 4b: File-flag trigger protocol ──────────────────────────────────────

def write_trigger(record: SpheroidRecord, trigger_dir: Path,
                  nd2_out_dir: Path) -> Path:
    """
    Write trigger.ini for the NIS-E macro.
    INI-style key=value, one per line (no JSON — NIS-E macro can't parse it).
    """
    trigger_dir.mkdir(parents=True, exist_ok=True)

    nd2_out = Path(nd2_out_dir) / f"{record.spheroid_id}.nd2"
    trigger_path = trigger_dir / TRIGGER_FILENAME

    lines = [
        f"[{INI_SECTION}]",
        f"rank={record.rank}",
        f"spheroid_id={record.spheroid_id}",
        f"stage_x={record.verified_x_um}",
        f"stage_y={record.verified_y_um}",
        f"z_centre={record.z_centre_um}",
        f"z_half={record.z_half_um}",
        f"z_step={record.z_step_um}",
        f"bin_path={record.bin_path}",
        f"nd2_out={nd2_out}",
        f"timestamp={datetime.now().isoformat(timespec='seconds')}",
    ]
    _atomic_write_crlf(trigger_path, lines)

    record.nd2_out_path = str(nd2_out)
    record.status       = Status.QUEUED
    return trigger_path


def parse_ini(text: str) -> dict:
    out: dict[str, str] = {}
    for line in text.splitlines():
        line = line.strip()
        if "=" in line and not line.startswith("#"):
            k, _, v = line.partition("=")
            out[k.strip()] = v.strip()
    return out


def wait_for_done(trigger_dir: Path, timeout_s: float = 300.0,
                  poll_s: float = 2.0) -> dict | None:
    """
    Poll for done.ini written by NIS-E macro.
    Returns parsed dict on success, None on timeout.
    """
    done_path  = trigger_dir / DONE_FILENAME
    deadline   = time.monotonic() + timeout_s
    while time.monotonic() < deadline:
        if done_path.exists():
            data = parse_ini(done_path.read_text(encoding="utf-8"))
            done_path.unlink(missing_ok=True)
            return data
        time.sleep(poll_s)
    return None


# ── State persistence ─────────────────────────────────────────────────────────

def save_state(state: PipelineState) -> None:
    out_dir = Path(state.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    d = asdict(state)
    (out_dir / STATE_FILENAME).write_text(
        json.dumps(d, indent=2), encoding="utf-8"
    )


def load_state(out_dir: Path) -> PipelineState | None:
    p = out_dir / STATE_FILENAME
    if not p.exists():
        return None
    d = json.loads(p.read_text(encoding="utf-8"))
    records = [SpheroidRecord(**r) for r in d.pop("records", [])]
    return PipelineState(**d, records=records)


# ── CSV export of pipeline state ──────────────────────────────────────────────

_STATE_COLS = [
    "spheroid_id", "rank", "well", "status",
    "mosaic_x_um", "mosaic_y_um", "mosaic_z_um",
    "verified_x_um", "verified_y_um",
    "z_centre_um", "z_half_um", "z_step_um",
    "mosaic_diam_um", "mosaic_score",
    "bin_path", "nd2_out_path",
]


def export_state_csv(records: list[SpheroidRecord], out_path: Path) -> None:
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=_STATE_COLS, extrasaction="ignore")
        w.writeheader()
        w.writerows(asdict(r) for r in records)


# ── Pipeline Dashboard ─────────────────────────────────────────────────────────

# Catppuccin Mocha palette
_C = {
    "bg":      "#1e1e2e",
    "bg2":     "#181825",
    "surf":    "#313244",
    "surf2":   "#45475a",
    "sub":     "#6c7086",
    "text":    "#cdd6f4",
    "text2":   "#bac2de",
    "green":   "#a6e3a1",
    "red":     "#f38ba8",
    "blue":    "#89b4fa",
    "sky":     "#89dceb",
    "yellow":  "#f9e2af",
    "peach":   "#fab387",
    "mauve":   "#cba6f7",
    "lavender":"#b4befe",
}

_STATUS_STYLE: dict[str, dict] = {
    # status  → facecolor, alpha, label
    Status.DETECTED:  {"fc": _C["sub"],     "alpha": 0.25, "lbl": "Detected"},
    Status.VERIFIED:  {"fc": _C["sky"],     "alpha": 0.50, "lbl": "Verified"},
    Status.Z_KNOWN:   {"fc": _C["blue"],    "alpha": 0.65, "lbl": "Z known"},
    Status.BIN_READY: {"fc": _C["yellow"],  "alpha": 0.75, "lbl": "Bin ready"},
    Status.QUEUED:    {"fc": _C["peach"],   "alpha": 0.85, "lbl": "Queued"},
    Status.IMAGING:   {"fc": _C["mauve"],   "alpha": 0.92, "lbl": "Imaging"},
    Status.IMAGED:    {"fc": _C["green"],   "alpha": 1.00, "lbl": "Imaged"},
    Status.FAILED:    {"fc": _C["red"],     "alpha": 1.00, "lbl": "Failed"},
}

_ROW_BG: dict[str, str] = {
    Status.DETECTED:  "#222236",
    Status.VERIFIED:  "#1d2b33",
    Status.Z_KNOWN:   "#1c2a3a",
    Status.BIN_READY: "#2b2b1e",
    Status.QUEUED:    "#2b2419",
    Status.IMAGING:   "#26213a",
    Status.IMAGED:    "#1b2b1e",
    Status.FAILED:    "#2b1b1e",
}


def _load_mosaic_img(nd2_path: Path, max_px: int = 1500
                     ) -> tuple[np.ndarray, float, float, float, int, int]:
    """
    Load mosaic nd2, MIP position-0 channel-0, downsample.
    Returns (display_img, display_pum, orig_pum, cx_um, cy_um, W_orig, H_orig).
    """
    import nd2 as _nd2lib
    from scipy.ndimage import zoom as _zoom

    with _nd2lib.ND2File(nd2_path) as f:
        pum_orig = float(f.voxel_size().x)
        try:
            fm  = f.frame_metadata(0)
            pos = fm.channels[0].position.stagePositionUm
            cx, cy = float(pos.x), float(pos.y)
        except Exception:
            cx, cy = 0.0, 0.0

        arr  = f.asarray()
        dims = list(f.sizes.keys())

        for dim_key, idx in [("P", 0), ("C", 0)]:
            if dim_key in dims:
                ax = dims.index(dim_key)
                sl = [slice(None)] * arr.ndim; sl[ax] = idx
                arr  = arr[tuple(sl)]
                dims = [d for d in dims if d != dim_key]
        if "Z" in dims:
            arr = arr.max(axis=dims.index("Z"))

    H_orig, W_orig = arr.shape
    arr = arr.astype(np.float32)

    ds = max(H_orig, W_orig) / max_px
    if ds > 1.0:
        arr = _zoom(arr, 1.0 / ds, order=1)
        pum_disp = pum_orig * ds
    else:
        pum_disp = pum_orig

    lo = float(np.percentile(arr, 1.0))
    hi = float(np.percentile(arr, 99.5))
    arr = np.clip((arr - lo) / (hi - lo + 1e-8), 0, 1).astype(np.float32)
    return arr, pum_disp, pum_orig, cx, cy, W_orig, H_orig


class PipelineDashboard:
    """
    Live SpheroidPA pipeline dashboard.

    Top panel  — 10X mosaic with spheroid circles and rank-order zig-zag path.
                 Circles start 25% opaque (DETECTED) and fill to 100% (IMAGED).
    Bottom panel — status table: one row per spheroid, color-coded by status.

    Usage (non-blocking updates during pipeline):
        dash = PipelineDashboard(mosaic_nd2=Path("..."), out_dir=Path("..."))
        dash.update(records)   # call after each status change → saves PNG
        dash.close()           # clean up when done

    Standalone replay:
        python spheroid_pipeline.py demo_data/pipeline_dryrun/pipeline_state.json
    """

    def __init__(self, mosaic_nd2: Path | None = None,
                 out_dir: Path | None = None,
                 interactive: bool = False,
                 figsize: tuple = (15, 10),
                 fig=None):
        if not _MPL_OK:
            raise ImportError("matplotlib and numpy required for PipelineDashboard")

        self._out_dir     = out_dir
        self._interactive = interactive
        self._mos_img     = None
        self._mos_extent  = None    # [xmin, xmax, ymax, ymin] stage µm for imshow

        # Load mosaic once
        if mosaic_nd2 and mosaic_nd2.exists():
            try:
                img, pum_d, pum_o, cx, cy, W, H = _load_mosaic_img(mosaic_nd2)
                self._mos_img = img
                hw  = W / 2 * pum_o
                hh  = H / 2 * pum_o
                # imshow extent: [left, right, bottom, top]
                # stage Y increases upward; image row-0 = stage y_max
                self._mos_extent = [cx - hw, cx + hw, cy + hh, cy - hh]
                print(f"[Dashboard] Mosaic loaded: {img.shape[1]}x{img.shape[0]} px "
                      f"(downsampled), stage extent "
                      f"X=[{cx-hw:.0f}, {cx+hw:.0f}] "
                      f"Y=[{cy-hh:.0f}, {cy+hh:.0f}] um")
            except Exception as exc:
                print(f"[Dashboard] Mosaic load failed: {exc}")

        # Figure
        if fig is not None:
            self._fig = fig
            self._fig.set_facecolor(_C["bg"])
        else:
            if interactive:
                matplotlib.use("TkAgg")
                plt.ion()
            else:
                matplotlib.use("Agg")
            self._fig = plt.figure(figsize=figsize, facecolor=_C["bg"])

        gs = mgridspec.GridSpec(
            2, 1, figure=self._fig,
            height_ratios=[3, 1], hspace=0.06,
            left=0.06, right=0.97, top=0.96, bottom=0.03,
        )
        self._ax_mos = self._fig.add_subplot(gs[0])
        self._ax_tbl = self._fig.add_subplot(gs[1])
        for ax in (self._ax_mos, self._ax_tbl):
            ax.set_facecolor(_C["bg2"])
            for sp in ax.spines.values():
                sp.set_edgecolor(_C["surf2"])

    # ── Public API ────────────────────────────────────────────────────────────

    def update(self, records: list[SpheroidRecord],
               save_path: Path | None = None) -> None:
        """Redraw the dashboard with current record states."""
        self._draw_mosaic(records)
        self._draw_table(records)
        self._fig.canvas.draw_idle()
        if self._interactive:
            plt.pause(0.05)

        # Auto-save
        dst = save_path or (self._out_dir / "pipeline_dashboard.png"
                            if self._out_dir else None)
        if dst:
            dst.parent.mkdir(parents=True, exist_ok=True)
            self._fig.savefig(dst, dpi=150, bbox_inches="tight",
                              facecolor=_C["bg"])

    def save(self, path: Path) -> None:
        self._fig.savefig(path, dpi=150, bbox_inches="tight", facecolor=_C["bg"])

    def show(self, block: bool = True) -> None:
        plt.show(block=block)

    def close(self) -> None:
        plt.close(self._fig)

    # ── Mosaic panel ──────────────────────────────────────────────────────────

    def _draw_mosaic(self, records: list[SpheroidRecord]) -> None:
        ax = self._ax_mos
        ax.cla()
        ax.set_facecolor(_C["bg2"])
        for sp in ax.spines.values():
            sp.set_edgecolor(_C["surf2"])

        # Background image
        if self._mos_img is not None:
            ax.imshow(self._mos_img, cmap="bone",
                      extent=self._mos_extent,
                      aspect="equal", interpolation="bilinear",
                      alpha=0.70, zorder=0)

        if not records:
            ax.text(0.5, 0.5, "No records yet",
                    ha="center", va="center", color=_C["sub"],
                    transform=ax.transAxes, fontsize=11)
            return

        recs = sorted(records, key=lambda r: r.rank)

        # ── Zig-zag capture path ──────────────────────────────────────────────
        # Use verified XY if available, else mosaic XY
        def _xy(r: SpheroidRecord) -> tuple[float, float]:
            # Always use mosaic XY so circles stay aligned with the mosaic image.
            # Verified XY is a different stage frame (sub-10X encoder zero) and
            # is shown in the table only.
            return r.mosaic_x_um, r.mosaic_y_um

        xs = [_xy(r)[0] for r in recs]
        ys = [_xy(r)[1] for r in recs]

        # Full path: faint dashed
        ax.plot(xs, ys, color=_C["surf2"], linewidth=0.9,
                linestyle=(0, (5, 4)), alpha=0.55, zorder=1)

        # Completed legs: solid green
        for i in range(len(recs) - 1):
            a, b = recs[i], recs[i + 1]
            if a.status == Status.IMAGED:
                ax.plot([_xy(a)[0], _xy(b)[0]],
                        [_xy(a)[1], _xy(b)[1]],
                        color=_C["green"], linewidth=1.4,
                        linestyle="-", alpha=0.75, zorder=2)

        # Currently imaging: animated-looking segment
        for r in recs:
            if r.status == Status.IMAGING:
                ax.plot(*_xy(r), "o", color=_C["mauve"],
                        markersize=20, alpha=0.25, zorder=3)

        # ── Spheroid circles ──────────────────────────────────────────────────
        for r in recs:
            x, y = _xy(r)
            radius_um  = r.mosaic_diam_um / 2.0
            st = _STATUS_STYLE.get(r.status, _STATUS_STYLE[Status.DETECTED])
            fc    = st["fc"]
            alpha = st["alpha"]

            # Filled disc (opacity encodes pipeline progress)
            disc = mpatches.Circle((x, y), radius_um,
                                   facecolor=fc, edgecolor="none",
                                   alpha=alpha, zorder=4)
            ax.add_patch(disc)

            # Edge ring — always visible
            ring = mpatches.Circle((x, y), radius_um,
                                   facecolor="none",
                                   edgecolor=fc,
                                   linewidth=1.4,
                                   alpha=min(alpha + 0.25, 1.0),
                                   zorder=5)
            ax.add_patch(ring)

            # Rank number below the circle
            ax.text(x, y - radius_um - 35,
                    f"#{r.rank}",
                    ha="center", va="top",
                    fontsize=7.5, color=fc,
                    fontfamily="monospace",
                    fontweight="bold",
                    alpha=min(alpha + 0.35, 1.0),
                    zorder=6,
                    path_effects=[mpe.withStroke(linewidth=1.5,
                                                 foreground=_C["bg"])])

            # Check / cross badge inside circle for terminal states
            if r.status == Status.IMAGED:
                ax.text(x, y, "✓", ha="center", va="center",
                        fontsize=10, color=_C["bg"],
                        fontweight="bold", zorder=7)
            elif r.status == Status.FAILED:
                ax.text(x, y, "✗", ha="center", va="center",
                        fontsize=10, color=_C["bg"],
                        fontweight="bold", zorder=7)
            elif r.status == Status.IMAGING:
                ax.text(x, y, "▶", ha="center", va="center",
                        fontsize=8, color=_C["bg"],
                        fontweight="bold", zorder=7)

        # ── Axis styling ──────────────────────────────────────────────────────
        ax.set_xlabel("Stage X (µm)", color=_C["text2"], fontsize=8,
                      labelpad=2)
        ax.set_ylabel("Stage Y (µm)", color=_C["text2"], fontsize=8,
                      labelpad=2)
        ax.tick_params(colors=_C["sub"], labelsize=7)
        ax.set_aspect("equal", adjustable="datalim")

        # Title with timestamp and counts
        n_done = sum(1 for r in records if r.status == Status.IMAGED)
        n_fail = sum(1 for r in records if r.status == Status.FAILED)
        ts     = datetime.now().strftime("%H:%M:%S")
        ax.set_title(
            f"SpheroidPA Pipeline  ·  {n_done}/{len(records)} imaged"
            + (f"  ·  {n_fail} failed" if n_fail else "")
            + f"  ·  {ts}",
            color=_C["text"], fontsize=9, pad=5, fontfamily="monospace")

        # Legend — status color key
        handles = [
            mpatches.Patch(facecolor=v["fc"], alpha=v["alpha"],
                           edgecolor=v["fc"], linewidth=0.8, label=v["lbl"])
            for v in _STATUS_STYLE.values()
        ]
        ax.legend(handles=handles, loc="upper right",
                  fontsize=6, ncol=2, framealpha=0.85,
                  facecolor=_C["surf"], edgecolor=_C["surf2"],
                  labelcolor=_C["text2"])

    # ── Status table panel ────────────────────────────────────────────────────

    def _draw_table(self, records: list[SpheroidRecord]) -> None:
        ax = self._ax_tbl
        ax.cla()
        ax.axis("off")
        ax.set_facecolor(_C["bg2"])

        if not records:
            return

        recs = sorted(records, key=lambda r: r.rank)

        COLS = ["Rank", "Spheroid ID", "Status",
                "Mosaic XY (µm)", "Verified XY (µm)", "Z-centre (µm)",
                "Diam (µm)", "Score", "Bin file"]
        COL_W = [0.04, 0.14, 0.08, 0.13, 0.13, 0.10, 0.07, 0.06, 0.25]

        cell_text = []
        for r in recs:
            mos_xy = f"({r.mosaic_x_um:.0f}, {r.mosaic_y_um:.0f})"
            ver_xy = (f"({r.verified_x_um:.0f}, {r.verified_y_um:.0f})"
                      if r.verified_x_um else "—")
            zc     = f"{r.z_centre_um:.1f}" if r.z_centre_um else "—"
            bn     = Path(r.bin_path).name if r.bin_path else "—"
            cell_text.append([
                str(r.rank), r.spheroid_id, r.status,
                mos_xy, ver_xy, zc,
                f"{r.mosaic_diam_um:.1f}", f"{r.mosaic_score:.4f}", bn,
            ])

        tbl = ax.table(
            cellText=cell_text,
            colLabels=COLS,
            colWidths=COL_W,
            cellLoc="center",
            loc="center",
            bbox=[0, 0, 1, 1],
        )
        tbl.auto_set_font_size(False)
        tbl.set_fontsize(6.5)

        # Header row
        for ci in range(len(COLS)):
            cell = tbl[(0, ci)]
            cell.set_facecolor(_C["surf"])
            cell.set_edgecolor(_C["surf2"])
            cell.set_text_props(color=_C["lavender"], fontweight="bold",
                                fontfamily="monospace")

        # Data rows — color by status
        for ri, r in enumerate(recs):
            row_bg   = _ROW_BG.get(r.status, _C["bg2"])
            text_col = _STATUS_STYLE.get(
                r.status, _STATUS_STYLE[Status.DETECTED])["fc"]
            for ci in range(len(COLS)):
                cell = tbl[(ri + 1, ci)]
                cell.set_facecolor(row_bg)
                cell.set_edgecolor(_C["surf"])
                cell.set_text_props(color=text_col, fontfamily="monospace")


# ── Standalone CLI ─────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import sys
    sys.stdout.reconfigure(encoding="utf-8")

    if len(sys.argv) < 2:
        print("Usage: python spheroid_pipeline.py <pipeline_state.json>")
        sys.exit(1)

    state = load_state(Path(sys.argv[1]).parent)
    if state is None:
        print("No pipeline_state.json found.")
        sys.exit(1)

    mosaic = Path(state.mosaic_nd2) if state.mosaic_nd2 else None
    dash = PipelineDashboard(
        mosaic_nd2  = mosaic,
        out_dir     = Path(state.out_dir),
        interactive = True,
    )
    dash.update(state.records)
    dash.show(block=True)
