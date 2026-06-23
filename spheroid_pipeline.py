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


def apply_anchor_transform(records: list[SpheroidRecord], corrs: list) -> dict:
    """Map every record's mosaic XY -> verified (stage) XY from anchor correspondences.

    corrs: list of ((mosaic_x, mosaic_y), (stage_x, stage_y)) pairs, one per verified
    anchor (stage = the anchor's matched mosaic point + its measured dx/dy).

    With >= 2 anchors, fit a SimilarityTransform (scale + rotation + translation),
    which captures a mosaic<->stage axis flip (a ~180 deg rotation) that a single
    averaged translation cannot -- that flip is why captures landed off the spheroids.
    With 1 anchor, fall back to a pure translation. Returns a summary dict for the GUI.
    """
    if len(corrs) >= 2:
        import math
        import numpy as np
        from skimage.transform import SimilarityTransform
        src = np.array([c[0] for c in corrs], dtype=float)
        dst = np.array([c[1] for c in corrs], dtype=float)
        tform = SimilarityTransform()
        if not tform.estimate(src, dst):
            raise ValueError("anchor transform fit failed (degenerate anchors?)")
        pts = np.array([(r.mosaic_x_um, r.mosaic_y_um) for r in records], dtype=float)
        mapped = tform(pts)
        for r, xy in zip(records, mapped):
            r.verified_x_um = round(float(xy[0]), 2)
            r.verified_y_um = round(float(xy[1]), 2)
            if r.status == Status.DETECTED:
                r.status = Status.VERIFIED
        return {"mode": "similarity", "n": len(corrs), "scale": float(tform.scale),
                "rotation_deg": math.degrees(tform.rotation),
                "tx": float(tform.translation[0]), "ty": float(tform.translation[1])}

    # 1 anchor: a pure translation can't represent the characterized mosaic<->stage
    # ~180 deg axis flip, so only the matched spheroid would land. Apply the known
    # flip (R = -I) and pin the translation from this one anchor:
    #   stage = -mosaic + t,  with  t = stage_a + mosaic_a.
    # This is a rough placement -- a 2nd good anchor gives the precise fitted
    # rotation/scale (the >=2-anchor branch above).
    (mx, my), (sx, sy) = corrs[0]
    tx, ty = sx + mx, sy + my
    for r in records:
        r.verified_x_um = round(tx - r.mosaic_x_um, 2)
        r.verified_y_um = round(ty - r.mosaic_y_um, 2)
        if r.status == Status.DETECTED:
            r.status = Status.VERIFIED
    return {"mode": "flip1", "n": 1, "tx": tx, "ty": ty}


# ── Step 3: Automated autofocus trigger / global registration ─────────────────

AF_TRIGGER_PREFIX = "af_trigger_"
AF_DONE_PREFIX    = "af_done_"

# NIS-E macros read these files with Int_GetKeyString/Int_GetKeyValue, which
# require a Windows-style [section] header. parse_ini() ignores the header line
# (no '='), so the same files round-trip on the Python side.
INI_SECTION = "spheroid"

# Fixed, space-free pointer file the NIS-E capture daemon reads to find this run's
# folders (proven by macro_selftest/test_spacedir.mac). trigger_autofocus_all writes
# it so the daemon follows the GUI's Step 1 save dir instead of a hardcoded path.
SESSION_INI = Path("C:/SpheroidPA/session.ini")


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


def trigger_autofocus_all(records: list[SpheroidRecord], work_dir: Path,
                          z_centre: float | None = None,
                          z_half: float | None = None,
                          z_step: float | None = None) -> int:
    """Write per-rank trigger files for the NIS-E capture daemon. Returns count written.

    z_centre/z_half/z_step (the Z-stack geometry from GUI Step 3) are written into
    every trigger so nis_macro_capture_zstack.mac reads them instead of hardcoding.
    They are omitted when None (the macro then falls back to its built-in defaults).

    Two safeguards (added after a capture run imaged duplicated coordinates — the
    daemon had consumed a folder holding a mix of fresh and stale triggers from
    earlier re-trigger cycles):

      1. Verify-distinct (pre-check): every record must have a distinct
         (stage_x, stage_y). A duplicate means a corrupted records list, and
         capturing it would image the same spot twice — raise before writing
         anything or touching the instrument.
      2. Clear-stale: remove any existing af_trigger_* / af_done_* in work_dir
         first, so the daemon only ever sees this run's triggers (no leftovers
         from a prior cycle).
    """
    # 1. verify-distinct, before we write or delete anything
    seen: dict[tuple, int] = {}
    for r in records:
        key = (round(r.verified_x_um, 2), round(r.verified_y_um, 2))
        if key in seen:
            raise ValueError(
                f"Duplicate trigger coordinate {key} for rank {r.rank} and rank "
                f"{seen[key]} — records list is corrupted; aborting before capture."
            )
        seen[key] = r.rank

    work_dir.mkdir(parents=True, exist_ok=True)

    # 2. clear stale triggers/done from any previous re-trigger cycle
    for stale in (list(work_dir.glob(f"{AF_TRIGGER_PREFIX}*.ini"))
                  + list(work_dir.glob(f"{AF_DONE_PREFIX}*.ini"))):
        stale.unlink(missing_ok=True)

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
        if z_centre is not None:
            lines.append(f"z_centre={z_centre}")
        if z_half is not None:
            lines.append(f"z_half={z_half}")
        if z_step is not None:
            lines.append(f"z_step={z_step}")
        _atomic_write_crlf(p, lines)
        n += 1

    # Point the NIS-E capture daemon at THIS run's folders (no hardcoded path in the
    # macro): write the fixed pointer file it reads. work_dir is <save-dir>/autofocus,
    # so the nd2 captures go to the sibling <save-dir>/nd2.
    nd2_dir = work_dir.parent / "nd2"
    nd2_dir.mkdir(parents=True, exist_ok=True)
    SESSION_INI.parent.mkdir(parents=True, exist_ok=True)
    _atomic_write_crlf(SESSION_INI, [
        "[paths]",
        f"work_dir={work_dir.as_posix()}",
        f"nd2_dir={nd2_dir.as_posix()}",
    ])
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


def recenter_from_captures(records: list[SpheroidRecord], nd2_dir: Path,
                           flip_x: bool = False, flip_y: bool = False,
                           suffix: str = "_zstack.nd2") -> list[tuple]:
    """Second-pass re-centering. For each record, read its capture
    <spheroid_id><suffix> in nd2_dir, measure the spheroid's intensity-weighted
    centroid offset from the frame centre, and nudge verified_x/y so the NEXT
    capture lands centred -- exact per-spheroid, independent of anchor quality.

    flip_x/flip_y set the camera-pixel -> stage-um axis sign (verify on one
    spheroid; toggle a flip if it nudges the wrong way). Pixel size is read from
    each nd2 (so it works at any magnification). Returns
    [(rank, dcol_px, drow_px, dx_um, dy_um), ...] for review.
    """
    import numpy as np
    import nd2 as _nd2
    sx = -1.0 if flip_x else 1.0
    sy = -1.0 if flip_y else 1.0
    out: list[tuple] = []
    for r in records:
        p = Path(nd2_dir) / f"{r.spheroid_id}{suffix}"
        if not p.exists():
            continue
        with _nd2.ND2File(p) as f:
            a = np.squeeze(np.asarray(f.asarray())).astype(np.float32)
            try:
                px = float(f.voxel_size().x)
            except Exception:
                px = 0.0
        if px <= 0.0:
            continue
        while a.ndim > 3:
            a = a[a.shape[0] // 2]
        proj = a.max(0) if a.ndim == 3 else a          # max over Z
        h, w = proj.shape
        bg = float(np.median(proj))
        mad = float(np.median(np.abs(proj - bg))) * 1.4826 + 1e-6
        mask = proj > bg + max(60.0, 6.0 * mad)        # bright spheroid pixels
        if int(mask.sum()) < 100:
            continue
        # TODO(re-center): two spheroids close together in one FOV -> this whole-frame
        # intensity centroid lands BETWEEN them, so the nudge moves to the midpoint
        # instead of the target (e.g. adjacent ranks #3/#4). Fix: label connected components
        # in `mask` and use the component nearest the frame centre (the targeted
        # spheroid), not the global centroid.
        ys, xs = np.nonzero(mask)
        wgt = proj[ys, xs] - bg
        cx = float((xs * wgt).sum() / wgt.sum())
        cy = float((ys * wgt).sum() / wgt.sum())
        dcol = cx - w / 2.0
        drow = cy - h / 2.0
        # move the stage opposite the image-centroid offset to bring it to centre
        dx = round(-dcol * px * sx, 2)
        dy = round(-drow * px * sy, 2)
        r.verified_x_um = round(r.verified_x_um + dx, 2)
        r.verified_y_um = round(r.verified_y_um + dy, 2)
        out.append((r.rank, round(dcol, 1), round(drow, 1), dx, dy))
    return out


def focus_scores(stack, method: str = "variance") -> list:
    """Per-plane sharpness score for a (Z,Y,X) stack -- higher = sharper, background-
    subtracted. `method`: variance (normalized), tenengrad, laplacian, brenner, fft.
    (Trialled in v1.4.3; none matched the eye reliably, so they're a manual tool here.)"""
    import numpy as np
    a = np.asarray(stack, dtype=np.float32)
    if a.ndim == 2:
        a = a[None, ...]
    while a.ndim > 3:
        a = a[a.shape[0] // 2]
    m = (method or "variance").lower()
    scores = []
    for p in a:
        bg = float(np.median(p))
        q = np.clip(p - bg, 0.0, None)
        if "tenengrad" in m:
            gy, gx = np.gradient(q)
            s = float((gx * gx + gy * gy).mean())
        elif "laplacian" in m:
            lap = (-4.0 * q + np.roll(q, 1, 0) + np.roll(q, -1, 0)
                   + np.roll(q, 1, 1) + np.roll(q, -1, 1))
            s = float(lap.var())
        elif "brenner" in m:
            d = q[:, 2:] - q[:, :-2]
            s = float((d * d).mean())
        elif "fft" in m or "freq" in m:
            F = np.abs(np.fft.fftshift(np.fft.fft2(q)))
            h, w = F.shape
            yy, xx = np.ogrid[:h, :w]
            rr = (min(h, w) / 8.0) ** 2
            hi = ((yy - h / 2.0) ** 2 + (xx - w / 2.0) ** 2) > rr
            s = float(F[hi].sum() / (F.sum() + 1e-6))
        else:  # normalized variance
            mu = float(q.mean()) + 1e-6
            s = float(q.var()) / (mu * mu)
        scores.append(s)
    return scores


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


# A rig-exported reference .bin supplies the detector/camera identity
# (hwUnit_Name, DetectorType, ChannelBits) and the per-item device fields
# (e.g. CH2PMTHighVoltage, iShowFlags, eItemType). Without them NIS-E rejects the
# bin as "Incompatible Z Correction and Camera" and the panel stays empty. We
# template off the reference and substitute only ZStack + the laser-power ramp.

_LP_FIELD_CANDIDATES = (
    "CH2LaserPower", "CH1LaserPower", "CH3LaserPower",
    "CH4LaserPower", "CH5LaserPower",
)


def load_bin_template(reference_bin: Path) -> tuple:
    """Parse a rig reference .bin -> (metadata, template_items).

    metadata: root-level scalars (hwUnit_Name, DetectorType, ChannelBits, ...).
    template_items: per-item field dicts in index order (keeps HV / flags).
    """
    import nis_bin_to_csv as _bindec
    parsed = _bindec.NisBinParser(Path(reference_bin).read_bytes()).parse()
    root = next(iter(parsed.values()))
    item_keys = sorted(
        (k for k in root if k.startswith("Item") and k[4:].isdigit()),
        key=lambda k: int(k[4:]),
    )
    template_items = [dict(root[k]) for k in item_keys]
    if not template_items:
        raise ValueError(f"No Item* entries in reference bin: {reference_bin}")
    metadata = {k: v for k, v in root.items()
                if not (k.startswith("Item") and k[4:].isdigit())}
    return metadata, template_items


def _detect_lp_field(template_items: list[dict]) -> str:
    """Pick the laser-power field the reference actually uses (first candidate
    present and non-zero in any item)."""
    for cand in _LP_FIELD_CANDIDATES:
        if any(abs(float(it.get(cand, 0) or 0)) > 1e-9 for it in template_items):
            return cand
    return "CH2LaserPower"


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
                 bin_metadata: dict | None = None,
                 reference_bin: Path | None = None) -> Path:
    """
    Generate a per-spheroid Z Intensity Correction .bin file.

    If reference_bin is given (recommended), the bin inherits that rig export's
    detector/camera metadata + per-item device fields, and only ZStack and the
    laser-power ramp are substituted -- this is what NIS-E actually loads. Without
    it, a bare bin is written (legacy path; NIS-E flags it as "Incompatible Z
    Correction and Camera").
    Returns the path to the written bin file.
    """
    if record.z_centre_um == 0.0:
        raise ValueError(f"Rank {record.rank}: z_centre_um not set — run Step 3 first")

    z_start = round(record.z_centre_um - record.z_half_um, 3)
    z_end   = round(record.z_centre_um + record.z_half_um, 3)

    bin_dir.mkdir(parents=True, exist_ok=True)
    out_path = bin_dir / f"{record.spheroid_id}.bin"

    if reference_bin:
        metadata, template = load_bin_template(reference_bin)
        lp_field = _detect_lp_field(template)
        n = len(template)
        # Use exactly the reference's item count (sidesteps the 19-item
        # ROOT_TAIL_U64S constraint); span z_start..z_end inclusive.
        if n > 1:
            zs = [round(z_start + i * (z_end - z_start) / (n - 1), 3) for i in range(n)]
        else:
            zs = [round((z_start + z_end) / 2.0, 3)]
        items = []
        for i, zv in enumerate(zs):
            it = dict(template[i])          # keep HV / flags / detector fields
            it["ZStack"] = zv
            it[lp_field] = round(_beer_lambert(P0_pct, L_um, max(0.0, zv - z_start)), 1)
            items.append(it)
        from csv_to_nis_bin import build_bin
        out_path.write_bytes(build_bin(metadata, items))
    else:
        items = _build_bin_items(z_start, z_end, record.z_step_um, P0_pct, L_um, ch_field)
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
