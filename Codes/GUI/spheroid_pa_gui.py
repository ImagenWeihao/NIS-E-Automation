"""
spheroid_pa_gui.py  v1.17.2
NIS-E Spheroid PA Pipeline — ND2-native I/O

Pipeline:  Job A ND2 → [parse metadata + detect spheroids]
                      → [Bridge: Beer-Lambert correction]
                      → [Job B: z-centre / PSF filter]
                      → [Job C: PA-ready targets]

ND2 metadata auto-populates: pixel size, Z step, objective,
channel wavelengths, stage positions per tile.

Hard dependencies (install once):
    pip install nd2 numpy scipy

Output: optional final CSV export (ND2 write not supported by nd2 library).
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import csv
import math
import threading
from pathlib import Path
from datetime import datetime

# ── Optional scientific stack ─────────────────────────────────────────────────

_MISSING: list = []

try:
    import nd2 as _nd2lib
except ImportError:
    _MISSING.append("nd2          →  pip install nd2")

try:
    import numpy as np
except ImportError:
    np = None
    _MISSING.append("numpy        →  pip install numpy")

try:
    from scipy import ndimage as _ndi
except ImportError:
    _ndi = None
    _MISSING.append("scipy        →  pip install scipy")

_DEPS_OK = len(_MISSING) == 0

try:
    from matplotlib.figure import Figure as _MplFigure
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg as _FigCanvas
    _MPL_TK_OK = True
except ImportError:
    _MPL_TK_OK = False

try:
    from PIL import Image as _PILImage, ImageTk as _PILImageTk
    _PIL_OK = True
except ImportError:
    _PIL_OK = False


# ── Defaults (from SLIM045 / SLIM033) ────────────────────────────────────────

DEFAULT_LASER_PARAMS = [
    {"wavelength": 850,  "P0_pct": 15.0, "L_um": 165.0},
    {"wavelength": 890,  "P0_pct": 20.0, "L_um": 160.0},
    {"wavelength": 940,  "P0_pct": 30.0, "L_um": 170.0},
    {"wavelength": 1050, "P0_pct": 45.0, "L_um": 194.0},
]

DEFAULT_PA_WAVELENGTH  = 850
DEFAULT_MIN_DIAM_UM    = 50
DEFAULT_MAX_DIAM_UM    = 400
DEFAULT_PSF_CAUTION_UM = 100
DEFAULT_BASE_LOOPS     = 25
DEFAULT_MAX_LOOPS      = 60
DEFAULT_DETECT_SIGMA   = 2.0    # Gaussian blur radius (pixels)
DEFAULT_DETECT_K       = 1.5    # threshold = mean + K * std
MAX_POWER_PCT          = 100.0
# PA activation power hard cap. 80% visibly damaged a spheroid (2026-07-07, well
# A02 sph#9 -- a sharply saturated hot spot appeared in the identical location
# across all 3 independent Pre/Post-PA viz channels, mean intensity fell while
# saturated-pixel count rose 3-14x, consistent with localized burn damage).
# NOTE: pa_trigger.ini's power_pct is NOT read by the step3_zstack_PA JOB itself
# (it's set manually in NIS-E's own Job Wizard) -- this GUI-side cap is a strong
# default/reminder, not a technical enforcement of the real laser power.
MAX_PA_ACTIVATION_POWER_PCT = 30.0

# Sharpness-peak prominence floor for per-spheroid Find Z-Center: a captured wide/coarse
# 1050 stack must show a focus peak at least this fraction above its baseline to be
# trusted as a real focal plane (find_zcenter_all confidence gate). Below it the result
# falls back to the geometric middle -- see the ±25 µm through-focus / edge-snap note.
ZC_PROMINENCE_MIN = 0.35

# NIS-E "mGold" optical-config group. Exact names as configured on the rig (Optical
# Configuration list) -- see NISE010/NISE011 for the activation-vs-visualization
# wavelength rationale (IMPA004, SLIM025/026/031/043).
#
# Activation range (~750-850 nm): IMPA004's 10-nm wavelength sweep to find the best
# 2P photoconversion wavelength per PA-protein/dye; "_PA"/"_PA1"/"_PA2" suffixed
# configs are purpose-built for firing. 850 nm/30% (mGold/PAmKate) is the current
# pipeline default; SLIM031 activated PA-JF646 at 800 nm instead -- the optimum is
# protein/dye-specific, hence the full sweep is offered rather than one fixed value.
PA_ACTIVATION_OC = "850 nm power loop full reso2"    # mGold/PAmKate 2P activation (IMPA004)

# LOCKED TO ONE ENTRY, deliberately. This was an 19-item sweep of 2P wavelengths sitting in
# a free-text dropdown next to the laser-firing button. On 2026-08-07 it held
# "750 nm power loop full reso" when a PA was launched: the run went full-depth across the
# spheroid (37 planes, 100% tissue coverage) before it could be stopped, consuming an
# irreplaceable sample -- and the confirmation dialog said "850 nm" as hardcoded text, so it
# actively concealed the mismatch. The other 18 entries had only ever been selected by
# accident. To sweep wavelengths again, edit PA_ACTIVATION_OC here: a code change is the
# right amount of friction for something that fires a laser.
PA_ACTIVATION_OC_LIST = [PA_ACTIVATION_OC]

# Visualization range (~890-1050 nm): sits above the activation band, so imaging here
# does not trigger photoconversion. Each config targets a specific marker/channel --
# 890 nm = mBeRFP (constitutive T-cell identity marker, SLIM025/026/043); 940 nm =
# PAsfGFP (the photoactivatable T-cell readout, imaged before AND after PA to confirm
# conversion); 1050 nm = spheroid-depth / faded-square re-imaging (SLIM031, pa_validate
# default). 950/970/980/1000 nm are Galvo/Resonance alternates for the same channels.
PA_VIZ_OC_LIST = [
    "1050nm_Galvo_561nm_NDD2_WC",       # default: spheroid depth / PA-validate re-image
    "1050nm_Reso_561nm_NDD2_JL1",
    "1000nm_Galvo_561nm_NDD2_JL1",
    "1000nm_Reso_561nm_NDD2_JL2",
    "940nm_Galvo_488nm_NDD2_JL",         # PAsfGFP before/after (T-cell activation state) --
                                          # confirmed 2026-07-09: this (not the "..._JL2" sibling
                                          # below) has 488 live w/ sane gain; JL2 leaves 640 active
    "940nm_Galvo_488nm_NDD2_JL2",
    "940nm_Reso_488nm_NDD1_JL",
    "940nm_Galvo_561nm_NDD2_BT",
    "890nm_Galvo_600nm_NDD2_BT",         # mBeRFP (T-cell identity marker)
    "950nm_Galvo_488nm_NDD2_JL1",
    "950 nm power loop 32lines Resonance1",
    "970nm_Galvo_488nm_NDD2_JL1",
    "970 nm power loop 128lines Galvo",
    "970nm_Resonance_400Hz",
    "980nm_Galvo_488nm_NDD2_JL1",
    "980nm_Galvo_488nm_NDD2_JL2_band",
    "1040nm_Resonance_440Hz_488nm_NDD",
]

# Catppuccin Mocha
BG       = "#1e1e2e"
BG2      = "#181825"
SURFACE  = "#313244"
SURFACE2 = "#45475a"
SUBTEXT  = "#6c7086"
TEXT     = "#cdd6f4"
TEXT2    = "#bac2de"
LAVENDER = "#b4befe"
MAUVE    = "#cba6f7"
BLUE     = "#89b4fa"
GREEN    = "#a6e3a1"
YELLOW   = "#f9e2af"
PEACH    = "#fab387"
RED      = "#f38ba8"


# ── ND2 loading & metadata ────────────────────────────────────────────────────

def _safe_attr(obj, *attrs, default=None):
    """Safely traverse a chain of attributes."""
    for a in attrs:
        try:
            obj = getattr(obj, a)
        except (AttributeError, TypeError, IndexError):
            return default
    return obj if obj is not None else default


def extract_nd2_metadata(f) -> dict:
    """
    Pull hardware parameters out of an open nd2.ND2File.
    Returns a flat dict of everything useful for populating the GUI.
    """
    meta = {}

    # Voxel / pixel size
    try:
        vox = f.voxel_size()
        meta["pixel_um"]  = round(float(vox.x), 4)
        meta["z_step_um"] = round(float(vox.z), 4) if vox.z and vox.z > 0 else None
    except Exception:
        meta["pixel_um"]  = None
        meta["z_step_um"] = None

    # Image dimensions
    meta["sizes"] = dict(f.sizes)
    meta["shape"]  = list(f.shape)
    meta["n_positions"] = f.sizes.get("P", 1)
    meta["n_z"]         = f.sizes.get("Z", 1)
    meta["n_channels"]  = f.sizes.get("C", 1)

    # Channels
    channels = []
    try:
        for ch in f.metadata.channels:
            info = {}
            info["name"]           = _safe_attr(ch, "channel", "name",           default="")
            info["excitation_nm"]  = _safe_attr(ch, "channel", "excitationLambdaNm", default=None)
            info["emission_nm"]    = _safe_attr(ch, "channel", "emissionLambdaNm",   default=None)
            # Objective (same for all channels, just read from first)
            if not channels:
                meta["objective"]  = _safe_attr(ch, "microscope", "objectiveName",              default="")
                meta["mag"]        = _safe_attr(ch, "microscope", "objectiveMagnification",      default=None)
                meta["na"]         = _safe_attr(ch, "microscope", "objectiveNumericalAperture",  default=None)
            channels.append(info)
    except Exception:
        pass

    meta["channels"] = channels
    return meta


# Fluorescence channel display colors -- matches NIS-E's own per-channel LUT (the
# A1plus Pad detector display), NOT a false-color heatmap. Applied to every nd2-
# derived figure/thumbnail (Captured Z-Stacks viewer, comparison figures) instead
# of plain grayscale, per explicit instruction.
NM_CHANNEL_RGB = {
    "405": (64, 115, 255),
    "488": (26, 255, 26),
    "561": (217, 217, 26),
    "640": (255, 38, 38),
}


def _tint_channel(gray_u8, chan_name):
    """gray_u8: 2D uint8 array already contrast-stretched to 0-255. Returns an RGB
    uint8 array tinted with this channel's NIS-E LUT color (white if unknown)."""
    import numpy as _np
    color = NM_CHANNEL_RGB.get(chan_name, (255, 255, 255))
    rgb = _np.empty(gray_u8.shape + (3,), dtype=_np.uint8)
    g = gray_u8.astype(_np.float32)
    for k, c in enumerate(color):
        rgb[..., k] = (g * (c / 255.0)).astype(_np.uint8)
    return rgb


def _center_crop_to_fov(arr, pixel_um: float, target_fov_um: float):
    """Center-crop a (..., Y, X) stack's last two dims so its physical extent
    matches target_fov_um, given this stack's own pixel_um. No-op if the stack's
    native FOV is already <= target (e.g. 1050nm's OC is natively ~2x more zoomed
    than 890/940nm -- an intentional per-channel difference, not a bug -- so the
    wider channels get cropped down to it for a like-for-like comparison, not the
    other way around)."""
    if not pixel_um or not target_fov_um:
        return arr
    h, w = arr.shape[-2], arr.shape[-1]
    crop_px = int(round(target_fov_um / pixel_um))
    crop_px = max(1, min(crop_px, h, w))
    if crop_px >= h and crop_px >= w:
        return arr
    y0 = (h - crop_px) // 2
    x0 = (w - crop_px) // 2
    return arr[..., y0:y0 + crop_px, x0:x0 + crop_px]


def _get_detection_image(f, pos_idx: int, ch_idx: int) -> "np.ndarray":
    """
    Return a 2-D float32 image for one (position, channel) pair.
    If Z-stack, returns max-projection.
    """
    sizes = f.sizes
    arr   = f.asarray()

    # Build an index dict for the dimension order
    dims  = list(f.sizes.keys())   # e.g. ['P','Z','C','Y','X']

    # Squeeze to (Y, X) for a given P and C
    # We handle the two most common NIS-E output shapes:
    #   (P, Z, C, Y, X)  →  arr[p, :, c, :, :]  then max over Z
    #   (P, C, Y, X)     →  arr[p, c, :, :]
    #   (Z, C, Y, X)     →  arr[:, c, :, :]  then max over Z  (single position)
    #   (C, Y, X)        →  arr[c, :, :]
    def _axis(key):
        return dims.index(key) if key in dims else None

    p_ax = _axis("P")
    z_ax = _axis("Z")
    c_ax = _axis("C")

    img = arr

    # Select position
    if p_ax is not None:
        sl = [slice(None)] * img.ndim
        sl[p_ax] = pos_idx
        img = img[tuple(sl)]
        # axes shift after indexing
        dims2 = [d for d in dims if d != "P"]
        z_ax  = dims2.index("Z") if "Z" in dims2 else None
        c_ax  = dims2.index("C") if "C" in dims2 else None
    else:
        dims2 = dims

    # Select channel
    if c_ax is not None:
        sl = [slice(None)] * img.ndim
        sl[c_ax] = ch_idx
        img = img[tuple(sl)]
        dims2 = [d for d in dims2 if d != "C"]
        z_ax  = dims2.index("Z") if "Z" in dims2 else None

    # Max-project over Z
    if z_ax is not None:
        img = img.max(axis=z_ax)

    return img.astype(np.float32)


def detect_spheroids(image: "np.ndarray", pixel_um: float,
                     min_diam_um: float, max_diam_um: float,
                     sigma: float = DEFAULT_DETECT_SIGMA,
                     threshold_k: float = DEFAULT_DETECT_K) -> list:
    """
    Detect circular bright objects (spheroids) in a 2-D fluorescence image.
    Returns list of dicts: {cx_px, cy_px, diameter_um, area_px}.
    """
    blurred = _ndi.gaussian_filter(image, sigma=sigma)

    mean_v, std_v = blurred.mean(), blurred.std()
    binary  = blurred > (mean_v + threshold_k * std_v)
    filled  = _ndi.binary_fill_holes(binary)

    labeled, n_obj = _ndi.label(filled)

    min_area = math.pi * (min_diam_um / pixel_um / 2) ** 2
    max_area = math.pi * (max_diam_um / pixel_um / 2) ** 2

    results = []
    for lbl in range(1, n_obj + 1):
        mask    = labeled == lbl
        area_px = int(mask.sum())
        if not (min_area <= area_px <= max_area):
            continue
        cy, cx    = _ndi.center_of_mass(mask)
        diam_um   = 2.0 * math.sqrt(area_px / math.pi) * pixel_um
        results.append({
            "cx_px":       cx,
            "cy_px":       cy,
            "diameter_um": round(diam_um, 2),
            "area_px":     area_px,
        })
    return results


def _stage_position(f, seq_idx: int) -> tuple:
    """Return (x_um, y_um, z_um) stage position for frame seq_idx."""
    try:
        fm = f.frame_metadata(seq_idx)
        pos = fm.channels[0].position.stagePositionUm
        return float(pos.x), float(pos.y), float(pos.z)
    except Exception:
        return 0.0, 0.0, 0.0


def _well_label(f, pos_idx: int) -> str:
    """Try to recover a well label from NIS-E XY position names."""
    try:
        exp = f.experiment
        for loop in exp:
            if hasattr(loop, "parameters") and hasattr(loop.parameters, "points"):
                pts = loop.parameters.points
                if pos_idx < len(pts):
                    name = getattr(pts[pos_idx], "name", "")
                    if name:
                        return name
    except Exception:
        pass
    return f"P{pos_idx+1:02d}"


def load_nd2_discovery(nd2_path: Path, ch_idx: int,
                        min_diam_um: float, max_diam_um: float,
                        sigma: float, threshold_k: float) -> tuple:
    """
    Load a NIS-E Job A ND2 file.
    Returns (spheroid_records, hw_metadata).

    spheroid_records: list of dicts ready for run_bridge()
      keys: well, x_um, y_um, z_surface_um, diameter_um
    hw_metadata: dict from extract_nd2_metadata()
    """
    with _nd2lib.ND2File(nd2_path) as f:
        hw = extract_nd2_metadata(f)
        pixel_um = hw["pixel_um"]
        if not pixel_um:
            raise ValueError("Could not read pixel size from ND2 metadata.")

        n_pos = hw["n_positions"]
        records = []

        for p_idx in range(n_pos):
            # Sequential index for this position (assuming single Z plane or we max-project)
            # For (P, Z, C, Y, X) the seq index for position p is p * n_z * n_ch
            # frame_metadata is 0-indexed over all sequences
            seq_idx = p_idx * hw["n_z"] * hw["n_channels"]

            sx, sy, sz = _stage_position(f, seq_idx)
            well       = _well_label(f, p_idx)

            img     = _get_detection_image(f, p_idx, ch_idx)
            objects = detect_spheroids(img, pixel_um, min_diam_um, max_diam_um,
                                       sigma, threshold_k)

            img_h, img_w = img.shape
            for obj in objects:
                # Convert pixel centroid → stage coordinate
                # (0,0) pixel = top-left of tile; stage XY offset from tile centre
                dx_um = (obj["cx_px"] - img_w / 2.0) * pixel_um
                dy_um = (obj["cy_px"] - img_h / 2.0) * pixel_um
                records.append({
                    "well":         well,
                    "x_um":         round(sx + dx_um, 2),
                    "y_um":         round(sy + dy_um, 2),
                    "z_surface_um": round(sz, 2),
                    "diameter_um":  obj["diameter_um"],
                })

    return records, hw


# ── Demo / synthetic data ─────────────────────────────────────────────────────

def generate_demo_spheroids() -> tuple:
    """
    Generate realistic synthetic spheroid records and fake HW metadata.
    Used for testing the pipeline without a real ND2 file.

    Simulates a 96-well plate scan with 2 spheroids per well across 8 wells,
    spanning a range of diameters so that some are PSF-limited, some mid-depth,
    and some fully clear.
    """
    import random
    random.seed(42)

    # Fake hardware metadata (mimics extract_nd2_metadata output)
    hw = {
        "pixel_um":     0.65,
        "z_step_um":    2.0,
        "objective":    "Plan Apo Lambda 20x",
        "mag":          20,
        "na":           0.75,
        "n_positions":  8,
        "n_z":          10,
        "n_channels":   2,
        "sizes":        {"P": 8, "Z": 10, "C": 2, "Y": 512, "X": 512},
        "shape":        [8, 10, 2, 512, 512],
        "channels": [
            {"name": "GFP",  "excitation_nm": 488,  "emission_nm": 525},
            {"name": "DAPI", "excitation_nm": 405,  "emission_nm": 460},
        ],
    }

    # Well layout (2 rows × 4 cols of a 96-well plate, 9 mm pitch)
    wells = ["A1","A2","A3","A4","B1","B2","B3","B4"]
    pitch_um = 9000.0

    records = []
    for w_idx, well in enumerate(wells):
        col = w_idx % 4
        row = w_idx // 4
        stage_x = col * pitch_um + random.uniform(-200, 200)
        stage_y = row * pitch_um + random.uniform(-200, 200)
        stage_z = 380.0 + random.uniform(-30, 30)

        # Two spheroids per well, spread across diameter range
        for offset, diam in enumerate([
            random.uniform(60, 140),    # small/medium — likely OK
            random.uniform(180, 380),   # large — may be PSF-limited
        ]):
            dx = random.uniform(-1500, 1500)
            dy = random.uniform(-1500, 1500)
            records.append({
                "well":         well,
                "x_um":         round(stage_x + dx, 2),
                "y_um":         round(stage_y + dy, 2),
                "z_surface_um": round(stage_z, 2),
                "diameter_um":  round(diam, 2),
            })

    return records, hw


# ── Pipeline stages ───────────────────────────────────────────────────────────

def beer_lambert_power(P0: float, L: float, depth_um: float) -> float:
    return min(P0 * math.exp(depth_um / L), MAX_POWER_PCT)


def run_bridge(raw: list, laser_params: dict, pa_wavelength: int,
               psf_caution_depth: float, base_loops: int, max_loops: int,
               z_step_um: float = 2.0) -> list:
    """
    For each spheroid, emit one record per Z plane so that PA power is
    Beer-Lambert corrected for the exact depth of that plane, not just
    the equatorial depth.

    depth_at_plane  = z_plane_um - z_surface_um   (distance light travels
                                                    through tissue to reach
                                                    that plane from the top)
    pa_power_pct    = P0 * exp(depth_at_plane / L)
    """
    import math as _math
    out = []
    p = laser_params[pa_wavelength]
    for s in raw:
        radius     = s["diameter_um"] / 2.0
        z_middle   = round(s["z_surface_um"] + radius, 2)
        z_start    = round(z_middle - radius, 2)   # = z_surface_um
        z_end      = round(z_middle + radius, 2)
        n_z_planes = _math.ceil(s["diameter_um"] / z_step_um)
        fov_um     = round(s["diameter_um"] * 1.2, 2)

        # Equatorial depth used only for PSF warning and loop scaling
        eq_depth   = radius
        psf_warn   = eq_depth > psf_caution_depth
        loops      = int(base_loops + (eq_depth / 200.0) * (max_loops - base_loops))
        loops      = max(base_loops, min(loops, max_loops))

        for plane_idx in range(n_z_planes):
            z_plane  = round(z_start + plane_idx * z_step_um, 2)
            # Depth from coverslip surface to this plane (Beer-Lambert input)
            depth_here = max(round(z_plane - s["z_surface_um"], 2), 0.0)
            pa_pwr     = beer_lambert_power(p["P0_pct"], p["L_um"], depth_here)
            out.append({
                **s,
                "plane_idx":        plane_idx,
                "z_plane_um":       z_plane,
                "depth_um":         depth_here,
                "z_middle_um":      z_middle,
                "z_start_um":       z_start,
                "z_end_um":         z_end,
                "z_step_um":        z_step_um,
                "n_z_planes":       n_z_planes,
                "fov_um":           fov_um,
                "pa_power_pct":     round(pa_pwr, 1),
                "pa_loops":         loops,
                "pa_wavelength_nm": pa_wavelength,
                "psf_warning":      psf_warn,
            })
    return out


def run_z_center(bridged: list, min_diam: float, max_diam: float) -> list:
    return [t for t in bridged
            if min_diam <= t["diameter_um"] <= max_diam and not t["psf_warning"]]


def run_photoactivation(b_targets: list) -> list:
    return [{**t, "pa_ready": True} for t in b_targets]


# ── CSV helpers ───────────────────────────────────────────────────────────────

# Columns written to the NIS-E JOBS target CSV — in the order a macro reads them.
# Internal pipeline fields (psf_warning, pa_ready, depth_um, z_surface_um) are excluded.
_NISELEMENTS_COLS = [
    "well",
    "x_um", "y_um",
    "z_middle_um", "z_start_um", "z_end_um", "z_step_um", "n_z_planes",
    "plane_idx", "z_plane_um",
    "fov_um",
    "pa_power_pct", "pa_loops", "pa_wavelength_nm",
]


def write_niselements_csv(records: list, path: Path) -> None:
    """Write only the columns a NIS-E JOBS macro needs, in macro-friendly order."""
    if not records:
        return
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=_NISELEMENTS_COLS, extrasaction="ignore")
        w.writeheader()
        w.writerows(records)


def write_csv(records: list, path: Path) -> None:
    """Write all columns (used internally / for debugging)."""
    if not records:
        return
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(records[0].keys()))
        w.writeheader()
        w.writerows(records)


# ── GUI ───────────────────────────────────────────────────────────────────────

class App(tk.Tk):

    PHASES = [("A", "Job A  ND2"), ("Bridge", "Bridge"),
              ("B", "Z-Centre"), ("C", "Photo-act.")]

    def __init__(self):
        super().__init__()
        self.title("SpheroidPA  v1.17.2 — NIS-E Spheroid PA Pipeline")
        self.geometry("1680x880")
        self.minsize(1000, 660)
        self.configure(bg=BG)

        self._raw:       list = []
        self._bridged:   list = []
        self._b_targets: list = []
        self._c_targets: list = []
        self._hw_meta:   dict = {}
        self._laser_vars: list = []
        self._ch_names:  list = []   # populated after ND2 load

        # Pipeline tab state
        self._pl_zc_found: dict = {}   # {sid: find_zcenter result} for the viewer highlight
        self._pl_records:  list = []   # list[SpheroidRecord]
        self._pl_state     = None      # PipelineState | None
        self._pl_anchors:  list = []   # list of anchor result dicts
        self._pl_offset    = (0.0, 0.0)
        self._pl_dashboard = None      # PipelineDashboard | None
        self._pl_mpl_canvas = None     # FigureCanvasTkAgg | None
        self._pl_dash_mosaic = ""      # mosaic path used to build current dashboard

        self._build_ui()

        if not _DEPS_OK:
            self.after(200, self._show_dep_warning)

    # ── Top-level layout ──────────────────────────────────────────────────────

    def _build_ui(self):
        hdr = tk.Frame(self, bg=BG2)
        hdr.pack(fill="x")
        tk.Label(hdr, text="  SpheroidPA  v1.17.2",
                 bg=BG2, fg=MAUVE, font=("Segoe UI", 13, "bold"), pady=8
                 ).pack(side="left")
        tk.Label(hdr, text="Screen → Anchor → Autofocus → Capture  ",
                 bg=BG2, fg=SUBTEXT, font=("Segoe UI", 8)).pack(side="right")

        self._phase_lbls = {}

        style = ttk.Style(self)
        style.theme_use("clam")
        style.configure("TNotebook",     background=BG, borderwidth=0)
        style.configure("TNotebook.Tab", background=SURFACE, foreground=TEXT2,
                        padding=[14, 4], font=("Segoe UI", 9))
        style.map("TNotebook.Tab",
                  background=[("selected", SURFACE2)],
                  foreground=[("selected", MAUVE)])

        nb = ttk.Notebook(self)
        nb.pack(fill="both", expand=True, padx=8, pady=(4, 8))
        self._nb = nb

        tabs = {
            "pipe":    ("  PA Workflow  ",       tk.Frame(nb, bg=BG)),
            "zviewer": ("  Captured Z-Stacks  ", tk.Frame(nb, bg=BG)),
            "log":     ("  Log  ",               tk.Frame(nb, bg=BG)),
        }
        for key, (label, frame) in tabs.items():
            nb.add(frame, text=label)
            setattr(self, f"_tab_{key}", frame)

        self._build_pipeline_tab(self._tab_pipe)
        self._pl_build_zviewer(self._tab_zviewer)
        self._build_log_tab(self._tab_log)

    def _build_status_strip(self):
        strip = tk.Frame(self, bg=BG2, pady=5)
        strip.pack(fill="x", padx=8)
        self._phase_lbls = {}
        tk.Label(strip, text="Pipeline:", bg=BG2, fg=SUBTEXT,
                 font=("Segoe UI", 8)).pack(side="left", padx=(4, 6))
        for i, (key, label) in enumerate(self.PHASES):
            if i:
                tk.Label(strip, text="→", bg=BG2, fg=SUBTEXT,
                         font=("Segoe UI", 10)).pack(side="left")
            lbl = tk.Label(strip, text=f"  {label}  ",
                           bg=SURFACE, fg=SUBTEXT,
                           font=("Segoe UI", 8, "bold"), padx=5, pady=2)
            lbl.pack(side="left", padx=2)
            self._phase_lbls[key] = lbl

    # ── Config tab ────────────────────────────────────────────────────────────

    def _build_config_tab(self, parent):
        canvas = tk.Canvas(parent, bg=BG, highlightthickness=0)
        sb     = ttk.Scrollbar(parent, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=sb.set)
        sb.pack(side="right", fill="y")
        canvas.pack(side="left", fill="both", expand=True)

        inner = tk.Frame(canvas, bg=BG)
        win   = canvas.create_window((0, 0), window=inner, anchor="nw")
        canvas.bind("<Configure>",
                    lambda e: canvas.itemconfigure(win, width=e.width))
        inner.bind("<Configure>",
                   lambda e: canvas.configure(scrollregion=canvas.bbox("all")))

        p = dict(padx=12, pady=3)

        # ── ND2 input ─────────────────────────────────────────────────────────
        self._section(inner, "Input — Job A ND2 File")
        fp = tk.Frame(inner, bg=BG); fp.pack(fill="x", **p)
        self._nd2_path = tk.StringVar()
        self._out_path = tk.StringVar()
        self._nd2_file_row(fp, "Job A ND2 file (input):",    self._nd2_path)
        self._file_row(fp,     "Final output CSV (optional):", self._out_path, save=True)

        # ── ND2 Metadata display ──────────────────────────────────────────────
        self._section(inner, "ND2 Hardware Metadata  (auto-populated on load)")
        self._meta_frame = tk.Frame(inner, bg=BG2, relief="flat", padx=10, pady=6)
        self._meta_frame.pack(fill="x", padx=12, pady=4)
        self._meta_lbl = tk.Label(self._meta_frame,
                                   text="No ND2 loaded — browse and click  Load Metadata",
                                   bg=BG2, fg=SUBTEXT, font=("Consolas", 8),
                                   justify="left", anchor="w")
        self._meta_lbl.pack(fill="x")

        # Load metadata button
        btn_meta = tk.Frame(inner, bg=BG); btn_meta.pack(fill="x", padx=12, pady=(2, 0))
        self._btn(btn_meta, "Load Metadata from ND2",
                  self._load_metadata_thread, BLUE, "#1e1e2e")

        # ── Channel selector ──────────────────────────────────────────────────
        self._section(inner, "Detection Channel")
        ch_row = tk.Frame(inner, bg=BG); ch_row.pack(fill="x", **p)
        tk.Label(ch_row, text="Channel for spheroid detection:",
                 bg=BG, fg=TEXT2, font=("Segoe UI", 9)).pack(side="left", padx=(0, 8))
        self._ch_var = tk.StringVar(value="0 — (load ND2 first)")
        self._ch_menu = ttk.Combobox(ch_row, textvariable=self._ch_var,
                                      state="readonly", width=32,
                                      font=("Segoe UI", 9))
        self._ch_menu.pack(side="left")

        # Detection sensitivity
        det_row = tk.Frame(inner, bg=BG); det_row.pack(fill="x", **p)
        self._sigma_var = tk.StringVar(value=str(DEFAULT_DETECT_SIGMA))
        self._k_var     = tk.StringVar(value=str(DEFAULT_DETECT_K))
        for label, var, tip in [
            ("Gaussian blur σ (px):", self._sigma_var, "Higher = smoother, less noise"),
            ("Threshold k (mean+k·σ):", self._k_var,  "Higher = only brighter objects"),
        ]:
            tk.Label(det_row, text=label, bg=BG, fg=TEXT2,
                     font=("Segoe UI", 9)).pack(side="left", padx=(0, 4))
            tk.Entry(det_row, textvariable=var, width=6,
                     bg=SURFACE, fg=TEXT, insertbackground=TEXT,
                     relief="flat", font=("Segoe UI", 9)).pack(side="left", padx=(0, 16))
            tk.Label(det_row, text=tip, bg=BG, fg=SUBTEXT,
                     font=("Segoe UI", 8)).pack(side="left", padx=(0, 20))

        # ── Beer-Lambert table ────────────────────────────────────────────────
        self._section(inner, "Beer-Lambert Parameters  (editable, auto-fills from ND2)")
        self._build_laser_table(inner)

        # ── Thresholds ────────────────────────────────────────────────────────
        self._section(inner, "Analysis Thresholds")
        tf = tk.Frame(inner, bg=BG); tf.pack(fill="x", **p)
        self._min_d   = tk.StringVar(value=str(DEFAULT_MIN_DIAM_UM))
        self._max_d   = tk.StringVar(value=str(DEFAULT_MAX_DIAM_UM))
        self._psf_d   = tk.StringVar(value=str(DEFAULT_PSF_CAUTION_UM))
        self._pa_wl   = tk.StringVar(value=str(DEFAULT_PA_WAVELENGTH))
        self._base_lp = tk.StringVar(value=str(DEFAULT_BASE_LOOPS))
        self._max_lp  = tk.StringVar(value=str(DEFAULT_MAX_LOOPS))
        fields = [
            ("Min spheroid diameter (µm):", self._min_d),
            ("Max spheroid diameter (µm):", self._max_d),
            ("PSF caution depth (µm):",      self._psf_d),
            ("PA wavelength (nm):",           self._pa_wl),
            ("Base PA loops (surface):",      self._base_lp),
            ("Max PA loops (deep):",          self._max_lp),
        ]
        for i, (lbl, var) in enumerate(fields):
            r, c = divmod(i, 2)
            tk.Label(tf, text=lbl, bg=BG, fg=TEXT2,
                     font=("Segoe UI", 9), anchor="w").grid(
                row=r, column=c*2, sticky="w", padx=(0, 4), pady=2)
            tk.Entry(tf, textvariable=var, width=9, bg=SURFACE, fg=TEXT,
                     insertbackground=TEXT, relief="flat",
                     font=("Segoe UI", 9)).grid(
                row=r, column=c*2+1, sticky="w", padx=(0, 24), pady=2)

        # ── Actions ───────────────────────────────────────────────────────────
        self._section(inner, "Actions")
        btn = tk.Frame(inner, bg=BG); btn.pack(fill="x", padx=12, pady=8)
        self._btn(btn, "▶  Run Full Pipeline",    self._run_pipeline_thread, MAUVE,    "#1e1e2e", bold=True)
        self._btn(btn, "Demo (synthetic data)",  self._run_demo_thread,    BLUE,     "#1e1e2e")
        self._btn(btn, "Export PA Targets CSV",  self._export_final,       GREEN,    "#1e1e2e")
        self._btn(btn, "Clear Log",              self._clear_log,           SURFACE2, TEXT)

    # ── Helper widgets ────────────────────────────────────────────────────────

    def _section(self, parent, text):
        f = tk.Frame(parent, bg=BG); f.pack(fill="x", padx=12, pady=(12, 2))
        tk.Label(f, text=text, bg=BG, fg=MAUVE,
                 font=("Segoe UI", 9, "bold")).pack(side="left")
        tk.Frame(f, bg=SURFACE2, height=1).pack(
            side="left", fill="x", expand=True, padx=(8, 0), pady=5)

    def _nd2_file_row(self, parent, label, var):
        row = tk.Frame(parent, bg=BG); row.pack(fill="x", pady=2)
        tk.Label(row, text=label, width=26, anchor="w",
                 bg=BG, fg=TEXT2, font=("Segoe UI", 9)).pack(side="left")
        tk.Entry(row, textvariable=var, bg=SURFACE, fg=TEXT,
                 insertbackground=TEXT, relief="flat",
                 font=("Segoe UI", 9)).pack(side="left", fill="x",
                                             expand=True, padx=(0, 6))
        self._btn(row, "Browse", lambda v=var: self._browse_nd2(v), SURFACE2, TEXT, side="left")

    def _file_row(self, parent, label, var, save=False):
        row = tk.Frame(parent, bg=BG); row.pack(fill="x", pady=2)
        tk.Label(row, text=label, width=26, anchor="w",
                 bg=BG, fg=TEXT2, font=("Segoe UI", 9)).pack(side="left")
        tk.Entry(row, textvariable=var, bg=SURFACE, fg=TEXT,
                 insertbackground=TEXT, relief="flat",
                 font=("Segoe UI", 9)).pack(side="left", fill="x",
                                             expand=True, padx=(0, 6))
        cmd = (lambda v=var: self._browse_save_csv(v)) if save \
              else (lambda v=var: self._browse_open_any(v))
        self._btn(row, "Browse", cmd, SURFACE2, TEXT, side="left")

    def _browse_nd2(self, var):
        p = filedialog.askopenfilename(
            title="Select Job A ND2 file",
            filetypes=[("Nikon ND2 files", "*.nd2"), ("All files", "*.*")])
        if p:
            var.set(p)

    def _browse_open_any(self, var):
        p = filedialog.askopenfilename(filetypes=[("All files", "*.*")])
        if p: var.set(p)

    def _browse_save_csv(self, var):
        p = filedialog.asksaveasfilename(
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv")])
        if p: var.set(p)

    def _btn(self, parent, text, cmd, bg, fg, bold=False, side="left"):
        tk.Button(parent, text=text, command=cmd,
                  bg=bg, fg=fg, relief="flat",
                  font=("Segoe UI", 9, "bold" if bold else "normal"),
                  padx=12, pady=5, cursor="hand2").pack(side=side, padx=(0, 6))

    def _build_laser_table(self, parent):
        frame = tk.Frame(parent, bg=BG)
        frame.pack(fill="x", padx=12, pady=4)
        headers = ["λ (nm)", "P₀ (%)", "L (µm)", "Power @ 100µm"]
        for j, h in enumerate(headers):
            tk.Label(frame, text=h, bg=SURFACE, fg=LAVENDER,
                     font=("Segoe UI", 8, "bold"),
                     width=18, pady=3, relief="flat").grid(
                row=0, column=j, padx=1, pady=1, sticky="ew")

        self._laser_vars = []
        for i, lp in enumerate(DEFAULT_LASER_PARAMS):
            rv = {k: tk.StringVar(value=str(lp[k]))
                  for k in ("wavelength", "P0_pct", "L_um")}
            row_bg = BG if i % 2 == 0 else BG2

            for j, key in enumerate(("wavelength", "P0_pct", "L_um")):
                state  = "disabled" if key == "wavelength" else "normal"
                fg_clr = BLUE if key == "wavelength" else TEXT
                e = tk.Entry(frame, textvariable=rv[key], width=18,
                             bg=row_bg, fg=fg_clr,
                             disabledforeground=BLUE,
                             insertbackground=TEXT,
                             state=state, relief="flat",
                             font=("Segoe UI", 9), justify="center")
                e.grid(row=i+1, column=j, padx=1, pady=1, sticky="ew")

            pv = tk.StringVar()
            tk.Label(frame, textvariable=pv, bg=row_bg, fg=YELLOW,
                     font=("Segoe UI", 9), justify="center").grid(
                row=i+1, column=3, padx=1, pady=1, sticky="ew")

            def _upd(rv=rv, pv=pv):
                try:
                    pw = beer_lambert_power(float(rv["P0_pct"].get()),
                                            float(rv["L_um"].get()), 100.0)
                    pv.set(f"{pw:.1f}%")
                except ValueError:
                    pv.set("—")

            for key in ("P0_pct", "L_um"):
                rv[key].trace_add("write", lambda *_, f=_upd: f())
            _upd()
            self._laser_vars.append(rv)

    # ── Spheroid result tabs ──────────────────────────────────────────────────

    def _build_all_tab(self, parent):
        tb = tk.Frame(parent, bg=BG); tb.pack(fill="x", padx=8, pady=6)
        self._all_count = tk.Label(tb, text="No data", bg=BG, fg=SUBTEXT,
                                   font=("Segoe UI", 8))
        self._all_count.pack(side="right", padx=6)
        self._build_legend(parent)
        COLS = ("well","x","y","z_mid","depth","diam","fov","power","loops","psf","stage")
        HDRS = ("Well","X (µm)","Y (µm)","Z mid","Depth (µm)","∅ (µm)","FOV (µm)",
                "Power %","Loops","PSF","Stage")
        WIDS = (55, 80, 80, 75, 80, 70, 75, 70, 55, 45, 80)
        self._all_tree = self._make_tree(parent, COLS, HDRS, WIDS)

    def _build_ready_tab(self, parent):
        tb = tk.Frame(parent, bg=BG); tb.pack(fill="x", padx=8, pady=6)
        self._btn(tb, "Remove Selected", self._remove_ready_sel, SURFACE2, TEXT)
        self._rdy_count = tk.Label(tb, text="No data", bg=BG, fg=SUBTEXT,
                                   font=("Segoe UI", 8))
        self._rdy_count.pack(side="right", padx=6)
        COLS = ("well","x","y","z_mid","depth","diam","fov","power","loops","pa_wl")
        HDRS = ("Well","X (µm)","Y (µm)","Z mid","Depth (µm)","∅ (µm)","FOV (µm)",
                "Power %","Loops","λ PA (nm)")
        WIDS = (55, 85, 85, 80, 85, 70, 75, 70, 55, 70)
        self._rdy_tree = self._make_tree(parent, COLS, HDRS, WIDS)

    def _build_legend(self, parent):
        leg = tk.Frame(parent, bg=BG); leg.pack(fill="x", padx=8, pady=(0, 4))
        for clr, txt in ((GREEN,  "✓ OK (<60µm)"),
                          (YELLOW, "~ Mid depth (60–100µm)"),
                          (RED,    "⚠ PSF-limited (>100µm, hardware ceiling per SLIM045)")):
            tk.Label(leg, text="█", fg=clr, bg=BG,
                     font=("Segoe UI", 10)).pack(side="left")
            tk.Label(leg, text=txt + "   ", fg=SUBTEXT, bg=BG,
                     font=("Segoe UI", 8)).pack(side="left")

    def _make_tree(self, parent, cols, hdrs, wids):
        frame = tk.Frame(parent, bg=BG)
        frame.pack(fill="both", expand=True, padx=8, pady=(0, 8))
        tree = ttk.Treeview(frame, columns=cols, show="headings",
                             style="T.Treeview", selectmode="extended")
        style = ttk.Style()
        style.configure("T.Treeview",
                        background=BG, foreground=TEXT,
                        fieldbackground=BG, rowheight=24,
                        font=("Segoe UI", 9))
        style.configure("T.Treeview.Heading",
                        background=SURFACE, foreground=LAVENDER,
                        font=("Segoe UI", 8, "bold"))
        style.map("T.Treeview", background=[("selected", SURFACE2)])
        for col, hdr, w in zip(cols, hdrs, wids):
            tree.heading(col, text=hdr)
            tree.column(col, width=w, anchor="center", minwidth=36)
        tree.tag_configure("ok",  foreground=GREEN)
        tree.tag_configure("mid", foreground=YELLOW)
        tree.tag_configure("psf", foreground=RED, background="#2a1020")
        vsb = ttk.Scrollbar(frame, orient="vertical", command=tree.yview)
        tree.configure(yscrollcommand=vsb.set)
        vsb.pack(side="right", fill="y")
        tree.pack(fill="both", expand=True)
        return tree

    def _build_log_tab(self, parent):
        # ── NIS-E macro dispatcher: drive macros from the GUI ─────────────────
        # Start nis_macro_dispatcher.mac ONCE in NIS-E; these buttons then write a
        # cmd.ini the dispatcher runs in-process (RunMacro), so no .mac has to be
        # hand-loaded per step.
        disp = tk.LabelFrame(parent, text=" NIS-E Macro Dispatcher ", bg=BG, fg=MAUVE,
                             font=("Segoe UI", 9, "bold"), bd=1, relief="groove")
        disp.pack(fill="x", padx=8, pady=(8, 0))
        tk.Label(disp, text="Start nis_macro_dispatcher.mac once in NIS-E (Macro > Run, or flag it "
                            "Run-on-Startup). These run the capture macros in-process — no per-step "
                            "macro loading. (The PA macros have their own cards in Step 4.)",
                 bg=BG, fg=SUBTEXT, font=("Segoe UI", 8), justify="left",
                 wraplength=920).pack(anchor="w", padx=8, pady=(4, 2))
        row = tk.Frame(disp, bg=BG); row.pack(fill="x", padx=8, pady=(0, 4))
        for label, action, aid, color in [
            ("Autofocus", "autofocus", 1, BLUE),
            ("Z-Stack", "zstack", 2, BLUE),
            ("Z-Corr Capture", "zcorrected", 3, BLUE),
        ]:
            self._btn(row, label, lambda a=action, i=aid: self._pl_send_command(a, i),
                      color, "#1e1e2e", side="left")
        # Global abort: drops a one-shot abort.ini every long-running macro polls, and
        # releases the GUI dispatcher lock. Macros stop at their next loop boundary --
        # an in-flight Z-series/Capture still needs Ctrl+Pause in NIS-E for a hard stop.
        self._btn(row, "Abort All", self._pl_abort_all, RED, "#1e1e2e", side="right")
        self._pl_dispatch_status = tk.Label(disp, text="Dispatcher: idle.", bg=BG, fg=SUBTEXT,
                                            font=("Segoe UI", 9), anchor="w", justify="left",
                                            wraplength=920)
        self._pl_dispatch_status.pack(fill="x", padx=8, pady=(0, 6))

        self._log = scrolledtext.ScrolledText(
            parent, bg=BG2, fg=TEXT, font=("Consolas", 9),
            relief="flat", insertbackground=TEXT, state="disabled")
        self._log.pack(fill="both", expand=True, padx=8, pady=8)
        self._log.tag_configure("info",  foreground=BLUE)
        self._log.tag_configure("ok",    foreground=GREEN)
        self._log.tag_configure("warn",  foreground=YELLOW)
        self._log.tag_configure("error", foreground=RED)
        self._log.tag_configure("dim",   foreground=SUBTEXT)

    # ── Metadata loading ──────────────────────────────────────────────────────

    def _load_metadata_thread(self):
        threading.Thread(target=self._load_metadata, daemon=True).start()

    def _load_metadata(self):
        if not _DEPS_OK:
            self._log_line("Missing dependencies — see warning above.", "error"); return

        path_str = self._nd2_path.get().strip()
        if not path_str:
            self._log_line("Set the ND2 file path first.", "warn"); return

        path = Path(path_str)
        if not path.exists():
            self._log_line(f"File not found: {path}", "error"); return

        try:
            with _nd2lib.ND2File(path) as f:
                hw = extract_nd2_metadata(f)
            self._hw_meta = hw
            self.after(0, lambda: self._populate_metadata_ui(hw))
            self._log_line(f"Metadata loaded from: {path.name}", "ok")
        except Exception as e:
            self._log_line(f"ERROR reading ND2 metadata: {e}", "error")

    def _populate_metadata_ui(self, hw: dict):
        """Fill in metadata panel and auto-update channel selector + BL table."""
        # Build human-readable metadata string
        lines = []
        if hw.get("pixel_um"):
            lines.append(f"Pixel size:   {hw['pixel_um']} µm/px")
        if hw.get("z_step_um"):
            lines.append(f"Z step:       {hw['z_step_um']} µm")
        if hw.get("objective"):
            na  = f"  NA {hw['na']}" if hw.get("na") else ""
            mag = f"{hw['mag']}x " if hw.get("mag") else ""
            lines.append(f"Objective:    {mag}{hw['objective']}{na}")

        sz = hw.get("sizes", {})
        lines.append(f"Positions:    {hw.get('n_positions', 1)}  |  "
                     f"Z planes: {hw.get('n_z', 1)}  |  "
                     f"Channels: {hw.get('n_channels', 1)}")

        ch_lines = []
        self._ch_names = []
        for i, ch in enumerate(hw.get("channels", [])):
            name = ch.get("name", f"Ch{i}")
            exc  = ch.get("excitation_nm")
            em   = ch.get("emission_nm")
            exc_s = f"  ex={exc}nm" if exc else ""
            em_s  = f"  em={em}nm"  if em  else ""
            label = f"[{i}]  {name}{exc_s}{em_s}"
            ch_lines.append(label)
            self._ch_names.append(label)
        if ch_lines:
            lines.append("Channels:")
            lines.extend(f"  {c}" for c in ch_lines)

        self._meta_lbl.configure(text="\n".join(lines), fg=TEXT)

        # Update channel selector
        if self._ch_names:
            self._ch_menu["values"] = self._ch_names
            self._ch_menu.current(0)

        # Auto-populate Beer-Lambert wavelengths from channel excitation
        excitations = [ch.get("excitation_nm") for ch in hw.get("channels", [])
                       if ch.get("excitation_nm")]
        if excitations:
            self._log_line(
                f"Detected excitation wavelengths from ND2: {excitations} nm  "
                "— verify Beer-Lambert table below.", "info")

        # Only update PA wavelength if a channel matches an existing Beer-Lambert entry.
        # Imaging channels (488, 405, 561 …) are never PA lasers — don't override the
        # default 850 nm just because they're absent from the Beer-Lambert table.
        bl_wls = {int(rv["wavelength"].get()) for rv in self._laser_vars}
        matched = [wl for wl in excitations if int(wl) in bl_wls]
        if matched:
            self._pa_wl.set(str(int(matched[0])))

    # ── Pipeline execution ────────────────────────────────────────────────────

    def _run_pipeline_thread(self):
        threading.Thread(target=self._run_pipeline, daemon=True).start()

    def _run_pipeline(self):
        if not _DEPS_OK:
            self._log_line(
                "Cannot run: missing packages.\n" +
                "\n".join(f"  {m}" for m in _MISSING), "error")
            return

        self._log_line("─" * 64, "dim")
        self._log_line(f"Pipeline started  {datetime.now():%Y-%m-%d %H:%M:%S}", "info")

        nd2_path = Path(self._nd2_path.get().strip())
        if not nd2_path.exists():
            self._log_line(f"ERROR: ND2 file not found:\n  {nd2_path}", "error")
            self._set_phase("A", "error"); return

        # Parse parameters
        try:
            laser_params = {}
            for rv in self._laser_vars:
                wl = int(rv["wavelength"].get())
                laser_params[wl] = {"P0_pct": float(rv["P0_pct"].get()),
                                    "L_um":   float(rv["L_um"].get())}
            pa_wl       = int(self._pa_wl.get())
            psf_caution = float(self._psf_d.get())
            base_loops  = int(self._base_lp.get())
            max_loops   = int(self._max_lp.get())
            min_diam    = float(self._min_d.get())
            max_diam    = float(self._max_d.get())
            sigma       = float(self._sigma_var.get())
            thresh_k    = float(self._k_var.get())
        except ValueError as e:
            self._log_line(f"ERROR: bad parameter value — {e}", "error"); return

        if pa_wl not in laser_params:
            self._log_line(f"ERROR: PA wavelength {pa_wl} nm not in Beer-Lambert table", "error"); return

        # Resolve selected channel index
        ch_sel = self._ch_var.get()
        try:
            ch_idx = int(ch_sel.split("]")[0].strip("[")) if "]" in ch_sel else 0
        except (ValueError, IndexError):
            ch_idx = 0
        self._log_line(f"Detection channel index: {ch_idx}  ({ch_sel})", "info")

        try:
            # ── Job A: load ND2 + detect spheroids ───────────────────────────
            self._set_phase("A", "running")
            self._log_line(f"Loading ND2: {nd2_path.name}", "info")
            raw, hw = load_nd2_discovery(
                nd2_path, ch_idx, min_diam, max_diam, sigma, thresh_k)
            self._hw_meta = hw
            self.after(0, lambda h=hw: self._populate_metadata_ui(h))
            self._log_line(
                f"Job A: detected {len(raw)} spheroids across "
                f"{hw.get('n_positions', '?')} positions", "ok")
            for r in raw:
                self._log_line(
                    f"  {r['well']:6s}  XY=({r['x_um']:8.1f}, {r['y_um']:8.1f})  "
                    f"Z={r['z_surface_um']:.1f}µm  ∅={r['diameter_um']:.1f}µm", "info")
            self._set_phase("A", "done")

            # ── Bridge ───────────────────────────────────────────────────────
            self._set_phase("Bridge", "running")
            z_step = hw.get("z_step_um") or 2.0
            bridged = run_bridge(raw, laser_params, pa_wl,
                                 psf_caution, base_loops, max_loops, z_step)
            n_psf = sum(1 for t in bridged if t["psf_warning"])
            self._log_line(
                f"Bridge: {len(bridged)} spheroids — "
                f"{len(bridged)-n_psf} clear, {n_psf} PSF-limited", "ok")
            for t in bridged:
                flag = "⚠ PSF" if t["psf_warning"] else "✓"
                self._log_line(
                    f"  {t['well']:6s}  depth={t['depth_um']:6.1f}µm  "
                    f"PA={t['pa_power_pct']:5.1f}% × {t['pa_loops']:2d} loops  {flag}",
                    "warn" if t["psf_warning"] else "ok")
            self._set_phase("Bridge", "done")

            # ── Job B: z-centre filter + auto-export ─────────────────────────
            self._set_phase("B", "running")
            b_targets = run_z_center(bridged, min_diam, max_diam)
            excl = len(bridged) - len(b_targets)
            self._log_line(
                f"Job B (z-centre): {len(b_targets)} candidates "
                f"({excl} excluded — PSF-limited / diameter out of range)", "ok")
            saved_b = self._auto_export_job_b(b_targets, nd2_path)
            self._log_line(f"Job B CSV saved → {saved_b}", "ok")
            self._set_phase("B", "done")

            # ── Job C: PA-ready stamp ─────────────────────────────────────────
            self._set_phase("C", "running")
            c_targets = run_photoactivation(b_targets)
            self._log_line(
                f"Job C (photo-act.): {len(c_targets)} targets PA-ready", "ok")
            self._set_phase("C", "done")

            self._raw       = raw
            self._bridged   = bridged
            self._b_targets = b_targets
            self._c_targets = c_targets

            self.after(0, self._refresh_all_tree)
            self.after(0, self._refresh_ready_tree)
            self._log_line(
                f"Done — {len(c_targets)}/{len(raw)} targets ready for PA", "ok")
            self.after(0, lambda: self._nb.select(2))

        except Exception as e:
            import traceback
            self._log_line(f"ERROR: {e}", "error")
            self._log_line(traceback.format_exc(), "dim")
            for key in ("A", "Bridge", "B", "C"):
                if self._phase_lbls[key].cget("bg") == PEACH:
                    self._set_phase(key, "error")

    # ── Job B export helpers ──────────────────────────────────────────────────

    def _resolve_job_b_path(self, nd2_path: Path = None) -> Path:
        """Return the save path for the Job B CSV.

        Priority:
          1. Whatever is already in the output-path field.
          2. nd2_path's directory with an auto-generated filename.
          3. Current working directory as last fallback.
        """
        out = self._out_path.get().strip()
        if out:
            return Path(out)

        stem = nd2_path.stem if nd2_path else "spheroid_demo"
        base = nd2_path.parent if nd2_path else Path.cwd()
        return base / f"{stem}_jobB_targets.csv"

    def _auto_export_job_b(self, b_targets: list, nd2_path: Path = None) -> Path:
        """Write the NIS-E JOBS target CSV and update the output-path field."""
        save_path = self._resolve_job_b_path(nd2_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        write_niselements_csv(b_targets, save_path)
        self.after(0, lambda p=str(save_path): self._out_path.set(p))
        return save_path

    # ── Tree refresh ──────────────────────────────────────────────────────────

    def _tag_for(self, t: dict) -> str:
        if t["psf_warning"]:    return "psf"
        if t["depth_um"] > 60:  return "mid"
        return "ok"

    def _refresh_all_tree(self):
        self._all_tree.delete(*self._all_tree.get_children())
        b_ids = {id(t) for t in self._b_targets}
        for t in self._bridged:
            stage = ("B+C ✓" if id(t) in b_ids else
                     "Bridge ✓" if not t["psf_warning"] else "PSF ⚠")
            self._all_tree.insert("", "end", tags=(self._tag_for(t),), values=(
                t["well"], t["x_um"], t["y_um"], t["z_middle_um"],
                t["depth_um"], t["diameter_um"], t["fov_um"],
                t["pa_power_pct"], t["pa_loops"],
                "⚠" if t["psf_warning"] else "✓", stage))
        n     = len(self._bridged)
        n_psf = sum(1 for t in self._bridged if t["psf_warning"])
        self._all_count.configure(
            text=f"{n} total  |  {n_psf} PSF-limited  |  {n - n_psf} clear")

    def _refresh_ready_tree(self):
        self._rdy_tree.delete(*self._rdy_tree.get_children())
        for t in self._c_targets:
            self._rdy_tree.insert("", "end", tags=(self._tag_for(t),), values=(
                t["well"], t["x_um"], t["y_um"], t["z_middle_um"],
                t["depth_um"], t["diameter_um"], t["fov_um"],
                t["pa_power_pct"], t["pa_loops"], t["pa_wavelength_nm"]))
        self._rdy_count.configure(
            text=f"{len(self._c_targets)} PA-ready targets")

    def _remove_ready_sel(self):
        sel = self._rdy_tree.selection()
        if not sel: return
        all_items = self._rdy_tree.get_children()
        indices   = {all_items.index(i) for i in sel}
        self._c_targets = [t for i, t in enumerate(self._c_targets)
                           if i not in indices]
        self._refresh_ready_tree()
        self._log_line(f"Removed {len(indices)} selected targets", "warn")

    # ── Demo mode ─────────────────────────────────────────────────────────────

    def _run_demo_thread(self):
        threading.Thread(target=self._run_demo, daemon=True).start()

    def _run_demo(self):
        self._log_line("─" * 64, "dim")
        self._log_line("DEMO MODE — synthetic spheroid data (no ND2 file needed)", "info")

        try:
            pa_wl       = int(self._pa_wl.get())
            psf_caution = float(self._psf_d.get())
            base_loops  = int(self._base_lp.get())
            max_loops   = int(self._max_lp.get())
            min_diam    = float(self._min_d.get())
            max_diam    = float(self._max_d.get())
            laser_params = {}
            for rv in self._laser_vars:
                wl = int(rv["wavelength"].get())
                laser_params[wl] = {"P0_pct": float(rv["P0_pct"].get()),
                                    "L_um":   float(rv["L_um"].get())}
        except ValueError as e:
            self._log_line(f"ERROR: bad parameter — {e}", "error"); return

        if pa_wl not in laser_params:
            self._log_line(f"ERROR: PA wavelength {pa_wl} nm not in Beer-Lambert table", "error"); return

        self._set_phase("A", "running")
        raw, hw = generate_demo_spheroids()
        self._hw_meta = hw
        self.after(0, lambda h=hw: self._populate_metadata_ui(h))
        self._log_line(f"Demo: generated {len(raw)} synthetic spheroids across 8 wells", "ok")
        for r in raw:
            self._log_line(
                f"  {r['well']:4s}  XY=({r['x_um']:8.1f}, {r['y_um']:8.1f})  "
                f"Z={r['z_surface_um']:.1f}µm  ∅={r['diameter_um']:.1f}µm", "info")
        self._set_phase("A", "done")

        self._set_phase("Bridge", "running")
        z_step = hw.get("z_step_um") or 2.0
        bridged = run_bridge(raw, laser_params, pa_wl, psf_caution, base_loops, max_loops, z_step)
        n_psf = sum(1 for t in bridged if t["psf_warning"])
        self._log_line(
            f"Bridge: {len(bridged)} spheroids — {len(bridged)-n_psf} clear, {n_psf} PSF-limited", "ok")
        for t in bridged:
            flag = "⚠ PSF" if t["psf_warning"] else "✓"
            self._log_line(
                f"  {t['well']:4s}  depth={t['depth_um']:6.1f}µm  "
                f"PA={t['pa_power_pct']:5.1f}% × {t['pa_loops']:2d} loops  {flag}",
                "warn" if t["psf_warning"] else "ok")
        self._set_phase("Bridge", "done")

        self._set_phase("B", "running")
        b_targets = run_z_center(bridged, min_diam, max_diam)
        excl = len(bridged) - len(b_targets)
        self._log_line(
            f"Job B (z-centre): {len(b_targets)} candidates "
            f"({excl} excluded — PSF-limited or out of diameter range)", "ok")
        nd2_str = self._nd2_path.get().strip()
        nd2_ref = Path(nd2_str) if nd2_str else None
        saved_b = self._auto_export_job_b(b_targets, nd2_ref)
        self._log_line(f"Job B CSV saved → {saved_b}", "ok")
        self._set_phase("B", "done")

        self._set_phase("C", "running")
        c_targets = run_photoactivation(b_targets)
        self._log_line(f"Job C: {len(c_targets)} targets PA-ready", "ok")
        self._set_phase("C", "done")

        self._raw       = raw
        self._bridged   = bridged
        self._b_targets = b_targets
        self._c_targets = c_targets

        self.after(0, self._refresh_all_tree)
        self.after(0, self._refresh_ready_tree)
        self._log_line(
            f"Demo complete — {len(c_targets)}/{len(raw)} targets ready for PA", "ok")
        self.after(0, lambda: self._nb.select(2))

    # ── Misc ──────────────────────────────────────────────────────────────────

    def _export_final(self):
        if not self._c_targets:
            messagebox.showwarning("No targets", "Run the pipeline first."); return
        path = self._out_path.get().strip() or filedialog.asksaveasfilename(
            defaultextension=".csv", filetypes=[("CSV files", "*.csv")])
        if path:
            write_csv(self._c_targets, Path(path))
            self._log_line(f"Exported {len(self._c_targets)} PA targets → {path}", "ok")

    def _clear_log(self):
        self._log.configure(state="normal")
        self._log.delete("1.0", "end")
        self._log.configure(state="disabled")

    def _set_phase(self, key: str, state: str):
        if key not in self._phase_lbls:
            return
        colors = {"idle": (SURFACE, SUBTEXT), "running": (PEACH, "#1e1e2e"),
                  "done": (GREEN, "#1e1e2e"),  "error":   (RED, "#1e1e2e")}
        bg, fg = colors.get(state, colors["idle"])
        self.after(0, lambda: self._phase_lbls[key].configure(bg=bg, fg=fg))

    def _log_line(self, text: str, tag: str = "info"):
        now = datetime.now()
        stamped = f"[{now.strftime('%H:%M:%S')}] {text}"
        # Persist every log line to a tailing file (full date + tag) so errors can be
        # read headless -- without the GUI -- for troubleshooting. Best-effort: logging
        # must never break the UI, so any failure here is swallowed. Tail it at
        # C:\SpheroidPA\logs\spheroidpa_gui.log.
        try:
            _lp = Path(r"C:\SpheroidPA\logs\spheroidpa_gui.log")
            _lp.parent.mkdir(parents=True, exist_ok=True)
            with _lp.open("a", encoding="utf-8") as _f:
                _f.write(f"[{now.strftime('%Y-%m-%d %H:%M:%S')}] [{tag}] {text}\n")
        except Exception:
            pass
        def _ins():
            self._log.configure(state="normal")
            self._log.insert("end", stamped + "\n", tag)
            self._log.see("end")
            self._log.configure(state="disabled")
        self.after(0, _ins)

    def _show_dep_warning(self):
        msg = ("Missing Python packages — install them once then restart:\n\n" +
               "\n".join(f"  {m}" for m in _MISSING) +
               "\n\nRun in a terminal:\n  pip install nd2 numpy scipy")
        self._log_line(msg, "error")
        messagebox.showwarning("Missing dependencies", msg)


    # ── Pipeline tab ──────────────────────────────────────────────────────────

    def _build_pipeline_tab(self, parent):
        # ── Shared instance variables (all handlers reference these) ─────────
        self._pl_mosaic_path = tk.StringVar()
        # Default to the macro WORK_DIR so Step 3 writes .../work/autofocus and
        # Step 4 writes .../work, matching the daemons out of the box.
        self._pl_out_dir     = tk.StringVar(value=r"S:\Images\Weihao\NISeA\NIS-E-Automation\work")
        self._pl_well_id     = tk.StringVar(value="")
        self._pl_plate_calib = None   # 3-well affine well->stage-XY calib (lazy-loaded from PLATE_CALIB_INI)
        self._pl_n_spheroids = tk.StringVar(value="")
        self._pl_z_centre    = tk.StringVar(value="7680.0")
        self._pl_z_half      = tk.StringVar(value="90.0")
        self._pl_z_step      = tk.StringVar(value="5.0")
        # Find Z-Center (per spheroid): a WIDE/COARSE 1050 search stack (the ±90/5 Capture
        # slab is too thin to see a ~300 um spheroid's through-focus), then peak-detect each
        # spheroid's focal plane. zc_offset is a one-time OC->PA focal ΔZ added to a confident
        # centre when it is written into z_centre.
        self._pl_zc_half     = tk.StringVar(value="150.0")
        self._pl_zc_step     = tk.StringVar(value="15.0")
        self._pl_zc_oc       = tk.StringVar(value="1050nm_Galvo_561nm_NDD2_WC")
        self._pl_zc_offset   = tk.StringVar(value="0.0")
        self._pl_zc_apply    = tk.BooleanVar(value=True)
        self._pl_recenter_flipx = tk.BooleanVar(value=False)
        self._pl_recenter_flipy = tk.BooleanVar(value=False)
        self._pl_recenter_n  = tk.StringVar(value="all")
        self._pl_recentered  = {}
        self._pl_locate_only = tk.BooleanVar(value=False)
        # Step 3 z-stack captures auto-close in NIS-E by default (mirrors Step 4's
        # _pl_keep_open_in_nise): close_after=1 -> the capture macro CloseCurrentDocument(2)
        # after ImageSaveAs, so a multi-spheroid pass doesn't pile up windows. Tick the
        # Step-3 checkbox to KEEP them open (to eyeball each capture).
        self._pl_s3_keep_open = tk.BooleanVar(value=False)
        self._pl_z_rank      = tk.StringVar()
        self._pl_z_nd2_path  = tk.StringVar()
        self._pl_trigger_dir = tk.StringVar(value=r"S:\Images\Weihao\NISeA\NIS-E-Automation\work")
        self._pl_nd2_out_dir = tk.StringVar(value=r"S:\Images\Weihao\NISeA\NIS-E-Automation\work\nd2")
        self._pl_ch_field    = tk.StringVar(value="CH2LaserPower")
        self._pl_P0          = tk.StringVar(value="15.0")
        self._pl_L_um        = tk.StringVar(value="165.0")
        # Beer-Lambert depth compensation: ON -> depth-adaptive .bin ramp (P0/L);
        # OFF (default) -> flat .bin at the fixed Photoactivation Power %.
        self._pl_beer_lambert = tk.BooleanVar(value=False)
        self._pl_ref_bin     = tk.StringVar()
        # Photoactivation (Step 4) -- parameters for the NIS-E JOB step3_zstack_PA.
        self._pl_pa_job      = tk.StringVar(value="step3_zstack_PA_WC")
        self._pl_pa_oc       = tk.StringVar(value=PA_ACTIVATION_OC)
        self._pl_pa_power    = tk.StringVar(value="30")
        self._pl_pa_power.trace_add("write", self._pl_pa_power_clamp)
        # No separate well var: PA Setup's and Job3's Well fields both bind to
        # _pl_well_id (Step 1's "working well"), so setting it once syncs everywhere.
        self._pl_pa_loops    = tk.StringVar(value="70")
        self._pl_pa_zoom     = tk.StringVar(value="8")
        self._pl_pa_dichroic = tk.BooleanVar(value=True)   # True = dichroic OUT
        self._pl_pa_interlock = tk.BooleanVar(value=True)  # True = remove A1 interlock first
        self._pl_pa_a1on     = tk.BooleanVar(value=False)  # True = A1 confirmed powered ON; pa_setup/pa_validate guard aborts if unchecked
        # Where the PA stacks land. The step3_zstack_PA JOB saves THROUGH the global ND
        # Acquisition experiment, so the ND dialog's "Save to File" Path is the real
        # destination -- and it persists across sessions/users. On 2026-08-07 it was still
        # pointing at a previous user's D:\Images S636B\Brandon\... folder, which is why
        # every GUI-created pa/ folder from 0623..0727 is EMPTY. Blank = <out_dir>/pa.
        self._pl_pa_save_dir    = tk.StringVar(value="")
        self._pl_pa_save_prefix = tk.StringVar(value="PA")
        # Path last prepared by _pl_pa_prepare_save_dir (folder created + recorded in
        # pa_trigger.ini), so the log can report where the PA output is expected.
        self._pl_pa_savepath_applied = None
        # NIS-E JOB launch (Jobs_RunJobByName via nis_macro_run_job.mac, action 10).
        self._pl_job_project = tk.StringVar(value="IMAGEN")
        # Neither job has a StringVar / dropdown any more: each Run button passes a
        # LITERAL job name -- Job1 "Step1_Locate_via_scan" (imaging only, ungated) and
        # Job3 "step3_zstack_PA" (gated, fires the 850 nm laser). Hard-wiring the names
        # means no control can be pointed at an unintended job. project stays a var
        # (IMAGEN, editable) because both jobs are launched under the same project.
        # Per-macro "include in Run Pipeline" selections.
        # pa_setup has no tick of its own since the Photoactivation card merged setup and
        # fire (2026-08-08). It mirrors _pl_pa_sel_pajob: pa_setup exists only to arm the
        # rig for the JOB, and the JOB cannot fire without it, so "laser in the sequence,
        # prep not" was never a state worth being able to express -- it fires through the
        # dichroic with the interlock still in. Kept as a separate var so the sequence
        # builder and the plan log line stay readable.
        self._pl_pa_sel_setup    = tk.BooleanVar(value=True)
        # _pl_pa_sel_points removed with the PA Points card (2026-08-07) -- see the removal
        # note in the Step-4 card builder.
        self._pl_dispatch_busy   = False   # one dispatcher command in flight at a time

        # Pre-PA card: baseline viz capture(s) before PA Setup, via the SAME
        # dispatcher z-stack action (id 2) -- one pass per checked OC, each routed
        # into its own nd2/prePA_<tag> subfolder. See SLIM025/026/031/043.
        self._pl_pa_sel_prepa   = tk.BooleanVar(value=True)    # include in Run Pipeline
        # Every card ships ticked (operator request, 2026-08-08): the whole point of Step 4
        # is the full pre-PA -> PA -> post-PA sequence, and a card silently left off reads
        # as a hang rather than a choice. The laser is still gated twice and NOT by this
        # tick: 'A1 powered ON' must be set, and _pl_pa_run_pipeline_thread shows a modal
        # naming the optical config and every point before the worker thread starts.
        self._pl_pa_sel_pajob   = tk.BooleanVar(value=True)    # launch step3_zstack_PA_WC
        self._pl_pa_sel_pajob.trace_add(
            "write", lambda *_: self._pl_pa_sel_setup.set(self._pl_pa_sel_pajob.get()))
        # ON by default: post-PA is imaging only (no laser), and the pipeline already
        # refuses to run it unless the PA genuinely produced output -- so a stray tick
        # cannot manufacture a bogus comparison. Left OFF, a "Run All" ends silently after
        # the PA and looks like a hang, which is how it read on 2026-08-07.
        self._pl_pa_sel_postpa  = tk.BooleanVar(value=True)    # After-PA captures
        # Re-centre off the pre-PA captures before firing. ON by default: the pipeline
        # already takes those captures, they are the exact input the manual Step-3
        # Re-center uses, and forgetting the manual click silently leaves the PA square
        # off-target. 2026-08-07 measured 21.7 um of centring error on a run where it was
        # skipped, versus 2.7-4.1 um on re-centred spheroids -- with a 39.8 um square
        # half-width against an 80 um spheroid radius, that is the difference between
        # activating the middle and clipping the rim.
        self._pl_pa_sel_recenter = tk.BooleanVar(value=True)
        self._pl_prepa_oc_890   = tk.BooleanVar(value=False)   # mBeRFP (T-cell identity)
        self._pl_prepa_oc_940   = tk.BooleanVar(value=True)    # PAsfGFP (before readout)
        self._pl_prepa_oc_1050  = tk.BooleanVar(value=False)   # spheroid depth / faded-square
        self._pl_prepa_locate_only = tk.BooleanVar(value=False)  # z_half=0 -> 1-plane (mirrors Step 3)
        # After-PA Validation: independent channel selection, saves to postPA_<tag>/
        self._pl_postpa_oc_890  = tk.BooleanVar(value=False)
        self._pl_postpa_oc_940  = tk.BooleanVar(value=True)
        self._pl_postpa_oc_1050 = tk.BooleanVar(value=False)
        self._pl_postpa_locate_only = tk.BooleanVar(value=False)
        # Before/after comparability. pa_setup arms the rig for PA (dichroic OUT, zoom to
        # the PA square) and NOTHING in the After-PA path puts it back -- that card
        # dispatches only capture_zstack (action 2), which never touches either. On
        # 2026-08-07 that left zoom 8 + dichroic OUT against a zoom 2 + dichroic IN
        # baseline: a 79.5 um field vs a 318.2 um one, where the difference would read as
        # photoactivation. ON by default: the pre-PA pass snapshots the rig's laser powers
        # and the post-PA pass restores them along with the dichroic and zoom, so the two
        # halves match by construction instead of by the operator remembering two numbers.
        #
        # No longer exposed in the UI (2026-08-08). There is no legitimate reason to
        # capture the post-PA set through the PA light path, and the mismatch is invisible
        # in the images -- it just looks like signal change. Kept as a var so the
        # snapshot/restore code path stays readable and testable.
        self._pl_lightpath_restore = tk.BooleanVar(value=True)
        # Keep each captured ND2 open in NIS-E after it is saved. OFF by default: a
        # multi-spheroid x multi-OC x z-stack run opens dozens of image windows, which
        # clutters NIS-E and can stall the acquisition itself. The GUI's "Captured
        # Z-Stacks" tab is where captures get reviewed, so there is no need to also
        # hold them all open on the rig. Tick only for live debugging.
        self._pl_keep_open_in_nise = tk.BooleanVar(value=False)

        # Initialization card (Step 4, first card): rig known-good defaults dispatched
        # as init_rig (action 8) -> init_trigger.ini. z_centre/z_half/z_step reuse the
        # Step 3 vars and a1_on reuses _pl_pa_a1on. Run Pipeline dispatches it immediately
        # before the pre-PA captures AND again before the post-PA captures, so both halves
        # of the comparison start from the same rig state; the button stays for standalone
        # use.
        self._pl_init_zoom         = tk.StringVar(value="2.0")
        self._pl_init_dichroic_out = tk.BooleanVar(value=False)   # unchecked = dichroic IN
        self._pl_init_save_repoint = tk.BooleanVar(value=True)
        # Not exposed (2026-08-08). init_trigger.ini's key set is the macro's contract, so
        # both keys are still written -- with values that are not the operator's to pick:
        #   power_pct  the 850 nm activation power lives in the OC, not here. The
        #              A1PlusPad wiring that would let the GUI set it is not in place, so
        #              a field for it only invited belief in a control that does nothing.
        #   pfs_on     PFS is always off on this rig; it cannot hold through a dichroic-OUT
        #              PA anyway, and _pl_init_dichroic_toggle already force-cleared it.
        INIT_PFS_ON = False
        self._pl_init_power        = tk.StringVar(value=str(MAX_PA_ACTIVATION_POWER_PCT))
        self._pl_init_pfs          = tk.BooleanVar(value=INIT_PFS_ON)

        # ── Outer layout: sidebar | [ middle (steps+dashboard) | table ]
        # A horizontal paned window holds two draggable panes that never overlap:
        # the step content + live dashboard, and the merged Spheroid State &
        # Dashboard table. The captured-spheroid Z-stack viewer is now its own
        # "Captured Z-Stacks" tab. Initial split is set in _pl_init_sash once the
        # window has a width.
        outer = tk.Frame(parent, bg=BG)
        outer.pack(fill="both", expand=True)

        # Left sidebar (fixed width)
        sidebar = tk.Frame(outer, bg=SURFACE, width=148)
        sidebar.pack(side="left", fill="y")
        sidebar.pack_propagate(False)

        paned = ttk.PanedWindow(outer, orient="horizontal")
        paned.pack(side="left", fill="both", expand=True)
        self._pl_paned = paned

        # Pane 1: step content (top) + live dashboard (fills the rest)
        middle = tk.Frame(paned, bg=BG)
        # Pane 2: merged Spheroid State & Dashboard table
        side_panel = tk.Frame(paned, bg=BG2)
        paned.add(middle, weight=3)
        paned.add(side_panel, weight=4)

        content_host = tk.Frame(middle, bg=BG)
        content_host.pack(side="top", fill="x")

        dash_host = tk.Frame(middle, bg=BG)
        dash_host.pack(side="top", fill="both", expand=True)

        # ── Sidebar step cards ────────────────────────────────────────────────
        self._pl_step_cards = {}
        self._pl_step_dots  = {}

        _step_defs = [
            ("s1", "1", "Screen\nMosaic"),
            ("s2", "2", "Anchor\nOffset"),
            ("s3", "3", "Autofocus\n+ Reg"),
            ("s4", "4", "Generate\n+ Capture"),
        ]

        for key, num, name in _step_defs:
            card = tk.Frame(sidebar, bg=BG2, cursor="hand2", padx=6, pady=10)
            card.pack(fill="x", padx=4, pady=(4, 0))

            num_lbl = tk.Label(card, text=num, bg=BG2, fg=MAUVE,
                               font=("Segoe UI", 20, "bold"))
            num_lbl.pack()
            name_lbl = tk.Label(card, text=name, bg=BG2, fg=TEXT2,
                                font=("Segoe UI", 8))
            name_lbl.pack()
            dot = tk.Label(card, text="*", bg=BG2, fg=SUBTEXT,
                           font=("Segoe UI", 9))
            dot.pack()

            self._pl_step_cards[key] = card
            self._pl_step_dots[key]  = dot

            for w in (card, num_lbl, name_lbl, dot):
                w.bind("<Button-1>", lambda e, k=key: self._pl_show_step(k))

        # ── Step content frames (one per step, swapped by _pl_show_step) ─────
        self._pl_step_frames = {}
        p = dict(padx=12, pady=3)

        # Step 1 content
        f_s1 = tk.Frame(content_host, bg=BG)
        self._pl_step_frames["s1"] = f_s1

        f1 = tk.Frame(f_s1, bg=BG); f1.pack(fill="x", **p)

        # Job1 (Step1_Locate_via_scan, the 10X whole-well mosaic JOB) settings.
        # Well ID here is the CANONICAL working well: PA Setup's and the Job3
        # card's Well fields are bound to this SAME StringVar, so setting it once
        # (e.g. "A02") automatically updates every other job's well field too.
        job1 = tk.LabelFrame(f1, text=" Job1 (X10 Mosaic -- Step1_Locate_via_scan_WC) ",
                             bg=BG, fg=MAUVE, font=("Segoe UI", 9, "bold"),
                             bd=1, relief="groove")
        job1.pack(fill="x", pady=(0, 4))
        row_w = tk.Frame(job1, bg=BG); row_w.pack(fill="x", pady=(2, 4), padx=4)
        tk.Label(row_w, text="Well ID (working well -- shared by all jobs):", width=34, anchor="w",
                 bg=BG, fg=TEXT2, font=("Segoe UI", 9)).pack(side="left")
        tk.Entry(row_w, textvariable=self._pl_well_id, width=14,
                 bg=SURFACE, fg=TEXT, insertbackground=TEXT, relief="flat",
                 font=("Segoe UI", 9)).pack(side="left")

        # NIS-E JOB launch: Job1 (the 10X mosaic) is fired directly via Jobs_RunJobByName
        # under the project below. The job NAME is hard-wired -- no dropdown / jobs-dir
        # picker (see row_j2). Project stays editable (defaults IMAGEN).
        row_j = tk.Frame(job1, bg=BG); row_j.pack(fill="x", pady=(2, 2), padx=4)
        tk.Label(row_j, text="Project:", bg=BG, fg=TEXT2,
                 font=("Segoe UI", 9)).pack(side="left", padx=(0, 3))
        tk.Entry(row_j, textvariable=self._pl_job_project, width=10, bg=SURFACE, fg=TEXT,
                 insertbackground=TEXT, relief="flat", font=("Segoe UI", 9)).pack(side="left")
        row_j2 = tk.Frame(job1, bg=BG); row_j2.pack(fill="x", pady=(0, 4), padx=4)
        # Job1 is HARD-WIRED to the imaging-only mosaic JOB. There is deliberately NO
        # dropdown here: a selectable list could offer step3_zstack_PA, whose JOB fires
        # the 850 nm PA activation laser -- and this Run button is UNGATED (gated=False).
        # We pass the literal job name (not a mutable StringVar) so nothing on this
        # control can ever be pointed at the PA job / fire the laser. Laser firing is
        # confined to the Step 4 Job3 card, which is gated (A1 confirm + confirm dialog).
        self._btn(row_j2, "Run Job1 (10X mosaic)",
                  self._pl_run_job1_at_well, BLUE, "#1e1e2e", side="left")
        tk.Label(row_j2, text="Fires: Step1_Locate_via_scan_WC  (fixed -- imaging only, no laser)",
                 bg=BG, fg=TEXT2, font=("Segoe UI", 9)).pack(side="left", padx=(8, 0))
        # 10X mosaic ND2 -- moved here, directly UNDER Run Job1. After the job runs it
        # auto-fills with the mosaic just saved (newest folder under the Step1 jobs dir);
        # or Browse to a file manually.
        self._nd2_file_row(job1, "10X mosaic ND2:", self._pl_mosaic_path)

        # ── Plate calibration ────────────────────────────────────────────────
        # Run Job1 drives the mosaic by stage COORDINATE (NIS-E WellSelection is not
        # param-drivable -- Index/Position are rebuilt from the UI selection at run time).
        # That needs a well -> stage affine, fit here from >=3 reference wells and persisted
        # to PLATE_CALIB_INI so it survives a GUI restart. Re-run this after any plate re-align.
        self._pl_build_plate_calib_panel(f1)

        row_d = tk.Frame(f1, bg=BG); row_d.pack(fill="x", pady=2)
        tk.Label(row_d, text="Output directory:", width=26, anchor="w",
                 bg=BG, fg=TEXT2, font=("Segoe UI", 9)).pack(side="left")
        tk.Entry(row_d, textvariable=self._pl_out_dir, bg=SURFACE, fg=TEXT,
                 insertbackground=TEXT, relief="flat",
                 font=("Segoe UI", 9)).pack(side="left", fill="x", expand=True, padx=(0, 6))
        self._btn(row_d, "Browse",
                  lambda: self._pl_browse_dir(self._pl_out_dir), SURFACE2, TEXT, side="left")
        row_n = tk.Frame(f1, bg=BG); row_n.pack(fill="x", pady=2)
        tk.Label(row_n, text="Top N spheroids:", width=26, anchor="w",
                 bg=BG, fg=TEXT2, font=("Segoe UI", 9)).pack(side="left")
        tk.Entry(row_n, textvariable=self._pl_n_spheroids, width=6, bg=SURFACE, fg=TEXT,
                 insertbackground=TEXT, relief="flat",
                 font=("Segoe UI", 9)).pack(side="left")
        tk.Label(row_n, text="(leave blank for all)",
                 bg=BG, fg=SUBTEXT, font=("Segoe UI", 8)).pack(side="left", padx=(6, 0))
        btn1 = tk.Frame(f_s1, bg=BG); btn1.pack(fill="x", padx=12, pady=(4, 0))
        self._btn(btn1, "Run Screener", self._pl_run_screener_thread, BLUE, "#1e1e2e")
        # Steps 1-3 are re-entered from disk, not re-run: a snapshot is written after
        # every step that changes the record set, and startup restores it automatically.
        # This button is for restoring by hand after pointing Output dir at an older run.
        self._btn(btn1, "Restore Session", self._pl_restore_pipeline_state, SURFACE2, TEXT)
        self._pl_screen_lbl = tk.Label(f_s1, text="No spheroids screened yet.",
                                        bg=BG, fg=SUBTEXT, font=("Segoe UI", 9), anchor="w")
        self._pl_screen_lbl.pack(fill="x", padx=12)

        # Step 2 content
        # TODO(step2): let the user remove unwanted spheroid IDs from the detected
        # list here (debris / edge / merged or close-pair detections) before
        # anchoring -- e.g. a multi-select or per-row "exclude" so they drop out of
        # records + the dashboard + the triggers, instead of carrying every
        # detection through to capture/PA.
        f_s2 = tk.Frame(content_host, bg=BG)
        self._pl_step_frames["s2"] = f_s2

        tk.Label(f_s2,
                 text="Navigate stage to a spheroid in NIS-E, capture an nd2 in the mosaic's channel (10X/wide field matches best); use 2-3 strong, well-separated spheroids (a 3rd anchor tightens the fit -- see NISE010). Type its rank, browse the nd2.\n"
                      "Typing a rank FORCES that cell as the anchor (NCC is only a cross-check — a disagreement is shown as a warning, never overridden). Leave rank blank to auto-match by NCC.",
                 bg=BG, fg=SUBTEXT, font=("Segoe UI", 8), justify="left").pack(anchor="w", padx=12, pady=(8, 2))

        self._pl_anchor_frames = []
        self._pl_anchor_paths  = []
        self._pl_anchor_ranks  = []
        anchor_outer = tk.Frame(f_s2, bg=BG)
        anchor_outer.pack(fill="x", padx=12, pady=2)
        _anchor_labels = ["Anchor 1 (required):", "Anchor 2 (optional):", "Anchor 3 (optional):"]
        for i in range(len(_anchor_labels)):
            af = tk.Frame(anchor_outer, bg=BG2, relief="flat", pady=4, padx=6)
            af.pack(fill="x", pady=2)
            tk.Label(af, text=_anchor_labels[i], bg=BG2, fg=LAVENDER,
                     font=("Segoe UI", 9, "bold")).grid(row=0, column=0, sticky="w")
            tk.Label(af, text="Rank:", bg=BG2, fg=TEXT2,
                     font=("Segoe UI", 8)).grid(row=0, column=1, padx=(8, 2), sticky="e")
            rk_var = tk.StringVar(value="")
            tk.Entry(af, textvariable=rk_var, width=5,
                     bg=SURFACE, fg=TEXT, insertbackground=TEXT, relief="flat",
                     font=("Segoe UI", 9)).grid(row=0, column=2, padx=(0, 6), sticky="w")
            nd2_var = tk.StringVar()
            tk.Entry(af, textvariable=nd2_var, bg=SURFACE, fg=TEXT,
                     insertbackground=TEXT, relief="flat",
                     font=("Segoe UI", 9)).grid(row=0, column=3, padx=6, sticky="ew")
            af.columnconfigure(3, weight=1)
            self._btn(af, "Browse ND2",
                      lambda v=nd2_var: self._browse_nd2(v), SURFACE2, TEXT, side=None,
                      grid=(0, 4))
            status_lbl = tk.Label(af, text="not verified",
                                   bg=BG2, fg=SUBTEXT, font=("Segoe UI", 8))
            status_lbl.grid(row=1, column=0, columnspan=5, sticky="w", pady=(2, 0))
            self._pl_anchor_paths.append(nd2_var)
            self._pl_anchor_ranks.append(rk_var)
            self._pl_anchor_frames.append((af, status_lbl))

        btn2 = tk.Frame(f_s2, bg=BG); btn2.pack(fill="x", padx=12, pady=(4, 0))
        self._btn(btn2, "Verify Anchors + Apply Offset",
                  self._pl_verify_anchors_thread, MAUVE, "#1e1e2e")
        self._pl_offset_lbl = tk.Label(f_s2, text="Offset: not estimated",
                                        bg=BG, fg=SUBTEXT, font=("Segoe UI", 9), anchor="w")
        self._pl_offset_lbl.pack(fill="x", padx=12)

        # Step 3 content (Z-stack capture triggers)
        f_s3 = tk.Frame(content_host, bg=BG)
        self._pl_step_frames["s3"] = f_s3

        tk.Label(f_s3,
                 text=("NIS-E macro nis_macro_capture_zstack.mac navigates to each spheroid\n"
                       "and runs a Z-stack centred on Middle plane Z, spanning +/- Z half-range\n"
                       "at Z step. These three values are written into each trigger and drive\n"
                       "the macro's stack."),
                 bg=BG, fg=SUBTEXT, font=("Segoe UI", 8), justify="left"
                 ).pack(anchor="w", padx=12, pady=(8, 4))

        z_params = tk.Frame(f_s3, bg=BG); z_params.pack(fill="x", padx=12, pady=2)
        for lbl, var in [("Middle plane Z (um):", self._pl_z_centre),
                          ("Z half-range (um):",   self._pl_z_half),
                          ("Z step (um):",         self._pl_z_step)]:
            tk.Label(z_params, text=lbl, bg=BG, fg=TEXT2,
                     font=("Segoe UI", 9)).pack(side="left", padx=(0, 4))
            tk.Entry(z_params, textvariable=var, width=7, bg=SURFACE, fg=TEXT,
                     insertbackground=TEXT, relief="flat",
                     font=("Segoe UI", 9)).pack(side="left", padx=(0, 16))

        btn3a = tk.Frame(f_s3, bg=BG); btn3a.pack(fill="x", padx=12, pady=(4, 0))
        # Step 3 has exactly TWO actions: Capture (here) and Recenter (btn3b). "Capture"
        # regenerates triggers from the checked rows, pins the 1050 nm depth OC, and runs
        # the z-stack macro -- a full stack by default (ND_RunZSeriesExp, the confocal
        # path), or a single plane if "Locate only" is ticked. Captures land in nd2/.
        self._btn(btn3a, "Capture", self._pl_capture_thread, GREEN, "#1e1e2e")
        tk.Checkbutton(btn3a, text="Locate only (1 plane)", variable=self._pl_locate_only,
                       bg=BG, fg=TEXT2, selectcolor=SURFACE, activebackground=BG,
                       activeforeground=TEXT, font=("Segoe UI", 8)).pack(side="left", padx=(8, 0))
        tk.Checkbutton(btn3a, text="Keep captures open in NIS-E", variable=self._pl_s3_keep_open,
                       bg=BG, fg=TEXT2, selectcolor=SURFACE, activebackground=BG,
                       activeforeground=TEXT, font=("Segoe UI", 8)).pack(side="left", padx=(8, 0))
        # Capture is 2-photon (1050 confocal): it removes the A1 interlock, so it requires
        # this confirmation (shared with Step 4's PA cards -- same _pl_pa_a1on var).
        tk.Checkbutton(btn3a, text="A1 powered ON (2P)", variable=self._pl_pa_a1on,
                       bg=BG, fg="#f7768e", selectcolor=SURFACE, activebackground=BG,
                       activeforeground=TEXT, font=("Segoe UI", 8)).pack(side="left", padx=(8, 0))
        self._pl_af_status_lbl = tk.Label(f_s3, text="Autofocus: not started.",
                                           bg=BG, fg=SUBTEXT, font=("Segoe UI", 9), anchor="w",
                                           justify="left", wraplength=480)
        self._pl_af_status_lbl.pack(fill="x", padx=12)

        btn3b = tk.Frame(f_s3, bg=BG); btn3b.pack(fill="x", padx=12, pady=(6, 0))
        # Recenter: reads the previous 1050 nm captures, re-centers each spheroid, and
        # AUTO-UPDATES the triggers with the corrected coords (then Capture again to
        # verify, or proceed to Step 4). '# to use' limits how many ranks; flip X/Y set
        # the camera-pixel -> stage-um axis sign.
        self._btn(btn3b, "Recenter", self._pl_recenter_thread, MAUVE, "#1e1e2e")
        tk.Label(btn3b, text="  # to use:", bg=BG, fg=TEXT2, font=("Segoe UI", 8)).pack(side="left")
        tk.Entry(btn3b, textvariable=self._pl_recenter_n, width=5, bg=SURFACE, fg=TEXT,
                 insertbackground=TEXT, relief="flat", font=("Segoe UI", 9)).pack(side="left", padx=(2, 6))
        tk.Checkbutton(btn3b, text="flip X", variable=self._pl_recenter_flipx, bg=BG, fg=TEXT2,
                       selectcolor=SURFACE, activebackground=BG, activeforeground=TEXT,
                       font=("Segoe UI", 8)).pack(side="left", padx=(8, 0))
        tk.Checkbutton(btn3b, text="flip Y", variable=self._pl_recenter_flipy, bg=BG, fg=TEXT2,
                       selectcolor=SURFACE, activebackground=BG, activeforeground=TEXT,
                       font=("Segoe UI", 8)).pack(side="left", padx=(4, 0))

        # Find Z-Center card: capture ONE wide/coarse 1050 stack per spheroid (deliberately
        # wider than the ±90/5 Capture slab, whose 50 um span can't see a ~300 um spheroid's
        # through-focus), then peak-detect each spheroid's focal plane and write it to
        # z_centre. Two actions: "Capture Z-Search" (runs on the rig, same 2P/1050 idiom as
        # Capture but wide/coarse) and "Find Centers" (offline peak detection on the captures).
        zc_card = tk.Frame(f_s3, bg=BG2, bd=1, relief="solid")
        zc_card.pack(fill="x", padx=12, pady=(10, 2))
        tk.Label(zc_card, text="Find Z-Center (per spheroid)", bg=BG2, fg=MAUVE,
                 font=("Segoe UI", 9, "bold")).pack(anchor="w", padx=8, pady=(4, 0))
        zc_p = tk.Frame(zc_card, bg=BG2); zc_p.pack(fill="x", padx=8, pady=2)
        for lbl, var in [("Search half-range (um):", self._pl_zc_half),
                         ("Search step (um):",       self._pl_zc_step),
                         ("OC->PA dZ offset (um):",   self._pl_zc_offset)]:
            tk.Label(zc_p, text=lbl, bg=BG2, fg=TEXT2,
                     font=("Segoe UI", 9)).pack(side="left", padx=(0, 4))
            tk.Entry(zc_p, textvariable=var, width=7, bg=SURFACE, fg=TEXT,
                     insertbackground=TEXT, relief="flat",
                     font=("Segoe UI", 9)).pack(side="left", padx=(0, 12))
        zc_oc = tk.Frame(zc_card, bg=BG2); zc_oc.pack(fill="x", padx=8, pady=2)
        tk.Label(zc_oc, text="Search OC:", bg=BG2, fg=TEXT2,
                 font=("Segoe UI", 9)).pack(side="left", padx=(0, 4))
        tk.Entry(zc_oc, textvariable=self._pl_zc_oc, bg=SURFACE, fg=TEXT,
                 insertbackground=TEXT, relief="flat",
                 font=("Segoe UI", 9)).pack(side="left", fill="x", expand=True, padx=(0, 8))
        zc_btn = tk.Frame(zc_card, bg=BG2); zc_btn.pack(fill="x", padx=8, pady=(2, 4))
        self._btn(zc_btn, "Capture Z-Search", self._pl_zc_capture_thread, GREEN, "#1e1e2e")
        self._btn(zc_btn, "Find Centers", self._pl_zc_find_thread, MAUVE, "#1e1e2e")
        tk.Checkbutton(zc_btn, text="Apply found centres to z_centre", variable=self._pl_zc_apply,
                       bg=BG2, fg=TEXT2, selectcolor=SURFACE, activebackground=BG2,
                       activeforeground=TEXT, font=("Segoe UI", 8)).pack(side="left", padx=(8, 0))
        self._pl_zc_status_lbl = tk.Label(zc_card, text="Find Z-Center: idle.", bg=BG2, fg=SUBTEXT,
                                          font=("Segoe UI", 9), anchor="w", justify="left",
                                          wraplength=480)
        self._pl_zc_status_lbl.pack(fill="x", padx=8, pady=(0, 4))

        # Z rank/nd2 widgets kept for fallback _pl_record_z handler (no button exposed)
        self._pl_z_rank_combo = ttk.Combobox(f_s3, textvariable=self._pl_z_rank, width=6,
                                              state="readonly", font=("Segoe UI", 9))
        self._pl_z_status_lbl = tk.Label(f_s3, text="", bg=BG, fg=SUBTEXT,
                                          font=("Segoe UI", 9), anchor="w")

        # Step 4 content
        f_s4 = tk.Frame(content_host, bg=BG)
        self._pl_step_frames["s4"] = f_s4

        s4 = tk.Frame(f_s4, bg=BG); s4.pack(fill="x", padx=12, pady=3)
        for lbl, var in [("Trigger dir:", self._pl_trigger_dir),
                          ("ND2 output dir:", self._pl_nd2_out_dir)]:
            row_ = tk.Frame(s4, bg=BG); row_.pack(fill="x", pady=2)
            tk.Label(row_, text=lbl, width=20, anchor="w",
                     bg=BG, fg=TEXT2, font=("Segoe UI", 9)).pack(side="left")
            tk.Entry(row_, textvariable=var, bg=SURFACE, fg=TEXT,
                     insertbackground=TEXT, relief="flat",
                     font=("Segoe UI", 9)).pack(side="left", fill="x", expand=True, padx=(0, 6))
            self._btn(row_, "Browse",
                      lambda v=var: self._pl_browse_dir(v), SURFACE2, TEXT, side="left")
        ref_row = tk.Frame(s4, bg=BG); ref_row.pack(fill="x", pady=2)
        tk.Label(ref_row, text="Reference .bin:", width=20, anchor="w",
                 bg=BG, fg=TEXT2, font=("Segoe UI", 9)).pack(side="left")
        tk.Entry(ref_row, textvariable=self._pl_ref_bin, bg=SURFACE, fg=TEXT,
                 insertbackground=TEXT, relief="flat",
                 font=("Segoe UI", 9)).pack(side="left", fill="x", expand=True, padx=(0, 6))
        self._btn(ref_row, "Browse",
                  lambda: self._pl_browse_bin(self._pl_ref_bin), SURFACE2, TEXT, side="left")
        tk.Label(s4, text="(rig-exported Z-correction .bin — supplies detector/camera identity "
                          "so NIS-E accepts generated bins)",
                 bg=BG, fg=SUBTEXT, font=("Segoe UI", 8),
                 justify="left", wraplength=480).pack(anchor="w", padx=2)
        # Beer-Lambert depth compensation toggle. ON -> generate the depth-adaptive
        # P0*exp(depth/L) ramp below. OFF (default) -> flat bin at the fixed
        # Photoactivation Power %, and the P0 / L fields are disabled.
        bl_row = tk.Frame(s4, bg=BG); bl_row.pack(fill="x", pady=(4, 0))
        tk.Checkbutton(bl_row, text="Beer-Lambert Intensity Compensation",
                       variable=self._pl_beer_lambert,
                       command=self._pl_toggle_beer_lambert,
                       bg=BG, fg=TEXT, selectcolor=SURFACE, activebackground=BG,
                       activeforeground=TEXT, font=("Segoe UI", 9)).pack(side="left")
        tk.Label(s4, text="OFF: flat bin at P0 (%) (fixed, all planes).  "
                          "ON: depth-adaptive P0*exp(depth/L) ramp (uses P0 and L below).",
                 bg=BG, fg=SUBTEXT, font=("Segoe UI", 8),
                 justify="left", wraplength=480).pack(anchor="w", padx=2)
        lp_row = tk.Frame(s4, bg=BG); lp_row.pack(fill="x", pady=2)
        self._pl_bl_fields = {}
        for key, lbl, var in [("ch", "Laser channel field:", self._pl_ch_field),
                              ("p0", "P0 (%):", self._pl_P0),
                              ("l",  "L (um):", self._pl_L_um)]:
            tk.Label(lp_row, text=lbl, bg=BG, fg=TEXT2,
                     font=("Segoe UI", 9)).pack(side="left", padx=(0, 4))
            e = tk.Entry(lp_row, textvariable=var, width=14, bg=SURFACE, fg=TEXT,
                         insertbackground=TEXT, relief="flat", font=("Segoe UI", 9))
            e.pack(side="left", padx=(0, 14))
            self._pl_bl_fields[key] = e
        self._pl_toggle_beer_lambert()   # apply initial enabled/disabled state
        btn4 = tk.Frame(f_s4, bg=BG); btn4.pack(fill="x", padx=12, pady=(6, 2))
        self._btn(btn4, "Generate All Bins",   self._pl_generate_bins_thread, GREEN, "#1e1e2e")
        self._btn(btn4, "Start Capture Queue", self._pl_capture_queue_thread, MAUVE, "#1e1e2e")
        self._pl_capture_lbl = tk.Label(f_s4, text="Capture: idle",
                                         bg=BG, fg=SUBTEXT, font=("Segoe UI", 9), anchor="w",
                                         justify="left", wraplength=480)
        self._pl_capture_lbl.pack(fill="x", padx=12)

        # ── Right pane: per-step content ──────────────────────────────────────
        # Steps 1-3 show the merged Spheroid State & Dashboard table; Step 4 swaps
        # it for the Photoactivation macro cards (toggled in _pl_show_step). The
        # live dashboard stays in the middle pane for every step.
        self._pl_table_host = tk.Frame(side_panel, bg=BG2)
        self._pl_pa_host    = tk.Frame(side_panel, bg=BG)

        tk.Label(self._pl_table_host, text="Spheroid State & Dashboard", bg=BG2, fg=MAUVE,
                 font=("Segoe UI", 9, "bold")).pack(anchor="w", padx=8, pady=(6, 2))
        tk.Label(self._pl_table_host, text="Click the Use box to include/exclude a spheroid from the "
                              "Step 3 triggers (default: all checked).",
                 bg=BG2, fg=SUBTEXT, font=("Segoe UI", 8), justify="left",
                 wraplength=520).pack(anchor="w", padx=8, pady=(0, 2))

        selrow = tk.Frame(self._pl_table_host, bg=BG2)
        selrow.pack(fill="x", padx=8, pady=(0, 2))
        self._btn(selrow, "Select All / None", self._pl_toggle_all_use, SURFACE2, TEXT)

        tbl_frame = tk.Frame(self._pl_table_host, bg=BG2)
        tbl_frame.pack(fill="both", expand=True, padx=4, pady=(0, 6))

        # Dashboard-styled Treeview: dark field, lavender monospace headings.
        dstyle = ttk.Style(self)
        dstyle.configure("Dash.Treeview",
                         background=BG2, foreground=TEXT, fieldbackground=BG2,
                         rowheight=22, borderwidth=0, font=("Consolas", 9))
        dstyle.configure("Dash.Treeview.Heading",
                         background=SURFACE, foreground=LAVENDER, relief="flat",
                         font=("Consolas", 9, "bold"))
        dstyle.map("Dash.Treeview",
                   background=[("selected", SURFACE2)], foreground=[("selected", TEXT)])
        dstyle.map("Dash.Treeview.Heading", background=[("active", SURFACE2)])

        # Ranks the operator unchecked in the Use column -> excluded from Step 3.
        self._pl_excluded = set()
        cols = ("use", "rank", "id", "status", "mosaic_xy", "verified_xy",
                "z_centre", "diam", "score", "bin")
        self._pl_tree = ttk.Treeview(tbl_frame, columns=cols, show="headings",
                                     height=12, style="Dash.Treeview")
        # Use = click-to-toggle checkbox ([x]/[ ]); the rest mirror the dashboard
        # table. Widths sum to all-visible at the default table-pane width; the
        # horizontal scrollbar covers overflow from long IDs / bin names.
        hdrs = {"use": ("Use", 40), "rank": ("Rank", 44), "id": ("Spheroid ID", 108),
                "status": ("Status", 84), "mosaic_xy": ("Mosaic XY (µm)", 96),
                "verified_xy": ("Verified XY (µm)", 112), "z_centre": ("Z-centre (µm)", 80),
                "diam": ("Diam (µm)", 60), "score": ("Score", 64), "bin": ("Bin file", 120)}
        for c, (h, w_) in hdrs.items():
            self._pl_tree.heading(c, text=h)
            self._pl_tree.column(c, width=w_, minwidth=32, anchor="center", stretch=False)
        # Per-status row colors, reusing the dashboard's _ROW_BG / _STATUS_STYLE.
        try:
            import spheroid_pipeline as _pl
            for st, sty in _pl._STATUS_STYLE.items():
                self._pl_tree.tag_configure(
                    str(st), background=_pl._ROW_BG.get(st, BG2), foreground=sty["fc"])
        except Exception:
            pass
        self._pl_tree.tag_configure("recentered", background="#3a2c52", foreground="#e6d8ff")
        self._pl_tree.bind("<Button-1>", self._pl_toggle_use)
        vsb_tree = ttk.Scrollbar(tbl_frame, orient="vertical", command=self._pl_tree.yview)
        hsb_tree = ttk.Scrollbar(tbl_frame, orient="horizontal", command=self._pl_tree.xview)
        self._pl_tree.configure(yscrollcommand=vsb_tree.set, xscrollcommand=hsb_tree.set)
        vsb_tree.pack(side="right", fill="y")
        hsb_tree.pack(side="bottom", fill="x")
        self._pl_tree.pack(side="left", fill="both", expand=True)

        # ── Right pane (Step 4 only): Photoactivation per-macro dispatcher cards
        #    Shown instead of the table when Step 4 is selected (see _pl_show_step).
        #    Each PA macro is its own card: edit params -> Reload writes them to
        #    pa_trigger.ini, Run dispatches that macro via nis_macro_dispatcher.mac;
        #    tick a card to include it in "Run Pipeline" (checked macros, in order).
        # Scrollable container: the Step-4 cards (Validation, PA Setup, Job3,
        # After-PA Validation) plus the Run Pipeline bar are taller than the pane, so the
        # lower cards sit below the fold. Wrap them in a canvas + vertical scrollbar
        # (mirrors the Config-tab pattern) so everything stays reachable at any window
        # height; the wheel scrolls only while the cursor is over this panel.
        _pa_canvas = tk.Canvas(self._pl_pa_host, bg=BG, highlightthickness=0)
        _pa_sb = ttk.Scrollbar(self._pl_pa_host, orient="vertical", command=_pa_canvas.yview)
        _pa_canvas.configure(yscrollcommand=_pa_sb.set)
        _pa_sb.pack(side="right", fill="y")
        _pa_canvas.pack(side="left", fill="both", expand=True)
        _pa_inner = tk.Frame(_pa_canvas, bg=BG)
        _pa_win = _pa_canvas.create_window((0, 0), window=_pa_inner, anchor="nw")
        _pa_canvas.bind("<Configure>", lambda e: _pa_canvas.itemconfigure(_pa_win, width=e.width))
        _pa_inner.bind("<Configure>",
                       lambda e: _pa_canvas.configure(scrollregion=_pa_canvas.bbox("all")))
        _pa_canvas.bind("<Enter>", lambda e: _pa_canvas.bind_all(
            "<MouseWheel>", lambda ev: _pa_canvas.yview_scroll(int(-ev.delta / 120), "units")))
        _pa_canvas.bind("<Leave>", lambda e: _pa_canvas.unbind_all("<MouseWheel>"))

        pa = tk.LabelFrame(_pa_inner, text=" Photoactivation - NIS-E Macro Dispatcher ",
                           bg=BG, fg=MAUVE, font=("Segoe UI", 9, "bold"),
                           bd=1, relief="groove")
        pa.pack(fill="both", expand=True, padx=4, pady=(6, 6))
        # Kept short and CURRENT. The previous text still described the pre-v1.16 flow --
        # "Run Pipeline ... does NOT fire the PA itself, manually Run step3_zstack_PA in
        # JOBS Explorer" -- which stopped being true when the PA job and post-PA captures
        # joined the sequence. A stale instruction on the laser step is worse than none.
        _hdr = tk.Frame(pa, bg=BG); _hdr.pack(fill="x", padx=8, pady=(3, 4))
        self._info(
            _hdr,
            "Run Pipeline: init -> pre-PA -> re-centre -> PA setup -> 850 nm JOB -> "
            "init -> post-PA. Needs nis_macro_dispatcher.mac running in NIS-E.",
            "Every card runs in the order shown, waiting for each to report back. Untick a "
            "card's header to leave it out.\n\n"
            "The 850 nm JOB fires the laser. It is gated twice: 'A1 powered ON' must be "
            "ticked, and a confirm dialog naming the optical config and every point appears "
            "before the sequence starts. Post-PA captures run only if the JOB actually "
            "produced output -- a post-PA set with no PA in between is a mislabelled copy "
            "of the baseline.\n\n"
            "Initialization runs twice, once before each capture phase, so the before and "
            "after sets are acquired from the same rig state.\n\n"
            "Viz OCs: 890 = mBeRFP (T-cell identity), 940 = PAsfGFP (the before/after "
            "readout), 1050 = spheroid depth / faded-square re-image. The activation OC is "
            "locked to 850 nm.",
            bg=BG, wrap=520)

        # Card -1 -- Initialization (dispatcher init_rig action, id 8): normalize the rig
        # + GUI to known-good defaults (confocal zoom=2, z-defaults, dichroic IN,
        # Save-to-File -> pipeline nd2_dir). Writes init_trigger.ini and dispatches
        # nis_macro_init_rig.mac. First card in Step 4.
        #
        # Run Pipeline dispatches this immediately before the pre-PA captures AND again
        # before the post-PA captures (2026-08-08) -- so it is NOT a standalone card any
        # more, and there is no tick to include it: both halves of a before/after set must
        # start from the same rig state or the comparison is not one. The button stays for
        # running it on its own.
        #
        # PA power % and PFS ON were removed: the activation power lives in the 850 nm OC
        # (the A1PlusPad wiring that would let the GUI set it is not in place), and PFS is
        # always off on this rig. Both keys are still written to init_trigger.ini -- see
        # the var block for why.
        init = tk.Frame(pa, bg=BG2, bd=1, relief="solid"); init.pack(fill="x", padx=4, pady=3)
        hdri = tk.Frame(init, bg=BG2); hdri.pack(fill="x", padx=6, pady=(3, 0))
        tk.Label(hdri, text="Initialization  (rig known-good defaults)", bg=BG2, fg=MAUVE,
                 font=("Segoe UI", 9, "bold")).pack(side="left")
        tk.Label(hdri, text="  runs automatically before pre-PA and post-PA", bg=BG2,
                 fg=GREEN, font=("Segoe UI", 8)).pack(side="left")
        bodyi = tk.Frame(init, bg=BG2); bodyi.pack(fill="x", padx=24, pady=(0, 4))
        ri = tk.Frame(bodyi, bg=BG2); ri.pack(fill="x", pady=(2, 0))
        for lbl, var, w in [("Confocal zoom:", self._pl_init_zoom, 5),
                            ("z_centre:", self._pl_z_centre, 8), ("z_half:", self._pl_z_half, 6),
                            ("z_step:", self._pl_z_step, 6)]:
            tk.Label(ri, text=lbl, bg=BG2, fg=TEXT2, font=("Segoe UI", 9)).pack(side="left", padx=(0, 3))
            tk.Entry(ri, textvariable=var, width=w, bg=SURFACE, fg=TEXT, insertbackground=TEXT,
                     relief="flat", font=("Segoe UI", 9)).pack(side="left", padx=(0, 8))
        ri3 = tk.Frame(bodyi, bg=BG2); ri3.pack(fill="x", pady=(2, 0))
        tk.Checkbutton(ri3, text="Dichroic OUT", variable=self._pl_init_dichroic_out, bg=BG2, fg=TEXT2,
                       selectcolor=SURFACE, activebackground=BG2, activeforeground=TEXT,
                       font=("Segoe UI", 9)).pack(side="left")
        tk.Checkbutton(ri3, text="Re-point Save-to-File to pipeline nd2_dir",
                       variable=self._pl_init_save_repoint, bg=BG2, fg=TEXT2, selectcolor=SURFACE,
                       activebackground=BG2, activeforeground=TEXT,
                       font=("Segoe UI", 9)).pack(side="left", padx=(10, 0))
        tk.Checkbutton(ri3, text="A1 powered ON", variable=self._pl_pa_a1on, bg=BG2, fg="#f7768e",
                       selectcolor=SURFACE, activebackground=BG2, activeforeground=TEXT,
                       font=("Segoe UI", 9, "bold")).pack(side="left", padx=(10, 0))
        btn_init = tk.Frame(init, bg=BG2); btn_init.pack(fill="x", padx=24, pady=(2, 4))
        self._btn(btn_init, "Initialize Rig", self._pl_init_rig, BLUE, "#1e1e2e", side="left")
        self._info(bodyi,
                   "Sets zoom / save-path / z-defaults / dichroic to known-good.",
                   "Dispatched as init_rig (action 8) via init_trigger.ini.\n\n"
                   "Run Pipeline runs it twice -- once before the pre-PA captures and again "
                   "before the post-PA captures -- so both halves of the before/after set "
                   "are acquired from the same rig state. PA Setup leaves the rig at "
                   "dichroic OUT and the PA zoom; without this the post-PA pass would "
                   "image a 79.5 um field against a 318.2 um baseline, and the difference "
                   "would read as photoactivation.\n\n"
                   "PFS stays off: it cannot hold focus through a dichroic-OUT PA. PA power "
                   "is not set here -- it lives in the 850 nm optical config.\n\n"
                   "Detector gains and per-OC re-save are rig-side, handled in the macro.")

        # Card 0 -- Validation (dispatcher z-stack action, id 2): before- AND after-PA
        # multi-channel capture, one pass per checked OC, each routed into its own
        # nd2/prePA_<tag> subfolder (SLIM025/026/031/043 before/after protocol). Run it
        # once before PA Setup (baseline) and again after the JOB (post-PA) -- same
        # card, same mechanism; replaces the old single-channel 1050nm-only PA Validate
        # card (removed 2026-07-07). Internal names keep the "_pl_prepa_*" prefix.
        card, body = self._pl_pa_card(pa, "Before-PA Validation  (baseline multi-channel capture -> nd2/prePA_<tag>)",
                                      self._pl_pa_sel_prepa)
        self._info(body,
                   "Baseline z-stack per checked channel. Positions come from Step 3.",
                   "Channel selection only -- position and Z come from whatever "
                   "af_trigger_NN.ini files already exist (generate them in Step 3, or write "
                   "them by hand).\n\n"
                   "Each checked OC injects oc= into every existing trigger and runs one "
                   "full z-stack pass via the dispatcher (nis_macro_capture_zstack.mac, "
                   "action 2) into its own nd2/prePA_<tag> subfolder.\n\n"
                   "The After-PA Validation card captures the matching post-JOB set into "
                   "nd2/postPA_<tag>, at the same positions and Z.")
        r = tk.Frame(body, bg=BG2); r.pack(fill="x", pady=(2, 0))
        tk.Checkbutton(r, text="890nm (mBeRFP - T-cell identity)", variable=self._pl_prepa_oc_890,
                       bg=BG2, fg=TEXT2, selectcolor=SURFACE, activebackground=BG2,
                       activeforeground=TEXT, font=("Segoe UI", 9)).pack(anchor="w")
        tk.Checkbutton(r, text="940nm (PAsfGFP - before readout)", variable=self._pl_prepa_oc_940,
                       bg=BG2, fg=TEXT2, selectcolor=SURFACE, activebackground=BG2,
                       activeforeground=TEXT, font=("Segoe UI", 9)).pack(anchor="w")
        tk.Checkbutton(r, text="1050nm (spheroid depth / faded-square)", variable=self._pl_prepa_oc_1050,
                       bg=BG2, fg=TEXT2, selectcolor=SURFACE, activebackground=BG2,
                       activeforeground=TEXT, font=("Segoe UI", 9)).pack(anchor="w")
        btn_prepa = tk.Frame(card, bg=BG2); btn_prepa.pack(fill="x", padx=24, pady=(2, 6))
        self._btn(btn_prepa, "Run Validation Captures", self._pl_prepa_run_thread, BLUE, "#1e1e2e",
                  side="left")
        tk.Checkbutton(btn_prepa, text="Locate only (1 plane)", variable=self._pl_prepa_locate_only,
                       bg=BG2, fg=TEXT2, selectcolor=SURFACE, activebackground=BG2,
                       activeforeground=TEXT, font=("Segoe UI", 8)).pack(side="left", padx=(8, 0))
        # Applies to EVERY capture dispatched from Step 4 (both Validation phases), not
        # just this card -- the macro closes each ND2 right after ImageSaveAs unless this
        # is ticked. Left unticked, NIS-E stays clean during long multi-channel runs.
        tk.Checkbutton(btn_prepa, text="Keep captured images open in NIS-E (debug; may stall acquisition)",
                       variable=self._pl_keep_open_in_nise,
                       bg=BG2, fg=TEXT2, selectcolor=SURFACE, activebackground=BG2,
                       activeforeground=TEXT, font=("Segoe UI", 8)).pack(side="left", padx=(12, 0))

        # Card 1 -- Photoactivation: rig prep (pa_setup, action 4) AND the 850 nm JOB
        # launch (run_job, action 10), merged into one card 2026-08-08.
        #
        # They were two cards with two ticks, but never two decisions: pa_setup exists only
        # to arm the rig for the JOB, and the JOB cannot fire without it. Two ticks meant a
        # state where the laser was in the sequence and its prep was not -- which fires
        # through the dichroic with the interlock still in. The single header tick drives
        # both (see the _pl_pa_sel_setup mirror below).
        #
        # Trimmed to what actually reaches the hardware. Power %, Zoom and Loops were
        # removed: the A1PlusPad wiring that would let the GUI drive the 850 nm scan is not
        # in place, so those three fields never changed the OC -- the JOB runs entirely on
        # "850 nm power loop full reso2"'s own settings. Leaving editable boxes in front of
        # an operator is worse than having none: it invites setting a power and believing
        # it took. Removing A1 interlock and Dichroic OUT DO work, so they stay.
        #
        # RED throughout, unlike every other card: ticking this one puts the laser in the
        # sequence.
        card = tk.Frame(pa, bg=BG2, bd=1, relief="solid"); card.pack(fill="x", padx=4, pady=3)
        hdrp = tk.Frame(card, bg=BG2); hdrp.pack(fill="x", padx=6, pady=(3, 0))
        tk.Checkbutton(hdrp, text="Photoactivation  (setup + launch step3_zstack_PA_WC - FIRES 850 nm LASER)",
                       variable=self._pl_pa_sel_pajob, bg=BG2, fg=RED,
                       selectcolor=SURFACE, activebackground=BG2, activeforeground=TEXT,
                       font=("Segoe UI", 9, "bold")).pack(side="left")
        body = tk.Frame(card, bg=BG2); body.pack(fill="x", padx=24, pady=(0, 2))
        r = tk.Frame(body, bg=BG2); r.pack(fill="x")
        tk.Label(r, text="Activation OC:", bg=BG2, fg=TEXT2, font=("Segoe UI", 9)).pack(side="left", padx=(0, 3))
        # readonly: the list holds ONE entry, and this also blocks free-text entry -- an
        # editable Combobox would let any OC name be typed straight into a laser-firing
        # parameter, which is the hole the single-entry list is meant to close.
        ttk.Combobox(r, textvariable=self._pl_pa_oc, width=30, font=("Segoe UI", 9),
                     state="readonly", values=PA_ACTIVATION_OC_LIST).pack(side="left")
        tk.Label(r, text="(locked)", bg=BG2, fg=SUBTEXT,
                 font=("Segoe UI", 8)).pack(side="left", padx=(6, 0))
        tk.Label(r, text="  Well:", bg=BG2, fg=TEXT2, font=("Segoe UI", 9)).pack(side="left", padx=(8, 3))
        tk.Entry(r, textvariable=self._pl_well_id, width=5, bg=SURFACE, fg=TEXT,
                 insertbackground=TEXT, relief="flat", font=("Segoe UI", 9)).pack(side="left")
        r = tk.Frame(body, bg=BG2); r.pack(fill="x", pady=(2, 0))
        tk.Checkbutton(r, text="Dichroic OUT", variable=self._pl_pa_dichroic, bg=BG2, fg=TEXT2,
                       selectcolor=SURFACE, activebackground=BG2, activeforeground=TEXT,
                       font=("Segoe UI", 9)).pack(side="left")
        tk.Checkbutton(r, text="Remove A1 interlock", variable=self._pl_pa_interlock, bg=BG2, fg=TEXT2,
                       selectcolor=SURFACE, activebackground=BG2, activeforeground=TEXT,
                       font=("Segoe UI", 9)).pack(side="left", padx=(10, 0))
        # A1-present guard: pa_setup/pa_validate abort (no A1 calls) unless this is ticked.
        tk.Checkbutton(r, text="A1 powered ON", variable=self._pl_pa_a1on, bg=BG2, fg="#f7768e",
                       selectcolor=SURFACE, activebackground=BG2, activeforeground=TEXT,
                       font=("Segoe UI", 9, "bold")).pack(side="left", padx=(10, 0))
        self._info(body,
                   "Arms the rig, then fires the 850 nm laser at Step 3's trigger point.",
                   "SETUP. Removing the A1 interlock and swinging the dichroic OUT are the "
                   "two things the setup step does to the hardware; both are required "
                   "before the JOB can fire.\n\n"
                   "It does NOT set the activation power, zoom or loop count. Those live "
                   "in the '850 nm power loop full reso2' optical config and are edited in "
                   "NIS-E. The GUI fields for them were removed because the A1PlusPad "
                   "wiring that would make them effective is not in place -- they read as "
                   "controls but changed nothing.\n\n"
                   "JOB. Hard-wired to step3_zstack_PA_WC -- deliberately no dropdown, so "
                   "this gated button can never be pointed at a different job. The point "
                   "and well are set at launch via Jobs_RunJobInitParam from the current "
                   "af_trigger, so it fires on exactly the supplied coordinates rather "
                   "than whatever PredefinedPoint the job has stored (which has been seen "
                   "15 mm off target, on plastic between wells).\n\n"
                   "The JOB's own z-stack depth comes from the job, not from here -- the "
                   "GUI cannot read ZStackDefinition, so check it in NIS-E before firing. "
                   "The completion line reports the actual file count.\n\n"
                   "Firing is gated twice: 'A1 powered ON' must be ticked, and a confirm "
                   "dialog shows the optical config and every point before anything "
                   "happens.\n\n"
                   "Well is the canonical working well, shared with Step 1: set it once "
                   "and it follows everywhere.",
                   fg=RED)

        # ── PA save folder ────────────────────────────────────────────────────
        rsd = tk.Frame(body, bg=BG2); rsd.pack(fill="x", pady=(6, 0))
        tk.Label(rsd, text="PA folder:", bg=BG2, fg=TEXT2,
                 font=("Segoe UI", 9)).pack(side="left", padx=(0, 4))
        tk.Entry(rsd, textvariable=self._pl_pa_save_dir, bg=SURFACE, fg=TEXT,
                 insertbackground=TEXT, relief="flat", font=("Segoe UI", 9)
                 ).pack(side="left", fill="x", expand=True, padx=(0, 6))
        self._btn(rsd, "Browse", lambda: self._pl_browse_dir(self._pl_pa_save_dir),
                  SURFACE2, TEXT, side="left")
        self._info(
            body,
            "Defaults to <Output directory>/pa, created automatically and sent to the job.",
            "The folder is created if missing and follows the Output directory, so it moves "
            "with the run folder. A path you type yourself is never overwritten.\n\n"
            "It is sent to the JOB at launch as StoreToFsOnly.Folder -- the 'Alternative "
            "Storage Location' task -- so this box IS the destination; there is no second "
            "place to keep in sync. Previously that path lived only on the job and kept "
            "whatever the last session set, which is how PA output once landed in another "
            "user's folder.\n\n"
            "Note this is NOT NIS-E's ND 'Save to File'. The JOB ignores that dialog; an "
            "earlier build pushed the path there and merely made every validation capture "
            "drop a duplicate copy into pa/.")
        self._pl_pa_save_lbl = tk.Label(body, text="", bg=BG2, fg=YELLOW,
                                        font=("Segoe UI", 8), anchor="w",
                                        justify="left", wraplength=460)
        self._pl_pa_save_lbl.pack(fill="x")

        # Setup and fire, in the order the pipeline runs them. Reload writes
        # pa_trigger.ini; Run Setup dispatches pa_setup alone (useful when re-arming after
        # a manual NIS-E poke); Run Photoactivation is the gated laser button.
        rjb = tk.Frame(card, bg=BG2); rjb.pack(fill="x", padx=24, pady=(6, 6))
        self._btn(rjb, "Reload", lambda: self._pl_pa_reload("pa_setup"),
                  SURFACE2, TEXT, side="left")
        self._btn(rjb, "Run Setup", lambda: self._pl_pa_run_macro("pa_setup", 4),
                  BLUE, "#1e1e2e", side="left")
        self._btn(rjb, "Run Photoactivation",
                  self._pl_run_pa_job, RED, "#1e1e2e", side="left")

        # Job3 card REMOVED 2026-08-08. It held a duplicate Well box (already on PA Setup,
        # same _pl_well_id) and a "Watch for PA done" button that polled pa_done.ini for an
        # hour. Both are obsolete: the JOB is launched by the pipeline, not by hand in JOBS
        # Explorer, and the pipeline waits on the PA's own output and logs the result --
        # "PA job -> complete (N file(s), output quiet for 90s)". pa_done.ini is not even
        # written on the dispatcher path, so the watcher would have timed out on a run that
        # succeeded. See _pl_pa_wait_for_job.

        # PA Points card REMOVED 2026-08-07 -- discontinued, not merely hidden.
        #
        # It built an ND multipoint (ND_ResetMultipointExp / ND_AppendMultipointPoint) for
        # step3_zstack_PA to "Import Point Set from ND". That route is gone: the job's
        # NDToPointSet task was deleted, and the PA point now comes from PredefinedPoints
        # set by Jobs_RunJobInitParam at launch (see _pl_pa_point_params). The job's POINT
        # LOOP iterates PredefinedPoints.Positions and never consults the ND multipoint,
        # so the card had no effect on where the laser fires.
        #
        # It was also actively dangerous to leave: its completion popup told the operator
        # "In step3_zstack_PA: Import Point Set from ND, then Run" -- and following that
        # instruction OVERWRITES the parameterised point with the ND list, which is exactly
        # the failure deleting NDToPointSet was meant to prevent.
        #
        # nis_macro_pa_points.mac and dispatcher action 5 are left in place (harmless, and
        # still useful for a manual ND multipoint), but nothing in the GUI dispatches them.

        # The standalone 'Photoactivation JOB' card was merged into the Photoactivation
        # card above on 2026-08-08 -- setup and fire were never two decisions. Its PA
        # folder box, its (i) text and its Run Photoactivation button all live there now.

        # Card 3 -- After-PA Validation (dispatcher z-stack action, id 2): the SAME
        # multi-channel capture as the before-PA Validation card, but its passes save
        # to nd2/postPA_<tag> (not prePA_<tag>) so the after-PA images no longer collide
        # with the baseline -- auto-produces the Before_pa/After_pa layout the comparison
        # figures expect. Run this AFTER the step3_zstack_PA JOB. Manual card (not in the
        # Run Pipeline sequence, which stops before the laser fires).
        apa = tk.Frame(pa, bg=BG2, bd=1, relief="solid"); apa.pack(fill="x", padx=4, pady=3)
        hdra = tk.Frame(apa, bg=BG2); hdra.pack(fill="x", padx=6, pady=(3, 0))
        tk.Checkbutton(hdra, text="After-PA Validation  (post-Job capture -> nd2/postPA_<tag>)",
                       variable=self._pl_pa_sel_postpa, bg=BG2, fg=MAUVE,
                       selectcolor=SURFACE, activebackground=BG2, activeforeground=TEXT,
                       font=("Segoe UI", 9, "bold")).pack(side="left")
        bodya = tk.Frame(apa, bg=BG2); bodya.pack(fill="x", padx=24, pady=(0, 4))
        self._info(bodya,
                   "Same positions and Z as the baseline, saved to nd2/postPA_<tag>.",
                   "Runs AFTER the 850 nm JOB, from Step 3's checked rows -- the same "
                   "positions and Z geometry as the Before-PA Validation, into a separate "
                   "subfolder so the two sets never collide.\n\n"
                   "In a pipeline these captures run ONLY if the PA job actually produced "
                   "output. Without a PA in between they would be a second copy of the "
                   "baseline wearing a post-PA label, which silently corrupts the "
                   "comparison rather than failing visibly.")
        ra = tk.Frame(bodya, bg=BG2); ra.pack(fill="x", pady=(2, 0))
        tk.Checkbutton(ra, text="890nm (mBeRFP - T-cell identity)", variable=self._pl_postpa_oc_890,
                       bg=BG2, fg=TEXT2, selectcolor=SURFACE, activebackground=BG2,
                       activeforeground=TEXT, font=("Segoe UI", 9)).pack(anchor="w")
        tk.Checkbutton(ra, text="940nm (PAsfGFP - after readout)", variable=self._pl_postpa_oc_940,
                       bg=BG2, fg=TEXT2, selectcolor=SURFACE, activebackground=BG2,
                       activeforeground=TEXT, font=("Segoe UI", 9)).pack(anchor="w")
        tk.Checkbutton(ra, text="1050nm (spheroid depth / faded-square)", variable=self._pl_postpa_oc_1050,
                       bg=BG2, fg=TEXT2, selectcolor=SURFACE, activebackground=BG2,
                       activeforeground=TEXT, font=("Segoe UI", 9)).pack(anchor="w")
        btn_apa = tk.Frame(apa, bg=BG2); btn_apa.pack(fill="x", padx=24, pady=(2, 6))
        self._btn(btn_apa, "Run After-PA Captures", lambda: self._pl_prepa_run_thread("postPA"),
                  "#f9e2af", "#1e1e2e", side="left")
        tk.Checkbutton(btn_apa, text="Locate only (1 plane)", variable=self._pl_postpa_locate_only,
                       bg=BG2, fg=TEXT2, selectcolor=SURFACE, activebackground=BG2,
                       activeforeground=TEXT, font=("Segoe UI", 8)).pack(side="left", padx=(8, 0))
        # "include in Run Pipeline" for this card lives in its header (top-left), matching
        # every other Step-4 card. In a pipeline these captures run ONLY if the PA job
        # actually produced output -- a post-PA set with no PA in between is a mislabelled
        # copy of the baseline.
        # The "Restore pre-PA light path" checkbox was removed 2026-08-08 -- the restore is
        # now unconditional (_pl_lightpath_restore stays True and is never unset). There is
        # no legitimate reason to capture the post-PA set through the PA optics: pa_setup
        # leaves dichroic OUT and the PA zoom, which on 2026-08-07 meant a 79.5 um field
        # against a 318.2 um baseline. That mismatch is invisible in the images -- it reads
        # as photoactivation -- so it must not be one tick away.
        tk.Label(apa, text="Pre-PA light path (dichroic IN + zoom + laser powers) is "
                           "restored automatically before these captures.",
                 bg=BG2, fg=GREEN, font=("Segoe UI", 8), justify="left",
                 wraplength=460).pack(anchor="w", padx=24, pady=(0, 6))

        # Bottom bar: run the checked cards in order.
        bottom = tk.Frame(pa, bg=BG); bottom.pack(fill="x", padx=8, pady=(4, 4))
        self._btn(bottom, "Run Pipeline (checked, in order)",
                  self._pl_pa_run_pipeline_thread, MAUVE, "#1e1e2e")
        self._pl_pa_status = tk.Label(pa, text="Photoactivation: idle", bg=BG, fg=SUBTEXT,
                                      font=("Segoe UI", 9), anchor="w", justify="left", wraplength=520)
        self._pl_pa_status.pack(fill="x", padx=8, pady=(0, 6))

        # ── Live dashboard: always visible in the middle pane below the steps ──
        tk.Label(dash_host, text="Live Dashboard", bg=BG, fg=MAUVE,
                 font=("Segoe UI", 9, "bold")).pack(anchor="w", padx=12, pady=(8, 2))
        if not _MPL_TK_OK:
            tk.Label(dash_host, text="matplotlib not installed — dashboard unavailable.",
                     bg=BG, fg=SUBTEXT, font=("Segoe UI", 9)).pack(anchor="w", padx=12)
        else:
            self._pl_dash_frame = tk.Frame(dash_host, bg=BG2, relief="flat")
            self._pl_dash_frame.pack(fill="both", expand=True, padx=12, pady=(0, 10))
            tk.Label(self._pl_dash_frame,
                     text="Dashboard will appear after Step 1 completes.",
                     bg=BG2, fg=SUBTEXT, font=("Segoe UI", 9)).pack(pady=4)

        # Select Step 1 by default
        self._pl_show_step("s1")
        # Set the initial dividers, then scan for any already-captured ND2s.
        self.after(150, self._pl_init_sash)
        self.after(220, self._pl_zv_refresh)
        # Resume the previous session's Steps 1-3 if this run folder has a snapshot.
        # Deferred so the log/table widgets exist; quiet=True so a first-ever launch
        # (no state file anywhere) doesn't open with a scary "nothing found" line.
        self.after(300, lambda: self._pl_restore_pipeline_state(quiet=True))
        # Confirm on close, and hold the window open until any in-flight macro reports
        # back -- closing mid-macro abandons a running capture (or a firing PA laser)
        # with nothing polling for its result. See _pl_on_close.
        self.protocol("WM_DELETE_WINDOW", self._pl_on_close)
        # Keep the PA folder tracking the Output directory, and create it. Deferred so the
        # log widget exists; the trace then follows every later change of run folder.
        self.after(320, self._pl_pa_sync_save_dir)
        self._pl_out_dir.trace_add("write", self._pl_pa_sync_save_dir)

    def _pl_init_sash(self):
        """Place the paned-window divider so the merged table pane gets the bulk of
        the width (all ten columns visible at full width) while the step/dashboard
        pane keeps ~500 px. Retries until the paned window has been laid out."""
        try:
            total = self._pl_paned.winfo_width()
        except Exception:
            return
        if total <= 100:
            self.after(80, self._pl_init_sash)
            return
        # Two panes: middle (steps+dashboard) | table; the table gets the rest.
        middle_w = 500   # enough for Step 4's rows
        if total - middle_w < 560:
            middle_w = max(460, total - 560)
        try:
            self._pl_paned.sashpos(0, middle_w)
        except Exception:
            pass

    # ── Step 3 captured-spheroid Z-stack viewer ───────────────────────────────

    def _pl_build_zviewer(self, parent):
        self._pl_zv_files  = {}
        self._pl_zv_stack  = None
        self._pl_zv_step   = None
        self._pl_zv_b2t    = True
        self._pl_zv_mid    = 0.0  # geometric middle plane index
        self._pl_zv_basez  = 0.0  # GUI Middle-plane Z (um) at the middle plane
        self._pl_zv_zabs   = None # true per-plane Z from ND2 events (or None)
        self._pl_zv_center = 0    # focus/centre plane index
        self._pl_zv_thumbs = []   # keep PhotoImage refs alive

        tk.Label(parent, text="Captured Spheroid Z-Stacks", bg=BG, fg=MAUVE,
                 font=("Segoe UI", 9, "bold")).pack(anchor="w", padx=12, pady=(8, 2))

        if not _PIL_OK:
            tk.Label(parent, text="Pillow not installed - image viewer unavailable "
                                  "(pip install pillow).",
                     bg=BG, fg=SUBTEXT, font=("Segoe UI", 9)).pack(anchor="w", padx=12)
            return

        ctl = tk.Frame(parent, bg=BG); ctl.pack(fill="x", padx=12, pady=(0, 2))
        tk.Label(ctl, text="Spheroid:", bg=BG, fg=TEXT2,
                 font=("Segoe UI", 9)).pack(side="left", padx=(0, 4))
        self._pl_zv_sel = tk.StringVar()
        self._pl_zv_combo = ttk.Combobox(ctl, textvariable=self._pl_zv_sel, width=18,
                                          state="readonly", font=("Segoe UI", 9))
        self._pl_zv_combo.pack(side="left", padx=(0, 6))
        self._pl_zv_combo.bind("<<ComboboxSelected>>", self._pl_zv_on_select)
        self._btn(ctl, "Auto Load", self._pl_zv_autoload, MAUVE, "#1e1e2e", side="left")
        self._btn(ctl, "Add...", self._pl_zv_add_file, BLUE, "#1e1e2e", side="left")
        self._btn(ctl, "Refresh", self._pl_zv_refresh, SURFACE2, TEXT, side="left")

        # Summary line: spheroid, middle-plane Z, plane count, geometry check.
        self._pl_zv_info = tk.Label(parent, text="No Z-stack loaded.",
                                    bg=BG, fg=TEXT, font=("Segoe UI", 9), anchor="w",
                                    justify="left", wraplength=520)
        self._pl_zv_info.pack(fill="x", padx=12, pady=(2, 2))

        tk.Label(parent, text="All focus planes — caption shows each plane's Z depth (from the "
                              "ND2 events Z Coord when recorded, else GUI Middle-plane Z + step); "
                              "focus-centre plane highlighted. Scroll horizontally.",
                 bg=BG, fg=SUBTEXT, font=("Segoe UI", 8)).pack(anchor="w", padx=12)

        # Horizontal filmstrip: one thumbnail per Z plane, captioned with its Z
        # depth. The middle plane (z-centre) is highlighted for orientation.
        strip_wrap = tk.Frame(parent, bg=BG2)
        strip_wrap.pack(fill="both", expand=True, padx=12, pady=(2, 10))
        self._pl_zv_strip_canvas = tk.Canvas(strip_wrap, bg=BG2, highlightthickness=0)
        hsb = ttk.Scrollbar(strip_wrap, orient="horizontal",
                            command=self._pl_zv_strip_canvas.xview)
        vsb = ttk.Scrollbar(strip_wrap, orient="vertical",
                            command=self._pl_zv_strip_canvas.yview)
        self._pl_zv_strip_canvas.configure(xscrollcommand=hsb.set, yscrollcommand=vsb.set)
        hsb.pack(side="bottom", fill="x")
        vsb.pack(side="right", fill="y")
        self._pl_zv_strip_canvas.pack(side="top", fill="both", expand=True)
        self._pl_zv_strip = tk.Frame(self._pl_zv_strip_canvas, bg=BG2)
        self._pl_zv_strip_canvas.create_window((0, 0), window=self._pl_zv_strip, anchor="nw")
        self._pl_zv_strip.bind(
            "<Configure>",
            lambda e: self._pl_zv_strip_canvas.configure(
                scrollregion=self._pl_zv_strip_canvas.bbox("all")))
        self._pl_zv_strip_canvas.bind(
            "<Enter>",
            lambda e: self._pl_zv_strip_canvas.bind_all("<MouseWheel>", self._pl_zv_wheel))
        self._pl_zv_strip_canvas.bind(
            "<Leave>", lambda e: self._pl_zv_strip_canvas.unbind_all("<MouseWheel>"))
        tk.Label(self._pl_zv_strip, text="Select or Add a captured Z-stack ND2.",
                 bg=BG2, fg="#888888", font=("Segoe UI", 9)).pack(padx=20, pady=20)

    def _pl_zv_refresh(self):
        if not _PIL_OK or not hasattr(self, "_pl_zv_combo"):
            return
        out = self._pl_out_dir.get().strip()
        cand = []
        if out:
            cand += [Path(out) / "nd2", Path(out) / "autofocus", Path(out)]
        nd2_out = self._pl_nd2_out_dir.get().strip()
        if nd2_out:
            cand.append(Path(nd2_out))
        files: dict = {}
        for d in cand:
            try:
                if d.is_dir():
                    for p in sorted(d.glob("*.nd2")):
                        files.setdefault(p.name, str(p))
            except Exception:
                pass
        self._pl_zv_files = files
        names = list(files.keys())
        self._pl_zv_combo.configure(values=names)
        cur = self._pl_zv_sel.get()
        if names and cur not in names:
            self._pl_zv_sel.set(names[0])
            self._pl_zv_load(files[names[0]])
        elif not names:
            self._pl_zv_info.configure(text="No captured ND2 found in nd2/ or autofocus/.")

    def _pl_zv_autoload(self):
        """Auto-discover the latest capture output dir (from C:/SpheroidPA/session.ini,
        which the capture daemon follows) and load all captured ND2 stacks into the
        viewer, newest first -- no manual path needed."""
        if not _PIL_OK or not hasattr(self, "_pl_zv_combo"):
            return
        import configparser
        import spheroid_pipeline as _pl
        cand = []
        try:
            cp = configparser.ConfigParser()
            cp.read(_pl.SESSION_INI)
            nd = cp.get("paths", "nd2_dir", fallback="").strip()
            if nd:
                cand.append(Path(nd))
            wd = cp.get("paths", "work_dir", fallback="").strip()
            if wd:
                cand.append(Path(wd).parent / "nd2")
        except Exception:
            pass
        nd2_out = self._pl_nd2_out_dir.get().strip()
        if nd2_out:
            cand.append(Path(nd2_out))
        out = self._pl_out_dir.get().strip()
        if out:
            cand.append(Path(out) / "nd2")
        files, chosen = {}, None
        for d in cand:
            try:
                if d.is_dir():
                    nd2s = sorted(d.glob("*.nd2"),
                                  key=lambda p: p.stat().st_mtime, reverse=True)
                    if nd2s:
                        chosen = d
                        for p in nd2s:
                            files.setdefault(p.name, str(p))
                        break
            except Exception:
                pass
        if not files:
            self._pl_zv_info.configure(
                text="Auto Load: no captured ND2 found (run a capture first).")
            return
        self._pl_zv_files = files
        names = list(files.keys())
        self._pl_zv_combo.configure(values=names)
        self._pl_zv_sel.set(names[0])
        if chosen is not None:
            self._pl_nd2_out_dir.set(str(chosen))
        self._pl_log(f"Auto Load: {len(names)} ND2 from {chosen} (newest first)")
        self._pl_zv_info.configure(text=f"Auto Load: reading {len(names)} stack(s) ...")
        items = [(n, files[n]) for n in names]
        threading.Thread(target=self._pl_zv_load_all, args=(items,), daemon=True).start()

    def _pl_zv_load_captured(self, items, plabel="Validation"):
        """Auto-load the ND2s a Validation run JUST captured into the "Captured
        Z-Stacks" tab -- called from _pl_prepa_capture_ocs on success instead of
        requiring a manual Auto Load/Refresh click. `items` is [(display_name,
        path), ...], one entry per (checked OC, spheroid).

        MERGES into the existing list (never rebinds): the After-PA card would
        otherwise wipe every prePA entry, since the display names are phase-
        qualified but the dict was replaced wholesale -- defeating the whole point
        of the before/after split. Also skips paths the capture didn't actually
        produce, so a partial run can't leave dead entries in the dropdown."""
        if not _PIL_OK or not hasattr(self, "_pl_zv_combo") or not items:
            return
        live = [(n, p) for n, p in items if Path(p).exists()]
        missing = len(items) - len(live)
        if missing:
            self._pl_log(f"{plabel}: {missing} captured stack(s) not found on disk -- skipped")
        if not live:
            return
        if not isinstance(getattr(self, "_pl_zv_files", None), dict):
            self._pl_zv_files = {}
        had_rows = bool(self._pl_zv_files) and getattr(self, "_pl_zv_list_mode", False)
        self._pl_zv_files.update({n: p for n, p in live})
        names = list(self._pl_zv_files.keys())
        self._pl_zv_combo.configure(values=names)
        self._pl_zv_sel.set(live[0][0])
        self._pl_log(f"{plabel}: auto-loading {len(live)} captured stack(s) into Captured Z-Stacks "
                     f"({len(names)} total listed)")
        self._pl_zv_info.configure(text=f"{plabel}: reading {len(live)} stack(s) ...")
        # append when rows are already rendered (e.g. After-PA following Pre-PA), so the
        # before/after rows sit side by side instead of the new run clearing the strip.
        threading.Thread(target=self._pl_zv_load_all, args=(live,),
                         kwargs={"append": had_rows}, daemon=True).start()

    def _pl_zv_add_file(self):
        if not _PIL_OK or not hasattr(self, "_pl_zv_combo"):
            return
        p = filedialog.askopenfilename(
            title="Select an ND2 file to view",
            filetypes=[("Nikon ND2 files", "*.nd2"), ("All files", "*.*")])
        if not p:
            return
        name = Path(p).name
        key = name
        if key in self._pl_zv_files and self._pl_zv_files[key] != p:
            key = f"{name}  [{Path(p).parent.name}]"
        self._pl_zv_files[key] = p
        self._pl_zv_combo.configure(values=list(self._pl_zv_files.keys()))
        self._pl_zv_sel.set(key)
        self._pl_zv_load(p)

    def _pl_zv_on_select(self, event=None):
        p = self._pl_zv_files.get(self._pl_zv_sel.get())
        if p:
            self._pl_zv_load(p)

    def _pl_zv_load(self, path):
        self._pl_zv_info.configure(text=f"Loading {Path(path).name} ...")
        threading.Thread(target=self._pl_zv_load_worker, args=(path,), daemon=True).start()

    def _pl_zv_load_worker(self, path):
        try:
            import nd2 as _nd2
            import numpy as _np
            step = None; home = 0; b2t = True
            with _nd2.ND2File(path) as f:
                sizes = dict(f.sizes)
                chan_name = None
                try:
                    chan_name = _safe_attr(f.metadata.channels[0], "channel", "name", default=None)
                except Exception:
                    chan_name = None
                arr = _np.asarray(f.asarray())
                try:
                    for lp in f.experiment:
                        if type(lp).__name__ == "ZStackLoop":
                            pr = lp.parameters
                            step = float(getattr(pr, "stepUm"))
                            home = int(getattr(pr, "homeIndex", 0))
                            b2t  = bool(getattr(pr, "bottomToTop", True))
                            break
                except Exception:
                    pass
                # True per-plane absolute Z from the acquisition events log
                # (the scanned Ti ZDrive / Z Coord per frame), plus the home/focus
                # plane where the Z-series offset is 0. This is the real recorded
                # Z; stagePositionUm.z is only the constant coarse-stage snapshot.
                zabs = None; ev_home = None
                try:
                    rows = f.events()
                    if rows:
                        keys = list(rows[0].keys())
                        def _fk(*subs):
                            for k in keys:
                                kl = k.lower()
                                if all(s in kl for s in subs):
                                    return k
                            return None
                        zi_key = _fk("z", "index")
                        zc_key = _fk("z", "coord") or _fk("zdrive") or _fk("ti", "z")
                        zs_key = _fk("series")
                        nz = int(sizes.get("Z", 1))
                        if zi_key and zc_key:
                            zmap = {}; smap = {}
                            for r in rows:
                                iv = r.get(zi_key)
                                if iv is None:
                                    continue
                                iv = int(iv)
                                cv = r.get(zc_key)
                                if cv is not None:
                                    zmap[iv] = float(cv)
                                sv = r.get(zs_key) if zs_key else None
                                if sv is not None:
                                    smap[iv] = float(sv)
                            if len(zmap) >= nz and all(j in zmap for j in range(nz)):
                                zabs = [zmap[j] for j in range(nz)]
                            if smap:
                                ev_home = min(smap, key=lambda k: abs(smap[k]))
                except Exception:
                    zabs = None; ev_home = None
            order = list(sizes.keys())
            sel = [slice(None) if ax in ("Z", "Y", "X") else 0 for ax in order]
            arr = arr[tuple(sel)]
            kept = [ax for ax in order if ax in ("Z", "Y", "X")]
            if "Z" in kept:
                arr = _np.moveaxis(arr, kept.index("Z"), 0)
            else:
                arr = arr[None, ...]
            arr = _np.ascontiguousarray(arr)
            vmin = float(_np.percentile(arr, 1)); vmax = float(_np.percentile(arr, 99.5))
            if vmax <= vmin:
                vmax = vmin + 1.0
        except Exception as exc:
            self.after(0, lambda e=exc: self._pl_zv_info.configure(text=f"Load error: {e}"))
            return
        n_planes = arr.shape[0]
        # GUI Middle-plane Z / half-range / step (fallback anchor + geometry check).
        gui_centre = gui_half = gui_step = None
        try: gui_centre = float(self._pl_z_centre.get())
        except Exception: pass
        try: gui_half = float(self._pl_z_half.get())
        except Exception: pass
        try: gui_step = float(self._pl_z_step.get())
        except Exception: pass
        if step is None:
            step = gui_step if (gui_step and gui_step > 0) else 1.0
        # Prefer the TRUE per-plane Z recorded by NIS-E in the events log
        # (Z Coord / Ti ZDrive). Only trust it when it actually varies across the
        # stack; otherwise fall back to GUI Middle-plane Z + step reconstruction.
        if zabs is not None and (max(zabs) - min(zabs)) <= 1e-6:
            zabs = None
        self._pl_zv_stack  = arr
        self._pl_zv_chan   = chan_name
        self._pl_zv_vmin, self._pl_zv_vmax = vmin, vmax
        self._pl_zv_step   = step
        self._pl_zv_b2t    = b2t
        self._pl_zv_mid    = (n_planes - 1) / 2.0          # geometric middle index
        self._pl_zv_basez  = gui_centre if gui_centre is not None else 0.0
        self._pl_zv_zabs   = zabs                          # true per-plane Z, or None
        if zabs is not None and ev_home is not None:
            self._pl_zv_center = int(ev_home)              # events home (Z-Series = 0)
        else:
            self._pl_zv_center = int(round(self._pl_zv_mid))
        self._pl_zv_check  = self._pl_zv_geom_check(n_planes, home, gui_half, gui_step)
        self._pl_zv_list_mode = False
        self.after(0, self._pl_zv_render_filmstrip)

    def _pl_zv_geom_check(self, n, home_idx, gui_half, gui_step):
        """Double-check the stack against the GUI Z geometry: does the plane count
        match 2*half/step+1 (so the middle really is z_centre)?"""
        msgs = []
        if gui_half and gui_step and gui_step > 0:
            exp = int(round(2.0 * gui_half / gui_step)) + 1
            if exp == n:
                msgs.append(f"{n} planes match +/-{gui_half:g}/{gui_step:g} um")
            else:
                msgs.append(f"{n} planes vs expected {exp} for +/-{gui_half:g}/{gui_step:g} um")
        return "  |  ".join(msgs)

    def _pl_zv_render_filmstrip(self):
        import numpy as _np
        for w in self._pl_zv_strip.winfo_children():
            w.destroy()
        self._pl_zv_thumbs = []
        arr = self._pl_zv_stack
        if arr is None:
            return
        z = arr.shape[0]
        vmin, vmax = self._pl_zv_vmin, self._pl_zv_vmax
        zabs   = self._pl_zv_zabs
        center = self._pl_zv_center
        mid    = self._pl_zv_mid
        sign   = 1.0 if self._pl_zv_b2t else -1.0
        chan   = getattr(self, "_pl_zv_chan", None)
        # Find Z-Center override: if this spheroid has a found result, highlight its focal
        # plane (nearest by absolute Z when zabs is populated, else its plane_index) in GREEN
        # when confident / AMBER on the geometric-middle fallback, keeping the blue geometric
        # middle as a faint secondary marker.
        geom_center = center
        found = None; found_idx = None
        if isinstance(getattr(self, "_pl_zc_found", None), dict) and self._pl_zc_found:
            found = self._pl_zc_found.get(self._pl_zc_sid_from_key(self._pl_zv_sel.get()))
        if found:
            if zabs is not None:
                found_idx = int(_np.argmin(_np.abs(
                    _np.asarray(zabs, dtype=float) - float(found["z_centre_um"]))))
            else:
                found_idx = int(found.get("plane_index", geom_center))
            found_idx = max(0, min(found_idx, z - 1))
            self._pl_zv_center = found_idx
            center = found_idx
        THUMB = 150
        for i in range(z):
            a = _np.clip((arr[i].astype(_np.float32) - vmin) / (vmax - vmin), 0.0, 1.0)
            a = (a * 255.0).astype("uint8")
            im = _PILImage.fromarray(_tint_channel(a, chan), mode="RGB")
            im.thumbnail((THUMB, THUMB))
            photo = _PILImageTk.PhotoImage(im)
            self._pl_zv_thumbs.append(photo)
            if zabs is not None:
                absz = zabs[i]
            else:
                absz = self._pl_zv_basez + (i - mid) * (self._pl_zv_step or 0.0) * sign
            is_found = (found is not None and i == found_idx)
            is_geom  = (i == geom_center)
            if is_found:
                conf = bool(found.get("confidence"))
                edge = GREEN if conf else PEACH
                cellbg = SURFACE2
                note = (f"found center  Zc={float(found['z_centre_um']):.1f} um" if conf
                        else f"middle (fallback: {found.get('reason')})")
                note_fg = edge
            elif is_geom and found is not None:
                edge = LAVENDER; cellbg = BG2; note = "middle"; note_fg = SUBTEXT
            elif is_geom:
                edge = BLUE; cellbg = SURFACE2; note = "focus centre"; note_fg = BLUE
            else:
                edge = SURFACE; cellbg = BG2; note = ""; note_fg = TEXT2
            cell = tk.Frame(self._pl_zv_strip, bg=cellbg, padx=3, pady=3,
                            highlightthickness=2, highlightbackground=edge,
                            highlightcolor=edge)
            cell.pack(side="left", padx=3, pady=4)
            tk.Label(cell, image=photo, bg=cellbg).pack()
            cap = f"Z {absz:.1f} um\n" + note
            tk.Label(cell, text=cap, bg=cellbg, fg=note_fg,
                     font=("Segoe UI", 8), justify="center").pack()
        self._pl_zv_strip_canvas.update_idletasks()
        self._pl_zv_strip_canvas.configure(
            scrollregion=self._pl_zv_strip_canvas.bbox("all"))
        chk = getattr(self, "_pl_zv_check", "")
        if zabs is not None:
            centre_z = zabs[center]; src = "from ND2 events (Z Coord)"
        else:
            centre_z = self._pl_zv_basez; src = "from GUI Middle-plane Z"
        self._pl_zv_info.configure(text=(
            f"Spheroid: {self._pl_zv_sel.get()}     "
            f"centre plane Z {centre_z:.1f} um ({src})     {z} planes, step {self._pl_zv_step:g} um     "
            f"{arr.shape[2]}x{arr.shape[1]} px"
            + (f"\n[check] {chk}" if chk else "")))

    def _pl_zv_wheel(self, event):
        try:
            if getattr(self, "_pl_zv_list_mode", False) and not (event.state & 0x0001):
                self._pl_zv_strip_canvas.yview_scroll(int(-event.delta / 120), "units")
            else:
                self._pl_zv_strip_canvas.xview_scroll(int(-event.delta / 120), "units")
        except Exception:
            pass

    def _pl_zv_begin_list(self):
        for w in self._pl_zv_strip.winfo_children():
            w.destroy()
        self._pl_zv_thumbs = []
        self._pl_zv_list_mode = True
        self._pl_zv_strip_canvas.xview_moveto(0.0)
        self._pl_zv_strip_canvas.yview_moveto(0.0)

    def _pl_zv_load_all(self, items, append=False):
        """List mode: render every (name, path) as its own filmstrip row, newest first.
        Lighter than the single-stack view (smaller thumbs, no per-plane caption).

        append=True keeps the rows already on screen and adds these below them --
        used by the After-PA card so its rows join the pre-PA ones instead of
        replacing them (_pl_zv_begin_list destroys every child of the strip, which
        would otherwise wipe the before/after comparison at the render layer even
        though _pl_zv_files was correctly merged)."""
        import nd2 as _nd2
        import numpy as _np
        if not append:
            self.after(0, self._pl_zv_begin_list)
        gui_centre = gui_step = None
        try: gui_centre = float(self._pl_z_centre.get())
        except Exception: pass
        try: gui_step = float(self._pl_z_step.get())
        except Exception: pass

        # Metadata-only first pass (voxel_size + frame size, no array read): different
        # OCs can carry different native zoom (1050nm ~318 um FOV vs 890/940nm ~636 um --
        # see MD/CLAUDE.md "940/1050nm scale difference" note), so wider-FOV channels
        # get center-cropped down to the narrowest one before display, for a like-for-
        # like comparison across rows of the same spheroid.
        pix_um = {}
        chan_of = {}
        for name, path in items:
            try:
                with _nd2.ND2File(path) as f:
                    vox = f.voxel_size()
                    sz = dict(f.sizes)
                    pix_um[name] = (float(vox.x), sz.get("X", 0))
                    try:
                        chan_of[name] = _safe_attr(f.metadata.channels[0], "channel", "name", default=None)
                    except Exception:
                        chan_of[name] = None
            except Exception:
                pix_um[name] = (None, 0)
        fovs = [p * n for p, n in pix_um.values() if p and n]
        target_fov_um = min(fovs) if fovs else None

        done = 0
        for name, path in items:
            try:
                with _nd2.ND2File(path) as f:
                    sizes = dict(f.sizes)
                    arr = _np.asarray(f.asarray())
                    b2t = True
                    try:
                        for lp in f.experiment:
                            if type(lp).__name__ == "ZStackLoop":
                                b2t = bool(getattr(lp.parameters, "bottomToTop", True))
                                break
                    except Exception:
                        pass
                order = list(sizes.keys())
                sel = [slice(None) if ax in ("Z", "Y", "X") else 0 for ax in order]
                arr = arr[tuple(sel)]
                kept = [ax for ax in order if ax in ("Z", "Y", "X")]
                arr = _np.moveaxis(arr, kept.index("Z"), 0) if "Z" in kept else arr[None, ...]
                arr = _np.ascontiguousarray(arr)

                p_um, _ = pix_um.get(name, (None, 0))
                cropped = False
                if p_um and target_fov_um:
                    arr2 = _center_crop_to_fov(arr, p_um, target_fov_um)
                    if arr2.shape != arr.shape:
                        arr = arr2
                        cropped = True

                vmin = float(_np.percentile(arr, 1)); vmax = float(_np.percentile(arr, 99.5))
                if vmax <= vmin:
                    vmax = vmin + 1.0
                n = arr.shape[0]
                st = dict(arr=arr, vmin=vmin, vmax=vmax, n_planes=n,
                          center=int(round((n - 1) / 2.0)),
                          basez=(gui_centre if gui_centre is not None else 0.0),
                          step=(gui_step or 0.0), b2t=b2t, cropped=cropped,
                          chan=chan_of.get(name))
                self.after(0, lambda nm=name, s=st: self._pl_zv_append_row(nm, s))
                done += 1
            except Exception:
                continue
        self.after(0, lambda d=done: self._pl_zv_info.configure(
            text=(f"Auto Load: {d} stack(s) listed, newest first — one row per spheroid "
                  "(focus-centre plane outlined). Rows marked (FOV-matched crop) are "
                  "center-cropped to the narrowest channel's physical field of view for "
                  "like-for-like comparison. Wheel scrolls the list; pick one in the "
                  "dropdown for the full-size filmstrip.")))

    def _pl_zv_append_row(self, name, st):
        import numpy as _np
        arr = st["arr"]; vmin = st["vmin"]; vmax = st["vmax"]; center = st["center"]
        # Find Z-Center override for the list row (no per-plane zabs here, so map via
        # plane_index); keep the geometric middle as a faint secondary marker.
        geom_center = center
        found = None; found_idx = None
        if isinstance(getattr(self, "_pl_zc_found", None), dict) and self._pl_zc_found:
            found = self._pl_zc_found.get(self._pl_zc_sid_from_key(name))
        if found:
            found_idx = max(0, min(int(found.get("plane_index", geom_center)), arr.shape[0] - 1))
            center = found_idx
        row = tk.Frame(self._pl_zv_strip, bg=BG2)
        row.pack(side="top", fill="x", anchor="w", pady=(2, 8))
        tag = "  (FOV-matched crop)" if st.get("cropped") else ""
        note = ""
        if found is not None:
            note = ("   [found Zc={:.1f} um]".format(float(found["z_centre_um"]))
                    if found.get("confidence")
                    else "   [middle fallback: {}]".format(found.get("reason")))
        tk.Label(row,
                 text=f"{name}     centre Z {st['basez']:.1f} um     {st['n_planes']} planes{tag}{note}",
                 bg=BG2, fg=MAUVE, font=("Segoe UI", 9, "bold"), anchor="w"
                 ).pack(side="top", anchor="w", padx=2)
        sub = tk.Frame(row, bg=BG2); sub.pack(side="top", anchor="w")
        chan = st.get("chan")
        for i in range(arr.shape[0]):
            a = _np.clip((arr[i].astype(_np.float32) - vmin) / (vmax - vmin), 0.0, 1.0)
            a = (a * 255.0).astype("uint8")
            im = _PILImage.fromarray(_tint_channel(a, chan), mode="RGB")
            im.thumbnail((96, 96))
            photo = _PILImageTk.PhotoImage(im)
            self._pl_zv_thumbs.append(photo)
            is_found = (found is not None and i == found_idx)
            is_geom  = (i == geom_center)
            if is_found:
                edge = GREEN if found.get("confidence") else PEACH
                cbg = SURFACE2
            elif is_geom and found is not None:
                edge = LAVENDER; cbg = BG2
            elif is_geom:
                edge = BLUE; cbg = SURFACE2
            else:
                edge = SURFACE; cbg = BG2
            cell = tk.Frame(sub, bg=cbg, padx=1, pady=1, highlightthickness=2,
                            highlightbackground=edge, highlightcolor=edge)
            cell.pack(side="left", padx=1, pady=1)
            tk.Label(cell, image=photo, bg=cbg).pack()
        self._pl_zv_strip_canvas.update_idletasks()
        self._pl_zv_strip_canvas.configure(
            scrollregion=self._pl_zv_strip_canvas.bbox("all"))

    # ── Pipeline tab helpers ──────────────────────────────────────────────────

    def _pl_browse_dir(self, var: tk.StringVar):
        d = filedialog.askdirectory(title="Select directory")
        if d:
            var.set(d)

    def _pl_browse_bin(self, var: tk.StringVar):
        p = filedialog.askopenfilename(
            title="Select rig reference .bin",
            filetypes=[("NIS-E Z-correction bin", "*.bin"), ("All files", "*.*")])
        if p:
            var.set(p)

    def _pl_show_step(self, key: str):
        for k, frame in self._pl_step_frames.items():
            frame.pack_forget()
        self._pl_step_frames[key].pack(side="top", fill="x")
        # Right pane: Step 4 shows the Photoactivation macro cards; steps 1-3 show
        # the merged Spheroid State & Dashboard table.
        if key == "s4":
            self._pl_table_host.pack_forget()
            self._pl_pa_host.pack(fill="both", expand=True)
        else:
            self._pl_pa_host.pack_forget()
            self._pl_table_host.pack(fill="both", expand=True)
        for k, card in self._pl_step_cards.items():
            bg = SURFACE2 if k == key else BG2
            card.configure(bg=bg)
            for w in card.winfo_children():
                w.configure(bg=bg)

    def _pl_update_step_dots(self):
        s1_fg = GREEN if self._pl_records else SUBTEXT
        s2_fg = GREEN if any(r.status not in ("DETECTED",) for r in self._pl_records) else SUBTEXT
        s3_fg = GREEN if any(r.status in ("Z_KNOWN", "BIN_READY", "QUEUED", "IMAGING", "IMAGED")
                             for r in self._pl_records) else SUBTEXT
        s4_fg = GREEN if any(r.status in ("BIN_READY", "QUEUED", "IMAGING", "IMAGED")
                             for r in self._pl_records) else SUBTEXT
        for k, fg in [("s1", s1_fg), ("s2", s2_fg), ("s3", s3_fg), ("s4", s4_fg)]:
            if k in self._pl_step_dots:
                self.after(0, lambda lbl=self._pl_step_dots[k], c=fg: lbl.configure(fg=c))

    def _pl_refresh_dashboard(self):
        if not _MPL_TK_OK or not self._pl_records:
            return
        try:
            import spheroid_pipeline as _pl
        except ImportError:
            return

        mosaic_path = self._pl_mosaic_path.get().strip()
        out_dir     = self._pl_out_dir.get().strip()

        if self._pl_mpl_canvas is None or self._pl_dash_mosaic != mosaic_path:
            # Clear placeholder label
            for w in self._pl_dash_frame.winfo_children():
                w.destroy()
            if self._pl_dashboard is not None:
                try:
                    self._pl_dashboard.close()
                except Exception:
                    pass

            fig = _MplFigure(figsize=(13, 8), facecolor="#1e1e2e")
            canvas = _FigCanvas(fig, master=self._pl_dash_frame)
            canvas.get_tk_widget().pack(fill="both", expand=True)
            self._pl_mpl_canvas = canvas
            self._pl_dashboard  = _pl.PipelineDashboard(
                mosaic_nd2 = Path(mosaic_path) if mosaic_path else None,
                out_dir    = Path(out_dir) if out_dir else None,
                fig        = fig,
            )
            self._pl_dash_mosaic = mosaic_path

        self._pl_dashboard.update(self._pl_records)
        self._pl_mpl_canvas.draw()

    def _pl_update_table(self):
        for row in self._pl_tree.get_children():
            self._pl_tree.delete(row)
        rc   = getattr(self, "_pl_recentered", {})
        excl = getattr(self, "_pl_excluded", set())
        for r in self._pl_records:
            mos = (f"({r.mosaic_x_um:.0f}, {r.mosaic_y_um:.0f})"
                   if r.mosaic_x_um else "—")
            xy = (f"({r.verified_x_um:.1f}, {r.verified_y_um:.1f})"
                  if r.verified_x_um else "—")
            tags = (str(r.status),)
            if r.rank in rc:
                dx, dy = rc[r.rank]
                xy = f"{xy}  Δ{dx:+.0f},{dy:+.0f}"
                tags = ("recentered",)
            zc = f"{r.z_centre_um:.1f}" if r.z_centre_um else "—"
            dm = f"{r.mosaic_diam_um:.1f}" if r.mosaic_diam_um else "—"
            sc = f"{r.mosaic_score:.4f}" if r.mosaic_score else "—"
            bn = Path(r.bin_path).name if r.bin_path else "—"
            use = "[ ]" if r.rank in excl else "[x]"
            self._pl_tree.insert("", "end", iid=f"r{r.rank}", tags=tags, values=(
                use, r.rank, r.spheroid_id, r.status, mos, xy, zc, dm, sc, bn))

    def _pl_toggle_use(self, event):
        """Toggle a spheroid's Use checkbox when its Use cell ([x]/[ ]) is clicked."""
        if self._pl_tree.identify_region(event.x, event.y) != "cell":
            return
        if self._pl_tree.identify_column(event.x) != "#1":   # the "use" column
            return
        row = self._pl_tree.identify_row(event.y)
        if not row or not row.startswith("r"):
            return
        try:
            rank = int(row[1:])
        except ValueError:
            return
        excl = self._pl_excluded
        if rank in excl:
            excl.discard(rank)
        else:
            excl.add(rank)
        self._pl_tree.set(row, "use", "[ ]" if rank in excl else "[x]")

    def _pl_toggle_all_use(self):
        """Select-all / none for the Use column: include every spheroid, or -- if all
        are already included -- exclude every one. Same effect as clicking each Use box."""
        if not self._pl_records:
            return
        if getattr(self, "_pl_excluded", set()):    # some excluded -> include all
            self._pl_excluded = set()
        else:                                        # all included -> exclude all
            self._pl_excluded = {r.rank for r in self._pl_records}
        self._pl_update_table()

    def _pl_log(self, msg: str):
        self._log_line(f"[Pipeline] {msg}", "info")

    # ── Step 1 handlers ───────────────────────────────────────────────────────

    def _pl_run_screener_thread(self):
        threading.Thread(target=self._pl_run_screener, daemon=True).start()

    # ── Session persistence ───────────────────────────────────────────────────
    # Steps 1-3 are expensive (a 121-tile mosaic + screening + anchoring), and closing
    # the GUI -- deliberately or by a crash mid-troubleshooting -- used to throw all of
    # it away because _pl_records lived only in memory. save_state/load_state already
    # existed in spheroid_pipeline but nothing ever called them. Now every step that
    # changes the record set writes pipeline_state.json into the Output dir, and startup
    # restores it, so a restart resumes at Step 4 instead of re-imaging the plate.

    # ── Safe shutdown ─────────────────────────────────────────────────────────
    # Closing mid-macro is not a clean no-op: the NIS-E side keeps running, so the GUI
    # stops polling cmd_done, the af_trigger/cmd handshake is left half-consumed, and a
    # capture or PA can complete with nothing recording where it went. Worse, closing
    # during a step3_zstack_PA run abandons a firing 850 nm laser with no watcher.

    CLOSE_WAIT_TIMEOUT_S = 1800        # re-ask after 30 min rather than wait forever

    def _pl_inflight_macro(self):
        """Describe the macro currently in flight, or None if the rig is idle.

        Two independent signals, because either alone can miss one:
          - _pl_dispatch_busy: the GUI's own single-in-flight lock. Misses a macro whose
            cmd_done poller thread died (then the flag is stuck False or True wrongly).
          - cmd.ini present on disk with no cmd_done.ini beside it: the dispatcher has
            been handed work and has not reported finishing. Survives GUI-side thread
            death, and is what actually reflects NIS-E's state.

        The always-on nis_macro_dispatcher.mac is deliberately NOT counted -- it idles
        waiting for cmd.ini and never "finishes", so waiting on it would hang forever."""
        import configparser, spheroid_pipeline as _pl
        work = ""
        try:
            cp = configparser.ConfigParser()
            if _pl.SESSION_INI.exists():
                cp.read(_pl.SESSION_INI)
            work = cp.get("paths", "work_dir", fallback="").strip()
        except Exception:
            pass
        if not work:
            work = self._pl_trigger_dir.get().strip()
        action = None
        if work:
            wd = Path(work)
            cmd, done = wd / "cmd.ini", wd / "cmd_done.ini"
            try:
                if cmd.exists() and not done.exists():
                    cp = configparser.ConfigParser()
                    cp.read(cmd)
                    action = cp.get("command", "action", fallback="a macro")
            except Exception:
                action = "a macro"
        if action is None and getattr(self, "_pl_dispatch_busy", False):
            action = "a macro (GUI dispatcher lock held)"
        return action

    def _pl_on_close(self):
        """WM_DELETE_WINDOW handler: confirm, and hold the window open until any in-flight
        macro reports back. Never force-closes silently."""
        busy = self._pl_inflight_macro()
        if busy is None:
            if messagebox.askyesno("Close SpheroidPA",
                                   "Close SpheroidPA?\n\nNo macro is currently running."):
                self._pl_shutdown()
            return
        if not messagebox.askyesno(
                "Macro still running",
                f"'{busy}' is still running in NIS-E.\n\n"
                f"Closing now would stop the GUI polling for its result: the macro keeps "
                f"running on the rig, and nothing records where its output went.\n\n"
                f"Yes  — wait for it to finish, then close (recommended)\n"
                f"No   — stay open, do not close"):
            return
        self._pl_log(f"Shutdown: waiting for '{busy}' to finish before closing...")
        self._pl_wait_then_close(busy, 0)

    def _pl_wait_then_close(self, busy, waited):
        """Poll once a second until the macro reports back, then close. Re-asks rather than
        waiting unboundedly, so a macro that never reports cannot wedge the window shut."""
        still = self._pl_inflight_macro()
        if still is None:
            self._pl_log("Shutdown: macro finished — closing.")
            self._pl_shutdown()
            return
        if waited and waited % 15 == 0:
            try:
                self._pl_dispatch_status.configure(
                    text=f"Closing after '{still}' finishes… {waited}s", fg=YELLOW)
            except Exception:
                pass
        if waited >= self.CLOSE_WAIT_TIMEOUT_S:
            if messagebox.askyesno(
                    "Still running",
                    f"'{still}' has not reported back after "
                    f"{self.CLOSE_WAIT_TIMEOUT_S // 60} minutes.\n\n"
                    f"Yes — close anyway (the macro keeps running in NIS-E)\n"
                    f"No  — keep waiting"):
                self._pl_log(f"Shutdown: closing with '{still}' still unreported "
                             f"(operator override).")
                self._pl_shutdown()
                return
            waited = 0
        self.after(1000, lambda: self._pl_wait_then_close(busy, waited + 1))

    def _pl_shutdown(self):
        """Persist the session, then tear the window down. The state save is best-effort:
        losing it costs a re-run of Steps 1-3, whereas raising here would leave the window
        un-closable."""
        try:
            self._pl_save_pipeline_state()
        except Exception as e:
            self._pl_log(f"Shutdown: session save failed ({e}) — closing anyway.")
        try:
            self.destroy()
        except Exception:
            pass

    def _pl_save_pipeline_state(self):
        """Persist records + the Step 1-3 UI state to <out_dir>/pipeline_state.json.
        Best-effort: this must never break a pipeline step, so all failures are logged
        and swallowed (a lost snapshot costs a re-run; a raised exception here would
        abort the step that just succeeded)."""
        try:
            import spheroid_pipeline as _pl
            out = self._pl_out_dir.get().strip()
            if not out or not self._pl_records:
                return
            def _f(var, default=0.0):
                try:
                    return float(var.get())
                except (ValueError, AttributeError):
                    return default
            st = _pl.PipelineState(
                well_id     = self._pl_well_id.get().strip(),
                mosaic_nd2  = self._pl_mosaic_path.get().strip(),
                out_dir     = out,
                trigger_dir = self._pl_trigger_dir.get().strip(),
                nd2_out_dir = self._pl_nd2_out_dir.get().strip(),
                records     = self._pl_records,
                excluded_ranks = sorted(getattr(self, "_pl_excluded", set())),
                z_centre    = _f(self._pl_z_centre),
                z_half      = _f(self._pl_z_half),
                z_step      = _f(self._pl_z_step),
            )
            _pl.save_state(st)
            self._pl_state = st
        except Exception as e:
            self._pl_log(f"Session save failed (continuing): {e}")

    def _pl_restore_pipeline_state(self, out_dir=None, quiet=False):
        """Restore a previous session's Steps 1-3 from <out_dir>/pipeline_state.json.

        out_dir defaults to the GUI's Output dir, falling back to the parent of
        session.ini's work_dir -- which is the run folder the dispatcher is already
        pointed at, so a fresh GUI (whose Output dir is still the generic default)
        finds the session that was actually in progress. Returns True if restored."""
        try:
            import spheroid_pipeline as _pl, configparser
            cand = []
            if out_dir:
                cand.append(Path(out_dir))
            cur = self._pl_out_dir.get().strip()
            if cur:
                cand.append(Path(cur))
            try:
                cp = configparser.ConfigParser()
                if _pl.SESSION_INI.exists():
                    cp.read(_pl.SESSION_INI)
                wd = cp.get("paths", "work_dir", fallback="").strip()
                if wd:
                    cand.append(Path(wd).parent)
            except Exception:
                pass
            st = None
            for c in cand:
                try:
                    st = _pl.load_state(c)
                except Exception as e:
                    self._pl_log(f"Session restore: {c} unreadable ({e})")
                    continue
                if st and st.records:
                    break
                st = None
            if st is None:
                if not quiet:
                    self._pl_log("Session restore: no pipeline_state.json with records found.")
                return False

            self._pl_state   = st
            self._pl_records = st.records
            self._pl_excluded = set(st.excluded_ranks or [])
            for var, val in ((self._pl_well_id, st.well_id),
                             (self._pl_mosaic_path, st.mosaic_nd2),
                             (self._pl_out_dir, st.out_dir),
                             (self._pl_trigger_dir, st.trigger_dir),
                             (self._pl_nd2_out_dir, st.nd2_out_dir)):
                if val:
                    var.set(val)
            for var, val in ((self._pl_z_centre, st.z_centre),
                             (self._pl_z_half, st.z_half),
                             (self._pl_z_step, st.z_step)):
                if val:
                    var.set(str(val))
            ranks = [str(r.rank) for r in st.records]
            try:
                self._pl_z_rank_combo.configure(values=ranks)
                if ranks and not self._pl_z_rank.get().strip():
                    self._pl_z_rank.set(ranks[0])
            except Exception:
                pass
            self._pl_update_table()
            self._pl_refresh_dashboard()
            self._pl_update_step_dots()
            n_ex = len(self._pl_excluded)
            msg = (f"Session restored: {len(st.records)} spheroid(s) for well "
                   f"{st.well_id or '?'}"
                   + (f", {n_ex} excluded" if n_ex else "")
                   + f" (from {Path(st.out_dir).as_posix()}). Steps 1-3 do NOT need re-running.")
            self._pl_log(msg)
            try:
                self._pl_screen_lbl.configure(text=msg, fg=GREEN)
            except Exception:
                pass
            return True
        except Exception as e:
            self._pl_log(f"Session restore failed: {e}")
            return False

    def _pl_run_screener(self):
        try:
            import spheroid_pipeline as _pl
        except ImportError as e:
            messagebox.showerror("Import error", f"spheroid_pipeline.py not found:\n{e}")
            return
        nd2 = self._pl_mosaic_path.get().strip()
        if not nd2:
            messagebox.showwarning("No ND2", "Select a 10X mosaic ND2 first."); return
        out = self._pl_out_dir.get().strip()
        if not out:
            out = str(Path(nd2).parent / "screener_out")
            self.after(0, lambda: self._pl_out_dir.set(out))
        well = self._pl_well_id.get().strip() or "Well"
        self._pl_log(f"Screening {Path(nd2).name} ...")
        try:
            records = _pl.screen_mosaic(Path(nd2), Path(out), well)
            n_str = self._pl_n_spheroids.get().strip()
            if n_str:
                try:
                    records = records[:max(1, int(n_str))]
                except ValueError:
                    pass
            self._pl_records = records
            # Suggest the top-ranked spheroids as editable anchor defaults
            # (records are score-ordered; rank 1 is highest). The operator can
            # overwrite these with whichever spheroid they actually navigate to.
            for i in range(min(2, len(records))):
                if i < len(self._pl_anchor_ranks):
                    self.after(0, lambda v=records[i].rank, var=self._pl_anchor_ranks[i]:
                               var.set(str(v)))
            # Populate rank combobox for Z step
            rank_vals = [str(r.rank) for r in records]
            self.after(0, lambda: self._pl_z_rank_combo.configure(values=rank_vals))
            if rank_vals:
                self.after(0, lambda: self._pl_z_rank.set(rank_vals[0]))
            self.after(0, self._pl_update_table)
            self.after(0, self._pl_refresh_dashboard)
            msg = f"{len(records)} spheroid(s) detected and ranked."
            self.after(0, lambda: self._pl_screen_lbl.configure(
                text=msg, fg=GREEN))
            self._pl_log(msg)
            self._pl_update_step_dots()
            self._pl_save_pipeline_state()
        except Exception as exc:
            self.after(0, lambda: self._pl_screen_lbl.configure(
                text=f"Error: {exc}", fg=RED))
            self._pl_log(f"Screener error: {exc}")

    # ── Step 2 handlers ───────────────────────────────────────────────────────

    def _pl_verify_anchors_thread(self):
        threading.Thread(target=self._pl_verify_anchors, daemon=True).start()

    def _pl_verify_anchors(self):
        try:
            import spheroid_pipeline as _pl
        except ImportError as e:
            messagebox.showerror("Import error", str(e)); return
        if not self._pl_records:
            messagebox.showwarning("No records", "Run Step 1 first."); return
        mosaic = self._pl_mosaic_path.get().strip()
        out    = self._pl_out_dir.get().strip()
        if not mosaic or not out:
            messagebox.showwarning("Missing paths", "Set mosaic ND2 and output dir."); return
        screen_csv = Path(out) / "spheroid_screen_latest.csv"

        corrs = []
        for i, (nd2_var, rk_var, (_, lbl)) in enumerate(
                zip(self._pl_anchor_paths, self._pl_anchor_ranks, self._pl_anchor_frames)):
            nd2_path = nd2_var.get().strip()
            if not nd2_path:
                self.after(0, lambda l=lbl: l.configure(
                    text="skipped — no nd2 selected", fg=YELLOW)); continue
            try:
                exp_rank = int(rk_var.get().strip())
            except ValueError:
                exp_rank = None
            self._pl_log(f"Verifying anchor {i+1}: {Path(nd2_path).name}"
                         + (f" (expected rank {exp_rank})" if exp_rank else " (auto-match)"))
            try:
                dx, dy, ncc, matched_rank, auto_rank = _pl.estimate_offset_from_nd2(
                    Path(mosaic), Path(nd2_path), self._pl_records, Path(out), exp_rank)
                # matched_rank == exp_rank when the operator typed one (forced);
                # the correspondence uses THAT cell, not the NCC best-match.
                rec_m = next((r for r in self._pl_records if r.rank == matched_rank), None)
                if rec_m is not None:
                    mxy = (rec_m.mosaic_x_um, rec_m.mosaic_y_um)
                    corrs.append((mxy, (mxy[0] + dx, mxy[1] + dy)))
                self.after(0, lambda v=matched_rank, var=rk_var: var.set(str(v)))
                # NCC is only a cross-check when a rank is forced: warn on disagreement,
                # never override.
                mism = (exp_rank is not None and auto_rank != exp_rank)
                if exp_rank is not None:
                    txt = (f"forced rank {matched_rank}  dx={dx:+.1f} um  dy={dy:+.1f} um"
                           + (f"   ⚠ NCC best-match is rank {auto_rank} "
                              f"(NCC={ncc:.3f}) — check this anchor!"
                              if mism else f"   (NCC agrees, {ncc:.3f})"))
                else:
                    txt = (f"auto-matched rank {matched_rank}  NCC={ncc:.4f}  "
                           f"dx={dx:+.1f} um  dy={dy:+.1f} um")
                self.after(0, lambda l=lbl, t=txt, mm=mism:
                           l.configure(text=t, fg=(YELLOW if mm else GREEN)))
                self._pl_log(f"  Anchor {i+1}: "
                             + (f"FORCED rank {matched_rank}" if exp_rank is not None
                                else f"auto rank {matched_rank}")
                             + f"; NCC best-match rank {auto_rank} ({ncc:.4f}); "
                             + f"dx={dx:+.1f} dy={dy:+.1f} um"
                             + ("   *** MISMATCH ***" if mism else ""))
            except ValueError as exc:
                self.after(0, lambda l=lbl, e=exc: l.configure(
                    text=f"FAILED: {e}", fg=RED))
                self._pl_log(f"  Anchor {i+1} failed: {exc}")

        if not corrs:
            self.after(0, lambda: self._pl_offset_lbl.configure(
                text="Offset: FAILED — no successful anchor", fg=RED))
            self._pl_log("All anchors failed — cannot proceed.")
            return
        info = _pl.apply_anchor_transform(self._pl_records, corrs)
        if info["mode"] == "similarity":
            self._pl_offset = (info["tx"], info["ty"])
            msg = (f"Transform applied ({info['n']} anchors): scale={info['scale']:.3f}  "
                   f"rot={info['rotation_deg']:+.1f} deg  (flip-aware similarity)")
        elif info["mode"] == "flip1":
            self._pl_offset = (info["tx"], info["ty"])
            msg = ("1 anchor + assumed 180-deg mosaic/stage flip applied (rough). "
                   "Add a 2nd strong anchor for the precise fitted transform.")
        else:
            self._pl_offset = (info.get("dx", 0.0), info.get("dy", 0.0))
            msg = (f"Offset applied: dx={info.get('dx', 0.0):+.1f} um  "
                   f"dy={info.get('dy', 0.0):+.1f} um  (translation)")
        self.after(0, self._pl_update_table)
        self.after(0, self._pl_refresh_dashboard)
        self.after(0, lambda m=msg: self._pl_offset_lbl.configure(text=m, fg=GREEN))
        self._pl_log(msg)
        self._pl_update_step_dots()
        self._pl_save_pipeline_state()

    # ── Step 3 handlers ───────────────────────────────────────────────────────

    def _pl_record_z_thread(self):
        threading.Thread(target=self._pl_record_z, daemon=True).start()

    def _pl_record_z(self):
        try:
            import spheroid_pipeline as _pl
        except ImportError as e:
            messagebox.showerror("Import error", str(e)); return
        rank_s = self._pl_z_rank.get().strip()
        nd2_s  = self._pl_z_nd2_path.get().strip()
        if not rank_s or not nd2_s:
            messagebox.showwarning("Missing input", "Select rank and ND2 file."); return
        try:
            rank   = int(rank_s)
            z_half = float(self._pl_z_half.get())
            z_step = float(self._pl_z_step.get())
        except ValueError:
            messagebox.showerror("Bad value", "Rank, Z-half and Z-step must be numbers."); return
        self._pl_log(f"Reading Z from {Path(nd2_s).name} for rank {rank} ...")
        try:
            _pl.record_z_centre(self._pl_records, rank, Path(nd2_s), z_half, z_step)
            self.after(0, self._pl_update_table)
            self.after(0, self._pl_refresh_dashboard)
            rec = next(r for r in self._pl_records if r.rank == rank)
            msg = f"Rank {rank}: z_centre={rec.z_centre_um:.2f} um"
            self.after(0, lambda: self._pl_z_status_lbl.configure(text=msg, fg=GREEN))
            self._pl_log(msg)
        except Exception as exc:
            self.after(0, lambda: self._pl_z_status_lbl.configure(
                text=f"Error: {exc}", fg=RED))
            self._pl_log(f"Z record error: {exc}")

    def _pl_trigger_autofocus(self):
        import spheroid_pipeline as _pl
        if not self._pl_records:
            messagebox.showwarning("No records", "Run Steps 1-2 first."); return False
        out = self._pl_out_dir.get().strip()
        if not out:
            messagebox.showwarning("No output dir", "Set output directory in Step 1."); return False
        # Writing a fresh trigger set means the operator is starting new work -- drop any
        # leftover abort flag so the macro they run next isn't killed at rank 1.
        self._pl_abort_clear()
        try:
            z_centre = float(self._pl_z_centre.get())
            z_half   = float(self._pl_z_half.get())
            z_step   = float(self._pl_z_step.get())
        except ValueError:
            messagebox.showerror(
                "Bad value", "Middle plane Z, Z half-range, and Z step must be numbers."); return False
        locate = self._pl_locate_only.get()
        if locate:
            z_half = 0.0          # single plane at z_centre -- fast locate pass for re-centering
        excl = getattr(self, "_pl_excluded", set())
        recs = [r for r in self._pl_records if r.rank not in excl]
        if not recs:
            messagebox.showwarning(
                "None selected",
                "Every spheroid is unchecked in the Use column — nothing to trigger."); return False
        # Auto-close each capture in NIS-E by default (Step 3 checkbox); "0" keeps open.
        close_after = "0" if self._pl_s3_keep_open.get() else "1"
        try:
            n = _pl.trigger_autofocus_all(recs, Path(out) / "autofocus",
                                          z_centre=z_centre, z_half=z_half, z_step=z_step,
                                          close_after=close_after)
        except Exception as e:
            # Same defect class as _pl_regen_triggers: this mkdirs, writes N trigger
            # files and session.ini. An OSError (share drop / read-only run folder) or
            # the duplicate-coordinate ValueError previously escaped uncaught, and the
            # stale-trigger unlink has ALREADY run -- leaving a partial set on disk.
            messagebox.showerror("Trigger write failed",
                                 f"Could not write triggers:\n\n{type(e).__name__}: {e}")
            self.after(0, lambda m=str(e): self._pl_af_status_lbl.configure(
                text=f"Trigger write FAILED — {m}", fg=RED))
            self._pl_log(f"Step 3: trigger write FAILED ({type(e).__name__}): {e}")
            return False
        mode = "1-plane LOCATE" if locate else f"Z-stack +/-{z_half:g}/{z_step:g} um"
        skipped = len(self._pl_records) - len(recs)
        extra = f" ({skipped} unchecked skipped)" if skipped else ""
        self._pl_log(f"Step 3: {n} {mode} triggers written to autofocus/ (centre={z_centre}){extra}")
        self.after(0, lambda m=mode, e=extra: self._pl_af_status_lbl.configure(
            text=f"{n} {m} triggers written{e}.", fg=YELLOW))
        self._pl_save_pipeline_state()
        return True

    def _pl_poll_af_thread(self):
        threading.Thread(target=self._pl_poll_af, daemon=True).start()

    def _pl_poll_af(self):
        import spheroid_pipeline as _pl
        out = self._pl_out_dir.get().strip()
        if not out or not self._pl_records:
            return
        done = _pl.poll_autofocus_done(self._pl_records, Path(out) / "autofocus")
        try:
            z_half = float(self._pl_z_half.get())
            z_step = float(self._pl_z_step.get())
        except ValueError:
            z_half, z_step = 90.0, 10.0
        _pl.apply_autofocus_results(self._pl_records, done, z_half, z_step)
        n_done  = len(done)
        n_total = len(self._pl_records)
        msg = f"{n_done}/{n_total} autofocus captures complete."
        fg = GREEN if n_done == n_total else (YELLOW if n_done > 0 else SUBTEXT)
        self.after(0, lambda: self._pl_af_status_lbl.configure(text=msg, fg=fg))
        self.after(0, self._pl_update_table)
        self.after(0, self._pl_refresh_dashboard)
        self._pl_log(f"Step 3 poll: {msg}")
        self._pl_update_step_dots()

    def _pl_recenter_thread(self):
        # Recenter now REWRITES the on-disk trigger set (auto-retrigger), so it must NOT
        # run while a Capture's macro is mid-sweep -- clearing/rewriting the triggers the
        # running daemon is consuming would corrupt the batch (remaining ranks captured
        # oc-less / wrong-channel, good captures overwritten). Hold the same single-in-
        # flight lock Capture uses, and release it in _pl_recenter_run's finally.
        if getattr(self, "_pl_dispatch_busy", False):
            messagebox.showinfo("Dispatcher busy",
                                "A capture is running via the dispatcher — wait for it to finish before re-centering."); return
        self._pl_dispatch_busy = True
        self._pl_af_status_lbl.configure(text="Re-alignment in process …", fg=YELLOW)
        try:
            threading.Thread(target=self._pl_recenter_run, daemon=True).start()
        except Exception as e:
            # Same guard as _pl_capture_thread: if the worker never starts, its finally
            # (the only clearer) never runs -- release the lock so the GUI isn't wedged.
            self._pl_dispatch_busy = False
            self._pl_log(f"Recenter: failed to start worker thread: {e}")
            self._pl_af_status_lbl.configure(text="Recenter: could not start (see log).", fg=RED)

    def _pl_recenter_run(self):
        try:
            self._pl_recenter()
        finally:
            self._pl_dispatch_busy = False

    def _pl_recenter(self):
        import spheroid_pipeline as _pl
        if not self._pl_records:
            messagebox.showwarning("No records", "Run Steps 1-3 (capture) first."); return
        # Find the folder that actually holds the *_zstack.nd2 captures, robustly, so
        # the user never has to hand-set a dir: session.ini nd2_dir -> <Output dir>/nd2
        # -> the Step-4 ND2 field (or its /nd2 child).
        def _caps(d):
            try:
                return bool(d) and Path(d).is_dir() and any(Path(d).glob("*_zstack.nd2"))
            except Exception:
                return False
        nd2_dir = ""
        try:
            import configparser
            cp = configparser.ConfigParser(); cp.read(_pl.SESSION_INI)
            c = cp.get("paths", "nd2_dir", fallback="").strip()
            if _caps(c):
                nd2_dir = c
        except Exception:
            pass
        if not nd2_dir:
            out = self._pl_out_dir.get().strip()
            if out and _caps(str(Path(out) / "nd2")):
                nd2_dir = str(Path(out) / "nd2")
        if not nd2_dir:
            c = self._pl_nd2_out_dir.get().strip()
            if _caps(c):
                nd2_dir = c
            elif c and _caps(str(Path(c) / "nd2")):
                nd2_dir = str(Path(c) / "nd2")
        if not nd2_dir:
            messagebox.showwarning(
                "No captures",
                "Couldn't find any *_zstack.nd2. Set the Step 1 Output dir to your run "
                "folder (e.g. ...\\0622) and capture first."); return
        excl = getattr(self, "_pl_excluded", set())
        targets = sorted((r for r in self._pl_records if r.rank not in excl),
                         key=lambda r: r.rank)
        n_s = self._pl_recenter_n.get().strip().lower()
        if n_s not in ("", "all", "0"):
            try:
                n = int(n_s)
            except ValueError:
                messagebox.showwarning("Bad number", "'# to use' must be a whole number or 'all'."); return
            if n > 0:
                targets = targets[:n]
        try:
            res, skipped = _pl.recenter_from_captures(targets, Path(nd2_dir),
                                             flip_x=self._pl_recenter_flipx.get(),
                                             flip_y=self._pl_recenter_flipy.get())
        except Exception as exc:
            self._pl_log(f"Re-center error: {exc}")
            self.after(0, lambda e=exc: self._pl_af_status_lbl.configure(
                text=f"Re-center error: {e}", fg=RED)); return
        if not res and not skipped:
            self.after(0, lambda: self._pl_af_status_lbl.configure(
                text="Re-center: no captures found in the ND2 dir (capture first).", fg=RED)); return
        self._pl_recentered = {rank: (dx, dy) for rank, dcol, drow, dx, dy in res}
        for rank, dcol, drow, dx, dy in res:
            self._pl_log(f"  rank {rank}: centroid ({dcol:+.0f},{drow:+.0f})px -> nudge ({dx:+.1f},{dy:+.1f}) um")
        # Skipped ranks keep their PRE-recenter verified_x/y -- surface this loudly,
        # since silently leaving them uncorrected is exactly what produced the
        # "still off-center after 3 recenters" bug (7/10 spheroids were silently
        # skipped every round by the old MAD-based threshold).
        for rank, reason in skipped:
            self._pl_log(f"  rank {rank}: SKIPPED (not corrected) -- {reason}")
        self.after(0, self._pl_update_table)
        self.after(0, self._pl_refresh_dashboard)
        deltas = "   ".join(f"#{rank} Δ({dx:+.0f},{dy:+.0f})" for rank, dcol, drow, dx, dy in res)
        warn = f"   ⚠ {len(skipped)} SKIPPED (not corrected): {[r for r, _ in skipped]}" if skipped else ""
        # Auto-update the on-disk triggers with the corrected coords so the operator
        # doesn't have to click a separate Trigger: recenter_from_captures already nudged
        # verified_x/y in-place, so regen writes the corrected positions (honoring Locate-
        # only). The next Capture then re-images at the corrected centres.
        # _pl_regen_triggers returns True BOTH when it wrote triggers AND when it SKIPPED
        # (no Output dir -> kept existing on-disk triggers). Distinguish them: without an
        # Output dir the correction lives only in memory and is NOT persisted, so don't
        # claim it was written.
        out_set = bool(self._pl_out_dir.get().strip())
        retrig = self._pl_regen_triggers()
        persisted = retrig and out_set
        failed = not retrig
        if persisted:
            self._pl_log("Re-center: triggers updated with corrected coordinates.")
            upd = "  Triggers updated -- click Capture to re-image."
        elif retrig and not out_set:
            self._pl_log("Re-center: correction is IN MEMORY only -- no Output dir set, triggers NOT written.")
            upd = "  Correction in MEMORY only -- set Step-1 Output dir + Capture to persist it."
        else:
            # Genuine regen failure (bad Step-3 Z values / write error). _pl_regen_triggers
            # routes its detailed reason to the Step-4 label, so log a Step-3 pointer here too.
            self._pl_log("Re-center: trigger update FAILED (bad Step-3 Z values or write error) -- correction NOT persisted.")
            upd = "  Trigger update FAILED (see log) -- correction NOT persisted."
        self.after(0, lambda n=len(res), d=deltas, w=warn, u=upd, ok=persisted, bad=failed: self._pl_af_status_lbl.configure(
            text=(f"Re-aligned {n} from captures (Δum):  {d}.{u}"
                  f"  (toggle flip X/Y if more off-center).{w}"),
            fg=(RED if bad else (GREEN if (ok and not warn) else YELLOW))))
        self._pl_log(f"Re-center: adjusted {len(res)} spheroids from their captures"
                     + (f"; {len(skipped)} SKIPPED (see above)." if skipped else "."))

    def _pl_capture_thread(self):
        """Step 3 "Capture": write triggers from the checked rows on the MAIN thread
        (honoring Locate-only + Keep-open; _pl_trigger_autofocus shows modal dialogs on
        failure, so it must not run off-thread), then on a worker thread inject the 1050 nm
        depth OC and run the z-stack macro (action 2). Full stack by default
        (ND_RunZSeriesExp, the confocal path); a single plane only if Locate-only is ticked.
        Captures land in the plain nd2/ dir that 'Recenter' reads."""
        if getattr(self, "_pl_dispatch_busy", False):
            messagebox.showinfo("Dispatcher busy",
                                "A macro is already running via the dispatcher — wait for it to finish."); return
        # The 1050 nm Capture is 2-PHOTON: the macro removes the A1 interlock and selects a
        # confocal/2P OC. Both crash the grabber if the A1 is OFF, and without the interlock
        # cleared + NIS-E in 2-Photon mode the frame comes off the widefield camera (a DIA
        # capture, not the 1050 confocal). Require the A1-on confirmation, like Step 4.
        if not self._pl_pa_a1on.get():
            messagebox.showwarning(
                "A1 not confirmed",
                "The 1050 nm Capture is 2-photon: it removes the A1 interlock and needs the "
                "A1/2P powered ON.\n\nTick 'A1 powered ON (2P)' first, and make sure NIS-E is "
                "switched to the 2-Photon acquisition mode (not Flash 4 + BF)."); return
        if not self._pl_trigger_autofocus():
            return   # nothing selected / bad value / write failed -- dialog already shown
        self._pl_abort_clear()          # a stale flag must not kill the run we're starting
        self._pl_dispatch_busy = True
        self.after(0, lambda: self._pl_af_status_lbl.configure(
            text="Capture: 1050 OC + z-stack starting…", fg=YELLOW))
        try:
            threading.Thread(target=self._pl_capture_1050, daemon=True).start()
        except Exception as e:
            # If the worker never starts, its finally (the only clearer) never runs --
            # release the lock here so the GUI isn't left wedged.
            self._pl_dispatch_busy = False
            self._pl_log(f"Capture: failed to start worker thread: {e}")
            self.after(0, lambda: self._pl_af_status_lbl.configure(
                text="Capture: could not start (see log).", fg=RED))

    def _pl_capture_1050(self, z_half=None, z_step=None, oc_name=None, tag=None,
                         load_captured=False):
        """Worker: inject the 1050 nm depth OC (the SAME OC Step 4's Validation uses --
        see _pl_prepa_checked_ocs) into the triggers _pl_trigger_autofocus just wrote, then
        dispatch the z-stack macro (action 2). Preserves the z geometry the triggers carry,
        so a full stack (Locate-only OFF) runs ND_RunZSeriesExp -- the confocal path that
        actually images the 1050 channel; the single-plane path (Locate-only ON) now also
        images the OC's channel via the macro's 1-plane-ND-series fix. Captures land in the
        plain nd2/ dir 'Recenter' reads. Guarded trigger restore + busy release in finally
        (models _pl_prepa_capture_ocs)."""
        import configparser, time, spheroid_pipeline as _pl
        # Defaults reproduce the Step-3 Capture exactly; the Find Z-Center "Capture Z-Search"
        # passes a wide/coarse (z_half, z_step) + search oc/tag and load_captured=True to
        # auto-show the filmstrips. When z_half/z_step are given they OVERRIDE the trigger's
        # geometry in the injection loop below (z_centre + close_after still come from it).
        if oc_name is None:
            oc_name, tag = "1050nm_Galvo_561nm_NDD2_WC", "1050nm_depth"   # == _pl_prepa_checked_ocs 1050
        elif tag is None:
            tag = "1050nm_zsearch"
        light_path = "2-Photon"   # NIS-E light path for the A1/2P acquisition (macro SelectLightPath)
        nosepiece_20x = 1         # 20X nosepiece POSITION for the macro's Stg_SetNosepiecePosition (0-based,
                                  # same as the nd2 NosPosition: 10X=0, 20X=1 -- confirmed from the 10X mosaic
                                  # AND a 20X capture). The earlier "10X" result was the WRONG function name
                                  # (NosepiecePosition, not callable -> "Cannot Evaluate"), not an index issue.
        wd = None
        triggers = []
        status = None
        try:
            # work_dir was just synced by _pl_trigger_autofocus -> trigger_autofocus_all.
            cp = configparser.ConfigParser()
            try:
                if _pl.SESSION_INI.exists():
                    cp.read(_pl.SESSION_INI)
            except Exception:
                pass
            work = cp.get("paths", "work_dir", fallback="").strip() or self._pl_trigger_dir.get().strip()
            if not work:
                o = self._pl_out_dir.get().strip()
                work = str(Path(o) / "autofocus") if o else ""
            if not work:
                self.after(0, lambda: self._pl_af_status_lbl.configure(
                    text="Capture: no work dir -- set the Output dir in Step 1.", fg=RED)); return
            wd = Path(work)
            for p in sorted(wd.glob(f"{_pl.AF_TRIGGER_PREFIX}*.ini")):
                d = _pl.parse_ini(p.read_text(encoding="utf-8"))
                if {"rank", "stage_x", "stage_y"} <= d.keys():
                    triggers.append(d)
            if not triggers:
                self.after(0, lambda: self._pl_af_status_lbl.configure(
                    text="Capture: no af_trigger_*.ini were written.", fg=RED)); return
            # Inject the OC into each trigger, PRESERVING the z geometry + close_after
            # _pl_trigger_autofocus already wrote (a full stack stays a full stack; a
            # Locate-only single plane stays single) -- unless the caller overrides
            # z_half/z_step (the Z-Search wide/coarse sweep), in which case those win.
            z_override = {"z_half": z_half, "z_step": z_step}
            for d in triggers:
                lines = ["[spheroid]", f"rank={d['rank']}",
                         f"spheroid_id={d.get('spheroid_id', '')}",
                         f"stage_x={d['stage_x']}", f"stage_y={d['stage_y']}"]
                for k in ("z_centre", "z_half", "z_step", "close_after"):
                    if z_override.get(k) is not None:
                        lines.append(f"{k}={z_override[k]}")
                    elif k in d:
                        lines.append(f"{k}={d[k]}")
                lines.append(f"oc={oc_name}")
                # 2-photon prep (macro does these ONCE, before the first SelectOptConf, gated on
                # a1_on): switch the light path to 2-Photon (SelectOptConf alone leaves the camera
                # on Flash 4 + BF -> a DIA capture) and clear the A1 interlock so the 2P scanner
                # can image. a1_on was confirmed in _pl_capture_thread; mirrors pa_setup. Step-4
                # validation writes NONE of these, so it is unchanged.
                lines.append(f"light_path={light_path}")
                if nosepiece_20x > 0:
                    lines.append(f"nosepiece_pos={nosepiece_20x}")   # switch to 20X before imaging
                lines.append("a1_on=1")
                lines.append("remove_interlock=1")
                _pl._atomic_write_crlf(wd / f"af_trigger_{int(d['rank']):02d}.ini", lines)
            # Dispatch the z-stack once (the macro loops all triggers). We deliberately do
            # NOT redirect session.ini nd2_dir: the captures land in the plain nd2/ dir the
            # trigger write set, which is where 'Recenter' looks for *_zstack.nd2.
            self.after(0, lambda n=len(triggers): self._pl_af_status_lbl.configure(
                text=f"Capture: {tag} for {n} spheroid(s)…", fg=YELLOW))
            self._pl_log(f"Capture: dispatched 1050 z-stack ({oc_name}) for {len(triggers)} spheroid(s)")
            # Same no-overwrite guarantee for the plain nd2/ dir these land in. A Step-3
            # Locate pass writes the same <id>_zstack.nd2 name as a full stack, so without
            # this a 1-plane locate silently destroys a 61-plane capture.
            try:
                _cp = configparser.ConfigParser(); _cp.read(_pl.SESSION_INI)
                _nd = _cp.get("paths", "nd2_dir", fallback="").strip()
            except Exception:
                _nd = ""
            if _nd:
                self._pl_archive_existing_nd2(
                    _nd, [d.get("spheroid_id", "") for d in triggers], f"Capture {tag}")
            done = wd / "cmd_done.ini"
            try:
                done.unlink()
            except FileNotFoundError:
                pass
            _pl._atomic_write_crlf(wd / "cmd.ini",
                                   ["[command]", f"action=zstack_{tag}", "action_id=2"])
            status = None
            for _ in range(1800):
                if done.exists():
                    cp2 = configparser.ConfigParser()
                    try:
                        cp2.read(done)
                        status = cp2.get("command", "status",
                                         fallback=cp2.get("spheroid", "status", fallback="?"))
                    except Exception:
                        time.sleep(1.0); continue
                    if status == "?":
                        # cmd_done.ini is built key-by-key with status written LAST, so a
                        # file that parses without one is still being written -- not a
                        # finished run with an unknown result. Breaking here reported "?"
                        # and aborted a post-PA sequence whose capture had actually
                        # succeeded (2026-08-07, cost the 1050 nm pass).
                        time.sleep(1.0); continue
                    break
                if self._pl_abort_requested(wd):
                    status = "aborted"; break
                time.sleep(1.0)
            self._pl_log(f"Capture: 1050 z-stack -> {status or 'timeout'}")
            if status != "ok":
                msg = (f"Capture: 1050 z-stack -> {status}" if status
                       else "Capture: timed out (is nis_macro_dispatcher.mac running?)")
                self.after(0, lambda m=msg: self._pl_af_status_lbl.configure(text=m, fg=RED)); return
            self.after(0, lambda: self._pl_af_status_lbl.configure(
                text="Capture: 1050 stacks captured -- click 'Recenter'.", fg=GREEN))
            self.after(0, self._pl_poll_af_thread)   # refresh the table/dashboard from af_done
            if load_captured:
                # Auto-show the just-captured filmstrips so Find Z-Center's highlights land
                # somewhere visible. Captures kept the plain nd2/ dir (nd2_dir was never
                # redirected), named {sid}_zstack.nd2 exactly as Recenter/Validation expect.
                nd2_cfg = cp.get("paths", "nd2_dir", fallback="").strip()
                base_nd2 = Path(nd2_cfg) if nd2_cfg else wd.parent / "nd2"
                items = []
                for d in triggers:
                    sid = (str(d.get("spheroid_id") or "").strip()
                           or f"rank{int(d['rank']):02d}")
                    items.append((f"{sid}_zstack.nd2", str(base_nd2 / f"{sid}_zstack.nd2")))
                self.after(0, lambda it=items: self._pl_zv_load_captured(it, "Z-Search"))
        except Exception as e:
            self._pl_log(f"Capture: FAILED ({type(e).__name__}): {e}")
            self.after(0, lambda m=str(e): self._pl_af_status_lbl.configure(
                text=f"Capture: FAILED -- {m}", fg=RED))
        finally:
            # Restore the plain (oc-less) triggers the macro consumed after each capture,
            # so 'Recenter' / PA Points don't find an empty autofocus/ folder. Re-persists
            # the SAME coords already written -- no staleness.
            #
            # ONLY when the macro actually REPORTED (status set): on a poll timeout (status
            # stays None) the daemon may still be mid-sweep, and rewriting its triggers
            # oc-less now would corrupt the remaining ranks (same class as the Recenter
            # race). status is set on ok/aborted, None on timeout / an early exception.
            try:
                if wd is not None and triggers and status is not None:
                    for d in triggers:
                        lines = ["[spheroid]", f"rank={d['rank']}",
                                 f"spheroid_id={d.get('spheroid_id', '')}",
                                 f"stage_x={d['stage_x']}", f"stage_y={d['stage_y']}"]
                        for k in ("z_centre", "z_half", "z_step"):
                            if k in d:
                                lines.append(f"{k}={d[k]}")
                        _pl._atomic_write_crlf(wd / f"af_trigger_{int(d['rank']):02d}.ini", lines)
            except Exception as e:
                self._pl_log(f"Capture: WARNING - could not restore triggers: {e}")
            self._pl_dispatch_busy = False

    # ── Find Z-Center handlers ────────────────────────────────────────────────

    def _pl_zc_capture_thread(self):
        """Find Z-Center "Capture Z-Search": the EXACT Step-3 Capture idiom
        (_pl_capture_thread) but drives a WIDE/COARSE stack from _pl_zc_half/_pl_zc_step
        and the search OC, so each spheroid's full through-focus is imaged. Same A1/2P
        interlock guards; non-blocking; auto-loads the captures into the viewer on success."""
        if getattr(self, "_pl_dispatch_busy", False):
            messagebox.showinfo("Dispatcher busy",
                                "A macro is already running via the dispatcher — wait for it to finish."); return
        if not self._pl_pa_a1on.get():
            messagebox.showwarning(
                "A1 not confirmed",
                "The Z-Search capture is 2-photon: it removes the A1 interlock and needs the "
                "A1/2P powered ON.\n\nTick 'A1 powered ON (2P)' first, and make sure NIS-E is "
                "switched to the 2-Photon acquisition mode (not Flash 4 + BF)."); return
        try:
            z_half = float(self._pl_zc_half.get()); z_step = float(self._pl_zc_step.get())
        except ValueError:
            messagebox.showerror("Bad value", "Search half-range and step must be numbers."); return
        oc = self._pl_zc_oc.get().strip()
        if not oc:
            messagebox.showwarning("No OC", "Set a search OC for the Z-Search capture."); return
        if not self._pl_trigger_autofocus():
            return   # nothing selected / bad value / write failed -- dialog already shown
        self._pl_abort_clear()          # a stale flag must not kill the run we're starting
        self._pl_dispatch_busy = True
        self.after(0, lambda: self._pl_zc_status_lbl.configure(
            text=f"Capture Z-Search: wide +/-{z_half:g}/{z_step:g} um 1050 stack starting…", fg=YELLOW))
        try:
            threading.Thread(target=self._pl_capture_1050,
                             kwargs=dict(z_half=z_half, z_step=z_step, oc_name=oc,
                                         tag="1050nm_zsearch", load_captured=True),
                             daemon=True).start()
        except Exception as e:
            self._pl_dispatch_busy = False
            self._pl_log(f"Capture Z-Search: failed to start worker thread: {e}")
            self.after(0, lambda: self._pl_zc_status_lbl.configure(
                text="Capture Z-Search: could not start (see log).", fg=RED))

    def _pl_zc_find_thread(self):
        threading.Thread(target=self._pl_zc_find, daemon=True).start()

    def _pl_zc_search_dir(self) -> str:
        """Resolve the folder holding the *_zstack.nd2 search captures the same way
        Recenter does: session.ini nd2_dir -> work_dir/nd2 -> Output dir/nd2 -> the
        Step-4 ND2 field (or its /nd2 child). Returns "" if none carry captures."""
        import configparser, spheroid_pipeline as _pl
        def _caps(d) -> bool:
            try:
                return bool(d) and Path(d).is_dir() and any(Path(d).glob("*_zstack.nd2"))
            except Exception:
                return False
        try:
            cp = configparser.ConfigParser(); cp.read(_pl.SESSION_INI)
            c = cp.get("paths", "nd2_dir", fallback="").strip()
            if _caps(c):
                return c
            wd = cp.get("paths", "work_dir", fallback="").strip()
            if wd and _caps(str(Path(wd).parent / "nd2")):
                return str(Path(wd).parent / "nd2")
        except Exception:
            pass
        out = self._pl_out_dir.get().strip()
        if out and _caps(str(Path(out) / "nd2")):
            return str(Path(out) / "nd2")
        c = self._pl_nd2_out_dir.get().strip()
        if _caps(c):
            return c
        if c and _caps(str(Path(c) / "nd2")):
            return str(Path(c) / "nd2")
        return ""

    def _pl_zc_find(self):
        """Peak-detect each spheroid's focal plane from the wide search captures
        (find_zcenter_all), store the results for the viewer highlight, apply confident
        centres (+ the OC->PA dZ offset) to z_centre when the toggle is on, and auto-load
        the captures so the found planes are highlighted."""
        import spheroid_pipeline as _pl
        if not self._pl_records:
            self.after(0, lambda: self._pl_zc_status_lbl.configure(
                text="Find Centers: run Steps 1-3 (capture) first.", fg=RED)); return
        search_dir = self._pl_zc_search_dir()
        if not search_dir:
            self.after(0, lambda: self._pl_zc_status_lbl.configure(
                text="Find Centers: no *_zstack.nd2 found -- Capture Z-Search first.", fg=RED)); return
        try:
            offset = float(self._pl_zc_offset.get())
        except ValueError:
            self.after(0, lambda: self._pl_zc_status_lbl.configure(
                text="Find Centers: OC->PA dZ offset must be a number.", fg=RED)); return
        try:
            res = _pl.find_zcenter_all(self._pl_records, search_dir,
                                       prominence_min=ZC_PROMINENCE_MIN)
        except Exception as e:
            self._pl_log(f"Find Centers: FAILED ({type(e).__name__}): {e}")
            self.after(0, lambda m=str(e): self._pl_zc_status_lbl.configure(
                text=f"Find Centers: FAILED -- {m}", fg=RED)); return
        self._pl_zc_found = res
        apply = self._pl_zc_apply.get()
        by_sid = {str(r.spheroid_id): r for r in self._pl_records}
        n_conf = 0
        for sid, result in res.items():
            conf = bool(result.get("confidence"))
            n_conf += int(conf)
            self._pl_log(
                f"Find Z-Center: sid={sid} plane={result.get('plane_index')} "
                f"z_centre_um={float(result.get('z_centre_um', 0.0)):.1f} conf={conf} "
                f"prom={float(result.get('prominence', 0.0)):.2f} ({result.get('reason')})")
            if apply and conf:
                r = by_sid.get(str(sid))
                if r is not None:
                    r.z_centre_um = float(result["z_centre_um"]) + offset
        self.after(0, self._pl_update_table)
        self.after(0, self._pl_refresh_dashboard)
        applied = f" ({n_conf} applied +{offset:g} um)" if apply else " (not applied)"
        self.after(0, lambda n=len(res), c=n_conf, a=applied: self._pl_zc_status_lbl.configure(
            text=f"Find Centers: {c}/{n} confident{a}. Highlights shown in Captured Z-Stacks.",
            fg=(GREEN if c else YELLOW)))
        self.after(0, self._pl_zv_autoload)   # re-render so found-centre highlights appear

    def _pl_zc_sid_from_key(self, key: str) -> str:
        """Extract a spheroid id from a viewer entry key, matching how _pl_zv entries are
        keyed: "{sid}_zstack.nd2", the "Add..." "name  [parent]" variant, or the phase-
        qualified "{phase}/{tag}/{sid}" Validation display name."""
        if not key:
            return ""
        s = key.split("  [")[0].strip()
        if "/" in s:
            s = s.rsplit("/", 1)[-1]
        if s.endswith(".nd2"):
            s = s[:-4]
        if s.endswith("_zstack"):
            s = s[:-7]
        return s.strip()

    # ── Step 4 handlers ───────────────────────────────────────────────────────

    def _info(self, parent, short, detail, bg=None, fg=None, wrap=460):
        """A one-line label plus a clickable (i) that toggles the full explanation below it.

        Cards accumulate caveats -- which job fires, why a folder matters, what a failure
        means -- and printing all of it permanently turns the panel into a wall of small
        text that gets skimmed and then ignored. Keep the operating line visible, put the
        reasoning one click away. Nothing is deleted, only folded."""
        bg = bg or BG2
        row = tk.Frame(parent, bg=bg); row.pack(fill="x", anchor="w")
        tk.Label(row, text=short, bg=bg, fg=fg or SUBTEXT, font=("Segoe UI", 8),
                 justify="left", wraplength=wrap).pack(side="left")
        btn = tk.Label(row, text=" ⓘ", bg=bg, fg=BLUE,
                       font=("Segoe UI", 10, "bold"), cursor="hand2")
        btn.pack(side="left")
        det = tk.Label(parent, text=detail, bg=bg, fg=SUBTEXT, font=("Segoe UI", 8),
                       justify="left", wraplength=wrap, anchor="w")

        def toggle(_=None):
            if det.winfo_ismapped():
                det.pack_forget()
                btn.configure(fg=BLUE)
            else:
                # after=row keeps the detail attached to ITS line; without it Tk appends
                # to the end of the parent and the text lands under an unrelated widget.
                det.pack(fill="x", anchor="w", padx=(14, 0), after=row)
                btn.configure(fg=YELLOW)
        btn.bind("<Button-1>", toggle)
        return row

    def _pl_pa_card(self, parent, title, sel_var):
        """One PA-macro card: a [select] checkbox header + a body frame for params.
        Returns (card, body); the caller fills `body` then calls _pl_pa_card_buttons."""
        card = tk.Frame(parent, bg=BG2, bd=1, relief="solid")
        card.pack(fill="x", padx=4, pady=3)
        hdr = tk.Frame(card, bg=BG2); hdr.pack(fill="x", padx=6, pady=(3, 0))
        tk.Checkbutton(hdr, text=title, variable=sel_var, bg=BG2, fg=MAUVE,
                       selectcolor=SURFACE, activebackground=BG2, activeforeground=TEXT,
                       font=("Segoe UI", 9, "bold")).pack(side="left")
        body = tk.Frame(card, bg=BG2); body.pack(fill="x", padx=24, pady=(0, 2))
        return card, body

    def _pl_pa_card_buttons(self, card, action, action_id, with_reload=True):
        btns = tk.Frame(card, bg=BG2); btns.pack(fill="x", padx=24, pady=(2, 6))
        if with_reload:
            self._btn(btns, "Reload", lambda a=action: self._pl_pa_reload(a),
                      SURFACE2, TEXT, side="left")
        self._btn(btns, "Run", lambda a=action, i=action_id: self._pl_pa_run_macro(a, i),
                  BLUE, "#1e1e2e", side="left")

    def _pl_pa_power_clamp(self, *_):
        """Hard-cap PA Setup's Power % at MAX_PA_ACTIVATION_POWER_PCT (30%) -- 80%
        visibly damaged a spheroid on 2026-07-07. This clamps the GUI field only;
        it can't enforce the real laser power, since step3_zstack_PA is a manual
        JOB run in NIS-E and doesn't read pa_trigger.ini's power_pct at all."""
        s = self._pl_pa_power.get().strip()
        if not s:
            return
        try:
            v = float(s)
        except ValueError:
            return   # mid-edit (e.g. a lone "-" or "."); let them keep typing
        if v > MAX_PA_ACTIVATION_POWER_PCT:
            self._pl_pa_power.set(str(MAX_PA_ACTIVATION_POWER_PCT))

    def _pl_pa_sync_save_dir(self, *_):
        """Keep the PA folder pointed at <Output directory>/pa and make sure it exists.

        Auto-fills whenever the box is empty or still holds a previous run's auto value,
        so moving to a new run folder carries the PA folder with it -- the 0807/pa vs
        0807_2/pa mismatch on 2026-08-07 came from a stale hand-set path. A path the
        operator typed themselves is left alone: `_pl_pa_save_dir_auto` records what this
        method last set, so anything else is treated as deliberate."""
        out = self._pl_out_dir.get().strip()
        if not out:
            return
        want = (Path(out) / "pa")
        cur = self._pl_pa_save_dir.get().strip()
        auto = getattr(self, "_pl_pa_save_dir_auto", "")
        if cur and cur != auto:
            return                      # operator-chosen -- do not overwrite
        self._pl_pa_save_dir.set(want.as_posix())
        self._pl_pa_save_dir_auto = want.as_posix()
        try:
            want.mkdir(parents=True, exist_ok=True)
        except Exception as e:
            self._pl_log(f"PA folder: could not create {want.as_posix()} ({e})")

    def _pl_pa_effective_save_dir(self):
        """The folder the PA stacks should land in: the Step 4 'PA save folder' box when
        set, else <out_dir>/pa. Returns a Path, or None when neither is available."""
        explicit = self._pl_pa_save_dir.get().strip()
        if explicit:
            return Path(explicit)
        base = self._pl_out_dir.get().strip()
        return (Path(base) / "pa") if base else None

    def _pl_pa_prepare_save_dir(self, quiet=True):
        """Ensure the PA folder exists and is recorded in pa_trigger.ini.

        Called automatically before every PA launch -- there is no button for it. It used
        to be a manual 'Prepare folder' click, which meant the one run where it was
        forgotten wrote into a stale path from a previous session.

        It DELIBERATELY does not re-point NIS-E's ND Acquisition 'Save to File' (an earlier
        build dispatched action 11 to do that). The step3_zstack_PA JOB writes via its own
        'Alternative Storage Location' task and ignores the ND dialog, so pointing the ND
        target at pa/ achieved nothing for the JOB while making EVERY ND acquisition drop a
        second copy there -- on 2026-08-07 each validation capture left a byte-identical
        duplicate as PA.nd2, PA001.nd2, ... in the PA output folder.

        Returns True when the folder is ready."""
        save_dir = self._pl_pa_effective_save_dir()
        if save_dir is None:
            if not quiet:
                messagebox.showwarning("PA folder",
                                       "Set an Output directory in Step 1 first.")
            self._pl_log("PA folder: no Output directory -- cannot prepare the pa/ folder.")
            return False
        try:
            save_dir.mkdir(parents=True, exist_ok=True)
        except Exception as exc:
            self._pl_log(f"PA folder: could not create {save_dir.as_posix()} ({exc})")
            if not quiet:
                messagebox.showerror("PA folder", f"Could not create {save_dir}:\n{exc}")
            return False
        if self._pl_pa_write_trigger() is None:      # records save_dir/save_prefix
            return False
        self._pl_pa_savepath_applied = save_dir.as_posix()
        try:
            self._pl_pa_save_lbl.configure(
                text=f"Folder ready: {save_dir.as_posix()}   "
                     f"(set the SAME path on the job's 'Alternative Storage Location')",
                fg=GREEN)
        except Exception:
            pass
        self._pl_log(f"PA folder ready: {save_dir.as_posix()} "
                     f"(prefix '{self._pl_pa_save_prefix.get().strip() or 'PA'}')")
        return True

    def _pl_pa_write_trigger(self):
        """Write pa_trigger.ini (all PA-macro params) into the session work dir and
        refresh session.ini (work_dir + macro_dir for the dispatcher). Returns the
        work_dir Path, or None on error (after showing a message)."""
        import configparser, spheroid_pipeline as _pl
        save_dir = self._pl_pa_effective_save_dir()
        cp0 = configparser.ConfigParser()
        try:
            if _pl.SESSION_INI.exists():
                cp0.read(_pl.SESSION_INI)
        except Exception:
            pass
        work = cp0.get("paths", "work_dir", fallback="").strip() or self._pl_trigger_dir.get().strip()
        if not work:
            messagebox.showwarning("No work dir",
                                   "Run Step 3 (Trigger) first, or set the Trigger dir."); return None
        wd = Path(work)
        try:
            wd.mkdir(parents=True, exist_ok=True)
            if save_dir:
                save_dir.mkdir(parents=True, exist_ok=True)
        except Exception as exc:
            messagebox.showerror("Dir error", str(exc)); return None
        # Atomic CRLF write (NIS-E Int_GetKeyString needs Windows line endings, and
        # the dispatcher may poll mid-write) -- same guarantee as cmd.ini/session.ini.
        pa_lines = [
            "[photoactivation]",
            f"job={self._pl_pa_job.get().strip()}",
            f"optical_config={self._pl_pa_oc.get().strip()}",
            f"power_pct={self._pl_pa_power.get().strip()}",
            f"well={self._pl_well_id.get().strip()}",
            f"loops={self._pl_pa_loops.get().strip()}",
            f"zoom={self._pl_pa_zoom.get().strip()}",
            f"dichroic_out={'1' if self._pl_pa_dichroic.get() else '0'}",
            f"remove_interlock={'1' if self._pl_pa_interlock.get() else '0'}",
            f"a1_on={'1' if self._pl_pa_a1on.get() else '0'}",
            f"save_dir={save_dir.as_posix() if save_dir else ''}",
            f"save_prefix={self._pl_pa_save_prefix.get().strip() or 'PA'}",
        ]
        try:
            _pl._atomic_write_crlf(wd / "pa_trigger.ini", pa_lines)
        except Exception as exc:
            messagebox.showerror("Write failed", str(exc)); return None
        try:
            nd2 = cp0.get("paths", "nd2_dir", fallback="").strip()
            lines = ["[paths]", f"work_dir={wd.as_posix()}"]
            if nd2:
                lines.append(f"nd2_dir={nd2}")
            lines.append(f"macro_dir={_pl.MACRO_DIR.as_posix()}")
            _pl.SESSION_INI.parent.mkdir(parents=True, exist_ok=True)
            _pl._atomic_write_crlf(_pl.SESSION_INI, lines)
        except Exception:
            pass
        return wd

    def _pl_pa_reload(self, action):
        wd = self._pl_pa_write_trigger()
        if wd is None:
            return
        self._pl_log(f"PA: params written to {wd / 'pa_trigger.ini'} (for '{action}')")
        self._pl_pa_status.configure(
            text=f"Params reloaded to pa_trigger.ini. '{action}' will use them on Run.", fg=YELLOW)

    def _pl_regen_triggers(self, honor_locate=True):
        """Rewrite af_trigger_*.ini from Step 3's currently-checked selection before a
        Step 4 Run, so every action uses the current coords (incl. re-centered) and a
        fresh, de-staled trigger set -- never whatever happened to be left on disk from
        a prior rescan / archived set / capture-consumed cycle. Delegates to
        trigger_autofocus_all, which also clears stale af_trigger/af_done and re-syncs
        session.ini's work_dir/nd2_dir to this run. Returns True to proceed, False to
        abort (status shown). No records (hand-placed triggers, Steps 1-3 skipped) ->
        True: leave the on-disk triggers untouched.
        """
        import spheroid_pipeline as _pl
        if not self._pl_records:
            return True
        out = self._pl_out_dir.get().strip()
        if not out:
            # No output dir -> nothing to regenerate into; the action falls back to
            # whatever triggers are already on disk. Say so rather than pretending
            # the regen happened (this is the exact staleness this helper exists to kill).
            self._pl_log("Step 4: no Output dir set -- skipped trigger regen, "
                         "using existing on-disk af_trigger_*.ini")
            return True
        try:
            z_centre = float(self._pl_z_centre.get())
            z_half   = float(self._pl_z_half.get())
            z_step   = float(self._pl_z_step.get())
        except ValueError:
            self.after(0, lambda: self._pl_pa_status.configure(
                text="Step 4: Step 3's Z centre/half/step must be numbers.", fg=RED))
            return False
        # Honor Step 3's "Locate only (1 plane)" exactly as _pl_trigger_autofocus does --
        # otherwise a Step 4 Run silently upgrades a 1-plane locate set into a full
        # Z-stack, i.e. the instrument does something the operator didn't ask for.
        #
        # EXCEPT for the Validation path (honor_locate=False): it has its OWN per-phase
        # "Locate only" checkbox and decides per pass via
        #   z_half = 0.0 if locate_only else d.get('z_half')
        # If we baked Step 3's zero onto disk here, that d.get() fallthrough would read
        # 0.0 back and the per-phase box could only ever force locate ON, never OFF --
        # so ticking Step 3's box would silently pin every Validation pass to 1 plane.
        # Write the full z_half instead and let the per-phase flag be authoritative.
        if (honor_locate
                and getattr(self, "_pl_locate_only", None) is not None
                and self._pl_locate_only.get()):
            z_half = 0.0
        excl = getattr(self, "_pl_excluded", set())
        recs = [r for r in self._pl_records if r.rank not in excl]
        if not recs:
            self.after(0, lambda: self._pl_pa_status.configure(
                text="Step 4: every spheroid is unchecked in Step 3's Use column — nothing to run.",
                fg=RED))
            return False
        try:
            n = _pl.trigger_autofocus_all(recs, Path(out) / "autofocus",
                                          z_centre=z_centre, z_half=z_half, z_step=z_step)
        except Exception as e:
            # MUST be broad: trigger_autofocus_all mkdirs, writes N trigger files and
            # session.ini. An OSError/PermissionError (S: share drop, read-only run
            # folder, C:/SpheroidPA ACL) previously escaped -- and on the worker-thread
            # call sites it died silently in the daemon thread with no dialog, after the
            # stale-trigger unlink had ALREADY run, leaving a partial trigger set the
            # next dispatch would happily image.
            self.after(0, lambda m=str(e): self._pl_pa_status.configure(
                text=f"Step 4: trigger regen FAILED — {m}", fg=RED))
            self._pl_log(f"Step 4: trigger regen FAILED ({type(e).__name__}): {e}")
            return False
        skipped = len(self._pl_records) - len(recs)
        self._pl_log(f"Step 4: regenerated {n} af_trigger from Step 3 selection"
                     + (" [locate-only, z_half=0]" if z_half == 0.0 else "")
                     + (f" ({skipped} unchecked skipped)" if skipped else ""))
        return True

    # ── Global macro abort ────────────────────────────────────────────────────
    def _pl_abort_path(self, work_dir=None):
        """Path of the one-shot abort flag the running macros poll."""
        import spheroid_pipeline as _pl
        if work_dir is None:
            import configparser
            cp = configparser.ConfigParser()
            try:
                if _pl.SESSION_INI.exists():
                    cp.read(_pl.SESSION_INI)
            except Exception:
                pass
            w = cp.get("paths", "work_dir", fallback="").strip() or self._pl_trigger_dir.get().strip()
            if not w:
                o = self._pl_out_dir.get().strip()
                w = str(Path(o) / "autofocus") if o else ""
            if not w:
                return None
            work_dir = Path(w)
        return Path(work_dir) / "abort.ini"

    def _pl_abort_requested(self, work_dir=None):
        """True if an abort flag is present (checked by the GUI-side capture loops)."""
        p = self._pl_abort_path(work_dir)
        try:
            return bool(p and p.exists())
        except Exception:
            return False

    def _pl_abort_clear(self, work_dir=None):
        p = self._pl_abort_path(work_dir)
        try:
            if p and p.exists():
                p.unlink()
        except Exception:
            pass

    def _pl_abort_all(self):
        """Abort All: drop a one-shot abort.ini that every long-running macro polls,
        and release the GUI-side dispatcher lock so the UI isn't stuck behind a macro
        that may never report back.

        Note: a macro already inside a blocking NIS-E call (mid ND Z-series, mid
        Capture) can only honor this at its NEXT loop iteration -- a true immediate
        hard-stop still requires Ctrl+Pause in NIS-E. Say so plainly rather than implying
        the button kills the scanner mid-frame."""
        import spheroid_pipeline as _pl
        if not getattr(self, "_pl_dispatch_busy", False):
            # Nothing is running -- planting the flag now would just be a landmine the
            # NEXT dispatch trips over. Clear any leftover instead.
            self._pl_abort_clear()
            self._pl_dispatch_status.configure(text="Dispatcher: idle — nothing to abort "
                                                    "(cleared any stale abort flag).", fg=SUBTEXT)
            messagebox.showinfo("Abort", "No macro is currently running via the dispatcher.\n\n"
                                         "Any stale abort flag has been cleared.")
            return
        p = self._pl_abort_path()
        if p is None:
            messagebox.showwarning("Abort", "No work dir known yet — run Step 3 (Trigger) first.")
            return
        try:
            p.parent.mkdir(parents=True, exist_ok=True)
            _pl._atomic_write_crlf(p, ["[abort]", "requested=1"])
        except Exception as e:
            messagebox.showerror("Abort", f"Could not write the abort flag:\n{e}")
            return
        self._pl_dispatch_busy = False
        self._pl_log(f"ABORT requested -- wrote {p}; GUI dispatcher lock released")
        self._pl_dispatch_status.configure(text="ABORT requested — macros stop at their next loop step.",
                                           fg=RED)
        messagebox.showinfo(
            "Abort requested",
            "Abort flag written.\n\nRunning macros stop at their next loop step (per-spheroid / "
            "per-pass boundary).\n\nA macro already inside a Z-series or Capture cannot be "
            "interrupted from here — press Ctrl+Pause in NIS-E for an immediate hard stop.")

    def _pl_pa_run_macro(self, action, action_id):
        if getattr(self, "_pl_dispatch_busy", False):
            messagebox.showinfo("Dispatcher busy",
                                "A macro is already running via the dispatcher — wait for it to finish."); return
        if not self._pl_regen_triggers():
            return
        if self._pl_pa_write_trigger() is None:
            return
        hint = {
            "pa_points": "PA Points: sent to NIS-E -- building ND multipoint from the active "
                         "triggers. Watch NIS-E for the 'Built ND multipoint with N spheroid(s)' popup.",
            "pa_setup":  "PA Setup: sent to NIS-E -- prepping the rig (interlock / OC / dichroic / "
                         "zoom). Watch NIS-E for the 'PA setup done' popup.",
        }.get(action, f"'{action}': sent to NIS-E via the dispatcher.")
        self._pl_pa_status.configure(text=hint, fg=YELLOW)
        self._pl_log(f"PA: dispatched '{action}' (id {action_id}) -- awaiting NIS-E")
        self._pl_send_command(action, action_id)

    def _pl_prepa_checked_ocs(self, phase="prePA"):
        """Return [(oc_name, tag), ...] for the checked Validation checkboxes of the
        given phase ('prePA' before PA / 'postPA' after), in a fixed 890->940->1050 order."""
        v890, v940, v1050 = (
            (self._pl_prepa_oc_890, self._pl_prepa_oc_940, self._pl_prepa_oc_1050)
            if phase == "prePA" else
            (self._pl_postpa_oc_890, self._pl_postpa_oc_940, self._pl_postpa_oc_1050))
        ocs = []
        if v890.get():
            ocs.append(("890nm_Galvo_600nm_NDD2_BT", "890nm_mBeRFP"))
        if v940.get():
            ocs.append(("940nm_Galvo_488nm_NDD2_JL", "940nm_PAsfGFP"))
        if v1050.get():
            ocs.append(("1050nm_Galvo_561nm_NDD2_WC", "1050nm_depth"))
        return ocs

    def _pl_archive_existing_nd2(self, nd2_dir, spheroid_ids, label="Capture"):
        """Move any capture this pass is about to overwrite into <nd2_dir>/_archive/.

        nis_macro_capture_zstack.mac saves to <nd2_dir>/<spheroid_id>_zstack.nd2 with
        ImageSaveAs, which overwrites in place, without warning, and with no undo. On
        2026-08-08 a re-run against the same Output directory with 'Locate only' ticked
        replaced three 61-plane baselines (32.6 MB each) with 1-plane captures (884 KB).
        The share (\\\\10.21.16.98) has no snapshots, so that raw data is simply gone --
        and it was the baseline half of the only complete before/after set of the night.

        Archiving rather than versioning in place is deliberate. Every analysis script
        does glob('*.nd2')[0] on these folders, and glob does not sort -- a sibling copy
        would make which file gets analysed depend on OS directory order. Moving the old
        one under _archive/ keeps exactly one canonical capture per folder, so nothing
        downstream has to change.

        The archived name carries the ORIGINAL capture's mtime, not the time of the
        overwrite: that is the timestamp identifying which run it belongs to.

        Never raises. Failing to archive must not abort a capture the operator is
        waiting on -- it degrades to the old overwrite behaviour with a loud log line.
        """
        try:
            d = Path(nd2_dir)
            if not d.is_dir():
                return
            wanted = {str(s).strip() for s in spheroid_ids if str(s).strip()}
            if not wanted:
                return
            arch = d / "_archive"
            moved = []
            for f in d.glob("*_zstack.nd2"):
                if f.stem[:-len("_zstack")] not in wanted:
                    continue                      # a different spheroid; this pass won't touch it
                try:
                    ts = time.strftime("%Y%m%d_%H%M%S", time.localtime(f.stat().st_mtime))
                    mb = f.stat().st_size / 1e6
                except OSError:
                    ts, mb = "unknown", 0.0
                arch.mkdir(parents=True, exist_ok=True)
                dest = arch / f"{f.stem}__{ts}.nd2"
                n = 1
                while dest.exists():              # same second, or a repeat archive
                    dest = arch / f"{f.stem}__{ts}_{n}.nd2"
                    n += 1
                f.rename(dest)                    # same filesystem -> atomic, no copy cost
                moved.append(f"{dest.name} ({mb:.1f} MB)")
            if moved:
                self._pl_log(f"{label}: archived {len(moved)} existing capture(s) rather than "
                             f"overwriting -> {d.name}/_archive/  [{', '.join(moved)}]")
        except Exception as exc:
            self._pl_log(f"{label}: WARNING -- could not archive existing captures in "
                         f"{nd2_dir} ({exc}). The capture will OVERWRITE them.")

    def _pl_prepa_capture_ocs(self, phase="prePA"):
        """Capture a full Z-stack at each checked Pre-PA OC, via the dispatcher's
        z-stack action (id 2) -- same mechanism as every other dispatcher card,
        just looped once per checked OC.

        First regenerates af_trigger_NN.ini from Step 3's currently-checked selection
        (via _pl_regen_triggers -> trigger_autofocus_all: clears stale triggers, writes
        current/re-centered coords, re-syncs session.ini) so this pass never images a
        stale position. Then injects the oc= key into each trigger (and optionally
        forces z_half=0 for a locate-only pass). With no records (hand-placed triggers,
        Steps 1-3 skipped) it leaves whatever is on disk untouched.

        Each pass is routed into its own nd2/prePA_<tag> subfolder (atomic session.ini
        rewrite before each dispatch) so multiple OC passes don't overwrite each
        other or the plain Step-3 capture. Must be called with _pl_dispatch_busy
        already held by the caller (mirrors _pl_pa_run_pipeline's locking contract).
        Returns True if every checked OC completed 'ok' (or none were checked --
        a no-op is not an error); False on the first failure (status already shown).
        """
        import configparser, time, spheroid_pipeline as _pl
        plabel = "Pre-PA" if phase == "prePA" else "After-PA"
        ocs = self._pl_prepa_checked_ocs(phase)
        if not ocs:
            return True
        # honor_locate=False: this path has its own per-phase "Locate only" checkbox and
        # applies it per pass below -- see the note in _pl_regen_triggers.
        if not self._pl_regen_triggers(honor_locate=False):
            return False

        # Resolve work_dir/nd2_dir the SAME way every other dispatcher action does
        # (session.ini's [paths] first -- it reflects whatever Step 3/the dispatcher
        # last set, e.g. .../work/0706/autofocus -- THEN _pl_trigger_dir, THEN
        # _pl_out_dir/autofocus as a last resort). Using _pl_out_dir alone (Step 1's
        # Output-dir field, which defaults to .../work with no run subfolder) looked
        # in the wrong folder and always found zero triggers.
        cp = configparser.ConfigParser()
        try:
            if _pl.SESSION_INI.exists():
                cp.read(_pl.SESSION_INI)
        except Exception:
            pass
        work = cp.get("paths", "work_dir", fallback="").strip() or self._pl_trigger_dir.get().strip()
        if not work:
            o = self._pl_out_dir.get().strip()
            work = str(Path(o) / "autofocus") if o else ""
        if not work:
            self.after(0, lambda: self._pl_pa_status.configure(
                text=f"{plabel}: no work dir -- run Step 3 (Trigger) first, or set the Trigger dir.",
                fg=RED)); return False
        wd = Path(work)
        nd2_dir_cfg = cp.get("paths", "nd2_dir", fallback="").strip()
        base_nd2_dir = Path(nd2_dir_cfg) if nd2_dir_cfg else wd.parent / "nd2"
        trig_files = sorted(wd.glob(f"{_pl.AF_TRIGGER_PREFIX}*.ini"))
        triggers = []
        for p in trig_files:
            d = _pl.parse_ini(p.read_text(encoding="utf-8"))
            if {"rank", "stage_x", "stage_y"} <= d.keys():
                triggers.append(d)
        if not triggers:
            self.after(0, lambda: self._pl_pa_status.configure(
                text=f"{plabel}: no af_trigger_*.ini found -- generate them in Step 3 first.",
                fg=RED)); return False
        locate_only = (self._pl_prepa_locate_only if phase == "prePA"
                       else self._pl_postpa_locate_only).get()
        # close_after=1 -> the capture macro closes each ND2 right after ImageSaveAs, so a
        # long multi-spheroid x multi-OC run doesn't leave dozens of windows open in NIS-E
        # (which clutters it and can stall the acquisition).
        close_after = "0" if self._pl_keep_open_in_nise.get() else "1"
        captured: list[tuple[str, str]] = []   # (display name, path) -> auto-loaded into "Captured Z-Stacks" on success

        def _point_nd2_dir(d):
            _pl._atomic_write_crlf(_pl.SESSION_INI, [
                "[paths]", f"work_dir={wd.as_posix()}",
                f"nd2_dir={d.as_posix()}", f"macro_dir={_pl.MACRO_DIR.as_posix()}"])

        # try/finally: session.ini's nd2_dir is redirected into the per-phase subfolder
        # for each pass. If ANYTHING throws mid-loop (share drop, mkdir PermissionError,
        # bad rank) the daemon thread dies and nd2_dir would stay pointed at
        # nd2/<phase>_<tag>/ -- silently poisoning every later capture. Always restore.
        try:
            # Light path FIRST, before any capture: snapshot on the pre-PA pass so the
            # baseline laser powers are read off the rig; restore on the post-PA pass so
            # the dichroic/zoom pa_setup changed are put back. Without this the two halves
            # of the comparison can be acquired through different optics -- see
            # _pl_lightpath_step.
            if self._pl_lightpath_restore.get():
                self._pl_lightpath_step(wd, "snapshot" if phase == "prePA" else "restore",
                                        plabel)
            for oc_name, tag in ocs:
                if self._pl_abort_requested(wd):
                    self._pl_log(f"{plabel}: ABORT requested -- stopping before '{tag}'")
                    self.after(0, lambda: self._pl_pa_status.configure(
                        text=f"{plabel}: aborted by operator.", fg=RED))
                    return False
                self._pl_log(f"{plabel}: capturing {tag} ({oc_name}) for {len(triggers)} spheroid(s)...")
                self.after(0, lambda t=tag: self._pl_pa_status.configure(
                    text=f"{plabel}: running {t}...", fg=YELLOW))
                for d in triggers:
                    lines = ["[spheroid]", f"rank={d['rank']}",
                             f"spheroid_id={d.get('spheroid_id', '')}",
                             f"stage_x={d['stage_x']}", f"stage_y={d['stage_y']}"]
                    if "z_centre" in d:
                        lines.append(f"z_centre={d['z_centre']}")
                    lines.append(f"z_half={'0.0' if locate_only else d.get('z_half', '0.0')}")
                    if "z_step" in d:
                        lines.append(f"z_step={d['z_step']}")
                    lines.append(f"oc={oc_name}")
                    lines.append(f"close_after={close_after}")
                    _pl._atomic_write_crlf(wd / f"af_trigger_{int(d['rank']):02d}.ini", lines)
                nd2_sub = base_nd2_dir / f"{phase}_{tag}"
                nd2_sub.mkdir(parents=True, exist_ok=True)
                # Never overwrite a capture that is already there -- archive it first.
                # This folder is exactly where the 2026-08-08 baseline loss happened.
                self._pl_archive_existing_nd2(
                    nd2_sub, [d.get("spheroid_id", "") for d in triggers],
                    f"{plabel} {tag}")
                _point_nd2_dir(nd2_sub)

                done = wd / "cmd_done.ini"
                try:
                    done.unlink()
                except FileNotFoundError:
                    pass
                _pl._atomic_write_crlf(wd / "cmd.ini",
                                       ["[command]", f"action=zstack_{tag}", "action_id=2"])
                self._pl_log(f"{plabel}: dispatched '{tag}' (id 2)")
                status = None
                for _ in range(1800):
                    if done.exists():
                        cp = configparser.ConfigParser()
                        try:
                            cp.read(done)
                            status = cp.get("command", "status",
                                            fallback=cp.get("spheroid", "status", fallback="?"))
                        except Exception:
                            time.sleep(1.0); continue
                        if status == "?":
                            # cmd_done.ini is built key-by-key with status written LAST, so a
                            # file that parses without one is still being written -- not a
                            # finished run with an unknown result. Breaking here reported "?"
                            # and aborted a post-PA sequence whose capture had actually
                            # succeeded (2026-08-07, cost the 1050 nm pass).
                            time.sleep(1.0); continue
                        break
                    if self._pl_abort_requested(wd):
                        status = "aborted"
                        break
                    time.sleep(1.0)
                self._pl_log(f"{plabel}: '{tag}' -> {status or 'timeout'}")
                if status != "ok":
                    msg = (f"{plabel} stopped at '{tag}' -> {status}" if status
                           else f"{plabel}: '{tag}' timed out (is nis_macro_dispatcher.mac running?)")
                    self.after(0, lambda m=msg: self._pl_pa_status.configure(text=m, fg=RED))
                    return False
                for d in triggers:
                    sid = (d.get("spheroid_id") or "").strip() or f"rank{int(d['rank']):02d}"
                    p = nd2_sub / f"{sid}_zstack.nd2"
                    # phase-qualify the display key: tags are identical across phases, so
                    # without the prefix a postPA entry silently overwrites its prePA twin.
                    captured.append((f"{phase}/{tag}/{sid}", str(p)))
        finally:
            # Guarded: this writes session.ini, the very operation that can fail on a
            # share drop. Unguarded it would (a) skip the trigger restore below --
            # reinstating the PA-Points-strands-on-empty-folder bug -- and (b) replace
            # any in-flight exception with its own, hiding the real cause.
            try:
                _point_nd2_dir(base_nd2_dir)   # restore the plain nd2/ dir for normal captures
            except Exception as e:
                self._pl_log(f"{plabel}: WARNING - could not restore session.ini nd2_dir: {e}")
            # Restore the plain (oc-less) triggers PA Points depends on. Each OC pass
            # above rewrote af_trigger_NN.ini with oc= set, and nis_macro_capture_zstack.mac
            # DELETES it after every successful capture -- so by the time we leave this
            # function, every trigger touched here is gone. Run Pipeline's next step is
            # PA Points, which only ever reads whatever's currently in autofocus/: without
            # this restore it finds an EMPTY folder and silently builds a 0-point
            # multipoint. Re-persist using the SAME coordinates already used above -- no
            # staleness introduced, just putting back what was already valid.
            #
            # MUST be in the finally, not after it: a mid-run bail (abort, non-ok status,
            # or any exception) leaves triggers already consumed by the macro, so the
            # non-restoring paths were exactly the ones that stranded PA Points.
            try:
                for d in triggers:
                    lines = ["[spheroid]", f"rank={d['rank']}",
                             f"spheroid_id={d.get('spheroid_id', '')}",
                             f"stage_x={d['stage_x']}", f"stage_y={d['stage_y']}"]
                    if "z_centre" in d:
                        lines.append(f"z_centre={d['z_centre']}")
                    lines.append(f"z_half={'0.0' if locate_only else d.get('z_half', '0.0')}")
                    if "z_step" in d:
                        lines.append(f"z_step={d['z_step']}")
                    _pl._atomic_write_crlf(wd / f"af_trigger_{int(d['rank']):02d}.ini", lines)
            except Exception as e:
                self._pl_log(f"{plabel}: WARNING - could not restore triggers for PA Points: {e}")

        self.after(0, lambda n=len(ocs): self._pl_pa_status.configure(
            text=f"{plabel}: {n} OC pass(es) complete.", fg=GREEN))
        if captured:
            self.after(0, lambda items=captured, lb=plabel: self._pl_zv_load_captured(items, lb))
        return True

    def _pl_prepa_run_thread(self, phase="prePA"):
        if getattr(self, "_pl_dispatch_busy", False):
            messagebox.showinfo("Dispatcher busy",
                                "A macro is already running via the dispatcher — wait for it to finish."); return
        if not self._pl_prepa_checked_ocs(phase):
            lbl = "Pre-PA" if phase == "prePA" else "After-PA"
            messagebox.showwarning("No OC checked", f"Tick at least one {lbl} OC checkbox."); return
        self._pl_abort_clear()       # a stale flag must not kill the run we're starting
        self._pl_dispatch_busy = True
        threading.Thread(target=lambda: self._pl_prepa_run(phase), daemon=True).start()

    # Zoom used for the validation captures. The pre-PA baseline is acquired at the Init
    # card's confocal zoom; the post-PA restore must reproduce THAT, not whatever zoom
    # pa_setup left behind for the PA square.
    LIGHTPATH_ZOOM_FALLBACK = 2.0

    def _pl_lightpath_step(self, wd, mode, plabel):
        """Snapshot (pre-PA) or restore (post-PA) the imaging light path via dispatcher
        action 12. Returns True to continue.

        Restore is what makes a before/after comparison mean anything: pa_setup extracts
        the dichroic and drives the zoom to the PA square, and the capture path (action 2)
        never puts either back. Snapshot runs on the pre-PA pass so the laser powers come
        off the rig rather than being assumed.

        A failure here is NOT fatal to the capture -- it is reported loudly and the caller
        proceeds, because aborting a post-PA run (which cannot be repeated on that sample)
        over a restore hiccup is worse than capturing with a warning attached."""
        import configparser, time, spheroid_pipeline as _pl
        try:
            zoom = float(self._pl_init_zoom.get())
        except (ValueError, AttributeError):
            zoom = self.LIGHTPATH_ZOOM_FALLBACK
        state = (Path(wd) / "lightpath_state.ini").as_posix()
        lines = ["[lightpath]", f"mode={mode}",
                 f"dichroic_in={'1' if mode == 'restore' else '0'}",
                 f"zoom={zoom if mode == 'restore' else 0.0}",
                 f"a1_on={'1' if self._pl_pa_a1on.get() else '0'}",
                 f"state={state}"]
        try:
            _pl._atomic_write_crlf(Path(wd) / "lightpath_trigger.ini", lines)
        except Exception as e:
            self._pl_log(f"{plabel}: light-path {mode} SKIPPED (could not write trigger: {e})")
            return True
        done = Path(wd) / "cmd_done.ini"
        try:
            done.unlink()
        except FileNotFoundError:
            pass
        _pl._atomic_write_crlf(Path(wd) / "cmd.ini",
                               ["[command]", f"action=lightpath_{mode}", "action_id=12"])
        if mode == "restore":
            self._pl_log(f"{plabel}: restoring pre-PA light path (dichroic IN, zoom {zoom:g}, "
                         f"laser powers from {Path(state).name}) before capture")
        else:
            self._pl_log(f"{plabel}: snapshotting laser powers -> {Path(state).name}")
        status = None
        for _ in range(180):                       # 3 min; this is a fast macro
            if done.exists():
                cp = configparser.ConfigParser()
                try:
                    cp.read(done)
                    status = cp.get("command", "status",
                                    fallback=cp.get("spheroid", "status", fallback="?"))
                except Exception:
                    time.sleep(1.0); continue
                if status == "?":
                    # cmd_done.ini is built key-by-key with status written LAST, so a
                    # file that parses without one is still being written -- not a
                    # finished run with an unknown result. Breaking here reported "?"
                    # and aborted a post-PA sequence whose capture had actually
                    # succeeded (2026-08-07, cost the 1050 nm pass).
                    time.sleep(1.0); continue
                break
            time.sleep(1.0)
        if status != "ok":
            # Loud, and surfaced in the UI: a post-PA set captured without the restore is
            # not comparable to its baseline, and that is invisible in the images.
            self._pl_log(f"{plabel}: light-path {mode} -> {status or 'timeout'} "
                         f"-- CHECK the dichroic/zoom by hand before trusting this data.")
            self.after(0, lambda: self._pl_pa_status.configure(
                text=f"{plabel}: light-path {mode} FAILED — verify dichroic IN and zoom "
                     f"{zoom:g} manually.", fg=RED))
        return True

    def _pl_prepa_run(self, phase="prePA"):
        try:
            self._pl_prepa_capture_ocs(phase)
        finally:
            self._pl_dispatch_busy = False

    # Waiting on the PA JOB. The job is asynchronous: cmd_done returns when the LAUNCH
    # MACRO returns, not when the job finishes, and the dispatcher writes status=ok
    # unconditionally -- so neither proves the laser fired. New nd2 files appearing is the
    # only real evidence, and the run is over when they stop appearing.
    # Only files the PA JOB itself produces count as evidence it fired. The job names its
    # output "<timestamp>__Point0000_ZStack00NN_Channel...nd2"; anything else in that folder
    # is not PA output. This matters because NIS-E's ND "Save to File" has repeatedly dumped
    # copies of the VALIDATION captures into the same folder (PA.nd2, PA001.nd2, ...), and
    # counting those made the gate below report success for a laser that never fired --
    # which would then wave the run on to post-PA on an unactivated spheroid
    # (2026-08-07: three such duplicates sat in work/0807_4/pa with no PA run at all).
    PA_OUTPUT_GLOB = "*ZStack*.nd2"
    PA_JOB_LAUNCH_TIMEOUT_S = 300     # no output at all by now -> the job never started
    PA_JOB_QUIET_S          = 90      # output stopped growing this long -> finished
    PA_JOB_MAX_S            = 5400    # hard cap so the pipeline can never hang forever

    def _pl_pa_output_dirs(self):
        """Folders a PA run might write into. The JOB's destination lives on its own
        'Alternative Storage Location' task, which this GUI cannot read -- so watch the
        folder the operator configured AND the shared work/PA that the job has actually
        been writing to. Whichever fills is reported."""
        dirs = []
        d = self._pl_pa_effective_save_dir()
        if d:
            dirs.append(Path(d))
        try:
            base = Path(self._pl_out_dir.get().strip()).parent / "PA"
            if base not in dirs:
                dirs.append(base)
        except Exception:
            pass
        return [x for x in dirs if x]

    def _pl_pa_wait_for_job(self, before, work_dir=None):
        """Block until the PA job has produced output and gone quiet.

        Returns (ok, message). ok=False means NO output ever appeared -- the job did not
        start -- and the caller must NOT go on to post-PA validation: a post-PA set with no
        PA in between is not a control, it is a mislabelled duplicate of the baseline.

        work_dir is where abort.ini lives. It must be the WORK DIR: the first version
        derived it from `before`, which is a set of nd2 FILE paths, so the abort check
        looked for the flag next to an image and could never fire (found 2026-08-07 on the
        first real pipeline run -- the operator's Abort All was ignored by this loop)."""
        import time
        dirs = self._pl_pa_output_dirs()
        t0 = time.time()
        seen = 0
        last_change = None
        while time.time() - t0 < self.PA_JOB_MAX_S:
            if self._pl_abort_requested(work_dir):
                return False, "aborted by operator"
            now = 0
            for d in dirs:
                try:
                    now += sum(1 for p in d.glob(self.PA_OUTPUT_GLOB)
                               if str(p) not in before)
                except Exception:
                    pass
            if now != seen:
                seen, last_change = now, time.time()
                self._pl_log(f"PA job: {seen} output file(s) so far")
            if seen == 0:
                if time.time() - t0 > self.PA_JOB_LAUNCH_TIMEOUT_S:
                    return False, (f"no PA output after "
                                   f"{self.PA_JOB_LAUNCH_TIMEOUT_S // 60} min in "
                                   f"{', '.join(d.as_posix() for d in dirs)} -- the job did "
                                   f"not start (check its Alternative Storage Location, and "
                                   f"that it validates)")
            elif last_change and time.time() - last_change > self.PA_JOB_QUIET_S:
                return True, f"{seen} file(s), output quiet for {self.PA_JOB_QUIET_S}s"
            time.sleep(5)
        return False, f"PA job still writing after {self.PA_JOB_MAX_S // 60} min -- gave up waiting"

    # Preference order for which pre-PA capture to re-centre from: the structural /
    # brightest channel first. 890 mBeRFP is routinely near the noise floor on these
    # samples (p99.9 of 62-79 counts), so segmenting it is the least reliable option and
    # is only used when nothing better was captured.
    RECENTER_TAG_ORDER = ["1050nm_depth", "940nm_PAsfGFP", "890nm_mBeRFP"]

    def _pl_pipeline_recenter(self):
        """Re-centre the triggers off the pre-PA captures this pipeline just took, then
        rewrite them so the PA fires on the corrected coordinates.

        Returns True to continue regardless of outcome: a failed re-centre leaves the
        pre-existing (anchored) coordinates in place, which is exactly what a run without
        re-centring would have used anyway. Aborting the sequence over it would be a
        bigger loss than proceeding slightly off-centre."""
        import configparser, spheroid_pipeline as _pl
        if not self._pl_records:
            return True
        cp = configparser.ConfigParser()
        try:
            if _pl.SESSION_INI.exists():
                cp.read(_pl.SESSION_INI)
        except Exception:
            pass
        base = cp.get("paths", "nd2_dir", fallback="").strip()
        if not base:
            self._pl_log("Pipeline re-center: no nd2_dir -- SKIPPED (triggers unchanged).")
            return True
        # A single-plane ("Locate only") capture is NOT a valid re-centre input: the
        # centroid is measured on whatever one plane happened to be taken, so the nudge
        # comes out near zero and looks like "already centred". 2026-08-07: a locate-only
        # run produced +0.3/-0.8 um on a spheroid the operator's manual pass had corrected
        # by 36.4 um. Refuse loudly rather than return a meaningless correction.
        MIN_PLANES = 5
        src, thin = None, []
        for tag in self.RECENTER_TAG_ORDER:
            d = Path(base) / f"prePA_{tag}"
            if not (d.is_dir() and any(d.glob("*_zstack.nd2"))):
                continue
            try:
                import nd2 as _nd2
                with _nd2.ND2File(str(next(d.glob("*_zstack.nd2")))) as f:
                    nz = f.sizes.get("Z", 1)
            except Exception:
                nz = 0
            if nz >= MIN_PLANES:
                src = d
                break
            thin.append(f"{tag}({nz}p)")
        if src is None:
            why = (f"only single/thin stacks {thin} -- 'Locate only' is on, or z_half=0"
                   if thin else "no prePA_* captures found")
            self._pl_log(f"Pipeline re-center: SKIPPED -- {why}. Triggers UNCHANGED; the PA "
                         f"will fire on the un-recentred coordinates.")
            self.after(0, lambda w=why: self._pl_pa_status.configure(
                text=f"Re-center SKIPPED — {w}", fg=RED))
            return True
        excl = getattr(self, "_pl_excluded", set())
        targets = sorted((r for r in self._pl_records if r.rank not in excl),
                         key=lambda r: r.rank)
        try:
            res, skipped = _pl.recenter_from_captures(
                targets, src,
                flip_x=self._pl_recenter_flipx.get(),
                flip_y=self._pl_recenter_flipy.get())
        except Exception as e:
            self._pl_log(f"Pipeline re-center: FAILED ({e}) -- triggers unchanged, "
                         f"continuing on the anchored coordinates.")
            return True
        if not res and not skipped:
            self._pl_log(f"Pipeline re-center: nothing measurable in {src.name} -- "
                         f"triggers unchanged.")
            return True
        self._pl_recentered = {rank: (dx, dy) for rank, dcol, drow, dx, dy in res}
        self._pl_log(f"Pipeline re-center: from {src.name}")
        for rank, dcol, drow, dx, dy in res:
            self._pl_log(f"   rank {rank}: nudge ({dx:+.1f}, {dy:+.1f}) um")
        # Skipped ranks keep their PRE-recenter coords -- say so, since silently leaving
        # them uncorrected is the failure this whole step exists to prevent.
        for rank, reason in skipped:
            self._pl_log(f"   rank {rank}: NOT corrected -- {reason}")
        self.after(0, self._pl_update_table)
        self.after(0, self._pl_refresh_dashboard)
        self._pl_regen_triggers()      # write the corrected coords for the PA launch
        return True

    def _pl_pa_run_pipeline_thread(self):
        if getattr(self, "_pl_dispatch_busy", False):
            messagebox.showinfo("Dispatcher busy",
                                "A macro is already running via the dispatcher — wait for it to finish."); return
        # Gate the laser HERE, on the main thread, before the worker starts. Tk modals are
        # not thread-safe, and an unattended sequence must not stop halfway waiting for a
        # dialog nobody is watching. Once confirmed, the run proceeds without prompting.
        if self._pl_pa_sel_pajob.get():
            if not self._pl_pa_a1on.get():
                messagebox.showwarning("A1 not confirmed",
                                       "Tick 'A1 powered ON' before running a pipeline that "
                                       "fires the PA laser."); return
            pts, _ = self._pl_pa_points_from_triggers()
            n = len(pts) if pts else 0
            oc = self._pl_pa_oc.get().strip()
            listing = "\n".join(f"   [{i}] ({x:.1f}, {y:.1f}) um  z={z:.1f}"
                                for i, (_n, x, y, z) in enumerate(pts or []))
            if not messagebox.askyesno(
                    "Run Pipeline will FIRE THE PA LASER",
                    f"This pipeline includes the photoactivation JOB.\n\n"
                    f"Optical config : {oc}\n"
                    f"Points ({n}):\n{listing or '   (none — job would use its STORED point)'}\n\n"
                    f"Post-PA captures will run only if the PA actually produces output.\n\n"
                    f"Proceed?"):
                self._pl_log("Pipeline: cancelled at the PA-laser confirmation.")
                return
            self._pl_log(f"Pipeline: PA laser confirmed by operator (OC '{oc}', {n} point(s)).")
        self._pl_abort_clear()       # a stale flag must not kill the run we're starting
        self._pl_dispatch_busy = True
        threading.Thread(target=self._pl_pa_run_pipeline, daemon=True).start()

    def _pl_pipeline_dispatch(self, wd, action, aid, timeout_s=1800):
        """Dispatch one dispatcher action from inside Run Pipeline and wait for its
        cmd_done.ini. Returns the status string ('ok', 'aborted', an error, or None on
        timeout).

        Dispatches DIRECTLY rather than through _pl_send_command: that path refuses when
        _pl_dispatch_busy is set, and the pipeline holds that lock for its whole sequence.
        Routing a pipeline step through it produces a modal on a worker thread and a step
        that never runs while the log still says it was sent (2026-08-07, three PA runs)."""
        import configparser, time, spheroid_pipeline as _pl
        if self._pl_abort_requested(wd):
            return "aborted"
        done = wd / "cmd_done.ini"
        try:
            done.unlink()
        except FileNotFoundError:
            pass
        _pl._atomic_write_crlf(wd / "cmd.ini",
                               ["[command]", f"action={action}", f"action_id={aid}"])
        self.after(0, lambda a=action: self._pl_pa_status.configure(
            text=f"Pipeline: running '{a}'...", fg=YELLOW))
        self._pl_log(f"PA pipeline: dispatched '{action}'")
        for _ in range(timeout_s):
            if done.exists():
                cp = configparser.ConfigParser()
                try:
                    cp.read(done)
                    status = cp.get("command", "status",
                                    fallback=cp.get("spheroid", "status", fallback="?"))
                except Exception:
                    time.sleep(1.0); continue
                if status == "?":
                    # cmd_done.ini is built key-by-key with status written LAST, so a file
                    # that parses without one is still being written -- not a finished run
                    # with an unknown result. Breaking here reported "?" and aborted a
                    # post-PA sequence whose capture had actually succeeded (2026-08-07,
                    # cost the 1050 nm pass).
                    time.sleep(1.0); continue
                return status
            if self._pl_abort_requested(wd):
                return "aborted"
            time.sleep(1.0)
        return None

    def _pl_pipeline_init_rig(self, wd, when):
        """Normalize the rig (init_rig, action 8) inside Run Pipeline. Called once before
        the pre-PA captures and again before the post-PA captures so both halves of the
        before/after set start from the same state -- PA Setup leaves the rig at dichroic
        OUT and the PA zoom, and nothing else puts it back. Returns True to continue.

        A failure here is NOT fatal to the run: the post-PA path also restores the light
        path explicitly (_pl_lightpath_step), so a failed init is a warning, not a reason
        to throw away captures that are about to succeed."""
        try:
            desc = self._pl_write_init_trigger(wd)
        except Exception as exc:
            self._pl_log(f"PA pipeline: init_rig SKIPPED before {when} -- "
                         f"could not write init_trigger.ini ({exc})")
            return True
        self._pl_log(f"PA pipeline: init_rig before {when} ({desc})")
        status = self._pl_pipeline_dispatch(wd, "init_rig", 8, timeout_s=300)
        if status == "aborted":
            self._pl_log("PA pipeline: ABORT requested during init_rig")
            return False
        if status != "ok":
            self._pl_log(f"PA pipeline: init_rig -> {status or 'timeout'} "
                         f"(continuing; the light-path restore still applies)")
        else:
            self._pl_log("PA pipeline: init_rig -> ok")
        return True

    def _pl_pa_run_pipeline(self):
        """Run the checked PA cards in order, each via the dispatcher, waiting for
        cmd_done.ini between steps. Stops on the first non-ok result. Holds the
        single-in-flight _pl_dispatch_busy lock for the whole sequence."""
        import configparser, time, spheroid_pipeline as _pl
        try:
            if not self._pl_regen_triggers():
                return
            # Declare the plan before executing it. Without this there is no way to tell,
            # from the log alone, whether a step was SKIPPED because it was unticked or
            # because something failed -- which is exactly the ambiguity that made the
            # 2026-08-07 run unreadable.
            self._pl_log(
                "PA pipeline: selected -> init_rig (always, before pre-PA and post-PA), "
                f"pre-PA={self._pl_pa_sel_prepa.get()} "
                f"(ocs={[t for _, t in self._pl_prepa_checked_ocs()]}), "
                f"recenter={self._pl_pa_sel_recenter.get()}, "
                f"setup={self._pl_pa_sel_setup.get()}, "
                f"PA-job={self._pl_pa_sel_pajob.get()}, "
                f"post-PA={self._pl_pa_sel_postpa.get()} "
                f"(ocs={[t for _, t in self._pl_prepa_checked_ocs('postPA')]})")
            # Z geometry and the locate-only flags, because THOSE decide whether the
            # baseline is a usable stack or a single plane -- and a single-plane baseline
            # makes the whole before/after run worthless while every step still reports ok.
            _lo_pre = self._pl_prepa_locate_only.get()
            _lo_post = self._pl_postpa_locate_only.get()
            self._pl_log(
                f"PA pipeline: Z -> centre={self._pl_z_centre.get()} "
                f"half={self._pl_z_half.get()} step={self._pl_z_step.get()}   "
                f"locate-only: pre-PA={_lo_pre} post-PA={_lo_post}")
            if _lo_pre or _lo_post or str(self._pl_z_half.get()).strip() in ("0", "0.0", ""):
                self._pl_log("PA pipeline: *** WARNING -- captures will be SINGLE-PLANE. "
                             "The before/after comparison needs a z-stack; untick "
                             "'Locate only' and set a non-zero z_half. ***")
                self.after(0, lambda: self._pl_pa_status.configure(
                    text="WARNING: locate-only / z_half=0 — captures will be single-plane",
                    fg=RED))
            do_prepa = self._pl_pa_sel_prepa.get() and bool(self._pl_prepa_checked_ocs())
            seq = []
            if self._pl_pa_sel_setup.get():    seq.append(("pa_setup", 4))
            # pa_points deliberately NOT sequenced -- the card is gone and the ND
            # multipoint it built is not what positions the PA. See the removal note in
            # the Step-4 card builder.
            if not seq and not do_prepa:
                self.after(0, lambda: self._pl_pa_status.configure(
                    text="Pipeline: no macros checked.", fg=RED)); return
            # The work dir is needed BEFORE the pre-PA captures now, because init_rig runs
            # first. Resolve it without a dialog (worker thread) and fall back to the
            # trigger writer's own resolution if it is not set up yet.
            wd_early = self._pl_resolve_work_dir()
            if do_prepa:
                if wd_early is not None and not self._pl_pipeline_init_rig(wd_early, "pre-PA"):
                    self.after(0, lambda: self._pl_pa_status.configure(
                        text="Pipeline: aborted during init.", fg=RED)); return
                self._pl_log("PA pipeline: Pre-PA viz capture(s) first")
                if not self._pl_prepa_capture_ocs():
                    return   # status already shown by the helper
                # Re-centre off those captures BEFORE anything downstream uses the
                # coordinates. The manual Step-3 Re-center is easy to forget, and the
                # pipeline has the exact same input sitting on disk by this point.
                if self._pl_pa_sel_recenter.get():
                    self.after(0, lambda: self._pl_pa_status.configure(
                        text="Pipeline: re-centering from the pre-PA captures...", fg=YELLOW))
                    self._pl_pipeline_recenter()
            if not seq and not (self._pl_pa_sel_pajob.get() or self._pl_pa_sel_postpa.get()):
                self.after(0, lambda: self._pl_pa_status.configure(
                    text="Pipeline complete: Pre-PA only.", fg=GREEN)); return
            wd = self._pl_pa_write_trigger()
            if wd is None:
                return
            self._pl_log(f"PA pipeline: {[a for a, _ in seq]}")
            for action, aid in seq:
                status = self._pl_pipeline_dispatch(wd, action, aid)
                if status == "aborted":
                    self._pl_log(f"PA pipeline: ABORT requested -- stopping at '{action}'")
                    self.after(0, lambda: self._pl_pa_status.configure(
                        text="Pipeline: aborted by operator.", fg=RED)); return
                self._pl_log(f"PA pipeline: '{action}' -> {status or 'timeout'}")
                if status != "ok":
                    msg = (f"Pipeline stopped at '{action}' -> {status}" if status
                           else f"Pipeline: '{action}' timed out (is nis_macro_dispatcher.mac running?)")
                    self.after(0, lambda m=msg: self._pl_pa_status.configure(text=m, fg=RED)); return
            # ── PA JOB (fires the 850 nm laser) ───────────────────────────────
            # Confirmed up front on the main thread in _pl_pa_run_pipeline_thread, so no
            # modal here: an unattended sequence must not stall on a dialog nobody sees.
            pa_ran = False
            if self._pl_pa_sel_pajob.get():
                if self._pl_abort_requested(wd):
                    self.after(0, lambda: self._pl_pa_status.configure(
                        text="Pipeline: aborted before the PA job.", fg=RED)); return
                before = set()
                for d in self._pl_pa_output_dirs():
                    try:
                        before |= {str(p) for p in d.glob("*.nd2")}
                    except Exception:
                        pass
                self.after(0, lambda: self._pl_pa_status.configure(
                    text="Pipeline: PA JOB running (850 nm)...", fg=RED))
                # Same automatic folder prep as the standalone button.
                self._pl_pa_sync_save_dir()
                self._pl_pa_prepare_save_dir()
                self._pl_log("PA pipeline: launching step3_zstack_PA_WC (laser)")
                # Dispatch DIRECTLY, exactly as this pipeline dispatches its macro steps.
                # It must NOT go through _pl_run_job -> _pl_send_command: that path's first
                # act is to refuse when _pl_dispatch_busy is set, and this pipeline holds
                # that lock for its whole sequence. The result was a "Dispatcher busy"
                # modal on a worker thread and a PA that never launched, while
                # _pl_run_job still logged "launched" because it logs after dispatching
                # regardless of outcome (2026-08-07: three runs lost to this).
                pa_params = self._pl_pa_point_params()
                if not self._pl_pa_a1on.get():
                    self._pl_log("PA pipeline: A1 not confirmed -- PA step SKIPPED.")
                    self.after(0, lambda: self._pl_pa_status.configure(
                        text="Pipeline: A1 not confirmed — PA skipped.", fg=RED)); return
                try:
                    _pl._atomic_write_crlf(wd / "job_trigger.ini", [
                        "[job]", f"project={self._pl_job_project.get().strip()}",
                        "name=step3_zstack_PA_WC"])
                    pj = wd / "job_params.json"
                    if pa_params:
                        import json as _json
                        pj.write_text(_json.dumps(pa_params, ensure_ascii=False),
                                      encoding="utf-16")
                    else:
                        pj.unlink(missing_ok=True)
                    (wd / "cmd_done.ini").unlink(missing_ok=True)
                    _pl._atomic_write_crlf(wd / "cmd.ini",
                                           ["[command]", "action=run_job", "action_id=10"])
                except Exception as exc:
                    self._pl_log(f"PA pipeline: could not write the PA launch files ({exc})")
                    self.after(0, lambda e=exc: self._pl_pa_status.configure(
                        text=f"Pipeline: PA launch write failed — {e}", fg=RED)); return
                self._pl_log(f"PA pipeline: dispatched run_job (id 10) for "
                             f"step3_zstack_PA_WC with {len(pa_params or {})} param(s)")
                ok, why = self._pl_pa_wait_for_job(before, wd)
                pa_ran = ok
                self._pl_log(f"PA pipeline: PA job -> {'complete' if ok else 'NOT RUN'} ({why})")
                if not ok:
                    # Hard stop. A post-PA set with no PA in between is not a control -- it
                    # is a mislabelled second copy of the baseline, and it silently corrupts
                    # the comparison. The operator explicitly asked for this behaviour.
                    self.after(0, lambda w=why: self._pl_pa_status.configure(
                        text=f"Pipeline STOPPED: PA did not run ({w}). Post-PA SKIPPED.",
                        fg=RED))
                    messagebox.showerror(
                        "PA did not run — post-PA skipped",
                        f"The photoactivation job produced no output:\n\n{why}\n\n"
                        f"After-PA captures were NOT run: without a PA in between they "
                        f"would just duplicate the baseline and corrupt the comparison.")
                    return

            # ── Post-PA validation (only reachable if the PA actually ran) ────
            if not self._pl_pa_sel_postpa.get():
                # Say so. Ending the sequence silently after a successful PA looks
                # identical to the pipeline hanging, which is exactly how it read on
                # 2026-08-07: the PA finished, nothing followed, and there was no line
                # explaining that post-PA was simply not selected.
                self._pl_log("PA pipeline: post-PA NOT selected -- stopping here. Tick "
                             "'include in Run Pipeline' on the After-PA card, or click "
                             "'Run After-PA Captures' now.")
            if self._pl_pa_sel_postpa.get():
                if self._pl_pa_sel_pajob.get() and not pa_ran:
                    return                                   # unreachable, kept explicit
                if not self._pl_pa_sel_pajob.get():
                    self._pl_log("PA pipeline: post-PA requested WITHOUT the PA job in this "
                                 "run -- assuming the PA was fired separately.")
                # Re-normalize the rig first. PA Setup left it at dichroic OUT and the PA
                # zoom; without this the post-PA set is captured through the PA optics and
                # the mismatch reads as photoactivation. The explicit light-path restore in
                # _pl_prepa_capture_ocs still runs after it -- init_rig sets the rig's
                # known-good state, the restore reproduces THIS run's pre-PA state.
                if not self._pl_pipeline_init_rig(wd, "post-PA"):
                    self.after(0, lambda: self._pl_pa_status.configure(
                        text="Pipeline: aborted before the post-PA captures.", fg=RED)); return
                self.after(0, lambda: self._pl_pa_status.configure(
                    text="Pipeline: After-PA captures...", fg=YELLOW))
                self._pl_log("PA pipeline: After-PA viz capture(s)")
                if not self._pl_prepa_capture_ocs("postPA"):
                    return                                   # helper showed the status
            done_bits = [f"{len(seq)} macro(s)"]
            if self._pl_pa_sel_pajob.get():
                done_bits.append("PA job")
            if self._pl_pa_sel_postpa.get():
                done_bits.append("post-PA")
            self._pl_log(f"PA pipeline: COMPLETE -- {', '.join(done_bits)}")
            self.after(0, lambda b=", ".join(done_bits): self._pl_pa_status.configure(
                text=f"Pipeline complete: {b}.", fg=GREEN))
        except Exception as exc:
            # This runs on a daemon thread. Without this handler an exception unwinds out
            # of the thread and is DISCARDED -- no console, no log, no dialog: the sequence
            # simply stops mid-run and looks like a hang. On 2026-08-07 a Run All ended
            # after the PA with no post-PA and no error line, and there was no way to tell
            # a crash from a skipped step. Log the full traceback so there is always
            # something to read.
            import traceback
            self._pl_log(f"PA pipeline: CRASHED -- {type(exc).__name__}: {exc}")
            for ln in traceback.format_exc().splitlines():
                self._pl_log(f"    {ln}")
            self.after(0, lambda e=exc: self._pl_pa_status.configure(
                text=f"Pipeline CRASHED — {type(e).__name__}: {e} (see log)", fg=RED))
        finally:
            self._pl_dispatch_busy = False

    def _pl_toggle_beer_lambert(self):
        """P0 is the base power in both modes (the flat fixed power when off), so it
        stays editable; only the depth-decay length L is Beer-Lambert-only."""
        st = "normal" if self._pl_beer_lambert.get() else "disabled"
        try:
            self._pl_bl_fields["l"].configure(state=st)
        except Exception:
            pass

    def _pl_generate_bins_thread(self):
        threading.Thread(target=self._pl_generate_bins, daemon=True).start()

    def _pl_generate_bins(self):
        try:
            import spheroid_pipeline as _pl
        except ImportError as e:
            messagebox.showerror("Import error", str(e)); return
        if not self._pl_records:
            messagebox.showwarning("No records", "Run Steps 1-3 first."); return
        out = self._pl_out_dir.get().strip()
        if not out:
            messagebox.showwarning("No output dir", "Set output directory."); return
        bin_dir = Path(out) / "bins"
        ch_field = self._pl_ch_field.get().strip() or "CH2LaserPower"
        if self._pl_beer_lambert.get():
            try:
                P0    = float(self._pl_P0.get())
                L_um  = float(self._pl_L_um.get())
            except ValueError:
                messagebox.showerror("Bad value", "P0 and L must be numbers."); return
            mode_txt = f"Beer-Lambert depth ramp (P0={P0:g}%, L={L_um:g}um)"
        else:
            # Fixed power: flat bin at P0 (%). L -> inf makes exp(depth/L)=1, so
            # every plane gets the same power.
            try:
                P0 = float(self._pl_P0.get())
            except ValueError:
                messagebox.showerror("Bad value", "P0 (%) must be a number for fixed-power bins."); return
            L_um = math.inf
            mode_txt = f"fixed power {P0:g}% (flat, P0)"
        self._pl_log(f"Generate bins: {mode_txt}")
        ref_s   = self._pl_ref_bin.get().strip()
        ref_bin = Path(ref_s) if ref_s else None
        if ref_bin and not ref_bin.exists():
            self._pl_log(f"Reference .bin not found: {ref_bin}"); ref_bin = None
        if ref_bin is None:
            self._pl_log("WARNING: no reference .bin set — generated bins will be flagged "
                         "'Incompatible Z Correction and Camera' by NIS-E. Browse a rig export.")
        # Fall back to the Step-3 Middle plane Z / half / step when a record has no
        # recorded z_centre (e.g. no capture has recorded one yet -- Capture auto-refreshes
        # it) -- the spheroids are all captured centred on that plane, so it IS the correct
        # bin z-centre.
        try: gui_z = float(self._pl_z_centre.get())
        except Exception: gui_z = 0.0
        try: gui_h = float(self._pl_z_half.get())
        except Exception: gui_h = 90.0
        try: gui_s = float(self._pl_z_step.get())
        except Exception: gui_s = 10.0
        n_ok = 0
        for r in self._pl_records:
            if r.z_centre_um == 0.0 and gui_z > 0.0:
                r.z_centre_um = gui_z
                if not r.z_half_um: r.z_half_um = gui_h
                if not r.z_step_um: r.z_step_um = gui_s
            if r.z_centre_um == 0.0:
                self._pl_log(f"Rank {r.rank}: skipping bin — Z not recorded "
                             "(set Step-3 Middle plane Z, or Capture to record it)"); continue
            try:
                _pl.generate_bin(r, P0, L_um, ch_field, bin_dir, reference_bin=ref_bin)
                n_ok += 1
                self._pl_log(f"Rank {r.rank}: bin written → {Path(r.bin_path).name}")
            except Exception as exc:
                self._pl_log(f"Rank {r.rank}: bin error — {exc}")
        self.after(0, self._pl_update_table)
        self.after(0, self._pl_refresh_dashboard)
        msg = f"Bins generated: {n_ok}/{len(self._pl_records)} — {mode_txt}"
        self.after(0, lambda: self._pl_capture_lbl.configure(text=msg, fg=GREEN))
        self._pl_log(msg)
        self._pl_update_step_dots()

    def _pl_capture_queue_thread(self):
        threading.Thread(target=self._pl_capture_queue, daemon=True).start()

    def _pl_capture_queue(self):
        try:
            import spheroid_pipeline as _pl
        except ImportError as e:
            messagebox.showerror("Import error", str(e)); return
        # Prefer the session.ini work_dir (where nis_macro_capture_zcorrected.mac polls)
        # so the GUI and macro always agree on the trigger location, matching Step 3.
        trigger_dir = self._pl_trigger_dir.get().strip()
        try:
            import configparser
            cp = configparser.ConfigParser(); cp.read(_pl.SESSION_INI)
            wd = cp.get("paths", "work_dir", fallback="").strip()
            if wd:
                trigger_dir = wd
        except Exception:
            pass
        nd2_out_dir = self._pl_nd2_out_dir.get().strip()
        if not trigger_dir or not nd2_out_dir:
            messagebox.showwarning("Missing dirs", "Set trigger dir and ND2 output dir."); return
        ready = [r for r in self._pl_records if r.status == "BIN_READY"]
        if not ready:
            messagebox.showwarning("No bins", "Generate bins first (Step 4)."); return

        total = len(ready)
        for i, rec in enumerate(ready):
            self._pl_log(f"Triggering rank {rec.rank} ({i+1}/{total}) ...")
            self.after(0, lambda r=rec, n=i+1, t=total:
                       self._pl_capture_lbl.configure(
                           text=f"Imaging rank {r.rank} ({n}/{t}) — waiting for NIS-E ...",
                           fg=PEACH))
            trigger_path = _pl.write_trigger(rec, Path(trigger_dir), Path(nd2_out_dir))
            self._pl_log(f"  trigger written: {trigger_path.name}")
            self.after(0, self._pl_update_table)

            done = _pl.wait_for_done(Path(trigger_dir), timeout_s=600.0)
            if done is None:
                self._pl_log(f"  Rank {rec.rank}: TIMEOUT — no done.ini in 10 min")
                rec.status = "FAILED"
            elif done.get("status", "ok") != "ok":
                rec.status = "FAILED"
                self._pl_log(f"  Rank {rec.rank}: macro reported status={done.get('status')}")
            else:
                rec.status     = "IMAGED"
                rec.nd2_out_path = done.get("nd2_path", rec.nd2_out_path)
                self._pl_log(f"  Rank {rec.rank}: IMAGED -> {done.get('nd2_path', '?')}")
            self.after(0, self._pl_update_table)
            self.after(0, self._pl_refresh_dashboard)

        done_count = sum(1 for r in ready if r.status == "IMAGED")
        msg = f"Capture queue done: {done_count}/{total} imaged."
        self.after(0, lambda: self._pl_capture_lbl.configure(text=msg, fg=GREEN))
        self._pl_log(msg)

        # Save state CSV for post-hoc reference
        out = self._pl_out_dir.get().strip()
        if out:
            import spheroid_pipeline as _pl2
            _pl2.export_state_csv(self._pl_records, Path(out) / "pipeline_state.csv")
            self._pl_log(f"State CSV written → {out}/pipeline_state.csv")

    # ── _btn grid-mode overload ────────────────────────────────────────────────

    # The WellSelection task's JSON parameter key for Jobs_RunJobInitParam is UNDOCUMENTED.
    # Discover it ONCE on the rig with the NIS-E Debug task (insert it after the "Select Wells"
    # task in Step1_Locate_via_scan, pick the selection parameter, run with a known well, read
    # the printed "<Task>.<Param>" key + the value FORMAT it shows), then set both below.
    #   WELL_PARAM_KEY -- the JSON key. The Debug task shows the read path as
    #     "Job.WellSelection.Selection.Wells[0].Name"; the InitParam JSON key drops the "Job."
    #     prefix (doc keys like "OCSel.OptConf" have none). So we set Wells[0].Name -- the
    #     first (and only) selected well's name.
    #   _pl_well_param_value() -- the value is just the well name (e.g. "B02").
    # If setting Name alone does not MOVE the scan (task also keys off Index/Position), the
    # selection is too structured for InitParam and we pivot to the point-set approach.
    # Empty key -> the job launches WITHOUT the well (plain Jobs_RunJobByName) -- safe no-op.
    WELL_PARAM_KEY = "WellSelection.Selection.Wells[0].Name"      # Step1 mosaic: RELABELS output file only
    WELL_INDEX_KEY = "WellSelection.Selection.Wells[0].Index"     # NOT settable -- rebuilt from the UI selection at run time
    WELL_NCOLS     = 12   # plate columns (Cellvis 96-well) for the well<->index numbering
    PLATE_CALIB_INI = r"C:/SpheroidPA/plate_calib.ini"           # 3-well affine calib: well -> stage XY (survives GUI restart)
    # step3_zstack_PA_WC (PA) is parameterized too. The ORIGINAL step3_zstack_PA had an
    # "Import points from ND acquisition" (NDToPointSet, Replace Points) task that overwrote
    # PredefinedPoints at run time from the LAST ND acquisition -> the 850 nm laser fired on a
    # STALE spot (observed 2026-07-27: both PA nd2 landed at the last-imaged coord regardless of
    # the stored point; the in-job point was even seen at -50 mm). step3_zstack_PA_WC has that
    # task DELETED, so the manual PredefinedPoints is authoritative and we set the point (and
    # well) DIRECTLY at launch via Jobs_RunJobInitParam. Keys confirmed via the job's Debug task (drop the "Job."
    # prefix, per Example 18/19). X/Y are native stage um (Position.Stage.x/y); Z is um (Position.z,
    # a sibling of Stage -- not Stage.z). Single point (Positions[0]) for now.
    PA_POINT_X_KEY = "PredefinedPoints.Positions.Positions[0].Position.Stage.x"
    PA_POINT_Y_KEY = "PredefinedPoints.Positions.Positions[0].Position.Stage.y"
    PA_POINT_Z_KEY = "PredefinedPoints.Positions.Positions[0].Position.z"
    # Per-point templates for multi-spheroid PA -- the [0] constants above are these at i=0.
    # Note z is Position.z, a SIBLING of Position.Stage, NOT Position.Stage.z: rig-confirmed
    # 2026-08-07 off the Debug task's parameter tree ('.Stage.z' throws Invalid Expression).
    PA_POINT_X_FMT = "PredefinedPoints.Positions.Positions[{i}].Position.Stage.x"
    PA_POINT_Y_FMT = "PredefinedPoints.Positions.Positions[{i}].Position.Stage.y"
    PA_POINT_Z_FMT = "PredefinedPoints.Positions.Positions[{i}].Position.z"
    # Point-array length. In the Debug parameter tree PredefinedPoints' ArraySize renders
    # NORMAL while WellSelection's renders RED (red = read-only output), which suggests this
    # one is settable -- UNVERIFIED on the rig. Sent so the job's point count matches what we
    # parameterize. If NIS ignores it, the job keeps whatever count is defined by hand, and
    # any surplus points fire the 850 nm laser at their STORED coordinates -- which is why
    # _pl_run_pa_job hard-confirms the count with the operator before firing.
    PA_POINT_COUNT_KEY = "PredefinedPoints.Positions.Positions.ArraySize"
    # Where the JOB writes its output. Read off the Debug task's parameter tree as
    # Job.StoreToFsOnly.Folder (the 'Alternative Storage Location' task); the Job. prefix is
    # dropped for Jobs_RunJobInitParam, exactly like the PredefinedPoints keys.
    #
    # This closes the last stale-path hole. The folder used to live ONLY on the job, so it
    # silently kept whatever the previous session set: PA output landed in another user's
    # D:\Images S636B\Brandon\... folder on 2026-08-07, and later in a shared work/PA while
    # the GUI said work/0807_4/pa. Sending it per run means the GUI's folder IS the
    # destination, with no second place to keep in sync.
    #
    # Kept on C: is NOT an option -- removing the storage task entirely would send ~11 GB
    # per run to C:\ProgramData\...\jobsdb (180 GB free, already 17 GB of Step1 mosaics),
    # versus 14 TB on the share.
    PA_STORE_FOLDER_KEY = "StoreToFsOnly.Folder"

    def _pl_well_param_value(self, well):
        """Format the GUI well id into the value WellSelection Wells[0].Name expects. The plate
        names wells WITHOUT a leading zero (the Debug task showed 'B1', not 'B01'), so normalize
        e.g. 'B02' -> 'B2', 'a01' -> 'A1'. Unrecognised formats pass through unchanged."""
        import re
        m = re.match(r"^\s*([A-Za-z]+)\s*0*(\d+)\s*$", well)
        return f"{m.group(1).upper()}{int(m.group(2))}" if m else well

    def _pl_well_index_value(self, well):
        """Row-major 0-based well INDEX for WellSelection.Selection.Wells[0].Index -- the field
        that actually repositions the Step 1 mosaic (Name only relabels). Rig-verified 2026-08-07:
        B1 -> 12, i.e. Index = row_idx*WELL_NCOLS + (col-1), row A=0..H=7. So A1=0, A2=1, B2=13,
        H12=95. Returns None if the well can't be parsed (single-letter rows only)."""
        import re
        m = re.match(r"^\s*([A-Za-z])\s*0*(\d+)\s*$", well)
        if not m:
            return None
        return (ord(m.group(1).upper()) - ord("A")) * self.WELL_NCOLS + (int(m.group(2)) - 1)

    # ── Plate coordinate calibration (well -> stage XY) ────────────────────────
    # WellSelection can't be parameterized: its Wells[] are rebuilt from the task's internal UI
    # selection at run time (rig-proven -- setting Index/Position is ignored; only .Name relabels).
    # So the mosaic is driven by a stage COORDINATE instead -- an affine map (origin + column/row
    # pitch vectors, rotation-aware) fit from 3 reference wells, set as the job's PredefinedPoint
    # (the SAME param proven settable on the PA job). Persisted to PLATE_CALIB_INI.
    @staticmethod
    def _plate_well_rc(well):
        """'A02'/'a2' -> (row_idx 0-based, col 1-based); None if unparseable (single-letter rows)."""
        import re
        m = re.match(r"^\s*([A-Za-z])\s*0*(\d+)\s*$", str(well))
        return (ord(m.group(1).upper()) - ord("A"), int(m.group(2))) if m else None

    @staticmethod
    def _solve3(A, b):
        """Solve 3x3  A x = b  by Cramer's rule (numpy is optional in this GUI). None if singular."""
        def d3(m):
            return (m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])
                  - m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])
                  + m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]))
        D = d3(A)
        if abs(D) < 1e-9:
            return None
        res = []
        for j in range(3):
            Aj = [row[:] for row in A]
            for i in range(3):
                Aj[i][j] = b[i]
            res.append(d3(Aj) / D)
        return res

    @classmethod
    def _plate_fit(cls, refs):
        """Least-squares affine fit over ALL provided wells (>=3). Returns {bx,cx,rx,by,cy,ry} or
        None. Model: x = bx + (col-1)*cx + row*rx ; y = by + (col-1)*cy + row*ry. Solved via the
        normal equations (AtA is 3x3 -> Cramer), so no numpy dependency; using >3 wells averages
        out plate non-idealities (corner-to-corner the grid isn't perfectly affine)."""
        pts = []
        for w, x, y in refs:
            rc = cls._plate_well_rc(w)
            if rc is None:
                return None
            pts.append(([1.0, rc[1] - 1, rc[0]], float(x), float(y)))
        if len(pts) < 3:
            return None
        AtA = [[0.0]*3 for _ in range(3)]
        Atx = [0.0]*3; Aty = [0.0]*3
        for a, x, y in pts:
            for i in range(3):
                Atx[i] += a[i]*x; Aty[i] += a[i]*y
                for j in range(3):
                    AtA[i][j] += a[i]*a[j]
        sx = cls._solve3([r[:] for r in AtA], Atx)
        sy = cls._solve3([r[:] for r in AtA], Aty)
        if sx is None or sy is None:
            return None
        return {"bx": sx[0], "cx": sx[1], "rx": sx[2], "by": sy[0], "cy": sy[1], "ry": sy[2]}

    @classmethod
    def _plate_well_xy(cls, calib, well):
        """Stage (x,y) for a well from a calib dict; None if no calib / unparseable well."""
        if not calib:
            return None
        rc = cls._plate_well_rc(well)
        if rc is None:
            return None
        row, col = rc
        return (calib["bx"] + (col - 1)*calib["cx"] + row*calib["rx"],
                calib["by"] + (col - 1)*calib["cy"] + row*calib["ry"])

    def _pl_load_plate_calib(self):
        """Load the persisted 3-well affine calib from PLATE_CALIB_INI. Caches on
        self._pl_plate_calib and returns it (or None if absent/invalid)."""
        import configparser
        try:
            p = Path(self.PLATE_CALIB_INI)
            if p.exists():
                cp = configparser.ConfigParser(); cp.read(p)
                refs = [(cp.get(s, "well"), cp.getfloat(s, "x"), cp.getfloat(s, "y"))
                        for s in (f"ref{i}" for i in range(1, 10)) if cp.has_section(s)]
                if len(refs) >= 3:
                    self._pl_plate_calib = self._plate_fit(refs)
                    return self._pl_plate_calib
        except Exception as e:
            self._pl_log(f"Plate calib: load failed: {e}")
        self._pl_plate_calib = None
        return None

    def _pl_save_plate_calib(self, refs):
        """Fit + persist 3 reference wells. refs: [(well,x,y)]*3. Returns the calib dict or None."""
        calib = self._plate_fit(refs)
        if calib is None:
            self._pl_log("Plate calib: fit failed -- pick 3 wells that span BOTH rows and columns.")
            return None
        import configparser
        cp = configparser.ConfigParser()
        for i, (w, x, y) in enumerate(refs, 1):
            cp[f"ref{i}"] = {"well": str(w).strip().upper(), "x": f"{float(x):.3f}", "y": f"{float(y):.3f}"}
        try:
            p = Path(self.PLATE_CALIB_INI); p.parent.mkdir(parents=True, exist_ok=True)
            with p.open("w") as f:
                cp.write(f)
        except Exception as e:
            self._pl_log(f"Plate calib: save failed: {e}"); return None
        self._pl_plate_calib = calib
        cpit = math.hypot(calib["cx"], calib["cy"]); rpit = math.hypot(calib["rx"], calib["ry"])
        self._pl_log(f"Plate calib saved: A1 origin=({calib['bx']:.1f},{calib['by']:.1f}) um, "
                     f"col-pitch={cpit:.1f}, row-pitch={rpit:.1f} um (expect ~9000).")
        return calib

    # Number of reference-well rows offered in the calibration panel. 3 is the minimum for an
    # affine fit; the extra rows are optional and simply average out plate non-idealities.
    PLATE_CALIB_ROWS = 4

    def _pl_build_plate_calib_panel(self, parent):
        """Build the 'Plate Calibration' panel: N (well, stage x, stage y) rows + Compute & Save.

        How the operator fills it: in NIS-E, click a well in the job's WellSelection task and read
        that well's `...Position.Stage.x/y` off the Debug task (or move the stage to the well
        centre and read the stage readout). Enter >=3 wells that span BOTH rows and columns --
        corner wells (A1 / A12 / H1 / H12) give the best-conditioned fit. Rows left blank are
        ignored. 'Compute & Save' fits the affine and writes PLATE_CALIB_INI."""
        box = tk.LabelFrame(parent, text=" Plate Calibration (well -> stage XY) ",
                            bg=BG, fg=MAUVE, font=("Segoe UI", 9, "bold"),
                            bd=1, relief="groove")
        box.pack(fill="x", pady=(0, 4))
        tk.Label(box, text="Reference wells: enter >=3 spanning both rows and columns "
                           "(corners are best). Stage XY in um.",
                 bg=BG, fg=SUBTEXT, font=("Segoe UI", 8), anchor="w").pack(fill="x", padx=6, pady=(2, 0))

        hdr = tk.Frame(box, bg=BG); hdr.pack(fill="x", padx=6)
        for txt, w in (("Well", 8), ("Stage X (um)", 14), ("Stage Y (um)", 14)):
            tk.Label(hdr, text=txt, width=w, anchor="w", bg=BG, fg=TEXT2,
                     font=("Segoe UI", 8, "bold")).pack(side="left", padx=(0, 4))

        self._pl_calib_rows = []
        for _ in range(self.PLATE_CALIB_ROWS):
            r = tk.Frame(box, bg=BG); r.pack(fill="x", padx=6, pady=1)
            vs = []
            for w in (8, 14, 14):
                v = tk.StringVar()
                tk.Entry(r, textvariable=v, width=w, bg=SURFACE, fg=TEXT,
                         insertbackground=TEXT, relief="flat",
                         font=("Segoe UI", 9)).pack(side="left", padx=(0, 4))
                vs.append(v)
            self._pl_calib_rows.append(tuple(vs))

        btn = tk.Frame(box, bg=BG); btn.pack(fill="x", padx=6, pady=(3, 2))
        self._btn(btn, "Compute & Save", self._pl_calib_compute_save, BLUE, "#1e1e2e", side="left")
        self._btn(btn, "Reload from file", self._pl_calib_fill_from_file, SURFACE2, TEXT, side="left")
        self._btn(btn, "Preview Well ID", self._pl_calib_preview, SURFACE2, TEXT, side="left")

        self._pl_calib_lbl = tk.Label(box, text="", bg=BG, fg=SUBTEXT,
                                      font=("Segoe UI", 8), anchor="w", justify="left")
        self._pl_calib_lbl.pack(fill="x", padx=6, pady=(0, 4))
        self._pl_calib_fill_from_file()

    def _pl_calib_read_rows(self):
        """Rows -> [(well, x, y)], skipping blank rows. Raises ValueError naming the bad row."""
        refs = []
        for i, (wv, xv, yv) in enumerate(self._pl_calib_rows, 1):
            w, x, y = wv.get().strip(), xv.get().strip(), yv.get().strip()
            if not (w or x or y):
                continue
            if not (w and x and y):
                raise ValueError(f"row {i} is incomplete -- fill well, X and Y (or clear it).")
            if self._plate_well_rc(w) is None:
                raise ValueError(f"row {i}: '{w}' is not a well ID (expected e.g. A1, A01, H12).")
            try:
                refs.append((w.upper(), float(x), float(y)))
            except ValueError:
                raise ValueError(f"row {i}: X/Y must be numbers (got '{x}', '{y}').")
        return refs

    def _pl_calib_fill_from_file(self):
        """Prefill the rows from PLATE_CALIB_INI and show the fit that is currently in force."""
        import configparser
        try:
            p = Path(self.PLATE_CALIB_INI)
            refs = []
            if p.exists():
                cp = configparser.ConfigParser(); cp.read(p)
                for s in (f"ref{i}" for i in range(1, 10)):
                    if cp.has_section(s):
                        refs.append((cp.get(s, "well"), cp.get(s, "x"), cp.get(s, "y")))
            for (wv, xv, yv), val in zip(self._pl_calib_rows,
                                         refs + [("", "", "")] * len(self._pl_calib_rows)):
                wv.set(val[0]); xv.set(val[1]); yv.set(val[2])
            if len(refs) > len(self._pl_calib_rows):
                # The ini holds more refs than the panel can show. Say so loudly: a later
                # 'Compute & Save' rewrites the file from the ROWS, which would silently
                # discard the ones that never made it on screen.
                self._pl_log(f"Plate calib: {self.PLATE_CALIB_INI} has {len(refs)} reference "
                             f"wells but the panel shows {len(self._pl_calib_rows)} -- "
                             f"'Compute & Save' would DROP the rest. Raise PLATE_CALIB_ROWS "
                             f"or edit the ini directly.")
        except Exception as e:
            self._pl_calib_lbl.configure(text=f"Could not read {self.PLATE_CALIB_INI}: {e}", fg=RED)
            return
        self._pl_calib_show(self._pl_load_plate_calib(), loaded=True)

    def _pl_calib_show(self, calib, loaded=False):
        """Render the fit summary: pitches (expect ~9000 um on a 96-well plate), rotation, and the
        worst reference-well residual -- the number that says whether the fit can be trusted."""
        if not calib:
            self._pl_calib_lbl.configure(
                text="No calibration -- Run Job1 cannot drive by coordinate (it will image the "
                     "job's stored point). Enter >=3 wells and Compute & Save.", fg=YELLOW)
            return
        cpit = math.hypot(calib["cx"], calib["cy"]); rpit = math.hypot(calib["rx"], calib["ry"])
        rot  = math.degrees(math.atan2(calib["cy"], -calib["cx"]))
        worst = ""
        try:
            refs = self._pl_calib_read_rows()
            if refs:
                res = max(math.hypot(*(a - b for a, b in
                                       zip(self._plate_well_xy(calib, w), (x, y))))
                          for w, x, y in refs)
                worst = f" | worst residual {res:.1f} um"
        except ValueError:
            pass
        self._pl_calib_lbl.configure(
            text=(f"{'Loaded' if loaded else 'Saved'}: A1 origin "
                  f"({calib['bx']:.1f}, {calib['by']:.1f}) um | col-pitch {cpit:.1f} | "
                  f"row-pitch {rpit:.1f} um (expect ~9000) | rotation {rot:+.3f} deg{worst}"),
            fg=GREEN)

    def _pl_calib_compute_save(self):
        """Fit + persist the entered reference wells, then refresh the summary."""
        try:
            refs = self._pl_calib_read_rows()
        except ValueError as e:
            messagebox.showwarning("Plate calibration", str(e)); return
        if len(refs) < 3:
            messagebox.showwarning("Plate calibration",
                                   "Need at least 3 reference wells spanning both rows and "
                                   "columns."); return
        calib = self._pl_save_plate_calib(refs)
        if calib is None:
            self._pl_calib_lbl.configure(
                text="Fit failed -- the wells are collinear (all in one row or one column). "
                     "Pick wells that span BOTH.", fg=RED)
            return
        self._pl_calib_show(calib)

    def _pl_calib_preview(self):
        """Show where the CURRENT Well ID maps to, so the operator can sanity-check the fit
        against the stage readout BEFORE firing a mosaic at it."""
        well  = self._pl_well_id.get().strip()
        calib = getattr(self, "_pl_plate_calib", None) or self._pl_load_plate_calib()
        if not calib:
            self._pl_calib_lbl.configure(text="No calibration loaded -- Compute & Save first.",
                                         fg=YELLOW); return
        xy = self._plate_well_xy(calib, well)
        if xy is None:
            self._pl_calib_lbl.configure(
                text=f"Well ID '{well or '(blank)'}' is not a well (expected e.g. A02).", fg=RED)
            return
        self._pl_calib_lbl.configure(
            text=f"Well {well.upper()} -> stage ({xy[0]:.1f}, {xy[1]:.1f}) um  "
                 f"-- Run Job1 will centre the mosaic here.", fg=BLUE)
        self._pl_log(f"Plate calib preview: {well.upper()} -> ({xy[0]:.1f}, {xy[1]:.1f}) um")

    def _pl_well_params(self, key, tag):
        """Build {key: value} (the WellSelection Name) from the shared _pl_well_id. NOTE: this only
        RELABELS the mosaic's output file -- WellSelection can't be repositioned by param (Index and
        Position are rebuilt from the UI selection at run time). Coordinate drive lives in
        _pl_run_job1_at_well. None when no well/key."""
        well = self._pl_well_id.get().strip()
        if well and key:
            return {key: self._pl_well_param_value(well)}
        return None

    def _pl_pa_points_from_triggers(self):
        """Read EVERY af_trigger_*.ini into [(name, x, y, z)], ordered by rank.

        Split out from _pl_pa_point_params so the point list can be counted and shown to the
        operator (and unit-tested) without building the param dict. Triggers with unreadable
        coords are DROPPED with a log line rather than defaulting to 0 -- a 0 would be stage
        origin, i.e. the centre of the plate."""
        import configparser, spheroid_pipeline as _pl
        cp0 = configparser.ConfigParser()
        try:
            if _pl.SESSION_INI.exists():
                cp0.read(_pl.SESSION_INI)
        except Exception:
            pass
        work = cp0.get("paths", "work_dir", fallback="").strip() or self._pl_trigger_dir.get().strip()
        if not work:
            return None, "no work dir"
        pts = []
        for t in sorted(Path(work).glob("af_trigger_*.ini")):
            cp = configparser.ConfigParser()
            try:
                cp.read(t)
                pts.append((t.name,
                            float(cp.get("spheroid", "stage_x")),
                            float(cp.get("spheroid", "stage_y")),
                            float(cp.get("spheroid", "z_centre"))))
            except Exception as e:
                self._pl_log(f"PA: SKIPPING {t.name} -- unreadable coords ({e}).")
        if not pts:
            return None, "no readable af_trigger_*.ini"
        return pts, None

    def _pl_pa_point_params(self):
        """Build the JOB parameters that set step3_zstack_PA_WC's PredefinedPoints to the
        CURRENT spheroid selection, so the 850 nm PA fires on exactly the supplied points
        instead of the job's stale stored ones.

        MULTI-POINT: every af_trigger_*.ini becomes Positions[i], each carrying its own
        Stage.x/y AND its own z. Per-point z matters because ZStackDefinition is
        Relative-to-Home and 'Move to Point' (Move to Z ticked) re-anchors Home at each
        point -- so every spheroid gets its stack centred on its own depth. Today all
        triggers carry the same z_centre; they diverge once Step 3 writes per-record depths
        (trigger_autofocus_all(per_record_z=True)), and this code needs no change then.

        Returns the params dict, or None -- and None means the job runs on its STORED points,
        which is why the caller confirms before firing."""
        pts, err = self._pl_pa_points_from_triggers()
        if pts is None:
            self._pl_log(f"PA: {err} -- launching WITHOUT point params (the job fires on its "
                         f"STORED PredefinedPoints, which may be stale).")
            return None
        params = {}
        for i, (name, x, y, z) in enumerate(pts):
            params[self.PA_POINT_X_FMT.format(i=i)] = x
            params[self.PA_POINT_Y_FMT.format(i=i)] = y
            params[self.PA_POINT_Z_FMT.format(i=i)] = z
        # ArraySize is sent ONLY when growing the array beyond the single point the job
        # already holds. Rig evidence 2026-08-07: three launches sending x/y/z + well all
        # produced output; the two that additionally sent ArraySize produced NOTHING -- a
        # silent no-op, no jobrun, no nd2, while cmd_done still said ok. An unrecognised key
        # appears to make Jobs_RunJobInitParam reject the WHOLE parameter set, which then
        # falls back to the job's stored point. So for the single-point case (the norm) send
        # exactly the payload that is proven to work, and accept the extra risk only when
        # multi-point genuinely needs the array grown.
        if len(pts) > 1:
            params[self.PA_POINT_COUNT_KEY] = len(pts)
            self._pl_log(f"PA: sending {self.PA_POINT_COUNT_KEY}={len(pts)} -- UNPROVEN key; "
                         f"if the job produces no output, this is the first thing to drop.")
        well = self._pl_well_id.get().strip()
        if well:
            params[self.WELL_PARAM_KEY] = self._pl_well_param_value(well)

        # Drive the JOB's output folder from the GUI. Sent as a NATIVE Windows path
        # (backslashes) to match what the job's own dialog shows -- the value goes straight
        # into a NIS-E path field, not through Python.
        store = self._pl_pa_effective_save_dir()
        if store is not None:
            params[self.PA_STORE_FOLDER_KEY] = str(Path(store))

        # Log EVERY point, not just a count: this is the last human-readable record of what
        # the laser was told to hit, and the Debug breakpoint that used to show it is gone.
        self._pl_log(f"PA point params: {len(pts)} point(s), well={well or '(none)'}, "
                     f"store={str(Path(store)) if store else '(job default)'}")
        for i, (name, x, y, z) in enumerate(pts):
            self._pl_log(f"   [{i}] {name}: Stage.x/y=({x:.1f}, {y:.1f}) um  z={z:.1f} um")
        if len({round(p[3], 1) for p in pts}) == 1 and len(pts) > 1:
            self._pl_log(f"   NOTE: all {len(pts)} points share z={pts[0][3]:.1f} um. Enable "
                         f"per-record Z (trigger_autofocus_all(per_record_z=True)) once the "
                         f"auto z-plane is validated, or every stack is centred on one depth.")
        return params

    def _pl_run_pa_job(self):
        """Step 4 'Run Photoactivation Job': regenerate the current triggers, then launch
        step3_zstack_PA_WC with its PredefinedPoint[0] + well set to the CURRENT spheroid via
        Jobs_RunJobInitParam -- the activation fires on the exact supplied point, never the job's
        stale in-job point. Still gated (A1 confirm + fire-confirm) inside _pl_run_job."""
        if not self._pl_regen_triggers():
            return
        # Prepare the PA folder automatically -- no button, no click to forget. Creates
        # <Output directory>/pa if absent and records it in pa_trigger.ini for the macros.
        self._pl_pa_sync_save_dir()
        self._pl_pa_prepare_save_dir()
        params = self._pl_pa_point_params()
        # Multi-point confirmation. We cannot read the job's PredefinedPoints count back, and
        # the mismatch is dangerous in one direction: if the job holds MORE points than we
        # parameterise, the surplus keep their STORED coordinates and the 850 nm laser fires
        # at them too -- off-target, on a sample that cannot be re-run. ArraySize is sent to
        # size the array, but that key is unverified, so a human confirms the count.
        pts, _ = self._pl_pa_points_from_triggers()
        n = len(pts) if pts else 0
        if n > 1:
            lines = "\n".join(f"   [{i}] ({x:.1f}, {y:.1f}) um   z={z:.1f}"
                              for i, (_nm, x, y, z) in enumerate(pts))
            if not messagebox.askyesno(
                    "Multi-spheroid PA",
                    f"This will photoactivate {n} points:\n\n{lines}\n\n"
                    f"step3_zstack_PA_WC's PredefinedPoints task must hold EXACTLY {n} "
                    f"point(s).\n\nIf the job holds more, the extra ones fire at their STORED "
                    f"coordinates — off target, and not repeatable on this sample.\n\n"
                    f"Confirmed the job has {n} point(s)?"):
                self._pl_log(f"PA: launch cancelled -- operator did not confirm the "
                             f"{n}-point job configuration.")
                return
        self._pl_run_job("step3_zstack_PA_WC", True, params=params)

    # Where the Step1_Locate_via_scan_WC JOB writes its mosaic -- a NEW timestamped folder per run,
    # under the NIS-E jobs DB (NOT the pipeline nd2_dir). See the [[nise-step1-locate-scan-save-path]] note.
    MOSAIC_JOB_DIR = r"C:/ProgramData/Laboratory Imaging/Jobs/jobsdb Projects/IMAGEN/Step1_Locate_via_scan_WC"

    # The widefield optical config the Step-1 mosaic must be taken on. Selected explicitly
    # before every Run Job1 so the mosaic is correct from ANY starting state -- brightfield,
    # 2-photon, whatever the previous step left behind. The job's own CaptureDefinition
    # names this OC, but it does NOT switch to it: on 2026-08-07 a mosaic fired straight
    # after a 1050 nm 2P capture came back unilluminated (max 356 vs 7287, and the nd2 had
    # no MultiLaser(EPI) entry at all), because the 555 line was never turned on. The
    # 18:11 mosaic on the same job worked purely because the rig happened to already be on
    # widefield. Asserting the OC here removes that dependency on prior state.
    MOSAIC_OC = "Flash 4.0 - EPI + Tucam(In) + Twincam(Out):Spectra3_555"

    def _pl_run_job1_at_well(self):
        """Step 1 'Run Job1': launch Step1_Locate_via_scan_WC at the GUI Well ID. WellSelection
        can't be repositioned by param, so when a plate calibration exists we compute the well's
        stage (x,y) and drive the mosaic via the job's PredefinedPoint (Jobs_RunJobInitParam);
        otherwise we fall back to the Name param (RELABELS only -- won't move) and warn. After
        launching, watches the JOB's output dir and auto-fills the mosaic ND2 field."""
        well  = self._pl_well_id.get().strip()
        calib = getattr(self, "_pl_plate_calib", None) or self._pl_load_plate_calib()
        xy    = self._plate_well_xy(calib, well) if (calib and well) else None
        if xy is not None:
            # PredefinedPoint ONLY. Do NOT add WELL_PARAM_KEY here: the WellSelection task was
            # DELETED from Step1_Locate_via_scan_WC (that deletion is the whole point of the
            # coordinate drive), so a "WellSelection.*" key now names a task the job no longer
            # has. The mosaic's filename therefore loses its well label -- harmless, because
            # _pl_watch_mosaic picks the LARGEST nd2 in the run folder, not one matched by name.
            params = {self.PA_POINT_X_KEY: round(xy[0], 3),
                      self.PA_POINT_Y_KEY: round(xy[1], 3)}
            self._pl_log(f"Run Job1: well {well} -> stage ({xy[0]:.1f}, {xy[1]:.1f}) um "
                         f"via plate calib (PredefinedPoint drive).")
        else:
            # No calib -> nothing we can send. The Name param can't help either (no WellSelection
            # task to relabel), so launch un-parameterized and say plainly where it will land.
            params = None
            self._pl_log(f"Run Job1: NO plate calibration -- launching UN-PARAMETERIZED. The "
                         f"mosaic will image the job's STORED PredefinedPoint, NOT well "
                         f"{well or '(none)'}. Calibrate the plate to enable coordinate drive.")
        # Put the rig on the widefield mosaic OC FIRST. Run Job1 has to work from any
        # starting state -- brightfield, 2-photon, mid-PA-sequence -- and the job does not
        # switch to its own CaptureDefinition OC on its own (see MOSAIC_OC). Dispatched
        # synchronously on a worker so the job launch waits for it; the launch still
        # happens if the OC select fails, with the failure logged loudly, because an
        # unilluminated mosaic is obvious on inspection whereas a silently skipped Run
        # Job1 is not.
        threading.Thread(target=self._pl_run_job1_after_oc, args=(params,),
                         daemon=True).start()

    def _pl_run_job1_after_oc(self, params):
        import configparser, time, spheroid_pipeline as _pl
        cp = configparser.ConfigParser()
        try:
            if _pl.SESSION_INI.exists():
                cp.read(_pl.SESSION_INI)
        except Exception:
            pass
        work = cp.get("paths", "work_dir", fallback="").strip() or self._pl_trigger_dir.get().strip()
        if work:
            wd = Path(work)
            try:
                # state= points at a file that does not exist ON PURPOSE: selecting the OC
                # already sets the light source and its powers, and replaying the snapshot
                # over the top would undo exactly what we just asked for.
                _pl._atomic_write_crlf(wd / "lightpath_trigger.ini", [
                    "[lightpath]", "mode=restore", "dichroic_in=0", "zoom=0.0",
                    f"a1_on={'1' if self._pl_pa_a1on.get() else '0'}",
                    f"oc={self.MOSAIC_OC}",
                    f"state={(wd / '_no_such_state.ini').as_posix()}"])
                (wd / "cmd_done.ini").unlink(missing_ok=True)
                _pl._atomic_write_crlf(wd / "cmd.ini",
                                       ["[command]", "action=mosaic_oc", "action_id=12"])
                self._pl_log(f"Run Job1: selecting mosaic OC '{self.MOSAIC_OC}'")
                done, status = wd / "cmd_done.ini", None
                for _ in range(120):
                    if done.exists():
                        c = configparser.ConfigParser()
                        try:
                            c.read(done)
                            status = c.get("command", "status", fallback="?")
                        except Exception:
                            time.sleep(1.0); continue
                        if status == "?":
                            # cmd_done.ini is built key-by-key with status written LAST, so a
                            # file that parses without one is still being written -- not a
                            # finished run with an unknown result. Breaking here reported "?"
                            # and aborted a post-PA sequence whose capture had actually
                            # succeeded (2026-08-07, cost the 1050 nm pass).
                            time.sleep(1.0); continue
                        break
                    time.sleep(1.0)
                if status != "ok":
                    self._pl_log(f"Run Job1: mosaic OC select -> {status or 'timeout'} "
                                 f"-- the mosaic may come back UNILLUMINATED; check the nd2 "
                                 f"has a MultiLaser(EPI) 555 entry.")
            except Exception as e:
                self._pl_log(f"Run Job1: could not select the mosaic OC ({e}) -- continuing.")
        self._pl_run_job("Step1_Locate_via_scan_WC", False, params=params)
        # Snapshot existing run folders BEFORE the job creates its new one, so the watcher
        # picks only the folder generated by THIS click.
        try:
            base = Path(self.MOSAIC_JOB_DIR)
            before = {str(p) for p in base.iterdir() if p.is_dir()} if base.is_dir() else set()
            self.after(0, lambda: self._pl_screen_lbl.configure(
                text="Mosaic: waiting for the JOB to save its nd2…", fg=YELLOW))
            threading.Thread(target=self._pl_watch_mosaic, args=(base, before), daemon=True).start()
        except Exception as e:
            self._pl_log(f"Run Job1: could not start mosaic auto-load watcher: {e}")

    def _pl_watch_mosaic(self, base, before):
        """Poll the Step1 mosaic JOB output dir for the NEW timestamped folder (not in `before`);
        once it holds a STABLE nd2, set _pl_mosaic_path to the LARGEST nd2 in it (the stitched
        whole-well mosaic). ~15 min, then gives up (the operator can Browse manually)."""
        import time
        last_size = -1
        stable = 0
        for _ in range(300):          # ~15 min at 3 s/poll
            time.sleep(3)
            try:
                if not base.is_dir():
                    continue
                new = [p for p in base.iterdir() if p.is_dir() and str(p) not in before]
                if not new:
                    continue
                folder = max(new, key=lambda p: p.stat().st_mtime)   # newest of the new folders
                nd2s = list(folder.rglob("*.nd2"))
                if not nd2s:
                    continue
                nd2 = max(nd2s, key=lambda p: p.stat().st_size)      # largest = the stitched mosaic
                sz = nd2.stat().st_size
                if sz > 0 and sz == last_size:
                    stable += 1
                else:
                    stable = 0; last_size = sz
                if stable >= 3:        # size unchanged ~9 s -> the write has finished
                    p = str(nd2)
                    self.after(0, lambda pp=p: self._pl_mosaic_path.set(pp))
                    self.after(0, lambda pp=p: self._pl_screen_lbl.configure(
                        text=f"Mosaic auto-loaded: {Path(pp).name} — click Run Screener.", fg=GREEN))
                    self._pl_log(f"Run Job1: mosaic auto-loaded -> {p}")
                    return
            except Exception:
                continue
        self.after(0, lambda: self._pl_screen_lbl.configure(
            text="Mosaic auto-load timed out — Browse to the nd2 manually.", fg=YELLOW))
        self._pl_log("Run Job1: mosaic auto-load watcher timed out (no new nd2) — Browse manually.")

    def _pl_run_job(self, name_var, gated, params=None, preconfirmed=False):
        """Launch a NIS-E JOB (nis_macro_run_job.mac, action 10). Writes job_trigger.ini,
        then dispatches. When gated, requires the A1-powered confirmation AND an explicit
        confirm click (the JOB fires the 850 nm PA activation laser).

        params: optional dict of {"<Task>.<Param>": value} (e.g. a WellSelection well). When
        given, the JOB is launched via Jobs_RunJobInitParam (sets the params AND launches) --
        the GUI writes job_params.json next to the trigger as UTF-16 (the docs REQUIRE the
        JSON be Unicode), and the macro ReadFiles it. When None, a stale job_params.json is
        removed so the macro falls back to plain Jobs_RunJobByName."""
        import configparser, spheroid_pipeline as _pl
        cp = configparser.ConfigParser()
        try:
            if _pl.SESSION_INI.exists():
                cp.read(_pl.SESSION_INI)
        except Exception:
            pass
        work = cp.get("paths", "work_dir", fallback="").strip() or self._pl_trigger_dir.get().strip()
        if not work:
            o = self._pl_out_dir.get().strip()
            work = str(Path(o) / "autofocus") if o else ""
        if not work:
            messagebox.showwarning("No work dir",
                                   "Run Step 3 (Trigger) first, or set the Trigger dir."); return
        # name_var may be a StringVar (Job3 dropdown) OR a plain literal string (Job1,
        # hard-wired to the imaging-only mosaic so it can never point at the PA job).
        name = (name_var.get() if hasattr(name_var, "get") else str(name_var)).strip()
        if not name:
            self._pl_log("Run Job: no job name set -- aborted"); return
        if gated:
            if not self._pl_pa_a1on.get():
                messagebox.showwarning("A1 not confirmed",
                                       "Tick 'A1 powered ON' before firing the PA laser."); return
            # preconfirmed: Run Pipeline already took BOTH gates on the main thread before
            # spawning its worker (Tk modals are not thread-safe, and an unattended run must
            # not stall on a dialog nobody is watching). The A1 check above still applies.
            if preconfirmed:
                self._pl_log(f"Run Job: '{name_var if isinstance(name_var, str) else ''}' "
                             f"firing under a pipeline pre-confirmation.")
            if not preconfirmed and not messagebox.askyesno(
                    "Fire PA laser?",
                    f"This launches the '{name}' JOB, which FIRES the PA activation laser.\n\n"
                    f"Optical config: {self._pl_pa_oc.get().strip() or '(unset)'}\n"
                    f"A1 confirmed ON.\n\n"
                    f"Proceed?"):
                return
        project = self._pl_job_project.get().strip()
        wd = Path(work)
        try:
            wd.mkdir(parents=True, exist_ok=True)
            _pl._atomic_write_crlf(wd / "job_trigger.ini",
                                   ["[job]", f"project={project}", f"name={name}"])
            # Job parameters (e.g. a WellSelection well): write a Unicode/UTF-16 JSON that the
            # macro ReadFiles into Jobs_RunJobInitParam. Remove any stale file when there are
            # no params, so a plain launch uses Jobs_RunJobByName (the macro keys off its
            # presence). UTF-16 (with BOM) is Windows "Unicode" -- required by the docs.
            pj = wd / "job_params.json"
            if params:
                import json as _json
                pj.write_text(_json.dumps(params, ensure_ascii=False), encoding="utf-16")
            else:
                pj.unlink(missing_ok=True)
        except Exception as exc:
            messagebox.showerror("Run Job", f"Failed to write job trigger/params: {exc}"); return
        self._pl_send_command("run_job", 10)
        pmsg = f" with params {list(params)}" if params else ""
        self._pl_log(f"Run Job: launched '{name}' (project {project}){pmsg} via action run_job (id 10)")

    def _pl_send_command(self, action, action_id):
        """Write a dispatcher command (cmd.ini) for the always-on nis_macro_dispatcher.mac
        to RunMacro the matching step macro in-process; poll cmd_done.ini for the result.
        Reuses session.ini's work_dir so the dispatcher and the step daemons agree."""
        if getattr(self, "_pl_dispatch_busy", False):
            messagebox.showinfo("Dispatcher busy",
                                "A macro is already running via the dispatcher — wait for it to finish."); return
        # A leftover abort.ini would make the macro we are about to dispatch bail at
        # rank 1 having captured nothing. Nothing else clears it on this path, and it
        # otherwise persists in the run folder indefinitely.
        self._pl_abort_clear()
        import configparser, spheroid_pipeline as _pl
        cp = configparser.ConfigParser()
        try:
            if _pl.SESSION_INI.exists():
                cp.read(_pl.SESSION_INI)
        except Exception:
            pass
        work = cp.get("paths", "work_dir", fallback="").strip() or self._pl_trigger_dir.get().strip()
        if not work:
            o = self._pl_out_dir.get().strip()
            work = str(Path(o) / "autofocus") if o else ""
        if not work:
            messagebox.showwarning("No work dir",
                                   "Run Step 3 (Trigger) first, or set the Trigger dir."); return
        nd2 = cp.get("paths", "nd2_dir", fallback="").strip()
        macro_dir = _pl.MACRO_DIR.as_posix()
        wd = Path(work)
        try:
            wd.mkdir(parents=True, exist_ok=True)
            _pl.SESSION_INI.parent.mkdir(parents=True, exist_ok=True)
            lines = ["[paths]", f"work_dir={wd.as_posix()}"]
            if nd2:
                lines.append(f"nd2_dir={nd2}")
            lines.append(f"macro_dir={macro_dir}")
            _pl._atomic_write_crlf(_pl.SESSION_INI, lines)
            try:
                (wd / "cmd_done.ini").unlink()
            except FileNotFoundError:
                pass
            _pl._atomic_write_crlf(wd / "cmd.ini",
                                   ["[command]", f"action={action}", f"action_id={action_id}"])
        except Exception as exc:
            messagebox.showerror("Dispatcher", f"Failed to write command: {exc}"); return
        self._pl_log(f"Dispatcher: sent '{action}' (id {action_id}) -> {wd / 'cmd.ini'}")
        self._pl_dispatch_status.configure(
            text=f"Sent '{action}'. If nothing runs, start nis_macro_dispatcher.mac once in NIS-E.",
            fg=YELLOW)
        self._pl_dispatch_busy = True
        threading.Thread(target=self._pl_poll_cmd_done, args=(wd, action), daemon=True).start()

    def _pl_poll_cmd_done(self, work_dir, action):
        import configparser, time
        done = Path(work_dir) / "cmd_done.ini"
        try:
            for _ in range(1800):      # up to ~30 min (a step macro may loop internally)
                if done.exists():
                    cp = configparser.ConfigParser()
                    try:
                        cp.read(done)
                        # status normally under [command]; tolerate an older dispatcher
                        # build that wrote it under [spheroid] so the GUI still shows ok.
                        status = cp.get("command", "status",
                                        fallback=cp.get("spheroid", "status", fallback="?"))
                        message = cp.get("command", "message",
                                         fallback=cp.get("spheroid", "message", fallback=""))
                    except Exception:
                        time.sleep(1.0); continue
                    if status == "?":
                        # Still being written -- the macro sets status LAST. Reporting "?"
                        # here would show the operator a failed-looking result for a run
                        # that is merely a fraction of a second from succeeding.
                        time.sleep(1.0); continue
                    fg = GREEN if status == "ok" else RED
                    self.after(0, lambda s=status, a=action: self._pl_dispatch_status.configure(
                        text=f"Dispatcher: '{a}' -> {s}", fg=fg))
                    self._pl_log(f"Dispatcher: '{action}' completed -> {status}")
                    if message:
                        self._pl_log(f"   ⚠ {message}")
                    return
                time.sleep(1.0)
            self.after(0, lambda a=action: self._pl_dispatch_status.configure(
                text=f"Dispatcher: '{a}' — no cmd_done.ini yet (is nis_macro_dispatcher.mac running?)",
                fg=YELLOW))
        finally:
            self._pl_dispatch_busy = False

    # _pl_init_dichroic_toggle removed 2026-08-08 with the PFS ON checkbox. It existed to
    # force PFS off (and grey the box) whenever Dichroic OUT was ticked, because the focus
    # lock cannot hold through a dichroic-OUT PA. PFS is now unconditionally off on this
    # rig, so there is nothing left to interlock -- see the INIT_PFS_ON note in __init__.
    def _pl_resolve_work_dir(self):
        """The session work dir (where every trigger/handshake file lives), resolved the
        same way everywhere: session.ini first, then the Trigger dir field, then
        <output>/autofocus. Returns a Path or None -- and never shows a dialog, so the
        pipeline worker thread can call it too."""
        import configparser, spheroid_pipeline as _pl
        cp = configparser.ConfigParser()
        try:
            if _pl.SESSION_INI.exists():
                cp.read(_pl.SESSION_INI)
        except Exception:
            pass
        work = cp.get("paths", "work_dir", fallback="").strip() or self._pl_trigger_dir.get().strip()
        if not work:
            o = self._pl_out_dir.get().strip()
            work = str(Path(o) / "autofocus") if o else ""
        return Path(work) if work else None

    def _pl_write_init_trigger(self, wd):
        """Write init_trigger.ini (rig known-good defaults) into `wd`. Returns a short
        description of what was written, or raises. Split out of _pl_init_rig so the Run
        Pipeline worker can write it without touching a Tk dialog.

        power_pct and pfs_on are no longer operator-settable (see the var block): the
        activation power lives in the 850 nm OC, and PFS is always off on this rig. Both
        keys stay in the file because the key set is nis_macro_init_rig.mac's contract --
        Int_GetKeyString on a missing key is not a case the macro handles."""
        import spheroid_pipeline as _pl
        zoom     = float(self._pl_init_zoom.get().strip())
        z_centre = float(self._pl_z_centre.get().strip())
        z_half   = float(self._pl_z_half.get().strip())
        z_step   = float(self._pl_z_step.get().strip())
        power    = min(float(self._pl_init_power.get().strip()), MAX_PA_ACTIVATION_POWER_PCT)
        dic      = '1' if self._pl_init_dichroic_out.get() else '0'
        wd.mkdir(parents=True, exist_ok=True)
        _pl._atomic_write_crlf(wd / "init_trigger.ini", [
            "[init]",
            f"zoom={zoom}",
            f"power_pct={power}",
            f"z_centre={z_centre}",
            f"z_half={z_half}",
            f"z_step={z_step}",
            f"dichroic_out={dic}",
            f"pfs_on={'1' if self._pl_init_pfs.get() else '0'}",
            f"a1_on={'1' if self._pl_pa_a1on.get() else '0'}",
            f"save_repoint={'1' if self._pl_init_save_repoint.get() else '0'}",
        ])
        return (f"zoom={zoom}, z={z_centre}+/-{z_half}@{z_step}, dichroic_out={dic}, "
                f"pfs=off")

    def _pl_init_rig(self):
        """Write init_trigger.ini and dispatch nis_macro_init_rig.mac (action 8) to
        normalize the rig: confocal zoom, z-defaults, dichroic, and re-point Save-to-File
        at the pipeline nd2_dir. z_centre/z_half/z_step come from Step 3; the A1-powered
        guard reuses PA Setup's _pl_pa_a1on. Run Pipeline calls _pl_write_init_trigger
        directly instead of this -- see _pl_pipeline_init_rig."""
        wd = self._pl_resolve_work_dir()
        if wd is None:
            messagebox.showwarning("No work dir",
                                   "Run Step 3 (Trigger) first, or set the Trigger dir."); return
        try:
            desc = self._pl_write_init_trigger(wd)
        except ValueError:
            messagebox.showwarning("Bad number",
                                   "Zoom / z_centre / z_half / z_step must all be numeric."); return
        except Exception as exc:
            messagebox.showerror("Write failed", str(exc)); return
        self._pl_send_command("init_rig", 8)
        self._pl_log(f"Init: rig defaults written to {wd / 'init_trigger.ini'} ({desc}) "
                     f"-> dispatched init_rig (id 8)")

    # _pl_watch_pa_done / _pl_watch_pa_done_thread removed 2026-08-08 with the Job3 card.
    # They polled <work_dir>/pa_done.ini for ~1 h so the operator could see "PA JOB
    # complete" without watching NIS-E. Run Pipeline now waits on the PA itself
    # (_pl_pa_wait_for_job, which watches for *ZStack*.nd2 appearing in the PA folder and
    # then going quiet) and logs the result, so a separate watcher button was a second,
    # weaker answer to a question already answered -- and pa_done.ini is not written at all
    # when the job is launched through the dispatcher.

    def _btn(self, parent, text, command, bg, fg, bold=False, side="left", grid=None):
        b = tk.Button(parent, text=text, command=command,
                      bg=bg, fg=fg, relief="flat", cursor="hand2",
                      font=("Segoe UI", 9, "bold" if bold else "normal"),
                      padx=10, pady=4, activebackground=SURFACE2,
                      activeforeground=TEXT)
        if grid is not None:
            b.grid(row=grid[0], column=grid[1], padx=(4, 0), pady=2)
        elif side:
            b.pack(side=side, padx=(0, 6))
        else:
            b.pack(padx=(0, 6))
        return b


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    app = App()
    app.mainloop()
