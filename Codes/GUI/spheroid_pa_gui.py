"""
spheroid_pa_gui.py  v1.13.0
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

# NIS-E "mGold" optical-config group. Exact names as configured on the rig (Optical
# Configuration list) -- see NISE010/NISE011 for the activation-vs-visualization
# wavelength rationale (IMPA004, SLIM025/026/031/043).
#
# Activation range (~750-850 nm): IMPA004's 10-nm wavelength sweep to find the best
# 2P photoconversion wavelength per PA-protein/dye; "_PA"/"_PA1"/"_PA2" suffixed
# configs are purpose-built for firing. 850 nm/30% (mGold/PAmKate) is the current
# pipeline default; SLIM031 activated PA-JF646 at 800 nm instead -- the optimum is
# protein/dye-specific, hence the full sweep is offered rather than one fixed value.
PA_ACTIVATION_OC_LIST = [
    "850 nm power loop full reso2",     # default: mGold/PAmKate 2P activation (IMPA004)
    "750 nm power loop full reso",
    "760 nm power loop full reso1",
    "770 nm power loop full reso2",
    "780 nm power loop full reso1",
    "790 nm power loop full reso2",
    "800 nm power loop full reso1",
    "810 nm power loop full reso2",
    "820 nm power loop full reso1",
    "830 nm power loop full reso2",
    "840 nm power loop full reso1",
    "850 nm power loop full galvo_BT",
    "750nm_Galvo_405nm_NDD2_PA1",
    "750nm_Galvo_405nm_NDD2_PA_band",
    "760nm_Galvo_405nm_NDD2_PA1",
    "800nm_Galvo_405nm_NDD2_PA2",
    "800nm_Galvo_488nm_NDD2_PA",
    "800nm_Reso_405nm_NDD2_PA2",
    "880nm_Galvo_405nm_NDD2_PA1",
]

# Visualization range (~890-1050 nm): sits above the activation band, so imaging here
# does not trigger photoconversion. Each config targets a specific marker/channel --
# 890 nm = mBeRFP (constitutive T-cell identity marker, SLIM025/026/043); 940 nm =
# PAsfGFP (the photoactivatable T-cell readout, imaged before AND after PA to confirm
# conversion); 1050 nm = spheroid-depth / faded-square re-imaging (SLIM031, pa_validate
# default). 950/970/980/1000 nm are Galvo/Resonance alternates for the same channels.
PA_VIZ_OC_LIST = [
    "1050nm_Galvo_561nm_NDD2_JL2",       # default: spheroid depth / PA-validate re-image
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
        self.title("SpheroidPA  v1.13.0 — NIS-E Spheroid PA Pipeline")
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
        tk.Label(hdr, text="  SpheroidPA  v1.13.0",
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
        # an in-flight Z-series/Capture still needs Esc in NIS-E for a hard stop.
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
        stamped = f"[{datetime.now().strftime('%H:%M:%S')}] {text}"
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
        self._pl_n_spheroids = tk.StringVar(value="")
        self._pl_z_centre    = tk.StringVar(value="7680.0")
        self._pl_z_half      = tk.StringVar(value="90.0")
        self._pl_z_step      = tk.StringVar(value="5.0")
        self._pl_recenter_flipx = tk.BooleanVar(value=False)
        self._pl_recenter_flipy = tk.BooleanVar(value=False)
        self._pl_recenter_n  = tk.StringVar(value="all")
        self._pl_recentered  = {}
        self._pl_locate_only = tk.BooleanVar(value=False)
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
        self._pl_pa_job      = tk.StringVar(value="step3_zstack_PA")
        self._pl_pa_oc       = tk.StringVar(value="850 nm power loop full reso2")
        self._pl_pa_power    = tk.StringVar(value="30")
        self._pl_pa_power.trace_add("write", self._pl_pa_power_clamp)
        # No separate well var: PA Setup's and Job3's Well fields both bind to
        # _pl_well_id (Step 1's "working well"), so setting it once syncs everywhere.
        self._pl_pa_loops    = tk.StringVar(value="70")
        self._pl_pa_zoom     = tk.StringVar(value="8")
        self._pl_pa_dichroic = tk.BooleanVar(value=True)   # True = dichroic OUT
        self._pl_pa_interlock = tk.BooleanVar(value=True)  # True = remove A1 interlock first
        self._pl_pa_a1on     = tk.BooleanVar(value=False)  # True = A1 confirmed powered ON; pa_setup/pa_validate guard aborts if unchecked
        # NIS-E JOB launch (Jobs_RunJobByName via nis_macro_run_job.mac, action 10).
        # jobs_dir points at NIS-E's jobs/project folder so the name dropdowns list the
        # actual JOB names on the rig instead of a free-typed string the operator can mistype.
        self._pl_jobs_dir    = tk.StringVar(value="")
        self._pl_job_project = tk.StringVar(value="IMAGEN")
        self._pl_job1_name   = tk.StringVar(value="Step1_Locate_via_scan")
        self._pl_job3_name   = tk.StringVar(value="step3_zstack_PA")
        # Per-macro "include in Run Pipeline" selections.
        self._pl_pa_sel_setup    = tk.BooleanVar(value=True)
        self._pl_pa_sel_points   = tk.BooleanVar(value=True)
        self._pl_dispatch_busy   = False   # one dispatcher command in flight at a time

        # Pre-PA card: baseline viz capture(s) before PA Setup, via the SAME
        # dispatcher z-stack action (id 2) -- one pass per checked OC, each routed
        # into its own nd2/prePA_<tag> subfolder. See SLIM025/026/031/043.
        self._pl_pa_sel_prepa   = tk.BooleanVar(value=False)   # include in Run Pipeline
        self._pl_prepa_oc_890   = tk.BooleanVar(value=False)   # mBeRFP (T-cell identity)
        self._pl_prepa_oc_940   = tk.BooleanVar(value=True)    # PAsfGFP (before readout)
        self._pl_prepa_oc_1050  = tk.BooleanVar(value=False)   # spheroid depth / faded-square
        self._pl_prepa_locate_only = tk.BooleanVar(value=False)  # z_half=0 -> 1-plane (mirrors Step 3)
        # After-PA Validation: independent channel selection, saves to postPA_<tag>/
        self._pl_postpa_oc_890  = tk.BooleanVar(value=False)
        self._pl_postpa_oc_940  = tk.BooleanVar(value=True)
        self._pl_postpa_oc_1050 = tk.BooleanVar(value=False)
        self._pl_postpa_locate_only = tk.BooleanVar(value=False)
        # Keep each captured ND2 open in NIS-E after it is saved. OFF by default: a
        # multi-spheroid x multi-OC x z-stack run opens dozens of image windows, which
        # clutters NIS-E and can stall the acquisition itself. The GUI's "Captured
        # Z-Stacks" tab is where captures get reviewed, so there is no need to also
        # hold them all open on the rig. Tick only for live debugging.
        self._pl_keep_open_in_nise = tk.BooleanVar(value=False)

        # Initialization card (Step 4, first card): rig known-good defaults dispatched
        # as init_rig (action 8) -> init_trigger.ini. z_centre/z_half/z_step reuse the
        # Step 3 vars and a1_on reuses _pl_pa_a1on. Standalone (not in Run Pipeline).
        self._pl_init_zoom         = tk.StringVar(value="2.0")
        self._pl_init_power        = tk.StringVar(value="30")
        self._pl_init_dichroic_out = tk.BooleanVar(value=False)   # unchecked = dichroic IN
        self._pl_init_pfs          = tk.BooleanVar(value=True)
        self._pl_init_save_repoint = tk.BooleanVar(value=True)

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
        job1 = tk.LabelFrame(f1, text=" Job1 (X10 Mosaic -- Step1_Locate_via_scan) ",
                             bg=BG, fg=MAUVE, font=("Segoe UI", 9, "bold"),
                             bd=1, relief="groove")
        job1.pack(fill="x", pady=(0, 4))
        self._nd2_file_row(job1, "10X mosaic ND2:", self._pl_mosaic_path)
        row_w = tk.Frame(job1, bg=BG); row_w.pack(fill="x", pady=(2, 4), padx=4)
        tk.Label(row_w, text="Well ID (working well -- shared by all jobs):", width=34, anchor="w",
                 bg=BG, fg=TEXT2, font=("Segoe UI", 9)).pack(side="left")
        tk.Entry(row_w, textvariable=self._pl_well_id, width=14,
                 bg=SURFACE, fg=TEXT, insertbackground=TEXT, relief="flat",
                 font=("Segoe UI", 9)).pack(side="left")

        # NIS-E JOB launch: point at the rig's jobs/project folder so the name dropdown
        # lists real JOB names (Jobs_RunJobByName), then fire Job1 (the 10X mosaic) directly.
        row_j = tk.Frame(job1, bg=BG); row_j.pack(fill="x", pady=(2, 2), padx=4)
        tk.Label(row_j, text="NIS-E Jobs dir:", bg=BG, fg=TEXT2,
                 font=("Segoe UI", 9)).pack(side="left", padx=(0, 3))
        tk.Entry(row_j, textvariable=self._pl_jobs_dir, width=28, bg=SURFACE, fg=TEXT,
                 insertbackground=TEXT, relief="flat", font=("Segoe UI", 9)).pack(side="left")
        self._btn(row_j, "Browse",
                  lambda: self._pl_browse_dir(self._pl_jobs_dir), SURFACE2, TEXT, side="left")
        tk.Label(row_j, text="Project:", bg=BG, fg=TEXT2,
                 font=("Segoe UI", 9)).pack(side="left", padx=(6, 3))
        tk.Entry(row_j, textvariable=self._pl_job_project, width=10, bg=SURFACE, fg=TEXT,
                 insertbackground=TEXT, relief="flat", font=("Segoe UI", 9)).pack(side="left")
        row_j2 = tk.Frame(job1, bg=BG); row_j2.pack(fill="x", pady=(0, 4), padx=4)
        self._btn(row_j2, "Run Job1 (10X mosaic)",
                  lambda: self._pl_run_job(self._pl_job1_name, False), BLUE, "#1e1e2e", side="left")
        _cb_job1 = ttk.Combobox(row_j2, textvariable=self._pl_job1_name, width=26,
                                font=("Segoe UI", 9), values=self._pl_list_jobs())
        _cb_job1.configure(postcommand=lambda cb=_cb_job1: cb.configure(values=self._pl_list_jobs()))
        _cb_job1.pack(side="left", padx=(6, 0))

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
        self._btn(btn3a, "Trigger NIS-E Z-Stack Captures",
                  self._pl_trigger_autofocus, BLUE, "#1e1e2e")
        tk.Checkbutton(btn3a, text="Locate only (1 plane)", variable=self._pl_locate_only,
                       bg=BG, fg=TEXT2, selectcolor=SURFACE, activebackground=BG,
                       activeforeground=TEXT, font=("Segoe UI", 8)).pack(side="left", padx=(8, 0))
        self._pl_af_status_lbl = tk.Label(f_s3, text="Autofocus: not started.",
                                           bg=BG, fg=SUBTEXT, font=("Segoe UI", 9), anchor="w",
                                           justify="left", wraplength=480)
        self._pl_af_status_lbl.pack(fill="x", padx=12)

        btn3b = tk.Frame(f_s3, bg=BG); btn3b.pack(fill="x", padx=12, pady=(6, 0))
        self._btn(btn3b, "Refresh Status", self._pl_poll_af_thread, SURFACE2, TEXT)
        self._btn(btn3b, "Re-center from captures", self._pl_recenter_thread, MAUVE, "#1e1e2e")
        tk.Label(btn3b, text="  # to use:", bg=BG, fg=TEXT2, font=("Segoe UI", 8)).pack(side="left")
        tk.Entry(btn3b, textvariable=self._pl_recenter_n, width=5, bg=SURFACE, fg=TEXT,
                 insertbackground=TEXT, relief="flat", font=("Segoe UI", 9)).pack(side="left", padx=(2, 6))
        tk.Checkbutton(btn3b, text="flip X", variable=self._pl_recenter_flipx, bg=BG, fg=TEXT2,
                       selectcolor=SURFACE, activebackground=BG, activeforeground=TEXT,
                       font=("Segoe UI", 8)).pack(side="left", padx=(8, 0))
        tk.Checkbutton(btn3b, text="flip Y", variable=self._pl_recenter_flipy, bg=BG, fg=TEXT2,
                       selectcolor=SURFACE, activebackground=BG, activeforeground=TEXT,
                       font=("Segoe UI", 8)).pack(side="left", padx=(4, 0))

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
        # Scrollable container: the Step-4 cards (Validation, PA Setup, Job3, PA Points,
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
        tk.Label(pa, text="Start nis_macro_dispatcher.mac once in NIS-E. Per card: edit params, "
                          "Reload writes pa_trigger.ini, Run dispatches that macro. Tick a card to "
                          "include it in Run Pipeline (checked macros run in order, waiting for each).\n"
                          "Run Pipeline covers Validation -> Setup -> Points; it does NOT fire the PA "
                          "itself -- after Points, manually Run step3_zstack_PA in NIS-E JOBS Explorer, "
                          "then run Validation again (same card) to capture the after-PA images. OC "
                          "dropdowns/checkboxes list the full mGold profile set: Activation OC = "
                          "750-850nm (per-protein/dye 2P activation sweep, IMPA004); Viz OC = "
                          "890-1050nm (890=mBeRFP T-cell identity, 940=PAsfGFP before/after, "
                          "1050=spheroid depth/faded-square re-image).",
                 bg=BG, fg=SUBTEXT, font=("Segoe UI", 8), justify="left",
                 wraplength=520).pack(anchor="w", padx=8, pady=(3, 4))

        # Card -1 -- Initialization (dispatcher init_rig action, id 8): normalize the rig
        # + GUI to known-good defaults before an experiment (confocal zoom=2, PA power
        # cap, z-defaults, dichroic IN, PFS ON, Save-to-File -> pipeline nd2_dir). Writes
        # init_trigger.ini and dispatches nis_macro_init_rig.mac. Standalone card (not in
        # the Run Pipeline sequence). First card in Step 4. z_centre/z_half/z_step reuse
        # the Step 3 vars; A1-powered guard reuses PA Setup's _pl_pa_a1on.
        init = tk.Frame(pa, bg=BG2, bd=1, relief="solid"); init.pack(fill="x", padx=4, pady=3)
        hdri = tk.Frame(init, bg=BG2); hdri.pack(fill="x", padx=6, pady=(3, 0))
        tk.Label(hdri, text="Initialization  (rig known-good defaults)", bg=BG2, fg=MAUVE,
                 font=("Segoe UI", 9, "bold")).pack(side="left")
        bodyi = tk.Frame(init, bg=BG2); bodyi.pack(fill="x", padx=24, pady=(0, 4))
        ri = tk.Frame(bodyi, bg=BG2); ri.pack(fill="x", pady=(2, 0))
        for lbl, var, w in [("Confocal zoom:", self._pl_init_zoom, 5),
                            ("PA power %:", self._pl_init_power, 5)]:
            tk.Label(ri, text=lbl, bg=BG2, fg=TEXT2, font=("Segoe UI", 9)).pack(side="left", padx=(0, 3))
            tk.Entry(ri, textvariable=var, width=w, bg=SURFACE, fg=TEXT, insertbackground=TEXT,
                     relief="flat", font=("Segoe UI", 9)).pack(side="left", padx=(0, 8))
        ri2 = tk.Frame(bodyi, bg=BG2); ri2.pack(fill="x", pady=(2, 0))
        for lbl, var, w in [("z_centre:", self._pl_z_centre, 8), ("z_half:", self._pl_z_half, 6),
                            ("z_step:", self._pl_z_step, 6)]:
            tk.Label(ri2, text=lbl, bg=BG2, fg=TEXT2, font=("Segoe UI", 9)).pack(side="left", padx=(0, 3))
            tk.Entry(ri2, textvariable=var, width=w, bg=SURFACE, fg=TEXT, insertbackground=TEXT,
                     relief="flat", font=("Segoe UI", 9)).pack(side="left", padx=(0, 8))
        ri3 = tk.Frame(bodyi, bg=BG2); ri3.pack(fill="x", pady=(2, 0))
        tk.Checkbutton(ri3, text="Dichroic OUT", variable=self._pl_init_dichroic_out, bg=BG2, fg=TEXT2,
                       selectcolor=SURFACE, activebackground=BG2, activeforeground=TEXT,
                       font=("Segoe UI", 9), command=self._pl_init_dichroic_toggle).pack(side="left")
        self._pl_init_pfs_chk = tk.Checkbutton(ri3, text="PFS ON", variable=self._pl_init_pfs, bg=BG2, fg=TEXT2,
                       selectcolor=SURFACE, activebackground=BG2, activeforeground=TEXT,
                       font=("Segoe UI", 9))
        self._pl_init_pfs_chk.pack(side="left", padx=(10, 0))
        self._pl_init_dichroic_toggle()
        tk.Checkbutton(ri3, text="Re-point Save-to-File to pipeline nd2_dir",
                       variable=self._pl_init_save_repoint, bg=BG2, fg=TEXT2, selectcolor=SURFACE,
                       activebackground=BG2, activeforeground=TEXT,
                       font=("Segoe UI", 9)).pack(side="left", padx=(10, 0))
        ri4 = tk.Frame(bodyi, bg=BG2); ri4.pack(fill="x", pady=(2, 0))
        tk.Checkbutton(ri4, text="A1 powered ON", variable=self._pl_pa_a1on, bg=BG2, fg="#f7768e",
                       selectcolor=SURFACE, activebackground=BG2, activeforeground=TEXT,
                       font=("Segoe UI", 9, "bold")).pack(side="left")
        btn_init = tk.Frame(init, bg=BG2); btn_init.pack(fill="x", padx=24, pady=(2, 4))
        self._btn(btn_init, "Initialize Rig", self._pl_init_rig, BLUE, "#1e1e2e", side="left")
        tk.Label(init, text="Sets zoom=2 / power cap / save-path / z-defaults / dichroic / PFS to "
                            "known-good before an experiment; detector gains and per-OC re-save are "
                            "rig-side (handled in the macro / rig).",
                 bg=BG2, fg=SUBTEXT, font=("Segoe UI", 8), justify="left",
                 wraplength=460).pack(anchor="w", padx=24, pady=(0, 6))

        # Card 0 -- Validation (dispatcher z-stack action, id 2): before- AND after-PA
        # multi-channel capture, one pass per checked OC, each routed into its own
        # nd2/prePA_<tag> subfolder (SLIM025/026/031/043 before/after protocol). Run it
        # once before PA Setup (baseline) and again after the JOB (post-PA) -- same
        # card, same mechanism; replaces the old single-channel 1050nm-only PA Validate
        # card (removed 2026-07-07). Internal names keep the "_pl_prepa_*" prefix.
        card, body = self._pl_pa_card(pa, "Before-PA Validation  (baseline multi-channel capture -> nd2/prePA_<tag>)",
                                      self._pl_pa_sel_prepa)
        tk.Label(body, text="Channel selection only -- position/Z come from whatever "
                            "af_trigger_NN.ini files already exist (generate them in Step 3, or "
                            "write them by hand). Each checked OC injects oc= into every existing "
                            "trigger and runs one full z-stack pass via the dispatcher "
                            "(nis_macro_capture_zstack.mac, action 2) into its own "
                            "nd2/prePA_<tag> subfolder. Run this before PA Setup for the baseline; the "
                            "After-PA Validation card below captures the post-JOB result.",
                 bg=BG2, fg=SUBTEXT, font=("Segoe UI", 8), justify="left",
                 wraplength=460).pack(anchor="w")
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

        # Card 1 -- PA Setup (action 4): hardware prep + centred activation square.
        card, body = self._pl_pa_card(pa, "PA Setup  (interlock / OC / dichroic / centred square)",
                                      self._pl_pa_sel_setup)
        r = tk.Frame(body, bg=BG2); r.pack(fill="x")
        tk.Label(r, text="Activation OC:", bg=BG2, fg=TEXT2, font=("Segoe UI", 9)).pack(side="left", padx=(0, 3))
        ttk.Combobox(r, textvariable=self._pl_pa_oc, width=26, font=("Segoe UI", 9),
                     values=PA_ACTIVATION_OC_LIST).pack(side="left")
        r = tk.Frame(body, bg=BG2); r.pack(fill="x", pady=(2, 0))
        for lbl, var, w in [("Power %:", self._pl_pa_power, 5), ("Zoom:", self._pl_pa_zoom, 5),
                            ("Well:", self._pl_well_id, 5), ("Loops:", self._pl_pa_loops, 5)]:
            tk.Label(r, text=lbl, bg=BG2, fg=TEXT2, font=("Segoe UI", 9)).pack(side="left", padx=(0, 3))
            tk.Entry(r, textvariable=var, width=w, bg=SURFACE, fg=TEXT, insertbackground=TEXT,
                     relief="flat", font=("Segoe UI", 9)).pack(side="left", padx=(0, 8))
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
        self._pl_pa_card_buttons(card, "pa_setup", 4)

        # Job3 settings (step3_zstack_PA -- the manual JOB itself, not a dispatcher
        # macro; no pipeline checkbox / Reload / Run here -- nothing to dispatch).
        # Well is bound to the SAME _pl_well_id as Job1's settings, so it's always
        # the current working well -- set it once in Step 1, syncs everywhere.
        job3 = tk.Frame(pa, bg=BG2, bd=1, relief="solid")
        job3.pack(fill="x", padx=4, pady=3)
        hdr3 = tk.Frame(job3, bg=BG2); hdr3.pack(fill="x", padx=6, pady=(3, 0))
        tk.Label(hdr3, text="Job3 (step3_zstack_PA -- manual JOB run)", bg=BG2, fg=MAUVE,
                 font=("Segoe UI", 9, "bold")).pack(side="left")
        body3 = tk.Frame(job3, bg=BG2); body3.pack(fill="x", padx=24, pady=(0, 6))
        tk.Label(body3, text="Run this JOB manually in NIS-E JOBS Explorer after PA Points. "
                             "Well below is the SAME working well as Job1 -- set it once in "
                             "Step 1 and every job's well field updates together.",
                 bg=BG2, fg=SUBTEXT, font=("Segoe UI", 8), justify="left",
                 wraplength=460).pack(anchor="w")
        r = tk.Frame(body3, bg=BG2); r.pack(fill="x", pady=(2, 0))
        tk.Label(r, text="Well:", bg=BG2, fg=TEXT2, font=("Segoe UI", 9)).pack(side="left", padx=(0, 3))
        tk.Entry(r, textvariable=self._pl_well_id, width=8, bg=SURFACE, fg=TEXT,
                 insertbackground=TEXT, relief="flat", font=("Segoe UI", 9)).pack(side="left")
        # The JOB runs manually; a companion macro (nis_macro_pa_done.mac / the JOB's
        # post-capture hook) writes pa_done.ini into the work dir when it finishes. This
        # button polls for that flag so the operator gets a timestamped "PA JOB complete"
        # in the log instead of watching NIS-E.
        rw = tk.Frame(body3, bg=BG2); rw.pack(fill="x", pady=(4, 0))
        self._btn(rw, "Watch for PA done", self._pl_watch_pa_done_thread, GREEN, "#1e1e2e", side="left")
        self._pl_pa_done_status = tk.Label(rw, text="PA done: not watching", bg=BG2, fg=SUBTEXT,
                                           font=("Segoe UI", 8))
        self._pl_pa_done_status.pack(side="left", padx=(8, 0))

        # Card 2 -- PA Points (action 5): build the ND multipoint from ALL active
        # Step 3 triggers -- no count of its own; Step 3's checked-row table IS the
        # selection (its "Trigger" click clears stale triggers and writes fresh ones
        # for exactly the checked rows each time).
        card, body = self._pl_pa_card(pa, "PA Points  (build ND multipoint from ALL active Step 3 triggers)",
                                      self._pl_pa_sel_points)
        r = tk.Frame(body, bg=BG2); r.pack(fill="x")
        tk.Label(r, text="Uses whatever af_trigger_*.ini Step 3 last wrote for its checked rows "
                         "-- check/uncheck spheroids and click Step 3 Trigger again to change the set.",
                 bg=BG2, fg=TEXT2, font=("Segoe UI", 9), wraplength=520, justify="left").pack(side="left")
        self._pl_pa_card_buttons(card, "pa_points", 5)

        # Photoactivation JOB launch (Jobs_RunJobByName via nis_macro_run_job.mac,
        # action 10): fires the 850 nm PA activation laser, so it is gated on the A1
        # power confirmation AND an explicit confirm click. Standalone card (not in the
        # Run Pipeline sequence, which deliberately stops before the laser fires).
        paj = tk.Frame(pa, bg=BG2, bd=1, relief="solid"); paj.pack(fill="x", padx=4, pady=3)
        hdrj = tk.Frame(paj, bg=BG2); hdrj.pack(fill="x", padx=6, pady=(3, 0))
        tk.Label(hdrj, text="Photoactivation JOB  (launch step3_zstack_PA - FIRES 850 nm LASER)",
                 bg=BG2, fg=MAUVE, font=("Segoe UI", 9, "bold")).pack(side="left")
        bodyj = tk.Frame(paj, bg=BG2); bodyj.pack(fill="x", padx=24, pady=(0, 6))
        tk.Label(bodyj, text="Fires the PA activation laser. Requires 'A1 powered ON' + a confirm "
                             "click; run PA Points first so the ND multipoint is loaded.",
                 bg=BG2, fg=RED, font=("Segoe UI", 8), justify="left",
                 wraplength=460).pack(anchor="w")
        rj = tk.Frame(bodyj, bg=BG2); rj.pack(fill="x", pady=(2, 0))
        tk.Label(rj, text="Job:", bg=BG2, fg=TEXT2, font=("Segoe UI", 9)).pack(side="left", padx=(0, 3))
        _cb_job3 = ttk.Combobox(rj, textvariable=self._pl_job3_name, width=26,
                                font=("Segoe UI", 9), values=self._pl_list_jobs())
        _cb_job3.configure(postcommand=lambda cb=_cb_job3: cb.configure(values=self._pl_list_jobs()))
        _cb_job3.pack(side="left")
        rjb = tk.Frame(bodyj, bg=BG2); rjb.pack(fill="x", pady=(4, 0))
        self._btn(rjb, "Run Photoactivation",
                  lambda: self._pl_run_job(self._pl_job3_name, True), RED, "#1e1e2e", side="left")

        # Card 3 -- After-PA Validation (dispatcher z-stack action, id 2): the SAME
        # multi-channel capture as the before-PA Validation card, but its passes save
        # to nd2/postPA_<tag> (not prePA_<tag>) so the after-PA images no longer collide
        # with the baseline -- auto-produces the Before_pa/After_pa layout the comparison
        # figures expect. Run this AFTER the step3_zstack_PA JOB. Manual card (not in the
        # Run Pipeline sequence, which stops before the laser fires).
        apa = tk.Frame(pa, bg=BG2, bd=1, relief="solid"); apa.pack(fill="x", padx=4, pady=3)
        hdra = tk.Frame(apa, bg=BG2); hdra.pack(fill="x", padx=6, pady=(3, 0))
        tk.Label(hdra, text="After-PA Validation  (post-Job capture -> nd2/postPA_<tag>)",
                 bg=BG2, fg=MAUVE, font=("Segoe UI", 9, "bold")).pack(side="left")
        bodya = tk.Frame(apa, bg=BG2); bodya.pack(fill="x", padx=24, pady=(0, 4))
        tk.Label(bodya, text="Run AFTER the step3_zstack_PA JOB. Same positions/Z as the before-PA "
                             "Validation (from Step 3's checked rows), saved to nd2/postPA_<tag> so "
                             "the before/after sets stay separate.",
                 bg=BG2, fg=SUBTEXT, font=("Segoe UI", 8), justify="left",
                 wraplength=460).pack(anchor="w")
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
            is_centre = (i == center)
            cellbg = SURFACE2 if is_centre else BG2
            edge   = BLUE if is_centre else SURFACE
            cell = tk.Frame(self._pl_zv_strip, bg=cellbg, padx=3, pady=3,
                            highlightthickness=2, highlightbackground=edge,
                            highlightcolor=edge)
            cell.pack(side="left", padx=3, pady=4)
            tk.Label(cell, image=photo, bg=cellbg).pack()
            cap = f"Z {absz:.1f} um\n" + ("focus centre" if is_centre else "")
            tk.Label(cell, text=cap, bg=cellbg,
                     fg=(BLUE if is_centre else TEXT2),
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
        row = tk.Frame(self._pl_zv_strip, bg=BG2)
        row.pack(side="top", fill="x", anchor="w", pady=(2, 8))
        tag = "  (FOV-matched crop)" if st.get("cropped") else ""
        tk.Label(row,
                 text=f"{name}     centre Z {st['basez']:.1f} um     {st['n_planes']} planes{tag}",
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
            is_c = (i == center)
            cbg = SURFACE2 if is_c else BG2
            edge = BLUE if is_c else SURFACE
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
            messagebox.showwarning("No records", "Run Steps 1-2 first."); return
        out = self._pl_out_dir.get().strip()
        if not out:
            messagebox.showwarning("No output dir", "Set output directory in Step 1."); return
        # Writing a fresh trigger set means the operator is starting new work -- drop any
        # leftover abort flag so the macro they run next isn't killed at rank 1.
        self._pl_abort_clear()
        try:
            z_centre = float(self._pl_z_centre.get())
            z_half   = float(self._pl_z_half.get())
            z_step   = float(self._pl_z_step.get())
        except ValueError:
            messagebox.showerror(
                "Bad value", "Middle plane Z, Z half-range, and Z step must be numbers."); return
        locate = self._pl_locate_only.get()
        if locate:
            z_half = 0.0          # single plane at z_centre -- fast locate pass for re-centering
        excl = getattr(self, "_pl_excluded", set())
        recs = [r for r in self._pl_records if r.rank not in excl]
        if not recs:
            messagebox.showwarning(
                "None selected",
                "Every spheroid is unchecked in the Use column — nothing to trigger."); return
        try:
            n = _pl.trigger_autofocus_all(recs, Path(out) / "autofocus",
                                          z_centre=z_centre, z_half=z_half, z_step=z_step)
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
            return
        mode = "1-plane LOCATE" if locate else f"Z-stack +/-{z_half:g}/{z_step:g} um"
        skipped = len(self._pl_records) - len(recs)
        extra = f" ({skipped} unchecked skipped)" if skipped else ""
        self._pl_log(f"Step 3: {n} {mode} triggers written to autofocus/ (centre={z_centre}){extra}")
        self.after(0, lambda m=mode, e=extra: self._pl_af_status_lbl.configure(
            text=f"{n} {m} triggers written{e}. Run nis_macro_capture_zstack.mac in NIS-E.", fg=YELLOW))

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
        self._pl_af_status_lbl.configure(text="Re-alignment in process …", fg=YELLOW)
        threading.Thread(target=self._pl_recenter, daemon=True).start()

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
        self.after(0, lambda n=len(res), d=deltas, w=warn: self._pl_af_status_lbl.configure(
            text=(f"Re-aligned {n} from captures (Δum):  {d}.   Re-Trigger + re-run the macro "
                  f"(toggle flip X/Y if more off-center).{w}"),
            fg=(YELLOW if warn else GREEN)))
        self._pl_log(f"Re-center: adjusted {len(res)} spheroids from their captures"
                     + (f"; {len(skipped)} SKIPPED (see above)." if skipped else "."))

    # ── Step 4 handlers ───────────────────────────────────────────────────────

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

    def _pl_pa_write_trigger(self):
        """Write pa_trigger.ini (all PA-macro params) into the session work dir and
        refresh session.ini (work_dir + macro_dir for the dispatcher). Returns the
        work_dir Path, or None on error (after showing a message)."""
        import configparser, spheroid_pipeline as _pl
        base = self._pl_out_dir.get().strip()
        save_dir = (Path(base) / "pa") if base else None
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
        hard-stop still requires Esc in NIS-E. Say so plainly rather than implying
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
            "interrupted from here — press Esc in NIS-E for an immediate hard stop.")

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
            ocs.append(("1050nm_Galvo_561nm_NDD2_JL2", "1050nm_depth"))
        return ocs

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

    def _pl_prepa_run(self, phase="prePA"):
        try:
            self._pl_prepa_capture_ocs(phase)
        finally:
            self._pl_dispatch_busy = False

    def _pl_pa_run_pipeline_thread(self):
        if getattr(self, "_pl_dispatch_busy", False):
            messagebox.showinfo("Dispatcher busy",
                                "A macro is already running via the dispatcher — wait for it to finish."); return
        self._pl_abort_clear()       # a stale flag must not kill the run we're starting
        self._pl_dispatch_busy = True
        threading.Thread(target=self._pl_pa_run_pipeline, daemon=True).start()

    def _pl_pa_run_pipeline(self):
        """Run the checked PA cards in order, each via the dispatcher, waiting for
        cmd_done.ini between steps. Stops on the first non-ok result. Holds the
        single-in-flight _pl_dispatch_busy lock for the whole sequence."""
        import configparser, time, spheroid_pipeline as _pl
        try:
            if not self._pl_regen_triggers():
                return
            do_prepa = self._pl_pa_sel_prepa.get() and bool(self._pl_prepa_checked_ocs())
            seq = []
            if self._pl_pa_sel_setup.get():    seq.append(("pa_setup", 4))
            if self._pl_pa_sel_points.get():   seq.append(("pa_points", 5))
            if not seq and not do_prepa:
                self.after(0, lambda: self._pl_pa_status.configure(
                    text="Pipeline: no macros checked.", fg=RED)); return
            if do_prepa:
                self._pl_log("PA pipeline: Pre-PA viz capture(s) first")
                if not self._pl_prepa_capture_ocs():
                    return   # status already shown by the helper
            if not seq:
                self.after(0, lambda: self._pl_pa_status.configure(
                    text="Pipeline complete: Pre-PA only.", fg=GREEN)); return
            wd = self._pl_pa_write_trigger()
            if wd is None:
                return
            self._pl_log(f"PA pipeline: {[a for a, _ in seq]}")
            for action, aid in seq:
                if self._pl_abort_requested(wd):
                    self._pl_log(f"PA pipeline: ABORT requested -- stopping before '{action}'")
                    self.after(0, lambda: self._pl_pa_status.configure(
                        text="Pipeline: aborted by operator.", fg=RED)); return
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
                status = None
                for _ in range(1800):
                    if done.exists():
                        cp = configparser.ConfigParser()
                        try:
                            cp.read(done); status = cp.get("command", "status",
                                fallback=cp.get("spheroid", "status", fallback="?"))
                        except Exception:
                            time.sleep(1.0); continue
                        break
                    if self._pl_abort_requested(wd):
                        status = "aborted"
                        break
                    time.sleep(1.0)
                self._pl_log(f"PA pipeline: '{action}' -> {status or 'timeout'}")
                if status != "ok":
                    msg = (f"Pipeline stopped at '{action}' -> {status}" if status
                           else f"Pipeline: '{action}' timed out (is nis_macro_dispatcher.mac running?)")
                    self.after(0, lambda m=msg: self._pl_pa_status.configure(text=m, fg=RED)); return
            self.after(0, lambda n=len(seq): self._pl_pa_status.configure(
                text=f"Pipeline complete: {n} macro(s) ok.", fg=GREEN))
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
        # recorded z_centre (e.g. Refresh Status wasn't run) -- the spheroids are all
        # captured centred on that plane, so it IS the correct bin z-centre.
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
                             "(set Step-3 Middle plane Z or run Refresh Status)"); continue
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

    def _pl_list_jobs(self) -> list:
        """Enumerate JOB names from the configured NIS-E jobs dir: immediate subdirectory
        names + file stems found in it, sorted and de-duplicated. Falls back to the two
        pipeline defaults if the dir is unset/missing or on any error."""
        defaults = ["Step1_Locate_via_scan", "step3_zstack_PA"]
        try:
            d = self._pl_jobs_dir.get().strip()
            if not d:
                return defaults
            p = Path(d)
            if not p.is_dir():
                return defaults
            names = set()
            for entry in p.iterdir():
                names.add(entry.name if entry.is_dir() else entry.stem)
            return sorted(names)
        except Exception:
            return defaults

    def _pl_run_job(self, name_var, gated):
        """Launch a NIS-E JOB by name via Jobs_RunJobByName (nis_macro_run_job.mac,
        action 10). Resolves work_dir like _pl_send_command, writes job_trigger.ini there,
        then dispatches. When gated, requires the A1-powered confirmation AND an explicit
        confirm click because the JOB fires the 850 nm PA activation laser."""
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
        name = name_var.get().strip()
        if not name:
            self._pl_log("Run Job: no job name set -- aborted"); return
        if gated:
            if not self._pl_pa_a1on.get():
                messagebox.showwarning("A1 not confirmed",
                                       "Tick 'A1 powered ON' before firing the PA laser."); return
            if not messagebox.askyesno(
                    "Fire PA laser?",
                    f"This launches the '{name}' JOB, which FIRES THE 850 nm PA activation laser.\n\n"
                    f"A1 confirmed ON. Make sure PA Points (the ND multipoint import) has run.\n\n"
                    f"Proceed?"):
                return
        project = self._pl_job_project.get().strip()
        wd = Path(work)
        try:
            wd.mkdir(parents=True, exist_ok=True)
            _pl._atomic_write_crlf(wd / "job_trigger.ini",
                                   ["[job]", f"project={project}", f"name={name}"])
        except Exception as exc:
            messagebox.showerror("Run Job", f"Failed to write job trigger: {exc}"); return
        self._pl_send_command("run_job", 10)
        self._pl_log(f"Run Job: launched '{name}' (project {project}) via action run_job (id 10)")

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

    def _pl_init_dichroic_toggle(self):
        """Dichroic OUT is the higher-priority parameter: the PFS focus lock cannot
        engage with the dichroic mirror extracted, so selecting OUT forces PFS
        unlocked (off) and greys out the PFS checkbox; IN re-enables it. The macro
        enforces the same override on the rig side as a safety."""
        if self._pl_init_dichroic_out.get():
            self._pl_init_pfs.set(False)
            self._pl_init_pfs_chk.configure(state="disabled")
        else:
            self._pl_init_pfs_chk.configure(state="normal")

    def _pl_init_rig(self):
        """Write init_trigger.ini (rig known-good defaults) into the session work dir and
        dispatch nis_macro_init_rig.mac (action 8) to normalize the rig before an
        experiment: confocal zoom, PA power cap, z-defaults, dichroic, PFS, and re-point
        Save-to-File at the pipeline nd2_dir. z_centre/z_half/z_step come from Step 3; the
        A1-powered guard reuses PA Setup's _pl_pa_a1on."""
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
        try:
            zoom     = float(self._pl_init_zoom.get().strip())
            power    = float(self._pl_init_power.get().strip())
            z_centre = float(self._pl_z_centre.get().strip())
            z_half   = float(self._pl_z_half.get().strip())
            z_step   = float(self._pl_z_step.get().strip())
        except ValueError:
            messagebox.showwarning("Bad number",
                                   "Zoom / power / z_centre / z_half / z_step must all be numeric."); return
        if power > MAX_PA_ACTIVATION_POWER_PCT:
            power = MAX_PA_ACTIVATION_POWER_PCT
            self._pl_init_power.set(str(MAX_PA_ACTIVATION_POWER_PCT))
        wd = Path(work)
        init_lines = [
            "[init]",
            f"zoom={zoom}",
            f"power_pct={power}",
            f"z_centre={z_centre}",
            f"z_half={z_half}",
            f"z_step={z_step}",
            f"dichroic_out={'1' if self._pl_init_dichroic_out.get() else '0'}",
            f"pfs_on={'1' if self._pl_init_pfs.get() else '0'}",
            f"a1_on={'1' if self._pl_pa_a1on.get() else '0'}",
            f"save_repoint={'1' if self._pl_init_save_repoint.get() else '0'}",
        ]
        try:
            wd.mkdir(parents=True, exist_ok=True)
            _pl._atomic_write_crlf(wd / "init_trigger.ini", init_lines)
        except Exception as exc:
            messagebox.showerror("Write failed", str(exc)); return
        self._pl_send_command("init_rig", 8)
        self._pl_log(f"Init: rig defaults written to {wd / 'init_trigger.ini'} "
                     f"(zoom={zoom}, power={power}%, z={z_centre}+/-{z_half}@{z_step}, "
                     f"dichroic_out={'1' if self._pl_init_dichroic_out.get() else '0'}, "
                     f"pfs={'1' if self._pl_init_pfs.get() else '0'}) -> dispatched init_rig (id 8)")

    def _pl_watch_pa_done_thread(self):
        threading.Thread(target=self._pl_watch_pa_done, daemon=True).start()

    def _pl_watch_pa_done(self):
        """Poll <work_dir>/pa_done.ini (written by nis_macro_pa_done.mac / the JOB
        post-capture hook when the manual step3_zstack_PA JOB finishes) for ~1 hour.
        When it appears, log a timestamped completion line, flip the status label GREEN,
        and delete the flag so the next run starts clean."""
        import configparser, time, spheroid_pipeline as _pl
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
            self.after(0, lambda: self._pl_pa_done_status.configure(
                text="PA done: no work dir (run Step 3 Trigger first)", fg=YELLOW))
            return
        flag = Path(work) / "pa_done.ini"
        self.after(0, lambda: self._pl_pa_done_status.configure(text="PA done: watching...", fg=YELLOW))
        for _ in range(1800):      # ~1 hour at 2 s per poll
            if flag.exists():
                ts = time.strftime('%H:%M:%S')
                self.after(0, lambda t=ts: self._pl_log(
                    f"[{t}] PA JOB complete - run After-PA Validation"))
                self.after(0, lambda: self._pl_pa_done_status.configure(
                    text="PA done: complete - run After-PA Validation", fg=GREEN))
                try:
                    flag.unlink()
                except FileNotFoundError:
                    pass
                return
            time.sleep(2)
        self.after(0, lambda: self._pl_pa_done_status.configure(
            text="PA done: timed out (no pa_done.ini after ~1 h)", fg=YELLOW))

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
