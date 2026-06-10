"""
spheroid_pa_gui.py  v1.2
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
        self.title("SpheroidPA  v1.2 — NIS-E Spheroid PA Pipeline")
        self.geometry("980x740")
        self.minsize(800, 600)
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
        tk.Label(hdr, text="  SpheroidPA  v1.2",
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
            "pipe": ("  PA Workflow  ", tk.Frame(nb, bg=BG)),
            "log":  ("  Log  ",         tk.Frame(nb, bg=BG)),
        }
        for key, (label, frame) in tabs.items():
            nb.add(frame, text=label)
            setattr(self, f"_tab_{key}", frame)

        self._build_pipeline_tab(self._tab_pipe)
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
        def _ins():
            self._log.configure(state="normal")
            self._log.insert("end", text + "\n", tag)
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
        self._pl_z_half      = tk.StringVar(value="18.0")
        self._pl_z_step      = tk.StringVar(value="2.0")
        self._pl_z_rank      = tk.StringVar()
        self._pl_z_nd2_path  = tk.StringVar()
        self._pl_trigger_dir = tk.StringVar(value=r"S:\Images\Weihao\NISeA\NIS-E-Automation\work")
        self._pl_nd2_out_dir = tk.StringVar(value=r"S:\Images\Weihao\NISeA\NIS-E-Automation\work\nd2")
        self._pl_ch_field    = tk.StringVar(value="CH2LaserPower")
        self._pl_P0          = tk.StringVar(value="15.0")
        self._pl_L_um        = tk.StringVar(value="165.0")
        self._pl_ref_bin     = tk.StringVar()

        # ── Outer layout: sidebar | step content | data panel ─────────────────
        outer = tk.Frame(parent, bg=BG)
        outer.pack(fill="both", expand=True)

        # Left sidebar (fixed width)
        sidebar = tk.Frame(outer, bg=SURFACE, width=148)
        sidebar.pack(side="left", fill="y")
        sidebar.pack_propagate(False)

        # Right-hand panel: spheroid state table only (full height, fixed width)
        side_panel = tk.Frame(outer, bg=BG2, width=380)
        side_panel.pack(side="right", fill="y")
        side_panel.pack_propagate(False)

        # Middle column: step content at top, live dashboard fills the rest
        middle = tk.Frame(outer, bg=BG)
        middle.pack(side="left", fill="both", expand=True)

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
        self._nd2_file_row(f1, "10X mosaic ND2:", self._pl_mosaic_path)
        row_w = tk.Frame(f1, bg=BG); row_w.pack(fill="x", pady=2)
        tk.Label(row_w, text="Well ID:", width=26, anchor="w",
                 bg=BG, fg=TEXT2, font=("Segoe UI", 9)).pack(side="left")
        tk.Entry(row_w, textvariable=self._pl_well_id, width=14,
                 bg=SURFACE, fg=TEXT, insertbackground=TEXT, relief="flat",
                 font=("Segoe UI", 9)).pack(side="left")
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
        f_s2 = tk.Frame(content_host, bg=BG)
        self._pl_step_frames["s2"] = f_s2

        tk.Label(f_s2,
                 text="Navigate stage to a spheroid in NIS-E, capture sub-10X nd2, type its rank, browse the nd2.\n"
                      "Leave rank blank to auto-match by NCC. After Verify, the box shows the actually matched rank.",
                 bg=BG, fg=SUBTEXT, font=("Segoe UI", 8), justify="left").pack(anchor="w", padx=12, pady=(8, 2))

        self._pl_anchor_frames = []
        self._pl_anchor_paths  = []
        self._pl_anchor_ranks  = []
        anchor_outer = tk.Frame(f_s2, bg=BG)
        anchor_outer.pack(fill="x", padx=12, pady=2)
        _anchor_labels = ["Anchor 1 (required):", "Anchor 2 (optional):"]
        for i in range(2):
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

        # Step 3 content (automated autofocus + global registration)
        f_s3 = tk.Frame(content_host, bg=BG)
        self._pl_step_frames["s3"] = f_s3

        tk.Label(f_s3,
                 text=("NIS-E macro nis_macro_z_autofocus.mac navigates to each spheroid's\n"
                       "corrected position, autofocuses, and captures a single-plane 20X ND2.\n"
                       "After captures complete, global registration refines all coordinates."),
                 bg=BG, fg=SUBTEXT, font=("Segoe UI", 8), justify="left"
                 ).pack(anchor="w", padx=12, pady=(8, 4))

        z_params = tk.Frame(f_s3, bg=BG); z_params.pack(fill="x", padx=12, pady=2)
        for lbl, var in [("Z half-range (um):", self._pl_z_half),
                          ("Z step (um):",       self._pl_z_step)]:
            tk.Label(z_params, text=lbl, bg=BG, fg=TEXT2,
                     font=("Segoe UI", 9)).pack(side="left", padx=(0, 4))
            tk.Entry(z_params, textvariable=var, width=7, bg=SURFACE, fg=TEXT,
                     insertbackground=TEXT, relief="flat",
                     font=("Segoe UI", 9)).pack(side="left", padx=(0, 16))

        btn3a = tk.Frame(f_s3, bg=BG); btn3a.pack(fill="x", padx=12, pady=(4, 0))
        self._btn(btn3a, "Trigger NIS-E Autofocus Captures",
                  self._pl_trigger_autofocus, BLUE, "#1e1e2e")
        self._pl_af_status_lbl = tk.Label(f_s3, text="Autofocus: not started.",
                                           bg=BG, fg=SUBTEXT, font=("Segoe UI", 9), anchor="w")
        self._pl_af_status_lbl.pack(fill="x", padx=12)

        btn3b = tk.Frame(f_s3, bg=BG); btn3b.pack(fill="x", padx=12, pady=(6, 0))
        self._btn(btn3b, "Refresh Status", self._pl_poll_af_thread, SURFACE2, TEXT)
        self._btn(btn3b, "Apply Global Registration",
                  self._pl_global_reg_thread, MAUVE, "#1e1e2e")
        self._pl_reg_lbl = tk.Label(f_s3, text="Registration: not run.",
                                     bg=BG, fg=SUBTEXT, font=("Segoe UI", 9), anchor="w")
        self._pl_reg_lbl.pack(fill="x", padx=12)

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
                 bg=BG, fg=SUBTEXT, font=("Segoe UI", 8)).pack(anchor="w", padx=2)
        lp_row = tk.Frame(s4, bg=BG); lp_row.pack(fill="x", pady=2)
        for lbl, var in [("Laser channel field:", self._pl_ch_field),
                          ("P0 (%):", self._pl_P0),
                          ("L (um):", self._pl_L_um)]:
            tk.Label(lp_row, text=lbl, bg=BG, fg=TEXT2,
                     font=("Segoe UI", 9)).pack(side="left", padx=(0, 4))
            tk.Entry(lp_row, textvariable=var, width=14, bg=SURFACE, fg=TEXT,
                     insertbackground=TEXT, relief="flat",
                     font=("Segoe UI", 9)).pack(side="left", padx=(0, 14))
        btn4 = tk.Frame(f_s4, bg=BG); btn4.pack(fill="x", padx=12, pady=(6, 2))
        self._btn(btn4, "Generate All Bins",   self._pl_generate_bins_thread, GREEN, "#1e1e2e")
        self._btn(btn4, "Start Capture Queue", self._pl_capture_queue_thread, MAUVE, "#1e1e2e")
        self._pl_capture_lbl = tk.Label(f_s4, text="Capture: idle",
                                         bg=BG, fg=SUBTEXT, font=("Segoe UI", 9), anchor="w")
        self._pl_capture_lbl.pack(fill="x", padx=12)

        # ── Right panel: spheroid state table (fills the column height) ───────
        tk.Label(side_panel, text="Spheroid State", bg=BG2, fg=MAUVE,
                 font=("Segoe UI", 9, "bold")).pack(anchor="w", padx=8, pady=(6, 2))

        tbl_frame = tk.Frame(side_panel, bg=BG2)
        tbl_frame.pack(fill="both", expand=True, padx=4, pady=(0, 6))

        cols = ("rank", "id", "status", "verified_xy", "z_centre", "bin")
        self._pl_tree = ttk.Treeview(tbl_frame, columns=cols, show="headings", height=8)
        hdrs = {"rank": ("Rank", 44), "id": ("ID", 130), "status": ("Status", 86),
                "verified_xy": ("Verified XY (um)", 140),
                "z_centre": ("Z (um)", 70), "bin": ("Bin", 150)}
        for c, (h, w_) in hdrs.items():
            self._pl_tree.heading(c, text=h)
            self._pl_tree.column(c, width=w_, minwidth=40, anchor="center")
        vsb_tree = ttk.Scrollbar(tbl_frame, orient="vertical", command=self._pl_tree.yview)
        hsb_tree = ttk.Scrollbar(tbl_frame, orient="horizontal", command=self._pl_tree.xview)
        self._pl_tree.configure(yscrollcommand=vsb_tree.set, xscrollcommand=hsb_tree.set)
        vsb_tree.pack(side="right", fill="y")
        hsb_tree.pack(side="bottom", fill="x")
        self._pl_tree.pack(side="left", fill="both", expand=True)

        # ── Middle dashboard: fills the large central area below the step form ─
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
        for r in self._pl_records:
            xy = (f"({r.verified_x_um:.1f}, {r.verified_y_um:.1f})"
                  if r.verified_x_um else "—")
            zc = f"{r.z_centre_um:.2f}" if r.z_centre_um else "—"
            bn = Path(r.bin_path).name if r.bin_path else "—"
            self._pl_tree.insert("", "end", values=(
                r.rank, r.spheroid_id, r.status, xy, zc, bn))

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

        dxs, dys = [], []
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
                dx, dy, ncc, matched_rank = _pl.estimate_offset_from_nd2(
                    Path(mosaic), Path(nd2_path), self._pl_records, Path(out), exp_rank)
                dxs.append(dx); dys.append(dy)
                self.after(0, lambda v=matched_rank, var=rk_var: var.set(str(v)))
                mism = (exp_rank is not None and matched_rank != exp_rank)
                self.after(0, lambda l=lbl, dx_=dx, dy_=dy, n=ncc, mr=matched_rank, mm=mism:
                           l.configure(
                               text=(f"matched rank {mr}  NCC={n:.4f}  "
                                     f"dx={dx_:+.1f} um  dy={dy_:+.1f} um"
                                     + ("  (differs from typed rank!)" if mm else "")),
                               fg=(YELLOW if mm else GREEN)))
                self._pl_log(f"  Anchor {i+1}: matched rank {matched_rank} "
                             f"NCC={ncc:.4f} dx={dx:+.1f} dy={dy:+.1f} um")
            except ValueError as exc:
                self.after(0, lambda l=lbl, e=exc: l.configure(
                    text=f"FAILED: {e}", fg=RED))
                self._pl_log(f"  Anchor {i+1} failed: {exc}")

        if not dxs:
            self.after(0, lambda: self._pl_offset_lbl.configure(
                text="Offset: FAILED — no successful anchor", fg=RED))
            self._pl_log("All anchors failed — cannot proceed.")
            return
        mean_dx = sum(dxs) / len(dxs)
        mean_dy = sum(dys) / len(dys)
        _pl.apply_offset(self._pl_records, mean_dx, mean_dy)
        self._pl_offset = (mean_dx, mean_dy)
        self.after(0, self._pl_update_table)
        self.after(0, self._pl_refresh_dashboard)
        msg = f"Offset applied: dx={mean_dx:+.1f} um  dy={mean_dy:+.1f} um  (from {len(dxs)} anchor(s))"
        self.after(0, lambda: self._pl_offset_lbl.configure(text=msg, fg=GREEN))
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
        n = _pl.trigger_autofocus_all(self._pl_records, Path(out) / "autofocus")
        self._pl_log(f"Step 3: {n} autofocus trigger files written to autofocus/")
        self.after(0, lambda: self._pl_af_status_lbl.configure(
            text=f"Triggers written for {n} spheroids. Start nis_macro_z_autofocus.mac in NIS-E.",
            fg=YELLOW))

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
            z_half, z_step = 18.0, 2.0
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

    def _pl_global_reg_thread(self):
        threading.Thread(target=self._pl_global_reg, daemon=True).start()

    def _pl_global_reg(self):
        import spheroid_pipeline as _pl
        out = self._pl_out_dir.get().strip()
        if not out or not self._pl_records:
            return
        done = _pl.poll_autofocus_done(self._pl_records, Path(out) / "autofocus")
        if not done:
            self.after(0, lambda: self._pl_reg_lbl.configure(
                text="No autofocus done files found — run captures first.", fg=RED)); return
        done_nd2s = {rank: d["nd2_path"] for rank, d in done.items()}
        mosaic = self._pl_mosaic_path.get().strip()
        if not mosaic:
            self.after(0, lambda: self._pl_reg_lbl.configure(
                text="Set mosaic ND2 path (Step 1).", fg=RED)); return
        self._pl_log(f"Step 3: running global registration on {len(done_nd2s)} ND2s ...")
        try:
            result = _pl.global_registration(
                self._pl_records, Path(mosaic), Path(out), done_nd2s)
            n  = result["n_matched"]
            mr = result["mean_residual_um"]
            t  = result["transform"]
            msg = (f"Global reg: {n} matches, mean residual={mr:.1f} um, "
                   f"dx={t.get('dx', 0):.1f} dy={t.get('dy', 0):.1f} um")
            if "scale" in t:
                msg += f" scale={t['scale']:.4f} angle={t.get('angle_deg', 0):.2f}deg"
            self.after(0, lambda: self._pl_reg_lbl.configure(text=msg, fg=GREEN))
            self._pl_log(f"Step 3 global reg done: {msg}")
        except ValueError as e:
            self.after(0, lambda: self._pl_reg_lbl.configure(
                text=f"Failed: {e}", fg=RED))
            self._pl_log(f"Step 3 global reg failed: {e}")
        self.after(0, self._pl_update_table)
        self.after(0, self._pl_refresh_dashboard)
        self._pl_update_step_dots()

    # ── Step 4 handlers ───────────────────────────────────────────────────────

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
        try:
            P0    = float(self._pl_P0.get())
            L_um  = float(self._pl_L_um.get())
        except ValueError:
            messagebox.showerror("Bad value", "P0 and L must be numbers."); return
        ch_field = self._pl_ch_field.get().strip() or "CH2LaserPower"
        ref_s   = self._pl_ref_bin.get().strip()
        ref_bin = Path(ref_s) if ref_s else None
        if ref_bin and not ref_bin.exists():
            self._pl_log(f"Reference .bin not found: {ref_bin}"); ref_bin = None
        if ref_bin is None:
            self._pl_log("WARNING: no reference .bin set — generated bins will be flagged "
                         "'Incompatible Z Correction and Camera' by NIS-E. Browse a rig export.")
        n_ok = 0
        for r in self._pl_records:
            if r.z_centre_um == 0.0:
                self._pl_log(f"Rank {r.rank}: skipping bin — Z not recorded"); continue
            try:
                _pl.generate_bin(r, P0, L_um, ch_field, bin_dir, reference_bin=ref_bin)
                n_ok += 1
                self._pl_log(f"Rank {r.rank}: bin written → {Path(r.bin_path).name}")
            except Exception as exc:
                self._pl_log(f"Rank {r.rank}: bin error — {exc}")
        self.after(0, self._pl_update_table)
        self.after(0, self._pl_refresh_dashboard)
        msg = f"Bins generated: {n_ok}/{len(self._pl_records)}"
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
        trigger_dir = self._pl_trigger_dir.get().strip()
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
