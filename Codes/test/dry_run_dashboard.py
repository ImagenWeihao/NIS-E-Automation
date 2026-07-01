"""
dry_run_dashboard.py — Dashboard test across all pipeline stages.

Runs the full dry-run pipeline and saves a dashboard PNG at every stage
to confirm the progressive transparency → solid fill-in looks correct.
"""
import sys
sys.stdout.reconfigure(encoding="utf-8")
import copy
from pathlib import Path
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "GUI"))  # Codes/GUI has the pipeline modules
import spheroid_pipeline as pl
from spheroid_pipeline import PipelineDashboard, Status

DEMO    = Path(r"C:\Users\weiha\SpheroidPA\demo_data")
SLIM    = DEMO / "2026-06-03 SLIM050 + SLIM047 QC check"
MOSAIC  = SLIM / "WellD05_Channel555_spIII1_Seq0000.nd2"
SUB10X  = SLIM / "spheroid 1" / "sph 1 1P" / "spheroid 1 10x 1P 555nm.nd2"
ND2_20X = SLIM / "spheroid 1" / "sph 1 1P" / "20x spheroid 1 1P 555nm.nd2"
OUT_DIR = DEMO / "pipeline_dryrun"

dash = PipelineDashboard(mosaic_nd2=MOSAIC, out_dir=OUT_DIR, interactive=False)

# ── Stage 0: freshly screened, all DETECTED ───────────────────────────────────
print("Running screener ...")
records = pl.screen_mosaic(MOSAIC, OUT_DIR, "WellD05")
dash.update(records, OUT_DIR / "dash_s0_detected.png")
print(f"  Saved dash_s0_detected.png  ({len(records)} DETECTED circles, 25% fill)")

# ── Stage 1: offset applied, all VERIFIED ─────────────────────────────────────
print("Estimating offset ...")
dx, dy, ncc, *_ = pl.estimate_offset_from_nd2(MOSAIC, SUB10X, records, OUT_DIR)
pl.apply_offset(records, dx, dy)
dash.update(records, OUT_DIR / "dash_s1_verified.png")
print(f"  Saved dash_s1_verified.png  (all VERIFIED, 50% fill)")

# ── Stage 2: z_centre recorded for rank 9 only ────────────────────────────────
pl.record_z_centre(records, 9, ND2_20X, z_half_um=18.0, z_step_um=2.0)
dash.update(records, OUT_DIR / "dash_s2_z_known.png")
print(f"  Saved dash_s2_z_known.png   (rank 9 Z_KNOWN, 65% fill)")

# ── Stage 3: simulate all z_centres known (fake Z for remaining ranks) ────────
for r in records:
    if r.z_centre_um == 0.0:
        r.z_centre_um = 7686.0
        r.z_half_um   = 18.0
        r.z_step_um   = 2.0
        r.status      = Status.Z_KNOWN
dash.update(records, OUT_DIR / "dash_s3_all_z.png")
print(f"  Saved dash_s3_all_z.png     (all Z_KNOWN, 65% fill)")

# ── Stage 4: bins generated for all ──────────────────────────────────────────
BIN_DIR = OUT_DIR / "bins"
for r in records:
    try:
        pl.generate_bin(r, P0_pct=15.0, L_um=165.0,
                        ch_field="CH1LaserPower", bin_dir=BIN_DIR)
    except Exception as exc:
        print(f"    rank {r.rank} bin failed: {exc}")
dash.update(records, OUT_DIR / "dash_s4_bins.png")
print(f"  Saved dash_s4_bins.png      (all BIN_READY, 75% fill)")

# ── Stage 5: simulate progressive imaging ─────────────────────────────────────
for i, r in enumerate(sorted(records, key=lambda r: r.rank)):
    r.status = Status.IMAGING
    dash.update(records, OUT_DIR / f"dash_s5_imaging_rank{r.rank:02d}.png")
    r.status = Status.IMAGED
dash.update(records, OUT_DIR / "dash_s6_all_done.png")
print(f"  Saved dash_s6_all_done.png  (all IMAGED, 100% fill, green)")

dash.close()
print("\nDashboard dry run complete. Check demo_data/pipeline_dryrun/dash_*.png")
