"""
dry_run_pipeline.py  —  Offline pipeline test (no NIS-E interaction)

Uses local demo data to exercise Steps 1-3 and bin generation.
Step 4 trigger/done handshake is skipped.
"""
import sys
sys.stdout.reconfigure(encoding="utf-8")

from pathlib import Path
import spheroid_pipeline as pl

DEMO = Path(r"C:\Users\weiha\SpheroidPA\demo_data")
SLIM = DEMO / "2026-06-03 SLIM050 + SLIM047 QC check"

MOSAIC   = SLIM / "WellD05_Channel555_spIII1_Seq0000.nd2"
SUB10X   = SLIM / "spheroid 1" / "sph 1 1P" / "spheroid 1 10x 1P 555nm.nd2"
ND2_20X  = SLIM / "spheroid 1" / "sph 1 1P" / "20x spheroid 1 1P 555nm.nd2"
OUT_DIR  = DEMO / "pipeline_dryrun"

print("=" * 60)
print("STEP 1: Screen 10X mosaic")
print("=" * 60)
records = pl.screen_mosaic(MOSAIC, OUT_DIR, "WellD05")
print(f"  Detected {len(records)} spheroid(s):")
for r in records:
    print(f"    Rank {r.rank:2d}  id={r.spheroid_id}  "
          f"diam={r.mosaic_diam_um:.1f} um  score={r.mosaic_score:.4f}  "
          f"status={r.status}")

print()
print("=" * 60)
print("STEP 2: Anchor offset estimation")
print("=" * 60)
anchors = pl.pick_anchors(records, n=2, seed=42)
print(f"  Selected anchors: ranks {[a.rank for a in anchors]}")

print(f"  Running NCC verification with sub-10X: {SUB10X.name}")
try:
    dx, dy, ncc = pl.estimate_offset_from_nd2(MOSAIC, SUB10X, records, OUT_DIR)
    print(f"  NCC={ncc:.4f}  dx={dx:+.1f} um  dy={dy:+.1f} um")
    pl.apply_offset(records, dx, dy)
    print(f"  Offset applied to all {len(records)} records.")
except ValueError as exc:
    print(f"  OFFSET FAILED: {exc}")
    sys.exit(1)

print()
print("  Records after offset:")
for r in records:
    print(f"    Rank {r.rank:2d}  verified=({r.verified_x_um:.1f}, {r.verified_y_um:.1f}) um  "
          f"status={r.status}")

print()
print("=" * 60)
print("STEP 3: Record Z-centre from 20X autofocus nd2")
print("=" * 60)

# From previous session: the 20X nd2 corresponds to Rank 9 (NCC=0.991)
Z_RANK   = 9
Z_HALF   = 18.0   # 18 um -> 19 items at 2 um step -> satisfies ROOT_TAIL_U64S = 19
Z_STEP   = 2.0

print(f"  Reading Z from {ND2_20X.name} for rank {Z_RANK} ...")
try:
    pl.record_z_centre(records, Z_RANK, ND2_20X, z_half_um=Z_HALF, z_step_um=Z_STEP)
    rec9 = next(r for r in records if r.rank == Z_RANK)
    print(f"  Rank {Z_RANK}: z_centre={rec9.z_centre_um:.2f} um  "
          f"z_range=[{rec9.z_centre_um - Z_HALF:.2f}, {rec9.z_centre_um + Z_HALF:.2f}]  "
          f"status={rec9.status}")
except ValueError as exc:
    print(f"  Z read failed: {exc}")

print()
print("=" * 60)
print("STEP 4a: Generate .bin for Rank 9")
print("=" * 60)

BIN_DIR  = OUT_DIR / "bins"
P0_PCT   = 15.0
L_UM     = 165.0
CH_FIELD = "CH1LaserPower"

rec9 = next(r for r in records if r.rank == Z_RANK)
try:
    bin_path = pl.generate_bin(rec9, P0_PCT, L_UM, CH_FIELD, BIN_DIR)
    print(f"  Bin written: {bin_path}")
    print(f"  Size: {bin_path.stat().st_size} bytes")
    print(f"  Status: {rec9.status}")
except Exception as exc:
    print(f"  Bin generation failed: {exc}")

print()
print("=" * 60)
print("STEP 4b: Write trigger.ini (offline — no NIS-E)")
print("=" * 60)
TRIGGER_DIR = OUT_DIR / "trigger"
ND2_OUT_DIR = OUT_DIR / "nd2_out"

if rec9.status == pl.Status.BIN_READY:
    trig = pl.write_trigger(rec9, TRIGGER_DIR, ND2_OUT_DIR)
    print(f"  Trigger written: {trig}")
    print()
    print(trig.read_text(encoding="utf-8"))
    print(f"  Status: {rec9.status}")
else:
    print(f"  Skipped — rank 9 status is {rec9.status}, expected BIN_READY")

print()
print("=" * 60)
print("Pipeline state export")
print("=" * 60)
pl.export_state_csv(records, OUT_DIR / "pipeline_state.csv")
print(f"  Saved: {OUT_DIR / 'pipeline_state.csv'}")
pl.save_state(pl.PipelineState(
    well_id     = "WellD05",
    mosaic_nd2  = str(MOSAIC),
    out_dir     = str(OUT_DIR),
    trigger_dir = str(TRIGGER_DIR),
    nd2_out_dir = str(ND2_OUT_DIR),
    records     = records,
    offset_x    = dx,
    offset_y    = dy,
))
print(f"  JSON state saved: {OUT_DIR / 'pipeline_state.json'}")

print()
print("Dry run complete.")
