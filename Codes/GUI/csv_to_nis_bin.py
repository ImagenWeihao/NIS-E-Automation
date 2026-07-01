"""
csv_to_nis_bin.py
Inverse of nis_bin_to_csv.py — write a Nikon NIS-E "Z Intensity Correction"
.bin file from the CSV we generated.

Usage:
    python csv_to_nis_bin.py <input.csv> [output.bin]
"""

import csv
import struct
import sys
from pathlib import Path

# ── Display-name → internal-name reverse map (must match nis_bin_to_csv.py) ──

DISPLAY_TO_INTERNAL = {
    "Z [um]":         "ZStack",
    "Exposure (ms)":  "ExposureTime",
    "HV1":            "CH1PMTHighVoltage",
    "HV2":            "CH2PMTHighVoltage",
    "HV3":            "CH3PMTHighVoltage",
    "HV4":            "CH4PMTHighVoltage",
    "TD HV":          "TDPMTHighVoltage",
    "LP opt":         "OptionCHLaserPower",
    "LP1 (%)":        "CH1LaserPower",
    "LP2 (%)":        "CH2LaserPower",
    "LP3 (%)":        "CH3LaserPower",
    "LP opt2":        "OptionCH2LaserPower",
    "LP4 (%)":        "CH4LaserPower",
    "LP5 (%)":        "CH5LaserPower",
    "SD HV":          "SDHighVoltage",
    "Gain opt":       "OptionCHGain",
    "Gain VF1":       "VFCH1Gain",
    "Gain VF2":       "VFCH2Gain",
    "Gain opt2":      "OptionCH2Gain",
    "Gain VF3":       "VFCH3Gain",
    "Gain VF4":       "VFCH4Gain",
    "Gain VF5":       "VFCH5Gain",
    "VAAS HV1":       "VAASCH1NearPMTHighVoltage",
    "VAAS HV2":       "VAASCH2NearPMTHighVoltage",
    "VAAS HV3":       "VAASCH3NearPMTHighVoltage",
    "VAAS HV4":       "VAASCH4NearPMTHighVoltage",
    "ShowFlags":      "iShowFlags",
    "ItemType":       "eItemType",
    "IsReference":    "bReference",
}

# ── Item field schema (exact order, exact internal names, exact types) ───────

ITEM_FIELDS_SCHEMA = [
    # (name, type, default)
    ("ZStack",                       "f64", 0.0),
    ("ExposureTime",                 "f64", 0.0),
    ("CH1PMTHighVoltage",            "f64", 0.0),
    ("CH2PMTHighVoltage",            "f64", 0.0),
    ("CH3PMTHighVoltage",            "f64", 0.0),
    ("CH4PMTHighVoltage",            "f64", 0.0),
    ("TDPMTHighVoltage",             "f64", 0.0),
    ("OptionCHLaserPower",           "f64", 0.0),
    ("CH1LaserPower",                "f64", 0.0),
    ("CH2LaserPower",                "f64", 0.0),
    ("CH3LaserPower",                "f64", 0.0),
    ("OptionCH2LaserPower",          "f64", 0.0),
    ("CH4LaserPower",                "f64", 0.0),
    ("CH5LaserPower",                "f64", 0.0),
    ("SDHighVoltage",                "f64", 0.0),
    ("OptionCHGain",                 "f64", 0.0),
    ("VFCH1Gain",                    "f64", 0.0),
    ("VFCH2Gain",                    "f64", 0.0),
    ("OptionCH2Gain",                "f64", 0.0),
    ("VFCH3Gain",                    "f64", 0.0),
    ("VFCH4Gain",                    "f64", 0.0),
    ("VFCH5Gain",                    "f64", 0.0),
    ("VAASCH1NearPMTHighVoltage",    "f64", 0.0),
    ("VAASCH2NearPMTHighVoltage",    "f64", 0.0),
    ("VAASCH3NearPMTHighVoltage",    "f64", 0.0),
    ("VAASCH4NearPMTHighVoltage",    "f64", 0.0),
    ("iShowFlags",                   "i32", 0),
    ("eItemType",                    "i32", 0),
    ("bReference",                   "bool", False),
]

# Observed Item0 footer pattern (29 u64s). Items 0-9 use these values;
# items 10-18 use each value + 2 (matching the item size difference of 2).
ITEM0_FOOTER = (
    362, 86, 400, 132, 438, 178, 526, 224, 564, 50,
    736, 476, 640, 314, 602, 270, 864, 926, 988, 1050,
    676, 706, 774, 804, 834, 26, 1166, 1140, 1112,
)

# Container "size" field values observed
ITEM_SIZE_FIRST10 = 1191
ITEM_SIZE_REST    = 1193

# Root-level tail: root_count * 8 bytes (root_count u64s = root container footer)
# root_count = 6 metadata + n_items; for 19 items root_count=25 → 25 u64s.
ROOT_TAIL_U64S = (
    196, 226, 164, 293, 1716,
    14523, 15948, 17373, 18798, 20223, 21648, 23073, 24498, 25923,
    3139, 4562, 5985, 7408, 8831, 10254, 11677, 13100,
    267, 112, 56,
)


# ── Low-level writers ────────────────────────────────────────────────────────

def _write_entry_header(out: bytearray, type_byte: int, name: str) -> None:
    """Write: type, name_len, NUL-terminated UTF-16-LE name."""
    encoded = (name + "\x00").encode("utf-16-le")
    name_len = len(encoded) // 2          # UTF-16 char count incl. NUL
    out += bytes([type_byte, name_len])
    out += encoded


def _write_value(out: bytearray, type_str: str, value) -> None:
    if   type_str == "f64":   out += struct.pack("<d", float(value))
    elif type_str == "i32":   out += struct.pack("<i", int(value))
    elif type_str == "bool":  out += bytes([1 if value else 0])
    elif type_str == "str":
        out += (str(value) + "\x00").encode("utf-16-le")
    else:
        raise ValueError(f"Unknown type: {type_str}")


def _write_item_container(out: bytearray, item_idx: int, fields: dict,
                          size_override: int = None) -> None:
    """Write one ItemN container with all 29 named children + 29 footer u64s."""
    name = f"Item{item_idx}"

    # Header: type, len, name, count(u32), size(u64)
    _write_entry_header(out, 0x0b, name)
    out += struct.pack("<I", len(ITEM_FIELDS_SCHEMA))             # count = 29
    size = size_override if size_override is not None else ITEM_SIZE_FIRST10
    out += struct.pack("<Q", size)

    # Children (all 29 fields in canonical order)
    for fname, ftype, default in ITEM_FIELDS_SCHEMA:
        val = fields.get(fname, default)
        if   ftype == "f64":  _write_entry_header(out, 0x06, fname)
        elif ftype == "i32":  _write_entry_header(out, 0x02, fname)
        elif ftype == "bool": _write_entry_header(out, 0x01, fname)
        _write_value(out, ftype, val)

    # Footer (29 u64 values). Items 0-9 use ITEM0_FOOTER; items 10+ use +2.
    bump = 0 if item_idx < 10 else 2
    for v in ITEM0_FOOTER:
        out += struct.pack("<Q", v + bump)


# ── CSV reader ───────────────────────────────────────────────────────────────

def read_csv(csv_path: Path) -> tuple:
    """Return (metadata_dict, list_of_item_dicts)."""
    metadata = {}
    data_lines = []
    with open(csv_path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                stripped = line.lstrip("#").strip()
                if ":" in stripped:
                    k, v = stripped.split(":", 1)
                    metadata[k.strip()] = v.strip()
            else:
                data_lines.append(line)

    # Parse the data section as CSV
    import io
    reader = csv.DictReader(io.StringIO("".join(data_lines)))
    items = []
    for row in reader:
        item = {}
        for disp_col, val_str in row.items():
            if disp_col == "item" or not val_str.strip():
                continue
            internal = DISPLAY_TO_INTERNAL.get(disp_col, disp_col)
            # type from schema
            type_str = next(t for n, t, _ in ITEM_FIELDS_SCHEMA if n == internal)
            if   type_str == "f64":  item[internal] = float(val_str)
            elif type_str == "i32":  item[internal] = int(val_str)
            elif type_str == "bool": item[internal] = val_str.lower() in ("true","1","yes")
        items.append(item)
    return metadata, items


# ── Main builder ─────────────────────────────────────────────────────────────

def build_bin(metadata: dict, items: list) -> bytes:
    out = bytearray()
    out += b"LV00"

    # Root container header
    _write_entry_header(out, 0x0b, "CLxZIntensityControl")

    # Root has 6 metadata children + N item containers
    n_items = len(items)
    root_count = 6 + n_items
    out += struct.pack("<I", root_count)

    # Root size (u64) — we'll patch this after we know the actual content size
    size_field_offset = len(out)
    out += struct.pack("<Q", 0)
    root_content_start = len(out)

    # ── Root metadata children, in observed order ────────────────────────────
    _write_entry_header(out, 0x08, "hwUnit_Name")
    _write_value(out, "str", metadata.get("hwUnit_Name", ""))

    _write_entry_header(out, 0x08, "hwUnit_ConnectionString")
    _write_value(out, "str", metadata.get("hwUnit_ConnectionString", ""))

    _write_entry_header(out, 0x02, "DetectorType")
    _write_value(out, "i32", int(metadata.get("DetectorType", 0)))

    _write_entry_header(out, 0x02, "ChannelBits")
    _write_value(out, "i32", int(metadata.get("ChannelBits", 0)))

    _write_entry_header(out, 0x01, "CorrectionRelative")
    _write_value(out, "bool", str(metadata.get("CorrectionRelative", "False")).lower() == "true")

    _write_entry_header(out, 0x02, "ItemCount")
    _write_value(out, "i32", n_items)

    # ── Item containers ─────────────────────────────────────────────────────
    for i, item in enumerate(items):
        size = ITEM_SIZE_FIRST10 if i < 10 else ITEM_SIZE_REST
        _write_item_container(out, i, item, size_override=size)

    # ── Root footer: 25 u64s + 3 trailing u64s (28 total) ────────────────────
    for v in ROOT_TAIL_U64S:
        out += struct.pack("<Q", v)

    # root size = total_bytes - LV00(4) - root_tail(root_count*8)
    # The tail is the root container's own footer (child_count * 8 bytes).
    root_content_size = len(out) - root_count * 8 - 4
    struct.pack_into("<Q", out, size_field_offset, root_content_size)

    return bytes(out)


# ── CLI ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python csv_to_nis_bin.py <input.csv> [output.bin]")
        sys.exit(1)

    csv_in = Path(sys.argv[1])
    bin_out = Path(sys.argv[2]) if len(sys.argv) >= 3 \
              else csv_in.with_name(csv_in.stem + "_converted.bin")

    metadata, items = read_csv(csv_in)
    print(f"CSV:       {csv_in}")
    print(f"Metadata:  {metadata}")
    print(f"Items:     {len(items)}")
    print(f"First item: {items[0] if items else 'none'}")

    blob = build_bin(metadata, items)
    bin_out.write_bytes(blob)
    print(f"Wrote:     {bin_out}  ({len(blob)} bytes)")
