"""
nis_bin_to_csv.py
Parse Nikon NIS-Elements "Z Intensity Correction" .bin files and
convert them to a flat, human-readable CSV.

Binary format (reverse-engineered from observation):
  Magic:  "LV00"  (4 bytes)
  Then a tree of typed entries; each entry is:
     type:      u8
     name_len:  u8   (UTF-16 character count INCLUDING the trailing NUL)
     name:      UTF-16-LE  (name_len * 2 bytes)
     value:     depends on type
        0x01 bool     1 byte
        0x02 int32    4 bytes LE
        0x06 double   8 bytes LE IEEE-754
        0x08 string   UTF-16-LE, NUL-terminated
        0x0b container {child_count: u32 LE, byte_size: u64 LE, children...}
"""

import csv
import struct
import sys
from pathlib import Path


# ── Parser ────────────────────────────────────────────────────────────────────

class NisBinParser:

    TYPE_BOOL      = 0x01
    TYPE_INT32     = 0x02
    TYPE_DOUBLE    = 0x06
    TYPE_STRING    = 0x08
    TYPE_CONTAINER = 0x0b

    def __init__(self, data: bytes):
        self.data = data
        self.pos  = 0

    # ── primitive readers ────────────────────────────────────────────────────

    def _u8(self) -> int:
        v = self.data[self.pos]; self.pos += 1
        return v

    def _u32(self) -> int:
        v = struct.unpack_from("<I", self.data, self.pos)[0]
        self.pos += 4
        return v

    def _u64(self) -> int:
        v = struct.unpack_from("<Q", self.data, self.pos)[0]
        self.pos += 8
        return v

    def _f64(self) -> float:
        v = struct.unpack_from("<d", self.data, self.pos)[0]
        self.pos += 8
        return v

    def _utf16_fixed(self, char_count: int) -> str:
        raw = self.data[self.pos : self.pos + char_count * 2]
        self.pos += char_count * 2
        return raw.decode("utf-16-le").rstrip("\x00")

    def _utf16_nullterm(self) -> str:
        start = self.pos
        while self.pos + 1 < len(self.data) and not (
                self.data[self.pos] == 0 and self.data[self.pos + 1] == 0):
            self.pos += 2
        raw = self.data[start : self.pos]
        self.pos += 2  # skip NUL
        return raw.decode("utf-16-le")

    # ── entry parser ─────────────────────────────────────────────────────────

    def _parse_entry(self) -> tuple:
        type_byte = self._u8()
        name_len  = self._u8()
        name      = self._utf16_fixed(name_len)

        if type_byte == self.TYPE_BOOL:
            value = bool(self._u8())
        elif type_byte == self.TYPE_INT32:
            value = self._u32()
        elif type_byte == self.TYPE_DOUBLE:
            value = self._f64()
        elif type_byte == self.TYPE_STRING:
            value = self._utf16_nullterm()
        elif type_byte == self.TYPE_CONTAINER:
            child_count = self._u32()
            _byte_size  = self._u64()    # not used — size field is unreliable in this format
            children    = {}
            for _ in range(child_count):
                child_name, child_value = self._parse_entry()
                children[child_name] = child_value
            # NIS-E containers append a child_count * 8 byte u64 footer table
            # after the named children. Skip past it before returning.
            self.pos = min(self.pos + child_count * 8, len(self.data))
            value = children
        else:
            raise ValueError(
                f"Unknown type byte 0x{type_byte:02x} at offset {self.pos - 1}")

        return name, value

    # ── public entry point ───────────────────────────────────────────────────

    def parse(self) -> dict:
        magic = self.data[:4].decode("ascii", errors="replace")
        if magic != "LV00":
            raise ValueError(f"Bad magic: {magic!r}  (expected 'LV00')")
        self.pos = 4
        root_name, root_value = self._parse_entry()
        return {root_name: root_value}


# ── Conversion to CSV ─────────────────────────────────────────────────────────

def _is_meaningful(value) -> bool:
    """Drop fields that are zero/empty for ALL items (they carry no info)."""
    if isinstance(value, bool):    return value is True
    if isinstance(value, (int,)):  return value != 0
    if isinstance(value, float):   return value != 0.0
    if isinstance(value, str):     return value != ""
    return True


# Human-readable column names matching the NIS-E "Z Intensity Correction" panel.
# The panel shows fields as e.g. "(HV2: 5) (LP2: 53.8)", which map to internal
# names CH2PMTHighVoltage and CH2LaserPower.
COLUMN_RENAMES = {
    "ZStack":                       "Z [um]",
    "ExposureTime":                 "Exposure (ms)",
    "CH1PMTHighVoltage":            "HV1",
    "CH2PMTHighVoltage":            "HV2",
    "CH3PMTHighVoltage":            "HV3",
    "CH4PMTHighVoltage":            "HV4",
    "TDPMTHighVoltage":             "TD HV",
    "OptionCHLaserPower":           "LP opt",
    "CH1LaserPower":                "LP1 (%)",
    "CH2LaserPower":                "LP2 (%)",
    "CH3LaserPower":                "LP3 (%)",
    "OptionCH2LaserPower":          "LP opt2",
    "CH4LaserPower":                "LP4 (%)",
    "CH5LaserPower":                "LP5 (%)",
    "SDHighVoltage":                "SD HV",
    "OptionCHGain":                 "Gain opt",
    "VFCH1Gain":                    "Gain VF1",
    "VFCH2Gain":                    "Gain VF2",
    "OptionCH2Gain":                "Gain opt2",
    "VFCH3Gain":                    "Gain VF3",
    "VFCH4Gain":                    "Gain VF4",
    "VFCH5Gain":                    "Gain VF5",
    "VAASCH1NearPMTHighVoltage":    "VAAS HV1",
    "VAASCH2NearPMTHighVoltage":    "VAAS HV2",
    "VAASCH3NearPMTHighVoltage":    "VAAS HV3",
    "VAASCH4NearPMTHighVoltage":    "VAAS HV4",
    "iShowFlags":                   "ShowFlags",
    "eItemType":                    "ItemType",
    "bReference":                   "IsReference",
}

# Decimal places to display for each field, matching the NIS-E panel format.
# Z shows 2 dp (e.g. 7754.23); laser power shows 1 dp (e.g. 53.8); HV shows 1 dp.
COLUMN_DECIMALS = {
    "ZStack":                    2,
    "ExposureTime":              3,
    "CH1PMTHighVoltage":         1,
    "CH2PMTHighVoltage":         1,
    "CH3PMTHighVoltage":         1,
    "CH4PMTHighVoltage":         1,
    "TDPMTHighVoltage":          1,
    "SDHighVoltage":             1,
    "OptionCHLaserPower":        1,
    "CH1LaserPower":             1,
    "CH2LaserPower":             1,
    "CH3LaserPower":             1,
    "OptionCH2LaserPower":       1,
    "CH4LaserPower":             1,
    "CH5LaserPower":             1,
    "OptionCHGain":              1,
    "VFCH1Gain":                 1,
    "VFCH2Gain":                 1,
    "OptionCH2Gain":             1,
    "VFCH3Gain":                 1,
    "VFCH4Gain":                 1,
    "VFCH5Gain":                 1,
    "VAASCH1NearPMTHighVoltage": 1,
    "VAASCH2NearPMTHighVoltage": 1,
    "VAASCH3NearPMTHighVoltage": 1,
    "VAASCH4NearPMTHighVoltage": 1,
}


def bin_to_csv(bin_path: Path, csv_path: Path,
               include_zero_columns: bool = False) -> dict:
    """
    Parse a NIS-E z-intensity-correction .bin file and write a CSV.
    Returns a dict of top-level metadata (excluding the Item* rows).
    """
    data   = bin_path.read_bytes()
    parsed = NisBinParser(data).parse()

    # The single root container holds everything.
    root_name, root = next(iter(parsed.items()))

    # Separate scalar metadata from indexed Item* containers.
    items_by_idx = {}
    metadata     = {}
    for k, v in root.items():
        if k.startswith("Item") and k[4:].isdigit() and isinstance(v, dict):
            items_by_idx[int(k[4:])] = v
        else:
            metadata[k] = v

    items = [items_by_idx[i] for i in sorted(items_by_idx)]
    if not items:
        raise ValueError("No Item* entries found in file.")

    # Determine columns: union of all keys across all items.
    all_keys = list(items[0].keys())                          # preserve order
    seen     = set(all_keys)
    for it in items[1:]:
        for k in it:
            if k not in seen:
                all_keys.append(k); seen.add(k)

    # Optionally drop columns whose value is zero/empty in every item.
    if not include_zero_columns:
        kept = []
        for k in all_keys:
            if any(_is_meaningful(it.get(k)) for it in items):
                kept.append(k)
        all_keys = kept

    # Map internal field names to human-readable column headers
    display_names = [COLUMN_RENAMES.get(k, k) for k in all_keys]
    fieldnames    = ["item", *display_names]

    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        # Emit metadata as a small header so the table is self-describing
        f.write(f"# Source: {bin_path.name}\n")
        f.write(f"# Root:   {root_name}\n")
        for mk, mv in metadata.items():
            f.write(f"# {mk}: {mv}\n")
        f.write(f"# Items:  {len(items)}\n")
        f.write("#\n")

        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for i, it in enumerate(items):
            row = {"item": i}
            for orig_k, disp_k in zip(all_keys, display_names):
                v = it.get(orig_k, "")
                if isinstance(v, float) and orig_k in COLUMN_DECIMALS:
                    dp = COLUMN_DECIMALS[orig_k]
                    v = f"{v:.{dp}f}"
                row[disp_k] = v
            w.writerow(row)

    return {"root": root_name, "n_items": len(items), "metadata": metadata,
            "columns": display_names, "csv_path": str(csv_path)}


# ── CLI ───────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage:  python nis_bin_to_csv.py  <input.bin>  [output.csv]")
        sys.exit(1)

    bin_in = Path(sys.argv[1])
    if not bin_in.exists():
        print(f"ERROR: file not found: {bin_in}")
        sys.exit(2)

    csv_out = Path(sys.argv[2]) if len(sys.argv) >= 3 \
              else bin_in.with_suffix(".csv")

    info = bin_to_csv(bin_in, csv_out)

    print(f"Parsed:  {bin_in}")
    print(f"Root:    {info['root']}")
    print(f"Items:   {info['n_items']}")
    print(f"Columns: {len(info['columns'])}  ->  {info['columns']}")
    print()
    print("Top-level metadata:")
    for k, v in info["metadata"].items():
        print(f"  {k:24s}  {v!r}")
    print()
    print(f"Wrote: {csv_out}")
