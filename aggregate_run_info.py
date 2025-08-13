#!/usr/bin/env python3
"""
Summarize Kallisto run_info.json files into a single CSV.

Looks for:
  kallisto_quant/
    <sampleA>/run_info.json
    <sampleB>/run_info.json
    ...

Outputs (default): kallisto_run_info.csv with columns:
  sample, n_targets, n_processed, n_pseudoaligned, n_unique, p_pseudoaligned, p_unique
"""

import argparse
import json
from pathlib import Path
import pandas as pd

FIELDS = [
    "n_targets",
    "n_processed",
    "n_pseudoaligned",
    "n_unique",
    "p_pseudoaligned",
    "p_unique",
]

def main():
    ap = argparse.ArgumentParser(description="Parse Kallisto run_info.json into a stats table.")
    ap.add_argument("-i", "--input-dir", default="kallisto_quant",
                    help="Directory with sample subfolders (default: %(default)s)")
    ap.add_argument("-o", "--output-csv", default="kallisto_run_info.csv",
                    help="Output CSV path (default: %(default)s)")
    args = ap.parse_args()

    in_dir = Path(args.input_dir)
    if not in_dir.exists():
        raise SystemExit(f"Input directory not found: {in_dir}")

    rows = []
    for sd in sorted(p for p in in_dir.iterdir() if p.is_dir()):
        run_info_path = sd / "run_info.json"
        if not run_info_path.exists():
            continue
        with run_info_path.open("r", encoding="utf-8") as fh:
            data = json.load(fh)
        row = {"sample": sd.name}
        for key in FIELDS:
            row[key] = data.get(key, None)
        rows.append(row)

    if not rows:
        raise SystemExit(f"No run_info.json files found under: {in_dir}")

    df = pd.DataFrame(rows).set_index("sample").sort_index()

    # Ensure numeric columns are numeric (coerce any odd values to NaN)
    for c in FIELDS:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    df.to_csv(args.output_csv)
    print(f"Wrote {args.output_csv} with shape {df.shape}")

if __name__ == "__main__":
    main()
