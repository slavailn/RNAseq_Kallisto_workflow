#!/usr/bin/env python3
"""
Aggregate Kallisto transcript-level outputs (est_counts and tpm)
into two CSV matrices with transcripts as rows and samples as columns.

Assumptions:
- Each sample subdir under input_dir contains an 'abundance.tsv'
- No duplicate target_id rows per file
- 'est_counts' and 'tpm' are numeric

USAGE:
 python3 aggregate_kallisto.py \
  -i kallisto_quant \
  --counts-out transcript_counts.csv \
  --tpm-out transcript_TPM.csv


"""

import argparse
from pathlib import Path
import pandas as pd

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input-dir", default="kallisto_quant",
                    help="Directory with sample subfolders (default: %(default)s)")
    ap.add_argument("--counts-out", default="transcript_counts.csv",
                    help="Output CSV for counts (default: %(default)s)")
    ap.add_argument("--tpm-out", default="transcript_TPM.csv",
                    help="Output CSV for TPM (default: %(default)s)")
    ap.add_argument("--fillna", type=float, default=0.0,
                    help="Fill value for missing transcripts (default: %(default)s)")
    args = ap.parse_args()

    in_dir = Path(args.input_dir)
    sample_dirs = sorted(p for p in in_dir.iterdir() if p.is_dir() and (p / "abundance.tsv").exists())
    if not sample_dirs:
        raise SystemExit(f"No sample subdirectories with abundance.tsv found in: {in_dir}")

    counts_cols = []
    tpm_cols = []

    for sd in sample_dirs:
        sample = sd.name
        df = pd.read_csv(sd / "abundance.tsv", sep="\t",
                         usecols=["target_id", "est_counts", "tpm"])
        counts_cols.append(df.set_index("target_id")["est_counts"].rename(sample))
        tpm_cols.append(df.set_index("target_id")["tpm"].rename(sample))

    counts = pd.concat(counts_cols, axis=1, join="outer").fillna(args.fillna).sort_index()
    tpm = pd.concat(tpm_cols, axis=1, join="outer").fillna(args.fillna).sort_index()

    counts.to_csv(args.counts_out)
    tpm.to_csv(args.tpm_out)

    print(f"Wrote {args.counts_out} with shape {counts.shape}")
    print(f"Wrote {args.tpm_out} with shape {tpm.shape}")

if __name__ == "__main__":
    main()
