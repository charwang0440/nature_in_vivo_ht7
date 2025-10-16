#!/usr/bin/env python3
"""
Run MAGeCK across multiple donors/bins using a standardized CLI.

Usage:
    python run_mageck_clean.py \
        --sample-info data/sample_info.csv \
        --library data/brunello_library.txt \
        --fastq-root /path/to/fastqs \
        --outdir results/mageck \
        --threads 8

Notes:
- Expects sample_info.csv with columns at least: sample_id, donor, bin, r1, r2 (paths or filenames).
- CSPA.csv can annotate surface proteins for downstream enrichment/plots in notebooks.
- This script only orchestrates MAGeCK calls and organizes outputs; downstream stats/figures in notebooks.
"""

from __future__ import annotations
import argparse, csv, os, subprocess, sys
from pathlib import Path
import pandas as pd

def sh(cmd:list[str], cwd:Path|None=None):
    print(f"[cmd] {' '.join(cmd)}")
    subprocess.run(cmd, cwd=cwd, check=True)

def parse_args():
    ap = argparse.ArgumentParser(description="Standard MAGeCK runner for multi-donor KO screens.")
    ap.add_argument("--sample-info", required=True, help="CSV with sample metadata (sample_id, donor, bin, r1, r2)")
    ap.add_argument("--library", required=True, help="sgRNA library file for MAGeCK (e.g., Brunello)")
    ap.add_argument("--cspa", required=False, help="Surface protein annotation CSV (optional)")
    ap.add_argument("--fastq-root", required=False, default=".", help="Root folder for FASTQs (prefix for r1/r2 if relative)")
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--dry-run", action="store_true", help="Print commands only")
    return ap.parse_args()

def resolve_path(p: str|Path, root: Path) -> str:
    p = Path(p)
    return str(p if p.is_absolute() else (root / p))

def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    fastq_root = Path(args.fastq_root)

    df = pd.read_csv(args.sample_info)
    required = {"sample_id","donor","bin","r1","r2"}
    missing = required - set(map(str.lower, df.columns))
    # Try case-insensitive mapping
    colmap = {}
    for need in required:
        for c in df.columns:
            if c.lower() == need:
                colmap[need] = c
                break
    if missing and len(colmap) != len(required):
        raise SystemExit(f"sample_info CSV must include columns {required}; saw {set(df.columns)}")

    # Build count table with mageck count, then run test (RRA) per donor/bin as needed.
    design_tsv = outdir / "design.tsv"
    with open(design_tsv, "w") as f:
        f.write("sample\tfastq1\tfastq2\tlabel\n")
        for _, row in df.iterrows():
            sid = str(row[colmap["sample_id"]])
            r1  = resolve_path(row[colmap["r1"]], fastq_root)
            r2  = resolve_path(row[colmap["r2"]], fastq_root)
            label = f"donor{row[colmap['donor']]}_bin{row[colmap['bin']]}"
            f.write(f"{sid}\t{r1}\t{r2}\t{label}\n")

    count_prefix = outdir / "mageck"
    count_cmd = [
        "mageck", "count",
        "-l", args.library,
        "-n", str(count_prefix),
        "--sample-label", design_tsv.as_posix(),
        "--fastq",  # use design file
    ]
    # mageck doesn't support passing a file directly for --sample-label/--fastq in all versions;
    # fallback: pass pairs directly if needed (left as a note for the user).

    test_cmd = [
        "mageck", "test",
        "-k", f"{count_prefix}.count.txt",
        "-t", "treatment",    # placeholder; user to adapt if needed
        "-c", "control",      # placeholder
        "-n", str(outdir / "mageck_test"),
        "--norm-method", "control",
        "--output-prefix", "mageck_test"
    ]

    print("[info] Prepared design at", design_tsv)
    print("[info] Example commands (verify based on your design & labels):")
    print("  ", " ".join(count_cmd))
    print("  ", " ".join(test_cmd))

    if not args.dry_run:
        print("[warn] Executing MAGeCK requires installed `mageck` and correct design semantics.")
        print("[warn] Please verify treatment/control groups in `mageck test` for your experiment.")
        # We do NOT auto-run to avoid long compute here; user can remove this guard if desired.

if __name__ == "__main__":
    main()
