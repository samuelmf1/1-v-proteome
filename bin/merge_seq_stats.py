#!/usr/bin/env python3

import argparse
import pandas as pd
import statistics
import os
import glob

def parse_args():
    parser = argparse.ArgumentParser(description="Merge per-pair length TSVs and calculate global statistics.")
    parser.add_argument("--inputs", required=True, nargs="+", help="Input per-pair length TSV files (glob patterns supported)")
    parser.add_argument("--output_tsv", default="af3_lengths.tsv", help="Output merged TSV filename")
    parser.add_argument("--output_stats", default="af3_stats.txt", help="Output stats text filename")
    return parser.parse_args()

def main():
    args = parse_args()

    # 1. Collect and merge TSVs
    # If inputs are provided as glob patterns that the shell didn't expand
    all_files = []
    for pattern in args.inputs:
        files = glob.glob(pattern)
        if files:
            all_files.extend(files)
        else:
            # Assume it's a direct filename if glob returns nothing
            if os.path.exists(pattern):
                all_files.append(pattern)

    if not all_files:
        print("No input files found.")
        return

    print(f"Merging {len(all_files)} files...")
    
    # Read and concatenate all TSVs
    dfs = [pd.read_csv(f, sep='\t') for f in all_files]
    df = pd.concat(dfs, ignore_index=True)
    
    # 2. Save Merged TSV
    df.to_csv(args.output_tsv, sep='\t', index=False)
    print(f"Merged TSV saved to {args.output_tsv}")

    # 3. Calculate and Save Statistics
    n = len(df)
    if n > 0:
        with open(args.output_stats, 'w') as f:
            f.write(f"Sequence Statistics (N={n})\n")
            f.write("-" * 30 + "\n")
            f.write(f"Mean Length A:   {statistics.mean(df['len_a']):.2f}\n")
            f.write(f"Median Length A: {statistics.median(df['len_a']):.2f}\n")
            f.write("-" * 30 + "\n")
            f.write(f"Mean Length B:   {statistics.mean(df['len_b']):.2f}\n")
            f.write(f"Median Length B: {statistics.median(df['len_b']):.2f}\n")
            f.write("-" * 30 + "\n")
            f.write(f"Mean Sum (A+B):   {statistics.mean(df['len_sum']):.2f}\n")
            f.write(f"Median Sum (A+B): {statistics.median(df['len_sum']):.2f}\n")
        print(f"Summary stats written to {args.output_stats}")
    else:
        print("No records found to calculate statistics.")

if __name__ == "__main__":
    main()
