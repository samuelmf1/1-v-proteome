#!/usr/bin/env python3
# filter_interactions.py
import pandas as pd
import argparse
from pathlib import Path
import os
import sys

def parse_args():
    parser = argparse.ArgumentParser(description="Filter interaction CSVs by scores and binding motifs.")
    parser.add_argument("--input", type=str, required=True, help="Input CSV file or directory containing CSVs")
    parser.add_argument("--outdir", type=str, default="filtered_interactions", help="Output directory")
    
    # Score thresholds
    parser.add_argument("--plddt_min", type=float, default=65.0, help="Minimum pLDDT for both residues")
    parser.add_argument("--pae_max", type=float, default=10.0, help="Maximum PAE for interaction")
    parser.add_argument("--iptm_min", type=float, default=0.65, help="Minimum global ipTM score")
    parser.add_argument("--ptm_min", type=float, default=0.65, help="Minimum global pTM score")
    
    # Motif parameters
    parser.add_argument("--motif_gap", type=int, default=10, help="Max sequence distance to cluster residues")
    parser.add_argument("--motif_min_res", type=int, default=2, help="Min residues in a cluster to valid motif")
    
    return parser.parse_args()

def find_motifs(residues, gap, min_res):
    """
    Cluster sorted residues. Returns True if any cluster has >= min_res.
    residues: sorted list of integers (residue numbers)
    gap: max distance to link residues
    min_res: minimum size of a valid cluster
    """
    if not residues:
        return False, []
        
    clusters = []
    current_cluster = [residues[0]]
    
    for r in residues[1:]:
        if r - current_cluster[-1] <= gap:
            current_cluster.append(r)
        else:
            clusters.append(current_cluster)
            current_cluster = [r]
    clusters.append(current_cluster)
    
    valid_clusters = [c for c in clusters if len(c) >= min_res]
    return len(valid_clusters) > 0, valid_clusters

def process_file(csv_path, args):
    try:
        df = pd.read_csv(csv_path)
    except Exception as e:
        return None, f"Error reading {csv_path}: {e}"

    # 1. Global score filter (iptm, ptm)
    # Check if columns exist
    if 'iptm' not in df.columns or 'ptm' not in df.columns:
        return None, f"Skipping {csv_path.name}: Missing iptm/ptm columns"
        
    # Assuming iptm/ptm are constant for the file, check first row
    if df.empty:
        return None, None
        
    global_iptm = df['iptm'].iloc[0]
    global_ptm = df['ptm'].iloc[0]
    
    if global_iptm < args.iptm_min or global_ptm < args.ptm_min:
        return None, None # File rejected based on global scores

    # 2. Row-wise filter (plddt, pae)
    # Ensure columns are numeric
    for col in ['plddt_A', 'plddt_B', 'pae']:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    
    df_filtered = df[
        (df['plddt_A'] >= args.plddt_min) &
        (df['plddt_B'] >= args.plddt_min) &
        (df['pae'] <= args.pae_max)
    ].copy()
    
    if df_filtered.empty:
        return None, None

    # 3. Motif Filter
    # Need to check chains separately? 
    # The file contains pairs. We want to see if the interface forms a motif on at least one side?
    # Or both sides? Usually a binding site on one protein is clustered.
    # Let's check if there is a valid motif on chain A OR chain B.
    
    # Get unique residues for A and B in the filtered set
    res_A = sorted(df_filtered['res_A'].unique())
    res_B = sorted(df_filtered['res_B'].unique())
    
    has_motif_A, clusters_A = find_motifs(res_A, args.motif_gap, args.motif_min_res)
    has_motif_B, clusters_B = find_motifs(res_B, args.motif_gap, args.motif_min_res)
    
    if has_motif_A or has_motif_B:
        # Create summary
        summary = {
            "file": csv_path.name,
            "iptm": global_iptm,
            "ptm": global_ptm,
            "contacts": len(df_filtered),
            "mean_plddt_A": df_filtered["plddt_A"].mean(),
            "mean_plddt_B": df_filtered["plddt_B"].mean(),
            "mean_pae": df_filtered["pae"].mean(),
            "motifs": []
        }
        if has_motif_A:
            for cls in clusters_A:
                if len(cls) >= args.motif_min_res:
                    summary["motifs"].append(f"Chain A: {cls[0]}-{cls[-1]} (n={len(cls)})")
        if has_motif_B:
            for cls in clusters_B:
                if len(cls) >= args.motif_min_res:
                    summary["motifs"].append(f"Chain B: {cls[0]}-{cls[-1]} (n={len(cls)})")
                    
        return df_filtered, summary
    else:
        return None, None

def main():
    args = parse_args()
    
    input_path = Path(args.input)
    os.makedirs(args.outdir, exist_ok=True)
    
    files_to_process = []
    if input_path.is_dir():
        files_to_process = sorted(list(input_path.glob("*.csv")))
    elif input_path.is_file():
        files_to_process = [input_path]
    else:
        print(f"Input not found: {input_path}")
        sys.exit(1)
        
    print(f"Found {len(files_to_process)} files to process.")
    print(f"Filters: pLDDT >= {args.plddt_min}, PAE <= {args.pae_max}, ipTM >= {args.iptm_min}, pTM >= {args.ptm_min}")
    print(f"Motif: >= {args.motif_min_res} residues within gap {args.motif_gap}")
    
    count_saved = 0
    for csv_file in files_to_process:
        result_df, summary = process_file(csv_file, args)
        if isinstance(summary, str): # Error message
             print(summary)
             continue
             
        if result_df is not None:
            out_name = csv_file.name.replace(".csv", "_filtered.csv")
            out_path = os.path.join(args.outdir, out_name)
            result_df.to_csv(out_path, index=False)
            
            print(f"--- {csv_file.name} ---")
            print(f"Global Scores: ipTM={summary['iptm']:.3f}, pTM={summary['ptm']:.3f}")
            print(f"Retained Contacts: {summary['contacts']}")
            print(f"Mean Interface Scores: pLDDT_A={summary['mean_plddt_A']:.1f}, pLDDT_B={summary['mean_plddt_B']:.1f}, PAE={summary['mean_pae']:.2f}")
            print("Motifs Found:")
            for m in summary['motifs']:
                print(f"  - {m}")
            print("-" * 40)
            
            count_saved += 1
            
    print(f"Finished. Saved {count_saved} filtered interaction files to {args.outdir}/")

if __name__ == "__main__":
    main()
