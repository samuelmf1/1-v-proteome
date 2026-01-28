#!/usr/bin/env python3
# interactions.py
import json
import numpy as np
import pandas as pd
from pathlib import Path
from itertools import combinations
from Bio.PDB import MMCIFParser, NeighborSearch
import os

def find_boltz_files(base_dir):
    """Locate the necessary files within the Boltz results structure."""
    base_path = Path(base_dir)
    # The structure: results/ID/boltz_results_ID/predictions/ID/
    pred_dirs = list(base_path.glob("**/predictions/*"))
    if not pred_dirs:
        raise FileNotFoundError(f"Could not find predictions folder inside {base_dir}")
    
    # Sort by modification time to get the latest prediction if multiple exist
    pred_path = sorted(pred_dirs, key=lambda x: x.stat().st_mtime)[-1]
    sample_id = pred_path.name
    
    try:
        files = {
            "cif": next(pred_path.glob(f"*{sample_id}*_model_0.cif")),
            "pae": next(pred_path.glob(f"pae_*{sample_id}*_model_0.npz")),
            "plddt": next(pred_path.glob(f"plddt_*{sample_id}*_model_0.npz")),
            "json": next(pred_path.glob(f"confidence_*{sample_id}*_model_0.json")),
            "outdir": pred_path
        }
    except StopIteration:
        raise FileNotFoundError(f"Missing one or more required files (.cif, .npz, .json) in {pred_path}")
        
    return files, sample_id

def get_npz_data(filepath):
    """Safely extract data from npz by looking for Boltz keys or fallback to first key."""
    data = np.load(filepath)
    # Boltz usually uses the filename prefix as the key (e.g., 'pae' or 'plddt')
    possible_keys = ['pae', 'plddt', 'arr_0']
    for k in possible_keys:
        if k in data.files:
            return data[k]
    # Final fallback: just give me the first thing in the archive
    return data[data.files[0]]

def load_scores(files):
    """Load and format confidence metrics with robust key detection."""
    with open(files["json"]) as f:
        summary = json.load(f)
    
    plddt_scores = get_npz_data(files["plddt"])
    # Auto-scale pLDDT if it's in [0, 1] range
    if plddt_scores.max() <= 1.0:
        plddt_scores = plddt_scores * 100.0

    return {
        "plddt": plddt_scores,
        "pae": get_npz_data(files["pae"]),
        "iptm": summary.get("iptm", 0.0),
        "ptm": summary.get("ptm", 0.0)
    }

def get_residue_mapping(model):
    """Map BioPython residues to the flat array indices used by Boltz."""
    res_list = list(model.get_residues())
    res_to_idx = {res.get_full_id(): i for i, res in enumerate(res_list)}
    return res_to_idx

def generate_pymol_script(df, cif_file, out_path):
    """Creates a PyMOL script to highlight interaction hotspots."""
    with open(out_path, 'w') as f:
        f.write(f"load {cif_file.name}\n")
        f.write("color gray80, all\n")
        f.write("show cartoon, all\n")
        
        # Color residues by hotspot score
        for _, row in df.iterrows():
            if row['hotspot_A'] > 2:
                f.write(f"color orange, chain {row['chain_A']} and resi {row['res_A']}\n")
                f.write(f"show sticks, chain {row['chain_A']} and resi {row['res_A']}\n")
            if row['hotspot_B'] > 2:
                f.write(f"color marine, chain {row['chain_B']} and resi {row['res_B']}\n")
                f.write(f"show sticks, chain {row['chain_B']} and resi {row['res_B']}\n")
        
        f.write("set stick_radius, 0.3\n")
        f.write("zoom\n")

def analyze_all_interfaces(files, scores, args):
    """Find contacts between all unique chain combinations."""
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("complex", files["cif"])
    model = structure[0]
    res_to_idx = get_residue_mapping(model)
    
    chains = [chain.id for chain in model]
    chain_pairs = list(combinations(chains, 2))
    
    ns = NeighborSearch(list(model.get_atoms()))
    all_close_pairs = ns.search_all(args.dist)
    
    all_data = []

    for chA, chB in chain_pairs:
        seen_pairs = set()
        pair_data = []
        
        for atom1, atom2 in all_close_pairs:
            res1, res2 = atom1.get_parent(), atom2.get_parent()
            id1, id2 = res1.get_parent().id, res2.get_parent().id
            
            if {id1, id2} == {chA, chB}:
                rA, rB = (res1, res2) if id1 == chA else (res2, res1)
                if (rA.get_full_id(), rB.get_full_id()) in seen_pairs:
                    continue
                
                idx_A = res_to_idx[rA.get_full_id()]
                idx_B = res_to_idx[rB.get_full_id()]
                
                pA = scores["plddt"][idx_A] if scores["plddt"].ndim == 1 else scores["plddt"][idx_A].mean()
                pB = scores["plddt"][idx_B] if scores["plddt"].ndim == 1 else scores["plddt"][idx_B].mean()
                pae_val = scores["pae"][idx_A, idx_B]

                if pA >= args.plddt_min and pB >= args.plddt_min and pae_val <= args.pae_max:
                    pair_data.append({
                        "chain_A": chA, "res_A": rA.get_id()[1], "res_A_name": rA.get_resname(),
                        "chain_B": chB, "res_B": rB.get_id()[1], "res_B_name": rB.get_resname(),
                        "plddt_A": round(float(pA), 2), "plddt_B": round(float(pB), 2),
                        "pae": round(float(pae_val), 2), "dist": round(float(atom1 - atom2), 3)
                    })
                    seen_pairs.add((rA.get_full_id(), rB.get_full_id()))
        
        if pair_data:
            df_pair = pd.DataFrame(pair_data)
            df_pair['hotspot_A'] = df_pair['res_A'].map(df_pair['res_A'].value_counts())
            df_pair['hotspot_B'] = df_pair['res_B'].map(df_pair['res_B'].value_counts())
            df_pair['iptm'] = round(scores['iptm'], 2)
            df_pair['ptm'] = round(scores['ptm'], 2)
            all_data.append(df_pair)

    return pd.concat(all_data) if all_data else pd.DataFrame()

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", type=str, required=True, help="Top-level Boltz result directory")
    parser.add_argument("--dist", type=float, default=4.0)
    parser.add_argument("--plddt_min", type=float, default=70.0)
    parser.add_argument("--pae_max", type=float, default=10.0)
    parser.add_argument("--outdir", type=str, default="interactions")
    args = parser.parse_args()

    try:
        files, sample_id = find_boltz_files(args.dir)
        scores = load_scores(files)
        
        print(f"--- Analysis Summary for {sample_id} ---")
        print(f"Global ipTM: {scores['iptm']:.3f}")
        print("-" * 45)
        
        results = analyze_all_interfaces(files, scores, args)
        
        if not results.empty:
            summary = results.groupby(['chain_A', 'chain_B']).size().reset_index(name='Contact_Count')
            print(summary.to_string(index=False))
            print("-" * 45)
            
            os.makedirs(args.outdir, exist_ok=True)
            csv_path = os.path.join(args.outdir, f"{sample_id}_interactions.csv")
            results.sort_values(["chain_A", "chain_B", "pae"]).to_csv(csv_path, index=False)
            
            pml_path = os.path.join(args.outdir, f"view_{sample_id}_interactions.pml")
            generate_pymol_script(results, files["cif"], pml_path)
            
            print(f"CSV results: {os.path.basename(csv_path)}")
            print(f"PyMOL script: {os.path.basename(pml_path)}")
            print(f"Total High-Conf Contacts: {len(results)}")
        else:
            print("No high-confidence interactions found with current thresholds.")
            
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()