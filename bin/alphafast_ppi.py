#!/usr/bin/env python3
import argparse
import sys
import json
from pathlib import Path

try:
    import pandas as pd
    from Bio.PDB import MMCIFParser, NeighborSearch
except ImportError:
    print("Error: Missing dependencies. Please run: pip install pandas biopython")
    sys.exit(1)

def analyze_interfaces(cif_path, json_path, dist_cutoff, plddt_min, pae_max, iptm_min):
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("complex", str(cif_path))
    model = structure[0]

    data = {}
    iptm_val = None
    
    if json_path and json_path.exists():
        with open(json_path) as f:
            data = json.load(f)
        
        # AF3 often puts metrics in 'summary' or top-level
        iptm_val = data.get("iptm")
        if iptm_val is None and "summary" in data:
            iptm_val = data["summary"].get("iptm")
            
    # Filter by global ipTM if threshold is set and data exists
    if iptm_val is not None and iptm_val < iptm_min:
        print(f"Skipping: ipTM ({iptm_val}) is below threshold ({iptm_min})")
        return pd.DataFrame()

    atom_plddts = data.get("atom_plddts", [])
    pae_matrix  = data.get("pae", [])

    all_atoms = list(model.get_atoms())
    res_list  = list(model.get_residues())
    res_to_token_idx = {res.get_full_id(): i for i, res in enumerate(res_list)}

    ns = NeighborSearch(all_atoms)
    all_close_pairs = ns.search_all(dist_cutoff)

    pair_data = []
    seen_pairs = set()

    for atom1, atom2 in all_close_pairs:
        res1, res2 = atom1.get_parent(), atom2.get_parent()
        chain1, chain2 = res1.get_parent().id, res2.get_parent().id

        if chain1 != chain2:
            r1_fid, r2_fid = res1.get_full_id(), res2.get_full_id()
            rA, rB = (res1, res2) if r1_fid < r2_fid else (res2, res1)
            pair_key = (rA.get_full_id(), rB.get_full_id())

            if pair_key in seen_pairs:
                continue

            pae_val = None
            p1, p2 = None, None

            if pae_matrix and res_to_token_idx:
                idxA = res_to_token_idx.get(rA.get_full_id())
                idxB = res_to_token_idx.get(rB.get_full_id())
                if idxA is not None and idxB is not None:
                    pae_val = pae_matrix[idxA][idxB]

            if atom_plddts:
                try:
                    at_idx1 = all_atoms.index(atom1)
                    at_idx2 = all_atoms.index(atom2)
                    p1, p2 = atom_plddts[at_idx1], atom_plddts[at_idx2]
                except (ValueError, IndexError):
                    pass

            # Metric Filtering
            keep = True
            if p1 is not None and p2 is not None:
                if p1 < plddt_min or p2 < plddt_min:
                    keep = False
            if pae_val is not None and pae_val > pae_max:
                keep = False

            if keep:
                pair_data.append({
                    "chain_A": rA.get_parent().id,
                    "res_A": rA.get_id()[1],
                    "res_name_A": rA.get_resname(),
                    "chain_B": rB.get_parent().id,
                    "res_B": rB.get_id()[1],
                    "res_name_B": rB.get_resname(),
                    "dist": round(atom1 - atom2, 3),
                    "pae": round(pae_val, 2) if pae_val is not None else "N/A",
                    "plddt_atom1": round(p1, 2) if p1 is not None else "N/A",
                    "plddt_atom2": round(p2, 2) if p2 is not None else "N/A",
                    "iptm": iptm_val if iptm_val is not None else "N/A"
                })
                seen_pairs.add(pair_key)

    return pd.DataFrame(pair_data)

def main():
    parser = argparse.ArgumentParser(description="Analyze AF3 interfaces with ipTM filtering.")
    parser.add_argument("cif", type=str, help="Path to the AF3 .cif file")
    parser.add_argument("--dist", type=float, default=5.0, help="Distance cutoff (Å)")
    parser.add_argument("--plddt_min", type=float, default=60.0, help="Min per-atom pLDDT")
    parser.add_argument("--pae_max", type=float, default=20.0, help="Max per-residue PAE")
    parser.add_argument("--iptm_min", type=float, default=0.0, help="Min global ipTM to save output")
    parser.add_argument("--outdir", type=str, default=".", help="Output directory")

    args = parser.parse_args()
    
    cif_path = Path(args.cif)
    json_path = cif_path.parent / cif_path.name.replace("_model.cif", "_confidences.json")
    if not json_path.exists():
        json_path = cif_path.with_suffix('.json')
    
    df = analyze_interfaces(cif_path, json_path if json_path.exists() else None, 
                            args.dist, args.plddt_min, args.pae_max, args.iptm_min)

    if not df.empty:
        outdir = Path(args.outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        out_file = outdir / f"{cif_path.stem}_interfaces.csv"
        df.to_csv(out_file, index=False)
        print(f"Success: {len(df)} pairs saved to {out_file}")

if __name__ == "__main__":
    main()

