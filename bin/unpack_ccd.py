#!/usr/bin/env python3

import pickle
import os
import sys
from rdkit import Chem
from pathlib import Path

# Ensure properties (like atom names) are pickled!
try:
    Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)
except:
    # Older RDKit versions might use 0xFFFF
    Chem.SetDefaultPickleProperties(0xFFFF)

def unpack_ccd(ccd_path, mols_dir):
    ccd_path = Path(ccd_path)
    mols_dir = Path(mols_dir)

    if not ccd_path.exists():
        print(f"Error: CCD not found at {ccd_path}")
        sys.exit(1)

    if mols_dir.exists() and (mols_dir / "ALA.pkl").exists():
        print("Mols directory seems populated. Skipping unpack.")
        return

    print(f"Unpacking {ccd_path} to {mols_dir}")
    mols_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        with open(ccd_path, "rb") as f:
            data = pickle.load(f)
        
        if isinstance(data, dict):
            count = 0
            for name, mol in data.items():
                # Sanitize name just in case
                safe_name = name.replace("/", "_")
                with open(mols_dir / f"{safe_name}.pkl", "wb") as out:
                    pickle.dump(mol, out)
                count += 1
            print(f"Successfully unpacked {count} molecules.")
        else:
            print(f"Error: CCD is not a dictionary: {type(data)}")
            sys.exit(1)
    except Exception as e:
        print(f"Error: Failed to unpack CCD: {e}")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: unpack_ccd.py <ccd.pkl> <mols_dir>")
        sys.exit(1)
    
    unpack_ccd(sys.argv[1], sys.argv[2])
