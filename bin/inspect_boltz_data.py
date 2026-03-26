#!/usr/bin/env python3
import numpy as np
import json
import glob
import sys
from pathlib import Path

def inspect(base_dir):
    print(f"Inspecting {base_dir}")
    path = Path(base_dir)
    # Replicate finding logic
    pred_dirs = list(path.glob("**/predictions/*"))
    if not pred_dirs:
        print("No prediction dirs found")
        return
    
    pred_path = sorted(pred_dirs, key=lambda x: x.stat().st_mtime)[-1]
    print(f"Prediction path: {pred_path}")
    
    sample_id = pred_path.name
    
    plddt_file = next(pred_path.glob(f"plddt_*{sample_id}*_model_0.npz"), None)
    pae_file = next(pred_path.glob(f"pae_*{sample_id}*_model_0.npz"), None)
    json_file = next(pred_path.glob(f"confidence_*{sample_id}*_model_0.json"), None)
    
    if plddt_file:
        print(f"\nLoading pLDDT: {plddt_file}")
        data = np.load(plddt_file)
        print(f"  Keys: {list(data.files)}")
        for k in data.files:
            arr = data[k]
            print(f"  Key '{k}': shape={arr.shape}, min={arr.min():.2f}, max={arr.max():.2f}, mean={arr.mean():.2f}")
            if arr.ndim == 1:
                print(f"  First 10 values: {arr[:10]}")
    else:
        print("pLDDT file not found")
        
    if pae_file:
        print(f"\nLoading PAE: {pae_file}")
        data = np.load(pae_file)
        print(f"  Keys: {list(data.files)}")
        for k in data.files:
            arr = data[k]
            print(f"  Key '{k}': shape={arr.shape}, min={arr.min():.2f}, max={arr.max():.2f}, mean={arr.mean():.2f}")
    else:
        print("PAE file not found")
        
    if json_file:
        print(f"\nLoading JSON: {json_file}")
        with open(json_file) as f:
            j = json.load(f)
            print(f"  Content: {j}")
            
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: script.py <dir>")
        sys.exit(1)
    inspect(sys.argv[1])
