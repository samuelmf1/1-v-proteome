#!/usr/bin/env python3

import sys
import os
import yaml
import subprocess
import glob
from pathlib import Path

def process_yaml(yaml_path, out_dir):
    """
    Reads a YAML file, checks for mixed MSA usage, and modifies it if necessary.
    Runs boltz predict on the file.
    """
    try:
        with open(yaml_path, 'r') as f:
            data = yaml.safe_load(f)
    except Exception as e:
        print(f"Error reading {yaml_path}: {e}", file=sys.stderr)
        return False

    if 'sequences' not in data:
        print(f"Invalid YAML format in {yaml_path}: no 'sequences' key", file=sys.stderr)
        return False

    sequences = data['sequences']
    has_msa = []
    
    # Check MSA status for each sequence
    for seq_entry in sequences:
        # Each entry is a dict like {'protein': {...}} or just {...}
        # Based on observed YAML, it's a list of dicts where keys are types
        # e.g. - protein: {id: A, sequence: ..., msa: ...}
        
        # Flatten to get the inner dict
        inner = None
        for key, val in seq_entry.items():
            if isinstance(val, dict) and 'sequence' in val:
                inner = val
                break
        
        if inner:
            if 'msa' in inner:
                has_msa.append(True)
            else:
                has_msa.append(False)

    # Check for mixed state
    # Mixed if we have both True and False in the list
    is_mixed = any(has_msa) and not all(has_msa)
    
    modified_yaml_path = yaml_path
    
    if is_mixed:
        print(f"Detected mixed MSA usage in {yaml_path}. Stripping custom MSAs to force auto-generation.")
        
        # Strip MSAs
        new_sequences = []
        for seq_entry in sequences:
            new_entry = seq_entry.copy()
            for key, val in new_entry.items():
                if isinstance(val, dict) and 'msa' in val:
                    # Remove msa key
                    del val['msa']
            new_sequences.append(new_entry)
            
        data['sequences'] = new_sequences
        
        # Write to a temp file
        base_name = os.path.basename(yaml_path)
        modified_yaml_path = f"modified_{base_name}"
        
        with open(modified_yaml_path, 'w') as f:
            yaml.dump(data, f)

    # Prepare output directory
    base_name = os.path.basename(yaml_path).replace('.yaml', '')
    target_out_dir = os.path.join(out_dir, base_name)
    
    # Run Boltz
    # Command from predict_boltz.sh:
    # boltz predict "$yaml_file" --cache "home_dir/.boltz" --use_msa_server --num_workers 0 --out_dir "results/${base_name}"
    
    # We need to preserve the environment variables set by the wrapper script or assumed
    # The wrapper script sets HOME, etc.
    
    cmd = [
        "boltz", "predict", modified_yaml_path,
        "--cache", "home_dir/.boltz",
        "--use_msa_server",
        "--num_workers", "0",
        "--out_dir", target_out_dir
    ]
    
    print(f"Running: {' '.join(cmd)}")
    
    try:
        # Capture output to print it? Or let it flow to stdout/stderr
        # We want to see it in the logs
        subprocess.run(cmd, check=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Boltz prediction failed for {yaml_path} (exit code {e.returncode})", file=sys.stderr)
        return False
    finally:
        # Cleanup modified file if created
        if modified_yaml_path != yaml_path and os.path.exists(modified_yaml_path):
            os.remove(modified_yaml_path)

def main():
    if len(sys.argv) < 2:
        print("Usage: batch_predict_boltz.py <out_dir> <yaml_file1> [yaml_file2 ...]")
        sys.exit(1)
        
    out_dir = sys.argv[1]
    yaml_files = sys.argv[2:]
    
    print(f"Processing {len(yaml_files)} YAML files...")
    
    success_count = 0
    fail_count = 0
    
    for yf in yaml_files:
        if process_yaml(yf, out_dir):
            success_count += 1
        else:
            fail_count += 1
            
    print(f"Batch processing complete. Success: {success_count}, Failed: {fail_count}")
    
    if fail_count > 0:
        # Decide if we want to fail the whole batch or just report errors
        # If we fail, the pipeline stops. 
        # For now, let's exit with error if ANY fail, so we don't silently miss predictions.
        sys.exit(1)

if __name__ == "__main__":
    main()
