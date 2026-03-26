#!/usr/bin/env python3
import sys
import os
import yaml
import hashlib
import argparse
from pathlib import Path
import logging

# Add local ColabFold fork to path
COLABFOLD_PATH = "/gpfs/commons/home/sfriedman/projects/fork/ColabFold"
if COLABFOLD_PATH not in sys.path:
    sys.path.insert(0, COLABFOLD_PATH)

try:
    from colabfold.mmseqs.search import mmseqs_search_monomer
except ImportError as e:
    import traceback
    traceback.print_exc()
    print(f"Error: Could not import colabfold from {COLABFOLD_PATH}. Exception: {e}", file=sys.stderr)
    sys.exit(1)

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

DB_BASE = Path("/gpfs/commons/home/sfriedman/projects/fork/ColabFold/colabfold_db")
MSA_DIR = Path("vlp/msa")

def get_sequence_hash(sequence):
    """Returns a SHA256 hash of the sequence (stripped)."""
    clean_seq = sequence.strip().upper()
    return hashlib.sha256(clean_seq.encode('utf-8')).hexdigest()

def ensure_msa(sequence, msa_dir, db_base, temp_dir, threads=8, sensitivity=8.0, use_env=True, db_load_mode=2):
    """
    Checks if MSA exists for the sequence. If not, generates it using colabfold.
    Returns the path to the MSA file.
    """
    seq_hash = get_sequence_hash(sequence)
    msa_path = msa_dir / f"{seq_hash}.a3m"
    
    if msa_path.exists():
        logger.info(f"MSA for sequence (hash: {seq_hash[:8]}) already exists at {msa_path}. Reuse.")
        return str(msa_path.resolve())
    
    logger.info(f"Generating MSA for sequence (hash: {seq_hash[:8]}) with {threads} threads...")
    
    # Create a temporary input file for this sequence
    temp_query = temp_dir / f"{seq_hash}.fasta"
    with open(temp_query, 'w') as f:
        f.write(f">{seq_hash}\n{sequence}\n")
    
    # Create a temporary output dir for this search
    search_out = temp_dir / f"search_{seq_hash}"
    search_out.mkdir(exist_ok=True, parents=True)
    
    # Create qdb for mmseqs
    # Replicating logic from colabfold/mmseqs/search.py
    query_file = search_out / "query.fas"
    with open(query_file, 'w') as f:
        f.write(f">{seq_hash}\n{sequence}\n")
    
    # Run createdb
    import subprocess
    cmd = ["mmseqs", "createdb", str(query_file), str(search_out / "qdb"), "--shuffle", "0", "--dbtype", "1"]
    
    try:
        subprocess.check_call(cmd)
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to create qdb: {e}")
        return None

    try:
        # Run mmseqs search
        mmseqs_search_monomer(
            dbbase=db_base,
            base=search_out,
            mmseqs=Path("mmseqs"),
            threads=threads,
            s=sensitivity,
            use_env=use_env,
            use_templates=False,
            db_load_mode=db_load_mode
        )
        
        # Look for any .a3m file in search_out
        potential_a3ms = list(search_out.glob("*.a3m"))
        result_a3m = None
        for p in potential_a3ms:
            if p.name == f"{seq_hash}.a3m":
                result_a3m = p
                break
        
        if not result_a3m and potential_a3ms:
            result_a3m = potential_a3ms[0]
            
        if result_a3m and result_a3m.exists():
            # Atomic move to final MSA dir
            # Use a temporary name in the same filesystem for atomicity
            temp_final_path = msa_path.with_suffix(".tmp")
            os.replace(result_a3m, temp_final_path)
            os.replace(temp_final_path, msa_path)
            logger.info(f"MSA generated and saved to {msa_path}")
        else:
            logger.error(f"MSA generation failed. Expected output {result_a3m} or any .a3m not found in {search_out}")
            return None
            
    except Exception as e:
        logger.error(f"Error generating MSA: {e}")
        return None

    return str(msa_path.resolve())

def process_yaml(yaml_file, output_yaml, msa_dir, db_base, temp_dir, threads=8, sensitivity=8.0, use_env=True, db_load_mode=2):
    with open(yaml_file, 'r') as f:
        data = yaml.safe_load(f)
    
    if 'sequences' not in data:
        return
    
    modified = False
    
    for entry in data['sequences']:
        for key, val in entry.items():
            if key == 'protein' and 'sequence' in val:
                seq = val['sequence']
                msa_path = ensure_msa(seq, msa_dir, db_base, temp_dir, threads=threads, sensitivity=sensitivity, use_env=use_env, db_load_mode=db_load_mode)
                if msa_path:
                    val['msa'] = msa_path
                    modified = True
                else:
                    logger.error(f"Failed to get MSA for a sequence in {yaml_file}")
    
    # Always write output YAML, even if not modified (so prediction step has its input)
    with open(output_yaml, 'w') as f:
        yaml.dump(data, f)
    if modified:
        logger.info(f"Written updated YAML to {output_yaml}")

def main():
    parser = argparse.ArgumentParser(description="Generate MSAs for Boltz YAML inputs using local ColabFold DB.")
    parser.add_argument("yaml_files", nargs='+', help="Input YAML files")
    parser.add_argument("--output_dir", default="msa_yamls", help="Directory to save updated YAMLs")
    parser.add_argument("--msa_dir", default="vlp/msa", help="Directory to save/cache A3M files")
    parser.add_argument("--threads", type=int, default=8, help="Threads per mmseqs search")
    parser.add_argument("--skip_existing_dir", help="Check this directory for existing PPI results and skip if found.")
    parser.add_argument("--sensitivity", type=float, default=8.0, help="MMseqs2 sensitivity (-s). Default 8.0. Lower (e.g. 3.0-6.0) is faster.")
    parser.add_argument("--no_env", action="store_true", help="Skip searching the environmental database (faster).")
    parser.add_argument("--db_load_mode", type=int, default=2, help="Database load mode (0: auto, 1: fread, 2: mmap, 3: mmap+touch).")
    
    args = parser.parse_args()
    
    msa_dir = Path(args.msa_dir).resolve()
    msa_dir.mkdir(parents=True, exist_ok=True)
    
    out_dir = Path(args.output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    
    temp_dir = Path("temp_msa_gen").resolve()
    temp_dir.mkdir(exist_ok=True)
    
    for yaml_file in args.yaml_files:
        base_name = os.path.basename(yaml_file).replace('.yaml', '')
        
        # Skip logic
        if args.skip_existing_dir:
            if os.path.isdir(os.path.join(args.skip_existing_dir, base_name)):
                logger.info(f"SKIPPING: {base_name} (PPI result already exists in {args.skip_existing_dir})")
                # Fall through to process_yaml to ensure output file is created for Nextflow.
                # It will be very fast due to MSA caching.

        out_yaml = out_dir / f"{base_name}.yaml"
        logger.info(f"Processing {yaml_file} -> {out_yaml}")
        process_yaml(yaml_file, out_yaml, msa_dir, DB_BASE, temp_dir, 
                     threads=args.threads, 
                     sensitivity=args.sensitivity, 
                     use_env=not args.no_env,
                     db_load_mode=args.db_load_mode)

if __name__ == "__main__":
    main()
