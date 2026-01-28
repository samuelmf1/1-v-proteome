#!/usr/bin/env python3

import argparse
import gzip
import json
import os
import sys
import csv
import statistics
import yaml

def parse_args():
    parser = argparse.ArgumentParser(description="Generate AlphaFold 3 input JSONs for Proteome Screening.")
    parser.add_argument("--targets", help="Path to targets TSV (id\\ttype\\tsequence)")
    parser.add_argument("--target_id", help="Single target ID")
    parser.add_argument("--target_type", help="Single target type (protein, dna, rna)")
    parser.add_argument("--target_seq", help="Single target sequence")
    parser.add_argument("--seqs", help="Path to STRING sequences fasta.gz (optional)")
    parser.add_argument("--outdir", required=True, help="Output directory for JSON files")
    parser.add_argument("--filetype", choices=['json', 'yaml'], default='json', help="Output file format (default: json)")
    return parser.parse_args()

def load_targets(filepath):
    targets = []
    with open(filepath, 'r') as f:
        # Check for header
        first_line = f.readline()
        f.seek(0)
        
        has_header = False
        if 'id' in first_line.lower() and 'type' in first_line.lower() and 'sequence' in first_line.lower():
            has_header = True
            
        reader = csv.reader(f, delimiter='\t')
        if has_header:
            next(reader, None) # Skip header
            
        for row in reader:
            if len(row) >= 3:
                # Assuming order if no header or just taking first 3 columns
                # The user specified: name(id) type sequence
                targets.append({
                    "id": row[0],
                    "type": row[1],
                    "sequence": row[2]
                })
    return targets

def get_seq_template(name, target_obj, string_id, string_seq):
    # Entity A: The custom target (from TSV)
    # Entity B: The STRING protein
    
    # Determine type for Entity A
    type_a = target_obj['type'] # protein, dna, rna
    # Entity B is always protein (STRING)
    
    return {
        "name": name,
        "sequences": [
            {type_a: {"id": ["A"], "sequence": target_obj['sequence']}},
            {"protein": {"id": ["B"], "sequence": string_seq}}
        ],
        "modelSeeds": [1],
        "dialect": "alphafold3",
        "version": 1
    }

def main():
    args = parse_args()
    
    # 1. Load Targets
    if args.target_id and args.target_type and args.target_seq:
        targets = [{
            "id": args.target_id,
            "type": args.target_type,
            "sequence": args.target_seq
        }]
        print(f"Using single target: {args.target_id}")
    elif args.targets:
        targets = load_targets(args.targets)
        print(f"Loaded {len(targets)} targets from {args.targets}")
    else:
        print("Error: Either --targets or (--target_id, --target_type, --target_seq) must be provided.")
        sys.exit(1)
    
    # 2. Setup Output
    os.makedirs(args.outdir, exist_ok=True)
    
    
    # 3. Stream STRING sequences and generate JSONs IF seqs provided
    if args.seqs:
        lengths_file = os.path.join(args.outdir, "af3_lengths.tsv")
        print(f"Streaming sequences from {args.seqs}...")
        
        count = 0
        len_a_list = []
        len_b_list = []
        len_sum_list = []
        
        with open(lengths_file, 'w', newline='') as csvfile:
            fieldnames = ['target_id', 'string_id', 'len_a', 'len_b', 'len_sum']
            writer = csv.DictWriter(csvfile, delimiter='\t', fieldnames=fieldnames)
            writer.writeheader()
            
            # Precompute sanitized target IDs
            sanitized_targets = []
            for t in targets:
                t_id_clean = t['id'].replace('/', '_').replace(' ', '_')
                sanitized_targets.append((t, t_id_clean))

            xopen = gzip.open if args.seqs.endswith('.gz') else open
            with xopen(args.seqs, 'rt') as f:
                current_id = None
                current_seq = []
                
                def process_sequence(sid, sseq):
                    nonlocal count
                    if not sid or not sseq: return
                    
                    # Sanitize STRING ID once
                    s_id_clean = sid.replace('/', '_').replace(' ', '_')

                    # Pair this STRING protein with EVERY target
                    for t, t_id_clean in sanitized_targets:
                        if args.filetype == 'json':
                            filename = f"{t_id_clean}__{s_id_clean}.json"
                            filepath = os.path.join(args.outdir, filename)
                            
                            data = get_seq_template(f"{t_id_clean}__{s_id_clean}", t, sid, sseq)
                            
                            with open(filepath, 'w') as out:
                                json.dump(data, out, indent=2)
                        
                        elif args.filetype == 'yaml':
                            filename = f"{t_id_clean}__{s_id_clean}.yaml"
                            filepath = os.path.join(args.outdir, filename)
                            
                            # Construct YAML data structure (simplified)
                            # Scalar IDs, minimal fields
                            type_a = t['type']
                            
                            # Entity A (Target)
                            id_a = "A"
                            # The user example showed D for DNA. Let's try to match that if type is DNA?
                            if type_a == 'dna':
                                id_a = "D" # Trying to match user style
                            
                            user_entity = {
                                type_a: {
                                    "id": id_a,
                                    "sequence": t['sequence']
                                }
                            }
                            
                            # Entity B (STRING Protein)
                            string_entity = {
                                "protein": {
                                    "id": "B" if id_a != "B" else "C",
                                    "sequence": sseq
                                }
                            }
                            
                            # Combine
                            yaml_data = {
                                "version": 1,
                                "sequences": [user_entity, string_entity]
                            }
                            
                            # Determine Header
                            if type_a == 'protein':
                                header = "# ppi_complex.yaml\n"
                            elif type_a == 'dna':
                                header = "# protein_dna_complex.yaml\n"
                            elif type_a == 'rna':
                                header = "# protein_rna_complex.yaml\n"
                            else:
                                header = "# complex.yaml\n"
                                
                            with open(filepath, 'w') as out:
                                out.write(header)
                                yaml.dump(yaml_data, out, sort_keys=False)
                        
                        # Track lengths
                        len_a = len(t['sequence'])
                        len_b = len(sseq)
                        len_sum = len_a + len_b
                        
                        writer.writerow({
                            'target_id': t_id_clean,
                            'string_id': s_id_clean,
                            'len_a': len_a,
                            'len_b': len_b,
                            'len_sum': len_sum
                        })
                        
                        len_a_list.append(len_a)
                        len_b_list.append(len_b)
                        len_sum_list.append(len_sum)
                        
                        count += 1
                        
                        if count % 10000 == 0:
                            print(f"Generated {count} files...")

                for line in f:
                    if line.startswith('>'):
                        if current_id:
                            process_sequence(current_id, "".join(current_seq))
                        current_id = line.strip()[1:]
                        current_seq = []
                    else:
                        current_seq.append(line.strip())
                
                # Process last sequence
                if current_id:
                    process_sequence(current_id, "".join(current_seq))
                
        print(f"Done. Generated {count} input files in {args.outdir}.")

        # 4. Generate Summary Statistics
        if len_a_list:
            stats_file = os.path.join(args.outdir, "af3_stats.txt")
            with open(stats_file, 'w') as f:
                f.write(f"Sequence Statistics (N={count})\n")
                f.write("-" * 30 + "\n")
                f.write(f"Mean Length A:   {statistics.mean(len_a_list):.2f}\n")
                f.write(f"Median Length A: {statistics.median(len_a_list):.2f}\n")
                f.write("-" * 30 + "\n")
                f.write(f"Mean Length B:   {statistics.mean(len_b_list):.2f}\n")
                f.write(f"Median Length B: {statistics.median(len_b_list):.2f}\n")
                f.write("-" * 30 + "\n")
                f.write(f"Mean Sum (A+B):   {statistics.mean(len_sum_list):.2f}\n")
                f.write(f"Median Sum (A+B): {statistics.median(len_sum_list):.2f}\n")
            print(f"Summary stats written to {stats_file}")

    else:
        # Generate files for targets ONLY
        print("Generating input files for targets only (monomers)...")
        count = 0
        for t in targets:
            t_id_clean = t['id'].replace('/', '_').replace(' ', '_')
            
            if args.filetype == 'yaml':
                filename = f"{t_id_clean}.yaml"
                filepath = os.path.join(args.outdir, filename)
                
                type_a = t['type']
                id_a = "A"
                if type_a == 'dna':
                    id_a = "D"
                
                user_entity = {
                    type_a: {
                        "id": id_a,
                        "sequence": t['sequence']
                    }
                }
                
                yaml_data = {
                    "version": 1,
                    "sequences": [user_entity]
                }
                
                header = f"# {type_a}_monomer.yaml\n"
                with open(filepath, 'w') as out:
                    out.write(header)
                    yaml.dump(yaml_data, out, sort_keys=False)
                    
            elif args.filetype == 'json':
                # Simplified JSON for monomers? Not requested but good to have
                filename = f"{t_id_clean}.json"
                filepath = os.path.join(args.outdir, filename)
                
                type_a = t['type']
                data = {
                    "name": t_id_clean,
                    "sequences": [
                        {type_a: {"id": ["A"], "sequence": t['sequence']}}
                    ],
                    "modelSeeds": [1],
                    "dialect": "alphafold3",
                    "version": 1
                }
                with open(filepath, 'w') as out:
                    json.dump(data, out, indent=2)
            
            count += 1
            
        print(f"Done. Generated {count} monomer input files for targets in {args.outdir}.")

if __name__ == "__main__":
    main()
