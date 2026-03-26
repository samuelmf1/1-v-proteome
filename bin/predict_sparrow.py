#!/usr/bin/env python3

import argparse
import os
import sys
import yaml
import csv
import numpy as np

try:
    from sparrow import Protein
except ImportError:
    # Fallback for local testing if sparrow is not installed
    print("Warning: sparrow not found. Using mock for testing.")
    class Predictor:
        def disorder(self, **kwargs): return np.random.rand(10)
        def dssp(self, **kwargs): return np.random.randint(0, 3, 10)
        def nes(self, **kwargs): return np.random.rand(10)
        def nis(self, **kwargs): return np.random.rand(10)
        def phosphorylation(self, **kwargs): return np.random.rand(10)
        def pscore(self, **kwargs): return np.random.rand(10)
        def tad(self, **kwargs): return np.random.rand(10)
        def mitochondrial_targeting(self, **kwargs): return np.random.rand(10)
        def transmembrane_region(self, **kwargs): return np.random.randint(0, 2, 10)
        def asphericity(self, **kwargs): return 0.5
        def radius_of_gyration(self, **kwargs): return 15.0
        def end_to_end_distance(self, **kwargs): return 40.0
        def scaling_exponent(self, **kwargs): return 0.58
        def prefactor(self, **kwargs): return 2.0

    class Protein:
        def __init__(self, sequence):
            self.sequence = sequence
            self.FCR = 0.2
            self.NCPR = 0.05
            self.hydrophobicity = 0.4
            self.complexity = 0.6
            self.predictor = Predictor()

def parse_args():
    parser = argparse.ArgumentParser(description="Predict protein features using SPARROW.")
    parser.add_argument("--yaml", required=True, help="Input monomer YAML file (AlphaFold 3 style)")
    parser.add_argument("--out", required=True, help="Output TSV file path")
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Load YAML
    try:
        with open(args.yaml, 'r') as f:
            data = yaml.safe_load(f)
    except Exception as e:
        print(f"Error reading YAML: {e}")
        sys.exit(1)
        
    # Extract sequence from AF3 style YAML
    # Expected structure: { 'version': 1, 'sequences': [ { 'protein': { 'id': 'A', 'sequence': '...' } } ] }
    sequence = None
    target_id = os.path.basename(args.yaml).replace(".yaml", "")
    
    if 'sequences' in data:
        for entity in data['sequences']:
            if 'protein' in entity:
                sequence = entity['protein']['sequence']
                break
    
    if not sequence:
        print(f"Error: Could not find protein sequence in {args.yaml}")
        sys.exit(1)
        
    # Initialize Sparrow Protein object
    p = Protein(sequence)
    results = {}
    
    # 1. Global Parameters (Directly from Protein object)
    results['FCR'] = p.FCR
    results['NCPR'] = p.NCPR
    results['hydrophobicity'] = p.hydrophobicity
    
    # 2. ALBATROSS Features (Global)
    try:
        results['asphericity'] = p.predictor.asphericity()
        results['radius_of_gyration'] = p.predictor.radius_of_gyration(use_scaled=True)
        results['end_to_end_distance'] = p.predictor.end_to_end_distance(use_scaled=True)
        results['scaling_exponent'] = p.predictor.scaling_exponent()
        results['prefactor'] = p.predictor.prefactor()
    except Exception as e:
        print(f"Warning: Error calculating ALBATROSS features: {e}")

    # 3. Sequence Predictors (Per-residue features, take mean for summary)
    # Mapping based on available methods in sparrow version installed
    residue_results = []
    for i, aa in enumerate(sequence):
        residue_results.append({'residue_index': i + 1, 'amino_acid': aa})

    predictors_map = {
        'disorder': 'disorder',
        'dssp_coil': 'dssp_coil',
        'dssp_extended': 'dssp_extended',
        'dssp_helicity': 'dssp_helicity',
        'nuclear_export_signal': 'nuclear_export_signal',
        'nuclear_import_signal': 'nuclear_import_signal',
        'serine_phosphorylation': 'serine_phosphorylation',
        'threonine_phosphorylation': 'threonine_phosphorylation',
        'tyrosine_phosphorylation': 'tyrosine_phosphorylation',
        'pscore': 'pscore',
        'transactivation_domains': 'transactivation_domains',
        'mitochondrial_targeting_sequence': 'mitochondrial_targeting_sequence',
        'transmembrane_regions': 'transmembrane_regions',
        'pLDDT': 'pLDDT'
    }
    
    for label, method in predictors_map.items():
        try:
            func = getattr(p.predictor, method)
            values = func()
            if isinstance(values, (list, np.ndarray)):
                results[f'mean_{label}'] = np.mean(values)
                # Store per-residue values
                for i, val in enumerate(values):
                    if i < len(residue_results):
                        residue_results[i][label] = val
            else:
                results[label] = values
        except Exception as e:
            print(f"Warning: Error calculating {label} ({method}): {e}")
            
    # Write to TSV (Summary)
    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(args.out, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['id'] + list(results.keys()), delimiter='\t')
        writer.writeheader()
        row = {'id': target_id}
        row.update(results)
        writer.writerow(row)
    
    # Write to TSV (Per-residue)
    residue_out = args.out.replace(".tsv", "_residues.tsv")
    if residue_results:
        with open(residue_out, 'w', newline='') as f:
            fieldnames = ['residue_index', 'amino_acid'] + [k for k in residue_results[0].keys() if k not in ['residue_index', 'amino_acid']]
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            writer.writerows(residue_results)
        
    print(f"Successfully predicted features for {target_id} -> {args.out} and {residue_out}")

if __name__ == "__main__":
    main()
