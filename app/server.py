#!/usr/bin/env python3
import http.server
import socketserver
import os
import glob
import json
import urllib.parse
import gzip
import argparse
import csv

DEFAULT_PORT = 8001
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# Data root directory (absolute path)
# Update this path if moving the app away from the data directory
DATA_ROOT = "/gpfs/commons/projects/one-v-proteome"

# Data directories (absolute paths)
DATASETS_DIR = os.path.join(BASE_DIR, "datasets")
SEQ_FILES_DIR = os.path.join(DATA_ROOT, "ppi")
MONOMERS_DIR = os.path.join(DATA_ROOT, "monomers")
MANIFEST_FILE = os.path.join(DATA_ROOT, "ppi_manifest.csv")

# Load protein info map
PROTEIN_INFO_FILE = os.path.join(DATASETS_DIR, "9606.protein.info.v12.0.txt.gz")
PROTEIN_MAP = {}

# Load GO Data
GO_ONTOLOGY_FILE = os.path.join(DATASETS_DIR, "go-plus.json.gz")
GO_ANNOTATION_FILE = os.path.join(DATASETS_DIR, "goa_human.gaf.gz")
GO_DEFINITIONS = {} # GO_ID -> { name, def }
GO_ANNOTATIONS = {} # Gene Symbol -> [ { go_id, aspect, qualifier } ]

# Load STRING Data
STRING_FILE = os.path.join(DATASETS_DIR, "9606.protein.physical.links.full.v12.0.txt.gz")
STRING_INTERACTIONS = {} # ProteinID -> [ { partner, score, experiments, database } ]
MANIFEST_ENTRIES = {}
MANIFEST_ENTRIES_CASEFOLD = {}

def normalize_target_id(target_name):
    return (target_name or "").strip().lower()

def find_first_existing_path(candidates):
    for candidate in candidates:
        if candidate and os.path.exists(candidate):
            return candidate
    return None

def get_manifest_structure_dir(row):
    structure_dir = row.get("structure_dir", "").strip()
    if not structure_dir:
        return None
    if os.path.isabs(structure_dir):
        return structure_dir
    return os.path.abspath(os.path.join(DATA_ROOT, structure_dir))

def get_structure_file_candidates(entry):
    structure_name = entry["structure_name"]
    structure_dir = entry["structure_dir"]
    return {
        "cif": [
            os.path.join(structure_dir, f"{structure_name}_model.cif"),
        ],
        "interfaces": [
            os.path.join(structure_dir, f"{structure_name}_model_interfaces.csv"),
            os.path.join(structure_dir, "interfaces.csv"),
        ],
        "metrics": [
            os.path.join(structure_dir, f"{structure_name}_summary_confidences.json"),
            os.path.join(structure_dir, f"{structure_name}_confidences.json"),
            os.path.join(structure_dir, f"confidence_{structure_name}_model.json"),
            os.path.join(structure_dir, f"{structure_name}_confidence.json"),
            os.path.join(structure_dir, f"confidence_{structure_name}.json"),
        ],
    }

def get_manifest_entry(file_id):
    entry = MANIFEST_ENTRIES.get(file_id)
    if entry is None:
        entry = MANIFEST_ENTRIES_CASEFOLD.get(file_id.lower())
    return entry

def calculate_interface_confidence(csv_path, file_id):
    if not csv_path or not os.path.exists(csv_path):
        return 0.0

    try:
        with open(csv_path, newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            total_plddt = 0.0
            count = 0
            for row in reader:
                try:
                    total_plddt += float(row['plddt_atom1']) + float(row['plddt_atom2'])
                    count += 2
                except (KeyError, TypeError, ValueError):
                    continue

            if count > 0:
                return (total_plddt / count) / 100.0
    except Exception as e:
        print(f"Warning: Could not calculate confidence for {file_id}: {e}")

    return 0.0

def load_manifest():
    print("Loading PPI manifest...")
    if not os.path.exists(MANIFEST_FILE):
        print(f"Warning: PPI manifest not found at {MANIFEST_FILE}")
        return

    try:
        with open(MANIFEST_FILE, newline='') as manifest_file:
            reader = csv.DictReader(manifest_file)
            for row in reader:
                structure_name = row.get("structure_name", "").strip()
                structure_dir = get_manifest_structure_dir(row)
                if not structure_name or not structure_dir:
                    continue

                entry = dict(row)
                entry["structure_name"] = structure_name
                entry["structure_dir"] = structure_dir
                entry["target_protein"] = row.get("target_protein", "").strip()
                entry["interacting_protein"] = row.get("interacting_protein", "").strip()

                MANIFEST_ENTRIES[structure_name] = entry
                MANIFEST_ENTRIES_CASEFOLD[structure_name.lower()] = entry

        print(f"Loaded {len(MANIFEST_ENTRIES)} manifest entries.")
    except Exception as e:
        print(f"Error loading PPI manifest: {e}")

def load_string_data():
    print("Loading STRING Data...")
    if os.path.exists(STRING_FILE):
        try:
            with gzip.open(STRING_FILE, 'rt') as f:
                header = next(f) # Skip header
                # Header: protein1 protein2 homology experiments experiments_transferred database database_transferred textmining textmining_transferred combined_score
                for line in f:
                    parts = line.split()
                    if len(parts) >= 10:
                        p1 = parts[0]
                        p2 = parts[1]
                        
                        # Parse scores
                        experiments = int(parts[3])
                        database = int(parts[5])
                        combined_score = int(parts[9])
                        
                        # Filter for high-confidence interactions only (≥700)
                        if combined_score < 700: continue

                        if p1 not in STRING_INTERACTIONS: STRING_INTERACTIONS[p1] = []
                        # if p2 not in STRING_INTERACTIONS: STRING_INTERACTIONS[p2] = [] # User only cares about target interactions, usually symmetric but let's just load p1->p2
                        
                        # Store essential info
                        interaction = {
                            "partner": p2,
                            "score": combined_score,
                            "experiments": experiments,
                            "database": database
                        }
                        STRING_INTERACTIONS[p1].append(interaction)
                        
        except Exception as e:
            print(f"Error loading STRING data: {e}")
    else:
        print("Warning: STRING file not found.")

def load_go_data():
    print("Loading GO Ontology...")
    try:
        with gzip.open(GO_ONTOLOGY_FILE, 'rt') as f:
            data = json.load(f)
            # The structure is data['graphs'][0]['nodes']
            for graph in data.get('graphs', []):
                for node in graph.get('nodes', []):
                    go_id = node.get('id')
                    if go_id and go_id.startswith('http://purl.obolibrary.org/obo/GO_'):
                        # Convert URI to ID format: GO:XXXXXXX
                        short_id = go_id.split('/')[-1].replace('_', ':')
                        
                        label = node.get('lbl', 'No label')
                        definition = "No definition"
                        if 'meta' in node and 'definition' in node['meta']:
                            definition = node['meta']['definition'].get('val', 'No definition')
                        
                        GO_DEFINITIONS[short_id] = {
                            "name": label,
                            "def": definition
                        }
    except Exception as e:
        print(f"Error loading GO ontology: {e}")

    print("Loading GO Annotations...")
    try:
        if os.path.exists(GO_ANNOTATION_FILE):
             with gzip.open(GO_ANNOTATION_FILE, 'rt') as f:
                for line in f:
                    if line.startswith('!'): continue
                    parts = line.split('\t')
                    if len(parts) >= 15:
                        # GAF 2.1/2.2 Format:
                        # Col 2: DB Object Symbol (e.g. NUDT4B)
                        # Col 4: GO ID (e.g. GO:0003723)
                        # Col 8: Aspect (F, P, C)
                        symbol = parts[2]
                        go_id = parts[4]
                        aspect = parts[8]
                        
                        if symbol not in GO_ANNOTATIONS:
                            GO_ANNOTATIONS[symbol] = []
                        
                        # Avoid duplicates
                        existing = [x['id'] for x in GO_ANNOTATIONS[symbol]]
                        if go_id not in existing:
                            GO_ANNOTATIONS[symbol].append({
                                "id": go_id,
                                "aspect": aspect
                            })
    except Exception as e:
         print(f"Error loading GO annotations: {e}")

# Load Protein Info
if os.path.exists(PROTEIN_INFO_FILE):
    print("Loading protein info...")
    try:
        with gzip.open(PROTEIN_INFO_FILE, 'rt') as f:
            for line in f:
                if line.startswith("#"): continue
                parts = line.split('\t')
                # string_protein_id	preferred_name	protein_size	annotation
                if len(parts) >= 4:
                    protein_id = parts[0].strip()
                    PROTEIN_MAP[protein_id] = {
                        "name": parts[1].strip(),
                        "size": parts[2].strip(),
                        "annotation": parts[3].strip()
                    }
    except Exception as e:
        print(f"Error loading protein info: {e}")
else:
    print("Warning: Protein info file not found.")

# Load GO data on startup
load_go_data()
load_string_data()
load_manifest()

# Cache file list at startup
CACHED_FILES = []
CACHED_FILES_JSON = b"[]"
CACHED_TARGETS = []
CACHED_TARGETS_JSON = b"[]"

def build_targets_cache():
    """Build the target cache from manifest entries with monomer support."""
    global CACHED_TARGETS, CACHED_TARGETS_JSON
    print("Building targets cache from manifest...")
    
    af3_dir = os.path.join(MONOMERS_DIR, "af3")
    if not os.path.exists(af3_dir):
        print(f"Warning: {af3_dir} not found.")
        return

    available_monomers = {
        d for d in os.listdir(af3_dir)
        if os.path.isdir(os.path.join(af3_dir, d))
    }

    manifest_targets = set()
    for entry in MANIFEST_ENTRIES.values():
        target_name = entry.get("target_protein", "") or entry.get("structure_name", "")
        if "__" in target_name:
            target_name = target_name.split("__", 1)[0]
        normalized_target = normalize_target_id(target_name)
        if normalized_target:
            manifest_targets.add(normalized_target)

    targets = sorted(manifest_targets & available_monomers)
    CACHED_TARGETS = targets
    CACHED_TARGETS_JSON = json.dumps(targets).encode()

def build_files_cache():
    """Build the interaction cache from manifest-backed AF3 outputs."""
    global CACHED_FILES, CACHED_FILES_JSON
    print("Building file list cache...")

    results = []

    for entry in MANIFEST_ENTRIES.values():
        file_id = entry["structure_name"]
        parts = file_id.split("__")
        target_name = entry["target_protein"] or (parts[0] if parts else "Chain A")
        chainA_name = normalize_target_id(target_name) or "unknown"
        chainB_name = "Chain B"

        display_name = file_id
        protein_info = None

        protein_id = entry["interacting_protein"] or (parts[-1] if len(parts) > 1 else "")
        if protein_id:
            protein_id_upper = protein_id.upper()
            if protein_id_upper in PROTEIN_MAP:
                info = PROTEIN_MAP[protein_id_upper]
                preferred = info["name"]
                chainB_name = preferred
                display_name = f"{chainA_name} vs {preferred} ({protein_id_upper})"
                
                protein_info = {
                    "id": protein_id_upper,
                    "name": info["name"],
                    "size": info["size"],
                    "annotation": info["annotation"]
                }
            else:
                chainB_name = protein_id
                display_name = f"{chainA_name} vs {protein_id}"
                
                protein_info = {
                    "id": protein_id,
                    "name": protein_id,
                    "size": "N/A",
                    "annotation": "No annotation found."
                }

        candidates = get_structure_file_candidates(entry)
        cif_path = find_first_existing_path(candidates["cif"])
        interfaces_path = find_first_existing_path(candidates["interfaces"])
        metrics_path = find_first_existing_path(candidates["metrics"])

        mtime_source = cif_path or interfaces_path or metrics_path or entry["structure_dir"]
        try:
            mtime = os.path.getmtime(mtime_source)
        except OSError:
            mtime = 0.0

        confidence = calculate_interface_confidence(interfaces_path, file_id)
        if confidence == 0.0:
            try:
                confidence = float(entry.get("ranking_score", "") or 0.0)
            except ValueError:
                confidence = 0.0

        other_targets = entry.get("other_target_proteins", "").strip()

        results.append({
            "id": file_id,
            "name": display_name,
            "chainA": chainA_name,
            "chainB": chainB_name,
            "proteinInfo": protein_info,
            "mtime": mtime,
            "confidence": confidence,
            "other_target_proteins": other_targets
        })
    
    results.sort(key=lambda x: x["name"])
    CACHED_FILES = results
    CACHED_FILES_JSON = json.dumps(results).encode()
    print(f"Cached {len(results)} interaction files.")

build_targets_cache()
build_files_cache()

class Handler(http.server.SimpleHTTPRequestHandler):
    def get_pair_entry(self, file_id):
        return get_manifest_entry(file_id)

    def log_message(self, format, *args):
        """Override to log to file as well as stderr."""
        # Standard formatting
        message = "%s - - [%s] %s\n" % (
            self.client_address[0],
            self.log_date_time_string(),
            format % args
        )
        
        # Write to server.log
        try:
            with open(os.path.join(os.getcwd(), "server.log"), "a") as f:
                f.write(message)
        except Exception:
            pass # Don't crash on logging error
            
        # Default behavior (console/stderr)
        super().log_message(format, *args)

    def do_GET(self):
        parsed_path = urllib.parse.urlparse(self.path)
        path = parsed_path.path
        
        if path == '/api/files':
            self.handle_files_list()
        elif path == '/api/targets':
            self.handle_targets_list()
        elif path.startswith('/data/csv/'):
            file_id = path.split('/')[-1]
            self.serve_csv(file_id)
        elif path.startswith('/data/cif/'):
            file_id = path.split('/')[-1]
            self.serve_cif(file_id)
        elif path.startswith('/data/monomer_cif/'):
            target_id = path.split('/')[-1]
            self.serve_monomer_cif(target_id)
        elif path.startswith('/api/go/'):
            symbol = path.split('/')[-1]
            self.handle_go_lookup(symbol)
        elif path.startswith('/api/go_batch'):
            self.handle_go_batch_lookup()
        elif path.startswith('/api/string/'):
            protein_id = path.split('/')[-1]
            self.handle_string_lookup(protein_id)
        elif path.startswith('/api/metrics/'):
            file_id = path.split('/')[-1]
            self.handle_metrics_lookup(file_id)
        elif path.startswith('/api/disorder/'):
            target_id = path.split('/')[-1]
            self.handle_disorder_lookup(target_id)
        else:
            # Serve local files (index.html etc)
            super().do_GET()

    def do_POST(self):
        parsed_path = urllib.parse.urlparse(self.path)
        path = parsed_path.path
        if path == '/api/go_batch':
            self.handle_go_batch_lookup()
        else:
            self.send_error(404, "Not Found")

    def handle_files_list(self):
        """Serve the pre-cached file list instantly."""
        self.send_response(200)
        self.send_header('Content-type', 'application/json')
        self.send_header('Content-Length', str(len(CACHED_FILES_JSON)))
        self.end_headers()
        self.wfile.write(CACHED_FILES_JSON)

    def handle_targets_list(self):
        """Serve the pre-cached target list instantly."""
        self.send_response(200)
        self.send_header('Content-type', 'application/json')
        self.send_header('Content-Length', str(len(CACHED_TARGETS_JSON)))
        self.end_headers()
        self.wfile.write(CACHED_TARGETS_JSON)

    def handle_go_lookup(self, symbol):
        try:
            # Decode symbol (e.g. URL encoded)
            symbol = urllib.parse.unquote(symbol)
            
            results = self.get_go_data(symbol)
            
            self.send_response(200)
            self.send_header('Content-type', 'application/json')
            self.end_headers()
            self.wfile.write(json.dumps(results).encode())
            
        except Exception as e:
            self.send_error(500, str(e))

    def get_go_data(self, symbol):
        results = []
        if symbol in GO_ANNOTATIONS:
            for item in GO_ANNOTATIONS[symbol]:
                go_id = item['id']
                aspect = item['aspect']
                
                term_info = GO_DEFINITIONS.get(go_id, {"name": "Unknown", "def": "No definition found"})
                
                results.append({
                    "id": go_id,
                    "aspect": aspect,
                    "name": term_info["name"],
                    "definition": term_info["def"]
                })
        return results

    def handle_go_batch_lookup(self):
        try:
            content_length = int(self.headers['Content-Length'])
            post_data = self.rfile.read(content_length)
            symbols = json.loads(post_data.decode('utf-8'))
            
            results = {}
            for symbol in symbols:
                results[symbol] = self.get_go_data(symbol)
            
            self.send_response(200)
            self.send_header('Content-type', 'application/json')
            self.end_headers()
            self.wfile.write(json.dumps(results).encode())
        except Exception as e:
            self.send_error(500, str(e))

    def handle_string_lookup(self, protein_id):
        try:
             # Look for exact ID (e.g. 9606.ENSP00000480527)
             # Note: ID might come in without 9606 prefix or with it.
             # In protein list we see: 9606.ENSP...
             # Frontend usually has full ID.
             
             results = []
             # AF3 pipeline lowercases IDs; STRING keys are uppercase
             lookup_id = protein_id.upper() if protein_id not in STRING_INTERACTIONS else protein_id
             if lookup_id in STRING_INTERACTIONS:
                 interactions = STRING_INTERACTIONS[lookup_id]
                 
                 # Sort by score descending
                 interactions.sort(key=lambda x: x['score'], reverse=True)
                 
                 # Take top 50
                 top_interactions = interactions[:50]
                 
                 for item in top_interactions:
                     partner_id = item['partner']
                     # Try to get gene name
                     partner_name = partner_id
                     if partner_id in PROTEIN_MAP:
                         partner_name = PROTEIN_MAP[partner_id]['name']
                     elif partner_id.replace("9606.", "") in PROTEIN_MAP: # Try without prefix if map keys exclude it (Map keys seem to be full ID based on previous code)
                          # Code check: protein_id = parts[0].strip() from file. 
                          # Example file line: 9606.ENSP...
                          pass
                     
                     results.append({
                         "partner_id": partner_id,
                         "partner_name": partner_name,
                         "score": item['score'],
                         "experiments": item['experiments'],
                         "database": item['database']
                     })
             
             self.send_response(200)
             self.send_header('Content-type', 'application/json')
             self.end_headers()
             self.wfile.write(json.dumps(results).encode())

        except Exception as e:
            self.send_error(500, str(e))

    def serve_csv(self, file_id):
        # Prevent traversal
        if '..' in file_id or '/' in file_id:
            self.send_error(400, "Invalid ID")
            return

        entry = self.get_pair_entry(file_id)
        if entry is None:
            self.send_error(404, "Structure not found in manifest")
            return

        filepath = find_first_existing_path(get_structure_file_candidates(entry)["interfaces"])

        if filepath and os.path.exists(filepath):
            self.serve_file_content(filepath, 'text/csv')
        else:
            self.send_error(404, "CSV not found")

    def serve_monomer_cif(self, target_id):
        # Prevent traversal
        if '..' in target_id or '/' in target_id:
            self.send_error(400, "Invalid ID")
            return
            
        # Monomer Path: MONOMERS_DIR/af3/<target_id>/output/<target_id>/<target_id>_model.cif
        filename = f"{target_id}_model.cif"
        path = os.path.join(MONOMERS_DIR, "af3", target_id, "output", target_id, filename)
        
        if os.path.exists(path):
            self.serve_file_content(path, 'text/plain')
        else:
            print(f"DEBUG: Failed to find monomer at {path}")
            self.send_error(404, "Monomer CIF not found")

    def serve_cif(self, file_id):
        # Prevent traversal
        if '..' in file_id or '/' in file_id:
            self.send_error(400, "Invalid ID")
            return

        entry = self.get_pair_entry(file_id)
        if entry is None:
            self.send_error(404, "Structure not found in manifest")
            return

        full_path = find_first_existing_path(get_structure_file_candidates(entry)["cif"])

        if full_path and os.path.exists(full_path):
            self.serve_file_content(full_path, 'text/plain')
        else:
            print(f"DEBUG: Failed to find {full_path}")
            self.send_error(404, f"CIF not found at expected path {full_path}")

    def handle_metrics_lookup(self, file_id):
        entry = self.get_pair_entry(file_id)
        if entry is not None:
            full_path = find_first_existing_path(get_structure_file_candidates(entry)["metrics"])
            if full_path is not None:
                try:
                    with open(full_path, 'r') as f:
                        data = json.load(f)
                    
                    # Transform summary_confidences.json data to frontend format
                    result = {}
                    
                    # Map ranking_score to confidence_score
                    result['confidence_score'] = data.get('ranking_score', 0.0)
                    
                    # Direct mappings
                    result['ptm'] = data.get('ptm', 0.0)
                    result['iptm'] = data.get('iptm', 0.0)
                    result['fraction_disordered'] = data.get('fraction_disordered', 0.0)
                    result['has_clash'] = data.get('has_clash', 0.0)
                    
                    # Chain pTMs: convert array to object with chain labels
                    chain_ptm = data.get('chain_ptm', [])
                    if chain_ptm:
                        result['chains_ptm'] = {}
                        for i, val in enumerate(chain_ptm):
                            result['chains_ptm'][chr(65 + i)] = val  # A, B, C, etc.
                    
                    # Chain ipTMs
                    chain_iptm = data.get('chain_iptm', [])
                    if chain_iptm:
                        result['chains_iptm'] = {}
                        for i, val in enumerate(chain_iptm):
                            result['chains_iptm'][chr(65 + i)] = val  # A, B, C, etc.
                    
                    # Chain pair ipTM (if only 2 chains, this is the protein-protein interface score)
                    chain_pair_iptm = data.get('chain_pair_iptm', [])
                    if chain_pair_iptm and len(chain_pair_iptm) > 1:
                        # For a 2-chain complex, [0][1] is the A-B interaction
                        result['protein_iptm'] = chain_pair_iptm[0][1] if len(chain_pair_iptm[0]) > 1 else 0.0
                    
                    # Chain pair PAE min
                    chain_pair_pae = data.get('chain_pair_pae_min', [])
                    if chain_pair_pae and len(chain_pair_pae) > 1 and len(chain_pair_pae[0]) > 1:
                        result['interface_pae_min'] = chain_pair_pae[0][1]
                    
                    self.send_response(200)
                    self.send_header('Content-type', 'application/json')
                    content = json.dumps(result).encode()
                    self.send_header('Content-Length', str(len(content)))
                    self.end_headers()
                    self.wfile.write(content)
                    return
                except Exception as e:
                    print(f"Error parsing metrics for {file_id}: {e}")
        
        # No confidence file found — return empty JSON so frontend handles gracefully
        self.send_response(200)
        self.send_header('Content-type', 'application/json')
        self.end_headers()
        self.wfile.write(b"{}")

    def serve_file_content(self, filepath, content_type):
        try:
            with open(filepath, 'rb') as f:
                content = f.read()
            self.send_response(200)
            self.send_header('Content-type', content_type)
            self.send_header('Content-Length', str(len(content)))
            self.end_headers()
            self.wfile.write(content)
        except Exception as e:
            self.send_error(500, str(e))

    def handle_disorder_lookup(self, target_id):
        # Prevent traversal
        if '..' in target_id or '/' in target_id:
            self.send_error(400, "Invalid ID")
            return

        # Path: MONOMERS_DIR/sparrow/<target_id>_sparrow_features_residues.tsv
        # Use glob for case-insensitive matching as filenames like Cas9_... exist
        pattern = os.path.join(MONOMERS_DIR, "sparrow", f"{target_id}_sparrow_features_residues.tsv")
        matches = glob.glob(pattern, recursive=False)
        
        # If no direct match, try case-insensitive walk/search
        if not matches:
            sparrow_dir = os.path.join(MONOMERS_DIR, "sparrow")
            if os.path.exists(sparrow_dir):
                target_filename = f"{target_id}_sparrow_features_residues.tsv".lower()
                for f in os.listdir(sparrow_dir):
                    if f.lower() == target_filename:
                        matches = [os.path.join(sparrow_dir, f)]
                        break

        if not matches:
            # Returning empty JSON is safer for the frontend to handle "no data"
            self.send_response(200)
            self.send_header('Content-type', 'application/json')
            self.end_headers()
            self.wfile.write(b"{}")
            return

        filepath = matches[0]

        try:
            disorder_map = {}
            with open(filepath, 'r') as f:
                header = next(f).strip().split('\t')
                try:
                    # Find indices
                    idx_residue = header.index('residue_index')
                    idx_disorder = header.index('disorder')
                except ValueError:
                    print(f"Error: Missing columns in {filepath}")
                    self.send_error(500, "Invalid TSV format")
                    return

                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) > max(idx_residue, idx_disorder):
                        resi = parts[idx_residue]
                        try:
                            val = float(parts[idx_disorder])
                            disorder_map[resi] = val
                        except ValueError:
                            pass # Skip non-numeric

            self.send_response(200)
            self.send_header('Content-type', 'application/json')
            self.end_headers()
            self.wfile.write(json.dumps(disorder_map).encode())
            
        except Exception as e:
            print(f"Error serving disorder data: {e}")
            self.send_error(500, str(e))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PPI Visualization Server")
    parser.add_argument("--port", type=int, default=DEFAULT_PORT, help=f"Port to serve on (default: {DEFAULT_PORT})")
    args = parser.parse_args()

    # Enable address reuse to avoid "Address already in use" errors
    socketserver.TCPServer.allow_reuse_address = True
    
    with socketserver.TCPServer(("", args.port), Handler) as httpd:
        print(f"Serving at http://localhost:{args.port}")
        httpd.serve_forever()
