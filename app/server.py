#!/usr/bin/env python3
import http.server
import socketserver
import os
import glob
import json
import urllib.parse
import gzip
import argparse

DEFAULT_PORT = 8001
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
SEQ_FILES_DIR = os.path.abspath("/gpfs/commons/home/sfriedman/projects/vlp/seq_files")
FILTERED_INTERACTIONS_DIR = os.path.join(SEQ_FILES_DIR, "filtered_interactions")

# Load protein info map
PROTEIN_INFO_FILE = os.path.join(BASE_DIR, "9606.protein.info.v12.0.txt.gz")
PROTEIN_MAP = {}

# Load GO Data
GO_ONTOLOGY_FILE = os.path.join(BASE_DIR, "go-plus.json.gz")
GO_ANNOTATION_FILE = os.path.join(BASE_DIR, "goa_human.gaf.gz")
GO_DEFINITIONS = {} # GO_ID -> { name, def }
GO_ANNOTATIONS = {} # Gene Symbol -> [ { go_id, aspect, qualifier } ]

# Load STRING Data
STRING_FILE = os.path.join(BASE_DIR, "9606.protein.physical.links.full.v12.0.txt.gz")
STRING_INTERACTIONS = {} # ProteinID -> [ { partner, score, experiments, database } ]

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
                        
                        # Filter for reasonable confidence or just store all?
                        # User wants physical links. The file is "physical.links".
                        # Let's verify if we should filter. Combined score > 400 is medium confidence.
                        # Let's filter > 150 to remove very low confidence noise but keep most.
                        if combined_score < 150: continue

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

class Handler(http.server.SimpleHTTPRequestHandler):
    def get_target_results_dir(self, file_id):
        """Derive the results directory from the target name in the file_id."""
        # file_id format: target__proteinID
        parts = file_id.split("__")
        if len(parts) > 0:
            target = parts[0]
            return os.path.join(SEQ_FILES_DIR, f"{target}_results")
        return None # Should not happen with valid file_ids

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
        try:
            # Pattern: *_interactions_filtered.csv
            pattern = os.path.join(FILTERED_INTERACTIONS_DIR, "*_interactions_filtered.csv")
            files = glob.glob(pattern)
            results = []
            
            for f in files:
                basename = os.path.basename(f)
                # Extract ID: remove _interactions_filtered.csv
                # e.g. psPAX2_RT_ribonucleaseH__9606.ENSP00000322801_interactions_filtered.csv
                if basename.endswith("_interactions_filtered.csv"):
                    file_id = basename[:-len("_interactions_filtered.csv")]
                    
                    # Extract Chain Names
                    # Format: prefix__9606.ENSP...
                    # Chain A = prefix
                    # Chain B = Preferred Name (or ID)
                    
                    parts = file_id.split("__")
                    chainA_name = parts[0] if len(parts) > 0 else "Chain A"
                    chainB_name = "Chain B" # Default
                    
                    display_name = file_id # Fallback
                    protein_info = None

                    if len(parts) > 1:
                        protein_id = parts[-1] 
                        if protein_id in PROTEIN_MAP:
                            info = PROTEIN_MAP[protein_id]
                            preferred = info["name"]
                            chainB_name = preferred
                            display_name = f"{chainA_name} vs {preferred} ({protein_id})"
                            
                            # Attach full info
                            protein_info = {
                                "id": protein_id,
                                "name": info["name"],
                                "size": info["size"],
                                "annotation": info["annotation"]
                            }
                        else:
                            chainB_name = protein_id
                            display_name = f"{chainA_name} vs {protein_id}"
                            
                    
                            # Minimal info fallback
                            protein_info = {
                                "id": protein_id,
                                "name": protein_id,
                                "size": "N/A",
                                "annotation": "No annotation found."
                            }
                    
                    # Get modification time (based on CIF file)
                    target_results_dir = self.get_target_results_dir(file_id)
                    cif_filename = f"{file_id}_model_0.cif"
                    cif_path = os.path.join(target_results_dir, file_id, f"boltz_results_{file_id}", "predictions", file_id, cif_filename)
                    if os.path.exists(cif_path):
                        mtime = os.path.getmtime(cif_path)
                    else:
                        mtime = os.path.getmtime(f) # Fallback to CSV mtime
                    
                    # Get Confidence Score
                    confidence = 0
                    json_filename = f"confidence_{file_id}_model_0.json"
                    # Path matched to handle_metrics_lookup
                    json_path = os.path.join(target_results_dir, file_id, f"boltz_results_{file_id}", "predictions", file_id, json_filename)
                    
                    if os.path.exists(json_path):
                        try:
                            with open(json_path, 'r') as jf:
                                jdata = json.load(jf)
                                confidence = jdata.get('confidence_score', 0)
                        except:
                            pass

                    results.append({
                        "id": file_id,
                        "name": display_name,
                        "chainA": chainA_name,
                        "chainB": chainB_name,
                        "proteinInfo": protein_info,
                        "mtime": mtime,
                        "confidence": confidence
                    })
            
            # Sort by name
            results.sort(key=lambda x: x["name"])
            
            self.send_response(200)
            self.send_header('Content-type', 'application/json')
            self.end_headers()
            self.wfile.write(json.dumps(results).encode())
            
        except Exception as e:
            self.send_error(500, str(e))

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
             if protein_id in STRING_INTERACTIONS:
                 interactions = STRING_INTERACTIONS[protein_id]
                 
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

        filename = f"{file_id}_interactions_filtered.csv"
        filepath = os.path.join(FILTERED_INTERACTIONS_DIR, filename)
        
        if os.path.exists(filepath):
            self.serve_file_content(filepath, 'text/csv')
        else:
            self.send_error(404, "CSV not found")

    def serve_monomer_cif(self, target_id):
        # Prevent traversal
        if '..' in target_id or '/' in target_id:
            self.send_error(400, "Invalid ID")
            return
            
        # Monomer Path: SEQ_FILES_DIR/monomers/results/<target>/boltz_results_<target>/predictions/<target>/<target>_model_0.cif
        monomers_dir = os.path.join(SEQ_FILES_DIR, "monomers", "results")
        
        filename = f"{target_id}_model_0.cif"
        path = os.path.join(monomers_dir, target_id, f"boltz_results_{target_id}", "predictions", target_id, filename)
        
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

        # CIF path logic:
        # seq_files/{target}_results/<ID>/boltz_results_<ID>/predictions/<ID>/<ID>_model_0.cif
        
        # Check specific directory structure
        # NOTE: User request said: "find the corresponding cif file from seq_files/psPAX2_RT_ribonucleaseH_results"
        # Based on previous `list_dir`, the structure is:
        # psPAX2_RT_ribonucleaseH_results/<ID>/boltz_results_<ID>/predictions/<ID>/...
        
        cif_filename = f"{file_id}_model_0.cif"
        
        # Construct path deeply
        # Path segment 1: ID
        path_seg1 = file_id
        # Path segment 2: boltz_results_ID
        path_seg2 = f"boltz_results_{file_id}"
        # Path segment 3: predictions
        path_seg3 = "predictions"
        # Path segment 4: ID
        path_seg4 = file_id
        
        target_results_dir = self.get_target_results_dir(file_id)
        full_path = os.path.join(target_results_dir, path_seg1, path_seg2, path_seg3, path_seg4, cif_filename)
        
        if os.path.exists(full_path):
             self.serve_file_content(full_path, 'text/plain') # Serving CIF as text
        else:
             print(f"DEBUG: Failed to find {full_path}")
             self.send_error(404, f"CIF not found at expected path {full_path}")

    def handle_metrics_lookup(self, file_id):
        # Path logic: same as CIF
        json_filename = f"confidence_{file_id}_model_0.json"
        
        path_seg1 = file_id
        path_seg2 = f"boltz_results_{file_id}"
        path_seg3 = "predictions"
        path_seg4 = file_id
        
        target_results_dir = self.get_target_results_dir(file_id)
        full_path = os.path.join(target_results_dir, path_seg1, path_seg2, path_seg3, path_seg4, json_filename)
        
        if os.path.exists(full_path):
            self.serve_file_content(full_path, 'application/json')
        else:
            print(f"DEBUG: Failed to find {full_path}")
            self.send_error(404, f"Metrics JSON not found")

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

        # Path: SEQ_FILES_DIR/monomers/sparrow_results/<target>_sparrow_features_residues.tsv
        filename = f"{target_id}_sparrow_features_residues.tsv"
        filepath = os.path.join(SEQ_FILES_DIR, "monomers", "sparrow_results", filename)
        
        if not os.path.exists(filepath):
            # Try falling back to just returning empty if file doesn't exist, or 404
            # Returning empty JSON is safer for the frontend to handle "no data"
            self.send_response(200)
            self.send_header('Content-type', 'application/json')
            self.end_headers()
            self.wfile.write(b"{}")
            return

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
