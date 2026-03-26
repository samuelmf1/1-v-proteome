import os
import json

MONOMERS_DIR = "/gpfs/commons/projects/one-v-proteome/monomers"

def build_targets_cache():
    af3_dir = os.path.join(MONOMERS_DIR, "af3")
    if os.path.exists(af3_dir):
        targets = [d for d in os.listdir(af3_dir) if os.path.isdir(os.path.join(af3_dir, d))]
        targets.sort()
        return targets
    return []

targets = build_targets_cache()
print(f"Targets found: {len(targets)}")
print(json.dumps(targets, indent=2))
