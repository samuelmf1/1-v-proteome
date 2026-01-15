#!/bin/bash
# Setup Boltz cache by downloading weights and unpacking CCD
set -e

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <cache_dir>"
    exit 1
fi

CACHE_DIR=$1

# Create dummy input to trigger initial download/setup
cat <<EOF > dummy.yaml
version: 1
sequences:
  - protein:
      id: A
      sequence: A
EOF

mkdir -p home_dir

# We want to populate the HOST cache_dir.
CACHE_PATH=$(readlink -f "$CACHE_DIR")
ln -s "$CACHE_PATH" home_dir/.boltz

export HOME="$PWD/home_dir"

# Run dummy prediction to force download of CCD and (incorrectly) v1 weights
echo "Triggering Boltz download (CCD + user weights)..."
boltz predict dummy.yaml --cache "$CACHE_PATH" --use_msa_server --out_dir results || true

# Manually download CORRECT Boltz-2 weights
echo "Downloading correct Boltz-2 weights from Hugging Face..."
wget -q -O "$CACHE_PATH/boltz2_conf.ckpt" "https://huggingface.co/boltz-community/boltz-2/resolve/main/boltz2_conf.ckpt"
wget -q -O "$CACHE_PATH/boltz2_aff.ckpt" "https://huggingface.co/boltz-community/boltz-2/resolve/main/boltz2_aff.ckpt"

# UNPACK CCD manually using standalone script
echo "Unpacking CCD..."
unpack_ccd.py "home_dir/.boltz/ccd.pkl" "home_dir/.boltz/mols"

touch boltz_ready
