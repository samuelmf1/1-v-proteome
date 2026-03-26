#!/bin/bash

# Ensure we are in the diagram directory
cd "$(dirname "$0")"

# Create the outputs directory if it doesn't exist
mkdir -p outputs

echo "🎨 Beginning Mermaid-CLI Diagram Rendering..."

for file in *.mmd; do
    # Extract the base filename without the .mmd extension
    base="${file%.mmd}"
    
    echo "⏳ Rendering ${base}..."
    
    # 1. Generate standard vector graphic (SVG)
    npx -y @mermaid-js/mermaid-cli -i "$file" -o "outputs/${base}.svg"
    
    # 2. Generate high-resolution 4x scaled raster graphic (PNG)
    npx -y @mermaid-js/mermaid-cli -i "$file" -o "outputs/${base}.png" -s 4
    
done

echo "✅ All diagrams have been rendered into the diagram/outputs/ directory!"
