#!/bin/bash

# Define the root directory containing your subfolders
INPUT_DIR="$SCRATCH/alphafast"

# Check if directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Directory $INPUT_DIR not found."
    exit 1
fi

echo "Scanning subfolders in $INPUT_DIR..."
echo "---------------------------------------"

# Loop through each item in the input directory
for TARGET_DIR in "$INPUT_DIR"/*/; do
    # Ensure it's a directory
    [ -d "$TARGET_DIR" ] || continue

    GOAL=$(ls "$TARGET_DIR/inputs" 2>/dev/null | wc -l | tr -d ' ' )
    if [ -z "$GOAL" ] || [ "$GOAL" -eq 0 ]; then GOAL=19699; fi # Fallback

    # Count the occurrences of completion.txt TERMS_OF_USE.md
    COUNT=$(find "$TARGET_DIR" -name "completion.txt" | wc -l)

    # Calculate percentage using awk for floating point math
    PERCENT=$(awk "BEGIN {printf \"%.1f\", ($COUNT/$GOAL)*100}")

    # Generate the ASCII bar (1 # for every 10%)
    # Using printf to repeat the '#' character
    BAR_SIZE=$(( COUNT * 10 / GOAL ))
    BAR=$(printf "%${BAR_SIZE}s" | tr ' ' '#')

    # Output the result
    printf "%-20s %-10s %d/%d (%s%%)\n" "$(basename "$TARGET_DIR")" "$BAR" "$COUNT" "$GOAL" "$PERCENT"
done
