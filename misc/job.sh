#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage: $0 <csv_file> [output_name]"
    echo ""
    echo "Examples:"
    echo "  $0 domain_decomposition_results.csv"
    echo "  $0 results.csv my_plot.png"
    exit 1
fi

CSV_FILE="$1"

if [ ! -f "$CSV_FILE" ]; then
    echo "Error: File '$CSV_FILE' not found!"
    exit 1
fi

if [ $# -ge 2 ]; then
    OUTPUT_NAME="$2"
else
    # Auto-generate from input filename
    BASENAME=$(basename "$CSV_FILE" .csv)
    OUTPUT_NAME="${BASENAME}_plot.png"
fi

# Extract title from filename (replace underscores with spaces, capitalize)
TITLE=$(basename "$CSV_FILE" .csv | sed 's/_/ /g' | sed 's/\b\(.\)/\u\1/g')

# Create Python script
python3 << EOF
from plotting_utilities import plot_1d_results, compute_statistics, plot_2d_results
import pandas as pd
import sys

csv_file = "$CSV_FILE"
output_name = "$OUTPUT_NAME"
title = "$TITLE"

# Auto-detect 1D or 2D
try:
    df = pd.read_csv(csv_file)

    if 'y' in df.columns:
        print(f"Detected 2D data in {csv_file}")
        print(f"Plotting 2D results...")
        plot_2d_results(csv_file, title=title, save_name=output_name)
    else:
        print(f"Detected 1D data in {csv_file}")
        print(f"Plotting 1D results...")
        plot_1d_results(csv_file, title=title, figsize=(12, 8), save_name=output_name)
        print("")
        compute_statistics(csv_file)

    print(f"\nPlot saved to: {output_name}")

except Exception as e:
    print(f"Error: {e}", file=sys.stderr)
    sys.exit(1)
EOF

# Check if plotting succeeded
if [ $? -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "Success!"
    echo "=========================================="
    echo "Input:  $CSV_FILE"
    echo "Output: $OUTPUT_NAME"
else
    echo "Error: Plotting failed!" >&2
    exit 1
fi