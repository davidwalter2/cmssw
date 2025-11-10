#!/bin/bash

# Input file with one line per set of coordinates
INPUT_FILE=$1 # "coords_Run1.txt"
# INPUT_FILE="coords_Run2.txt"

# Output file
# OUTPUT_FILE="field_results_170812_sx5.txt"
# OUTPUT_FILE="field_results_170812_run1.txt"
# OUTPUT_FILE="field_results_170812_run2.txt"
# OUTPUT_FILE="field_results_160812_run1.txt"
# OUTPUT_FILE="field_results_160812_run2.txt"

# OUTPUT_FILE="field_results_120812_run1.txt"
# OUTPUT_FILE="field_results_120812_run2.txt"
# OUTPUT_FILE="field_results_130503_run1.txt"
# OUTPUT_FILE="field_results_130503_run2.txt"

# OUTPUT_FILE="field_results_polyFit2D.txt"
# OUTPUT_FILE="field_results_polyFit3D.txt"

OUTPUT_FILE="test.txt"

# Run the program, feed coordinates, and save output
cmsRun queryField.py < "$INPUT_FILE" > "$OUTPUT_FILE"

cat "$OUTPUT_FILE"