#!/bin/bash

# Input file with one line per set of coordinates
INPUT_FILE=$1
OUTPUT_FILE="test.txt"

# Run the program, feed coordinates, and save output
cmsRun queryField.py < "$INPUT_FILE" > "$OUTPUT_FILE"

cat "$OUTPUT_FILE"