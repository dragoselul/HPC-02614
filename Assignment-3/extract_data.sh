#!/bin/bash

# Extract performance data from mm_batch_gpu output file
# Output format: Memory(KB) MFLOPS for gnuplot consumption

INPUT="${1:-mm_batch_gpu_*.out}"
OUTPUT="data.txt"

# Find the most recent output file if glob pattern used
if [[ "$INPUT" == *"*"* ]]; then
    INPUT=$(ls -t mm_batch_gpu_*.out 2>/dev/null | head -1)
    if [ -z "$INPUT" ]; then
        echo "Error: No mm_batch_gpu_*.out file found"
        exit 1
    fi
fi

echo "Processing: $INPUT"

# Clear output file
> $OUTPUT

# Methods to extract (in order matching gnuplot index)
# index 0: mkn_omp, index 1: lib, index 2: mkn_offload, etc.
METHODS=("mkn_omp" "lib" "mkn_offload" "mnk_offload" "blk_offload" "asy_offload" "lib_offload")

for METHOD in "${METHODS[@]}"; do
    echo "# $METHOD" >> $OUTPUT

    # Extract data block for this method
    # Format in output file: Memory(KB) MFLOPS # matmult_xxx
    awk -v method="$METHOD" '
        /^=+$/ { next }
        /Testing: / {
            gsub(/Testing: /, "")
            current = $0
            next
        }
        current == method && /# matmult_/ {
            # $1 = Memory(KB), $2 = MFLOPS
            print $1, $2
        }
    ' "$INPUT" >> $OUTPUT

    # Add blank lines between datasets for gnuplot index
    echo "" >> $OUTPUT
    echo "" >> $OUTPUT
done

echo "Data extracted to $OUTPUT"
echo "Methods extracted: ${METHODS[*]}"
