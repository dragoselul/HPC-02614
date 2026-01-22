#!/bin/bash
# Extract benchmark data from output file
# Usage: ./extract_data.sh mm_batch_gpu_JOBID.out

if [ $# -eq 0 ]; then
    echo "Usage: $0 <output_file>"
    exit 1
fi

INPUT=$1
OUTPUT="data.txt"

# Clear output file
> $OUTPUT

# Problem sizes
SIZES=(100 200 500 1000 2000 5000 10000)

# Extract data for each method
for METHOD in "mkn_omp" "mkn_offload" "mnk_offload" "blk_offload" "asy_offload" "lib_offload"; do
    echo "# $METHOD" >> $OUTPUT
    
    # Find the testing block and extract performance data
    awk -v method="$METHOD" '
        /Testing: / { if ($2 == method) found=1; else found=0 }
        found && /# matmult_/ { print }
    ' $INPUT | while read -r line; do
        # Extract memory, MFLOPS
        mem=$(echo "$line" | awk '{print $1}')
        mflops=$(echo "$line" | awk '{print $2}')
        if [ ! -z "$mflops" ]; then
            # Determine size based on line number
            idx=$((count++))
            echo "${SIZES[$idx]} $mflops $mem"
        fi
    done | awk '{print $1, $2, $3}' >> $OUTPUT
    
    # Add blank line between datasets
    echo "" >> $OUTPUT
    echo "" >> $OUTPUT
done

echo "Data extracted to $OUTPUT"
echo "Run: gnuplot plot_performance.gnu"
