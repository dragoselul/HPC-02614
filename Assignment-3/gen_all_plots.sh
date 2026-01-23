#!/bin/bash

cd /zhome/b6/0/227603/HPC-02614/Assignment-3

# Generate all performance plots
echo "Extracting data..."
./extract_data.sh

echo "Generating plots..."
gnuplot plot_mkn_omp_lib.gnu
gnuplot plot_mkn_mnk_offload.gnu
gnuplot plot_blk_offload.gnu
gnuplot plot_asy_offload.gnu
gnuplot plot_lib_offload.gnu
gnuplot plot_all_comparison.gnu

echo "Done! Generated plots:"
ls -la *.png
