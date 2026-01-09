#!/usr/bin/gnuplot
# Gnuplot script for Assignment-1 matrix multiplication performance plots
# Usage: gnuplot plot.gnu

set terminal pngcairo size 1200,900 enhanced font 'Arial,12'

# ============ Plot 1: All permutations comparison ============
set output 'permutation_comparison.png'
set title "Matrix Multiplication Performance: Loop Order Comparison"
set xlabel "Matrix Size (N for NxN)"
set ylabel "MFLOPS"
set logscale x
set grid
set key right bottom

plot 'nat.dat' using 1:2 with linespoints lw 2 pt 7 title 'nat', \
     'mnk.dat' using 1:2 with linespoints lw 2 pt 5 title 'mnk', \
     'mkn.dat' using 1:2 with linespoints lw 2 pt 9 title 'mkn', \
     'kmn.dat' using 1:2 with linespoints lw 2 pt 11 title 'kmn', \
     'knm.dat' using 1:2 with linespoints lw 2 pt 13 title 'knm', \
     'nmk.dat' using 1:2 with linespoints lw 2 pt 6 title 'nmk', \
     'nkm.dat' using 1:2 with linespoints lw 2 pt 8 title 'nkm'

# ============ Plot 2: Best vs BLAS library ============
set output 'naive_vs_blas.png'
set title "Naive Implementation vs BLAS Library"
set xlabel "Matrix Size (N for NxN)"
set ylabel "MFLOPS"

plot 'mkn.dat' using 1:2 with linespoints lw 2 pt 7 title 'mkn (best naive)', \
     'lib.dat' using 1:2 with linespoints lw 2 pt 5 title 'BLAS (lib)'

# ============ Plot 3: Block size comparison ============
set output 'block_comparison.png'
set title "Blocked Matrix Multiplication: Block Size Effect"
set xlabel "Matrix Size (N for NxN)"
set ylabel "MFLOPS"

plot 'blk_8.dat' using 1:2 with linespoints lw 2 pt 7 title 'bs=8', \
     'blk_16.dat' using 1:2 with linespoints lw 2 pt 5 title 'bs=16', \
     'blk_32.dat' using 1:2 with linespoints lw 2 pt 9 title 'bs=32', \
     'blk_64.dat' using 1:2 with linespoints lw 2 pt 11 title 'bs=64', \
     'blk_128.dat' using 1:2 with linespoints lw 2 pt 13 title 'bs=128', \
     'blk_256.dat' using 1:2 with linespoints lw 2 pt 6 title 'bs=256'

# ============ Plot 4: All on one (multiplot) ============
set output 'all_performance.png'
set terminal pngcairo size 1400,1000 enhanced font 'Arial,11'
set multiplot layout 2,2 title "Matrix Multiplication Performance Analysis"

# Top left: permutations
set title "Loop Order Comparison"
set xlabel "Matrix Size"
set ylabel "MFLOPS"
set key right bottom font ",9"
plot 'nat.dat' using 1:2 with linespoints lw 2 pt 7 ps 0.7 title 'nat', \
     'mkn.dat' using 1:2 with linespoints lw 2 pt 5 ps 0.7 title 'mkn', \
     'kmn.dat' using 1:2 with linespoints lw 2 pt 9 ps 0.7 title 'kmn', \
     'knm.dat' using 1:2 with linespoints lw 2 pt 11 ps 0.7 title 'knm'

# Top right: lib comparison
set title "Naive vs BLAS"
set key left top
plot 'mkn.dat' using 1:2 with linespoints lw 2 pt 7 title 'mkn', \
     'lib.dat' using 1:2 with linespoints lw 2 pt 5 title 'BLAS'

# Bottom left: block sizes
set title "Block Size Effect"
set key right bottom
plot 'blk_16.dat' using 1:2 with linespoints lw 2 pt 7 title 'bs=16', \
     'blk_32.dat' using 1:2 with linespoints lw 2 pt 5 title 'bs=32', \
     'blk_64.dat' using 1:2 with linespoints lw 2 pt 9 title 'bs=64', \
     'blk_128.dat' using 1:2 with linespoints lw 2 pt 13 title 'bs=128', \
     'blk_256.dat' using 1:2 with linespoints lw 2 pt 6 title 'bs=256'

# Bottom right: best blocked vs lib
set title "Best Blocked vs BLAS"
set key left top
plot 'blk_64.dat' using 1:2 with linespoints lw 2 pt 7 title 'blk (bs=128)', \
     'lib.dat' using 1:2 with linespoints lw 2 pt 5 title 'BLAS'

unset multiplot

print "Plots generated: permutation_comparison.png, native_vs_blas.png, block_comparison.png, all_performance.png"

