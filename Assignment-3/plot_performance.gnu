#!/usr/bin/gnuplot

set terminal pngcairo size 1600,1000 enhanced font 'Verdana,12'
set output 'performance_comparison.png'

set title "Matrix Multiplication: All Methods Performance Comparison" font ",16"
set xlabel "Memory Footprint (MB) [3×N²×8 bytes / 1024]" font ",12"
set ylabel "Performance (GFLOPS)" font ",12"

set logscale x
set logscale y
set grid

set key right bottom box

L1_MB = 2
L2_MB = 64
L3_MB = 512
set arrow from L1_MB, graph 0 to L1_MB, graph 1 nohead lc rgb "gray" lw 1 dt 3
set arrow from L2_MB, graph 0 to L2_MB, graph 1 nohead lc rgb "gray" lw 1 dt 3
set arrow from L3_MB, graph 0 to L3_MB, graph 1 nohead lc rgb "gray" lw 1 dt 3
set label "L1d" at L1_MB, graph 0.02 center font ",8" tc rgb "gray"
set label "L2" at L2_MB, graph 0.02 center font ",8" tc rgb "gray"
set label "L3" at L3_MB, graph 0.02 center font ",8" tc rgb "gray"

plot 'data.txt' index 0 using ($1/1024):($2/1000) with linespoints title 'mkn\\_omp (CPU)' lw 2 pt 7 ps 1.2 lc rgb "#1f77b4", \
     ''         index 1 using ($1/1024):($2/1000) with linespoints title 'lib (CPU BLAS)' lw 2 pt 5 ps 1.2 lc rgb "#2ca02c", \
     ''         index 2 using ($1/1024):($2/1000) with linespoints title 'mkn\\_offload' lw 2 pt 7 ps 1.2 lc rgb "#d62728", \
     ''         index 3 using ($1/1024):($2/1000) with linespoints title 'mnk\\_offload' lw 2 pt 9 ps 1.2 lc rgb "#9467bd", \
     ''         index 4 using ($1/1024):($2/1000) with linespoints title 'blk\\_offload' lw 2 pt 11 ps 1.2 lc rgb "#8c564b", \
     ''         index 5 using ($1/1024):($2/1000) with linespoints title 'asy\\_offload' lw 2 pt 13 ps 1.2 lc rgb "#e377c2", \
     ''         index 6 using ($1/1024):($2/1000) with linespoints title 'lib\\_offload (cuBLAS)' lw 2 pt 15 ps 1.2 lc rgb "#17becf"
