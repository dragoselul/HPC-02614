#!/usr/bin/gnuplot

set terminal pngcairo size 1400,900 enhanced font 'Verdana,12'
set output 'plot_mkn_omp_lib.png'

set title "CPU Performance: mkn\_omp vs BLAS Library" font ",16"
set xlabel "Memory Footprint (MB) [3×N²×8 bytes / 1024]" font ",12"
set ylabel "Performance (GFLOPS)" font ",12"

set logscale x
set logscale y
set grid

# Move legend to left top to avoid cache line labels
set key left top box

# CPU Cache sizes (total system):
L1_MB = 2
L2_MB = 64
L3_MB = 512

set arrow from L1_MB, graph 0 to L1_MB, graph 1 nohead lc rgb "red" lw 2 dt 2
set arrow from L2_MB, graph 0 to L2_MB, graph 1 nohead lc rgb "orange" lw 2 dt 2
set arrow from L3_MB, graph 0 to L3_MB, graph 1 nohead lc rgb "blue" lw 2 dt 2

set label "L1d (2 MB)" at L1_MB, graph 0.05 center font ",10" tc rgb "red"
set label "L2 (64 MB)" at L2_MB, graph 0.05 center font ",10" tc rgb "orange"
set label "L3 (512 MB)" at L3_MB, graph 0.05 center font ",10" tc rgb "blue"

plot 'data.txt' index 0 using ($1/1024):($2/1000) with linespoints title 'mkn\_omp (OpenMP)' lw 2 pt 7 ps 1.5 lc rgb "#1f77b4", \
     ''         index 1 using ($1/1024):($2/1000) with linespoints title 'lib (CPU BLAS)' lw 2 pt 5 ps 1.5 lc rgb "#2ca02c"
