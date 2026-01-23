#!/usr/bin/gnuplot

set terminal pngcairo size 1400,900 enhanced font 'Verdana,12'
set output 'plot_mkn_mnk_offload.png'

set title "GPU Offload: Loop Order Comparison (mkn vs mnk)" font ",16"
set xlabel "Memory Footprint (MB) [3×N²×8 bytes / 1024]" font ",12"
set ylabel "Performance (GFLOPS)" font ",12"

set logscale x
set logscale y
set grid

set key left top box

# index 2 = mkn_offload, index 3 = mnk_offload

plot 'data.txt' index 2 using ($1/1024):($2/1000) with linespoints title 'mkn\_offload' lw 2 pt 7 ps 1.5 lc rgb "#d62728", \
     ''         index 3 using ($1/1024):($2/1000) with linespoints title 'mnk\_offload' lw 2 pt 9 ps 1.5 lc rgb "#9467bd"
