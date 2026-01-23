#!/usr/bin/gnuplot

set terminal pngcairo size 1400,900 enhanced font 'Verdana,12'
set output 'plot_blk_offload.png'

set title "GPU Offload: Blocked Algorithm Performance" font ",16"
set xlabel "Memory Footprint (MB) [3×N²×8 bytes / 1024]" font ",12"
set ylabel "Performance (GFLOPS)" font ",12"

set logscale x
set logscale y
set grid

set key left top box

# index 4 = blk_offload

plot 'data.txt' index 4 using ($1/1024):($2/1000) with linespoints title 'blk\_offload (block size=32)' lw 2 pt 11 ps 1.5 lc rgb "#8c564b"
