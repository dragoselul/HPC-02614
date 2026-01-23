#!/usr/bin/gnuplot

set terminal pngcairo size 1400,900 enhanced font 'Verdana,12'
set output 'plot_asy_offload.png'

set title "GPU Offload: Asynchronous Data Transfer Performance" font ",16"
set xlabel "Memory Footprint (MB) [3×N²×8 bytes / 1024]" font ",12"
set ylabel "Performance (GFLOPS)" font ",12"

set logscale x
set logscale y
set grid

set key left top box

# index 5 = asy_offload

plot 'data.txt' index 5 using ($1/1024):($2/1000) with linespoints title 'asy\_offload (async transfer)' lw 2 pt 13 ps 1.5 lc rgb "#e377c2"
