#!/usr/bin/gnuplot

set terminal pngcairo size 1200,800 enhanced font 'Verdana,12'
set output 'performance_comparison.png'

set title "Matrix Multiplication Performance Comparison" font ",14"
set xlabel "Matrix Size (m=n=k)"
set ylabel "Performance (GFLOPS)"

set logscale x
set logscale y
set grid

set key left top

# Convert MFLOPS to GFLOPS (divide by 1000)
plot 'data.txt' index 0 using 1:($2/1000) with linespoints title 'mkn\_omp' lw 2 pt 7 ps 1.5, \
     ''         index 1 using 1:($2/1000) with linespoints title 'mkn\_offload' lw 2 pt 5 ps 1.5, \
     ''         index 2 using 1:($2/1000) with linespoints title 'mnk\_offload' lw 2 pt 9 ps 1.5, \
     ''         index 3 using 1:($2/1000) with linespoints title 'blk\_offload' lw 2 pt 11 ps 1.5, \
     ''         index 4 using 1:($2/1000) with linespoints title 'asy\_offload' lw 2 pt 13 ps 1.5, \
     ''         index 5 using 1:($2/1000) with linespoints title 'lib\_offload' lw 2 pt 15 ps 1.5
