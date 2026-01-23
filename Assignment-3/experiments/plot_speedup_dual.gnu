# plot_speedup_dual.gnu
set terminal pngcairo size 1200,600 enhanced font 'Arial,12'
set output 'dual_gpu_speedup.png'
set title "Dual GPU Speed-up"
set xlabel "Grid Size (N)"
set ylabel "Speed-up (x)"
set grid
set key top left

plot 'dual_speedup.csv' using 1:2 with linespoints pt 7 ps 1.5 lw 2 lc rgb '#0060ad' title 'Map', \
     'dual_speedup.csv' using 1:3 with linespoints pt 7 ps 1.5 lw 2 lc rgb '#dd181f' title 'Memcpy', \
     2 with lines lw 2 lc rgb '#ff0000' dt 2 title 'Ideal 2x'
