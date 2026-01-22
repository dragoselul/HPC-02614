# plot_gflops_vs_grid.gnu
set terminal pngcairo size 1200,600 enhanced font 'Arial,12'
set output 'gflops_vs_grid.png'
set title "Jacobi GPU Solver: GFLOPS vs Grid Size"
set xlabel "Grid Size (N)"
set ylabel "GFLOPS"
set logscale y
set grid
set key top left

plot 'grid_scaling_gpu_map.dat' using 1:5 with linespoints pt 7 ps 1.5 lw 2 lc rgb '#0060ad' title 'Map', \
     'grid_scaling_gpu_memcpy.dat' using 1:5 with linespoints pt 7 ps 1.5 lw 2 lc rgb '#dd181f' title 'Memcpy', \
     'grid_scaling_gpu_dual.dat' using 1:5 with linespoints pt 7 ps 1.5 lw 2 lc rgb '#00aa00' title 'Dual GPU'
