# plot_bw_vs_grid.gnu
set terminal pngcairo size 1200,600 enhanced font 'Arial,12'
set output 'bandwidth_vs_grid.png'
set title "Jacobi GPU Solver: Effective Memory Bandwidth (H100)"
set xlabel "Grid Size (N)"
set ylabel "Bandwidth (GB/s)"
set logscale y
set grid
set key top right

# Peak memory bandwidth of NVIDIA H100 PCIe
PEAK_BW_1GPU = 3350
PEAK_BW_2GPU = 2 * 3350

plot 'grid_scaling_gpu_map.dat' using 1:6 with linespoints pt 7 ps 1.5 lw 2 lc rgb '#0060ad' title 'Map', \
     'grid_scaling_gpu_memcpy.dat' using 1:6 with linespoints pt 7 ps 1.5 lw 2 lc rgb '#dd181f' title 'Memcpy', \
     'grid_scaling_gpu_dual.dat' using 1:6 with linespoints pt 7 ps 1.5 lw 2 lc rgb '#00aa00' title 'Dual GPU', \
     PEAK_BW_1GPU with lines lw 2 lc rgb '#ff0000' dt 2 title 'Peak BW 1 GPU', \
     PEAK_BW_2GPU with lines lw 2 lc rgb '#ff9900' dt 2 title 'Peak BW 2 GPU'
