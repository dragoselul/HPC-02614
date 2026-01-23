# plot_gflops_vs_memory.gnu
set terminal pngcairo size 1200,600 enhanced font 'Arial,12'
set output 'gflops_vs_memory.png'
set title "Jacobi GPU Solver: GFLOPS vs Memory Footprint"
set xlabel "Memory Footprint (MB)"
set ylabel "GFLOPS"
set logscale y
set logscale x
set grid
set key top left

# Assuming your data files have columns:
# 1 N  2 Memory_MB  3 Time_s  4 MUpdates_per_s  5 GFLOPS  6 Bandwidth_GB_s

plot 'grid_scaling_GPU_map.dat'     using 2:5 with linespoints pt 7 ps 1.5 lw 2 lc rgb '#0060ad' title 'Map', \
     'grid_scaling_GPU_memcpy.dat' using 2:5 with linespoints pt 7 ps 1.5 lw 2 lc rgb '#dd181f' title 'Memcpy', \
     'grid_scaling_GPU_dual.dat'   using 2:5 with linespoints pt 7 ps 1.5 lw 2 lc rgb '#00aa00' title 'Dual GPU', \
     'grid_scaling_GPU_norm.dat'   using 2:5 with linespoints pt 7 ps 1.5 lw 2 lc rgb '#9933cc' title 'Map with norm'
