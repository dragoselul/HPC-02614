# plot_gflops_vs_memory.gnu
set terminal pngcairo size 1200,600 enhanced font 'Arial,12'
set output 'gflops_vs_memory.png'
set title "Jacobi GPU Solver: GFLOPS vs Memory Footprint"
set xlabel "Memory Footprint (GB)"
set ylabel "GFLOPS"
set logscale y
set logscale x
set grid
set key top left

# Compute memory footprint from N using: footprint = 3 * N^3 * 8 / (1024^3) GB
# Assuming your data files have first column N and fifth column GFLOPS
plot 'grid_scaling_GPU_map.dat' using ($1**3*3*8/1024.0/1024.0/1024.0):5 \
     with linespoints pt 7 ps 1.5 lw 2 lc rgb '#0060ad' title 'Map', \
     'grid_scaling_GPU_memcpy.dat' using ($1**3*3*8/1024.0/1024.0/1024.0):5 \
     with linespoints pt 7 ps 1.5 lw 2 lc rgb '#dd181f' title 'Memcpy', \
     'grid_scaling_GPU_dual.dat' using ($1**3*3*8/1024.0/1024.0/1024.0):5 \
     with linespoints pt 7 ps 1.5 lw 2 lc rgb '#00aa00' title 'Dual GPU', \
     'grid_scaling_GPU_norm.dat' using ($1**3*3*8/1024.0/1024.0/1024.0):5 \
     with linespoints pt 7 ps 1.5 lw 2 lc rgb '#9933cc' title 'Map with norm'
