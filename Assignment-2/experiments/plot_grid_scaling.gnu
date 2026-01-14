# plot_grid_scaling.gnu - Plot grid size scaling with cache boundaries
#
# Usage: gnuplot plot_grid_scaling.gnu
#

set terminal pngcairo size 1400,600 enhanced font 'Arial,12'
set output 'grid_scaling.png'

set multiplot layout 1,2 title "Poisson 3D Solver: Grid Size Scaling Analysis" font ',14'

# Cache sizes in MB (adjust these for your CPU)
L1_KB = 512        # L1 cache per core in KB
L2_KB = 8192       # L2 cache per core in KB
L3_MB = 32        # L3 cache total in MB

L1_MB = L1_KB / 1024.0
L2_MB = L2_KB / 1024.0

# ============================================
# Plot 1: GFLOPS vs Memory Size
# ============================================
set title "Performance vs Memory Footprint"
set xlabel "Memory Usage (MB)"
set ylabel "GFLOPS"
set logscale x
set grid
set key top right

# Draw vertical lines for cache boundaries
set arrow 1 from L1_MB, graph 0 to L1_MB, graph 1 nohead lc rgb '#ff0000' lw 2 dt 2
set arrow 2 from L2_MB, graph 0 to L2_MB, graph 1 nohead lc rgb '#00aa00' lw 2 dt 2
set arrow 3 from L3_MB, graph 0 to L3_MB, graph 1 nohead lc rgb '#0000ff' lw 2 dt 2

# Labels for cache boundaries
set label 1 "L1" at L1_MB, graph 0.95 center tc rgb '#ff0000' font ',10'
set label 2 "L2" at L2_MB, graph 0.95 center tc rgb '#00aa00' font ',10'
set label 3 "L3" at L3_MB, graph 0.95 center tc rgb '#0000ff' font ',10'

plot 'grid_scaling.dat' using 3:4 with linespoints pt 7 ps 1.5 lw 2 lc rgb '#0060ad' title 'Actual GFLOPS'

# Clear arrows and labels for next plot
unset arrow 1
unset arrow 2
unset arrow 3
unset label 1
unset label 2
unset label 3

# ============================================
# Plot 2: GFLOPS vs Grid Size N
# ============================================
set title "Performance vs Grid Size"
set xlabel "Grid Size (N)"
set ylabel "GFLOPS"
unset logscale x
set grid
set key top right

# Calculate memory for given N: 3 * (N+2)^3 * 8 bytes
# Find N values that correspond to cache boundaries
# Memory = 3 * (N+2)^3 * 8 => N = (Memory / 24)^(1/3) - 2

N_L1 = (L1_KB * 1024.0 / 24.0) ** (1.0/3.0) - 2
N_L2 = (L2_KB * 1024.0 / 24.0) ** (1.0/3.0) - 2
N_L3 = (L3_MB * 1024.0 * 1024.0 / 24.0) ** (1.0/3.0) - 2

set arrow 1 from N_L1, graph 0 to N_L1, graph 1 nohead lc rgb '#ff0000' lw 2 dt 2
set arrow 2 from N_L2, graph 0 to N_L2, graph 1 nohead lc rgb '#00aa00' lw 2 dt 2
set arrow 3 from N_L3, graph 0 to N_L3, graph 1 nohead lc rgb '#0000ff' lw 2 dt 2

set label 1 sprintf("L1 (N≈%.0f)", N_L1) at N_L1, graph 0.95 center tc rgb '#ff0000' font ',10'
set label 2 sprintf("L2 (N≈%.0f)", N_L2) at N_L2, graph 0.95 center tc rgb '#00aa00' font ',10'
set label 3 sprintf("L3 (N≈%.0f)", N_L3) at N_L3, graph 0.95 center tc rgb '#0000ff' font ',10'

plot 'grid_scaling.dat' using 1:4 with linespoints pt 7 ps 1.5 lw 2 lc rgb '#0060ad' title 'Actual GFLOPS'

unset multiplot

# ============================================
# Single detailed plot with memory on secondary axis
# ============================================
set terminal pngcairo size 1000,700 enhanced font 'Arial,12'
set output 'cache_effects.png'

set title "Cache Effects on Jacobi Solver Performance" font ',14'
set xlabel "Grid Size (N)"
set ylabel "GFLOPS" tc rgb '#0060ad'
set y2label "Memory (MB)" tc rgb '#dd181f'
set y2tics
set ytics nomirror
set grid

# Reset arrows
unset arrow 1
unset arrow 2
unset arrow 3
unset label 1
unset label 2
unset label 3

# Cache boundary lines
set arrow 1 from N_L1, graph 0 to N_L1, graph 1 nohead lc rgb '#ff0000' lw 2 dt 2
set arrow 2 from N_L2, graph 0 to N_L2, graph 1 nohead lc rgb '#00aa00' lw 2 dt 2
set arrow 3 from N_L3, graph 0 to N_L3, graph 1 nohead lc rgb '#0000ff' lw 2 dt 2

set label 1 "L1 Cache" at N_L1, graph 0.92 center tc rgb '#ff0000' font ',10'
set label 2 "L2 Cache" at N_L2, graph 0.92 center tc rgb '#00aa00' font ',10'
set label 3 "L3 Cache" at N_L3, graph 0.92 center tc rgb '#0000ff' font ',10'

set key top center

plot 'grid_scaling.dat' using 1:4 with linespoints pt 7 ps 1.5 lw 2 lc rgb '#0060ad' title 'GFLOPS' axes x1y1, \
     'grid_scaling.dat' using 1:3 with linespoints pt 5 ps 1.2 lw 2 lc rgb '#dd181f' title 'Memory (MB)' axes x1y2

print "Plots generated: grid_scaling.png, cache_effects.png"

